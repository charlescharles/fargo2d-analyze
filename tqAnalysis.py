
from argparse import ArgumentParser
import numpy as np

"""
return tuple of secondary r, theta
(r, theta)
"""
def getTrajectory():
    sec = np.loadtxt('bigplanet0.dat')
    date = sec[:, -2]
    _, ixs = np.unique(date, return_index=True)
    sec = sec[ixs]

    secx = sec[:, 1]
    secy = sec[:, 2]
    secr = np.sqrt(np.square(secx)+np.square(secy))[::20]
    sect = np.arctan2(secy, secx)[::20]

    return secr, sect

"""
turn an array of shape (nr) into array of shape (nr, ns)
"""
def azimuthalStack(a, ns):
    return np.array([a] * ns).transpose()


"""
calculate torque density dT/dr(r) for specified azimuthal modes.
parameters (dens, r_sup, r_inf, r_med, theta) have shape (nr, ns)
`modes` has shape (n_modes)
returns array of shape (len(modes), nr)
"""
def computeTorqueDensity(mb, secr, sect, dens, r_med, theta, modes, indirect_term):
    nr, ns = r_med.shape
    n_modes = len(modes)
    psi = theta - sect

    dist = np.sqrt(np.square(r_med) + np.square(secr) - 2. * r_med * secr * np.cos(psi))

    accel = 1./ np.power(dist, 3)
    if indirect_term:
        accel -= 1. / np.power(secr, 3)

    accel = mb * secr * accel

    # specific torque
    spec_tq = r_med * accel * np.sin(psi) * secr / dist

    r_dtheta = 2. * np.pi * r_med / ns

    # dT/dr per cell. shape (nr, ns)
    cell_tq = r_dtheta * dens * spec_tq

    # tile cell_tq along z axis n_modes times
    # resulting array has shape (n_modes, nr, ns)
    cell_tq_mat = np.tile(cell_tq, (n_modes, 1, 1))

    # promote theta to (n_modes, nr, ns)
    theta_mat = np.tile(theta, (n_modes, 1, 1))

    # promote modes to (n_modes, nr, ns)
    modes_mat = np.array([np.zeros((nr, ns), dtype='float')+m for m in modes])

    # summed along ns to give shape (n_modes, nr)
    sine = np.sum(cell_tq_mat * np.sin(modes_mat * theta_mat), axis=2)
    cosine = np.sum(cell_tq_mat * np.cos(modes_mat * theta_mat), axis=2)

    return np.sqrt(np.square(sine) + np.square(cosine))


def computeL(dens, vtheta, r_sup, r_inf, r_med):
    nr, ns = dens.shape

    area = np.pi * (np.square(r_sup) - np.square(r_inf)) / ns

    return np.sum(dens * area * vtheta * r_med)


def computeFargoTorque(mb, secr, sect, dens, r_sup, r_inf, r_med, theta):
    nr, ns = dens.shape

    surf = np.pi * (np.square(r_sup) - np.square(r_inf)) / ns

    mcell = surf * dens

    xcell = r_med * np.cos(theta)
    ycell = r_med * np.sin(theta)

    xb = secr * np.cos(sect)
    yb = secr * np.sin(sect)

    dx = xcell - xb
    dy = ycell - yb

    dist3 = np.power(np.square(dx) + np.square(dy), 1.5)

    fxi = mb * mcell * dx / dist3
    fyi = mb * mcell * dy / dist3

    tq = yb * fxi - xb * fyi

    return tq.sum()


def initvars():
    nr, ns = 438, 574
    radialEdges = np.loadtxt('used_rad.dat')
    n = len(radialEdges)
    r_inf = radialEdges[:n - 1]
    r_sup = radialEdges[1:]

    cube = lambda x: np.power(x, 3)
    square = np.square

    radIntervals = np.multiply((2.0 / 3.0), ((cube(r_sup) - cube(r_inf)) / (square(r_sup) - square(r_inf))))
    thetaIntervals = np.linspace(0, 2*np.pi, 574)
    dr = np.ediff1d(radialEdges)

    r_med = azimuthalStack(radIntervals, ns)
    r_sup = azimuthalStack(r_sup, ns)
    r_inf = azimuthalStack(r_inf, ns)

    theta, r = np.meshgrid(thetaIntervals, radIntervals)

    secr, sectheta = getTrajectory()

    return nr, ns, secr, sectheta, r_inf, r_sup, r_med, theta, dr

def deltaL(dens, vtheta, r_med, dm):
    # avg vtheta at edge
    avgVtheta = np.sum(dens[-1] * vtheta[-1])/np.sum(dens[-1])

    return dm * avgVtheta * r_med[-1][0]

def mass(dens, r_sup, r_inf):
    nr, ns = dens.shape
    surf = np.pi * (np.square(r_sup) - np.square(r_inf)) / ns

    return np.sum(surf * dens)


def main():
    parser = ArgumentParser()
    parser.add_argument('-m', '--binary-mass', nargs='?', default=0.2857, type=float)
    parser.add_argument('-c', '--computation', nargs='?', default='all', type=str)
    parser.add_argument('-e', '--end', nargs='?', default=-1, type=int)
    args = parser.parse_args()

    mb = args.binary_mass
    compute = args.computation
    end = args.end

    print 'using binary mass ' + str(mb)
    print 'computing ' + compute

    if end == -1:
        print 'for all orbits'
    else:
        print 'until orbit ' + str(end)

    nr, ns, secr, sectheta, r_inf, r_sup, r_med, theta, dr = initvars()

    if compute == 'all':
        i = 0
        tqDensityFourier = []
        totalTqDirect = []
        angularMomentum = []
        dL = []
        m0 = mass(np.fromfile('gasdens0.dat').reshape(nr, ns), r_sup, r_inf)
        while True:
            if end > 0 and i > end:
                break

            try:
                dens = np.fromfile('gasdens'+str(i)+'.dat').reshape(nr, ns)
                vtheta = np.fromfile('gasvtheta'+str(i)+'.dat').reshape(nr, ns)
            except IOError:
                break

            tqDensityFourier.append(computeTorqueDensity(mb, secr[i], sectheta[i], dens, r_med, theta, np.arange(11), True))
            directDens = computeTorqueDensity(mb, secr[i], sectheta[i], dens, r_med, theta, [0], False)[0]
            totalTqDirect.append(np.sum(directDens * dr))
            angularMomentum.append(computeL(dens, vtheta, r_sup, r_inf, r_med))

            m1 = mass(dens, r_sup, r_inf)
            dL.append(deltaL(dens, vtheta, r_med, m1-m0))
            m0 = m1

            if i%100 == 0:
                print i
            i += 1

        print 'finished at ' + str(i)
        print 'saving'
        np.save('parsedDiagnostics/tqFourier', tqDensityFourier)
        np.save('parsedDiagnostics/tqDirect', totalTqDirect)
        np.save('parsedDiagnostics/angularMomentum', angularMomentum)
        np.save('parsedDiagnostics/deltaL', dL)
        return

    elif compute == 'fargo':
        i = 0
        fargoTq = []
        while True:
            try:
                dens = np.fromfile('gasdens'+str(i)+'.dat').reshape(nr, ns)
            except IOError:
                print 'finished at ' + str(i)
                break

            fargoTq.append(computeFargoTorque(mb, secr[i], sectheta[i], dens, r_sup, r_inf, r_med, theta))

            if i%100 == 0:
                print i
            i += 1

        print 'saving'
        np.save('parsedDiagnostics/fargoTq', fargoTq)
        return

    elif compute == 'deltal':
        i = 0
        dL = []
        m0 = mass(np.fromfile('gasdens0.dat').reshape(nr, ns), r_sup, r_inf)
        while True:
            try:
                dens = np.fromfile('gasdens'+str(i)+'.dat').reshape(nr, ns)
                vtheta = np.fromfile('gasvtheta'+str(i)+'.dat').reshape(nr, ns)
            except IOError:
                print 'finished at ' + str(i)
                break

            m1 = mass(dens, r_sup, r_inf)
            dL.append(deltaL(dens, vtheta, r_med, m1-m0))
            m0 = m1

            if i%100 == 0:
                print i
            i += 1

        print 'saving'
        np.save('parsedDiagnostics/deltaL', dL)
        return

    elif compute == 'angularmom':
        i = 0
        angularMomentum = []
        while True:
            if end > 0 and i > end:
                break

            try:
                dens = np.fromfile('gasdens'+str(i)+'.dat').reshape(nr, ns)
                vtheta = np.fromfile('gasvtheta'+str(i)+'.dat').reshape(nr, ns)
            except IOError:
                break
            angularMomentum.append(computeL(dens, vtheta, r_sup, r_inf, r_med))

            if i%100 == 0:
                print i
            i += 1

        print 'finished at ' + str(i)
        print 'saving'
        np.save('parsedDiagnostics/angularMomentum', angularMomentum)
        return


if __name__ == '__main__':
    main()
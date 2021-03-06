import numpy as np

def main():
    nr, ns = 438, 574
    radialEdges = np.loadtxt('used_rad.dat')
    n = len(radialEdges)
    r0 = radialEdges[:n - 1]
    r1 = radialEdges[1:]

    cube = lambda x: np.power(x, 3)
    square = np.square

    radIntervals = np.multiply((2.0 / 3.0), ((cube(r1) - cube(r0)) / (square(r1) - square(r0))))
    thetaIntervals = np.linspace(0, 2*np.pi, 574)
    rdiff = np.ediff1d(radialEdges)

    rmat = np.array([radIntervals] * ns).transpose()
    drmat = np.array([rdiff]*ns).transpose()

    theta, r = np.meshgrid(thetaIntervals, radIntervals)

    sec = np.loadtxt('bigplanet0.dat')
    secx = sec[:, 1]
    secy = sec[:, 2]
    secr = np.sqrt(np.square(secx)+np.square(secy))[::20]
    sect = np.arctan2(secy, secx)[::20]

    masslost = sec[:, -3]
    masslost = [sum(masslost[current: current+5]) for current in xrange(0, len(masslost), 5)]

    i = 0
    momLost = []
    while True:
        try:
            dens = np.fromfile('gasdens'+str(i)+'.dat').reshape(nr, ns)
            vr = np.fromfile('gasvrad'+str(i)+'.dat').reshape(nr, ns)
            vtheta = np.fromfile('gasvtheta'+str(i)+'.dat').reshape(nr, ns)
        except IOError:
            print 'finished at ' + str(i)
            break

        specMom = specificAngMom(secr[i], sect[i], r, theta, vr, vtheta)
        avgSpecMom = np.sum(specMom[0]*dens[0]) / np.sum(dens[0])

        momLost.append( avgSpecMom * masslost[i] * radIntervals[0] * rdiff[0] )

        if i%100 == 0:
            print i
        i += 1

    print 'saving'
    np.save('parsedDiagnostics/momLostInner', momLost)

def specificAngMom(secr, sect, r, theta, vr, vtheta):
    m = 0.2857
    psi = theta - sect

    a = m/(1.+m)
    c = np.sqrt(np.square(a) + np.square(r) - 2. * a * r * np.cos(psi))
    alpha = np.arcsin(a * np.sin(psi) / c)

    L1 = vtheta * np.cos(alpha)
    L2 = -1. * vr * np.sin(alpha)

    return L1 + L2

if __name__ == '__main__':
    main()
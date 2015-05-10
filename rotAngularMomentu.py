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

    i = 0
    angularMom = []
    while True:
        try:
            dens = np.fromfile('gasdens'+str(i)+'.dat').reshape(nr, ns)
            vtheta = np.fromfile('gasvtheta'+str(i)+'.dat').reshape(nr, ns)
        except IOError:
            print 'finished at ' + str(i)
            break

        angularMom.append(np.sum(dens * vtheta * rmat * drmat))
        if i%100 == 0:
            print i
        i += 1

    print 'saving'
    np.save('parsedDiagnostics/angularMomRot', angularMom)

if __name__ == '__main__':
    main()
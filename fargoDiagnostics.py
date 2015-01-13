import numpy as np
import math


def _thetaBroadcast(row, numRows, numFrames):
    """
    broadcast a theta row into a numFrames x numRows x len(row) array
    """

    arr = np.vstack([row] * numRows)
    return np.array([arr] * numFrames)


def _radialBroadcast(row, numCols, numFrames):
    """
    broadcast a radial row into a numFrames x len(row) x numCols array
    """

    col = np.array(row).reshape([-1, 1])
    arr = np.hstack([col] * numCols)
    return np.array([arr] * numFrames)


def _azimuthalMassAverage(arr, density):
    """
    compute and return the azimuthal mass-weighted average of `arr`
    """
    weightedSum = np.einsum("abc,abc->ab", arr, density)
    radialDensity = np.sum(density, 2)
    return np.divide(weightedSum, radialDensity)


def _lubowDiagnostics(radialIntervals, thetaIntervals, dens, vr, vtheta):
    numRadialIntervals = len(radialIntervals)
    numThetaIntervals = len(thetaIntervals)
    numTimeIntervals = len(vr)

    thetaDeltas = np.zeros_like(thetaIntervals) + (2 * math.pi)/numThetaIntervals

    # numTimeIntervals x numRadialIntervals x numThetaIntervals
    r = _radialBroadcast(radialIntervals, numThetaIntervals, numTimeIntervals)
    theta = _thetaBroadcast(thetaIntervals, numRadialIntervals, numTimeIntervals)
    dtheta = _thetaBroadcast(thetaDeltas, numRadialIntervals, numTimeIntervals)

    # vtheta/r averaged azimuthally
    omega = _azimuthalMassAverage(np.divide(vtheta, r), dens)

    vsin = np.multiply(np.multiply(vtheta, dtheta), np.sin(theta)).sum(2) / math.pi
    vcos = np.multiply(np.multiply(vtheta, dtheta), np.cos(theta)).sum(2) / math.pi

    # numTimeIntervals x numRadialIntervals
    r2d = np.array([radialIntervals] * numTimeIntervals)

    e = (2.0 / np.multiply(r2d, omega)) * np.sqrt(np.add(np.square(vsin), np.square(vcos)))
    peri = np.arctan2(vsin, vcos)

    return {
        "radialEccLubow": e,
        "radialPeriLubow": peri
    }


def _computeCellDiagnostics(radialIntervals, thetaIntervals, vr, vtheta):
    numRadialIntervals = len(radialIntervals)
    numThetaIntervals = len(thetaIntervals)
    numTimeIntervals = len(vr)

    # numTimeIntervals x numRadialIntervals x numThetaIntervals
    r = _radialBroadcast(radialIntervals, numThetaIntervals, numTimeIntervals)
    
    theta = _thetaBroadcast(thetaIntervals, numRadialIntervals, numTimeIntervals)
    
    # r * v_theta
    r_vtheta = np.multiply(r, vtheta)

    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)

    e_x = np.subtract(
            np.multiply(
                r_vtheta,
                np.add(
                    np.multiply(vr, sin_theta),
                    np.multiply(vtheta, cos_theta))),
            cos_theta)

    e_y = np.subtract(
            np.multiply(
                r_vtheta,
                np.subtract(
                    np.multiply(vtheta, sin_theta),
                    np.multiply(vr, cos_theta))),
            sin_theta)

    cellEccentricity = np.sqrt(np.add(np.square(e_x), np.square(e_y)))
    cellPeriastron = np.arctan2(e_y, e_x)

    return {
        "cellEccentricity": cellEccentricity,
        "cellPeriastron": cellPeriastron
        }


def diskMassAverage(arr, density, radialIntervals, numThetaIntervals):
    """
    return average of `arr` weighted by density
    """
    if len(arr.shape) == 3:
        numUsedTimeIntervals = len(arr)
    else:
        arr = np.array([arr])
        numUsedTimeIntervals = 1

    arr = arr[:, 1:, :]
    density = density[:, 1:, :]

    delta_r = np.ediff1d(radialIntervals)

    delta_r_mat = _radialBroadcast(delta_r, numThetaIntervals, numUsedTimeIntervals)
    r_mat = _radialBroadcast(radialIntervals[1:], numThetaIntervals, numUsedTimeIntervals)

    r_delta_r = np.multiply(delta_r_mat, r_mat)
    weightedArr = np.multiply(arr, density)

    weightedSum = np.multiply(r_delta_r, weightedArr).sum(1).sum(1)
    totalMass = np.multiply(r_delta_r, density).sum(1).sum(1)

    return np.divide(weightedSum, totalMass)

def radialDiskMassAverage(arr, dens, radialEdges, radialIntervals, numThetaIntervals):
    numUsedTimeIntervals = len(dens)

    d_theta = 2.*math.pi / numThetaIntervals

    sumRadialDens = dens.sum(2) * d_theta
    sumRadialWeights = np.multiply(sumRadialDens, arr)

    delta_r = np.ediff1d(radialEdges)

    dr_mat = np.array([delta_r] * numUsedTimeIntervals)
    r_mat = np.array([radialIntervals] * numUsedTimeIntervals)

    r_dr = np.multiply(dr_mat, r_mat)

    return np.multiply(r_dr, sumRadialWeights).sum(1)

def computeTotalMass(dens, radialEdges, radialIntervals, numThetaIntervals):
    numUsedTimeIntervals = len(dens)
    
    d_theta = 2.*math.pi / numThetaIntervals
    
    sumRadialDens = dens.sum(2) * d_theta
    
    delta_r = np.ediff1d(radialEdges)
    
    dr_mat = np.array([delta_r] * numUsedTimeIntervals)
    r_mat = np.array([radialIntervals] * numUsedTimeIntervals)
    
    r_dr = np.multiply(dr_mat, r_mat)
    
    return np.multiply(r_dr, sumRadialDens).sum(1)

def computeDiagnostics(radialEdges, radialIntervals, thetaIntervals, dens, vr, vtheta):
    diags = _computeCellDiagnostics(radialIntervals, thetaIntervals, vr, vtheta)
    radialEccMK = _azimuthalMassAverage(diags['cellEccentricity'], dens)
    radialPeriMK = _azimuthalMassAverage(diags['cellPeriastron'], dens)

    lubowDiagnostics = _lubowDiagnostics(radialIntervals, thetaIntervals, dens, vr, vtheta)
    radialEccLubow = lubowDiagnostics['radialEccLubow']
    radialPeriLubow = lubowDiagnostics['radialPeriLubow']

    numThetaIntervals = len(thetaIntervals)
    diskEccMK = diskMassAverage(diags['cellEccentricity'], dens, radialIntervals, numThetaIntervals)
    diskPeriMK = diskMassAverage(diags['cellPeriastron'], dens, radialIntervals, numThetaIntervals)

    radialDens = 2.0 * math.pi / numThetaIntervals * dens.sum(2)

    totalMass = computeTotalMass(dens, radialEdges, radialIntervals, numThetaIntervals)

    diskEccLubow = radialDiskMassAverage(radialPeriLubow, dens, radialEdges, radialIntervals, numThetaIntervals)
    diskPeriLubow = radialDiskMassAverage(radialPeriLubow, dens, radialEdges, radialIntervals, numThetaIntervals)

    diskEccMK = radialDiskMassAverage(radialEccMK, dens, radialEdges, radialIntervals, numThetaIntervals)

    m = radialDiskMassAverage(1.0, dens, radialEdges, radialIntervals, numThetaIntervals)
    print "totalmass: " + str(totalMass)
    print "avged mass: " + str(m)

    return {
        "radialEccMK": radialEccMK,
        "radialPeriMK": radialPeriMK,
        "radialEccLubow": radialEccLubow,
        "radialPeriLubow": radialPeriLubow,
        "radialDens": radialDens,
        "diskEccMK": diskEccMK,
        "diskPeriMK": diskPeriMK,
        "diskEccLubow": diskEccLubow,
        "diskPeriLubow": diskPeriLubow,
        "totalMass": totalMass
    }
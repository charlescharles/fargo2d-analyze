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
    print arr.shape
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


def radialDiskMassAverage(arr, radialDens, radialIntervals, totalMass):
    arr = arr[:, 1:]
    radialDens = radialDens[:, 1:]

    numTimeIntervals = len(radialDens)

    delta_r = np.ediff1d(radialIntervals)
    delta_r_mat = np.array(delta_r * numTimeIntervals)
    r_mat = np.array(radialIntervals[1:] * numTimeIntervals)
    r_delta_r = np.multiply(delta_r_mat, r_mat)

    weightedArr = np.multiply(np.multiply(np.multiply(arr, radialDens), 2 * math.pi), r_mat)

    weightedSum = np.multiply(weightedArr, delta_r).sum(1)

    return np.divide(weightedSum, totalMass)


def computeTotalMass(dens, radialIntervals, numThetaIntervals):
    dens = dens[:, 1:, :]
    
    numUsedTimeIntervals = len(dens)

    delta_r = np.ediff1d(radialIntervals)

    delta_r_mat = _radialBroadcast(delta_r, numThetaIntervals, numUsedTimeIntervals)
    r_mat = _radialBroadcast(radialIntervals[1:], numThetaIntervals, numUsedTimeIntervals)

    r_delta_r = np.multiply(delta_r_mat, r_mat)

    totalMass = np.multiply(np.multiply(r_delta_r, dens).sum(1).sum(1), 2 * math.pi)

    return totalMass


def computeDiagnostics(radialIntervals, thetaIntervals, dens, vr, vtheta):
    diags = _computeCellDiagnostics(radialIntervals, thetaIntervals, vr, vtheta)
    radialEccMK = _azimuthalMassAverage(diags['cellEccentricity'], dens)
    radialPeriMK = _azimuthalMassAverage(diags['cellPeriastron'], dens)

    lubowDiagnostics = _lubowDiagnostics(radialIntervals, thetaIntervals, dens, vr, vtheta)
    radialEccLubow = lubowDiagnostics['radialEccLubow']
    radialPeriLubow = lubowDiagnostics['radialPeriLubow']

    diskEccMK = diskMassAverage(diags['cellEccentricity'], dens, radialIntervals, len(thetaIntervals))
    diskPeriMK = diskMassAverage(diags['cellPeriastron'], dens, radialIntervals, len(thetaIntervals))

    radialDens = np.average(dens, 2)

    totalMass = computeTotalMass(dens, radialIntervals, len(thetaIntervals))

    diskEccLubow = radialDiskMassAverage(radialEccLubow, radialDens, radialIntervals, totalMass)
    diskPeriLubow = radialDiskMassAverage(radialPeriLubow, radialDens, radialIntervals, totalMass)

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
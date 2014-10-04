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


def _lubowEccentricity(radialIntervals, thetaIntervals, dens, vr, vtheta):
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

    print vsin.shape

    e = (2.0 / np.multiply(r2d, omega)) * np.sqrt(np.add(np.square(vsin), np.square(vcos)))

    return e


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


def computeRadialDiagnostics(radialIntervals, thetaIntervals, dens, vr, vtheta):
    cellEccMK, cellPeri = _computeCellDiagnostics(radialIntervals, thetaIntervals, vr, vtheta)

    radialEccMK = _azimuthalMassAverage(cellEccMK, dens)
    radialPeriMK = _azimuthalMassAverage(cellPeri, dens)
    radialEccLubow = _lubowEccentricity(radialIntervals, thetaIntervals, dens, vr, vtheta)

    radialDens = np.average(dens, 2)

    return {
        "radialEccMK": radialEccMK,
        "radialPeriMK": radialPeriMK,
        "radialEccLubow": radialEccLubow,
        "radialDens": radialDens
    }
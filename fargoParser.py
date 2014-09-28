__author__ = 'cguo'

import numpy as np
import glob
import re
import itertools as it
import math

class FargoParser:
    """
    a parser for Fargo2D output files
    """

    def expandedCycle(self, iterable, n):
        """
        cycle the expanded iterable, where the expanded iterable is the iterable with each element repeated n times
        e.g., expandedCycle([A, B, C], 2) => [A, A, B, B, C, C, A, A, B, B, C, C, ...]
        """
        saved = []
        for item in iterable:
            saved.append(item)
            for _ in range(n):
                yield item

        while saved:
            for item in saved:
                for _ in range(n):
                    yield item


    def __init__(self, outputDir):
        if outputDir.endswith('/'):
            outputDir = outputDir[:-1]

        self.outputDir = outputDir
        self.gasdens = None
        self.gasvrad = None
        self.gasvtheta = None

        self.G = 1.0
        self.M = 1.0

        self.initParams()
        self.parseGasOutput()


    def initParams(self):

        dims = np.loadtxt(self.pathTo("dims.dat"))
        self.rmax = dims[4]
        self.numOutputs = dims[5]
        self.Nrad = int(dims[6])
        self.Nsec = int(dims[7])

        self.gasradii = np.loadtxt(self.pathTo("used_rad.dat"))[:self.Nrad]

        dtheta = 2 * math.pi/self.Nsec
        self.gasthetas = np.arange(0, 2*math.pi + (dtheta/2), dtheta)[:self.Nsec]


        planetData = np.loadtxt(self.pathTo("orbit0.dat"))
        self.times = planetData[:, 0]
        self.Ntime = len(self.times)


    def pathTo(self, endpoint):
        return self.outputDir + "/" + endpoint


    def extractFileIndex(self, filePath):
        name = filePath.split('/')[-1]
        return int(re.search('[0-9]+', name).group(0))


    def parseGasValue(self, varType):
        """
        :param varType: "dens", "vrad", or "vtheta"
        :return:
        """

        fileFormat = "gas{}*.dat".format(varType)

        filePaths = glob.glob(self.pathTo(fileFormat))
        sortedPaths = sorted(filePaths, key=self.extractFileIndex)

        arrays = []
        for path in sortedPaths:
            arr = np.fromfile(path, dtype='double')
            arr.shape = (self.Nrad, self.Nsec)
            arrays.append(arr)

        return np.array(arrays)


    def parseGasOutput(self):
        """
        read and parse gas density, vrad, vtheta and set values on self
        :return:
        """

        self.gasdens = self.parseGasValue("dens")
        self.gasvrad = self.parseGasValue("vrad")
        self.gasvtheta = self.parseGasValue("vtheta")

        # number of actual outputs
        self.NtimesUsed = len(self.gasvtheta)


    def computeCellEccentricities(self):
        """
        compute eccentricity matrix
        :return:
        """

        Nsec = self.Nsec
        Nrad = self.Nrad
        Ntime = self.Ntime

        G = self.G
        M = self.M

        sin = math.sin
        cos = math.cos
        sqrt = math.sqrt
        pow = math.pow

        # increment time every Nsec * Nrad values
        ts = self.expandedCycle(self.times, Nsec * Nrad)

        # increment radius every Nrad values
        rs = self.expandedCycle(self.gasradii, Nsec)

        # increment theta every value
        thetas = it.cycle(self.gasthetas)

        # iterate through linearly
        vrs = np.nditer(self.gasvrad)
        vthetas = np.nditer(self.gasvtheta)

        # x, y components of eccentricity (Mueller, Kley 2012)
        # such that e_vec = ex * xhat + ex * yhat
        exs = []
        eys = []

        # magnitude of eccentricity
        # e = sqrt(ex^2 + ey^2)
        es = []
        for (r, theta, t, vr, vtheta) in zip(rs, thetas, ts, vrs, vthetas):
            ex = ((r * vtheta) * (vr * sin(theta) + vtheta * cos(theta)) / (G * M)) - cos(theta)
            ey = ((r * vtheta) / (G * M)) * (vtheta * sin(theta) - vr * cos(theta)) - sin(theta)
            #ey = sin(theta) - (r * vtheta) * (vr * cos(theta) - vtheta * sin(theta)) / (G * M)
            e = sqrt(pow(ex, 2) + pow(ey, 2))

            exs.append(ex)
            eys.append(ey)
            es.append(e)

        # reshape eccentricity matrices
        NtimeCalculated = len(exs) / (Nrad * Nsec)
        shape = [NtimeCalculated, Nrad, Nsec]
        exs, eys, es = map(np.array, [exs, eys, es])
        exs.shape = eys.shape = es.shape = shape

        # self.exs = exs
        # self.eys = eys
        self.cellEccentricities = es

    def save(self, array, filename):
        """
        saves array to filename in working directory
        :param array:
        :param filename:
        :return:
        """

        np.save(self.pathTo(filename), array)


    def weightedAverage(self, arr, axis):
        """
        NOTE: if we implement logarithmic scaling of Nrad we'll have to change this!
        compute and return the density-weighted average of `arr` along axis = {'r', 'theta'}
        :param arr:
        :return:
        """
        density = self.gasdens
        if axis == 'theta':
            weightedSum = np.einsum("abc,abc->ab", arr, density)
            radialDensity = np.sum(density, 2)
            return np.divide(weightedSum, radialDensity)

        elif axis == 'r':
            Nrad = self.Nrad
            NtimesUsed = self.NtimesUsed

            # differences between consecutive elements
            delta_r = np.ediff1d(self.gasradii)

            col = delta_r.reshape(-1, 1)
            arr = np.hstack([col] * Nrad)
            # 3-D array of dimensions (NtimesUsed, Nsec, Nrad)
            # matching dimensions of gas variables
            delta_r = np.array([arr] * NtimesUsed)

            # element-wise multiply r-steps, the array, and density
            # and sum over axis 1 (radius)
            weightedSum = np.multiply(np.multiply(delta_r, arr), density).sum(1)
            azimuthalDensity = np.multiply(delta_r, density).sum(1)

            return np.divide(weightedSum, azimuthalDensity)
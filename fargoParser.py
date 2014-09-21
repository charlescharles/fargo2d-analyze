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


    def initParams(self):
        self.gasradii = np.loadtxt(self.pathTo("used_rad.dat"))

        dims = np.loadtxt(self.pathTo("dims.dat"))
        self.rmax = dims[4]
        self.numOutputs = dims[5]
        self.Nrad = dims[6]
        self.Nsec = dims[7]

        dtheta = math.pi/self.nrad
        self.gasthetas = np.arange(0, 2*math.pi + (dtheta/2), dtheta)


        planetData = np.loadtxt("orbit0.dat")
        self.times = planetData[:, 0]
        self.Ntime = len(self.times)


    def pathTo(self, endpoint):
        return self.outputDir + endpoint


    def extractFileIndex(self, filePath):
        name = filePath.split('/')[-1]
        return int(re.search('[0-9]+', name).group(0))


    def parseGasValue(self, type):
        """
        :param type: "dens", "vrad", or "vtheta"
        :return:
        """

        fileFormat = "gas{}[0-9]+.dat".format(type)

        filePaths = glob.glob(self.pathTo(fileFormat))
        sortedPaths = sorted(filePaths, key=self.extractFileIndex)

        arrays = []
        for path in sortedPaths:
            arrays.append(np.fromfile(path, dtype='double'))

        return np.dstack(arrays)


    def parseGasOutput(self):
        """
        read and parse gas density, vrad, vtheta and set values on self
        :return:
        """

        self.gasdens = self.parseGasValue("dens")
        self.gasvrad = self.parseGasValue("vrad")
        self.gasvtheta = self.parseGasValue("vtheta")


    def cellEccentricity(self):
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
        rs = self.expandedCycle(self.gasradii, Nrad)

        # increment theta every value
        thetas = it.cycle(self.gasthetas)

        vrs = self.gasvrad
        vthetas = self.gasvtheta

        # x, y components of eccentricity (Mueller, Kley 2012)
        # such that e_vec = ex * xhat + ex * yhat
        exs = []
        eys = []

        # magnitude of eccentricity
        # e = sqrt(ex^2 + ey^2)
        es = []
        for (r, theta, t, vr, vtheta) in zip(rs, thetas, ts, vrs, vthetas):
            ex = (r * vtheta) * (vr * sin(theta) + vtheta * cos(theta)) / (G * M) - cos(theta)
            ey = sin(theta) - (r * vtheta) * (vr * cos(theta) - vtheta * sin(theta)) / (G * M)
            e = sqrt(pow(ex, 2) + pow(ey, 2))

            exs.append(ex)
            eys.append(ey)
            es.append(e)

        # reshape eccentricity matrices
        shape = [Nrad, Nsec, Ntime]
        exs.shape = eys.shape = es.shape = shape


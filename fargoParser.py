__author__ = 'cguo'

import numpy as np
import glob
import re
import itertools as it
import math
import os
import sys
import pickle
import logging
import gc

class FargoParser:
    """
    a parser for Fargo2D output files
    """

    logging.basicConfig(level=logging.DEBUG, filename='parserDiagnostics.log', filemode='a')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    logging.getLogger('').addHandler(console)


    def __init__(self, outputDir, batchSize=100):
        if outputDir.endswith('/'):
            outputDir = outputDir[:-1]

        self.outputDir = outputDir

        self.batchSize = batchSize

        self.G = 1.0
        self.M = 1.0

        self.initParams()
        self.saveRunParameters()
        self.run()


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


    def initParams(self):
        logging.info("\n*** reading run parameters ***\n")

        dims = np.loadtxt(self.pathTo("dims.dat"))
        self.maxRadius = dims[4]

        # Nrad
        self.numRadialIntervals = int(dims[6])

        # Nsec
        self.numThetaIntervals = int(dims[7])

        self.radialIntervals = np.loadtxt(self.pathTo("used_rad.dat"))[:self.numRadialIntervals]

        dtheta = 2 * math.pi/self.numThetaIntervals
        self.thetaIntervals = np.arange(0, 2*math.pi + (dtheta/2), dtheta)[:self.numThetaIntervals]

        planetData = np.loadtxt(self.pathTo("orbit0.dat"))
        self.timeIntervals = planetData[:, 0]
        self.numTimeIntervals = len(self.timeIntervals)

        filePaths = glob.glob(self.pathTo('gasdens*.dat'))
        self.totalNumOutputs = len(filePaths)


    def pathTo(self, endpoint):
        return self.outputDir + "/" + endpoint


    def extractFileIndex(self, filePath):
        name = filePath.split('/')[-1]
        return int(re.search('[0-9]+', name).group(0))


    def parseGasValue(self, varType, startIndex, endIndex):
        """
        :param varType: "dens", "vrad", or "vtheta"
        :return:
        """

        logging.info("\n*** parsing values for gas" + varType + " ***\n")

        fileFormat = "gas" + varType + "*.dat"

        filePaths = glob.glob(self.pathTo(fileFormat))
        sortedPaths = sorted(filePaths, key=self.extractFileIndex)

        arrays = []
        for path in sortedPaths[startIndex : endIndex]:
            arr = np.fromfile(path, dtype='double')
            arr.shape = (self.numRadialIntervals, self.numThetaIntervals)
            arrays.append(arr)

            logging.info("parsed gas output file: " + path.split('/')[-1])

        ret = np.array(arrays)

        return ret


    def parseGasOutput(self, startIndex, endIndex):
        """
        read and parse gas density, vrad, vtheta and set values on self
        :return: tuple of (gasdens, gasvrad, gasvtheta)
        """

        varTypes = ['dens', 'vrad', 'vtheta']

        return (self.parseGasValue(varType, startIndex, endIndex) for varType in varTypes)


    def computeCellEccentricities(self, vrad, vtheta, startIndex, endIndex):
        """
        compute eccentricity matrix
        :return:
        """

        logging.info("\n*** computing cell eccentricities ***\n")

        numThetaIntervals = self.numThetaIntervals
        numRadialIntervals = self.numRadialIntervals

        G = self.G
        M = self.M

        sin = math.sin
        cos = math.cos
        sqrt = math.sqrt
        pow = math.pow

        radialIntervals = self.radialIntervals
        thetaIntervals = self.thetaIntervals

        # increment radius every Nrad values
        rs = self.expandedCycle(radialIntervals, numThetaIntervals)

        # increment theta every value
        thetas = it.cycle(thetaIntervals)

        # iterate through linearly
        vrs = np.nditer(vrad)
        vthetas = np.nditer(vtheta)

        # magnitude of eccentricity
        # e = sqrt(ex^2 + ey^2)
        es = []

        logging.info("\n*** cycling through r, theta, vr, vtheta ***\n")

        for (r, theta, vr, vtheta) in zip(rs, thetas, vrs, vthetas):
            ex = ((r * vtheta) * (vr * sin(theta) + vtheta * cos(theta)) / (G * M)) - cos(theta)
            ey = ((r * vtheta) / (G * M)) * (vtheta * sin(theta) - vr * cos(theta)) - sin(theta)
            e = sqrt(pow(ex, 2) + pow(ey, 2))

            es.append(e)

            logging.info("r = " + str(r) + "; theta = " + str(theta) +\
                  "; vr = " + str(vr) + "; vtheta = " + str(vtheta) + "; ecc = " + str(e))

        # reshape eccentricity matrix
        numUsedTimeIntervals = len(es) / (numRadialIntervals * numThetaIntervals)
        es = np.array(es).reshape([numUsedTimeIntervals, numRadialIntervals, numThetaIntervals])

        return es

    def saveRunParameters(self):
        # outputDir = os.path.dirname(self.pathTo('parsedOutput'))
        # if not os.path.exists(outputDir):
        #     os.mkdir(outputDir)

        logging.info("\n*** saving run parameters ***\n")

        paramNames = ['numRadialIntervals', 'numThetaIntervals', 'radialIntervals', 'thetaIntervals',\
                      'timeIntervals', 'maxRadius', 'totalNumOutputs']

        params = dict((paramName, getattr(self, paramName)) for paramName in paramNames)

        logging.info("run params:")
        logging.info(str(params))

        with open(self.pathTo('parsed_parameters.pickle'), 'w') as f:
            pickle.dump(params, f)


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
            numRadialIntervals = self.numRadialIntervals
            numUsedTimeIntervals = len(arr)

            # differences between consecutive elements
            delta_r = np.ediff1d(self.gasradii)

            col = delta_r.reshape(-1, 1)
            arr = np.hstack([col] * numRadialIntervals)

            # 3-D array of dimensions (NtimesUsed, Nsec, Nrad)
            # matching dimensions of gas variables
            delta_r = np.array([arr] * numUsedTimeIntervals)

            # element-wise multiply r-steps, the array, and density
            # and sum over axis 1 (radius)
            weightedSum = np.multiply(np.multiply(delta_r, arr), density).sum(1)
            azimuthalDensity = np.multiply(delta_r, density).sum(1)

            return np.divide(weightedSum, azimuthalDensity)


    def runBatch(self, startIndex, endIndex):
        logging.info("\nrunning batch from " + str(startIndex) + " to " + str(endIndex) + "\n")

        dens, vrad, vtheta = self.parseGasOutput(startIndex, endIndex)

        logging.info("\ncomputing cell ecc from " + str(startIndex) + " to " + str(endIndex) + "\n")

        cellEcc = self.computeCellEccentricities(vrad, vtheta, startIndex, endIndex)

        vars = ['dens', 'vrad', 'vtheta', 'cellEcc']

        for var in vars:
            logging.info("\nsaving gas var " + var + " from " + str(startIndex) + " to " + str(endIndex) + "\n")
            np.save(self.pathTo('parsed_' + var + str(startIndex) + '-' + str(endIndex)), eval(var))

        dens = vrad = vtheta = cellEcc = None

        gc.collect()


    def run(self):
        batchSize = self.batchSize
        totalNumOutputs = self.totalNumOutputs
        for start in range(0, totalNumOutputs, batchSize):
            self.runBatch(start, start + batchSize)


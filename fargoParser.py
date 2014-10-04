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

    def thetaBroadcast(self, row, numRows, numFrames):
        """
        broadcast a theta row into a numFrames x numRows x len(row) array
        """

        arr = np.vstack([row] * numRows)
        return np.array([arr] * numFrames)


    def radialBroadcast(self, row, numCols, numFrames):
        """
        broadcast a radial row into a numFrames x len(row) x numCols array
        """

        col = np.array(row).reshape([-1, 1])
        arr = np.hstack([col] * numCols)
        return np.array([arr] * numFrames)


    def diskMassAverage(self, arr, density, params):
        """
        return average of `arr` weighted by density
        """
        numThetaIntervals = params['numThetaIntervals']
        numUsedTimeIntervals = len(arr)

        arr = arr[:, 1:, :]
        density = density[:, 1:, :]

        radialIntervals = params['radialIntervals']

        delta_r = np.ediff1d(params['radialIntervals'])

        delta_r_mat = self.radialBroadcast(delta_r, numThetaIntervals, numUsedTimeIntervals)
        r_mat = self.radialBroadcast(radialIntervals[1:], numThetaIntervals, numUsedTimeIntervals)

        r_delta_r = np.multiply(delta_r_mat, r_mat)
        weightedArr = np.multiply(arr, density)

        weightedSum = np.multiply(r_delta_r, weightedArr).sum(1).sum(1)

        totalMass = np.multiply(r_delta_r, density).sum(1).sum(1)

        return np.divide(weightedSum, totalMass)


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


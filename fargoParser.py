__author__ = 'cguo'

import numpy as np
import glob
import re
import math
import logging

class FargoParser:
    """
    a parser for Fargo2D output files. Has the following properties:
    paramNames = ['numRadialIntervals', 'numThetaIntervals', 'radialIntervals', 'thetaIntervals',
                  'timeIntervals', 'maxRadius', 'totalNumOutputs']

    methods:
    FargoParser(outputDirectory, batchSize): creates parser, reads run parameters

    getParams(): returns a dict of params : param value, for each param in paramNames above

    getNextBatch(): returns a three-tuple of (density, vr, vtheta) for the next batch

    hasRemainingBatches(): returns True iff there are batches left
    """

    logging.basicConfig(level=logging.DEBUG, filename='parserDiagnostics.log', filemode='a')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    logging.getLogger('').addHandler(console)

    def __init__(self, outputDir, batchSize=100):
        if outputDir.endswith('/'):
            outputDir = outputDir[:-1]

        self.outputDir = outputDir
        self._readRunParams()

        self.batchSize = batchSize
        self.startIndex = 0


    def _readRunParams(self):
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

        filePaths = glob.glob(self.pathTo('gasdens*.dat'))
        self.totalNumOutputs = len(filePaths)


        paramNames = ['numRadialIntervals', 'numThetaIntervals', 'radialIntervals', 'thetaIntervals',
                      'timeIntervals', 'maxRadius', 'totalNumOutputs']

        self.params = dict((paramName, getattr(self, paramName)) for paramName in paramNames)


    def _extractFileIndex(self, filePath):
        name = filePath.split('/')[-1]
        return int(re.search('[0-9]+', name).group(0))


    def _pathTo(self, endpoint):
        return self.outputDir + "/" + endpoint


    def _parseGasValue(self, varType, startIndex, endIndex):
        """
        :param varType: "dens", "vrad", or "vtheta"
        :return:
        """

        logging.info("\n*** parsing values for gas" + varType + " ***\n")

        fileFormat = "gas" + varType + "*.dat"

        filePaths = glob.glob(self._pathTo(fileFormat))
        sortedPaths = sorted(filePaths, key=self._extractFileIndex)

        arrays = []
        for path in sortedPaths[startIndex : endIndex]:
            arr = np.fromfile(path, dtype='double')
            arr.shape = (self.numRadialIntervals, self.numThetaIntervals)
            arrays.append(arr)

            logging.info("parsed gas output file: " + path.split('/')[-1])

        ret = np.array(arrays)

        return ret


    def _parseGasOutput(self, startIndex, endIndex):
        """
        read and parse gas density, vrad, vtheta
        :return: tuple of (gasdens, gasvrad, gasvtheta)
        """

        varTypes = ['dens', 'vrad', 'vtheta']

        return (self._parseGasValue(varType, startIndex, endIndex) for varType in varTypes)


    def getParams(self):
        return self.params


    def hasRemainingBatches(self):
        return (self.totalNumOutputs - self.lastReadIndex) > 0


    def getNextBatch(self):
        # read files in [startIndex, endIndex)
        startIndex = self.startIndex
        endIndex = min(self.startIndex + self.batchSize, self.totalNumOutputs)
        self.startIndex = endIndex

        logging.info("\nreading batch from " + str(startIndex) + " to " + str(endIndex) + "\n")
        return self._parseGasOutput(startIndex, endIndex)

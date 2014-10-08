__author__ = 'cguo'

from fargoParser import FargoParser
from fargoPlotter import FargoPlotter
from optparse import OptionParser
import fargoDiagnostics as fd
import numpy as np


class FargoDiagnosticsRunner:

    def __init__(self, inputDir, outputDir, plotDir, batchSize):
        self.outputDir = outputDir

        self.parser = FargoParser(inputDir, batchSize)

        params = self.parser.getParams()
        radIntervals = params['radialIntervals'] * params['maxRadius']
        timeIntervals = params['timeIntervals']

        self.params = params

        self.outputDir = outputDir

        self.plotter = FargoPlotter(radIntervals, timeIntervals, plotDir, 'Radius, AU', 'Time, binary periods')

    def runBatches(self):
        i = 0
        while self.parser.hasRemainingBatches():
            dens, vrad, vtheta = self.parser.getNextBatch()
            calculations = fd.computeRadialDiagnostics(self.params['radialIntervals'],\
                                                             self.params['thetaIntervals'], dens, vrad, vtheta)
            for j in range(0, len(dens), 20):
                print 'plotting'
                print calculations['radialEccMK'][j]
                self.plotter.vsRadius(calculations['radialEccMK'][j], 'eccMK',\
                                      'Eccentricity (Mueller-Kley)', 'Eccentricity vs radius', i + j)

            np.save(self.outputDir + '/radialEccMK' + str(i), calculations['radialEccMK'])
            np.save(self.outputDir + '/radialPeri' + str(i), calculations['radialPeri'])
            np.save(self.outputDir + '/radialEccLubow' + str(i), calculations['radialEccLubow'])
            np.save(self.outputDir + '/radialDens' + str(i), calculations['radialDens'])
            np.save(self.outputDir + '/diskEccMK' + str(i), calculations['diskEccMK'])
            np.save(self.outputDir + '/diskPeri' + str(i), calculations['diskPeri'])

            i += len(dens)


def main():
    optParser = OptionParser()
    optParser.add_option('-i', '--inputdirectory', action='store',
                        type='string', dest='inputDirectory')

    optParser.add_option('-b', '--batchsize', action='store',
                         type='int', dest='batchSize', default=100)
    
    optParser.add_option('-o', '--outputdirectory', action='store',
                         type='string', dest='outputDirectory')

    optParser.add_option('-p', '--plotdirectory', action='store',
                         type='string', dest='plotDirectory')

    (options, args) = optParser.parse_args()

    if not options.inputDirectory:
        optParser.error('you must specify an input directory with -i or --inputdirectory')

    runner = FargoDiagnosticsRunner(options.inputDirectory, options.outputDirectory, options.plotDirectory, 100)
    runner.runBatches()

if __name__ == '__main__':
    main()
__author__ = 'cguo'

from fargoParser import FargoParser
from fargoPlotter import FargoPlotter
from optparse import OptionParser
import fargoDiagnostics as fd


class FargoDiagnosticsRunner:

    def __init__(self, inputDir, outputDir, plotDir, batchSize):
        self.outputDir = outputDir

        self.parser = FargoParser(inputDir, batchSize)

        params = self.parser.getParams()
        radIntervals = params['radialIntervals'] * params['maxRadius']
        timeIntervals = params['timeIntervals']

        self.params = params

        self.plotter = FargoPlotter(radIntervals, timeIntervals, plotDir, 'Radius, AU', 'Time, binary periods')

    def runBatches(self):
        while self.parser.hasRemainingBatches():
            dens, vrad, vtheta = self.parser.getNextBatch()
            radialMeasurements = fd.computeRadialDiagnostics(self.params['radialIntervals'],\
                                                             self.params['thetaIntervals'], dens, vrad, vtheta)
            


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


if __name__ == '__main__':
    main()
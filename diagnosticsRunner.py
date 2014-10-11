__author__ = 'cguo'

from fargoParser import FargoParser
from fargoPlotter import FargoPlotter
from optparse import OptionParser
import fargoDiagnostics as fd
import numpy as np
import glob


class FargoDiagnosticsRunner:

    def __init__(self, inputDir, outputDir, plotDir, batchSize):
        self.outputDir = outputDir

        self.parser = FargoParser(inputDir, batchSize)

        params = self.parser.getParams()
        radIntervals = params['radialIntervals'] * 20
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
                                      'Eccentricity (Mueller-Kley)', 'Eccentricity vs radius', i + j, [0, 0.35])

                self.plotter.vsRadius(calculations['radialEccLubow'][j], 'eccLubow',\
                                      'Eccentricity (Lubow)', 'Lubow eccentricity vs radius', i + j, [0, 0.35])

                self.plotter.vsRadius(calculations['radialPeri'][j], 'peri',\
                                      'Periastron angle', 'Periastron angle vs radius', i + j, [-3.1, 3.1])

            np.save(self.outputDir + '/radialEccMK' + str(i), calculations['radialEccMK'])
            np.save(self.outputDir + '/radialPeri' + str(i), calculations['radialPeri'])
            np.save(self.outputDir + '/radialEccLubow' + str(i), calculations['radialEccLubow'])
            np.save(self.outputDir + '/radialDens' + str(i), calculations['radialDens'])
            np.save(self.outputDir + '/diskEccMK' + str(i), calculations['diskEccMK'])
            np.save(self.outputDir + '/diskPeri' + str(i), calculations['diskPeri'])

            i += len(dens)


    def runDiskTime(self):
        diagnosticTypes = [
            {
                'fileFormat': '/diskEccMK*.npy',
                'arrayFilename': 'eccMKVsTime.npy',
                'yName': 'diskEccMK',
                'yLabel': 'Disk eccentricity (Mueller-Kley)',
                'title': 'Disk eccentricity vs time'
            },
            {
                'fileFormat': '/diskPeri*.npy',
                'arrayFilename': 'periVsTime.npy',
                'yName': 'diskPeri',
                'yLabel': 'Disk periastron angle',
                'title': 'Disk periastron vs time'
            }
        ]

        for type in diagnosticTypes:
            filePaths = glob.glob(self.outputDir + type['fileFormat'])
            sortedPaths = sorted(filePaths, key=self.parser._extractFileIndex)

            arrays = []
            for path in sortedPaths:
                arrays.append(np.load(path))
            vsTime = np.concatenate(arrays)
            arrays = None

            np.save(self.outputDir + '/' + type['arrayFilename'], vsTime)

            self.plotter.vsTime(vsTime, type['yName'], type['yLabel'], type['title'])


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

    optParser.add_option('-d', '--diskonly', action='store_true',
                         dest='diskOnly')

    (options, args) = optParser.parse_args()

    if not options.inputDirectory:
        optParser.error('you must specify an input directory with -i or --inputdirectory')

    runner = FargoDiagnosticsRunner(options.inputDirectory, options.outputDirectory, options.plotDirectory, 100)
    if not options.diskOnly:
        runner.runBatches()
    runner.runDiskTime()

if __name__ == '__main__':
    main()
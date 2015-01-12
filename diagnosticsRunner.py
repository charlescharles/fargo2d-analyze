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
        radIntervals = params['radialIntervals']
        numOutputs = params['totalNumOutputs']
        timeIntervals = np.linspace(0, numOutputs/5.0, num=numOutputs)

        self.params = params

        self.outputDir = outputDir

        self.plotter = FargoPlotter(radIntervals * 20.0, timeIntervals, plotDir, 'Radius, AU', 'Time, binary periods')

    def runBatches(self):
        i = 0
        while self.parser.hasRemainingBatches():
            dens, vrad, vtheta = self.parser.getNextBatch()
            calculations = fd.computeDiagnostics(self.params['radialIntervals'],
                                                             self.params['thetaIntervals'], dens, vrad, vtheta)
            for j in range(0, len(dens), 20):
                print 'plotting'
                # self.plotter.vsRadius(calculations['radialEccMK'][j], 'eccMK',
                #                       'Eccentricity (Mueller-Kley)', 'Eccentricity vs radius', i + j, [0, 0.35])
                #
                # self.plotter.vsRadius(calculations['radialEccLubow'][j], 'eccLubow',
                #                       'Eccentricity (Lubow)', 'Lubow eccentricity vs radius', i + j, [0, 0.35])
                #
                # self.plotter.vsRadius(calculations['radialPeriMK'][j], 'periMK',
                #                       'Periastron (MK)', 'Periastron vs radius', i + j, [-3.1, 3.1])
                #
                # self.plotter.vsRadius(calculations['radialPeriLubow'][j], 'periLubow',
                #                       'Periastron (Lubow)', 'Periastron vs radius', i + j, [-3.1, 3.1])

                print "length of radialDens: " + str(len(calculations['radialDens']))

                self.plotter.threePanelVsRadius(calculations['radialDens'][j],
                                                calculations['radialEccMK'][j], calculations['radialEccLubow'][j],
                                                calculations['radialPeriMK'][j], calculations['radialPeriLubow'][j],
                                                "%.1f" % ((i + j)/5.0), 'threePanel', i + j)

            np.save(self.outputDir + '/radialEccMK' + str(i), calculations['radialEccMK'])
            np.save(self.outputDir + '/radialEccLubow' + str(i), calculations['radialEccLubow'])
            np.save(self.outputDir + '/radialPeriMK' + str(i), calculations['radialPeriMK'])
            np.save(self.outputDir + '/radialPeriLubow' + str(i), calculations['radialPeriLubow'])
            np.save(self.outputDir + '/radialDens' + str(i), calculations['radialDens'])
            np.save(self.outputDir + '/diskEccMK' + str(i), calculations['diskEccMK'])
            np.save(self.outputDir + '/diskPeriMK' + str(i), calculations['diskPeriMK'])
            np.save(self.outputDir + '/diskEccLubow' + str(i), calculations['diskEccLubow'])
            np.save(self.outputDir + '/diskPeriLubow' + str(i), calculations['diskPeriLubow'])
            np.save(self.outputDir + '/totalMass' + str(i), calculations['totalMass'])

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
                'fileFormat': '/diskPeriMK*.npy',
                'arrayFilename': 'periMKVsTime.npy',
                'yName': 'diskPeriMK',
                'yLabel': 'Disk periastron angle (MK)',
                'title': 'Disk periastron vs time'
            },
            {
                'fileFormat': '/diskEccLubow*.npy',
                'arrayFilename': 'eccLubowVsTime.npy',
                'yName': 'diskEccLubow',
                'yLabel': 'Disk eccentricity (Lubow)',
                'title': 'Disk eccentricity vs time'
            },
            {
                'fileFormat': '/diskPeriLubow*.npy',
                'arrayFilename': 'periLubowVsTime.npy',
                'yName': 'diskPeriLubow',
                'yLabel': 'Disk periastron angle (Lubow)',
                'title': 'Disk periastron vs time'
            },
            {
                'fileFormat': '/totalMass*.npy',
                'arrayFilename': 'massVsTime.npy',
                'yName': 'totalMass',
                'yLabel': 'Disk mass (code units)',
                'title': 'Disk mass vs time'
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

    runner = FargoDiagnosticsRunner(options.inputDirectory, options.outputDirectory, options.plotDirectory, options.batchSize)
    if not options.diskOnly:
        runner.runBatches()
    runner.runDiskTime()

if __name__ == '__main__':
    main()
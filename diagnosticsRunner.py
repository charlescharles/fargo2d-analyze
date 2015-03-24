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

    def _getDiagnostic(self, fmt):
        filePaths = glob.glob(self.outputDir + fmt)
        sortedPaths = sorted(filePaths, key=self.parser._extractFileIndex)

        arrays = [np.load(path) for path in sortedPaths]
        if len(arrays) == 0:
            print "didn't find any arrays for format " + fmt
            return []

        return np.concatenate(arrays)

    def runBatches(self):
        i = 0
        while self.parser.hasRemainingBatches():
            dens, vrad, vtheta = self.parser.getNextBatch()
            calculations = fd.computeDiagnostics(self.params['radialEdges'], self.params['radialIntervals'],
                                                             self.params['thetaIntervals'], dens, vrad, vtheta)

            avgDens = np.average(dens, axis=2)

            for j in range(0, len(dens), 20):
                print 'plotting'
                print "length of radialDens: " + str(len(calculations['radialDens']))

                self.plotter.threePanelVsRadius(avgDens[j],
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

            np.save(self.outputDir + '/diskRadius90' + str(i), calculations['diskRad90'])
            np.save(self.outputDir + '/diskRadius95' + str(i), calculations['diskRad95'])

            np.save(self.outputDir + '/lubowVsin' + str(i), calculations['lubowVsin'])
            np.save(self.outputDir + '/lubowVcos' + str(i), calculations['lubowVcos'])

            i += len(dens)


    def runDiskTime(self):
        diagnosticTypes = [
            {
                'fileFormat': '/diskEccMK*.npy',
                'arrayFilename': 'eccMKVsTime.npy',
                'yName': 'diskEccMK',
                'yLabel': 'Disk eccentricity (Mueller-Kley)',
                'title': 'Disk eccentricity vs time',
                'plot': True
            },
            {
                'fileFormat': '/diskPeriMK*.npy',
                'arrayFilename': 'periMKVsTime.npy',
                'yName': 'diskPeriMK',
                'yLabel': 'Disk periastron angle (MK)',
                'title': 'Disk periastron vs time',
                'plot': True
            },
            {
                'fileFormat': '/diskEccLubow*.npy',
                'arrayFilename': 'eccLubowVsTime.npy',
                'yName': 'diskEccLubow',
                'yLabel': 'Disk eccentricity (Lubow)',
                'title': 'Disk eccentricity vs time',
                'plot': True
            },
            {
                'fileFormat': '/diskPeriLubow*.npy',
                'arrayFilename': 'periLubowVsTime.npy',
                'yName': 'diskPeriLubow',
                'yLabel': 'Disk periastron angle (Lubow)',
                'title': 'Disk periastron vs time',
                'plot': True
            },
            {
                'fileFormat': '/totalMass*.npy',
                'arrayFilename': 'massVsTime.npy',
                'yName': 'totalMass',
                'yLabel': 'Disk mass (code units)',
                'title': 'Disk mass vs time',
                'plot': True
            },
            {
                'fileFormat': '/diskRadius90*.npy',
                'arrayFilename': 'radius90VsTime.npy',
                'yName': 'diskRadius90',
                'yLabel': 'Disk radius (a, 90%)',
                'title': 'Disk radius vs time',
                'plot': True
            },
            {
                'fileFormat': '/diskRadius95*.npy',
                'arrayFilename': 'radius95VsTime.npy',
                'yName': 'diskRadius95',
                'yLabel': 'Disk radius (a, 95%)',
                'title': 'Disk radius vs time',
                'plot': True
            },
            {
                'fileFormat': '/lubowVsin*.npy',
                'arrayFilename': 'vsinVsTime.npy',
                'yName': 'lubowVsin',
                'yLabel': 'Lubow V_sin',
                'title': 'Lubow V_sin',
                'plot': False
            },
            {
                'fileFormat': '/lubowVcos*.npy',
                'arrayFilename': 'vcosVsTime.npy',
                'yName': 'lubowVcos',
                'yLabel': 'Lubow V_cos',
                'title': 'Lubow V_cos',
                'plot': False
            }
        ]

        diags = {}

        for type in diagnosticTypes:
            diag = self._getDiagnostic(type['fileFormat'])
            diags[type['yName']] = diag

            np.save(self.outputDir + '/' + type['arrayFilename'], diag)

        for type in diagnosticTypes:
            diag = diags[type['yName']]
            if type['plot']:
                print "plotting vs time " + type['yName']
                self.plotter.vsTime(diag, type['yName'], type['yLabel'], type['title'])

        eccMK = diags['diskEccMK']
        periMK = diags['diskPeriMK']
        print "plotting twopanel vs time"
        self.plotter.twoPanelVsTime(eccMK, periMK, "eccPeriMK_vs_time")


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
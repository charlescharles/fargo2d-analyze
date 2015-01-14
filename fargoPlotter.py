__author__ = 'cguo'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

class FargoPlotter:
    """
    Produces and saves plots against radius and time

    methods:

    vsRadius(array, yName='', yDisplayLabel='', title='')

    vsTime(array, yName='', yDisplayLabel='', title='')
    """
    def __init__(self, radialIntervals, timeIntervals, outputDir, radialLabel='', timeLabel=''):
        self.radialIntervals = radialIntervals
        self.timeIntervals = timeIntervals
        self.outputDir = outputDir[:-1] if outputDir.endswith('/') else outputDir
        self.radialLabel = radialLabel
        self.timeLabel = timeLabel

        plt.ioff()

    def _setFigsize(self, size):
        plt.rcParams['figure.figsize'] = size

    def _resetFigsize(self):
        self._setFigsize((8, 6))

    def _pathTo(self, fileName):
        if not fileName.startswith('/'):
            return self.outputDir + '/' + fileName
        return self.outputDir + fileName


    def threePanelVsRadius(self, density, eccMK, eccLubow, periMK, periLubow, time, fname, index):
        self._setFigsize((7, 11))

        logDensity = np.log10(density)
        fig = plt.figure()

        fig.suptitle(time + " binary periods")

        plt.subplot(3, 1, 1)
        plt.plot(self.radialIntervals, logDensity)
        plt.xlabel(self.radialLabel)
        plt.ylabel('log(density)')

        plt.subplot(3, 1, 2)
        plt.plot(self.radialIntervals, eccMK, 'k', label="Mueller-Kley")
        plt.plot(self.radialIntervals, eccLubow, 'k--', label="Lubow")
        plt.xlabel(self.radialLabel)
        plt.ylabel('Eccentricity')
        plt.legend(loc='upper right')

        plt.subplot(3, 1, 3)
        plt.plot(self.radialIntervals, periMK, 'k', label="Mueller-Kley")
        plt.plot(self.radialIntervals, periLubow, 'k--', label="Lubow")
        plt.xlabel(self.radialLabel)
        plt.ylabel('Periastron')
        plt.legend(loc='upper right')

        fname = self._pathTo(fname + str(index) + '.png')

        print "saving 3-panel figure with name " + fname
        plt.savefig(fname)
        plt.close(fig)

        self._resetFigsize()

    def twoPanelVsTime(self, ecc, peri, fname):
        self._setFigsize((8, 12))

        fig = plt.figure()

        usedTimes = self.timeIntervals[:len(ecc)]

        plt.subplot(2, 1, 1)
        plt.plot(usedTimes, ecc)
        plt.xlabel(self.timeLabel)
        self.ylabel("Eccentricity")

        plt.subplot(2, 1, 2)
        plt.plot(usedTimes, peri)
        plt.xlabel(self.timeLabel)
        self.ylabel("Periastron")

        fname = self._pathTo(fname + '.png')

        print "saving 2-panel figure with name " + fname
        plt.savefig(fname + '.png')
        plt.close(fig)

        self._resetFigsize()


    def vsRadius(self, array, yName='', yDisplayLabel='', title='', index='', ylim=None):
        fig = plt.figure()
        plt.plot(self.radialIntervals, array)
        plt.xlabel(self.radialLabel)
        plt.ylabel(yDisplayLabel)
        plt.title(title)
        if ylim:
            plt.ylim(ylim)

        fname = self._pathTo(yName + '_vs_radius' + str(index) + '.png')

        print "saving figure with name " + fname
        plt.savefig(fname)
        plt.close(fig)


    def vsTime(self, array, yName='', yDisplayLabel='', title=''):
        plt.rcParams['figure.figsize'] = 9, 6

        fig = plt.figure()
        print 'array len is ' + str(len(array))
        print 'timeinterval len is ' + str(len(self.timeIntervals))
        print 'yname is ' + yName
        plt.plot(self.timeIntervals[:len(array)], array)
        plt.xlabel(self.timeLabel)
        plt.ylabel(yDisplayLabel)
        plt.title(title)

        plt.savefig(self._pathTo(yName + '_vs_time.png'))
        plt.close(fig)


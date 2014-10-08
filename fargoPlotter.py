__author__ = 'cguo'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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


    def _pathTo(self, fileName):
        if not fileName.startswith('/'):
            return self.outputDir + '/' + fileName
        return self.outputDir + fileName


    def vsRadius(self, array, yName='', yDisplayLabel='', title='', index=''):
        print len(array)
        print len(self.radialIntervals)

        fig = plt.figure()
        plt.plot(self.radialIntervals, array)
        plt.xlabel(self.radialLabel)
        plt.ylabel(yDisplayLabel)
        plt.title(title)

        fname = self._pathTo(yName + '_vs_radius' + str(index) + '.png')

        print "saving figure with name " + fname
        plt.savefig(fname)
        plt.close(fig)


    def vsTime(self, array, yName='', yDisplayLabel='', title=''):
        fig = plt.figure()
        plt.plot(self.timeIntervals, array)
        plt.xlabel(self.timeLabel)
        plt.ylabel(yDisplayLabel)
        plt.title(title)

        plt.savefig(self._pathTo(yName + '_vs_time.png'))
        plt.close(fig)


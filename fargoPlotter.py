__author__ = 'cguo'

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
        return self.outputDir + fileName


    def vsRadius(self, array, yName='', yDisplayLabel='', title=''):
        fig = plt.figure()
        plt.plot(self.radialIntervals, array)
        plt.xlabel(self.radialLabel)
        plt.ylabel(yDisplayLabel)
        plt.title(title)

        plt.savefig(self._pathTo(yName + '_vs_radius.png'))
        plt.close(fig)


    def vsTime(self, array, yName='', yDisplayLabel='', title=''):
        fig = plt.figure()
        plt.plot(self.timeIntervals, array)
        plt.xlabel(self.timeLabel)
        plt.ylabel(yDisplayLabel)
        plt.title(title)

        plt.savefig(self._pathTo(yName + '_vs_time.png'))
        plt.close(fig)


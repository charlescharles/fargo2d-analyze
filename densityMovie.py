__author__ = 'cguo'

from fargoParser import FargoParser
from optparse import OptionParser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

class FargoMovieMaker:

    def __init__(self, inputDir, outputDir, batchSize):
        self.outputDir = outputDir
        self.batchSize = batchSize
        self.parser = FargoParser(inputDir, batchSize)

        secondaryOrbit = np.loadtxt(inputDir + "/planet0.dat")
        secondaryX = secondaryOrbit[:, 1]
        secondaryY = secondaryOrbit[:, 2]

        self.secondaryRadius = np.sqrt(np.add(np.square(secondaryX), np.square(secondaryY)))
        self.secondaryTheta = np.arctan2(secondaryY, secondaryX)

        params = self.parser.getParams()
        self.params = params
        self.outputDir = outputDir

    def go(self, start, end):

        start = (start / self.batchSize) * self.batchSize
        end = (end / self.batchSize + 1) * self.batchSize
        cur = 0

        # skip to the start
        for _ in range(start / self.batchSize):
            cur += self.batchSize
            self.parser.getNextBatch()

        for _ in range((end - start) / self.batchSize):
            dens, _, _ = self.parser.getNextBatch()

            r, theta = np.meshgrid(self.params['radialIntervals'], self.params['thetaIntervals'])
            plt.ioff()
            #-- Plot... ------------------------------------------------
            for i in range(len(dens)):
                fig = plt.figure()
                ax = plt.subplot(111, polar=True)
                ax.contourf(theta, r, np.log(dens[i]).transpose(), cmap=plt.cm.afmhot)
                ax.scatter([self.secondaryTheta[cur]], [self.secondaryRadius[cur]], s=150)
                ax.set_rmax(1.5)
                #ax.set_title(r"$\theta_{sec}=" + "{0:.2f}$ rad".format(self.secondaryTheta[cur] % 6.283), va='bottom')
                plt.savefig(self.outputDir + "/figs/dens" + str(cur) + ".png")
                plt.close(fig)

                cur += 1

    def finish(self):
        os.system("tar -zcvf " + self.outputDir + "/animation.tar.gz " + self.outputDir + "/figs")


def main():
    optParser = OptionParser()
    optParser.add_option('-i', '--inputdirectory', action='store',
                        type='string', dest='inputDirectory')

    optParser.add_option('-b', '--batchsize', action='store',
                         type='int', dest='batchSize', default=100)

    optParser.add_option('-s', '--start', action='store',
                         type='int', dest='startIndex', default=0)

    optParser.add_option('-e', '--end', action='store',
                         type='int', dest='endIndex', default=100)
    
    optParser.add_option('-o', '--outputdirectory', action='store',
                         type='string', dest='outputDirectory')

    (options, args) = optParser.parse_args()

    if not options.inputDirectory:
        optParser.error('you must specify an input directory with -i or --inputdirectory')

    movies = FargoMovieMaker(options.inputDirectory, options.outputDirectory, options.batchSize)
    movies.go(options.startIndex, options.endIndex)
    movies.finish()

if __name__ == '__main__':
    main()
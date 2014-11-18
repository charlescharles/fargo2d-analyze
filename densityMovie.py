__author__ = 'cguo'

from fargoParser import FargoParser
from optparse import OptionParser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class FargoMovieMaker:

    def __init__(self, inputDir, outputDir, batchSize):
        self.outputDir = outputDir

        self.parser = FargoParser(inputDir, batchSize)

        secondaryOrbit = np.loadtxt(inputDir + "/planet0.dat")
        secondaryX = secondaryOrbit[:, 1]
        secondaryY = secondaryOrbit[:, 2]

        self.secondaryRadius = np.sqrt(np.add(np.square(secondaryX), np.square(secondaryY)))
        self.secondaryTheta = np.arctan2(secondaryY, secondaryX)

        params = self.parser.getParams()
        self.params = params
        self.outputDir = outputDir

    def go(self):
        dens, _, _ = self.parser.getNextBatch()

        r, theta = np.meshgrid(self.params['radialIntervals'], self.params['thetaIntervals'])
        plt.ioff()
        #-- Plot... ------------------------------------------------
        for i in range(100):
            fig = plt.figure()
            ax = plt.subplot(111, polar=True)
            ax.contourf(theta, r, np.log(dens[i]).transpose())
            ax.scatter([self.secondaryTheta[i]], [self.secondaryRadius[i]], s=150)
            ax.set_rmax(1.5)
            ax.set_title(r'$\theta_{sec}=' + "{0:.2f}".format(self.secondaryTheta[i] % 6.283) + r'$ rad', va='bottom')
            plt.savefig(self.outputDir + "/figs/dens" + str(i) + ".png")
            plt.close(fig)



def main():
    optParser = OptionParser()
    optParser.add_option('-i', '--inputdirectory', action='store',
                        type='string', dest='inputDirectory')

    optParser.add_option('-b', '--batchsize', action='store',
                         type='int', dest='batchSize', default=100)
    
    optParser.add_option('-o', '--outputdirectory', action='store',
                         type='string', dest='outputDirectory')

    (options, args) = optParser.parse_args()

    if not options.inputDirectory:
        optParser.error('you must specify an input directory with -i or --inputdirectory')

    movies = FargoMovieMaker(options.inputDirectory, options.outputDirectory, options.batchSize)
    movies.go()

if __name__ == '__main__':
    main()
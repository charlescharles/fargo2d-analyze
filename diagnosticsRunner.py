__author__ = 'cguo'

from fargoParser import FargoParser
from optparse import OptionParser
import os

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

    fargoParser = FargoParser(options.inputDirectory, options.batchSize)

    # outputDir = os.path.dirname(options.inputDirectory + '/parsedOutput')
    #
    # if not os.path.exists(outputDir):
    #     os.mkdir(outputDir)
    #
    # outputProperties = ['dims', 'gasVariables', 'cellEccentricities']
    #
    # for prop in outputProperties:
    #     fargoParser.saveProperty(prop)


if __name__ == '__main__':
    main()
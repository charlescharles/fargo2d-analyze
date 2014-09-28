__author__ = 'cguo'

from fargoParser import FargoParser
from optparse import OptionParser

def main():
    optParser = OptionParser()
    optParser.add_option("-i", "--inputdirectory", action="store",
                      type="string", dest="inputDirectory")

    (options, args) = optParser.parse_args()

    if not options.inputDirectory:
        optParser.error("you must specify an input directory with -i or --inputdirectory")

    fargoParser = FargoParser(options.inputDirectory)
    fargoParser.computeCellEccentricities()

    fargoParser.save(fargoParser.cellEccentricities, "cellEccentricity")


if __name__ == "__main__":
    main()
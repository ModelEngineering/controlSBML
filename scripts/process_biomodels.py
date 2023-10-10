""" Iterates across all models in BioModels to do staircase and make a SISOClosedLoop"""

from controlSBML.siso_maker import SISOMaker

import argparse

def main(first, last):
    SISOMaker.runBiomodels(start=first, end=last)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process BioModels")
    parser.add_argument("--first", help="first model tp process (0)", default=0, type=int)
    parser.add_argument("--last", help="last model tp process (5000)", default=5000, type=int)
    args = parser.parse_args()
    main(args.first, args.last)
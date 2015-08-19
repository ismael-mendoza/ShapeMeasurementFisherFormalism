#!/usr/bin/env python
"""Allows the user to do N fittings of a galaxy produced in generate.py given a SNR and compare with the fisher formalism by displaying biases and 
correlation coefficients,etc. and also by displaying a triangle plot that summarizes this results and displays data points."""

import argparse

import os

def main():
    """Give the instructions to run a bunch of fits either here or in SLAC and to reproduce results and produce triangle plot."""

    #have to add a require argument: filename to write and extract galaxies from.

    parser = argparse.ArgumentParser(description = 'Displays a triangle plot comparing the expected noise bias calculation with the fisher matrix formalism for a given file of galaxies.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--wdir', default = 'output', metavar = 'DIRECTORY', type = str,
    help = 'Specify a directory name where the output will be inputted and files read from.')
    parser.add_argument('--galaxy-file', default = 'galaxies', metavar = 'FILENAME', type = str,
    help = 'Specify a file where galaxies will be registered or read from.')
    parser.add_argument('-rltsdir', default = 'results', metavar = 'RESULTS', type = str,
    help = 'Specify a directory where output data from fits will be produced')
    parser.add_argument('--snr', metavar = 'SNR', default = 30, type = int,
    help = 'Specify signal to noise ratio of the run. Should use same snr between run_fits and produce_results, else it does not work.')
    parser.add_argument('--number-fits', metavar = 'N', default = 1, type = int,
    help = 'Specify how many fits are run.')
    parser.add_argument('--run-fits', action = 'store_true',
    help = 'Run N different instantiations of noise of the added galaxies to produce a triangle plot. Fits all galaxies in given file N times.')
    parser.add_argument('--run-fits-slac', metavar = 'SLAC_COMPUTER',
    help = 'Same as above but have to be logged in a SLAC computer.')
    parser.add_argument('--produce-results', action = 'store_true',
    help = 'Read each of the files produces from a fits trial in rltsdir and create a .csv file in the same directory, triangle.png (the graph of results) and sumresults.csv (condensation ofr all files that are processed by triangle.png)')


    args = parser.parse_args()

    if args.run_fits:
        for i in range(args.number_fits):
            os.system("python runfits.py " + str(i) + " " +  str(args.snr) + " " +  str(args.wdir) + " " + str(args.galaxy_file) + " " + str(args.rltsdir))

    elif args.run_fits_slac:
        os.system("bsub -q " + args.run_fits_slac + " -J \"name[1-" + str(args.number_fits) + "]\" \"python runfits.py \$LSB_JOBINDEX " 
            + args.snr + " " + args.wdir + " " + args.galaxy_file + " " + args.rltsdir + "\"")

    elif args.produce_results:
        os.system("python readfits.py " +  str(args.wdir) + " " + str(args.galaxy_file) + " " + str(args.rltsdir) + " " + str(args.snr))

if __name__ == '__main__':
    main()
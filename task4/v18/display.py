#!/usr/bin/env python
"""This program allows the user to display various statistical results (like plots) for the galaxies generated in generate.py"""

import argparse

import fisher

import plotsfisher 

import defaults

import os

import galfun

def main():
    names = defaults.names()
    cts = defaults.constants()

    parser = argparse.ArgumentParser(description = 'Display different results and plots from the galaxies in a given file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v','--verbose', action = 'store_true',
    help = 'Prints parameters of the galaxy used in the fisher formalism, plus default parameters important for the analysis.')
    parser.add_argument('--info-file', action = 'store_true',
    help = 'Write relevant information to a .txt file for future reference')
    parser.add_argument('--wdir', default = names.wdir, metavar = 'DIRECTORY', type = str,
    help = 'Specify a directory name where the output will be inputted and files read from.')
    parser.add_argument('--galaxy-file', default = names.galaxy_file, metavar = 'FILENAME', type = str,
    help = 'Specify a file where galaxies will be registered or read from.')
    parser.add_argument('--plots-dir', default = names.plots_dir, metavar = 'PLOTDIRECTORY', type = str,
    help = 'Specify a directory name where the plots will be saved')
    parser.add_argument('--hide', action = 'store_true',
    help = 'Do not show plots produced')
    parser.add_argument('-d', '--draw-galaxy', action = 'store_true',
    help = 'Show original galaxies')
    parser.add_argument('--partials', action = 'store_true',
    help = 'Show partial derivative images.')
    parser.add_argument('--second-partials', action = 'store_true',
    help = 'Show second partial derivative images.')
    parser.add_argument('-f', '--fisher', action = 'store_true',
    help = 'Show fisher matrix images.')
    parser.add_argument('--fisher-chi2', action = 'store_true',
    help = 'Show fisher matrix images with chi2 calculated values on top.')
    parser.add_argument('--covariance', action = 'store_true',
    help = 'Show covariance matrix elements')
    parser.add_argument('--correlation', action = 'store_true',
    help = 'Show correlation matrix elements.')   
    parser.add_argument('--bias-matrix', action = 'store_true',
    help = 'Show bias matrix images.')
    parser.add_argument('-b', '--biases', action = 'store_true',
    help = 'Show bias images.')
    parser.add_argument('--values', action = 'store_true', 
    help = 'Show values on top of appropiate images (fisher, biases).')
    parser.add_argument('--sigma-n', metavar = 'SIGMA_N', default = cts.sigma_n,
    help = 'Value of noise bias (standard deviation) in each pixel.')


    args = parser.parse_args()

    galaxies = galfun.galaxies(args.wdir, args.galaxy_file)

    fish = fisher.fisher_analysis(galaxies = galaxies, sigma_n = args.sigma_n)
    fisherplots = plotsfisher.fisherplots(fisher = fish, wdir = args.wdir, plots_dir = args.plots_dir, hide = args.hide)

    if args.draw_galaxy: 
        fisherplots.galaxy()

    if args.partials: 
        fisherplots.derivatives()

    if args.fisher and args.values: 
        fisherplots.fisherMatrixValues()

    elif args.fisher: 
        fisherplots.fisherMatrix()

    if args.fisher_chi2: 
        fisherplots.fisherMatrixChi2()

    if args.covariance: 
        fisherplots.covarianceMatrix()

    if args.correlation: 
        fisherplots.correlationMatrix()

    if args.second_partials: 
        fisherplots.secondDerivatives()
 
    if args.bias_matrix and args.values: 
        fisherplots.biasMatrixValues()

    elif args.bias_matrix: 
        fisherplots.biasMatrix()

    if args.biases and args.values: 
        fisherplots.biasValues()

    elif args.biases: 
        fisherplots.biases()

    
    info = defaults.info(galaxies, fish)
    if args.verbose: 
        for line in info.galaxy + info.fisher:
            print line

    if args.info_file:
        with open(os.path.join(args.wdir, names.info + '.txt'), 'w') as txtfile:
            for line in info.galaxy + info.fisher:
                txtfile.write(line + '\n') #last character to skip lines.

if __name__ == '__main__':
    main()
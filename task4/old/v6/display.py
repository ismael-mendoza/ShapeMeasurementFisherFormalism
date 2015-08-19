import argparse
import defaults
import os
import matplotlib.pyplot as plt
from functions import *
import csv
import math
import numpy as np
from fisher import *
from plotfisher import *

def main():
    parser = argparse.ArgumentParser(description = 'Display different results and plots from the galaxies in a given file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--wdir', default = 'output', metavar = 'DIRECTORY', type = str,
    help = 'Specify a directory name where the output will be inputted and files read from.')
    parser.add_argument('--galaxy-file', default = 'galaxies', metavar = 'FILENAME', type = str,
    help = 'Specify a file where galaxies will be registered or read from.')
    parser.add_argument('--plots-dir', default = 'plots', metavar = 'PLOTDIRECTORY', type = str,
    help = 'Specify a directory name where the plots will be saved')
    parser.add_argument('--hide', action = 'store_true',
    help = 'Do not show plots produced')
    parser.add_argument('--draw-galaxy', action = 'store_true',
    help = 'Show original galaxies')
    parser.add_argument('--partials', action = 'store_true',
    help = 'Show partial derivative images.')
    parser.add_argument('--second-partials', action = 'store_true',
    help = 'Show second partial derivative images.')
    parser.add_argument('--fisher', action = 'store_true',
    help = 'Show fisher matrix images.')
    parser.add_argument('--fisher-chi2', action = 'store_true',
    help = 'Show fisher matrix images with chi2 calculated values on top.')
    parser.add_argument('--covariance', action = 'store_true',
    help = 'Show covariance matrix elements')
    parser.add_argument('--correlation', action = 'store_true',
    help = 'Show correlation matrix elements.')   
    parser.add_argument('--bias-matrix', action = 'store_true',
    help = 'Show bias matrix images.')
    parser.add_argument('--biases', action = 'store_true',
    help = 'Show bias images.')
    parser.add_argument('--values', action = 'store_true', 
    help = 'Show values on top of appropiate images (fisher, biases).')
    parser.add_argument('--sigma-n', metavar = 'SIGMA_N', default = 1.,
    help = 'Value of noise bias (standard deviation) in each pixel.')

    ##maybe change steps from default. 


    args = parser.parse_args()


    #some checks.
    if not os.path.isdir(args.wdir):
        print ('Directory does not exists')
        return -1

    filename = os.path.join(args.wdir, args.galaxy_file + '.csv')
    if not os.path.isfile(filename):
        print('Galaxies file does not exist')
        return -1

    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        lst_params = [] #possible here create dictionary with total number of params of both galaxies.
        for row in reader: 
            lst_params.append(row)

    params = lst_params[0] 

    #convert number values in params to floats.
    for key, value in params.iteritems():
        try:
            params[key] = float(value)
        except ValueError:
            pass

    #possible here create all possible galaxies with all params. 
    gal_image = drawGalaxy(params)

    fisher_analysis = fisher.fisher(params = params, gal_image = gal_image, sigma_n = args.sigma_n)
    fisherplot = fisherplot.plotfisher(fisher = fisher_analysis, wdir = args.wdir, plots_dir = args.plots_dir, hide = args.hide)



    if args.draw_galaxy: 
        fisherplot.galaxy()

    if args.partials: 
        fisherplot.derivatives()

    if args.fisher and args.values: 
        fisherplot.fisherMatrixValues()

    elif args.fisher: 
        fisherplot.fisherMatrix()

    if args.fisher_chi2: 
        fisherplot.fisherMatrixChi2()

    if args.covariance: 
        fisherplot.covarianceMatrix()

    if args.correlation: 
        fisherplot.correlationMatrix()

    if args.second_partials: 
        fisherplot.secondDerivatives()
 
    if args.bias_matrix and args.values: 
        fisherplot.biasMatrixValues()

    elif args.bias_matrix: 
        fisherplot.biasMatrix()

    if args.biases and args.values: 
        fisherplot.biasValues()

    elif args.biases: 
        fisherplot.biases()




if __name__ == '__main__':
    main()
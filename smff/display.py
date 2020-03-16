#!/usr/bin/env python
"""This program allows the user to display various statistical results (like
plots) for the galaxies generated in generate.py
"""

import argparse

import analysis.fisher as fisher

import analysis.draw as draw

import analysis.defaults as defaults

import analysis.galfun as galfun

#import pickle 

def main():
    parser = argparse.ArgumentParser(description=(
        'Display different results fisher results for a galaxy.'), formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--project', default=defaults.PROJECT,
                        type=str,
                        help=('Specify a directory name where the project will'
                              'be saved. In this fashion each individual'
                              'project represents one analysis.'))

    parser.add_argument('--hide', action='store_true',
                        help='Do not show plots produced')

    parser.add_argument('-d', '--draw-galaxy', action='store_true',
                        help='Show original galaxies')

    parser.add_argument('-fp', '--partials', action='store_true',
                        help='Show partial derivative images.')

    parser.add_argument('-sp', '--second-partials', action='store_true',
                        help='Show second partial derivative images.')

    parser.add_argument('-f', '--fisher', action='store_true',
                        help='Show fisher matrix images.')

    parser.add_argument('--covariance', action='store_true',
                        help='Show covariance matrix elements')

    parser.add_argument('-c', '--correlation', action='store_true',
                        help='Show correlation matrix elements.')

    parser.add_argument('-bm', '--bias-matrix', action='store_true',
                        help='Show bias matrix images.')

    parser.add_argument('-b', '--biases', action='store_true',
                        help='Show bias images.')

    parser.add_argument('--values', action='store_true',
                        help='Show values on top of appropiate images.')

    parser.add_argument('--snr', required=True, type=float,
                        help='Value of noise bias (standard deviation)')

    parser.add_argument('--verbose', action='store_true',
                        help=('Prints parameters of the galaxy used in the'
                              'fisher formalism plus default parameters'
                              'important for the analysis.'))

    parser.add_argument('--all', action='store_true',
                        help='Stores all fisher images in a pdf format.')


    args = parser.parse_args()

    if args.all:
        args.hide = True

    g_parameters = galfun.GParameters(args.project)
    image_renderer = galfun.ImageRenderer(pixel_scale=float(defaults.PIXEL_SCALE),
                                          nx=float(defaults.NX),ny=float(defaults.NY))
    fish = fisher.Fisher(g_parameters=g_parameters,image_renderer=image_renderer, 
                         snr=float(args.snr))
    #pickle.dump(fish.image,open("model_repo.p","wb"))
    plots = draw.Plots(fish=fish, project=args.project,
                       plots_dir=defaults.PLOTS_DIR,
                       hide=args.hide)

    if args.draw_galaxy or args.all:
        plots.galaxy()
    if args.partials or args.all:
        plots.derivatives()
    if args.fisher and args.values or args.all:
        plots.fisherMatrixValues()
    if args.fisher and not args.values or args.all:
        plots.fisherMatrix()
    if args.covariance or args.all:
        plots.covarianceMatrix()
    if args.correlation or args.all:
        plots.correlationMatrix()
    if args.second_partials or args.all:
        plots.secondDerivatives()
    if args.bias_matrix and args.values or args.all:
        plots.biasMatrixValues()
    if args.biases and args.values or args.all:
        plots.biasValues()
    if args.bias_matrix and not args.values or args.all:
        plots.biasMatrix()
    if args.biases and not args.values or args.all:
        plots.biases()

if __name__ == '__main__':
    main()

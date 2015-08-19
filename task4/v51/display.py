#!/usr/bin/env python
"""This program allows the user to display various statistical results (like
plots) for the galaxies generated in generate.py
"""

import argparse

import fisher

import plotsfisher

import defaults

import galfun

import info


def main():
    parser = argparse.ArgumentParser(description=(
        'Display different results and plots from the galaxies in a given'
        'file.'), formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

    parser.add_argument('--fisher-chi2', action='store_true',
                        help='Show fisher matrix images with chi2 values')

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
                        help='Stores all images in a pdf format. Not showing'
                        'them.')

    parser.add_argument('--info-file', action='store_true',
                        help='Write relevant information to a .txt file for'
                             'future reference')

    parser.add_argument('--bplot1', action='store_true',
                        help='Plot the bias of a given galaxy for a fixed snr'
                             'and varying size.')

    parser.add_argument('--bplot2', action='store_true',
                        help='Plot the bias of a given galaxy for a fixed sized and varying snr.')

    parser.add_argument('--bplot3', action='store_true',
                        help='Plot the bias of a given galaxy for a fixed snr and varying size.')

    parser.add_argument('--bplot4', action='store_true',
                        help='Vary distance between two given galaxies and plot their bias.')

    parser.add_argument('--bplot5', action='store_true',
                        help='Vary the beta of a galaxy (keeping e fixed) and plot bias.')

    parser.add_argument('--bplot6', action='store_true',
                        help='Keep ellipticies of galaxies fixed but change their relative orientation while plotting their bias in a color map.')

    args = parser.parse_args()

    if args.all:
        args.hide = True

    g_parameters = galfun.GParameters(args.project)
    fish = fisher.Fisher(g_parameters=g_parameters, snr=float(args.snr))
    fisherplots = plotsfisher.FisherPlots(fish=fish, project=args.project,
                                          plots_dir=defaults.PLOTS_DIR,
                                          hide=args.hide)

    if args.draw_galaxy or args.all:
        fisherplots.galaxy()
    if args.partials or args.all:
        fisherplots.derivatives()
    if args.fisher and args.values or args.all:
        fisherplots.fisherMatrixValues()
    if args.fisher and not args.values or args.all:
        fisherplots.fisherMatrix()
    if args.fisher_chi2 or args.all:
        fisherplots.fisherMatrixChi2()
    if args.covariance or args.all:
        fisherplots.covarianceMatrix()
    if args.correlation or args.all:
        fisherplots.correlationMatrix()
    if args.second_partials or args.all:
        fisherplots.secondDerivatives()
    if args.bias_matrix and args.values or args.all:
        fisherplots.biasMatrixValues()
    if args.biases and args.values or args.all:
        fisherplots.biasValues()
    if args.bias_matrix and not args.values or args.all:
        fisherplots.biasMatrix()
    if args.biases and not args.values or args.all:
        fisherplots.biases()
    if args.bplot1:
        pass
    if args.bplot2:
        pass
    if args.bplot3:
        pass

    if g_parameters.num_galaxies == 2:
        if args.bplot4:
            fisherplots.biasPlot4()
        if args.bplot5:
            fisherplots.biasPlot5()
        if args.bplot6:
            fisherplots.biasPlot6()
    else:
        if args.bplot4 or args.bplot5 or args.bplot6:
            raise ValueError('The bias plot 4,5,6 only work when you have'
                             'exactly two galaxies specified in your file.')

    information = info.Info(g_parameters, fish)
    if args.verbose:
        information.printInfo()

    if args.info_file or args.all:
        information.writeInfo(args.project)

if __name__ == '__main__':
    main()

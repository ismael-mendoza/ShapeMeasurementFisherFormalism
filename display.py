
#!/usr/bin/env python
"""This program allows the user to display various statistical results (like
plots) for the galaxies generated in generate.py
"""

import argparse

import fisher

import draw

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

    parser.add_argument('--error_bars', action='store_true',
                        help='Show error_bars on bplots.')

    parser.add_argument('--bias_sigma', action='store_true',
                        help='Show bias over sigma instead of bias in bplots')

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

    parser.add_argument('--bplot1', action='store_true',
                        help='Plot the bias of a given galaxy for a fixed snr'
                             'and varying size.')

    parser.add_argument('--bplot2', action='store_true',
                        help='Plot the bias*snr2 of a given galaxy for'
                             'varying snr, keeping sized fixed.')

    parser.add_argument('--bplot3', action='store_true',
                        help='Plot the bias*(snr_norm/snr)**2 of a given'
                        'galaxy as a function of gal_hlr/psf_fwhm.')

    parser.add_argument('--bplot4', action='store_true',
                            help='Plot the bias/sigma of a given galaxy for a'
                            'fixed snr and varying size.')

    parser.add_argument('--bplot5', action='store_true',
                        help='Vary distance between two given galaxies and'
                        'plot their bias.')

    parser.add_argument('--bplot6', action='store_true',
                        help='Vary the beta of a galaxy (keeping e fixed) and'
                        'plot bias.')

    parser.add_argument('--bplot7', action='store_true',
                        help='Keep ellipticies of galaxies fixed but change'
                        'their relative orientation while plotting their bias'
                        'in a color map.')

    parser.add_argument('--noisy_image', action='store_true',
                        help='Display galaxy(ies) with noise just for'
                        'just for visualization of noise.')

    parser.add_argument('--condition_number', action='store_true',
                        help='Print the condition number of the Fisher'
                             'Matrix.')

    parser.add_argument('--ring_test', type=float, metavar='SHEAR',
                        help='Print the bias of the given shear g using'
                              'a ring test on the provided single galaxy.')


    args = parser.parse_args()

    if args.all:
        args.hide = True

    g_parameters = galfun.GParameters(args.project)
    fish = fisher.Fisher(g_parameters=g_parameters, snr=float(args.snr))
    plots = draw.Plots(fish=fish, project=args.project,
                              plots_dir=defaults.PLOTS_DIR,
                              hide=args.hide,
                              error_bars=args.error_bars,
                              bias_sigma=args.bias_sigma)

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

    if g_parameters.num_galaxies == 1:
        if args.bplot1:
            plots.biasPlot1()
        if args.bplot2:
            plots.biasPlot2()
        if args.bplot3:
            plots.biasPlot3()
        if args.bplot4:
            pass
    else:
        if args.bplot1 or args.bplot2 or args.bplot3 or args.bplot4:
            raise ValueError(
                'The bias plots 1,2,3,4 only work with an individual galaxy.')

    if g_parameters.num_galaxies == 2:
        if args.bplot5:
            plots.biasPlot5()
        if args.bplot6:
            plots.biasPlot6()
        if args.bplot7:
            plots.biasPlot7()
    else:
        if args.bplot5 or args.bplot6 or args.bplot7:
            raise ValueError('The bias plot 5,6,7 only work when you have'
                             'exactly two galaxies specified in your file.')


    #extra visualizations for presentations,
    if args.noisy_image:
        plots.noisy_image()
    if args.condition_number:
        print fish.fisher_condition_number
    if args.ring_test:
        print fisher.ringTest(fish, args.ring_test)

    information = info.Info(g_parameters, fish)
    information.writeInfo(args.project)
    if args.verbose:
        information.printInfo()

if __name__ == '__main__':
    main()

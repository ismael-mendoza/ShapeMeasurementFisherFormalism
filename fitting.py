#!/usr/bin/env python
"""Interface that allows user to do N fittings of a galaxy produced in
generate.py given a SNR and compare with the fisher formalism by
displaying biases and correlation coefficients,etc. and also by displaying
a triangle plot that summarizes this results and displays data points.
"""

import argparse

import os

import info

import galfun

import fisher

import defaults

def main():

    parser = argparse.ArgumentParser(
        description=('Displays a triangle plot comparing the expected'
                     'noise bias calculation with the fisher matrix'
                     'formalism for a given file of galaxies.'),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-p', '--project', default=defaults.PROJECT,
    type=str,
    help=('Specify a directory name where the project will be saved.'
          'In this fashion each individual project represents an'
          'analysis.'))

    parser.add_argument('--snr', default=None,
    type=float,
    help=('Specify signal to noise ratio of the run. No need to specify'
    'it twice between doing fits and producing results.'))

    parser.add_argument('-n', '--number-fits', default=1,
    type=int,
    help='Specify how many fits to run.')

    parser.add_argument('-rf', '--run-fits', action='store_true',
    help = ('Run N different instantiations of noise of the added'
            'galaxies to produce a triangle plot. Fits all galaxies in'
            'given file N times.'))

    parser.add_argument('-rfs', '--run-fits-slac',
    metavar='SLAC_COMPUTER',
    help = 'Same as above but have to be logged in a SLAC computer.')

    parser.add_argument('-pr', '--produce-results', action='store_true',
    help=('Read each of the files produces from a fits trial in project'
          'directory and creates a triangle plot comparing fisher to'
          'residuals.'))

    parser.add_argument('--verbose', action='store_true',
    help=('Prints technical information about the galaxy and the defaults'
            'used for the fitting'))

    args = parser.parse_args()

    rltsdir = os.path.join(args.project, defaults.RESULTS_DIR)
    snr_file_name = os.path.join(args.project, defaults.SNR_FILE)

    #delete rltsdir if it exsits and if snr is specified.
    if args.snr:
        snr = args.snr
        if os.path.isdir(rltsdir):
            os.system('rm -r ' + rltsdir)

    elif not args.snr and os.path.isfile(snr_file_name):
        with open(snr_file_name, 'r') as snrfile:
            snr = float(snrfile.readline())

    else:
        raise ValueError('SNR is nowhere to be found and was not'
                         'specified.')

    if not os.path.isdir(rltsdir):
        os.mkdir(rltsdir)

    #count number of existing results in rltsdir.
    existing_fits = 0
    for filename in os.listdir(rltsdir):
        existing_fits += 1

    if args.run_fits:
        for i in range(args.number_fits):
            os.system("python runfits.py " + str(i+1) + " " +
                      str(snr) + " " +  str(args.project) + " "
                      + str(existing_fits))

        #write snr to file, so no confusion as to what snr we have later.
        if args.snr:
            with open(snr_file_name, 'w') as snrfile:
                snrfile.write(str(snr))

    elif args.run_fits_slac:
        os.system("bsub -q " + str(args.run_fits_slac) + " -J \"name[1"
                  + "-" + str(args.number_fits) +
                  "]\" \"python runfits.py \$LSB_JOBINDEX "
                  + str(snr) + " " + str(args.project) + " " +
                  str(existing_fits) + "\"")

        if args.snr:
            with open(snr_file_name, 'w') as snrfile:
                snrfile.write(str(snr))

    elif args.produce_results:
        os.system("python readfits.py " +  str(args.project) + " " +
                  str(snr) + " " + str(existing_fits))

    #produce info file every time we run this program.
    g_parameters = galfun.GParameters(args.project)
    fish = fisher.Fisher(g_parameters, snr)
    image = galfun.drawGalaxies(g_parameters=g_parameters, image=True)
    init_values = defaults.getInitialValuesFit(g_parameters)
    maximums = defaults.getMaximums(g_parameters, image)
    minimums = defaults.getMinimums(g_parameters, image)
    information = info.Info(g_parameters, fish,
                            number_fits = existing_fits, minimums=minimums, maximums=maximums,init_values=init_values)
    information.writeInfo(args.project)

    if args.verbose:
        information.printInfo()

if __name__ == '__main__':
    main()

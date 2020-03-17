"""Interface that allows user to do N fittings of a galaxy produced in
generate.py given a SNR and compare with the fisher formalism by
displaying biases and correlation coefficients,etc. It writes the results into a folder 
inside the project folder specified.
"""

import argparse
import os

import analysis.defaults as defaults


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
                        help=('Run N different instantiations of noise of the added'
                              'galaxies to produce a triangle plot. Fits all galaxies in'
                              'given file N times.'))

    parser.add_argument('-rfs', '--run-fits-slac',
                        metavar='SLAC_COMPUTER',
                        help='Same as above but have to be logged in a SLAC computer.')

    args = parser.parse_args()

    results_dir = os.path.join(args.project, defaults.RESULTS_DIR)
    snr_file_name = os.path.join(args.project, defaults.SNR_FILE)

    # delete results_dir if it exsits and if snr is specified.
    if args.snr:
        snr = args.snr
        if os.path.isdir(results_dir):
            os.system('rm -r ' + results_dir)

    elif not args.snr and os.path.isfile(snr_file_name):
        with open(snr_file_name, 'r') as snrfile:
            snr = float(snrfile.readline())

    else:
        raise ValueError('SNR was not specified.')

    if not os.path.isdir(results_dir):
        os.mkdir(results_dir)

    # count number of existing results in rltsdir.
    existing_fits = 0
    for filename in os.listdir(results_dir):
        existing_fits += 1

    # print(existing_fits)

    if args.run_fits:
        for i in range(args.number_fits):
            os.system("python runfits.py " + str(i + 1) + " " +
                      str(snr) + " " + str(args.project) + " "
                      + str(existing_fits))

        # write snr to file, so no confusion as to what snr we have later.
        if args.snr:
            with open(snr_file_name, 'w') as snrfile:
                snrfile.write(str(snr))

    elif args.run_fits_slac:
        os.system("bsub -o output1.txt -q " + str(args.run_fits_slac) + " -J \"name[1"
                  + "-" + str(args.number_fits) +
                  "]\" \"python runfits.py \$LSB_JOBINDEX "
                  + str(snr) + " " + str(args.project) + " " +
                  str(existing_fits) + "\"")

        if args.snr:
            with open(snr_file_name, 'w') as snrfile:
                snrfile.write(str(snr))


if __name__ == '__main__':
    main()

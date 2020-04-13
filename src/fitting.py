#!/usr/bin/env python3

"""Interface that allows user to do N fittings of a galaxy produced in
generate.py given a SNR and compare with the fisher formalism by
displaying biases and correlation coefficients,etc. It writes the results into a folder 
inside the project folder specified.
"""

import argparse
import shutil
import subprocess
from pathlib import Path

from . import defaults


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

    project_path = Path(args.project)
    assert project_path.exists(), "There should be a project folder with a galaxy in the args.project specified."

    results_dir = project_path.joinpath(defaults.RESULTS_DIR)
    snr_file = project_path.joinpath(defaults.SNR_FILE)

    # delete results_dir if it exists and if snr is specified.
    if args.snr:
        snr = args.snr
        if results_dir.exists():
            shutil.rmtree(results_dir.as_posix())

    elif not args.snr and snr_file.exists():
        with open(snr_file, 'r') as snrfile:
            snr = float(snrfile.readline())

    else:
        raise ValueError('SNR was not specified.')

    results_dir.mkdir(exist_ok=True)

    # count number of existing results in results_dir.
    existing_fits = 0
    for _ in results_dir.iterdir():
        existing_fits += 1

    if args.run_fits:
        for i in range(args.number_fits):
            subprocess.run(f"python -m src.runfits {i + 1} {snr} {project_path} {existing_fits}",
                           shell=True)

        # write snr to file, so no confusion as to what snr we have later.
        if args.snr:
            with open(snr_file, 'w') as snrfile:
                snrfile.write(str(snr))

    elif args.run_fits_slac:
        subprocess.run(f'bsub -o data/output.txt -q {args.run_fits_slac} -J "name[1-{args.number_fits}]"'
                       f' "python -m src.runfits \$LSB_JOBINDEX {snr} {args.project} {existing_fits}"',
                       shell=True)

        if args.snr:
            with open(snr_file, 'w') as snrfile:
                snrfile.write(str(snr))


if __name__ == '__main__':
    main()

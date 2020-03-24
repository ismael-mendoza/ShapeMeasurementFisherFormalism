"""Generate a galaxy(ies) as specified by the user and saves it to a csv file."""
import argparse
import csv
import os
import shutil
from pathlib import Path

import analysis.defaults as defaults
import analysis.models as models


def csv_is_empty(filename):
    """Check if a given csv file is empty by making sure each of the rows in
    the file does not have a string on it.
    """
    with open(filename, 'r') as f:
        for row in f:
            if len(row) != 0:
                return False
        return True


def main():
    parser = argparse.ArgumentParser(description=('Generate galaxies'
                                                  'specified by the user that'
                                                  'are added to a file.'),
                                     formatter_class=(
                                         argparse.ArgumentDefaultsHelpFormatter))

    parser.add_argument('-p', '--project', default=defaults.PROJECT,
                        type=str,
                        help=('Specify a directory name where the project will'
                              'be saved. In this fashion each individual'
                              'project represents one analysis.'))

    parser.add_argument('-gal', '--id', required=True,
                        type=int,
                        help=('Add a galaxy with given ID to the project'
                              'file.'))

    parser.add_argument('--galaxy-model', required=True,
                        type=str, choices=models.get_all_models(),
                        help='Change the galaxy\'s model.')

    parser.add_argument('--psf_model',
                        type=str, choices=models.get_all_psf_models(),
                        help='Change the psf model.')

    parser.add_argument('--snr', type=float,
                        help='Value of noise bias (standard deviation). If'
                             'given an info file with fisher analysis is created.')

    # add all parameter arguments to the parser.
    for name in models.get_all_parameters():
        parser.add_argument('--' + name, default=None,
                            type=float,
                            help='Add a value for the parameter ' + name + '.')

    args = parser.parse_args()

    project_path = Path(args.project)
    if project_path.exists():
        project_path.rmdir()
    project_path.mkdir()

    # create file for galaxies if it does not exist.
    galaxy_file = project_path.joinpath(defaults.GALAXY_FILE)
    temp_file = project_path.joinpath('temp.csv')

    if galaxy_file.exists():
        # creat a copy to read from and compare.
        shutil.copyfile(galaxy_file, temp_file)

    # write galaxy data to a filename.
    with open(galaxy_file, 'w') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=models.get_fieldnames())
        args_dict = vars(args)

        # extract appropriate entries from dictionary of args.
        row_to_write = {k: v for (k, v) in args_dict.items()
                        if k in models.get_fieldnames()}

        if not temp_file.exists():
            writer.writeheader()
            writer.writerow(row_to_write)

        # otherwise more complicated comparison with previous entry.
        else:
            with open(temp_file, 'r') as tempfile:
                reader = csv.DictReader(tempfile)
                writer.writeheader()
                for row in reader:
                    if int(row['id']) != args.id:
                        writer.writerow(row)
                writer.writerow(row_to_write)
            os.remove(temp_file)


if __name__ == '__main__':
    main()

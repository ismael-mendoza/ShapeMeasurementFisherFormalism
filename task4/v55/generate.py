#!/usr/bin/env python
"""Generate a galaxy as specified by the user and save it to a csv file"""
import argparse

import defaults

import os

import csv

import shutil

import names


def csvIsEmpty(filename):
    """Check if a given csv file is empty by making sure each of the rows in
    the file does not have a string on it.
    """
    with open(filename, 'r') as f:
        for row in f:
            if len(row) != 0:
                return False

        return True


def main():
    # NEED TO DISPLAY information about models available for both psf and
    # galaxy.

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
                        type=str, choices=names.gal_models,
                        help='Change the galaxy\'s model.')

    parser.add_argument('--psf_model', required=True,
                        type=str, choices=names.psf_models,
                        help='Change the psf model.')

    # add all parameter arguments to the parser.
    for name in names.parameters:
        parser.add_argument('--' + name, default=None,
                            type=float,
                            help='Add a value for the parameter ' + name + '.')

    args = parser.parse_args()

    if not os.path.isdir(args.project):
        os.mkdir(args.project)

    # create file for galaxies if it does not exist.
    filename = os.path.join(args.project, defaults.GALAXY_FILE)
    if not os.path.isfile(filename):
        f = open(filename, 'w+')
        f.close()

    tempname = os.path.join(args.project, 'temp' + '.csv')

    # creat a copy to read from and compare.
    shutil.copyfile(filename, tempname)

    # write galaxy data to a filename.
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=names.fieldnames)
        args_dict = vars(args)

        # extract appropiate entries from dictionary of args.
        row_to_write = {k: v for (k, v) in args_dict.iteritems()
                        if k in names.fieldnames}

        if csvIsEmpty(tempname):
            writer.writeheader()
            writer.writerow(row_to_write)

        # otherwise more complicated comparison with previous entry.
        else:
            with open(tempname, 'r') as tempfile:
                reader = csv.DictReader(tempfile)
                writer.writeheader()
                for row in reader:
                    if int(row['id']) != args.id:
                        writer.writerow(row)
                writer.writerow(row_to_write)

        os.remove(tempname)

if __name__ == '__main__':
    main()

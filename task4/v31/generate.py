#!/usr/bin/env python
"""Generate a galaxy as specified by the user and save it to a csv file"""
import argparse

import defaults

import os

import csv

import galfun 

import shutil

import names

import info

def csvIsEmpty(filename):
    """Checks if a given csv file is empty by making sure each of the rows in the file does not have a string on it.
    """
    with open(filename, 'r') as f:
        for row in f:
            if len(row) != 0:
                return False
        else: 
            return True

def main():

    parser = argparse.ArgumentParser(description=('Generate galaxies'
    'specified by the user that are added to a file.'), 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--project', default=defaults.PROJECT, 
    type=str,
    help=('Specify a directory name where the project will be saved. In this fashion each individual project represents one analysis.'))

    parser.add_argument('--id', required = True, default = None, 
    type = int,
    help = ('Add galaxy with id GALAXYID, if GALAXYID already exists in the'
            'table, you will change the parameters of the galaxy with that'
            'GALAXYID but not create a new entry.'))
    parser.add_argument('--model', default = names.galaxy_models[0], 
    type = str, choices = names.galaxy_models,
    help = 'Change the galaxy model.')
    parser.add_argument('--psf_model', default = names.psf_models[0], 
    type = str, choices = names.psf_models,
    help = 'Change the psf model.')

    #add all parameter arguments to the parser. 
    for name in names.parameters:
        parser.add_argument('--' + name, default = defaults.PARAMETERS[name], 
            metavar = name.upper(), type = float, 
            help = 'Change the parameter ' + name + ' from its default value.'
            ) 

    parser.add_argument('--verbose', action='store_true',
    help='Prints parameters of the galaxy created.')


    args = parser.parse_args()

    if not os.path.isdir(args.project):
        os.mkdir(args.project)

    filename = os.path.join(args.project, defaults.GALAXY_FILE)
    if not os.path.isfile(filename):
        f = open(filename, 'w+')
        f.close()
    tempname = os.path.join(args.project, 'temp' + '.csv')

    #creat a copy to read from and compare.
    shutil.copyfile(filename, tempname) 

    #write galaxy data to a filename.
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=names.fieldnames)
        args_dict = vars(args)
        #extract appropiate entries from args_dict
        row_to_write={k:v for (k,v) in args_dict.iteritems() 
                      if k in names.fieldnames}

        #if file is empty just write the new row. 
        if csvIsEmpty(tempname):      
            writer.writeheader()
            writer.writerow(row_to_write)
        #otherwise more complicated comparison.
        else:
            with open(tempname, 'r') as tempfile:
                reader = csv.DictReader(tempfile)
                writer.writeheader()        
                for row in reader:
                    if(int(row['id']) != args.id):
                        writer.writerow(row)
                writer.writerow(row_to_write)

        os.remove(tempname)

        if args.verbose:
            information = info.Info(galfun.GParameters(args.project))
            information.printInfo()
                        
                        
if __name__ == '__main__':
    main()
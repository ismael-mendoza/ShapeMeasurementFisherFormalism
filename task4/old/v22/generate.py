#!/usr/bin/env python
"""Generate a galaxy as specified by the user and save it to a csv file"""
import argparse

import defaults

import os

import csv

import galfun 

import shutil

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
    names = defaults.names()
    dflt_params = defaults.parameters().dict 

    parser = argparse.ArgumentParser(description=('Generate galaxies'
    'specified by the user that are added to a file.'), 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--verbose', action='store_true',
    help='Prints parameters of the galaxy created.')
    parser.add_argument('--wd', default=names.wdir, metavar='DIRECTORY', 
    type=str,
    help=('Specify a directory name where the output will be inputted and' 
          'files read from.'))
    parser.add_argument('--galaxy-file', default = names.galaxy_file, metavar = 'FILENAME', type = str,
    help = 'Specify a file where galaxies will be registered and read from, it is a .csv file. Do not need to append .csv at the end.')
    parser.add_argument('--id', required = True, default = None, metavar = 'GALAXYID', type = int,
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
        parser.add_argument('--' + name, default = dflt_params[name], 
            metavar = name.upper(), type = float, 
            help = 'Change the parameter ' + name + ' from its default value.'
            ) 

    args = parser.parse_args()

    if not os.path.isdir(args.wd):
        os.mkdir(args.wd)

    filename = os.path.join(args.wd, args.galaxy_file + '.csv')
    if not os.path.isfile(filename):
        f = open(filename, 'w+')
        f.close()
    tempname = os.path.join(args.wd, 'temp' + '.csv')

    #creat a copy to read from and compare.
    shutil.copyfile(filename, tempname) 

    #write galaxy data to a filename.
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=names.fieldnames)
        row_to_write = dict(
            id = args.id,
            model = args.model,
            x0 = args.x0,
            y0 = args.y0,
            flux = args.flux,
            hlr = args.hlr,
            e1 = args.e1,
            e2 = args.e2,
            psf_model = args.psf_model,
            psf_flux = args.psf_flux,
            psf_fwhm = args.psf_fwhm
        )

        #if file is empty just write the new row. 
        if csvIsEmpty(tempname):      
            writer.writeheader()
            writer.writerow(row_to_write)
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
            galaxies = galfun.galaxies(wdir= args.wdir, 
                                      galaxy_file= args.galaxy_file)
            info = defaults.info(galaxies)
            for line in info.galaxy:
                print line
                        
                        
if __name__ == '__main__':
    main()
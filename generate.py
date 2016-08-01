#!/usr/bin/env python
"""Generate a galaxy as specified by the user and save it to a csv file"""
import argparse
import defaults
import os
import csv
import shutil
import galfun
import fisher
import info
import models 

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
                        type=str, choices=models.getAllModels(),
                        help='Change the galaxy\'s model.')

    parser.add_argument('--psf_model',
                        type=str, choices=models.getAllPsfModels(),
                        help='Change the psf model.')

    parser.add_argument('--snr', type=float,
                        help='Value of noise bias (standard deviation). If'
                        'given an info file with fisher analysis is created.')

    parser.add_argument('--pixel_scale', type=float,
                        help='Pixel scale to arcsecs to use.')


    parser.add_argument('--nx', type=float,
                        help='Width of pixel stamp to generate.')


    parser.add_argument('--ny', type=float,
                        help='Height of pixel stamp to generate.')

    # add all parameter arguments to the parser.
    for name in models.getAllParameters():
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
        writer = csv.DictWriter(csvfile, fieldnames=models.getFieldnames())
        args_dict = vars(args)

        # extract appropiate entries from dictionary of args.
        row_to_write = {k: v for (k, v) in args_dict.iteritems()
                        if k in models.getFieldnames()}
        
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
    
    # if os.path.isfile(tempname):
    #     print 'yes'
    #     exit(1)
    # print 'hello'
    # print 'hello'


    # #write image data to a filename
    # image_dict = {'pixel_scale':float(args.pixel_scale),
    #               'nx':float(args.nx),
    #               'ny':float(args.ny)
    #               }

    # image_filename = os.path.join(args.project, defaults.IMAGE_FILENAME)
    # with open(image_filename,'w') as csvfile: 
    #     writer = csv.DictWriter(csvfile, fieldnames=image_dict.keys())
    #     writer.writeheader()
    #     writer.writerow(image_dict)


    # if(args.snr):
    #     g_parameters = galfun.GParameters(args.project)
    #     image_renderer = galfun.ImageRenderer(pixel_scale=float(args.pixel_scale),
    #                                           nx=float(args.nx),ny=float(args.ny))
    #     fish = fisher.Fisher(g_parameters=g_parameters,image_renderer=image_renderer, 
    #                          snr=float(args.snr))
    #     information = info.Info(g_parameters, image_renderer, fish)
    #     information.writeInfo(args.project)

if __name__ == '__main__':
    main()

#!/usr/bin/env python
"""Contains various functions used throughout the module, may need to look for a better way to organize them specially drawGalaxy will become 
quite large soon. """

import os

import defaults

cts = defaults.constants()
names = defaults.names()


def csvIsEmpty(filename):
    """checks each row and if any is not empty, then the file is not empty"""

    with open(filename, 'r') as f:
        for row in f:
            if len(row) != 0:
                return False
        else: return True

def SaveFigureToPdf(figure,file_name, dir_name, plotdir_name, hide = True): 

     #save and preview pdf.
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    if not os.path.isdir(os.path.join(dir_name,plotdir_name)):
        os.mkdir(os.path.join(dir_name,plotdir_name))
    file_name = os.path.join(dir_name, plotdir_name, file_name)  #puts slashes in between things.
    figure.savefig(file_name, bbox_inches='tight')

    if not hide: 
        os.system("open " + file_name)






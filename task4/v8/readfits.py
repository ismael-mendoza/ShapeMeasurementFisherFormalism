"""This gathers the residuals and computes the bias from multiple .csv files scattered in a directory. 
Each file will not necessarily have only one row of results (although it will most of the times"""

import csv

import numpy as np 

import os

import sys


def main(argv):
    wdir, galaxy_file, rltsdir = argv[1], argv[2], argv[3]

    ######this is repeated in 3 files. 
    if not os.path.isdir(wdir):
        print ('Directory does not exists')
        return -1

    filename = os.path.join(wdir, galaxy_file + '.csv')
    if not os.path.isfile(filename):
        print('Galaxies file does not exist')
        return -1

    if not os.path.isdir(os.path.join(wdir, rltsdir)):
        os.mkdir(os.path.join(wdir, rltsdir))

    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        lst_params = [] #possible here create dictionary with total number of params of both galaxies.
        for row in reader: 
            lst_params.append(row)

    params = lst_params[0] 
    for key, value in params.iteritems():
        try:
            params[key] = float(value)
        except ValueError:
            pass

    ######


    residuals = dict()
    path = os.path.join(wdir,rltsdir)
    for filename in os.listdir(path):
        with open(os.path.join(path,filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for param in row.keys():
                    #store every residual in a list and then can calculate mean. 
                    residuals[param] = [] 
                    #calculate residual for each row. 
                    residuals[param].append(float(params[param]) - float(row[param]))

    biases = {param:np.mean(residuals[param]) for param in residuals.keys()}
    print biases

if __name__ == "__main__":
    main(sys.argv)
"""This gathers the residuals and computes the bias from multiple .csv files scattered in a directory. 
Each file will not necessarily have only one row of results (although it will most of the times"""

import csv

import numpy as np 

import os

import sys

import functions as fns


def main(argv):
    wdir, galaxy_file, rltsdir = argv[1], argv[2], argv[3]

    params = fns.getParamsCsv(wdir, galaxy_file)

    if not os.path.isdir(os.path.join(wdir, rltsdir)):
        os.mkdir(os.path.join(wdir, rltsdir))

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
"""This file contains the sanity check for the fisher 
formalism in task3.py and it involves averaging over the same galaxy but different 
instantiations of the noise in order to calculate the bias from the residuals and compare to the bias given by fisher.
It appends the results from one fit (that includes the adjusted parameter and the corresponding stderr (StandardError)) into a .csv file, 
need to do this a lot of times.

"""

import sys
import os
import math
import logging
import galsim
from lmfit import minimize, Parameters, Parameter, fit_report
import numpy as np
from constantsTask3 import *
from functionsTask3 import *
import csv

def objfuncDif(params, data, psf_params):
    """Without noise the fitting does not really maek sense"""
    model = drawGalaxy(params.valuesdict(), psf_params= psf_params.valuesdict()).array.ravel()
    return model - data

def objfuncChi(params, data, variance_noise, psf_params): 
    model = drawGalaxy(params.valuesdict(), psf_params = psf_params.valuesdict()).array.ravel()
    return (model - data) / math.sqrt(variance_noise)

def csvIsEmpty(filename):
    #check if header is there so we do not write it again. File is not there if the file is empty. 
    with open(filename, 'r') as f:
        if(f.read() == ''):
            return True
        else:
            return False

def main(argv):

    if(len(argv) == 4):
        results_number = argv[1] #what result folder you want to add your calculations to, should be the same parameter as other trials in the same result folder. 
        trial_number = argv[2] #each trial represents exactly the same conditions but run again to produce multiple outputs.
        results_filename = 'data' + argv[3] + '.csv' #want to use different files each time not the same. so the user should input the different file names by using the numbers from slac.
        noise_seed = int(argv[3])
    else:
        print('you need exacty 3 argument when you run this file')
        print('usage: ./task3res.py, result_number, trial_number, results_filename_number')
        return -1

    results_dir = results_folders + results_number #each result directory accounts to the same parameters (change this when parameters change.)
    trial_dir = trial_folders + trial_number

    gal_params = Parameters()
    gal_params.add('flux', value = 100.)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    gal_params.add('HLR', value = 3.)   # arcsec 
    gal_params.add('e1', value = .5)          #ellipticity: e1 
    gal_params.add('e2', value = -.5)         #ellipticity: e2
    gal_params.add('x0', value = 2.)          #shift in x origin. 
    gal_params.add('y0', value = 2.)          #shift in y

    psf_params = Parameters()
    SNR = 40
    psf_params.add('flux', value = 1.)
    psf_params.add('sigma', value = 1.)

    #get image and variance from list returned
    gal_image, variance_noise = drawGalaxy(gal_params.valuesdict(), psf_params = psf_params.valuesdict(), SNR= SNR, noise_seed = noise_seed)
    print variance_noise

    # #draw resulting galaxy.
    # figure1, subplt= plt.subplots(1,1)
    # figure1.suptitle('Initial Galaxy', fontsize = 20)
    # drawPlot(subplt, gal_image.array)
    # SaveFigureToPdfAndOpen(figure1, 'figure1.png') #this will open for preview because it is the default defined in my mac

    #we initialize to random values and set bounds accordingly if necessary. only for the parameter we want to fit.  
    fit_params = Parameters()
    fit_params.add('sigma', value = 1., min = 0) #the min here is important because the program will crash if you have a negative size.
    fit_params.add('flux', value = 50., min = 0)
    fit_params.add('e1', value = 0.05, min = -.7, max = .7) #can be both positive or ngative 
    fit_params.add('e2', value = 0.1, min = -.7, max = .7) #does not necessarily work always but it is enough for the galaxies we are looking at now.
    fit_params.add('x0', value = 2., min = gal_image.getXMin(), max= gal_image.getXMax())
    fit_params.add('y0', value = 2., min = gal_image.getYMin(), max =gal_image.getYMax())

    #this is the function that obtains the fit. 
    minimize(objfuncChi, fit_params, args=(gal_image.array.ravel(), variance_noise, psf_params))


    ##printing fit.
    # print("")
    # print("Start fit_report:")
    # print(fit_report(params))



    #write results into a csv file for reading it later. not only the results from the fit but also the original parameters. 
    if not os.path.isdir(results_dir):
        os.mkdir(results_dir)
    if not os.path.isdir(os.path.join(results_dir, trial_dir)):
        os.mkdir(os.path.join(results_dir,trial_dir))
    results_filename = os.path.join(results_dir, trial_dir, results_filename)

    with open(results_filename, 'w') as csvfile: #'a' for appending and not 'w' for overwriting.
        row_to_write = dict()
        for parameter_name in fit_params.keys():
            row_to_write[parameter_name + '_value'] = fit_params[parameter_name].value
            row_to_write[parameter_name + '_stderr'] = fit_params[parameter_name].stderr #maybe not needed, it depends apparently, will keep it for now.
        writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
        if(csvIsEmpty(results_filename)):
            writer.writeheader()
        writer.writerow(row_to_write)

    #write original parameters and other useful data to a file, only once per result directory
    with open(os.path.join(results_dir, initial_data_filename), 'w') as csvfile:
        row_to_write= dict()
        for parameter_name in gal_param_names:
            row_to_write[parameter_name] = gal_params[parameter_name].value
        for psf_parameter_name in psf_param_names:
            row_to_write[psf_parameter_name] = psf_params[psf_parameter_name].value
        row_to_write['variance_noise'] = variance_noise
        writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
        writer.writeheader()                
        writer.writerow(row_to_write)

if __name__ == "__main__":
    main(sys.argv)

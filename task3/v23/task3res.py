"""This file contains the sanity check for the fisher formalism in task3.py and it involves averaging over the same galaxy but different instantiations of the noise in order
to calculate the bias from the residuals and compare to the bias given by fisher.
It appends the results from one fit (that includes the adjusted parameter and the corresponding stderr (StandardError) ) into a .csv file, need to do this a lot of times 
then do the residual average. 

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
    with open('fit.csv', 'r') as f:
        if(f.read() == ''):
            return True
        else:
            return False

def main(argv):
    initial_data_filename = 'initial_data.csv'
    results_filename = 'fit.csv'

    orig_params = Parameters()
    orig_params.add('gal_flux', value = 100.)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    orig_params.add('gal_sigma', value = 3.) # arcsec 
    orig_params.add('e1', value = .5) #ellipticity: e1 
    orig_params.add('e2', value = -.5)#ellipticity: e2
    orig_params.add('x0', value = 2.) #shift in x origin. 
    orig_params.add('y0', value = 2.)     #shift in y

    psf_params = Parameters()
    noiseSNR = 70
    psf_params.add('psf_flux', value = 1.)
    psf_params.add('psf_sigma', value = 1.)

    #get image and variance from list returned
    gal_image, variance_noise = drawGalaxy(orig_params.valuesdict(), psf_params = psf_params.valuesdict(), noiseSNR= noiseSNR)

    #write original parameters and other useful data to a file (only do it once)
    if(not os.path.isfile(initial_data_filename)):
        with open(initial_data_filename, 'w') as csvfile:
            row_to_write= dict()
            for parameter_name in orig_params.keys():
                row_to_write['true_' + parameter_name] = orig_params[parameter_name].value
            for psf_parameter_name in psf_params.keys():
                row_to_write['true_' + parameter_name] = orig_params[parameter_name].value
            row_to_write['variance_noise'] = variance_noise
            writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
            writer.writeheader()                
            writer.writerow(row_to_write)


    # #draw resulting galaxy.
    # figure1, subplt= plt.subplots(1,1)
    # figure1.suptitle('Initial Galaxy', fontsize = 20)
    # drawPlot(subplt, gal_image.array)
    # SaveFigureToPdfAndOpen(figure1, 'figure1.png') #this will open for preview because it is the default defined in my mac

    #we initialize to random values and set bounds accordingly if necessary. only for the parameter we want to fit.  
    params = Parameters()
    params.add('gal_sigma', value = 1., min = 0) #the min here is important because the program will crash if you have a negative size.
    params.add('gal_flux', value = 50., min = 0)
    params.add('e1', value = 0.05, min = -.7, max = .7) #can be both positive or ngative 
    params.add('e2', value = 0.1, min = -.7, max = .7) #does not necessarily work always but it is enough for the galaxies we are looking at now.
    params.add('x0', value = 2., min = gal_image.getXMin(), max= gal_image.getXMax())
    params.add('y0', value = 2., min = gal_image.getYMin(), max =gal_image.getYMax())

    #this is the function that obtains the fit. 
    minimize(objfuncChi, params, args=(gal_image.array.ravel(), variance_noise, psf_params))


    ##printing fit.
    # print("")
    # print("Start fit_report:")
    # print(fit_report(params))



    #write results into a csv file for reading it later. not only the results from the fit but also the original parameters. 
    with open(results_filename, 'a') as csvfile: #'a' for appending and not 'w' for overwriting.
        row_to_write = dict()
        for parameter_name in params.keys():
            row_to_write[parameter_name + '_value'] = params[parameter_name].value
            row_to_write[parameter_name + '_stderr'] = params[parameter_name].stderr
        writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
        if(csvIsEmpty(results_filename)):
            writer.writeheader()
        writer.writerow(row_to_write)

if __name__ == "__main__":
    main(sys.argv)

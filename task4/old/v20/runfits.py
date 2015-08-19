#!/usr/bin/env python
"""Runs a fit once in a given galaxy image from generate.py, and writes results to a csv file that can be read from readfits.py"""

import os

import csv

import defaults

import sys

import lmfit

import math

import galfun

class Minimums:
    """min values for fit, may add more as needed"""
    def __init__(self, gal_image):
        cts = defaults.constants()
        self.dict = dict()
        self.dict['x0'] = - gal_image.getXMax() * cts.pixel_scale / 2
        self.dict['y0'] = - gal_image.getYMax() * cts.pixel_scale / 2
        self.dict['flux'] = 0. 
        self.dict['hlr'] = 0.
        self.dict['e1'] = -.7 
        self.dict['e2'] = -.7

class Maximums:
    """max values for fit, may add more as needed"""
    def __init__(self, gal_image):
        cts = defaults.constants()
        self.dict = dict()
        self.dict['x0'] = gal_image.getXMax() * cts.pixel_scale / 2
        self.dict['y0'] = gal_image.getYMax() * cts.pixel_scale / 2
        self.dict['flux'] = None 
        self.dict['hlr'] = None
        self.dict['e1'] = .7 
        self.dict['e2'] = .7

def objFunc(fit_params, data, variance_noise, **kwargs):
    model = galfun.drawGalaxies(fit_params.valuesdict(), 
        **kwargs.itervalues().next()).ravel()
    return (model - data) / math.sqrt(variance_noise)

def main(argv):
    noise_seed, snr, wdir, galaxy_file, rltsdir = (
    int(argv[1]), float(argv[2]), argv[3], argv[4], argv[5])

    if not os.path.isdir(os.path.join(wdir, rltsdir)):
        os.mkdir(os.path.join(wdir, rltsdir))

    galaxies = galfun.galaxies(wdir = wdir, galaxy_file = galaxy_file)

    variance_noise = galaxies.addNoise(snr = snr, noise_seed = noise_seed)

    mins = Minimums(galaxies.image).dict
    maxs = Maximums(galaxies.image).dict
    
    #initialize params for fit with true params. 
    fit_params = lmfit.Parameters()
    for param in galaxies.model_params.keys():
        no_subscript_param = param[:-2]
        fit_params.add(param, value = galaxies.model_params[param], 
                       min = mins[no_subscript_param], 
                       max = maxs[no_subscript_param])

    # nfit_params = dict()
    # for param in galaxies.params: 
    #     if param not in galaxies.model_params.keys():
    #         nfit_params[param] = galaxies.model_params[param]

    #run fitting. 
    lmfit.minimize(objFunc, fit_params, kws=dict(data = 
                                           galaxies.image.array.ravel(), 
                                           variance_noise = variance_noise, 
                                           kwargs = galaxies.params))

    #for now lets only worry about the actual fit values reported. 
    result_filename = os.path.join(wdir,rltsdir,''.join(['results',
                                                        str(noise_seed),'.csv'
                                                        ]))

    #write results of fits into a file. 
    with open(result_filename, 'w') as csvfile:
        row_to_write = fit_params.valuesdict()
        writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
        writer.writeheader()
        writer.writerow(row_to_write)

if __name__ == '__main__':
    main(sys.argv)
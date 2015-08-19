#!/usr/bin/env python
"""Runs a fit once in a given galaxy image from generate.py, and writes results to a csv file that can be read from readfits.py"""

import os

import csv

import defaults

import sys

from lmfit import Parameters, minimize, fit_report

import math

import galfun

def objfunc(fit_params, data, variance_noise, **kwargs): 
    model = galfun.drawGalaxy(fit_params.valuesdict(), **kwargs).array.ravel()
    return (model - data) / math.sqrt(variance_noise)

def main(argv):
    noise_seed, snr, wdir, galaxy_file, rltsdir = int(argv[1]), float(argv[2]), argv[3], argv[4], argv[5]

    if not os.path.isdir(os.path.join(wdir, rltsdir)):
        os.mkdir(os.path.join(wdir, rltsdir))

    galaxies = galfun.galaxies(wdir = wdir, galaxy_file = galaxy_file)

    variance_noise = galaxies.addNoise(snr = snr, noise_seed = noise_seed)

    mins = defaults.min(galaxies.image).dict
    maxs = defaults.max(galaxies.image).dict
    
    #initializa with original parameters for fitting pat says its a good idea.
    fit_params = Parameters()
    for param in galaxies.model_params.keys():
        fit_params.add(param, value = galaxies.model_params[param], min = mins[param], max = maxs[param])

    # nfit_params = dict()
    # for param in galaxies.params: 
    #     if param not in galaxies.model_params.keys():
    #         nfit_params[param] = galaxies.model_params[param]

    minimize(objfunc, fit_params, kws=dict(data = galaxies.image.array.ravel(), variance_noise = variance_noise, kwargs = galaxies.params))


    # #printing fit.
    # print("")
    # print("Start fit_report:")
    # print(fit_report(fit_params))

    #for now lets only worry about the actual fit values reported. 
    result_filename = os.path.join(wdir,rltsdir,'results' +str(noise_seed) + '.csv')
    with open(result_filename, 'w') as csvfile: #'a' for appending and not 'w' for overwriting.
        row_to_write = fit_params.valuesdict()
        writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
        writer.writeheader()
        writer.writerow(row_to_write)

if __name__ == '__main__':
    main(sys.argv)
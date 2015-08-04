#!/usr/bin/env python
"""Runs a fit once in a given galaxy image from generate.py, and writes results to a csv file that can be read from readfits.py"""

import os

import csv

import defaults

import sys

import lmfit

import math

import galfun


def objFunc(fit_params, data, variance_noise, **kwargs):
    model = galfun.drawGalaxies(params = fit_params.valuesdict(), 
                                image = True, **kwargs)
    return ((model - data).array.ravel()) / math.sqrt(variance_noise)


def main(argv):
    noise_seed, snr, project = (int(argv[1]), float(argv[2]), argv[3])

    if not os.path.isdir(os.path.join(project, defaults.RESULTS_DIR)):
        os.mkdir(os.path.join(project, defaults.RESULTS_DIR))

    g_parameters = galfun.GParameters(project)
    orig_image = galfun.drawGalaxies(g_parameters=g_parameters, image = True)

    image,variance_noise = galfun.addNoise(orig_image, snr, noise_seed)

    mins = defaults.getMinimums(g_parameters,image)
    maxs = defaults.getMaximums(g_parameters,image)
    
    #initialize params for fit with true params.
    #separating in fit_params and nfit_params is crucial for some reason. 
    fit_params = lmfit.Parameters()
    for param in g_parameters.model_params.keys():
        no_subscript_param = param[:-2]
        fit_params.add(param, value = g_parameters.model_params[param], 
                       min = mins[no_subscript_param], 
                       max = maxs[no_subscript_param])

    nfit_params = dict()
    for param in g_parameters.params: 
        if param not in g_parameters.model_params.keys():
            nfit_params[param] = g_parameters.params[param]

    lmfit.minimize(objFunc, fit_params, kws=dict(data = image, 
                                            variance_noise = variance_noise, 
                                            **nfit_params))


    # print("")
    # print("Start fit_report:")
    # print(lmfit.fit_report(fit_params))

    #for now lets only worry about the actual fit values reported. 
    result_filename = os.path.join(project,defaults.RESULTS_DIR,
                                   ''.join(['results',str(noise_seed),'.csv'
                                                        ]))

    #write results of fits into a file. 
    with open(result_filename, 'w') as csvfile:
        row_to_write = fit_params.valuesdict()
        writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
        writer.writeheader()
        writer.writerow(row_to_write)

if __name__ == '__main__':
    main(sys.argv)
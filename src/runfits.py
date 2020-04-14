#!/usr/bin/env python3

"""Runs a fit once in a given galaxy image from generate.py, and
writes results to a csv file that can be read from using gparameters.py
"""
import csv
import math
import os
import sys

import lmfit
import numpy as np

from . import defaults
from .analysis import fisher
from .analysis import gparameters
from .analysis import images


def obj_func(fit_params, image_renderer, data, variance_noise, **kwargs):
    gal_model = gparameters.get_galaxies_models(fit_params=fit_params.valuesdict(), **kwargs)
    model = image_renderer.get_image(gal_model)
    return ((model - data).array.ravel()) / math.sqrt(variance_noise)


def perform_fit(g_parameters, image_renderer, snr=20., noise_seed=None, method='leastsq'):
    if noise_seed is None:
        noise_seed = np.random.randint(99999999999999)

    fish = fisher.Fisher(g_parameters=g_parameters, image_renderer=image_renderer, snr=snr)
    orig_image = fish.image

    mins = defaults.get_minimums(g_parameters, orig_image)
    maxs = defaults.get_maximums(g_parameters, orig_image)
    init_values = defaults.get_initial_values_fit(g_parameters)
    nfit_params = g_parameters.nfit_params
    noisy_image, variance_noise = images.add_noise(orig_image, snr, noise_seed)

    fit_params = lmfit.Parameters()
    for param in g_parameters.fit_params:
        fit_params.add(param,
                       value=init_values[param],
                       min=mins[param],
                       max=maxs[param])

    results = lmfit.minimize(obj_func, fit_params, method=method, kws=dict(image_renderer=image_renderer,
                                                                           data=noisy_image,
                                                                           variance_noise=variance_noise,
                                                                           **nfit_params))
    return results


def main(argv):
    current_fit_number, snr, project, existing_fits, slen = (
        int(argv[1]), float(argv[2]), argv[3], int(argv[4]), int(argv[5]))

    assert slen % 2 == 1, "slen should be odd otherwise fit will fail. "

    noise_seed = existing_fits + current_fit_number

    if not os.path.isdir(os.path.join(project, defaults.RESULTS_DIR)):
        os.mkdir(os.path.join(project, defaults.RESULTS_DIR))

    g_parameters = gparameters.GParameters(project)
    image_renderer = images.ImageRenderer(pixel_scale=defaults.PIXEL_SCALE, nx=slen, ny=slen)
    results = perform_fit(g_parameters, image_renderer, snr=snr)

    filename = ''.join([defaults.RESULTS_DIR, str(noise_seed), '.csv'])
    result_filename = os.path.join(project, defaults.RESULTS_DIR, filename)

    # obtain dictionary of the result values that can be written to the csv file.
    results_values = dict()
    for param in results.params:
        results_values[param] = results.params[param].value

        # write results of fits into a file,
    with open(result_filename, 'w') as csvfile:
        row_to_write = results_values
        row_to_write['chi2'] = results.chisqr
        row_to_write['success'] = results.success
        row_to_write['errorbars'] = results.errorbars
        row_to_write['nfev'] = results.nfev
        row_to_write['nvarys'] = results.nvarys
        row_to_write['ndata'] = results.ndata
        row_to_write['nfree'] = results.nfree
        row_to_write['redchi'] = results.redchi
        writer = csv.DictWriter(csvfile, fieldnames=list(row_to_write.keys()))
        writer.writeheader()
        writer.writerow(row_to_write)


if __name__ == '__main__':
    main(sys.argv)

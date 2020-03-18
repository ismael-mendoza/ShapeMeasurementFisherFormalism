"""Runs a fit once in a given galaxy image from generate.py, and
writes results to a csv file that can be read from using galfun.py
"""
# ToDo: move to analysis
import csv
import math
import os
import sys

import lmfit

import analysis.defaults as defaults
import analysis.fisher as fisher
import analysis.readfits as galfun


def objFunc(fit_params, image_renderer, data, variance_noise, **kwargs):
    gal_model = galfun.getGalaxiesModels(fit_params=fit_params.valuesdict(), **kwargs)
    model = image_renderer.getImage(gal_model)
    return ((model - data).array.ravel()) / math.sqrt(variance_noise)


def main(argv):
    # import ipdb; ipdb.set_trace()

    current_fit_number, snr, project, existing_fits = (
        int(argv[1]), float(argv[2]), argv[3], int(argv[4]))

    noise_seed = current_fit_number + existing_fits

    if not os.path.isdir(os.path.join(project, defaults.RESULTS_DIR)):
        os.mkdir(os.path.join(project, defaults.RESULTS_DIR))

    g_parameters = galfun.GParameters(project)
    image_renderer = galfun.ImageRenderer(pixel_scale=defaults.PIXEL_SCALE, nx=defaults.NX, ny=defaults.NY)
    fish = fisher.Fisher(g_parameters=g_parameters, image_renderer=image_renderer, snr=snr)
    # orig_image = copy.deepcopy(fish.image)
    orig_image = fish.model.drawImage(nx=defaults.NX, ny=defaults.NX, scale=defaults.PIXEL_SCALE)
    mins = defaults.get_minimums(g_parameters, orig_image)
    maxs = defaults.getMaximums(g_parameters, orig_image)
    init_values = defaults.get_initial_values_fit(g_parameters)
    nfit_params = g_parameters.nfit_params
    noisy_image, variance_noise = galfun.addNoise(orig_image, snr, noise_seed)

    fit_params = lmfit.Parameters()
    for param in g_parameters.fit_params:
        fit_params.add(param,
                       value=init_values[param],
                       min=mins[param],
                       max=maxs[param])

    results = lmfit.minimize(objFunc, fit_params, kws=dict(image_renderer=image_renderer,
                                                           data=noisy_image,
                                                           variance_noise=(variance_noise),
                                                           **nfit_params))

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

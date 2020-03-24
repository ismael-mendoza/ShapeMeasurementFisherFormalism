"""Multipurpose module that contains important functions ranging from managing parameters of 
generated galaxies to extracting information from relevant files.
"""
import csv
import math
import os
import numpy as np

from .. import defaults
from . import models


def get_omit_fit(id_params, omit):
    omit_fit = {}

    for gal_id in id_params:
        params_omit = omit.get(gal_id, [])
        params = id_params[gal_id]
        galaxy_model = params['galaxy_model']
        cls = models.get_model_cls(galaxy_model)
        obj = cls(params_omit=params_omit)
        omit_fit[gal_id] = obj.omit_fit

    return omit_fit


def read_results(project, g_parameters, fish, limit=None):
    orig_image = fish.image
    mins = defaults.get_minimums(g_parameters, orig_image)
    maxs = defaults.getMaximums(g_parameters, orig_image)

    residuals = {}
    pulls = {}
    redchis = []  # list containing values of reduced chi2 for each fit.
    results_dir = os.path.join(project, defaults.RESULTS_DIR)

    files = os.listdir(results_dir)

    if limit is not None:
        files = files[:limit]

    # read results from results_dir's files.
    for filename in files:
        with open(os.path.join(results_dir, filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for i, row in enumerate(reader):
                redchis.append(float(row['redchi']))
                for param in g_parameters.fit_params:
                    if param not in residuals:
                        residuals[param] = []
                    if param not in pulls:
                        pulls[param] = []
                    residual = (float(row[param]) -
                                float(g_parameters.params[param]))

                    pull = (residual /
                            math.sqrt(fish.covariance_matrix[param, param]))

                    residuals[param].append(residual)
                    pulls[param].append(pull)

    biases = {param: np.mean(residuals[param]) for param in residuals}
    pull_means = {param: np.mean(pulls[param]) for param in residuals}
    res_stds = {param: np.std(residuals[param]) for param in residuals}
    pull_mins = {param: ((mins[param] - float(g_parameters.params[param])) /
                         math.sqrt(fish.covariance_matrix[param, param])) for
                 param in residuals}
    pull_maxs = {param: ((maxs[param] - float(g_parameters.params[param])) /
                         math.sqrt(fish.covariance_matrix[param, param])) for
                 param in residuals}

    return pulls, residuals, biases, pull_means, res_stds, pull_mins, pull_maxs, redchis

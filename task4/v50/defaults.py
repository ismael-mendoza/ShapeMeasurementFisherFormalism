#!/usr/bin/env python
"""Some of the defaults that are used in the overall program."""


def getSteps(g_parameters):
    """Get steps for a given set of galaxy parameters."""
    steps = dict()
    fit_params = g_parameters.fit_params
    for param in fit_params:
        if 'flux' in param or 'hlr' in param:
            steps[param] = fit_params[param] * .01
        else:
            steps[param] = .01

    return steps


def getInitialValuesFit(g_parameters):
    """Get initial values for the fit."""
    initial_values = dict()
    fit_params = g_parameters.fit_params
    for param in fit_params:
        if fit_params[param] == 0:
            initial_values[param] = FIT_DEVIATION
        else:
            initial_values[param] = fit_params[param]

    return initial_values


def getMinimums(g_parameters, gal_image):
    """min values for fit, may add more as needed."""
    minimums = dict()
    for param in g_parameters.fit_params:
        if 'flux' in param:
            minimums[param] = 0.

        elif 'hlr' in param:
            minimums[param] = 0.

        elif 'x0' in param:
            minimums[param] = - gal_image.getXMax() * PIXEL_SCALE / 2

        elif 'y0' in param:
            minimums[param] = - gal_image.getYMax() * PIXEL_SCALE / 2

        elif 'eta1' in param:
            minimums[param] = g_parameters.params[param] - 2.5

        elif 'eta2' in param:
            minimums[param] = g_parameters.params[param] - 2.5

        elif 'e1' in param:
            minimums[param] = -.7

        elif 'e2' in param:
            minimums[param] = -.7

    return minimums


def getMaximums(g_parameters, gal_image):
    """Max values for fit, may add more as needed."""
    maximums = dict()
    for param in g_parameters.fit_params:

        if 'flux' in param:
            maximums[param] = float('Inf')

        elif 'hlr' in param:
            maximums[param] = float('Inf')

        elif 'x0' in param:
            maximums[param] = gal_image.getXMax() * PIXEL_SCALE / 2

        elif 'y0' in param:
            maximums[param] = gal_image.getYMax() * PIXEL_SCALE / 2

        elif 'eta1' in param:
            maximums[param] = g_parameters.params[param] + 2.5

        elif 'eta2' in param:
            maximums[param] = g_parameters.params[param] + 2.5

        elif 'e1' in param:
            maximums[param] = .7

        elif 'e2' in param:
            maximums[param] = .7

    return maximums


# general global(module-level) constants.
NX = 48  # pixels
NY = 48
PIXEL_SCALE = .2  # arcsec/pixel
DPI = 600.  # resolution for pdf saving.
FONTSIZE_LABEL = 8
FONTSIZE_VALUE = 4
FONTSIZE_TITLE = 14
EXTENT_PULL = (-3, 3)
BINS_PULL = 40
SIG_DIGITS = 4
FIT_DEVIATION = .00001


# some default names for argparse and i/0 file management.
PROJECT = 'project'
PLOTS_DIR = 'plots'
RESULTS_DIR = 'results'
RESULTS_BASENAME = 'result'
GALAXY_FILE = 'galaxies.csv'
PSF_FILE = 'psf.csv'
INFO_FILE = 'info.txt'
FIGURE_BASENAME = 'figure'
FIGURE_EXTENSION = '.pdf'
TRIANGLE_NAME = 'triangle.pdf'
REDCHI_HIST_NAME = 'redchi_hist.pdf'
SNR_FILE = 'snr.txt'
MODEL = 'gaussian'
REDCHI_HIST_TITLE = 'Reduced chi2 histogram from data.'

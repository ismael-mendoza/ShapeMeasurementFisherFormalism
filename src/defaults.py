"""Some of the defaults that are used in the overall program."""

import numpy as np


def get_steps(g_parameters, image_renderer):
    """Return a dictionary containing the steps to be used in the
    derivatives of each parameter.

    The dictionary is of the form: 'parameter_name:value_of_step'

    Some parameter variations were copied from David's code suggestions. 

    Args:
    g_parameters(:class:`analysis.galfun.GParameters`): An object containing different
        forms of the galaxy parameters.

    Returns:
        A dict.
    """
    steps = dict()
    fit_params = g_parameters.fit_params
    for param in fit_params:
        if 'flux' in param:
            steps[param] = fit_params[param] * .01

        elif 'hlr' in param:
            steps[param] = fit_params[param] * .05

        elif 'beta' in param:
            steps[param] = .1

        elif 'x0' in param or 'y0' in param:
            steps[param] = image_renderer.pixel_scale / 3.

        elif 'g1' in param or 'g2' in param:
            steps[param] = .03

        elif 'e' in param:
            steps[param] = .03

        else:
            steps[param] = .01

    return steps


def get_initial_values_fit(g_parameters):
    """Return a dictionary containing the initial values to be used in the
    in the fitting of the parameters.

    The dictionary is of the form: 'parameter_name:initial_value'

    Args:
    g_parameters(:class:`analysis.galfun.GParameters`): An object containing different
        forms of the galaxy parameters.

    Returns:
        A dict.
    """
    initial_values = dict()
    fit_params = g_parameters.fit_params
    for param in fit_params:
        initial_values[param] = fit_params[param] + abs(np.random.uniform()) * (fit_params[param] / 10 + 0.2)
    return initial_values


def get_minimums(g_parameters, gal_image):
    """Return a dictionary containing the minimum values to be used in the
    in the fitting of the parameters.

    The dictionary is of the form: 'parameter_name:initial_value'

    Args:
    g_parameters(:class:`analysis.galfun.GParameters`): An object containing different
        forms of the galaxy parameters.

    Returns:
        A dict.
    """
    minimums = dict()
    for param in g_parameters.fit_params:
        if 'flux' in param:
            minimums[param] = 0.

        elif 'hlr' in param:
            minimums[param] = 0.

        elif 'x0' in param:
            minimums[param] = - gal_image.xmax * PIXEL_SCALE / 2

        elif 'y0' in param:
            minimums[param] = - gal_image.ymax * PIXEL_SCALE / 2

        elif 'eta1' in param:
            minimums[param] = g_parameters.params[param] - 2.5

        elif 'eta2' in param:
            minimums[param] = g_parameters.params[param] - 2.5

        elif 'e1' in param or 'g1' in param:
            minimums[param] = -.7

        elif 'e2' in param or 'g2' in param:
            minimums[param] = -.7

    return minimums


def get_maximums(g_parameters, gal_image):
    """Return a dictionary containing the minimum values to be used in the
    in the fitting of the parameters.

    The dictionary is of the form: 'parameter_name:maximum_value'

    Args:
    g_parameters(:class:`galfun.GParameters`): An object containing different
        forms of the galaxy parameters.

    Returns:
        A dict.
    """
    maximums = dict()
    for param in g_parameters.fit_params:

        if 'flux' in param:
            maximums[param] = float('Inf')

        elif 'hlr' in param:
            maximums[param] = float('Inf')

        elif 'x0' in param:
            maximums[param] = gal_image.xmax * PIXEL_SCALE / 2

        elif 'y0' in param:
            maximums[param] = gal_image.ymax * PIXEL_SCALE / 2

        elif 'eta1' in param:
            maximums[param] = g_parameters.params[param] + 2.5

        elif 'eta2' in param:
            maximums[param] = g_parameters.params[param] + 2.5

        elif 'e1' in param or 'g1' in param:
            maximums[param] = .7

        elif 'e2' in param or 'g2' in param:
            maximums[param] = .7

    return maximums


# general global(module-level) constants.
FIT_DEVIATION = .00001
SNR_NORM = 20.
NX = 40  # number of pixels.
NY = 40
PIXEL_SCALE = .2
SIG_DIGITS = 4
DPI = 300
FONTSIZE_LABEL = 8
FONTSIZE_VALUE = 4
FONTSIZE_TITLE = 14

# some default names for argparse and i/0 file management.
PROJECT = 'project'
PLOTS_DIR = 'plots'
RESULTS_DIR = 'results'
GALAXY_FILE = 'galaxies.csv'
SNR_FILE = 'snr.txt'
MODEL = 'gaussian'
FIGURE_BASENAME = 'figure'
FIGURE_EXTENSION = '.pdf'

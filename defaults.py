
#!/usr/bin/env python
"""Some of the defaults that are used in the overall program."""


def getSteps(g_parameters):
    """Return a dictionary containing the steps to be used in the
    derivatives of each parameter.

    The dictionary is of the form: 'parameter_name:value_of_step'

    Some parameter variations were copied from David's code suggestions. 

    Args:
    g_parameters(:class:`GParameters`): An object containing different
                                        forms of the galaxy parameters.

    Returns:
        A :py:dict.
    """
    steps = dict()
    fit_params = g_parameters.fit_params
    for param in fit_params:
        if 'flux' in param:
            steps[param] = fit_params[param] * .01

        if 'hlr' in param:
            steps[param] = fit_params[param] * .05

        if 'x' in param or 'y' in param:
            steps[param] = defaults.PIXEL_SCALE/3.

        if 'g1' in param or 'g2' in param:
            steps[param] = .03

        else:
            steps[param] = .01

    return steps


def getInitialValuesFit(g_parameters):
    """Return a dictionary containing the initial values to be used in the
    in the fitting of the parameters.

    The dictionary is of the form: 'parameter_name:initial_value'

    Args:
    g_parameters(:class:`GParameters`): An object containing different
                                        forms of the galaxy parameters.

    Returns:
        A :py:dict.
    """
    initial_values = dict()
    fit_params = g_parameters.fit_params
    for param in fit_params:
        if fit_params[param] == 0:
            initial_values[param] = FIT_DEVIATION
        else:
            initial_values[param] = fit_params[param]

    return initial_values


def getMinimums(g_parameters, gal_image):
    """Return a dictionary containing the minimum values to be used in the
    in the fitting of the parameters.

    The dictionary is of the form: 'parameter_name:initial_value'

    Args:
    g_parameters(:class:`GParameters`): An object containing different
                                        forms of the galaxy parameters.

    Returns:
        A :py:dict.
    """
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
    """Return a dictionary containing the minimum values to be used in the
    in the fitting of the parameters.

    The dictionary is of the form: 'parameter_name:maximum_value'

    Args:
    g_parameters(:class:`GParameters`): An object containing different
                                        forms of the galaxy parameters.

    Returns:
        A :py:dict.
    """
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
DPI = 600.  # resolution for pdf saving.
FONTSIZE_LABEL = 8
FONTSIZE_VALUE = 4
FONTSIZE_TITLE = 14
EXTENT_PULL = (-3, 3)
BINS_PULL = 40
SIG_DIGITS = 4
FIT_DEVIATION = .00001
SNR_NORM = 20.
FANCY = False


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

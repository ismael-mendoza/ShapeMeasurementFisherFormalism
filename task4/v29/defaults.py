#!/usr/bin/env python
"""Some of the defaults that are used in the overall program."""

def getSteps(g_parameters): 
    """Get steps for a given set of galaxy parameters"""

    steps = dict()
    for param in g_parameters.model_params.keys():
        param_no_subscript = param[:-2]
        if param_no_subscript == 'flux' or param_no_subscript == 'hlr': 
            steps[param] = g_parameters.model_params[param] * .01
        else:
            steps[param] = .01

    return steps


def getMinimums(g_parameters, gal_image):
    """min values for fit, may add more as needed"""

    minimums = dict()
    for param in g_parameters.model_param_names:
        param_no_subscript = param[:-2]
        if param_no_subscript == 'flux':
            minimums['flux'] = 0. 

        elif param_no_subscript == 'hlr':
            minimums['hlr'] = 0. 

        elif param_no_subscript == 'x0':
            minimums['x0'] = - gal_image.getXMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'y0':
            minimums['y0'] = - gal_image.getYMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'e1':
            minimums['e1'] = -.7 

        elif param_no_subscript == 'e2':
            minimums['e2'] = -.7

    return minimums


def getMaximums(g_parameters, gal_image):
    """max values for fit, may add more as needed"""

    maximums = dict()
    for param in g_parameters.model_param_names:
        param_no_subscript = param[:-2]
        if param_no_subscript == 'flux':
            maximums['flux'] = None

        elif param_no_subscript == 'hlr':
            maximums['hlr'] = None 

        elif param_no_subscript == 'x0':
            maximums['x0'] = gal_image.getXMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'y0':
            maximums['y0'] = gal_image.getYMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'e1':
            maximums['e1'] = .7 

        elif param_no_subscript == 'e2':
            maximums['e2'] = .7

    return maximums


#general global(module-level) constants.
NX = 48 #pixels 
NY = 48
PIXEL_SCALE = .2 #arcsec/pixel 
DPI = 600. #resolution for pdf saving.
FONTSIZE_LABEL = 8
FONTSIZE_VALUE = 6
FONTSIZE_TITLE = 14
SNR = 60. 
EXTENT_PULL = (-2.5,2.5)
BINS_PULL = 40
SIG_DIGITS = 4


#some default names for argparse and i/0 file management.
PROJECT = 'project'
PLOTS_DIR = 'plots'
RESULTS_DIR = 'results'
GALAXY_FILE = 'galaxies.csv'
INFO_FILE = 'info.txt'
FIGURE_BASENAME = 'figure'
FIGURE_EXTENSION = '.pdf'
TRIANGLE_NAME = 'triangle.pdf'
SNR_FILE = 'snr.txt'


#Default values of some parameters that are given to argparse.
PARAMETERS = dict()
PARAMETERS['x0'] = 0. 
PARAMETERS['y0'] = 0. 
PARAMETERS['flux'] = 100. 
PARAMETERS['hlr'] = 1.
PARAMETERS['e1'] = 0. 
PARAMETERS['e2'] = 0. 
PARAMETERS['psf_flux'] = 1. 
PARAMETERS['psf_fwhm'] = .7 
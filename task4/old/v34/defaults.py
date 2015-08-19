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
            minimums[param] = 0. 

        elif param_no_subscript == 'hlr':
            minimums[param] = 0. 

        elif param_no_subscript == 'x0':
            minimums[param] = - gal_image.getXMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'y0':
            minimums[param] = - gal_image.getYMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'e1':
            minimums[param] = -1 

        #no bounds on e2 as it is constrained by the expr. 
        elif param_no_subscript == 'e2':
            minimums[param] = None

    return minimums


def getMaximums(g_parameters, gal_image):
    """max values for fit, may add more as needed"""

    maximums = dict()
    for param in g_parameters.model_param_names:
        param_no_subscript = param[:-2]
        if param_no_subscript == 'flux':
            maximums[param] = None

        elif param_no_subscript == 'hlr':
            maximums[param] = None 

        elif param_no_subscript == 'x0':
            maximums[param] = gal_image.getXMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'y0':
            maximums[param] = gal_image.getYMax() * PIXEL_SCALE / 2

        elif param_no_subscript == 'e1':
            maximums[param] = 1

        elif param_no_subscript == 'e2':
            maximums[param] = None

    return maximums

def getExprs(g_parameters):

    exprs = dict()
    #assume both galaxies have same params (always, which they use depends on 
    #their model.)
    first_galaxy_params = g_parameters.id_params.itervalues().next()
    for param in first_galaxy_params:
        for gal_id in g_parameters.id_params:
            subscriptParam = param + '_' + gal_id
            if param == 'flux':
                exprs[subscriptParam] = None

            elif param == 'hlr':
                exprs[subscriptParam] = None 

            elif param == 'x0':
                exprs[subscriptParam] = None

            elif param == 'y0':
                exprs[subscriptParam] = None

            elif param == 'e1':
                exprs[subscriptParam] = None

            elif param == 'e2':
                exprs[subscriptParam] = ('math.sqrt' + '(e_' + gal_id + 
                                         '**2-e1**2)')

    return exprs


#general global(module-level) constants.
NX = 48 #pixels 
NY = 48
PIXEL_SCALE = .2 #arcsec/pixel 
DPI = 600. #resolution for pdf saving.
FONTSIZE_LABEL = 8
FONTSIZE_VALUE = 4
FONTSIZE_TITLE = 14
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
PARAMETERS['x0'] = 1. 
PARAMETERS['y0'] = 1. 
PARAMETERS['flux'] = 100. 
PARAMETERS['hlr'] = 1.
PARAMETERS['e1'] = 0. 
PARAMETERS['e2'] = 0. 
PARAMETERS['psf_flux'] = 1. 
PARAMETERS['psf_fwhm'] = .7 
PARAMETERS['psf_e1'] = 0.
PARAMETERS['psf_e2'] = 0.
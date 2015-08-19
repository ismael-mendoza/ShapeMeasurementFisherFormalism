#!/usr/bin/env python
"""Some of the defaults that are used in the overall program."""

import csv

import os


def getSteps(g_parameters): 
    """Get steps for a given set of galaxy parameters"""

    steps = dict()
    for param in g_parameters.param.keys():
        param_no_subscript = param[:-2]
        if param_no_subscript == 'flux' or param_no_subscript == 'hlr': 
            steps[param] = [param] * .01
        else:
            steps[param] = .01

    return steps


def getMinimums(gal_image):
    """min values for fit, may add more as needed"""

    minimums = dict()
    minimums['x0'] = - gal_image.getXMax() * PIXEL_SCALE / 2
    minimums['y0'] = - gal_image.getYMax() * PIXEL_SCALE / 2
    minimums['flux'] = 0. 
    minimums['hlr'] = 0.
    minimums['e1'] = -.7 
    minimums['e2'] = -.7

    return minimums


def getMaximums(gal_image):
    """max values for fit, may add more as needed"""

    maximums = dict()
    maximums['x0'] = gal_image.getXMax() * PIXEL_SCALE/ 2
    maximums['y0'] = gal_image.getYMax() * PIXEL_SCALE / 2
    maximums['flux'] = None 
    maximums['hlr'] = None
    maximums['e1'] = .7 
    maximums['e2'] = .7
    return maximums


class Names:
    def __init__(self):

        self.galaxy_models = ['gaussian','exponential']
        self.psf_models = ['gaussian']

        self.gal_models_parameters = dict()
        self.gal_models_parameters['gaussian'] = ['x0', 'y0', 'flux', 'hlr', 
                                                  'e1', 'e2'
                                                 ]
        self.gal_models_parameters['exponential'] = ['x0', 'y0', 'flux', 
                                                     'hlr', 'e1', 'e2'
                                                    ]
                                                    
        #include all of the parameters in each of the models in the order 
        #you desire to display them in. This is done manually just
        #because the order we want is very specific.
        self.gal_parameters = ['x0', 'y0', 'flux', 'hlr', 'e1', 'e2']

        self.psf_models_parameters = dict()
        self.psf_models_parameters['gaussian'] = ['psf_flux','psf_fwhm']
        self.psf_parameters =  list(set(
        [elt for sublist in self.psf_models_parameters.values() for elt in sublist]))

        self.parameters = self.gal_parameters + self.psf_parameters
        self.extra = ['id', 'model', 'psf_model']

        self.fieldnames = (self.extra + self.gal_parameters + 
                          self.psf_parameters)


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

#some default names for argparse and i/0 file management.
PROJECT_DIR = 'project'
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


#declarations of objects.
names = Names()
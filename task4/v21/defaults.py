#!/usr/bin/env python
"""Some of the default values that are used in the overall. In general
defaults should not be able to be change by the user. But some
exceptions for steps or max & min may be implemented in the future.
"""

import math

def dictToStringList(dct): 
    """Returns a list of strings consisting of the form 'key:value' for a
    given dictionary."""
    if type(dct.keys()[0]) == str: 
        return [str(key) + ': ' + str(dct[key]) for key in dct.keys()]
    else:
        return [str(','.join(key)) + ': ' + str(dct[key]) 
                for key in dct.keys()
               ]


class constants:
    """General global constants."""
    def __init__(self):
        self.nx = 48 #pixels 
        self.ny = 48
        self.pixel_scale = .2 #arcsec/pixel 
        self.snr = 60. 
        self.dpi = 600. #resolution for pdf saving.

        #fontsizes for labels and values in plotsfisher.
        self.fontsize_label = 8
        self.fontsize_values = 6 
        self.fontsize_titles = 14
        self.figure_extension = '.pdf'
        self.figure_basename = 'figure'


class parameters:
    """Default values of some parameters that is given to argparse."""
    def __init__(self):
        self.dict = dict()
        self.dict['x0'] = 0. 
        self.dict['y0'] = 0. 
        self.dict['flux'] = 100. 
        self.dict['hlr'] = 1.
        self.dict['e1'] = 0. 
        self.dict['e2'] = 0. 
        self.dict['psf_flux'] = 1. 
        self.dict['psf_fwhm'] = .7 


class names:
    """Contains the different names for parameters, headers, etc. in lists
    that are useful in various situations.
    """
    def __init__(self):

        #some default names for argparse.
        self.wdir = 'output'
        self.galaxy_file = 'galaxies'
        self.plots_dir = 'plots'
        self.rltsdir = 'results'
        self.info = 'info'

        self.galaxy_models = ['gaussian','exponential']
        self.psf_models = ['gaussian']

        self.gal_models_parameters = dict()
        self.gal_models_parameters['gaussian'] = ['x0','y0',
                                                 'flux','hlr','e1', 'e2'
                                                 ]

        self.gal_models_parameters['exponential'] = ['x0','y0', 'flux','hlr',
                                                     'e1', 'e2'
                                                    ]

        #include all of the ones in each model in the order you desire to
        #display them in. This is done manually just because the order we 
        #want is very specific.
        self.galaxy_parameters = ['x0', 'y0', 'flux', 'hlr', 'e1', 'e2']

        self.psf_models_parameters = dict()
        self.psf_models_parameters['gaussian'] = ['psf_flux','psf_fwhm']
        self.psf_parameters =  list(set(
        [elt for sublist in self.psf_models_parameters.values() 
        for elt in sublist
        ]))

        self.parameters = self.galaxy_parameters + self.psf_parameters
        self.extra = ['id', 'model', 'psf_model']
        self.fieldnames = (self.extra + self.galaxy_parameters + 
                          self.psf_parameters)


class info: 
    """Contains the different things that are written in the information text
    file and printed out with verbose option.
    CHANGE STRING CONCATENATION EVENTUALLY TO NOT USE +
    """
    
    def __init__(self, g_parameters, fish = None, fits_biases = None, num_fits = None):
        cts = constants()
        params = g_parameters.params 

        #contains all the lines of galaxy info to be writen into the file. 
        self.galaxy = list(
            (
            'Default values used in the analysis:',
            'nx: ' + str(cts.nx),
            'ny: ' + str(cts.ny),
            'pixel_scale: ' + str(cts.pixel_scale),
            '',
            'Galaxies drawn have the following parameters:'
            ) 
            + 
            tuple(param + ': ' + str(params[param]) for param in params.keys()
                 ) 
        )
        # if g_parameters.num_galaxies == 2:
        #     separation = math.sqrt()
        #     separation_info = ('Separation between galaxies in arcsecs:'
        #                        str(params['x0_1'])

        if fish:
            self.fisher = list(
                (
                '',
                'Fisher matrix elements:',
                )
                + 
                tuple(dictToStringList(fish.fisher_matrix))
                + (
                '',
                'Covariance matrix elements:'
                )
                + 
                tuple(dictToStringList(fish.covariance_matrix))
                +
                (
                '',
                'Biases for each parameter'
                )
                + 
                tuple(dictToStringList(fish.biases))
                + 
                (
                '',
                'Fisher analysis used the following snr): ' + str(fish.snr),
                '',
                'Steps used for the derivatives: '
                ) 
                + 
                tuple(param + ': ' + str(fish.steps[param]) for param in fish.steps.keys())
            )

        #add covariances in the future? (unlikely)
        if fits_biases and num_fits:
            self.fits = list(
                (
                'Number of fits: ' + str(num_fits),
                'Biases were obtained by averaging over fits of noisy instantiations of the same image: '
                )
                +
                tuple(dictToStringList(fits_biases))
            )
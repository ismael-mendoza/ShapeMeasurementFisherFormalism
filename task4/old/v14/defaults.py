#!/usr/bin/env python
"""Some of the default values that are used in the overall.
In general defaults should not be able to be change by the user. But some exceptions for steps or max & min may be implemented in the future.

"""

def dictToStringList(dct): 
    """Returns a list of strings consisting of the form 'key:value' for a given dictionary"""
    if type(dct.keys()[0]) == str: 
        return [str(key) + ': ' + str(dct[key]) for key in dct.keys()]
    else:
        return [str(','.join(key)) + ': ' + str(dct[key]) for key in dct.keys()]

class constants:
    """General global constants"""
    def __init__(self):
        self.nx = 48
        self.ny = 48
        self.pixel_scale = .2
        self.sigma_n = 1.
        self.snr = 60.

class parameters:
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

class steps: 
    """Define the steps for derivatives of each individual parameter."""
    def __init__(self, params):
        self.dict = dict()
        self.dict['flux'] = params['flux'] * .01
        self.dict['hlr'] = params['hlr'] *.01
        self.dict['e1'] = .01
        self.dict['e2'] = .01
        self.dict['x0'] = .01
        self.dict['y0'] = .01
        #steps['q'] = orig_params['q']*.01
        #steps['beta'] = orig_params['beta']*.01

class min:
    """min values for fit, may add more as needed"""
    def __init__(self, gal_image):
        self.dict = dict()
        self.dict['x0'] = 0.
        self.dict['y0'] = 0.
        self.dict['flux'] = 0. 
        self.dict['hlr'] = 0.
        self.dict['e1'] = -.7 
        self.dict['e2'] = -.7

class max:
    """max values for fit, may add more as needed"""
    def __init__(self, gal_image):
        cts = constants()
        self.dict = dict()
        self.dict['x0'] = gal_image.getXMax() * cts.pixel_scale / 2
        self.dict['y0'] = gal_image.getYMax() * cts.pixel_scale / 2
        self.dict['flux'] = None 
        self.dict['hlr'] = None
        self.dict['e1'] = .7 
        self.dict['e2'] = .7

class names:
    """Contains the different names for parameters, headers, etc. in lists"""
    def __init__(self):
        self.galaxy_models = ['gaussian','exponential']
        self.psf_models = ['gaussian']

        self.galaxy_parameters = dict()
        self.psf_parameters = dict()
        self.galaxy_parameters['gaussian'] = ['x0','y0','flux','hlr','e1', 'e2']
        self.psf_parameters['gaussian'] = ['psf_flux','psf_fwhm']
        self.galaxy_parameters['exponential'] = ['x0','y0', 'flux','hlr','e1', 'e2']

        self.parameters = [] 
        for lst in self.galaxy_parameters.values() + self.psf_parameters.values():
            for element in lst:
                self.parameters.append(element)
        self.parameters = list(set(self.parameters)) #remove repeated elements.

        self.fieldnames = ['id', 'model', 'psf_model'] + self.parameters


        #some for argparse.
        self.wdir = 'output'
        self.galaxy_file = 'galaxies'
        self.plots_dir = 'plots'
        self.rltsdir = 'results'
        self.info = 'info'

class info: 
    """Contains the different things that are written in the information text file"""
    
    def __init__(self, params, fish = -1, fits_biases = -1, N = -1):
        cts = constants()


        #contains all the lines of galaxy info to be writen into the file. 
        self.galaxy = list(
            (
            'Default values used in the analysis:',
            'nx: ' + str(cts.nx),
            'ny: ' + str(cts.ny),
            'pixel_scale: ' + str(cts.pixel_scale),
            '',
            'Galaxy drawn has the following parameters:'
            ) 
            + 
            tuple(param + ': ' + str(params[param]) for param in params.keys()) 
        )

        if fish != -1:
            self.fisher = list(
                (
                '',
                'Fisher matrix elements:',
                )
                + 
                tuple(dictToStringList(fish.fisher_matrix))
                + 
                (
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
                'Fisher analysis used the following noise bias deviation (sigma_n): ' + str(fish.sigma_n),
                '',
                'Steps used for the derivatives: '
                ) 
                + 
                tuple(param + ': ' + str(fish.steps[param]) for param in fish.steps.keys())
            )

        #add covariances in the future? (unlikely)
        if fits_biases != -1 and N != -1:
            self.fits = list(
                (
                'Number of fits: ' + str(N),
                'Biases were obtained by averaging over fits of noisy instantiations of the same image: '
                )
                +
                tuple(dictToStringList(fits_biases))
            )


            

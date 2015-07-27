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
        self.pixel_scale = .2 #arcsec/pixel 
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

class min:
    """min values for fit, may add more as needed"""
    def __init__(self, gal_image):
        cts = constants()
        self.dict = dict()
        self.dict['x0'] = - gal_image.getXMax() * cts.pixel_scale / 2
        self.dict['y0'] = - gal_image.getYMax() * cts.pixel_scale / 2
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

        #some default names for argparse.
        self.wdir = 'output'
        self.galaxy_file = 'galaxies'
        self.plots_dir = 'plots'
        self.rltsdir = 'results'
        self.info = 'info'


        self.galaxy_models = ['gaussian','exponential']
        self.psf_models = ['gaussian']

        self.gal_models_parameters = dict()
        self.gal_models_parameters['gaussian'] = ['x0','y0','flux','hlr','e1', 'e2']
        self.gal_models_parameters['exponential'] = ['x0','y0', 'flux','hlr','e1', 'e2']
        self.galaxy_parameters = list(set([elt for sublist in self.gal_models_parameters.values()for elt in sublist])) #remove repeated elements.

        self.psf_models_parameters = dict()
        self.psf_models_parameters['gaussian'] = ['psf_flux','psf_fwhm']        
        self.psf_parameters =  list(set([elt for sublist in self.psf_models_parameters.values() for elt in sublist]))

        self.parameters = self.galaxy_parameters + self.psf_parameters

        self.extra = ['id', 'model', 'psf_model'] #elements probably not in fits, more often then not, the non-numerical values written to the file. 
        
        self.fieldnames = self.extra + self.galaxy_parameters + self.psf_parameters


class info: 
    """Contains the different things that are written in the information text file"""
    
    def __init__(self, galaxies, fish = -1, fits_biases = -1, N = -1):
        cts = constants()
        params = galaxies.params 

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


            

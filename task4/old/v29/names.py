"""Define names of different parameters and models used throughout the 
package
"""

galaxy_models = ['gaussian','exponential']
psf_models = ['gaussian']

gal_models_parameters = dict()
gal_models_parameters['gaussian'] = ['x0', 'y0', 'flux', 'hlr', 'e1', 'e2']
gal_models_parameters['exponential'] = ['x0', 'y0', 'flux', 'hlr', 'e1', 'e2']
                                                    
#specify order as you want it to appear on plots/images.
gal_parameters = ['x0', 'y0', 'flux', 'hlr', 'e1', 'e2']

psf_models_parameters = dict()
psf_models_parameters['gaussian'] = ['psf_flux','psf_fwhm']
psf_parameters =  list(set(
[elt for sublist in psf_models_parameters.values() for elt in sublist]))

parameters = gal_parameters + psf_parameters
extra = ['id', 'model', 'psf_model']
fieldnames = extra + gal_parameters + psf_parameters
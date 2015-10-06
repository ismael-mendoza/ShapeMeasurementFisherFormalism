"""Define names of different parameters and models used throughout the
package
"""

extra = ['id', 'galaxy_model', 'psf_model']

gal_models_parameters = dict()
gal_models_parameters['gaussian'] = [

    'x0', 'y0',

    'flux',

    'hlr',
    'fwhm',
    'sigma',

    'e1', 'e2',
    'eta1', 'eta2',
    'e', 'beta'

]

gal_models_parameters['exponential'] = [

    'x0', 'y0',

    'flux',

    'hlr',
    'fwhm',
    'sigma',

    'e1', 'e2',
    'eta1', 'eta2'

]

gal_models_parameters['deVaucouleurs'] = [

    'x0', 'y0',

    'flux',

    'hlr',
    'fwhm',
    'sigma',

    'e1', 'e2',
    'eta1', 'eta2'

]

gal_models_parameters['bulge+disk'] = [
    'x0', 'y0',

    'flux_b', 'flux_d', 'flux_b/flux_total',

    'hlr_d', 'hlr_b', 'R_r',

    'e1', 'e2',
    'eta1', 'eta2',

    'delta_e', 'delta_theta',

    'n_d', 'n_b',
]


gal_models = gal_models_parameters.keys()
gal_parameters = list(set(
    [elt for sublist in gal_models_parameters.values() for elt in sublist]))


psf_models_parameters = dict()
psf_models_parameters['gaussian'] = [
    'psf_flux',

    'psf_fwhm',

    'psf_e1', 'psf_e2'
]

psf_models_parameters['moffat'] = ['psf_flux',

                                   'psf_fwhm',

                                   'psf_beta',

                                   'psf_e1', 'psf_e2']

psf_parameters = list(set(
    [elt for sublist in psf_models_parameters.values() for elt in sublist]))
psf_models = psf_models_parameters.keys()

# specify parameters to omit in fit and Fisher analysis in a given galaxy model
omit_fit = dict()
omit_fit['gaussian'] = extra + psf_parameters
omit_fit['exponential'] = extra + psf_parameters
omit_fit['deVaucouleurs'] = extra + psf_parameters
omit_fit['bulge+disk'] = extra + psf_parameters + ['R_r', 'delta_e',
                                                   'delta_theta',
                                                   'n_d', 'n_b', 'hlr_b']

parameters = gal_parameters + psf_parameters
fieldnames = extra + gal_parameters + psf_parameters

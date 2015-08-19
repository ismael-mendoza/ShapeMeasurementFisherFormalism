"""Define names of different parameters and models used throughout the 
package
"""

gal_models_parameters = dict()
gal_models_parameters['gaussian'] = [[('x0', 'y0')], 
                                     ['flux'], 
                                     ['hlr','fwhm','sigma'],
                                     [('e1','e2'),('eta1','eta2')]
                                     ]
gal_models_parameters['exponetial'] = [[('x0', 'y0')], 
                                      ['flux'], 
                                      ['hlr','fwhm','sigma'],
                                      [('e1','e2'),('eta1','eta2')]
                                      ]

gal_models_parameters['bulge+disk'] = [[('x0', 'y0')], 
                                      [('flux_b','flub_d')], 
                                      ['hlr_d','fwhm_d','sigma_d'],
                                      [('e1','e2'),('eta1','eta2')],
                                      [('R_r', 'delta_e','delta_theta')],
                                      [('n_d','n_b')],
                                      ]

#gal_models_parameters['sersic'] = []

galaxy_models = gal_models_parameters.keys()
                                                    

psf_models_parameters = dict()
psf_models_parameters['gaussian'] = ['flux',
                                     'fwhm',
                                     [('e1', 'e2')]
                                    ]

psf_models_parameters['moffat'] = [['flux'], 
                                   ['fwhm'], 
                                   ['beta'], 
                                   [('e1', 'e2')]
                                  ]
#add psf_ preffix to all
# psf_parameters = list(set(
# [elt for sublist in psf_models_parameters.values() for elt in sublist]))
# psf_models = psf_models_parameters.keys()

#all_parameters = gal_parameters + psf_parameters
extra = ['id', 'model', 'psf_model']




#specify parameters to omit in fit and Fisher analysis in a given galaxy model
omit_fit = dict()
omit_fit['gaussian'] = extra + psf_parameters
omit_fit['exponential'] = extra + psf_parameters
omit_fit['bulge+disk'] = extra + psf_parameters + [] #

fieldnames = extra + gal_parameters + psf_parameters
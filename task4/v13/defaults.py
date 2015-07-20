#!/usr/bin/env python
"""Some of the default values that are used in the overall.
In general defaults should not be able to be change by the user. But some exceptions for steps or max & min may be implemented in the future.

"""

class constants:
	"""General global constants"""
	def __init__(self):
		self.nx = 48
		self.ny = 48
		self.pixel_scale = .2

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
		self.galaxy_parameters['exponential'] = ['x0','y0','hlr','e1', 'e2']

		self.parameters = [] 
		for lst in self.galaxy_parameters.values() + self.psf_parameters.values():
			for element in lst:
				self.parameters.append(element)
		self.parameters = list(set(self.parameters)) #remove repeated elements.

		self.fieldnames = ['id', 'model', 'psf_model'] + self.parameters



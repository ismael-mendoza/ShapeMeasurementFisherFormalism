#a submodule with some defaults the program uses.

class constants:
	"""General global constants"""
	def __init__(self):
		self.nx = 48
		self.ny = 48
		self.pixel_scale = .2

class parameters():
	def __init__(self):
		self.dflt_params = dict()
		self.dflt_params['x0'] = 0. 
		self.dflt_params['y0'] = 0. 
		self.dflt_params['flux'] = 100. 
		self.dflt_params['hlr'] = .25
		self.dflt_params['e1'] = 0. 
		self.dflt_params['e2'] = 0. 
		self.dflt_params['psf_flux'] = 1. 
		self.dflt_params['psf_fwhm'] = .7 

class min():
	"""min values for fit, may add more as needed"""
	def __init__(self, gal_image):
		self.min_values = dict()
		self.min_values['x0'] = gal_image.getXMin()
		self.min_values['y0'] = gal_image.getYMin() 
		self.min_values['flux'] = 0. 
		self.min_values['hlr'] = 0.
		self.min_values['e1'] = -.7 
		self.min_values['e2'] = -.7

class max():
	"""max values for fit, may add more as needed"""
	def __init__(self, gal_image):
		self.max_values = dict()
		self.max_values['x0'] = gal_image.getXMax() 
		self.max_values['y0'] = gal_image.getYMax()
		self.max_values['flux'] = None 
		self.max_values['hlr'] = .25
		self.max_values['e1'] = .7 
		self.max_values['e2'] = .7

class names:
	"""Contains the different names for parameters, headers, etc. in lists"""
	def __init__(self):
		self.galaxy_models = ['gaussian']
		self.psf_models = ['gaussian']
		self.galaxy_parameters = dict()
		self.psf_parameters = dict()
		self.galaxy_parameters['gaussian'] = ['x0','y0','flux','hlr','e1', 'e2']
		self.psf_parameters['gaussian'] = ['psf_flux','psf_fwhm']
		self.parameters = [] #also remove repeated elements when necessary.
		for lst in self.galaxy_parameters.values() + self.psf_parameters.values():
			for element in lst:
				self.parameters.append(element)
		self.header = ['id', 'model', 'psf_model'] + self.parameters #later use list(set(x)) to remove repeated elements.

class steps: 
	"""Define the steps for derivatives of each individual parameter."""
	def __init__(self, params):
		self.dct = dict()
		self.dct['flux'] = params['flux'] * .01
		self.dct['hlr'] = params['hlr'] *.01
		self.dct['e1'] = .01
		self.dct['e2'] = .01
		self.dct['x0'] = .01
		self.dct['y0'] = .01
	    #steps['q'] = orig_params['q']*.01
	    #steps['beta'] = orig_params['beta']*.01

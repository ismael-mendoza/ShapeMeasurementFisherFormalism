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


class names:
	"""Contains the different names for parameters, headers, etc. in lists"""
	def __init__(self):
		self.galaxy_choices = ['gaussian']
		self.gauss = ['x0','y0','flux','hlr','e1', 'e2']
		self.psf_gauss = ['psf_flux','psf_fwhm']
		self.parameters = self.gauss + self.psf_gauss #keep adding new ones if necessary. 
		self.header = ['id', 'type'] + self.gauss + self.psf_gauss #later use list(set(x)) to remove repeated elements.



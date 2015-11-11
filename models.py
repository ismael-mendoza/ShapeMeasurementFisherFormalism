###Todo-
#how to initialize a subclass in python, to not copy paste self.__init__ in all models??
#

import galsim

import sys 


#referencing own module
curr_module = sys.modules[__name__] 

def getExtra():
    return ['id', 'galaxy_model', 'psf_model']

#make sure names of model class is the same name as the one to generate.
#write down this class all galaxy modules desired.
class model(object):

    def __init__(self, params=None):
        self.parameters = []

        #specify any special parameters to ignore in the fitting.
        self.omit_specific = []

        if(params):
            self.omit_fit = self.getOmitFit()
            self.gal = self.getGal(params)

    def getOmitFit(self):
        return getExtra() + getPsfParameters() + self.omit_specific

    def getGal(self, params):
        gal = self.getProfile(params)
        gal = self.shear(gal, params)
        gal = self.shift(gal, params)
        return gal

    def getProfile(self, params):
        pass

    def shear(self, gal, params):
        if 'e1' in params and 'e2' in params:
            return gal.shear(e1=params['e1'], e2=params['e2'])
        elif 'eta1' in params and 'eta2' in params:
            return gal.shear(eta1=params['eta1'], eta2=params['eta2'])
        elif 'q' in params and 'beta' in params:
            return gal.shear(q=params['q'], beta=(params['beta']*
                                                  galsim.radians))
        elif 'e' in params and 'beta' in params:
            return gal.shear(e=params['e'], beta=(params['beta'] *
                                                  galsim.radians))
        else:
            raise ValueError('The shear for the galaxy was not specified.')

    def shift(self, gal, params):
        if 'x0' and 'y0' in params:
            return gal.shift(params['x0'], params['y0'])
        else:
            raise ValueError('The shift for the galaxy was not specified.')


class gaussian(model):

    def __init__(self, params=None):

        self.parameters = [
            'x0', 'y0',

            'flux',

            'hlr',
            'fwhm',
            'sigma',

            'e1', 'e2',
            'eta1', 'eta2',
            'e', 'beta'
        ]
        self.omit_specific = []

        if(params):
            self.omit_fit = self.getOmitFit()
            self.gal = self.getGal(params)

    def getProfile(self, params):
        if 'hlr' in params:
            return galsim.Gaussian(flux=params['flux'],
                                   half_light_radius=params['hlr'])
        elif 'sigma' in params:
            return galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

        else:
            raise ValueError('Did not specify the size.')



class exponential(model):

    def __init__(self, params=None):
        self.parameters = [
            'x0', 'y0',

            'flux',

            'hlr',
            'fwhm',
            'sigma',

            'e1', 'e2',
            'eta1', 'eta2'
        ]
        self.omit_specific = []
        if(params):
            self.omit_fit = self.getOmitFit()
            self.gal = self.getGal(params)

    def getProfile(self, params):
        return galsim.Exponential(flux=params['flux'],
                                  half_light_radius=params['hlr'])


class deVaucouleurs(model):

    def __init__(self, params=None):
        self.parameters = [
            'x0', 'y0',

            'flux',

            'hlr',
            'fwhm',
            'sigma',

            'e1', 'e2',
            'eta1', 'eta2'
        ]
        self.omit_specific = []

        if(params):
            self.omit_fit = self.getOmitFit()
            self.gal = self.getGal(params)

    def getProfile(self, params):
        return galsim.DeVaucouleurs(half_light_radius=params['hlr'],
                                    flux=params['flux'])


class bulgeDisk(model):

    def __init__(self,params=None):
        self.parameters = [
        'x0', 'y0',

        'flux_b', 'flux_d', 'flux_b/flux_total',

        'hlr_d', 'hlr_b', 'R_r',

        'e1', 'e2',
        'eta1', 'eta2',

        'delta_e', 'delta_theta',

        'n_d', 'n_b'
        ]
        self.omit_specific = []

        if(params):
            self.omit_fit = self.getOmitFit()
            self.gal = self.getGal(params)

    def getProfile(self, params):
        if 'flux_b' and 'flux_d' in params:
            flux_b = params['flux_b']
            flux_d = params['flux_d']

        elif 'flux_b' and 'flux_b/flux_total' in params:
            raise NotImplementedError('Need to implement flux_b/flux_total')

        if 'hlr_b' and 'hlr_d' in params:
            hlr_d = params['hlr_d']
            hlr_b = params['hlr_b']

        if 'hlr_d' and 'R_r' in params:
            hlr_d = params['hlr_d']
            hlr_b = params['R_r'] * hlr_d

        bulge = galsim.Sersic(n=params['n_b'], half_light_radius=hlr_b,
                              flux=flux_b)

        disk = galsim.Sersic(n=params['n_d'], half_light_radius=hlr_d,
                             flux=flux_d)

        if 'delta_e' in params or 'delta_theta' in params:
            raise NotImplementedError('Need to implement delta_e or '
                                      'delta_theta.')

        return (bulge+disk)

class psf_model(object):

    def __init__(self, params=None):
        self.parameters = []
        if(params):
            self.psf = self.getProfile(params) #ignore shear for now.

    def getProfile(self, params):
        pass

    def shearPsf(psf, params):
        return  psf.shear(e1=params.get('psf_e1', 0),
                          e2=params.get('psf_e2', 0))



class psf_gaussian(psf_model):
    def __init__(self, params=None):
        self.parameters = [
            'psf_flux',

            'psf_fwhm',

            'psf_e1', 'psf_e2'
        ]
        if(params):
            self.psf = self.getProfile(params)

    def getProfile(self, params):
        if 'psf_hlr' in params:
            psf = galsim.Gaussian(
                flux=params['psf_flux'],
                half_light_radius=params['psf_hlr'])

        elif 'psf_sigma' in params:
            psf = galsim.Gaussian(flux=params['psf_flux'],
                                  sigma=params['psf_sigma'])
        elif 'psf_fwhm' in params:
            psf = galsim.Gaussian(flux=params['psf_flux'],
                                  fwhm=params['psf_fwhm'])
        else:
            raise ValueError('Size of PSF was not specified.')
        return psf

class psf_moffat(psf_model):

    def __init__(self, params=None):
        self.parameters = [
            'psf_flux',

           'psf_fwhm',

           'psf_beta',

           'psf_e1', 'psf_e2'
        ]
        if(params):
            self.psf = self.getProfile(params)

    def getProfile(self, params):
        return galsim.Moffat(beta=params['psf_beta'],
                            fwhm=params['psf_fwhm'],
                            flux=params['psf_flux'])

# def psfProfile(params):
#     # each psf model.
#     if params['psf_model'] == 'gaussian':
#         return gaussianPsfProfile(params)

#     elif params['psf_model'] == 'moffat':
#         return moffatPsfProfile(params)

#     else:
#         raise ValueError('PSF model was not specified.')


# def get_model(params):

#     if params['galaxy_model'] == 'gaussian':
#         return gaussianProfile(params)


#     elif params['galaxy_model'] == 'exp2onential':
#         return exponentialProfile(params)


#     elif params['galaxy_model'] == 'deVaucouleurs':
#         return deVaucouleursProfile(params)


#     elif params['galaxy_model'] == 'bulge+disk':
#         return bulgeDiskProfile(params)

#     else:
#         raise ValueError('The galaxy model was not specified.')


#iterate over all subclasses to get fieldnames.
#instantiate each subclass and get self.parameters from each and add to
#fieldnames.

def getGalParameters():
    gal_parameters = []
    subclasses = [cls for cls in vars(curr_module)['model'].__subclasses__()]
    for cls in subclasses:
        obj = cls()
        gal_parameters += obj.parameters

    #remove duplicates
    gal_parameters = list(set(gal_parameters))

    return gal_parameters


def getPsfParameters():
    psf_parameters = []
    subclasses = [cls for cls in vars(curr_module)['psf_model'].__subclasses__()]
    for cls in subclasses:
        obj = cls()
        psf_parameters += obj.parameters

    psf_parameters = list(set(psf_parameters))
        
    return psf_parameters

def getAllParameters():
    return getGalParameters() + getPsfParameters()

def getFieldnames():
    return getExtra() + getGalParameters() + getPsfParameters()


#return the corresponding class to the model specified in params
def getModelCls(params):
    subclasses = [cls for cls in vars(curr_module)['model'].__subclasses__()]
    for cls in subclasses:
        if(cls.__name__ == params['galaxy_model']):
            return cls
    raise NotImplementedError('Have not implemented that galaxy model')
#make dictionary of classes to names. no need.

#print([cls.__name__ for cls in vars()['Foo'].__subclasses__()])
#print([cls for cls in vars()['Foo'].__subclasses__()])

#return the corresponding psf class specified in params.
def getPsfModelCls(params):
    subclasses = [cls for cls in vars(curr_module)['psf_model'].__subclasses__()]
    for cls in subclasses:
        if(cls.__name__ == params['psf_model']):
            return cls
    raise NotImplementedError('Have not implemented that psf model')


#just to display choices in generate.py
def getAllModels():
    gal_models = []
    subclasses = [cls for cls in vars(curr_module)['model'].__subclasses__()]
    for cls in subclasses:
        gal_models.append(cls.__name__)

    return gal_models

def getAllPsfModels():
    psf_models = []
    subclasses = [cls for cls in vars(curr_module)['psf_model'].__subclasses__()]
    for cls in subclasses:
        psf_models.append(cls.__name__)

    return psf_models



"""
This module contains the models that are used for galaxies and psfs for galsim, more 
information of how to add your own models for both galaxies and psfs can be found in tutorial
"""
import galsim
import sys 


#referencing own module
curr_module = sys.modules[__name__] 

def getExtra():
    return ['id', 'galaxy_model', 'psf_model']

#make sure names of model class is the same name as the one to generate.
#write down this class in all galaxy modules.
class model(object):
    
    parameters = []
    omit_general = [] #omit always for all instances. 

    def __init__(self, params=None, params_omit=None):

        self.omit_fit = self.getOmitFit()
        
        if(params_omit):
            self.setOmitSpecific(params_omit)
        if(params):
            self.gal = self.getGal(params)

    def getOmitFit(self):
        #remove redundancy. 
        omit = list(set(getExtra() + getPsfParameters() + self.omit_general)) 
        return omit

    #pass in a list if you want to omit specific parameters 
    #for that instance. 
    def setOmitSpecific(self, params_omit):
        self.omit_fit = list(set(params_omit + self.omit_fit))


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

        elif 'g1' in params and 'g2' in params:
            return gal.shear(g1=params['g1'],g2=params['g2'])
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

    parameters = [
        'x0', 'y0',

        'flux',

        'hlr',
        'fwhm',
        'sigma',

        'e1', 'e2',
        'eta1', 'eta2',
        'e', 'q', 'beta',
    ]
    omit_general = []

    def __init__(self, params=None, params_omit=None):
        model.__init__(self, params, params_omit)

    def getProfile(self, params):
        if 'flux' not in params:
            raise ValueError('Did not specify flux')
            
        if 'hlr' in params:
            return galsim.Gaussian(flux=params['flux'],
                                   half_light_radius=params['hlr'])
        elif 'sigma' in params:
            return galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

        else:
            raise ValueError('Did not specify the size.')


class exponential(model):

    parameters = [
        'x0', 'y0',

        'flux',

        'hlr',
        'fwhm',
        'sigma',

        'e1', 'e2', 
        'g1', 'g2',
        'eta1', 'eta2',
        'q','beta'
    ]
    omit_general = []
    
    def __init__(self, params=None, params_omit=None):
        model.__init__(self, params, params_omit)


    def getProfile(self, params):
        if 'flux' not in params:
            raise ValueError('Did not specify flux')

        if 'hlr' in params:
            return galsim.Exponential(flux=params['flux'],
                                      half_light_radius=params['hlr'])

        else:
            raise ValueError('Did not specify the size.')


class deVaucouleurs(model):

    parameters = [
        'x0', 'y0',

        'flux',

        'hlr',
        'fwhm',
        'sigma',

        'e1', 'e2',
        'eta1', 'eta2'
    ]

    omit_general = []


    def __init__(self, params=None, params_omit=None):
        model.__init__(self, params, params_omit)

    def getProfile(self, params):
        return galsim.DeVaucouleurs(half_light_radius=params['hlr'],
                                    flux=params['flux'])


class bulgeDisk(model):

    parameters = [
    'x0', 'y0',

    'flux_b', 'flux_d', 'flux_b/flux_total',

    'hlr_d', 'hlr_b', 'R_r',

    'e1', 'e2',
    'eta1', 'eta2',

    'delta_e', 'delta_theta',

    'n_d', 'n_b'
    ]

    omit_general = ['delta_e', 'delta_theta', 'n_d', 'n_b']

    def __init__(self, params=None, params_omit=None):
        model.__init__(self, params, params_omit)

    def getProfile(self, params):
        if 'flux_b' in params and 'flux_d' in params:
            flux_b = params['flux_b']
            flux_d = params['flux_d']

        elif 'flux_b' and 'flux_b/flux_total' in params:
            raise NotImplementedError('Need to implement flux_b/flux_total')

        else:
            raise ValueError('Flux was not specified.')

        if 'hlr_b' in params and 'hlr_d' in params:
            hlr_d = params['hlr_d']
            hlr_b = params['hlr_b']

        elif 'hlr_d' and 'R_r' in params:
            hlr_d = params['hlr_d']
            hlr_b = params['R_r'] * hlr_d

        else: 
            raise ValueError('Size was not specified.')

        bulge = galsim.Sersic(n=params['n_b'], half_light_radius=hlr_b,
                              flux=flux_b)

        disk = galsim.Sersic(n=params['n_d'], half_light_radius=hlr_d,
                             flux=flux_d)

        if 'delta_e' in params or 'delta_theta' in params:
            raise NotImplementedError('Need to implement delta_e or '
                                      'delta_theta.')

        return (bulge+disk)

class bulgeDisk6(model):

    parameters = [
    'x0', 'y0',

    'flux',

    'hlr',

    'e1', 'e2',
    'eta1', 'eta2',

    'n_d', 'n_b'
    ]

    omit_general = ['n_d', 'n_b']

    def __init__(self, params=None, params_omit=None):
        model.__init__(self, params, params_omit)

    def getProfile(self, params):
        if 'flux' in params:
            flux = params['flux']

        else:
            raise ValueError('Flux was not specified.')

        if 'hlr' in params:
            hlr = params['hlr']

        else: 
            raise ValueError('Size was not specified.')

        bulge = galsim.Sersic(n=params['n_b'], half_light_radius=hlr,
                              flux=flux)

        disk = galsim.Sersic(n=params['n_d'], half_light_radius=hlr,
                             flux=flux)

        return (bulge+disk)

class psf_model(object):

    parameters = []

    def __init__(self, params=None):
        if(params):
            self.psf = self.getProfile(params) #ignore shear for now.

    def getProfile(self, params):
        pass

    def shearPsf(psf, params):
        return  psf.shear(e1=params.get('psf_e1', 0),
                          e2=params.get('psf_e2', 0))



class psf_gaussian(psf_model):
    
    parameters = [
        'psf_flux',

        'psf_fwhm',

        'psf_e1', 'psf_e2'
    ]

    def __init__(self, params=None):
         psf_model.__init__(self, params)

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

    parameters = [
       'psf_flux',

       'psf_fwhm',
       'psf_hlr',

       'psf_beta',

       'psf_e1', 'psf_e2'
    ]

    def __init__(self, params=None):
        psf_model.__init__(self, params)

    def getProfile(self, params):
        if 'psf_fwhm' in params:
            return galsim.Moffat(beta=params['psf_beta'],
                                 fwhm=params['psf_fwhm'],
                                 flux=params['psf_flux'])
        elif 'psf_hlr' in params:
            return galsim.Moffat(beta=params['psf_beta'],
                                 half_light_radius=params['psf_hlr'],
                                 flux=params['psf_flux'])


#iterate over all subclasses to get fieldnames.
#instantiate each subclass and get self.parameters from each and add to
#fieldnames.
def getGalParameters():
    gal_parameters = []
    subclasses = [cls for cls in vars(curr_module)['model'].__subclasses__()]
    for cls in subclasses:
        gal_parameters += cls.parameters

    #remove duplicates
    gal_parameters = list(set(gal_parameters))

    return gal_parameters


def getPsfParameters():
    psf_parameters = []
    subclasses = [cls for cls in vars(curr_module)['psf_model'].__subclasses__()]
    for cls in subclasses:
        psf_parameters += cls.parameters

    psf_parameters = list(set(psf_parameters))
        
    return psf_parameters

def getAllParameters():
    return getGalParameters() + getPsfParameters()

def getFieldnames():
    return getExtra() + getGalParameters() + getPsfParameters()


#return the corresponding class to the model specified in params
def getModelCls(model):
    subclasses = [cls for cls in vars(curr_module)['model'].__subclasses__()]
    for cls in subclasses:
        if(cls.__name__ == model):
            return cls
    raise NotImplementedError('Have not implemented that galaxy model')


#return the corresponding psf class specified in params.
def getPsfModelCls(model):
    subclasses = [cls for cls in vars(curr_module)['psf_model'].__subclasses__()]
    for cls in subclasses:
        if(cls.__name__ == model):
            return cls
    raise NotImplementedError('Have not implemented that psf model')


#used to display choices in generate.py
def getAllModels():
    gal_models = []
    subclasses = [cls for cls in vars(curr_module)['model'].__subclasses__()]
    for cls in subclasses:
        gal_models.append(cls.__name__)

    return gal_models

#used to display choices in generate.py
def getAllPsfModels():
    psf_models = []
    subclasses = [cls for cls in vars(curr_module)['psf_model'].__subclasses__()]
    for cls in subclasses:
        psf_models.append(cls.__name__)

    return psf_models



"""
This module contains the models that are used for galaxies and psfs for galsim, more 
information of how to add your own models for both galaxies and psfs can be found in the corresponding tutorial. 
"""
import inspect
import sys

import galsim

# referencing itself.
curr_module = sys.modules[__name__]


def get_extra():
    return ['id', 'galaxy_model', 'psf_model']


# make sure names of model class is the same name as the one to generate.
class Model(object):
    parameters = []
    omit_general = []  # omit always for all instances of this class.

    def __init__(self, params=None, params_omit=None):

        self.omit_fit = self.get_omit_fit()

        if params_omit:
            self.set_omit_specific(params_omit)
        if params:
            self.gal = self.get_gal(params)

    def get_omit_fit(self):
        # remove redundancy.
        omit = list(set(get_extra() + get_psf_parameters() + self.omit_general))
        return omit

    # pass in a list if you want to omit specific parameters
    # for that instance.
    def set_omit_specific(self, params_omit):
        self.omit_fit = list(set(params_omit + self.omit_fit))

    def get_gal(self, params):
        gal = self.get_profile(params)
        gal = self.shear(gal, params)
        gal = self.shift(gal, params)
        return gal

    def get_profile(self, params):
        pass

    def shear(self, gal, params):
        if 'e1' in params and 'e2' in params:
            return gal.shear(e1=params['e1'], e2=params['e2'])

        elif 'g1' in params and 'g2' in params:
            return gal.shear(g1=params['g1'], g2=params['g2'])

        elif 'eta1' in params and 'eta2' in params:
            return gal.shear(eta1=params['eta1'], eta2=params['eta2'])

        elif 'q' in params and 'beta' in params:
            return gal.shear(q=params['q'], beta=(params['beta'] *
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


class Gaussian(Model):
    parameters = [
        'flux',

        'x0', 'y0',

        'hlr',
        'fwhm',
        'sigma',

        'e1', 'e2',
        'eta1', 'eta2',
        'e', 'q', 'beta',
        'g1', 'g2'
    ]
    omit_general = []

    def __init__(self, params=None, params_omit=None):
        Model.__init__(self, params, params_omit)

    def get_profile(self, params):
        if 'flux' not in params:
            raise ValueError('Did not specify flux')

        if 'hlr' in params:
            return galsim.Gaussian(flux=params['flux'],
                                   half_light_radius=params['hlr'])
        elif 'sigma' in params:
            return galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

        else:
            raise ValueError('Did not specify the size.')


class Exponential(Model):
    parameters = [
        'x0', 'y0',

        'flux',

        'hlr',
        'fwhm',
        'sigma',

        'e1', 'e2',
        'g1', 'g2',
        'eta1', 'eta2',
        'q', 'beta'
    ]
    omit_general = []

    def __init__(self, params=None, params_omit=None):
        Model.__init__(self, params, params_omit)

    def get_profile(self, params):
        if 'flux' not in params:
            raise ValueError('Did not specify flux')

        if 'hlr' in params:
            return galsim.Exponential(flux=params['flux'],
                                      half_light_radius=params['hlr'])

        else:
            raise ValueError('Did not specify the size.')


class BulgeDisk(Model):
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
        Model.__init__(self, params, params_omit)

    def get_profile(self, params):
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

        return bulge + disk


class BulgeDisk6(Model):
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
        Model.__init__(self, params, params_omit)

    def get_profile(self, params):
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

        return bulge + disk


class PsfModel(object):
    parameters = []

    def __init__(self, params=None):
        if params:
            self.psf = self.get_profile(params)  # ignore shear for now.

    def get_profile(self, params):
        pass

    def shear_psf(self, params):
        return self.psf.shear(e1=params.get('psf_e1', 0),
                              e2=params.get('psf_e2', 0))


class GaussianPsf(PsfModel):
    parameters = [
        'psf_flux',

        'psf_fwhm',

        'psf_e1', 'psf_e2'
    ]

    def __init__(self, params=None):
        PsfModel.__init__(self, params)

    def get_profile(self, params):
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


class MoffatPsf(PsfModel):
    parameters = [
        'psf_flux',

        'psf_fwhm',
        'psf_hlr',

        'psf_beta',

        'psf_e1', 'psf_e2'
    ]

    def __init__(self, params=None):
        PsfModel.__init__(self, params)

    def get_profile(self, params):
        if 'psf_fwhm' in params:
            return galsim.Moffat(beta=params['psf_beta'],
                                 fwhm=params['psf_fwhm'],
                                 flux=params['psf_flux'])
        elif 'psf_hlr' in params:
            return galsim.Moffat(beta=params['psf_beta'],
                                 half_light_radius=params['psf_hlr'],
                                 flux=params['psf_flux'])


def get_gal_parameters():
    """
    Iterate over all subclasses to get fieldnames. Instantiate each subclass and get self.parameters
    from each and add to fieldnames.
    Returns:

    """
    gal_parameters = []
    classes = [obj for name, obj in inspect.getmembers(sys.modules[__name__]) if inspect.isclass(obj)]
    for cls in classes:
        if 'Psf' not in cls.__name__:
            gal_parameters += cls.parameters

    # remove duplicates
    gal_parameters = list(set(gal_parameters))

    return gal_parameters


def get_psf_parameters():
    psf_parameters = []
    classes = [obj for name, obj in inspect.getmembers(sys.modules[__name__]) if inspect.isclass(obj)]
    for cls in classes:
        if 'Psf' in cls.__name__ and cls.__name__ != 'PsfModel':
            psf_parameters += cls.parameters

    psf_parameters = list(set(psf_parameters))

    return psf_parameters


def get_all_parameters():
    return get_gal_parameters() + get_psf_parameters()


def get_fieldnames():
    return get_extra() + get_gal_parameters() + get_psf_parameters()


def get_model_cls(model):
    """Return the corresponding class to the model specified in params"""
    classes = [obj for name, obj in inspect.getmembers(sys.modules[__name__]) if inspect.isclass(obj)]
    for cls in classes:
        if cls.__name__.lower() == model:
            return cls
    raise NotImplementedError('Have not implemented that galaxy model')


def get_all_models():
    """Used to display choices in generate.py"""
    classes = [obj for name, obj in inspect.getmembers(sys.modules[__name__]) if inspect.isclass(obj)]
    gal_models = []
    for cls in classes:
        if 'Psf' not in cls.__name__ and 'Model' != cls.__name__:
            gal_models.append(cls.__name__.lower())
    return gal_models


def get_all_psf_models():
    """Used to display choices in generate.py"""
    classes = [obj for name, obj in inspect.getmembers(sys.modules[__name__]) if inspect.isclass(obj)]
    psf_models = []
    for cls in classes:
        if 'Psf' in cls.__name__ and cls.__name__ != 'PsfModel':
            psf_models.append(cls.__name__.lower())

    return psf_models

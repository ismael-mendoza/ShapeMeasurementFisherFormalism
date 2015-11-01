import galsim


extra = ['id', 'galaxy_model', 'psf_model']



###doubts,
#how to initialize a subclass in python??


def shearPsf(psf, params):
    return  psf.shear(e1=params.get('psf_e1', 0), e2=params.get('psf_e2', 0))


class model(object):

    def __init__(self):
        self.parameters = []
        self.omit_fit = []
        self.profile = get_profile(params)

    def get_profile(self, params):
        pass

    def shear(self, params):
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

    def shift(self, params):
        if 'x0' and 'y0' in params:
            gal = gal.shift(params['x0'], params['y0'])
        else:
            raise ValueError('The shift for the galaxy was not specified.')


class gaussian(model):

    def __init__(self):

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

        self.omit_fit = []


    def get_profile(self, params):
        if 'hlr' in params:
            return galsim.Gaussian(flux=params['flux'],
                                   half_light_radius=params['hlr'])
        elif 'sigma' in params:
            return galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

        else:
            raise ValueError('Did not specify the size.')



class exponential(model):

    def __init__(self):
        self.parameters = []
        self.omit_fit = []
        self.profile = get_profile(params)

    def get_profile(self, params):
    return galsim.Exponential(flux=params['flux'],
                              half_light_radius=params['hlr'])


class deVaucouleurs(model):

    def __init__(self):
        self.parameters = []
        self.omit_fit = []
        self.profile = get_profile(params)

    def get_profile(self, params):
        return galsim.DeVaucouleurs(half_light_radius=params['hlr'],
                                    flux=params['flux'])

class bulgeDisk(model):

    def __init__(self):
        self.parameters = []
        self.omit_fit = []
        self.profile = get_profile(params)

    def get_profile(self, params):
        pass



def gaussianProfile(params):



def exponentialProfile(params):

def deVaucouleursProfile(params):


def bulgeDiskProfile(params):
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

def psfProfile(params):
    # each psf model.
    if params['psf_model'] == 'gaussian':
        return gaussianPsfProfile(params)

    elif params['psf_model'] == 'moffat':
        return moffatPsfProfile(params)

    else:
        raise ValueError('PSF model was not specified.')

def gaussianPsfProfile(params):
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

def moffatPsfProfile(params):
    psf = galsim.Moffat(beta=params['psf_beta'],
                        fwhm=params['psf_fwhm'],
                        flux=params['psf_flux'])



def get_model(params):

    if params['galaxy_model'] == 'gaussian':
        return gaussianProfile(params)


    elif params['galaxy_model'] == 'exponential':
        return exponentialProfile(params)


    elif params['galaxy_model'] == 'deVaucouleurs':
        return deVaucouleursProfile(params)


    elif params['galaxy_model'] == 'bulge+disk':
        return bulgeDiskProfile(params)

    else:
        raise ValueError('The galaxy model was not specified.')


def drawGalaxy(params):
    """Return the image of a single galaxy optionally drawn with a psf.

    Look at the :mod:`names.py` to figure out which galaxy models and psf
    models are supported as well as their corresponding implemented parameters.
    This function uses galsim extensively to draw the different models.

    Args:
        params(dict): Dictionary containing the information of a single
        galaxy where the keys is the name(str) of the parameter and the
        values are the values of the parametes.

    Returns:
        A galsim.Image object.
    """

    gal = profile(params)

    gal = shear(gal, params)

    gal = shift(gal, params)

    final = gal

    # generate psf.
    if params.get('psf_flux', 0) != 0:

        if params.get('psf_flux', 1) != 1:
            raise ValueError('I do not think you want a psf of flux not 1')

        psf = psfProfile(params)

        psf = shear_psf(psf, params)

        final = galsim.Convolve([gal, psf])

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    return final.drawImage(scale=defaults.PIXEL_SCALE, nx=defaults.NX,
                           ny=defaults.NY)

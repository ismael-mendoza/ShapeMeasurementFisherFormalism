#!/usr/bin/env python

import os

import csv

import defaults

import galsim

import copy

import names

import math


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

    ####
    #cretea a class that only draws galaxies and has different methods for each of the following steps: model selection, setting shear, shift, psf.
    ####
    # check each model of galaxy
    if params['galaxy_model'] == 'gaussian':
        if 'hlr' in params:
            gal = galsim.Gaussian(
                flux=params['flux'], half_light_radius=params['hlr'])
        elif 'sigma' in params:
            gal = galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

    elif params['galaxy_model'] == 'exponential':
        gal = galsim.Exponential(
            flux=params['flux'], half_light_radius=params['hlr'])

    elif params['galaxy_model'] == 'deVaucouleurs':
        gal = galsim.DeVaucouleurs(half_light_radius=params['hlr'],
                                   flux=params['flux'])

    elif params['galaxy_model'] == 'bulge+disk':
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

        gal = bulge + disk
    else:
        raise ValueError('The galaxy model was not specified.')

    # set shear.
    if 'e1' in params and 'e2' in params:
        gal = gal.shear(e1=params['e1'], e2=params['e2'])
    elif 'eta1' in params and 'eta2' in params:
        gal = gal.shear(eta1=params['eta1'], eta2=params['eta2'])
    elif 'q'in params and 'beta' in params:
        gal = gal.shear(q=params['q'], beta=params['beta'] * galsim.radians)
    elif 'e' in params and 'beta' in params:
        gal = gal.shear(e=params['e'], beta=params['beta'] * galsim.radians)
    else:
        raise ValueError('The shear for the galaxy was not specified.')

    # shift galaxy centroid.
    if 'x0' and 'y0' in params:
        gal = gal.shift(params['x0'], params['y0'])
    else:
        raise ValueError('The shift for the galaxy was not specified.')

    # generate psf.
    if params.get('psf_flux', 0) != 0:

        if params.get('psf_flux', 1) != 1:
            raise ValueError('I do not think you want a psf of flux not'
                             '1')

        # each psf model.
        if params['psf_model'] == 'gaussian':
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

        elif params['psf_model'] == 'moffat':
            psf = galsim.Moffat(beta=params['psf_beta'],
                                fwhm=params['psf_fwhm'],
                                flux=params['psf_flux'])

        else:
            raise ValueError('PSF model was not specified.')

        psf = psf.shear(e1=params.get('psf_e1', 0), e2=params.get('psf_e2', 0))

        final = galsim.Convolve([gal, psf])
    else:
        final = gal

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    return final.drawImage(scale=defaults.PIXEL_SCALE, nx=defaults.NX,
                           ny=defaults.NY)



def drawGalaxies(fit_params=None, id_params=None, image=False,
                 g_parameters=None, **kwargs):
    """Return the image of a set of galaxies.

    One of the the following must be specified:
        fit_params (and nfit_params as **kwargs. Used specially in
        :mod:`runfits.py`).
        id_params
        g_parameters (from which id_params is extracted.)
    This function draws each of the galaxies specified in id_params and then
    sums them together to get a final galaxy.

    Args:
        fit_params(dict): Partial form of id_params that only includes the
                          parameters to be used for the fit.
                          For details, :class:`GParameters`
        id_params(dict): Dictionary containing each of the galaxies
                         parameters. For details, :class:`GParameters`
        g_parameters(:class:`GParameters`): An object containing different
                                            forms of the galaxy parameters.
        image(bool): If :bool:True returns an galsim.Image otherwise it returns
                     a np.array

    Returns:
        A galsim.Image or a np.array
    """

    if id_params is None and g_parameters is None:
        fit_params.update(kwargs)
        id_params = GParameters.convertParams_Id(fit_params)

    if g_parameters is not None:
        id_params = g_parameters.id_params

    gals = []

    for gal_id in id_params:
        gals.append(drawGalaxy(id_params[gal_id]))
    final = sum(gals)

    if image is False:
        return final.array
    else:
        return final


def addNoise(image, snr, noise_seed=0):
    """Set gaussian noise to the given galsim.Image.

    Args:
        image(galsim.Image): Galaxy image that noise is going to be added
                             to.
        snr(float): Signal to noise ratio.
        noise_seed(int): Seed to set to galsim.BaseDeviate which
                         will create the galsim.noise instance.

    Returns:
        A galsim.Image, variance_noise tuple. The image is the noisy version
        of the original image and variance_noise is the noise variance on each
        pixel due to the added noise.
    """

    noisy_image = copy.deepcopy(image)  # do not alter original image.
    bd = galsim.BaseDeviate(noise_seed)
    noise = galsim.GaussianNoise(rng=bd)
    variance_noise = noisy_image.addNoiseSNR(noise, snr,
                                             preserve_flux=True)
    return noisy_image, variance_noise


class GParameters(object):

    """Class that manages galaxies parameters obtained from galaxies.csv

    This class reads a galaxies.csv file located in the specified project
    directory and extracts the parameters of each galaxy contained in it.

        Args:
            project(str): String point to the directory specified by the user.
            id_params(dict): See Attributes for details. Makes it possible to
                             create a :class:`GParameters` object without
                             galaxies.csv file.

        Attributes:
            omit_fit(dict): Dictionary defined in :mod:`names.py` containing
                            the parameters that should not be included in the
                            analysis for a particular galaxy model.
            id_params(dict): Dictionary whose keys are the ids of each of the
                             galaxies specified in galaxies.csv, and that map
                             to another dictionary that can be taken in by
                             :func:`drawGalaxy`
            params(dict): Dictionary that encodes the same information as
                          id_params but in a different form. Combines each of
                          the dictionaries contained in id_params into a
                          single dictionary that contains all parameters but in
                          the form param_#.
            fit_params(dict): Dictionary similar to the params attribute but
                              without the parameters specified in omit_fit.
            nfit_params(dict): Dictionary that contains all the parameters not
                               contained in fit_params. Usually used for
                               **kwargs in conjunction with fit_params to draw
                               Galaxies.
            ordered_fit_names(list): A list containing the keys of fit_params
                                     in a desirable order.
            num_galaxies(int): Number of galaxies specified.
    """

    def __init__(self, project=None, id_params=None):
        if project:
            if not os.path.isdir(project):
                raise OSError('Directory given does not exist.')

            filename = os.path.join(project, defaults.GALAXY_FILE)
            if not os.path.isfile(filename):
                raise OSError('The given file name is not in the directory:')

            # extract params from each of the rows in the given csvfile.
            # also remove empty params.
            with open(filename, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                id_params = {}
                for row in reader:
                    row_to_store = copy.deepcopy(row)
                    gal_id = row['id']
                    for param in row:
                        if not row[param]:
                            row_to_store.pop(param)
                    row_to_store.pop('id')  # avoid redundancy
                    id_params[gal_id] = row_to_store

            # convert all appropiate values to floats,
            for gal_id in id_params:
                for key, value in id_params[gal_id].iteritems():
                    try:
                        id_params[gal_id][key] = float(value)
                    except ValueError:
                        pass

        self.omit_fit = names.omit_fit
        self.id_params = id_params
        self.params = GParameters.convertId_Params(self.id_params)
        self.fit_params = GParameters.convertId_Params(self.id_params,
                                                       self.omit_fit)
        self.nfit_params = self.getNFitParams()
        self.ordered_fit_names = self.sortModelParamsNames()
        self.num_galaxies = len(self.id_params.keys())

    def getNFitParams(self):
        """Extract :attr:`nfit_params from :attr:`params` by noticing which
        parameters are in fit_params.
        """
        nfit_params = dict()
        for param in self.params:
            if param not in self.fit_params:
                nfit_params[param] = self.params[param]
        return nfit_params

    def sortModelParamsNames(self):
        """Return the keys of :attr:`params` in an ordered specified by
        :mod:`names.py`. And when having more than one galaxy, all the
        parameters from one of the galaxies are ordered together.
        """
        param_names = []
        for gal_id in self.id_params:
            galaxy_model = self.id_params[gal_id]['galaxy_model']
            for name in names.gal_models_parameters[galaxy_model]:
                for param in self.id_params[gal_id]:
                    if param not in self.omit_fit.get(galaxy_model, []):
                        if param == name:
                            param_names.append(param + '_' + str(gal_id))
        return param_names

    @staticmethod
    def convertId_Params(id_params, omit_fit={}):
        """Converts id_params to the format of :attr:`params`.

            Args:
                id_params(dict): Same as :attr:`id_params`
                omit_fit(dict): Dictionary that has the same purpose as
                                :attr:`omit_fit`

            Returns:
                A dictionary params(dict).
        """
        params = {}
        for gal_id in id_params:
            galaxy_model = id_params[gal_id]['galaxy_model']
            for param in id_params[gal_id]:
                if param not in omit_fit.get(galaxy_model, []):
                    params[param + '_' + str(gal_id)] = (
                        id_params[gal_id][param])
        return params

    @staticmethod
    def convertParams_Id(params):
        """Convert a dictionary params in the format of :attr:`params` to a
        dictionary in the format :attr:`id_params`
        """
        id_params = {}
        ids = []
        for param in params.keys():
            if param[-1] not in ids:
                ids.append(param[-1])  # appends last character of param

        for gal_id in ids:
            ID_params = {}
            for param in params.keys():
                if param[-1] == gal_id:
                    # slice last 2 characters to avoid '_1'
                    ID_params[param[:-2]] = params[param]
            id_params[gal_id] = ID_params

        return id_params

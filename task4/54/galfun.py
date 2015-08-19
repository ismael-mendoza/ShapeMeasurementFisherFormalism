#!/usr/bin/env python

import os

import csv

import defaults

import galsim

import copy

import names


def drawGalaxy(params):
    """Draw an image of a galaxy with noise and convoluted with a given psf if
    desired. Also sets
    GaussianNoise to the image if SNR is given.
    """
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

    # set shear
    if 'e1' and 'e2' in params:
        gal = gal.shear(e1=params['e1'], e2=params['e2'])
    elif 'eta1' and 'eta2' in params:
        gal = gal.shear(eta1=params['eta1'], eta2=params['eta2'])
    elif 'q' and 'beta' in params:
        gal = gal.shear(q=params['q'], beta=params['beta'] * galsim.radians)
    else:
        raise ValueError('The shear for the galaxy was not specified.')

    # shift galaxy.
    if 'x0' and 'y0' in params:
        gal = gal.shift(params['x0'], params['y0'])
    else:
        raise ValueError('The shift for the galaxy was not specified.')

    # generate psf.
    if params.get('psf_flux', 0) != 0:

        if params.get('psf_flux', 1) != 1:
            raise ValueError('I do not think you want a psf of flux not 1')

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
    """Draw each of the galaxy represented in gals_params, important to note:
    *Convolution can be done at the end since we will assume that blended galaxies are in the same field of view/stamp
    so that they are convoluted with the same psf and convolution is an additive operation (integrating is additive)
    *Convolution can also be done for each galaxy separately and then adding them, we will always assume galaxies have same psf.
    *Noise should only be added once to the image (otherwise two different noises combined), so it should go here and not in the original drawGalaxy function.
    """
    """Create image of galaxies given a set of parameters.

    Galaxies can be draw in different ways, most common is to pass in g_parameters, but also for the fits a combination of fit_params and nfit_params as a keyword argument is needed.

    Args:
        fit_params(dict):Dictionary containing the parameters of a galaxy that
        are relevant when doing a fit.

    Args:
        entry(astropy.table.Row): A single row from a galaxy :mod:`descwl.catalog`.
        dx_arcsecs(float): Horizontal offset of catalog entry's centroid from image center
            in arcseconds.
        dy_arcsecs(float): Vertical offset of catalog entry's centroid from image center
            in arcseconds.
        filter_band(str): The LSST filter band to use for calculating flux, which must
            be one of 'u','g','r','i','z','y'.

    Returns:
        :class:`Galaxy`: A newly created galaxy source model.

    Raises:
        SourceNotVisible: All of the galaxy's components are being ignored.
        RuntimeError: Catalog entry is missing AB flux value in requested filter band.
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


def addNoise(image, snr, noise_seed):
    # set gaussian noise with given noise seed.
    noisy_image = copy.deepcopy(image)  # do not alter original image.
    bd = galsim.BaseDeviate(noise_seed)
    noise = galsim.GaussianNoise(rng=bd)
    variance_noise = noisy_image.addNoiseSNR(noise, snr,
                                             preserve_flux=True)
    return noisy_image, variance_noise


class GParameters(object):

    """Class that manages given galaxy parameters."""

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
        nfit_params = dict()
        for param in self.params:
            if param not in self.fit_params:
                nfit_params[param] = self.params[param]
        return nfit_params

    def sortModelParamsNames(self):
        """Return param names in the order _1, _2,..., and in each
        subscript follow
        order of defaults.py, retuns list of ordered params names. With it we
        can change order just by changing order in names file gal_parameters.
        Avoids omit.
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
        """Uses a dictionary of the form id:params, where params is a dictionary of parameters and converts it to a dictionary of the form param_id:value(param)
        *omit (if given) is a list of strings that are not desired in the final dictionary that are in the original idDict,
        useful because to remove items later have to access it by param_1 for example.
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

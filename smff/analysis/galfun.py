"""Multipurpose module that contains important functions ranging from managing parameters of 
generated galaxies to extracing information from relevant files. 
"""
import os
import csv
import galsim
import copy
import math
import numpy as np 

import analysis.models as models
import analysis.defaults as defaults 


def get_galaxy_model(params):
    """Return the image of a single galaxy optionally drawn with a psf.

    Look at :mod:`analysis.models` to figure out which galaxy models and psf
    models are supported as well as their corresponding implemented 
    parameters. You can even add your own galaxies if desired. 

    Args:
        params(dict): Dictionary containing the information of a single
        galaxy where the keys is the name(str) of the parameter and the
        values are the values of the parameters.

    Returns:
        A :class:`galsim.GSObject`
        
    """

    galaxy_model = params['galaxy_model']
    gal_cls = models.getModelCls(galaxy_model)
    gal_model = gal_cls(params)

    final = gal_model.gal

    if params.get('psf_flux', 0) != 0:

        if params.get('psf_flux', 1) != 1:
            raise ValueError('I do not think you want a psf of flux not 1')

        psf_cls = models.getPsfModelCls(params['psf_model'])
        psf_model = psf_cls(params)

        final = galsim.Convolve([final, psf_model.psf])

    return final


def getGalaxiesModels(fit_params=None, id_params=None, g_parameters=None, **kwargs):
    """Return the model of a set of galaxies.

    One of the the following must be specified:
        fit_params (and nfit_params as **kwargs. Used specially in :mod:`runfits.py`).
        id_params
        g_parameters (from which id_params is extracted.)
    This function generates each of the galaxies specified in id_params and then
    sums them together to get a final galaxy.

    Args:
        fit_params(dict): Partial form of id_params that only includes the
                          parameters to be used for the fit.
                          For details, :class:`GParameters`
        id_params(dict): Dictionary containing each of the galaxies
                         parameters. For details, :class:`GParameters`
        g_parameters(:class:`GParameters`): An object containing different
                                            forms of the galaxy parameters.
        image(bool): If :bool:True returns an galsim.Image otherwise it 
        returns a np.array

    Returns:
        A :class:`galsim.GSObject`

    """

    if id_params is None and g_parameters is None:
        fit_params.update(kwargs)
        id_params = GParameters.convertParams_Id(fit_params)

    if g_parameters is not None:
        id_params = g_parameters.id_params

    gals = []

    for gal_id in id_params:
        gals.append(get_galaxy_model(id_params[gal_id]))
    
    return galsim.Add(gals)

def getOmitFit(id_params, omit):

    omit_fit = {} 

    for gal_id in id_params:
        params_omit = omit.get(gal_id,[])
        params = id_params[gal_id]
        galaxy_model = params['galaxy_model']
        cls = models.getModelCls(galaxy_model)
        obj = cls(params_omit=params_omit)
        omit_fit[gal_id] = obj.omit_fit

    return omit_fit

def addNoise(image, snr, noise_seed=0):
    """Set gaussian noise to the given galsim.Image.

    Args:
        image(:class:`galsim.Image`): Galaxy image that noise is going to be added
                             to.
        snr(float): Signal to noise ratio.
        noise_seed(int): Seed to set to galsim.BaseDeviate which
                         will create the galsim.noise instance.

    Returns:
        A :class:`galsim.Image`, variance_noise tuple. The image is the noisy version
        of the original image and variance_noise is the noise variance on each
        pixel due to the added noise.

    """

    noisy_image = copy.deepcopy(image)  # do not alter original image.
    bd = galsim.BaseDeviate(noise_seed)
    noise = galsim.GaussianNoise(rng=bd)
    variance_noise = noisy_image.addNoiseSNR(noise, snr, preserve_flux=True)
    return noisy_image, variance_noise

def read_results(project,g_parameters, fish,limit=None):
    orig_image = fish.image
    mins = defaults.getMinimums(g_parameters, orig_image)
    maxs = defaults.getMaximums(g_parameters, orig_image)
    
    residuals = {}
    pulls = {}
    redchis = [] #list containing values of reduced chi2 for each fit.
    rltsdir = os.path.join(project, defaults.RESULTS_DIR)

    files = os.listdir(rltsdir)
    if limit!=None: 
        files = files[:limit]

    # read results from rltsdir's files.
    for filename in files:
        with open(os.path.join(rltsdir, filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for i,row in enumerate(reader):
                redchis.append(float(row['redchi']))
                for param in g_parameters.fit_params:
                    if param not in residuals:
                        residuals[param] = []
                    if param not in pulls:
                        pulls[param] = []
                    residual = (float(row[param]) -
                                float(g_parameters.params[param]))

                    # if param in ['x0_1', 'y0_1']:
                    #     residual +=0.1 

                    pull = (residual /
                            math.sqrt(fish.covariance_matrix[param, param]))

                    residuals[param].append(residual)
                    pulls[param].append(pull)

    biases = {param: np.mean(residuals[param]) for param in residuals}
    pull_means = {param: np.mean(pulls[param]) for param in residuals}
    res_stds = {param: np.std(residuals[param]) for param in residuals}
    pull_mins = {param: ((mins[param] - float(g_parameters.params[param])) /
                         math.sqrt(fish.covariance_matrix[param, param])) for
                 param in residuals}
    pull_maxs = {param: ((maxs[param] - float(g_parameters.params[param])) /
                         math.sqrt(fish.covariance_matrix[param, param])) for
                 param in residuals}

    return pulls,residuals,biases,pull_means,res_stds,pull_mins,pull_maxs,redchis


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
            omit_fit(dict): Dictionary defined in containing
                the parameters that should not be included in the
                analysis for a particular galaxy model but could be 
                ecessary for drawing the galaxy. 
            id_params(dict): Dictionary whose keys are the ids of each of the
                galaxies specified in galaxies.csv, and that map
                to another dictionary that can be taken in by
                :func:`get_galaxy_model`
            params(dict): Dictionary that encodes the same information as
                id_params but in a different form. Combines each of
                the dictionaries contained in id_params into a
                single dictionary that contains all parameters but 
                in the form param_#.
            fit_params(dict): Dictionary similar to the params attribute but
                without the parameters specified in omit_fit.
            nfit_params(dict): Dictionary that contains all the parameters not
                contained in fit_params. Usually used for kwargs in conjunction with fit_params to 
                draw Galaxies.
            ordered_fit_names(list): A list containing the keys of fit_params
                in a desirable order.
            num_galaxies(int): Number of galaxies specified.

    """

    def __init__(self, project=None, id_params=None, omit={}):
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
                for key, value in id_params[gal_id].items():
                    try:
                        id_params[gal_id][key] = float(value)
                    except ValueError:
                        pass

        self.id_params = id_params
        self.params = GParameters.convertId_Params(self.id_params)
        self.omit_fit = getOmitFit(id_params, omit)
        self.fit_params = GParameters.convertId_Params(self.id_params,
                                                       self.omit_fit)
        self.nfit_params = self.getNFitParams()
        self.ordered_fit_names = self.sortModelParamsNames()
        self.num_galaxies = len(list(self.id_params.keys()))

    def getNFitParams(self):
        """Extract :attr:`nfit_params` from :attr:`params` by noticing which
        parameters are in fit_params.
        """
        nfit_params = dict()
        for param in self.params:
            if param not in self.fit_params:
                nfit_params[param] = self.params[param]
        return nfit_params

    def sortModelParamsNames(self):
        """Return the keys of :attr:`params` in an ordered specified by
        the class parameters. And when having more than one galaxy, all the
        parameters from one of the galaxies are ordered together.
        """
        param_names = []
        for gal_id in self.id_params:
            galaxy_model = self.id_params[gal_id]['galaxy_model']
            cls = models.getModelCls(galaxy_model)
            for name in cls.parameters:
                for param in self.id_params[gal_id]:
                    if param not in self.omit_fit.get(gal_id, []):
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
                A dictionary of params (dict).
        """
        params = {}
        for gal_id in id_params:
            for param in id_params[gal_id]:
                if param not in omit_fit.get(gal_id, []):
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

class ImageRenderer(object):
    """Object used to produce the image of a galaxy. 

    Args:
        pixel_scale(float): Pixel_scale to use in the image (ratio of pixels to arcsecs)
        nx(float): Width of the image in pixels.
        ny(float): height of the image in pixels. 
        stamp(int): optional, galsim.Image of appropiate dimensions to draw the image. does 
            not actually use whatever was originally in the stamp. 
        bounds(tuple): When drawn, the image will be clipped to these bounds. 
        mask(:class:`np.array`): the pixels selected in this mask will be set to 0. 
    
    One of the the following must be specified:
        * stamp
        * nx,ny,pixel_scale

    This object is made so it can be passsed in to a class :class:`analysis.fisher.Fisher` object. 
    """

    def __init__(self, pixel_scale=None,nx=None, ny=None,stamp=None, project=None,
                 bounds=None, mask=None):

        self.pixel_scale = pixel_scale
        self.nx = nx 
        self.ny = ny
        self.stamp = stamp 
        self.bounds = bounds 
        self.mask = mask

        if self.stamp is None:

            if self.nx is not None and self.ny is not None and self.pixel_scale is not None:
                self.stamp = galsim.Image(self.nx, self.ny, scale=self.pixel_scale)


        else: 
            self.nx = self.stamp.array.shape[0]
            self.ny = self.stamp.array.shape[1]
            self.pixel_scale = self.stamp.scale

        if self.bounds is not None:
            self.stamp = self.stamp[bounds]

    def getImage(self, galaxy):
        img = copy.deepcopy(self.stamp)
        galaxy.drawImage(image=img,use_true_center=False)

        if self.mask is None:
            return img
        else: 
            img.array[mask] = 0.
            return img
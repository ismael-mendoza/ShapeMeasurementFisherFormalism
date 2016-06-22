import os
import csv
import defaults
import galsim
import copy
import models
import math
import numpy as np 

#Todo - 
# remove drawGalaxies returning two different things. 
#Check if other functions have side effects,better to do only one thing. 
#drawGalaxies and Gparameters have different ways of initializing them? 
#confusing, better to only 
#have one way? 
#condense drawGalaxy and getGalaxyModel functions


def getGalaxyModel(params):
    """Return the image of a single galaxy optionally drawn with a psf.

    Look at :mod:`names.py` to figure out which galaxy models and psf
    models are supported as well as their corresponding implemented 
    parameters.

    Args:
        params(dict): Dictionary containing the information of a single
        galaxy where the keys is the name(str) of the parameter and the
        values are the values of the parametes.

    Returns:
        A galsim.GSObject
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
        image(bool): If :bool:True returns an galsim.Image otherwise it 
        returns a np.array

    Returns:
        A galsim.GSObject
    """

    if id_params is None and g_parameters is None:
        fit_params.update(kwargs)
        id_params = GParameters.convertParams_Id(fit_params)

    if g_parameters is not None:
        id_params = g_parameters.id_params

    gals = []

    for gal_id in id_params:
        gals.append(getGalaxyModel(id_params[gal_id]))
    
    return galsim.Add(gals)


#assume params_omit is a dictionary from gal_id to parameters to omit, 
#must to the same.
#for omit fit
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
                          single dictionary that contains all parameters but 
                          in the form param_#.
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
                for key, value in id_params[gal_id].iteritems():
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
                A dictionary params(dict).
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
    """
    Everything on how to produce the image of the galaxy. 
    """

    def __init__(self, pixel_scale=None,nx=None, ny=None,stamp=None, project=None,
                 mean_sky_level=None,min_snr=None,truncate_radius=None):

        #add possibility of obtaining nx, ny and pixel_scale from project. 

        self.pixel_scale = pixel_scale
        self.nx = nx 
        self.ny = ny
        self.stamp = stamp 
        self.mean_sky_level = mean_sky_level
        self.min_snr = min_snr
        self.truncate_radius = truncate_radius
        self.truncation_mask = None

        if self.stamp is None:

            if (self.pixel_scale is not None and self.mean_sky_level is not None 
                and self.min_snr is not None and self.truncate_radius is not None):
                self.truncation_mask = self.getTruncationMask()

            elif self.nx is not None and self.ny is not None and self.pixel_scale is not None:
                self.stamp = galsim.Image(self.nx, self.ny, scale=self.pixel_scale)


            else:
                raise OSError('Did not enough attributes for the image renderer.')


    def getTruncationMask(self):
        sky_noise = math.sqrt(self.mean_sky_level)
        self.pixel_cut = self.min_snr*sky_noise
        # We will render each source into a square stamp with width = height = 2*padding + 1.
        self.padding = int(math.ceil(self.truncate_radius/self.pixel_scale - 0.5))
        size = 2*self.padding + 1
        self.stamp = galsim.Image(size,size,scale = self.pixel_scale, dtype = np.float32)
        # Prepare a truncation mask.
        pixel_grid = np.arange(-self.padding,self.padding+1)*self.pixel_scale
        pixel_x,pixel_y = np.meshgrid(pixel_grid,pixel_grid)
        pixel_radius = np.sqrt(pixel_x**2 + pixel_y**2)
        return (pixel_radius <= self.truncate_radius)

    def get_image_coordinates(image_width,dx_arcsecs,dy_arcsecs):
        """Convert a physical offset from the image center into image coordinates.

        Args:
            dx_arcsecs(float): Offset from the image center in arcseconds.
            dy_arcsecs(float): Offset from the image center in arcseconds.

        Returns:
            tuple: Corresponding floating-point image coordinates (x_pixels,y_pixels)
                whose :func:`math.floor` value gives pixel indices and whose...
        """
        x_pixels = 0.5*self.image_width + dx_arcsecs/self.pixel_scale
        y_pixels = 0.5*self.image_height + dy_arcsecs/self.pixel_scale
        return x_pixels,y_pixels

    def getCroppedBounds(self):
        keep_mask = (self.stamp.array*self.truncation_mask > self.pixel_cut)
        if np.sum(keep_mask) == 0:
            raise SourceNotVisible
        self.stamp.array[np.logical_not(keep_mask)] = 0.

        x_projection = (np.sum(keep_mask,axis=0) > 0)
        y_projection = (np.sum(keep_mask,axis=1) > 0)
        x_min_inset = np.argmax(x_projection)
        x_max_inset = np.argmax(x_projection[::-1])
        y_min_inset = np.argmax(y_projection)
        y_max_inset = np.argmax(y_projection[::-1])
        return galsim.BoundsI(
               x_min+x_min_inset,x_max-x_max_inset,
               y_min+y_min_inset,y_max-y_max_inset)

    def getImage(self, galaxy):

        galaxy.drawImage(image=self.stamp)
        if self.truncation_mask is not None:
            bounds = self.getCroppedBounds()
            self.stamp = self.stamp[bounds]
            return self.stamp

        return copy.deepcopy(self.stamp)


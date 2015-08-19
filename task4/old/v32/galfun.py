#!/usr/bin/env python
import os

import csv

import defaults

import galsim

import plotsfisher

import matplotlib.pyplot as plt

import copy

import names


def sortParamsNames(id_params, omit):
    """Returns params in the order _1, _2,..., and in each subscript follow 
    order of defaults.py, retuns list of ordered params names. With it we 
    can change order just by changing order in defaults file. 
    """
    param_names = []
    for params_id in id_params.keys():
        for name in names.gal_parameters:
            for param in id_params[params_id].keys():
                if param not in omit:
                    if param == name:
                        param_names.append(param + '_'+ str(params_id))
    return param_names

def convertParamsId(id_params, omit = []): 
    """Uses a dictionary of the form id:params, where params is a dictionary of parameters and converts it to a dictionary of the form param_id:value(param)
    *omit (if given) is a list of strings that are not desired in the final dictionary that are in the original idDict,
    useful because to remove items later have to access it by param_1 for example.
    """

    params = {}
    for params_id in id_params:
        for param in id_params[params_id].keys():
            if param not in omit:
                params[param + '_'+ str(params_id)] = id_params[params_id][param]
    return params


def convertParams(params):
    id_params ={}
    ids = []
    for param in params.keys():
        if param[-1] not in ids:
            ids.append(param[-1]) #appends last character of param

    for ID in ids:
        ID_params = {}
        for param in params.keys():
            if param[-1] == ID:
                #slice last 2 characters to avoid '_1'
                ID_params[param[:-2]] = params[param]
        id_params[ID] = ID_params

    return id_params


def drawGalaxy(params): 
    """Draws an image of a galaxy with noise and convoluted with a given psf if desired. Also sets 
    GaussianNoise to the image if SNR is given.
    """

    #check each model of galaxy
    if(params['model'] == 'gaussian'):
        if('hlr' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], half_light_radius=params['hlr'])
        elif('sigma' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

    elif(params['model'] == 'exponential'):
        gal = galsim.Exponential(flux = params['flux'], half_light_radius = params['hlr'])

    #set shear
    if('e1' and 'e2' in params.keys()): 
        gal = gal.shear(e1=params['e1'], e2 = params['e2'])
    elif('q' and 'beta' in params.keys()):
        gal = gal.shear(q=params['q'], beta = params['beta'] * galsim.radians)

    #shift galaxy.
    gal = gal.shift(params['x0'],params['y0'])

    #add a  psf to the galaxy, if psf_flux is not 0
    if(params['psf_flux'] != 0):

        #each psf model. 
        if(params['psf_model'] == 'gaussian'):
            if('psf_hlr' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], half_light_radius=params['psf_hlr'])
            elif('psf_sigma' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], 
                                      sigma=params['psf_sigma'])
            elif('psf_fwhm' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], 
                                      fwhm=params['psf_fwhm'])

        psf = psf.shear(e1=params['psf_e1'],e2=params['psf_e2'])

        final = galsim.Convolve([gal, psf])
    else:
        final = gal

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    return final.drawImage(scale=defaults.PIXEL_SCALE, nx=defaults.NX, 
                            ny=defaults.NY)


def drawGalaxies(params=None, id_params=None, image=False, 
                 g_parameters=None, **kwargs):
    """Draw each of the galaxy represented in gals_params, important to note: 
    *Convolution can be done at the end since we will assume that blended galaxies are in the same field of view/stamp 
    so that they are convoluted with the same psf and convolution is an additive operation (integrating is additive)
    *Convolution can also be done for each galaxy separately and then adding them, we will always assume galaxies have same psf. 
    *Noise should only be added once to the image (otherwise two different noises combined), so it should go here and not in the original drawGalaxy function. 
    """

    """ Will return an image of the galaxies that are specified by params.
    Either one (but not both) of params and id_params must be specified. 

    Args
    ----
    @param params       Dictionary of parameters of the form: 'flux_1:#'
    @param id_params        A list of tuples representing the peak 
                        positions objects in the blend.
    @param interpolate  If at least one component of rot_center is not a half-integer, use GalSim
                        to rotate the image.  This currently doesn't work very well!!!
    @param force_interpolate   Use GalSim to rotate the image, even if rot_center components are
                               half-integer and rotation via numpy array operations is possible.
                               This currently doesn't work very well!!!
    
    @returns templates, template_fractions, children
    
    """

    if id_params is None and g_parameters is None:
        params.update(kwargs)
        id_params = convertParams(params)

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


def displayGalaxy(image):
    """Displays a galaxy given an image of a galaxy, mostly used just for 
    testing.
    """
    figure, subplt= plt.subplots(1,1)
    figure.suptitle('Galaxy', fontsize = defaults.FONTSIZE_TITLE)
    plotsfisher.drawImage(subplt, image.array)
    plt.show()


def addNoise(image, snr, noise_seed):
    #set gaussian noise with given noise seed.
    noisy_image = copy.deepcopy(image) # do not alter original image.
    bd = galsim.BaseDeviate(noise_seed)
    noise = galsim.GaussianNoise(rng=bd)
    variance_noise = noisy_image.addNoiseSNR(noise, snr,  
                                            preserve_flux =True) 
    return noisy_image, variance_noise


class GParameters:
    """Class that manages given galaxy parameters."""
    def __init__(self, project): 

        #some error handling.
        if not os.path.isdir(project):
            raise OSError('Directory given does not exist.')

        filename = os.path.join(project, defaults.GALAXY_FILE)
        if not os.path.isfile(filename):
            raise OSError('The given file name is not in the directory:')

        #extract params from each of the rows (galaxies) in the given csvfile.
        with open(filename, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            id_params = {} #dictionary that contains all parameters of each of the galaxies inputted by the user. 
            for row in reader:
                gal_id = row.pop('id') #remove id entry from row. 
                id_params[gal_id] = row

        #convert all appropiate values to floats and change names of each corresponding galaxy parameter.
        for gal_id in id_params: 
            for key, value in id_params[gal_id].iteritems():
                try:
                    id_params[gal_id][key] = float(value)
                except ValueError:
                    pass

        self.id_params = id_params 
        self.params = convertParamsId(self.id_params)
        self.model_params = convertParamsId(self.id_params, 
                                            names.extra + names.psf_parameters
                                            )
        #in order.
        self.model_param_names = sortParamsNames(self.id_params, 
                                           names.extra + names.psf_parameters)
        self.num_galaxies = len(self.id_params.keys())


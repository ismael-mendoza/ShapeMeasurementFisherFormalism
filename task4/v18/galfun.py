#!/usr/bin/env python
import os

import csv

import defaults

import galsim

import plotsfisher

import matplotlib.pyplot as plt

cts = defaults.constants()

def convertParamsId(id_params, omit = []): 
    """Uses a dictionary of the form id:params, where params is a dictionary of parameters and converts it to a dictionary of the form param_id:value(param)
    *omit (if given) is a list of strings that are not desired in the final dictionary that are in the original idDict,
    useful because to remove items later have to access it by param_1 for example."""

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

def drawGalaxy(params, **kwargs): 
    """Draws an image of a galaxy with noise and convoluted with a given psf if desired. Also sets 
    GaussianNoise to the image if SNR is given."""

    params.update(kwargs)

    if(params['model'] == 'gaussian'):
        if('hlr' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], half_light_radius=params['hlr'])
        elif('sigma' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

    elif(params['model'] == 'exponential'):
        gal = galsim.Exponential(flux = params['flux'], half_light_radius = params['hlr'])

    #set shear, if desired. note that shear generates a new sheared galaxy, should not try to change attribute of the galaxy directly
    if('e1' and 'e2' in params.keys()): 
        gal = gal.shear(e1=params['e1'], e2 = params['e2'])
    elif('q' and 'beta' in params.keys()):
        gal = gal.shear(q=params['q'], beta = params['beta'] * galsim.radians) ##galsim.radians is only useful when you draw it (sheart it) do not need it anywhereelse.


    #shift galaxy if desired. be careful to do shear and shift in this order, because otherwise might act weirdly (takes center as the original for shear)
    gal = gal.shift(params['x0'],params['y0'])

    #add a  psf if psf_flux is not 0
    if(params['psf_flux'] != 0):

        if(params['psf_model'] == 'gaussian'):
            if('psf_hlr' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], half_light_radius=params['psf_hlr'])
            elif('psf_sigma' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], sigma=params['psf_sigma'])
            elif('psf_fwhm' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], fwhm=params['psf_fwhm'])
            final = galsim.Convolve([gal, psf])
    else:
        final = gal

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    return final.drawImage(scale=cts.pixel_scale, nx = cts.nx, ny = cts.ny)


def drawGalaxies(params = None, id_params = None, image = False, galaxies = None, display = False, **kwargs):
    """Draw each of the galaxy represented in gals_params, important to note: 
    *Convolution can be done at the end since we will assume that blended galaxies are in the same field of view/stamp 
    so that they are convoluted with the same psf and convolution is an additive operation (integrating is additive)
    *Convolution can also be done for each galaxy separately and then adding them, we will always assume galaxies have same psf. 
    *Noise should only be added once to the image (otherwise two different noises combined), so it should go here and not in the original drawGalaxy function. 
    """

    if id_params is None and galaxies is None:
        id_params = convertParams(params)

    if galaxies is not None:
        id_params = galaxies.id_params

    gals = []
    
    for gal_id in id_params:
        gals.append(drawGalaxy(id_params[gal_id]))
    final = sum(gals) #will sum all galaxy images in gals.

    if image is False:
        return final.array
    else:
        return final

    if display:
        figure, subplt= plt.subplots(1,1)
        figure.suptitle('Galaxy', fontsize = 20)
        plotsfisher.drawImage(subplt, final.array)
        plotsfisher.SaveFigureToPdf(figure, 'test.png', 'output', 'plots')

class galaxies:
    def __init__(self, wdir = None, galaxy_file = None, id_params = None, params = None): 
        #initialize a the object in 3 ways, by reading the galaxies' paramaters info from the given file, by given id_params dictionary or parameters dictionary. 

        names = defaults.names()

        if wdir is not None and galaxy_file is not None:
            #some error handling.
            if not os.path.isdir(wdir):
                print ('Directory does not exists')
                return -1

            filename = os.path.join(wdir, galaxy_file + '.csv')
            if not os.path.isfile(filename):
                print('Galaxies file does not exist')
                return -1

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

        elif id_params is not None:
            self.id_params = id_params #dictionary of the form {id:params}, key is id, values are params of each galaxy, the params do not have id as a key anymore.
            self.params = convertParamsId(self.id_params)

        else:
            self.params = params
            self.id_params = convertParams(self.params)

        self.model_params = convertParamsId(self.id_params, omit = names.extra + names.psf_parameters)
        self.image = drawGalaxies(self.params, image = True)

    def addNoise(self, snr, noise_seed):
        #set gaussian noise with given noise seed.
        bd = galsim.BaseDeviate(noise_seed)
        noise = galsim.GaussianNoise(rng=bd)
        variance_noise = self.image.addNoiseSNR(noise, snr, preserve_flux =True) #this will presumably change the image of the galaxies stored in the object.
        return variance_noise
    





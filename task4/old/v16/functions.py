#!/usr/bin/env python
"""Contains various functions used throughout the module, may need to look for a better way to organize them specially drawGalaxy will become 
quite large soon. """

import os

import defaults

import matplotlib.cm as cm

from copy import deepcopy

import galsim

import csv

import matplotlib.pyplot as plt

from matplotlib.patches import Ellipse

# def errorEllipse(pos,cov,nstd):

cts = defaults.constants()
names = defaults.names()


class galsParameters:

    def __init__(self, wdir, galaxy_file):
        #initialize a the object by reading the galaxies' paramaters info from the given file. 

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
            gals_params = {} #dictionary that contains all parameters of each of the galaxies inputted by the user. 
            for row in reader:
                gal_id = row.pop('id') #remove id entry from row. 
                gals_params[gal_id] = row

        #convert all appropiate values to floats and change names of each corresponding galaxy parameter.
        for gal_id in gals_params: 
            for key, value in gals_params[gal_id].iteritems():
                try:
                    gals_params[gal_id][key] = float(value)
                except ValueError:
                    pass

        self.id_params = gals_params #dictionary of the form {id:params}, key is id, values are params of each galaxy, the params do not have id as a key anymore.
        self.params = self.convertParamsId()
        self.model_params = self.convertParamsId(omit = names.extra)

        def convertParamsId(self, omit = []): 
            """Uses a dictionary of the form id:params, where params is a dictionary of parameters and converts it to a dictionary of the form param_id:value(param)
            *omit (if given) is a list of strings that are not desired in the final dictionary that are in the original idDict,
            useful because to remove items later have to access it by param_1 for example."""

            params = {}
            for params_id in self.id_params:
                for param in self.id_params[params_id].keys():
                    if param not in omit:
                        params[param + '_'+ str(params_id)] = self.id_params[params_id][param]


def csvIsEmpty(filename):
    """checks each row and if any is not empty, then the file is not empty"""
    with open(filename, 'r') as f:
        for row in f:
            if len(row) != 0:
                return False
        else: return True

def SaveFigureToPdf(figure,file_name, dir_name, plotdir_name, hide = True): 

     #save and preview pdf.
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    if not os.path.isdir(os.path.join(dir_name,plotdir_name)):
        os.mkdir(os.path.join(dir_name,plotdir_name))
    file_name = os.path.join(dir_name, plotdir_name, file_name)  #puts slashes in between things.
    figure.savefig(file_name, bbox_inches='tight')

    if not hide: 
        os.system("open " + file_name)

def drawImage(axis, plot, title = "", xlabel = "", ylabel = ""): 
    """draws a given plot with default values using imshow in the given subplot
    if they are not subplots should just pass plt. It also adds a title with name. -- OWN"""

    axis.axes.get_xaxis().set_ticks([]) #remove numbers and ticks but can use label.
    axis.axes.get_yaxis().set_ticks([])

    axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)

    #be careful none must be in quotes. 
    #vmax and vmin useful to center around zero which looks nice.
    axis.imshow(plot, interpolation = 'None', cmap=cm.RdYlGn,
            origin='lower', extent=[-cts.nx,cts.nx,-cts.ny,cts.ny],
            vmax=abs(plot).max(), vmin=-abs(plot).max())

def partialDifferentiate(func, parameter, step, **kwargs):
    """Partially derive f with respect to a parameter with a certain step.
    We are assuming that the function has a certain structure, namely, one of its arguments is a dictionary of 
    variables that can be changed and other (**kwargs) arguments are requisites or extra variables
    the function needs to be evaluated. This is because we are assuming we can add step to params[parameter]."""

    def Dfunc(params):
        """Evaluate the partial derivative at params."""
        params_up = deepcopy(params) #avoids altering params later.
        params_up[parameter] += step #increment the value of the parameter by step. 

        params_down = deepcopy(params)
        params_down[parameter] -= step 

        return (func(params_up, **kwargs) - func(params_down, **kwargs)) / (2 * step)
    return Dfunc

def secondPartialDifferentiate(func, parameter1, parameter2, step1, step2, **kwargs): 
    Df = partialDifferentiate(func, parameter1, step1, **kwargs)
    return partialDifferentiate(Df, parameter2, step2)

def chi2(params, gal_image, sigma_n, **kwargs): 
    """Returns chi2 given the modified parameters and the original galaxy, assume sigma_n is the same for all pixels -- OWN"""
    return ((((gal_image- drawGalaxy(params, **kwargs)).array/ (sigma_n)))**2).sum()

def drawGalaxy(params, **kwargs): 
    """Draws an image of a galaxy with noise and convoluted with a given psf if desired. Also sets 
    GaussianNoise to the image if SNR is given."""

    params.update(kwargs)

    #have to do this for each type of galaxy
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

    #add a gaussian psf if its flux it not 0
    if(params['psf_flux'] != 0):
        #have to do this for each type of psf. 
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
    return final.drawImage(scale=cts.pixel_scale, nx = cts.nx, ny = cts.ny).array


def drawGalaxies(gals_params, snr = -1, noise_seed = 0, **kwargs):
    """Draw each of the galaxy represented in gals_params, important to note: 
    *Convolution can be done at the end since we will assume that blended galaxies are in the same field of view/stamp 
    so that they are convoluted with the same psf and convolution is an additive operation (integrating is additive)
    *Convolution can also be done for each galaxy separately and then adding them, we will always assume galaxies have same psf. 
    *Noise should only be added once to the image (otherwise two different noises combined), so it should go here and not in the original drawGalaxy function. 
    """
    gals = []
    id_params = gals_params.id_params
    
    for gal_id in id_params:
        gals.append(drawGalaxy(id_params[gal_id], **kwargs))
    final = sum(gals) #will sum all galaxy arrays in gals.

    #careful with noise, must be added only once. 
    #add noise if snr is specified.
    if(snr != -1):
        #set gaussian noise with given noise seed.
        bd = galsim.BaseDeviate(noise_seed)
        noise = galsim.GaussianNoise(rng=bd)
        variance_noise = final.addNoiseSNR(noise, snr, preserve_flux =True)
        return  final, variance_noise
    else: 
        return final

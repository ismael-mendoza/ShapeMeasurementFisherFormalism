import sys
import os
import math
import galsim
from lmfit import minimize, Parameters, Parameter, fit_report
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from constantsTask3 import *
from copy import deepcopy #useful to not change both variables when one is copied.
import gc  

def bias(paramFunc, params, param_names, biases_of_params, steps):
    pass

def variance(paramFunc1, paramFunc2, params, param_names, CovM_of_params, steps):
    """Uses the existing covariance matrix of the initial parameters to compute the variance of other parameters that were not included initially in the model"""
    var = 0 
    for k in range(num_params):
        for l in range(num_params):
            var += partialDifferentiate(paramFunc1, param_names[k], steps[param_names[k]])(params) * partialDifferentiate(paramFunc2, param_names[l], steps[param_names[l]])(params) * CovM_of_params[param_names[k],param_names[l]]
    return var

def radiusMean_func(params):
    return math.sqrt(a1_func(params)**2 + a2_func(params)**2)

def q_func(params): 
    if('q' in params.keys()): 
        return params['q']

    elif('e1' and 'e2' in params.keys()): 
        return math.sqrt((1-math.sqrt(params['e1']**2 + params['e2']**2))/(1 + math.sqrt(params['e1']**2 + params['e2']**2)))

def beta_func(params):
    """Calculates the galaxy parameter beta in terms of the other 6 initial known parameters"""

    if('beta' in params.keys()):
        return params['beta']

    elif('e1' and 'e2' in params.keys()):
        return .5 * math.atan(params['e2'] /params['e1']) #should be in radians. 

def a1_func(params):
    """Calculates the galaxy parameter a1 in terms of the other 6 initial known parameters"""
    return a2_func(params) / q_func(params)

def a2_func(params):
    """Calculates the galaxy parameter a2 in terms of the other 6 initial known parameters"""

    if('sigma' in params.keys()):
        return params['sigma'] * math.sqrt(q_func(params))
    elif('hlr' in params.keys()):
        pass

def amplitude_func(params): 
    """Calculates the galaxy parameter amplitude in terms of the other 6 initial known parameters"""

    return params['flux']/params['sigma']

def SaveFigureToPdfAndOpen(figure,file_name): 

     #os to save and preview pdf.
    if not os.path.isdir('figures'):
            os.mkdir('figures')

    file_name = os.path.join('figures', file_name)  #puts slashes in between things.
    figure.savefig(file_name, bbox_inches='tight')

    os.system("open " + file_name)


def drawPlot(figure, plot, title = "", xlabel = "", ylabel = ""): 
    """draws a given plot with default values using imshow in the given subplot
    if they are not subplots should just pass plt. It also adds a title with name. -- OWN"""

    figure.axes.get_xaxis().set_ticks([]) #remove numbers and ticks but can use label.
    figure.axes.get_yaxis().set_ticks([])

    figure.set_title(title)
    figure.set_xlabel(xlabel)
    figure.set_ylabel(ylabel)

    #be careful none must be in quotes. 
    #vmax and vmin useful to center around zero which looks nice.
    figure.imshow(plot, interpolation = 'None', cmap=cm.RdYlGn,
            origin='lower', extent=[-nx,nx,-ny,ny],
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

def drawGalaxy(gal_params, psf_params = {}, SNR = -1, noise_seed = 0): 
    """Draws an image of a galaxy with noise and convoluted with a given psf if desired."""

    # Define the galaxy profile
    if('hlr' in gal_params.keys()):
        gal = galsim.Gaussian(flux=gal_params['flux'], half_light_radius=gal_params['hlr'])
    elif('sigma' in gal_params.keys()):
        gal = galsim.Gaussian(flux=gal_params['flux'], sigma=gal_params['sigma'])

    #set shear, if desired. note that shear generates a new sheared galaxy, should not try to change attribute of the galaxy directly
    if('e1' and 'e2' in gal_params.keys()): 
        gal = gal.shear(e1=gal_params['e1'], e2 = gal_params['e2'])
    elif('q' and 'beta' in gal_params.keys()):
        gal = gal.shear(q=gal_params['q'], beta = gal_params['beta'] * galsim.radians) ##galsim.radians is only useful when you draw it (sheart it) do not need it anywhereelse.

    #shift galaxy if desired. be careful to do shear and shift in this order, because otherwise might act weirdly (takes center as the original for shear)
    gal = gal.shift(gal_params['x0'],gal_params['y0'])

    #add a gaussian psf.
    if(len(psf_params) != 0):
        if('hlr' in psf_params.keys()):
            psf = galsim.Gaussian(flux=psf_params['flux'], half_light_radius=psf_params['hlr'])
        elif('sigma' in psf_params.keys()):
            psf = galsim.Gaussian(flux=psf_params['flux'], sigma=psf_params['sigma'])
        elif('fwhm' in psf_params.keys()):
            psf = galsim.Gaussian(flux=psf_params['flux'], fwhm=psf_params['fwhm'])
        
        final = galsim.Convolve([gal, psf])
    else:
        final = gal

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    image = final.drawImage(scale=pixel_scale, nx = nx, ny = ny)

    #set noise 
    if(SNR != -1):
        #set noise seed 
        bd = galsim.BaseDeviate(noise_seed)
        noise = galsim.GaussianNoise(rng=bd)
        variance_noise = image.addNoiseSNR(noise, SNR, preserve_flux =True)
        return  image, variance_noise

    return image
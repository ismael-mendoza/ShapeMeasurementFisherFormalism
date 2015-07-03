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

        return (func(params = params_up, **kwargs) - func(params =params_down, **kwargs)) / (2 * step)

    return Dfunc

def secondPartialDifferentiate(func, parameter1, parameter2, step1, step2, **kwargs): 
    Df = partialDifferentiate(func, parameter1, step1, **kwargs)
    return partialDifferentiate(Df, parameter2, step2)


def chi2(params, gal_image, sigma_n): 
    """Returns chi2 given the modified parameters and the original galaxy, assume sigma_n is the same for all pixels -- OWN"""
    return ((((gal_image- drawGalaxy(params)).array/ (sigma_n)))**2).sum()

def drawGalaxy(params): 
    """Draws image with noise when are simulating the galaxy and without noise when we are creating the model. 
    This is a good way to recycle code. -- OWN"""
    # Define the galaxy profile
    gal = galsim.Gaussian(flux=params['gal_flux'], sigma=params['gal_sigma'])

    #set shear, if desired. note that shear generates a new sheared galaxy, should not try to change attribute of the galaxy directly
    gal = gal.shear(e1=params['e1'], e2 = params['e2'])

    #shift galaxy if desired. be careful to do shear and shift in this order, because otherwise might act weirdly (takes center as the original for shear)
    gal = gal.shift(params['x0'],params['y0'])

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    # The returned image has a member, added_flux, which is gives the total flux actually added to 
    # the image.  One could use this value to check if the image is large enough for some desired
    # accuracy level.  Here, we just ignore it.
    image = gal.drawImage(scale=pixel_scale, nx = nx, ny = ny)

    #variance_noise = image.addNoiseSNR(galsim.GaussianNoise(), 50, True) 

    return image

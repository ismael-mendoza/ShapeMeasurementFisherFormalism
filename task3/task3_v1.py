"""We want to draw a fisher matrix for a resulting galaxy and its parameters. We do everything with a Gaussian in this case."""

import sys
import os
import math
import logging
import galsim
from lmfit import minimize, Parameters, Parameter, fit_report
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



#global constants for image (never going to use in a fit): (do not put in dictionary)
pixel_scale = .4
nx = 40
ny = 40


def drawPlot(subplt, plot, name): 
    """draws a given plot with default values using imshow in the given subplot
    if they are not subplots should just pass plt. It also adds a title with name."""

    #be careful none must be in quotes. 
    #vmax and vmin useful to center around zero which looks nice. 
    subplt.imshow(plot, interpolation = 'None', cmap=cm.RdYlGn,
            origin='lower', extent=[-nx,nx,-ny,ny],
            vmax=abs(plot).max(), vmin=-abs(plot).max())

    #for removing y and x axis. 
    subplt.axis('off') 
    subplt.set_title(name)

def partialDeriveGal(gal_image, params, parameter, step):
    """Calculates the partial derivative of a galaxy with respect to one of its parameters.
    Uses average of adding step and subtracting step (for some reason this is CRITICAL)"""
   
   #used to change each parameter slightly and draw an image with that. 
    modParams = params

    initValue = params[parameter].value

    modParams.add(parameter, value = initValue + step)
    mod_gal_imageUp = drawGalaxy(modParams)

    modParams.add(parameter, value = initValue - step)
    mod_gal_imageDown = drawGalaxy(modParams)

    return (mod_gal_imageUp - mod_gal_imageDown).array / (2* step)

    #return (gal_image - mod_gal_image).array / step, had it before did not work. 

def drawGalaxy(Params): 
    """Draws image with noise when are simulating the galaxy and without noise when we are creating the model. 
    This is a good way to recycle code."""

    params = Params.valuesdict()

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

    return image

def main(argv):
    """
    First draw a galaxy with the following characteristics, 
      - Use a circular Gaussian profile for the galaxy.
      - Convolve it by a circular Gaussian PSF.
      - Add Gaussian noise to the image.
      - It is sheared so that e1 != 0 and e2!= 0  and shifted. 

      Then we change a parameter and obtain a new image for each parameter that is slightly changed so we can calculate a derivative.
    """

    #first we want to set here the true values. Or the ones we want our galaxy to have.    
    #possible parameters for the galaxy formation. 

    #to change and calculate derivatives. 
    params = Parameters()
    params.add('gal_sigma', value = 1.) # arcsec 
    params.add('gal_flux', value = 100.)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    params.add('e1', value = 0) #ellipticity: e1 
    params.add('e2', value = 0)#ellipticity: e2
    params.add('x0', value = 0.) #shift in x origin. 
    params.add('y0', value = 0.)     #shift in y

    #get image of the original galaxy
    gal_image = drawGalaxy(params)

    params = Parameters()
    params.add('gal_sigma', value = 3.) # arcsec 
    params.add('gal_flux', value = 100.)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    params.add('e1', value = .1) #ellipticity: e1 
    params.add('e2', value = -.5)#ellipticity: e2
    params.add('x0', value = 2.) #shift in x origin. 
    params.add('y0', value = 2.)     #shift in y

    gal_image2 = drawGalaxy(params)

    #can draw multiple images at once by calling subplots multiple times.
    f, subplt = plt.subplots(1,1)
    drawPlot(subplt,  (gal_image + gal_image+ gal_image2).array, 'Initial Galaxy')

  
    #first we want only derivatives of the galaxy with respect to each parameter numerically so 6 plots in a row
    #Derivatives
    #define subplots 
    f, subplts = plt.subplots(1,6) #first number is the rows and the second number is the columns.

    #name for derivatives. 
    Dgal_sigma = partialDeriveGal(gal_image, params, 'gal_sigma', params['gal_sigma'].value * .0001)
    Dgal_flux = partialDeriveGal(gal_image, params, 'gal_flux', params['gal_flux'].value * .0001)
    De1 = partialDeriveGal(gal_image, params, 'e1', 1 * .0001)
    De2 = partialDeriveGal(gal_image, params, 'e2', 1 * .0001)
    Dx0 = partialDeriveGal(gal_image, params, 'x0', nx * .0001)
    Dy0 = partialDeriveGal(gal_image, params, 'y0', ny * .0001)

    #list with derivatives and list with their names. 
    names = ['Dgal_sigma', 'Dgal_flux', 'De1', 'De2', 'Dx0', 'Dy0'] 
    Ds = [Dgal_sigma, Dgal_flux, De1, De2, Dx0, Dy0]

    for i, D, name in zip(range(len(Ds)), Ds,names):
        drawPlot(subplts[i], D, name)

    # #labels and titles 
    # plt.xlabel('Smarts')
    # plt.ylabel('Something')
    # plt.title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')


    #be carefule none must be in quotes.
    # plt.imshow(gal_plot, interpolation = 'None', cmap=cm.RdYlGn,
    #             origin='lower', extent=[-nx,nx,-ny,ny],
    #             vmax=abs(plot).max(), vmin=-abs(plot).max()) 


    # plt.colorbar()

    #display the plots
    plt.show()






    #after havine one derivative for each down... 
        #do in a triangular plot fashion like in the presentation 
        #do not forget about squaring it too and multiplying to obtain fisher matrix eleemnts
        #also integrate over all pixels. 
        #do second derivatives to calculate bias matrix 
        #want to check answers analitically with the formulas of the paper.



if __name__ == "__main__":
    main(sys.argv)

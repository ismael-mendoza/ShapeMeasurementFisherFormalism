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
    if they are not subplots should just pass plt. It also adds a title with name. -- OWN"""

    # #Remove x axis but leave possibility for a x labe. 
    # subplt.tick_params(
    #     axis='y',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom='off',      # ticks along the bottom edge are off
    #     top='off',         # ticks along the top edge are off
    #     labelbottom='off') # labels along the bottom edge are off

    # #remove y axis completely including label. 
    # subplt.get_xaxis().set_visible(False)
    subplt.set_title(name)
    subplt.axes.get_xaxis().set_ticks([])
    subplt.axes.get_yaxis().set_ticks([])

    #be careful none must be in quotes. 
    #vmax and vmin useful to center around zero which looks nice. 
    subplt.imshow(plot, interpolation = 'None', cmap=cm.RdYlGn,
            origin='lower', extent=[-nx,nx,-ny,ny],
            vmax=abs(plot).max(), vmin=-abs(plot).max())

def partialDeriveGal(gal_image, params, parameter, step):
    """Calculates the partial derivative of a galaxy with respect to one of its parameters.
    Uses average of adding step and subtracting step (for some reason this is CRITICAL) -- OWN"""
   
   #used to change each parameter slightly and draw an image with that. 
    modParams = params

    parameter_initValue = params[parameter].value

    modParams.add(parameter, value = parameter_initValue + step)
    mod_gal_imageUp = drawGalaxy(modParams)

    modParams.add(parameter, value = parameter_initValue - step)
    mod_gal_imageDown = drawGalaxy(modParams)

    return (mod_gal_imageUp - mod_gal_imageDown).array / (2* step)

def PartialDeriveChi2(chi2, params, parameter, step):
    """partial derive chi2 with respect to a parameter --OWN """
    modParams = params

    parameter_initValue = params[parameter].value

    modParams.add(parameter, value = parameter_initValue + step)
    chi2_up = Chi2(modParams)

    modParams.add(parameter, value = parameter_initValue - step)
    chi2_down = Chi2(modParams)

    return (chi2_up - chi2_down).array/ (2* step)

def Chi2(gal_image, params, sigma_n): 
    """Returns chi2 given the modified parameters and the original galaxy, assume sigma_n is the same for all pixels -- OWN"""
    return ((((gal_image - drawGalaxy(params)).array/ (sigma_n)))**2).sum()

def drawGalaxy(Params): 
    """Draws image with noise when are simulating the galaxy and without noise when we are creating the model. 
    This is a good way to recycle code. -- OWN"""

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
    orig_params = Parameters()
    orig_params.add('x0', value = 0.) #shift in x origin. 
    orig_params.add('y0', value = 0.)     #shift in y
    orig_params.add('gal_sigma', value = 3.) # arcsec 
    orig_params.add('gal_flux', value = 100.)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    orig_params.add('e1', value = .1) #ellipticity: e1 
    orig_params.add('e2', value = -.5)#ellipticity: e2


    #get image of the original galaxy
    gal_image = drawGalaxy(orig_params)

    #can draw multiple images at once separately by calling subplots multiple times.
    f, subplt = plt.subplots(1,1)
    drawPlot(subplt, gal_image.array, 'Initial Galaxy')

  
    #first we want only derivatives of the galaxy with respect to each parameter numerically so 6 plots in a row
    #Derivatives
    #define subplots 
    f, subplts = plt.subplots(1,6) #first number is the rows and the second number is the columns.

    #Calculating derivative of a galaxy with respect to each of each parameters. 
    Dgal_sigma = partialDeriveGal(gal_image, orig_params, 'gal_sigma', orig_params['gal_sigma'].value * .0001)
    Dgal_flux = partialDeriveGal(gal_image, orig_params, 'gal_flux', orig_params['gal_flux'].value * .0001)
    De1 = partialDeriveGal(gal_image, orig_params, 'e1', 1 * .0001)
    De2 = partialDeriveGal(gal_image, orig_params, 'e2', 1 * .0001)
    Dx0 = partialDeriveGal(gal_image, orig_params, 'x0', nx * .0001)
    Dy0 = partialDeriveGal(gal_image, orig_params, 'y0', ny * .0001)

    #list with derivatives and list with their names. 
    names = ['Dx0', 'Dy0','Dgal_sigma', 'Dgal_flux', 'De1', 'De2' ] 
    Ds = [Dgal_sigma, Dgal_flux, De1, De2, Dx0, Dy0]

    for i in range(len(Ds)):
        drawPlot(subplts[i], Ds[i], names[i])



    # #labels and titles 
    # plt.xlabel('Smarts')
    # plt.ylabel('Something')
    # plt.title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')


    #Find fihser matrix. 
    

    #or with the general definition of fisher matrix integrating (.sum()) over all pixels and put them in a numpy array. 
    #initialize fisher matrix as a numpy array of only zero s
    FisherM = np.zeros((6,6))

    #this is the jiggle in one pixel due to the noise, we define it to be 1 for now, but we can rescale later
    sigma_n  = 1

    for i in range(len(Ds)): 
        for j in range(len(Ds)):
            FisherM[i][j] = (Ds[i] * Ds[j]).sum() /(sigma_n**2)

    #after havine one derivative for each down... 
        #do in a triangular plot fashion like in the presentation 
        #do not forget about multiplying to obtain fisher matrix eleemnts

    # f = plt.figure()
    # for i in range(len(Ds)): 
    #     for j in range(len(Ds)):
    #         if(i >= j):
    #             ax = f.add_subplot(6,6, 6 * i + j + 1)
    #             ax.set_xlabel('FishM' + '(' + str(i) + ',' + str(j) + ') =' + str(round(FisherM[i][j],5)))
    #             drawPlot(ax, Ds[i] * Ds[j], names[i] + ',' + names[j]) 


    f = plt.figure()
    for i in range(len(Ds)): 
        for j in range(len(Ds)):
            if(i >= j):
                ax = f.add_subplot(6,6, 6 * i + j + 1)
                if(i == 5):
                    ax.set_xlabel(names[j])
                if(j == 0):
                    ax.set_ylabel(names[i])
                drawPlot(ax, Ds[i] * Ds[j], "") 



    # #or with chi2
    # FisherM_chi2 = np.zeros((6,6))







    #do second derivatives to calculate bias matrix 

    #want to check answers analitically with the formulas from the paper.


    #variables for adjusting space between plots. 

    left  = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = -.6  # the amount of width reserved for blank space between subplots
    hspace = .05 # the amount of height reserved for white space between subplots

    #command for adjusting space between plots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)


    #display the plots
    plt.show()



if __name__ == "__main__":
    main(sys.argv)

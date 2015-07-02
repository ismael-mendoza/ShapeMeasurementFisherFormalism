"""We want to draw a fisher matrix for a resulting galaxy and its parameters. We do everything with a Gaussian in this case."""

import sys
import os
import math
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


#variables for adjusting space between plots. 

left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = .4   # the amount of width reserved for blank space between subplots
hspace = 1.0   # the amount of height reserved for white space between subplots




def drawPlot(subplt, plot, name): 
    """draws a given plot with default values using imshow in the given subplot
    if they are not subplots should just pass plt. It also adds a title with name. -- OWN"""

    #Remove x axis but leave possibility for a x labe. 
    subplt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are off

    #remove y axis completely including label. 
    subplt.get_yaxis().set_visible(False)
    subplt.set_title(name)

    #be careful none must be in quotes. 
    #vmax and vmin useful to center around zero which looks nice. 
    subplt.imshow(plot, interpolation = 'None', cmap=cm.RdYlGn,
            origin='lower', extent=[-nx,nx,-ny,ny],
            vmax=abs(plot).max(), vmin=-abs(plot).max())

def partialDeriveGal(gal_image, params, parameter, step):
    """Calculates the partial derivative of a galaxy with respect to one of its parameters.
    Uses average of adding step and subtracting step (for some reason this is CRITICAL)
    Returns a galaxy image array not a number or an element of the fisher matrix. -- OWN"""
   
   #used to change each parameter slightly and draw an image with that. 
    modParams = params

    parameter_initValue = params[parameter].value

    modParams.add(parameter, value = parameter_initValue + step)
    gal_imageUp = drawGalaxy(modParams)

    modParams.add(parameter, value = parameter_initValue - step)
    gal_imageDown = drawGalaxy(modParams)

    return (gal_imageUp - gal_imageDown).array / (2 * step)

def secondPartialDeriveGal(gal_image, params, parameter1, parameter2, step1, step2):
   
    modParams2 = params

    parameter2_initValue = params[parameter2].value 

    modParams2.add(parameter2, value = parameter2_initValue + step2)
    Dgal_imageUp = partialDeriveGal(gal_image, modParams2, parameter1, step1)

    modParams2.add(parameter2, value = parameter2_initValue - step2)
    Dgal_imageUp = partialDeriveGal(gal_image, modParams2, parameter1, step1)

    return ( Dgal_imageUp - Dgal_imageUp) / (2 * step2)

def partialDeriveChi2(gal_image, params, parameter, step, sigma_n):
    """partial derive chi2 with respect to a parameter evaluated on params. --OWN """
    modParams = params

    parameter_initValue = params[parameter].value

    modParams.add(parameter, value = parameter_initValue + step)
    chi2_up = chi2(gal_image, modParams, sigma_n)

    modParams.add(parameter, value = parameter_initValue - step)
    chi2_down = chi2(gal_image, modParams, sigma_n)

    return (chi2_up - chi2_down)/ (2 * step)

def secondPartialDeriveChi2(gal_image, params, parameter1, parameter2, step1, step2, sigma_n):

    modParams2 = params

    parameter2_initValue = params[parameter2].value 

    modParams2.add(parameter2, value = parameter2_initValue + step2)
    Dchi2_up = partialDeriveChi2(gal_image, modParams2, parameter1, step1, sigma_n)

    modParams2.add(parameter2, value = parameter2_initValue - step2)
    Dchi2_down = partialDeriveChi2(gal_image, modParams2, parameter1, step1, sigma_n)

    return (Dchi2_up - Dchi2_down) / (2 * step2)


def chi2(gal_image, params, sigma_n): 
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

    #variance_noise = image.addNoiseSNR(galsim.GaussianNoise(), 50, True) 

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

    orig_params = Parameters()
    orig_params.add('gal_sigma', value = 1.) # arcsec 
    orig_params.add('gal_flux', value = 100.)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    orig_params.add('e1', value = .1) #ellipticity: e1 
    orig_params.add('e2', value = -.5)#ellipticity: e2
    orig_params.add('x0', value = 0.) #shift in x origin. 
    orig_params.add('y0', value = 0.)     #shift in y

    #get image of the original galaxy
    gal_image = drawGalaxy(orig_params)

    #can draw multiple images at once separately by calling subplots multiple times.
    figure1, subplt = plt.subplots(1,1)
    drawPlot(subplt, gal_image.array, 'Initial Galaxy')

  
    #first we want only derivatives of the galaxy with respect to each parameter numerically so 6 plots in a row
    #Derivatives
    #define subplots 
    figure2, subplts = plt.subplots(1,6) #first number is the rows and the second number is the columns.

    #Calculating derivative of a galaxy with respect to each of each parameters.

    steps = Parameters() 
    steps.add('gal_sigma', value = orig_params['gal_sigma'].value * .0001)
    steps.add('gal_flux', value = orig_params['gal_flux'].value * .0001)
    steps.add('e1', value = 1*.0001) 
    steps.add('e2', value = 1*.0001)
    steps.add('x0', value = nx* .0001)
    steps.add('y0', value = ny* .0001)

    Dgal_sigma = partialDeriveGal(gal_image, orig_params, 'gal_sigma', steps['gal_sigma'].value)
    Dgal_flux = partialDeriveGal(gal_image, orig_params, 'gal_flux', steps['gal_flux'].value)
    De1 = partialDeriveGal(gal_image, orig_params, 'e1', steps['e1'].value)
    De2 = partialDeriveGal(gal_image, orig_params, 'e2', steps['e2'].value)
    Dx0 = partialDeriveGal(gal_image, orig_params, 'x0', steps['x0'].value)
    Dy0 = partialDeriveGal(gal_image, orig_params, 'y0', steps['y0'].value)

    #list with derivatives and list with their names_Ds_gal. 
    names_Ds_gal = []
    for name in orig_params.keys(): 
        names_Ds_gal.append('D' + name)

    Ds_gal = [Dgal_sigma, Dgal_flux, De1, De2, Dx0, Dy0]

    for i in range(len(Ds_gal)):
        drawPlot(subplts[i], Ds_gal[i], names_Ds_gal[i])



    # #labels and titles 
    # plt.xlabel('Smarts')
    # plt.ylabel('Something')
    # plt.title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')


    #Find fihser matrix. 
    

    #or with the general definition of fisher matrix integrating (.sum()) over all pixels and put them in a numpy array. 
    #initialize fisher matrix as a numpy array of only zeros
    #do not forget about multiplying derivatives to obtain fisher matrix eleemnts

    FisherM = np.zeros((6,6))

    #this is the jiggle in one pixel due to the noise, we define it to be 1 for now, but we can rescale later
    sigma_n  = 1

    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)):
            FisherM[i][j] = (Ds_gal[i] * Ds_gal[j]).sum() /(sigma_n**2)


    #print FisherM - FisherM.T ##it is symmetric!! 

    #or with chi2
    #but need to average over many galaxies with different noises... 
    FisherM_chi2 = np.zeros((6,6))

    for i in range(len(orig_params.keys())):
        for j in range(len(orig_params.keys())): 
            FisherM_chi2[i][j] = .5 * secondPartialDeriveChi2(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value, sigma_n)

    # #just to check. if fisherchi2 works. 
    # for i in range(len(orig_params.keys())):
    #     for j in range(len(orig_params.keys())): 
    #         print orig_params.keys()[i] + ',' + orig_params.keys()[j] + ': ' + str(FisherM_chi2[i][j])



    #do in a triangular plot fashion like in the presentation to show plots and fisher matrix elements.

    figure3 = plt.figure()
    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)):
            if(i >= j):
                ax = figure3.add_subplot(6,6, 6 * i + j + 1)
                ax.set_xlabel(str(round(FisherM[i][j],5)))
                drawPlot(ax, Ds_gal[i] * Ds_gal[j], 'G:' + names_Ds_gal[i] + '*' + names_Ds_gal[j]) 

    figure4 = plt.figure()
    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)):
            if(i >= j):
                ax = figure4.add_subplot(6,6, 6 * i + j + 1)
                ax.set_xlabel(str(round(FisherM_chi2[i][j],5)))
                drawPlot(ax, Ds_gal[i] * Ds_gal[j], 'chi2:D' + orig_params.keys()[i] + 'D' + orig_params.keys()[j]) 





    #do second derivatives to calculate bias matrix 

    #SecondDs_gal = np.zeros((6,6))

    for i in range(len(orig_params.keys())):
        for j in range(len(orig_params.keys())): 
            #print secondPartialDeriveGal(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value)
            print secondPartialDeriveGal(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value)

    # # figure5 = plt.figure()
    # for i in range(len(Ds_gal)): 
    #     for j in range(len(Ds_gal)):
    #         if(i >= j):
    #             ax = figure5.add_subplot(6,6, 6 * i + j + 1)
    #             ax.set_xlabel(str(round(FisherM_chi2[i][j],5)))
    #             drawPlot(ax, SecondDs_gal[i][j], 'D' + orig_params.keys()[i] + 'D' + orig_params.keys()[j]) 


#     figure5 = plt.figure()
#     for i in range(len(Ds_gal)): 
#         for j in range(len(Ds_gal)):
#             ax = figure5.add_subplot(6,6, i + j + 1)
#             drawPlot(ax, secondPartialDeriveGal(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value)
# , 'G:D' + orig_params.keys()[i] + 'D' + orig_params.keys()[j]) 


        #want to check answers analitically with the formulas from the paper.


    #command for adjusting space between plots
    figure3.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    figure4.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)


    #display the plots
    figure1.show()
    figure3.show()
    figure4.show()
    figure5.show()



if __name__ == "__main__":
    main(sys.argv)

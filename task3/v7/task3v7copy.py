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
from constantsTask3 import *
from functionsTask3 import *



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
    #initialize dictionary (for ease of use) with parameters.

    orig_params = dict()
    orig_params['gal_sigma'] =  1.  #0; arcsec
    orig_params['gal_flux'] = 100.  #1 ; total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    orig_params['e1'] = 0.          #2 ; ellipticity: e1 
    orig_params['e2'] = 0.         #3; ellipticity: e2
    orig_params['x0'] = 0.          #4;shift in x origin. 
    orig_params['y0'] = 0.          #5;shift in y


    #get image of the original galaxy
    gal_image = drawGalaxy(orig_params)

    #can draw multiple images at once separately by calling subplots multiple times.
    #draws a plot of the galaxy.
    figure1, subplt = plt.subplots(1,1)
    drawPlot(subplt, gal_image.array, 'Initial Galaxy')

  
    #first we want only derivatives of the galaxy with respect to each parameter numerically so 6 plots in a row
    #Derivatives
    #define subplots for the derviatives
    figure2, subplts = plt.subplots(1,6) #first number is the rows and the second number is the columns.

    #Calculating derivative of a galaxy with respect to each of each parameters.
    #define the steps for derivatives of each individual parameter.
    steps = dict() 
    steps['gal_sigma'] = orig_params['gal_sigma'] * .01
    steps['gal_flux']  = orig_params['gal_flux'] * .01
    steps['e1'] = 1*.01
    steps['e2'] = 1*.01
    steps['x0'] = nx* .01
    steps['y0'] = ny* .01


    #create a vector(numpy array with 1 row and 6 columns) with the derivatives of the model with respect to each parameter
    #doubt if more convenient than a dictionary but remember the numbers from the first list and maintain the order because vector is more natural.
    Ds_gal = []

    for parameter in orig_params.keys():
        Ds_gal.append(partialDerive(func = drawGalaxy, params = orig_params, parameter =parameter, step = steps[parameter]).array)

    # #list with derivatives and list with their names_Ds_gal. 
    names_Ds_gal = []
    for name in orig_params.keys(): 
         names_Ds_gal.append('D' + name)

    #Ds_gal = [Dgal_sigma, Dgal_flux, De1, De2, Dx0, Dy0]

    for i in range(len(Ds_gal)):
        drawPlot(subplts[i], Ds_gal[i], 'D' + orig_params.keys()[i])


#     #Find fihser matrix. 
#     #or with the general definition of fisher matrix integrating (.sum()) over all pixels and put them in a numpy array. 
#     #initialize fisher matrix as a numpy array of only zeros
#     #do not forget about multiplying derivatives to obtain fisher matrix eleemnts

    FisherM = np.zeros((6,6))

    #this is the jiggle in one pixel due to the noise, we define it to be 1 for now, but we can rescale later
    sigma_n  = 1

    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)):
            FisherM[i][j] = (Ds_gal[i] * Ds_gal[j]).sum() /(sigma_n**2)



    #print FisherM - FisherM.T ##it is symmetric!! 

    #or with chi2
    # #but need to expected value (average over many galaxies with different noises)... 
    # FisherM_chi2 = np.zeros((6,6))

    # for i in range(len(orig_params.keys())):
    #     for j in range(len(orig_params.keys())): 
    #         FisherM_chi2[i][j] = .5 * secondPartialDeriveChi2(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value, sigma_n)

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

    # figure4 = plt.figure()
    # for i in range(len(Ds_gal)): 
    #     for j in range(len(Ds_gal)):
    #         if(i >= j):
    #             ax = figure4.add_subplot(6,6, 6 * i + j + 1)
    #             ax.set_xlabel(str(round(FisherM_chi2[i][j],5)))
    #             drawPlot(ax, Ds_gal[i] * Ds_gal[j], 'chi2:D' + orig_params.keys()[i] + 'D' + orig_params.keys()[j]) 


    #do second derivatives to calculate bias matrix later. 

#     SecondDs_gal = np.zeros((6,6))

#     for i in range(len(orig_params.keys())):
#         for j in range(len(orig_params.keys())): 
#             #print secondPartialDeriveGal(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value)
#             print secondPartialDeriveGal(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value)

#     # figure5 = plt.figure()
#     for i in range(len(Ds_gal)): 
#         for j in range(len(Ds_gal)):
#             if(i >= j):
#                 ax = figure5.add_subplot(6,6, 6 * i + j + 1)
#                 ax.set_xlabel(str(round(FisherM_chi2[i][j],5)))
#                 drawPlot(ax, SecondDs_gal[i][j], 'D' + orig_params.keys()[i] + 'D' + orig_params.keys()[j]) 


#     figure5 = plt.figure()
#     for i in range(len(Ds_gal)): 
#         for j in range(len(Ds_gal)):
#             ax = figure5.add_subplot(6,6, i + j + 1)
#             drawPlot(ax, secondPartialDeriveGal(gal_image, orig_params, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]].value, steps[orig_params.keys()[j]].value)
# , 'G:D' + orig_params.keys()[i] + 'D' + orig_params.keys()[j]) 


    #want to check answers analitically with the formulas from the paper.


    #command for adjusting space between plots, constants defined in constantsTask3.py
    figure3.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    # figure4.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)


    #display the plots
    figure1.show()
    figure2.show()
    figure3.savefig('figure.png')
    #figure4.show()
    # figure5.show()



if __name__ == "__main__":
    main(sys.argv)

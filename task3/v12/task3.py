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
import matplotlib
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

    if(len(argv) == 1): 
        argv.append('')

    #first we want to set here the true values. Or the ones we want our galaxy to have.    
    #possible parameters for the galaxy formation. 
    #initialize dictionary (for ease of use) with parameters.

    orig_params = dict()
    orig_params['gal_sigma'] =  3.  #0; arcsec
    orig_params['gal_flux'] = 100.  #1 ; total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    orig_params['e1'] = .1          #2 ; ellipticity: e1 
    orig_params['e2'] = -.5         #3; ellipticity: e2
    orig_params['x0'] = 0.          #4;shift in x origin. 
    orig_params['y0'] = 0.          #5;shift in y


    #get image of the original galaxy
    gal_image = drawGalaxy(orig_params)

    #first we want only derivatives of the galaxy with respect to each parameter numerically so 6 plots in a row
    #Derivatives

    #define the steps for derivatives of each individual parameter.
    steps = dict() 
    steps['gal_sigma'] = orig_params['gal_sigma'] * .0001
    steps['gal_flux']  = orig_params['gal_flux'] * .0001
    steps['e1'] = 1*.0001
    steps['e2'] = 1*.0001
    steps['x0'] = nx* .0001
    steps['y0'] = ny* .0001


    #create a vector(numpy array with 1 row and 6 columns) with the derivatives of the model with respect to each parameter
    #doubt if more convenient than a dictionary but remember the numbers from the first list and maintain the order because vector is more natural.
    Ds_gal = []

    for parameter in orig_params.keys():
        Ds_gal.append(partialDifferentiate(func = drawGalaxy, parameter =parameter, step = steps[parameter])(orig_params).array)

    #Ds_gal = [partialDifferentiate(func = drawGalaxy, parameter =parameter, step = steps[parameter])(orig_params).array for parameter in orig_params.keys() ]

    #Find fihser matrix. 
    #or with the general definition of fisher matrix integrating (.sum()) over all pixels and put them in a numpy array. 
    #obtain fisher matrix eleemnts by multiplying derivatives.

    #this is the jiggle in one pixel due to the noise, we define it to be 1 for now, but we can rescale later, uniform for all pixels for now too.
    sigma_n  = 1
    FisherM = np.zeros([6,6])

    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)): 
            FisherM[i][j] = (Ds_gal[i] * Ds_gal[j]).sum() /(sigma_n**2)

    #FisherM = np.array([[(Ds_gal[i] * Ds_gal[j]).sum() /(sigma_n**2) for j in range(len(Ds_gal))] for i in range(len(Ds_gal))]) ##list comprehension is just too convoluted but was good practice. 

    #print FisherM - FisherM.T ##it is symmetric!! 

    #or with chi2
    #but need to expected value (average over many galaxies with different noises)... although if the galaxy is noiseless, burchat argues it is the same or very close.
    FisherM_chi2 = np.zeros([6,6])
    for i in range(len(orig_params.keys())): 
        for j in range(len(orig_params.keys())): 
            FisherM_chi2[i][j] = .5 * secondPartialDifferentiate(chi2, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]], steps[orig_params.keys()[j]], sigma_n = sigma_n, gal_image = gal_image)(orig_params)

    #FisherM_chi2 = np.array([[.5 * secondPartialDifferentiate(chi2, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]], steps[orig_params.keys()[j]], sigma_n = sigma_n, gal_image = gal_image)(orig_params) for j in range(len(orig_params.keys()))] for i in range(len(orig_params.keys()))])

    #print FisherM_chi2 - FisherM_chi2.T ##it is symmetric


    #do second derivatives of galaxies. Bias matrix can be later. 
    SecondDs_gal = [[] for i in range(6)]

    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)):
            SecondDs_gal[i].append(secondPartialDifferentiate(drawGalaxy, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]], steps[orig_params.keys()[j]])(orig_params).array)

    #bias matrix. 



#     #want to check answers analitically with the formulas from the paper.

    

    if(argv[1] == 'plot'):

        #can draw multiple images at once separately by calling subplots multiple times.
        #draws a plot of the galaxy.
        figure1, subplt= plt.subplots(1,1)
        figure1.suptitle('Initial Galaxy', fontsize = 20)
        drawPlot(subplt, gal_image.array)


        #figure with derivatives. 
        figure2 = plt.figure() 
        figure2.suptitle('Derivatives of model with respect to each parameter', fontsize = 20)
        for i in range(len(Ds_gal)):
            ax = figure2.add_subplot(1,6,i+1)
            drawPlot(ax, Ds_gal[i], title = 'D' + orig_params.keys()[i])



        #do in a triangular plot fashion like in the presentation to show plots and fisher matrix elements.
        figure3 = plt.figure()
        figure3.suptitle('Fisher matrix elements ', fontsize=14, fontweight='bold')
        for i in range(len(Ds_gal)): 
            for j in range(len(Ds_gal)):
                if(i >= j):
                    ax = figure3.add_subplot(6,6, 6 * i + j + 1)
                    drawPlot(ax, Ds_gal[i] * Ds_gal[j])
                    if(j == 0):
                        ax.set_ylabel('D' + orig_params.keys()[i] )
                    if(i == len(Ds_gal) - 1):
                        ax.set_xlabel('D' + orig_params.keys()[j])


        #now want with values of fisher matrix in the middle of the plot.
        figure4 = plt.figure()
        figure4.suptitle('Fisher matrix elements with values ', fontsize=14, fontweight='bold')
        for i in range(len(Ds_gal)): 
            for j in range(len(Ds_gal)):
                if(i >= j):
                    ax = figure4.add_subplot(6,6, 6 * i + j + 1)
                    drawPlot(ax, Ds_gal[i] * Ds_gal[j])
                    if(j == 0):
                        ax.set_ylabel('D' + orig_params.keys()[i] )
                    if(i == len(Ds_gal) - 1):
                        ax.set_xlabel('D' + orig_params.keys()[j])
                    ax.text(-20,0, str(round(FisherM[i][j],5)))

        figure5 = plt.figure()
        figure5.suptitle('Fisher matrix elements with values of chi2 method ', fontsize=14, fontweight='bold')
        for i in range(len(Ds_gal)): 
            for j in range(len(Ds_gal)):
                if(i >= j):
                    ax = figure5.add_subplot(6,6, 6 * i + j + 1)
                    drawPlot(ax, Ds_gal[i] * Ds_gal[j])
                    if(j == 0):
                        ax.set_ylabel('D' + orig_params.keys()[i] )
                    if(i == len(Ds_gal) - 1):
                        ax.set_xlabel('D' + orig_params.keys()[j])
                    ax.text(-20,0, str(round(FisherM_chi2[i][j],5)))

        figure6 = plt.figure()
        figure6.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=14, fontweight='bold')
        for i in range(len(Ds_gal)): 
            for j in range(len(Ds_gal)):
                if(i >= j):
                    ax = figure6.add_subplot(6,6, 6 * i + j + 1)
                    drawPlot(ax, SecondDs_gal[i][j]) 
                    if(j == 0):
                        ax.set_ylabel('D' + orig_params.keys()[i] )
                    if(i == len(Ds_gal) - 1):
                        ax.set_xlabel('D' + orig_params.keys()[j])


        #display the plots
        figure1.show()
        figure2.show()
        figure3.show()
        figure4.show()
        figure5.show()
        figure6.show()



    if(argv[1] == 'clear'):
        #clear memory of figures. did not work inside a function.
        figures =[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        for figure in figures:    
           figure.clf()
        plt.close()
        gc.collect()




if __name__ == "__main__":
    main(sys.argv)

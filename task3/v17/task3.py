"""We want to draw a fisher matrix for a resulting galaxy and its parameters. We do everything with a 6 paramater-Gaussian in this case."""

import sys
import os
import inspect
import subprocess
import math
import galsim
import subprocess
from lmfit import minimize, Parameters, Parameter, fit_report
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib
from constantsTask3 import *
from functionsTask3 import *



def main(argv):

    if(len(argv) == 1): 
        argv.append('')
    if(len(argv) == 2): 
        argv.append('0')

    #first we want to set here the true values. Or the ones we want our galaxy to have.    
    #possible parameters for the galaxy formation. 
    #initialize dictionary (for ease of use) with parameters.

    orig_params = dict()
    orig_params['gal_flux'] = 100.  #1 ; total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    orig_params['gal_sigma'] =  3.  #0; arcsec
    orig_params['e1'] = .1          #2 ; ellipticity: e1 
    orig_params['e2'] = -.5         #3; ellipticity: e2
    orig_params['x0'] = 0.          #4;shift in x origin. 
    orig_params['y0'] = 0.          #5;shift in y


    #get image of the original galaxy
    gal_image = drawGalaxy(orig_params)

    #first we want only derivatives of the galaxy with respect to each parameter numerically so 6 plots in a row

    #define the steps for derivatives of each individual parameter.
    steps = dict() 
    steps['gal_flux']  = orig_params['gal_flux'] * .0001
    steps['gal_sigma'] = orig_params['gal_sigma'] * .0001
    steps['e1'] = 1*.0001
    steps['e2'] = 1*.0001
    steps['x0'] = nx* .0001
    steps['y0'] = ny* .0001


    #create a vector(numpy array with 1 row and 6 columns) with the derivatives of the model with respect to each parameter
    #doubt if more convenient than a dictionary but remember the numbers from the first list and maintain the order because vector is more natural.
    Ds_gal = []

    for parameter in orig_params.keys():
        Ds_gal.append(partialDifferentiate(func = drawGalaxy, parameter =parameter, step = steps[parameter])(orig_params).array)

    #Find fisher matrix. 
    #or with the general definition of fisher matrix integrating (.sum()) over all pixels and put them in a numpy array. 
    #obtain fisher matrix eleemnts by multiplying derivatives.

    #this is the jiggle in one pixel due to the noise, we define it to be 1 for now, but we can rescale later, uniform for all pixels for now too.
    sigma_n  = 1

    FisherM = np.zeros([6,6])
    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)): 
            FisherM[i][j] = (Ds_gal[i] * Ds_gal[j]).sum() /(sigma_n**2)


    #or with chi2
    #but need to expected value (average over many galaxies with different noises)... although if the galaxy is noiseless, burchat argues it is the same or very close.
    FisherM_chi2 = np.zeros([6,6])
    for i in range(len(orig_params.keys())): 
        for j in range(len(orig_params.keys())): 
            FisherM_chi2[i][j] = .5 * secondPartialDifferentiate(chi2, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]], steps[orig_params.keys()[j]], sigma_n = sigma_n, gal_image = gal_image)(orig_params)



    #do second derivatives of galaxies. Bias matrix can be later. 
    SecondDs_gal = [[] for i in range(6)]

    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)):
            SecondDs_gal[i].append(secondPartialDifferentiate(drawGalaxy, orig_params.keys()[i], orig_params.keys()[j], steps[orig_params.keys()[i]], steps[orig_params.keys()[j]])(orig_params).array)

    #bias matrix.
    BiasM_images = [[[] for j in range(6)] for i in range(6)]
    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)): 
            for k in range(len(Ds_gal)): 
                BiasM_images[i][j].append((Ds_gal[i] * SecondDs_gal[j][k]) / (sigma_n**2))

    #summing each element over all pixels get the numerical values of each element of the bias matrix. 
    BiasM = [[[BiasM_images[i][j][k].sum() for k in range(6)] for j in range(6)] for i in range(6)]


    #Covariance matrix is inverse of Fisher Matrix: 
    CovM = np.linalg.inv(FisherM)

    #print np.inner(FisherM,CovM)  ##they are inverses!! 

    #now we want bias of each parameter per pixel, so we can see how each parameter contributes. 
    bias_images = []
    for i in range(6):
        sumation = 0 
        for j in range(6):
            for k in range(6):
                for l in range(6):
                    sumation += CovM[i][j]*CovM[k][l]*BiasM_images[j][k][l]
        bias_images.append((-.5) *sumation)


    biases = [image.sum() for image in bias_images]


    print biases[0]
    print CovM[0][0]
    print (2 * biases[0] * orig_params['gal_flux'] / CovM[0][0])



    #use lmfit over a lot of noisy images. 
    #want to check answers analitically with the formulas from the paper.

    




########################################################
#figures. 

# plot: #display plots. 
# 1 ##shows the initial galaxy
# 2 ##shows first derivatives of the galaxy 
# 3 ##shows Fisher matrix elements of the galaxy
# 4 ##shows Fisher matrix elements with values summed over each pixel. 
# 5 ##shows Fisher matrix elements with values of the chi2 method for calculating fisher matrix elements. 
# 6 ##shows covariance matrix elements (just numbers)
# 7 ##shows plots of the second derivatives of the initial galaxy with respect to each parameter. 
# 8 #figures shows 6 plots for each i of the bias tensor. 
# 9 #figures shows 6 plots for each i of the bias tensor with the sum over pixels of each of them.    

    if(argv[1] == 'plot'):

        plt.rc('text', usetex=False)
        list_of_plots = np.array([int(n) for n in argv[2:]])

        #here we slice arg[2:] because we have to have a numpy array only with numbers for this to work and only the arguments after 'plot' are numbers.
        if((list_of_plots == 1).any()): 
            figure1, subplt= plt.subplots(1,1)
            figure1.suptitle('Initial Galaxy', fontsize = 20)
            drawPlot(subplt, gal_image.array)
            SaveFigureToPdfAndOpen(figure1, 'figure1.png') #this will open for preview because it is the default defined in my mac


        if((list_of_plots == 2).any()):  
            figure2 = plt.figure() 
            figure2.suptitle('Derivatives of model with respect to each parameter', fontsize = 20)
            for i in range(len(Ds_gal)):
                ax = figure2.add_subplot(1,6,i+1)
                drawPlot(ax, Ds_gal[i], title = 'D' + orig_params.keys()[i])

            SaveFigureToPdfAndOpen(figure2, 'figure2.png')


        if((list_of_plots == 3).any()):
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
            SaveFigureToPdfAndOpen(figure3, 'figure3.png')


        if((list_of_plots == 4).any()):
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
                        ax.text(-20,0, str(round(FisherM[i][j],5)), fontweight='bold')
            SaveFigureToPdfAndOpen(figure4, 'figure4.png')


        if((list_of_plots == 5).any()):
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
                        ax.text(-20,0, str(round(FisherM_chi2[i][j],5)), fontweight='bold')

            SaveFigureToPdfAndOpen(figure5, 'figure5.png')

        if((list_of_plots == 6).any()):
            figure6 = plt.figure()
            figure6.suptitle('Covariance matrix elements', fontsize=14, fontweight='bold')
            for i in range(6): 
                for j in range(6):
                    if(i >= j):
                        ax = figure6.add_subplot(6,6, 6 * i + j + 1)
                        drawPlot(ax, Ds_gal[i] * Ds_gal[j] * 0)
                        if(j == 0):
                            ax.set_ylabel(orig_params.keys()[i] )
                        if(i == 6 - 1):
                            ax.set_xlabel(orig_params.keys()[j])
                        ax.text(-20,0, str(round(CovM[i][j],5)), fontweight='bold')

            SaveFigureToPdfAndOpen(figure6, 'figure6.png')

        if((list_of_plots == 7).any()):
            figure7 = plt.figure()
            figure7.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=14, fontweight='bold')
            for i in range(len(Ds_gal)): 
                for j in range(len(Ds_gal)):
                    if(i >= j):
                        ax = figure7.add_subplot(6,6, 6 * i + j + 1)
                        drawPlot(ax, SecondDs_gal[i][j]) 
                        if(j == 0):
                            ax.set_ylabel('D' + orig_params.keys()[i] )
                        if(i == len(Ds_gal) - 1):
                            ax.set_xlabel('D' + orig_params.keys()[j])

            SaveFigureToPdfAndOpen(figure7, 'figure7.png')

        if((list_of_plots == 8).any()):
            figuresOfBiasMatrix = [] 
            for i in range(len(Ds_gal)):
                figure = plt.figure()
                figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + orig_params.keys()[i] + ')' , fontsize=14, fontweight='bold')
                for j in range(len(Ds_gal)): 
                    for k in range(len(Ds_gal)):
                        if(j >= k):
                            ax = figure.add_subplot(6,6, 6 * j + k + 1)
                            drawPlot(ax, BiasM_images[i][j][k]) 
                            if(k == 0):
                                ax.set_ylabel('D' + orig_params.keys()[j])
                            if(j == len(Ds_gal) - 1):
                                ax.set_xlabel('D' + orig_params.keys()[k])
                figuresOfBiasMatrix.append(figure)

            for i, figure in enumerate(figuresOfBiasMatrix): 
                SaveFigureToPdfAndOpen(figure, 'figure' + str(8) + '_' + str(i) + '.png') 

        if((list_of_plots == 9).any()):
            figuresOfBiasMatrixNumbers = []    
            for i in range(len(Ds_gal)):
                figure = plt.figure()
                figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + orig_params.keys()[i] + ')' , fontsize=14, fontweight='bold')
                for j in range(len(Ds_gal)): 
                    for k in range(len(Ds_gal)):
                        if(j >= k):
                            ax = figure.add_subplot(6,6, 6 * j + k + 1)
                            drawPlot(ax, BiasM_images[i][j][k]) 
                            if(k == 0):
                                ax.set_ylabel('D' + orig_params.keys()[j] )
                            if(j == len(Ds_gal) - 1):
                                ax.set_xlabel('D' + orig_params.keys()[k])
                            ax.text(-20,0, str(round(BiasM[i][j][k],5)), fontweight='bold')
                figuresOfBiasMatrixNumbers.append(figure)

            for i, figure in enumerate(figuresOfBiasMatrixNumbers): 
                SaveFigureToPdfAndOpen(figure, 'figure' + str(9) + '_' + str(i) + '.png') 

        if((list_of_plots == 10).any()):
	    	figure10 = plt.figure() 
	        figure10.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 15)
	        for i in range(len(Ds_gal)):
	        	ax = figure10.add_subplot(1,6,i+1)
	        	drawPlot(ax, bias_images[i], title = orig_params.keys()[i])

	        SaveFigureToPdfAndOpen(figure10, 'figure10.png')

        if((list_of_plots == 11).any()):
	    	figure11 = plt.figure() 
	        figure11.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 14)
	        for i in range(len(Ds_gal)):
	        	ax = figure11.add_subplot(1,6,i+1)
	        	drawPlot(ax, bias_images[i], title = orig_params.keys()[i])
	        	ax.text(-20,0, str(round(biases[i],5)), fontweight='bold')

	        SaveFigureToPdfAndOpen(figure11, 'figure11.png')



if __name__ == "__main__":
    main(sys.argv)

"""We want to draw a fisher matrix for a resulting galaxy and its parameters. We do everything with a 6 paramater-Gaussian in this case.

v18/19
*transformed all lists and matrices to dictionaries to avoid mistakes in indices/names.
*going to change parameters to the ones in the paper to compare the biases for sanity check.
*added (the much wanted) constant number of parameters.
*fractional bias instead of bias.

*investigate if the granular elements are actually 0 by rounding to the 3 or 2nd significant digit.
*print std of the arrays in the image to see how close they are to 0.

v20
*new parameters from refrigier paper. 
*made modifications to drawGalaxy function so that it can also take q,beta as parameters.
*added variance function and way to calculate other parameters not used initially in the model 

v21
*change to initial parameters, as noticed that not necessary to use other ones. 
*want to add a bias function which does the same as variance function. 

v22 
*so now results are consistent with refrigier's which means that my method is on its way to being correct. Need to average bias by residuals in different noisy images to be sure and compare. 
*Also want to experiment convolving a galaxy with a psf and see how it behaves with the fisher formalism.


v23
*probably will make changes so that this formalsim in this file can be applied to galaxy parameters in task3res.py, careful not to add extra functions
or anything if I decide to make a copy of this file, just add lines needed to initialize the other data from .csv file. 

"""


import sys
import os
import inspect
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

    #read parameters from file and initialize dictionary with that.
    with open('initial_data.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        init_params = dict() 
        row = reader.next()
        for parameter_name in row.keys():
            init_params[parameter_name] = row[parameter_name]

    #first we want to set here the true values. Or the ones we want our galaxy to have.    
    #possible parameters for the galaxy formation. 
    #initialize dictionary (for ease of use) with parameters.

    orig_params = dict()
    orig_params['gal_flux'] = init_params['gal_flux']            #total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    orig_params['gal_sigma'] =  init_params['gal_sigma']            #arcsec
    orig_params['e1'] = init_params['e1']                 #ellipticity: e1 
    orig_params['e2'] = init_params['e2']                #ellipticity: e2
    orig_params['x0'] = init_params['x0']                    #shift in x origin. 
    orig_params['y0'] = init_params['y0']                     #shift in y
    #orig_params['q'] = .5                     #minor to major axis ration
    #orig_params['beta'] = 1.75 * np.pi        #angle


    #define psf parameters, 
    psf_params = dict()
    psf_params['psf_flux'] = init_params['psf_flux']
    psf_params['psf_sigma'] = init_params['psf_sigma']




    #names of parameters for dictionaries.
    param_names = orig_params.keys()

    #get image of the original galaxy
    gal_image = drawGalaxy(params = orig_params, psf_params = psf_params)

    #first we want only derivatives of the galaxy with respect to each parameter numerically so 6 plots in a row

    #define the steps for derivatives of each individual parameter.
    steps = dict() 
    steps['gal_flux']  = orig_params['gal_flux'] * .01
    steps['gal_sigma'] = .01
    steps['e1'] = .01
    steps['e2'] = .01
    steps['x0'] = .01
    steps['y0'] = .01
    #steps['q'] = orig_params['q']*.01
    #steps['beta'] = orig_params['beta']*.01



    #create a dictionary with the derivatives of the model with respect to each parameter.
    Ds_gal = {param_names[i]:partialDifferentiate(func = drawGalaxy, parameter = param_names[i], step = steps[param_names[i]], psf_params = psf_params)(orig_params).array for i in range(num_params)}

    #Find fisher matrix. 
    #or with the general definition of fisher matrix integrating (.sum()) over all pixels and put them in a numpy array. 
    #obtain fisher matrix eleemnts by multiplying derivatives.

    #this is the jiggle in one pixel due to the noise, we define it to be 1 for now, but we can rescale later, uniform for all pixels for now too.
    sigma_n  = math.sqrt(init_params['variance_noise'])
    #sigma_n = 1.

    FisherM_images = {}
    for i in range(num_params): 
        for j in range(num_params): 
            FisherM_images[param_names[i],param_names[j]] = (Ds_gal[param_names[i]] * Ds_gal[param_names[j]]) /(sigma_n**2)

    FisherM = {}
    for i in range(num_params): 
        for j in range(num_params): 
            FisherM[param_names[i],param_names[j]] = FisherM_images[param_names[i],param_names[j]].sum() #sum over all pixels. 

    #or with chi2
    #but need to expected value (average over many galaxies with different noises)... although if the galaxy is noiseless, burchat argues it is the same or very close.
    FisherM_chi2 = {}
    for i in range(num_params): 
        for j in range(num_params): 
            FisherM_chi2[param_names[i],param_names[j]] = .5 * secondPartialDifferentiate(chi2, param_names[i], param_names[j], steps[param_names[i]], steps[param_names[j]], sigma_n = sigma_n, gal_image = gal_image, psf_params = psf_params)(orig_params)



    #do second derivatives of galaxies. Bias matrix can be later. 
    SecondDs_gal = {}
    for i in range(num_params): 
        for j in range(num_params):
            SecondDs_gal[param_names[i],param_names[j]] = (secondPartialDifferentiate(drawGalaxy, param_names[i], param_names[j], steps[param_names[i]], steps[param_names[j]], psf_params = psf_params)(orig_params).array)

    #bias matrix.
    BiasM_images = {}
    for i in range(len(Ds_gal)): 
        for j in range(len(Ds_gal)): 
            for k in range(len(Ds_gal)): 
                BiasM_images[param_names[i],param_names[j],param_names[k]] = (Ds_gal[param_names[i]] * SecondDs_gal[param_names[j],param_names[k]]) / (sigma_n**2)

    #summing each element over all pixels get the numerical values of each element of the bias matrix. 
    BiasM = {}

    for i in range(num_params):
        for  j in range(num_params):
            for k in range(num_params):
                BiasM[param_names[i],param_names[j],param_names[k]] = BiasM_images[param_names[i],param_names[j],param_names[k]].sum()


    #Covariance matrix is inverse of Fisher Matrix: 
    CovM = {}
    FisherM_array = np.array([[FisherM[param_names[i],param_names[j]] for i in range(num_params)] for j in range(num_params)])#convert to numpy array because it can be useful.
    CovM_array = np.linalg.inv(FisherM_array)

    for i in range(num_params):
        for j in range(num_params):
            CovM[param_names[i],param_names[j]] = CovM_array[i][j]

    #now we want bias of each parameter per pixel, so we can see how each parameter contributes. 
    bias_images = {}
    for i in range(6):
        sumation = 0 
        for j in range(6):
            for k in range(6):
                for l in range(6):
                    sumation += CovM[param_names[i],param_names[j]]*CovM[param_names[k],param_names[l]]*BiasM_images[param_names[j],param_names[k],param_names[l]]
        bias_images[param_names[i]] = (-.5) * sumation


    biases = {param_names[i]:image.sum() for i,image in enumerate(bias_images.values())}


    # #some sanity checks. 
    # #want to check answers analitically with the formulas from the paper.
    # rhoA =  amplitude_func(orig_params) / math.sqrt(variance(amplitude_func, amplitude_func, orig_params, param_names, CovM, steps))


    # print 'stat. error in a2: '+ str((math.sqrt(variance(a2_func, a2_func, orig_params, param_names, CovM, steps)) / a2_func(orig_params)) *  rhoA)
    # print 'stat. error in a1: '+ str((math.sqrt(variance(a1_func, a1_func, orig_params, param_names, CovM, steps)) / a1_func(orig_params)) *  rhoA)
    # print 'Flux bias: ' + str((biases['gal_flux']/orig_params['gal_flux']) * (rhoA)**2)
    


    #use lmfit over a lot of noisy images, combining parts from task2. 

 

    




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
# 10 ##shows the bias contribution per pixel 
# 11 ##shows the bias contribution per pixel and the sums over all pixels.   

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
            for i in range(num_params):
                ax = figure2.add_subplot(1,6,i+1)
                drawPlot(ax, Ds_gal[param_names[i]], title = 'D' + param_names[i])

            SaveFigureToPdfAndOpen(figure2, 'figure2.png')


        if((list_of_plots == 3).any()):
            figure3 = plt.figure()
            figure3.suptitle('Fisher matrix elements ', fontsize=14, fontweight='bold')
            for i in range(num_params): 
                for j in range(num_params):
                    if(i >= j):
                        ax = figure3.add_subplot(6,6, 6 * i + j + 1)
                        drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                        if(j == 0):
                            ax.set_ylabel('D' + param_names[i] )
                        if(i == len(Ds_gal) - 1):
                            ax.set_xlabel('D' + param_names[j])
            SaveFigureToPdfAndOpen(figure3, 'figure3.png')


        if((list_of_plots == 4).any()):
            figure4 = plt.figure()
            figure4.suptitle('Fisher matrix elements with values ', fontsize=14, fontweight='bold')
            for i in range(num_params): 
                for j in range(num_params):
                    if(i >= j):
                        ax = figure4.add_subplot(6,6, 6 * i + j + 1)
                        drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                        if(j == 0):
                            ax.set_ylabel('D' + param_names[i])
                        if(i == len(Ds_gal) - 1):
                            ax.set_xlabel('D' + param_names[j])
                        ax.text(-20,0, str(round(FisherM[param_names[i],param_names[j]],5)), fontweight='bold')
            SaveFigureToPdfAndOpen(figure4, 'figure4.png')


        if((list_of_plots == 5).any()):
            figure5 = plt.figure()
            figure5.suptitle('Fisher matrix elements with values of chi2 method ', fontsize=14, fontweight='bold')
            for i in range(num_params): 
                for j in range(num_params):
                    if(i >= j):
                        ax = figure5.add_subplot(6,6, 6 * i + j + 1)
                        drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                        if(j == 0):
                            ax.set_ylabel('D' + param_names[i] )
                        if(i == len(Ds_gal) - 1):
                            ax.set_xlabel('D' + param_names[j])
                        ax.text(-20,0, str(round(FisherM_chi2[param_names[i],param_names[j]],5)), fontweight='bold')

            SaveFigureToPdfAndOpen(figure5, 'figure5.png')

        if((list_of_plots == 6).any()):
            figure6 = plt.figure()
            figure6.suptitle('Covariance matrix elements', fontsize=14, fontweight='bold')
            for i in range(num_params): 
                for j in range(num_params):
                    if(i >= j):
                        ax = figure6.add_subplot(6,6, 6 * i + j + 1)
                        drawPlot(ax, FisherM_images[param_names[i],param_names[j]] * 0)
                        if(j == 0):
                            ax.set_ylabel(param_names[i])
                        if(i == 6 - 1):
                            ax.set_xlabel(param_names[j])
                        ax.text(-20,0, str(round(CovM[param_names[i],param_names[j]],5)), fontweight='bold')

            SaveFigureToPdfAndOpen(figure6, 'figure6.png')

        if((list_of_plots == 7).any()):
            figure7 = plt.figure()
            figure7.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=14, fontweight='bold')
            for i in range(num_params): 
                for j in range(num_params):
                    if(i >= j):
                        ax = figure7.add_subplot(6,6, 6 * i + j + 1)
                        drawPlot(ax, SecondDs_gal[param_names[i],param_names[j]])
                        ax.text(0,20, "std:" + str(SecondDs_gal[param_names[i],param_names[j]].std().round(4)), fontsize = 8, fontweight='bold') 
                        if(j == 0):
                            ax.set_ylabel('D' + param_names[i] )
                        if(i == len(Ds_gal) - 1):
                            ax.set_xlabel('D' + param_names[j])

            SaveFigureToPdfAndOpen(figure7, 'figure7.png')

        if((list_of_plots == 8).any()):
            figuresOfBiasMatrix = [] 
            for i in range(num_params):
                figure = plt.figure()
                figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
                for j in range(num_params): 
                    for k in range(num_params):
                        if(j >= k):
                            ax = figure.add_subplot(6,6, 6 * j + k + 1)
                            drawPlot(ax, BiasM_images[param_names[i],param_names[j],param_names[k]]) 
                            if(k == 0):
                                ax.set_ylabel('D' + param_names[j])
                            if(j == len(Ds_gal) - 1):
                                ax.set_xlabel('D' + param_names[k])
                figuresOfBiasMatrix.append(figure)

            for i, figure in enumerate(figuresOfBiasMatrix): 
                SaveFigureToPdfAndOpen(figure, 'figure' + str(8) + '_' + str(i) + '.png') 

        if((list_of_plots == 9).any()):
            figuresOfBiasMatrixNumbers = []    
            for i in range(num_params):
                figure = plt.figure()
                figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
                for j in range(num_params): 
                    for k in range(num_params):
                        if(j >= k):
                            ax = figure.add_subplot(6,6, 6 * j + k + 1)
                            drawPlot(ax, BiasM_images[param_names[i],param_names[j],param_names[k]]) 
                            if(k == 0):
                                ax.set_ylabel('D' + param_names[j] )
                            if(j == len(Ds_gal) - 1):
                                ax.set_xlabel('D' + param_names[k])
                            ax.text(-20,0, str(round(BiasM[param_names[i],param_names[j],param_names[k]],5)), fontweight='bold')
                figuresOfBiasMatrixNumbers.append(figure)

            for i, figure in enumerate(figuresOfBiasMatrixNumbers): 
                SaveFigureToPdfAndOpen(figure, 'figure' + str(9) + '_' + str(i) + '.png') 

        if((list_of_plots == 10).any()):
	    	figure10 = plt.figure() 
	        figure10.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 15)
	        for i in range(num_params):
	        	ax = figure10.add_subplot(1,6,i+1)
	        	drawPlot(ax, bias_images[param_names[i]], title = param_names[i])

	        SaveFigureToPdfAndOpen(figure10, 'figure10.png')

        if((list_of_plots == 11).any()):
	    	figure11 = plt.figure() 
	        figure11.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 14)
	        for i in range(num_params):
	        	ax = figure11.add_subplot(1,6,i+1)
	        	drawPlot(ax, bias_images[param_names[i]], title = param_names[i])
	        	ax.text(-20,0, str(round(biases[param_names[i]] /orig_params[param_names[i]],5)), fontweight='bold')

	        SaveFigureToPdfAndOpen(figure11, 'figure11.png')



if __name__ == "__main__":
    main(sys.argv)

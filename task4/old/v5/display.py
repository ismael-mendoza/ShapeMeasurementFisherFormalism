import argparse
import defaults
import os
import matplotlib.pyplot as plt
from functions import *
import csv
import math
import numpy as np
from fisher import *

def main():
    parser = argparse.ArgumentParser(description = 'Display different results and plots from the galaxies in a given file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--wdir', default = 'output', metavar = 'DIRECTORY', type = str,
    help = 'Specify a directory name where the output will be inputted and files read from.')
    parser.add_argument('--galaxy-file', default = 'galaxies', metavar = 'FILENAME', type = str,
    help = 'Specify a file where galaxies will be registered or read from.')
    parser.add_argument('--plots-dir', default = 'plots', metavar = 'PLOTDIRECTORY', type = str,
    help = 'Specify a directory name where the plots will be saved')
    parser.add_argument('--hide', action = 'store_true',
    help = 'Do not show plots produced')
    parser.add_argument('--draw-galaxy', action = 'store_true',
    help = 'Show original galaxies')
    parser.add_argument('--partials', action = 'store_true',
    help = 'Show partial derivative images.')
    parser.add_argument('--second-partials', action = 'store_true',
    help = 'Show second partial derivative images.')
    parser.add_argument('--fisher', action = 'store_true',
    help = 'Show fisher matrix images.')
    parser.add_argument('--fisher-chi2', action = 'store_true',
    help = 'Show fisher matrix images with chi2 calculated values on top.')
    parser.add_argument('--covariance', action = 'store_true',
    help = 'Show covariance matrix elements')
    parser.add_argument('--correlation', action = 'store_true',
    help = 'Show correlation matrix elements.')   
    parser.add_argument('--bias-matrix', action = 'store_true',
    help = 'Show bias matrix images.')
    parser.add_argument('--biases', action = 'store_true',
    help = 'Show bias images.')
    parser.add_argument('--values', action = 'store_true', 
    help = 'Show values on top of appropiate images (fisher, biases).')
    parser.add_argument('--sigma-n', metavar = 'SIGMA_N', default = 1.,
    help = 'Value of noise bias (standard deviation) in each pixel.')

    ##maybe change steps from default. 


    names = defaults.names() #initialize names from default that can be useful parameters,etc.

    args = parser.parse_args()
    #some checks.
    if not os.path.isdir(args.wdir):
        print ('Directory does not exists')
        return -1

    filename = os.path.join(args.wdir, args.galaxy_file + '.csv')
    if not os.path.isfile(filename):
        print('Galaxies file does not exist')
        return -1

    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        lst_params = [] #possible here create dictionary with total number of params of both galaxies.
        for row in reader: 
            lst_params.append(row)

    params = lst_params[0] 

    #convert number values in params to floats.
    for key, value in params.iteritems():
        try:
            params[key] = float(value)
        except ValueError:
            pass


    #this is the jiggle in one pixel due to the noise, uniform for all pixels for now too.
    sigma_n = args.sigma_n

    #possible here create all possible galaxies with all params. 
    gal_image = drawGalaxy(params)

    # #initialize steps to be used from defaults.py
    # steps = defaults.steps(params).dct 

    # #initialize name of the parameters that will be plotted and statistically analyzed.
    # param_names = names.galaxy_parameters[params['model']]
    # num_params = len(param_names)

    #create a dictionary with the derivatives of the model with respect to each parameter.
    Ds_gal = derivativesGalaxy(gal_image = gal_image, params = params)
    # Ds_gal = {param_names[i]:partialDifferentiate(func = drawGalaxy, parameter = param_names[i], step = steps[param_names[i]])(params).array for i in range(num_params)}




    ##produce images for analysis.
    FisherM_images = fisherMatrixImages(gal_image = gal_image, params = params)
    # for i in range(num_params): 
    #     for j in range(num_params): 
    #         FisherM_images[param_names[i],param_names[j]] = (Ds_gal[param_names[i]] * Ds_gal[param_names[j]]) /(sigma_n**2)

    FisherM = fisherMatrix(gal_image = gal_image, params = params)
    # for i in range(num_params): 
    #     for j in range(num_params): 
    #         FisherM[param_names[i],param_names[j]] = FisherM_images[param_names[i],param_names[j]].sum() #sum over all pixels. 


    #or with chi2
    #but need to expected value (average over many galaxies with different noises)... although if the galaxy is noiseless, burchat argues it is the same or very close.
    FisherM_chi2 = fisherMatrixChi2(gal_image,params)
    # for i in range(num_params): 
    #     for j in range(num_params): 
    #         FisherM_chi2[param_names[i],param_names[j]] = .5 * secondPartialDifferentiate(chi2, param_names[i], param_names[j], steps[param_names[i]], steps[param_names[j]], sigma_n = sigma_n, gal_image = gal_image)(params)

    #do second derivatives of galaxies.
    # SecondDs_gal = {}
    # for i in range(num_params): 
    #     for j in range(num_params):
    #         SecondDs_gal[param_names[i],param_names[j]] = (secondPartialDifferentiate(drawGalaxy, param_names[i], param_names[j], steps[param_names[i]], steps[param_names[j]])(params).array)

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
    CovM_array = np.linalg.inv(FisherM_array) #need to be an array to inverse.

    for i in range(num_params):
        for j in range(num_params):
            CovM[param_names[i],param_names[j]] = CovM_array[i][j] #dictionary with covarianc matrix elements for completeness

    #now we want bias of each parameter per pixel, so we can see how each parameter contributes. 
    bias_images = {}
    for i in range(num_params):
        sumation = 0 
        for j in range(num_params):
            for k in range(num_params):
                for l in range(num_params):
                    sumation += CovM[param_names[i],param_names[j]]*CovM[param_names[k],param_names[l]]*BiasM_images[param_names[j],param_names[k],param_names[l]]
        bias_images[param_names[i]] = (-.5) * sumation

    biases = {param_names[i]:bias_images[param_names[i]].sum() for i in range(num_params)}


    ##plotting.
    plt.rc('text', usetex=False)
    if args.draw_galaxy:
        figure1, subplt= plt.subplots(1,1)
        figure1.suptitle('Initial Galaxy', fontsize = 20)
        drawPlot(subplt, gal_image.array)
        SaveFigureToPdf(figure1, 'figure1.png', args.wdir, args.plots_dir, hide = args.hide)

    if args.partials:
        figure2 = plt.figure() 
        figure2.suptitle('Derivatives of model with respect to each parameter', fontsize = 20)
        for i in range(num_params):
            ax = figure2.add_subplot(1,num_params,i+1)
            drawPlot(ax, Ds_gal[param_names[i]], title = param_names[i])
        SaveFigureToPdf(figure2, 'figure2.png', args.wdir, args.plots_dir, hide = args.hide)

    if args.fisher and args.values:
        figure4 = plt.figure()
        figure4.suptitle('Fisher matrix elements with values ', fontsize=14, fontweight='bold')
        for i in range(num_params): 
            for j in range(num_params):
                if(i >= j):
                    ax = figure4.add_subplot(num_params,num_params, num_params * i + j + 1)
                    drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(param_names[i])
                    if(i == num_params - 1):
                        ax.set_xlabel(param_names[j])
                    ax.text(-20,0, str(round(FisherM[param_names[i],param_names[j]],5)), fontweight='bold')
        SaveFigureToPdf(figure4, 'figure4.png', args.wdir, args.plots_dir, hide = args.hide)

    elif args.fisher:
        figure3 = plt.figure()
        figure3.suptitle('Fisher matrix elements ', fontsize=14, fontweight='bold')
        for i in range(num_params): 
                for j in range(num_params):
                    if(i >= j):
                        ax = figure3.add_subplot(num_params,num_params, num_params * i + j + 1)
                        drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                        if(j == 0):
                            ax.set_ylabel(param_names[i] )
                        if(i == num_params - 1):
                            ax.set_xlabel(param_names[j])
        SaveFigureToPdf(figure3, 'figure3.png', args.wdir, args.plots_dir, hide = args.hide)



    if args.fisher_chi2:
        figure5 = plt.figure()
        figure5.suptitle('Fisher matrix elements with values of chi2 method ', fontsize=14, fontweight='bold')
        for i in range(num_params): 
            for j in range(num_params):
                if(i >= j):
                    ax = figure5.add_subplot(num_params,num_params, num_params * i + j + 1)
                    drawPlot(ax, FisherM_images[param_names[i],param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(param_names[i] )
                    if(i == num_params - 1):
                        ax.set_xlabel(param_names[j])
                    ax.text(-20,0, str(round(FisherM_chi2[param_names[i],param_names[j]],5)), fontweight='bold')
        SaveFigureToPdf(figure5, 'figure5.png', args.wdir, args.plots_dir, hide = args.hide)

    if args.covariance:
        figure6 = plt.figure()
        figure6.suptitle('Covariance matrix elements', fontsize=14, fontweight='bold')
        for i in range(num_params): 
            for j in range(num_params):
                if(i >= j):
                    ax = figure6.add_subplot(num_params,num_params, num_params * i + j + 1)
                    drawPlot(ax, FisherM_images[param_names[i],param_names[j]] * 0)
                    if(j == 0):
                        ax.set_ylabel(param_names[i])
                    if(i == num_params - 1):
                        ax.set_xlabel(param_names[j])
                    ax.text(-20,0, str(round(CovM[param_names[i],param_names[j]],5)), fontweight='bold')
        SaveFigureToPdf(figure6, 'figure6.png', args.wdir, args.plots_dir, hide = args.hide)


    ########################change this to correlation and do the color thing to show when more correlated like David.
    if args.correlation:
        figure7 = plt.figure()
        figure7.suptitle('Covariance matrix elements', fontsize=14, fontweight='bold')
        for i in range(num_params): 
            for j in range(num_params):
                if(i >= j):
                    ax = figure7.add_subplot(num_params,num_params, num_params * i + j + 1)
                    drawPlot(ax, FisherM_images[param_names[i],param_names[j]] * 0)
                    if(j == 0):
                        ax.set_ylabel(param_names[i])
                    if(i == num_params - 1):
                        ax.set_xlabel(param_names[j])
                    ax.text(-20,0, str(round(CovM[param_names[i],param_names[j]],5)), fontweight='bold') 
        SaveFigureToPdf(figure7, 'figure7.png', args.wdir, args.plots_dir, hide = args.hide)

    if args.second_partials:
        figure8 = plt.figure()
        figure8.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=14, fontweight='bold')
        for i in range(num_params): 
            for j in range(num_params):
                if(i >= j):
                    ax = figure8.add_subplot(num_params,num_params, num_params * i + j + 1)
                    drawPlot(ax, SecondDs_gal[param_names[i],param_names[j]])
                    ax.text(0,20, "std:" + str(SecondDs_gal[param_names[i],param_names[j]].std().round(4)), fontsize = 8, fontweight='bold') 
                    if(j == 0):
                        ax.set_ylabel(param_names[i])
                    if(i == num_params - 1):
                        ax.set_xlabel(param_names[j])
        SaveFigureToPdf(figure8, 'figure8.png', args.wdir, args.plots_dir, hide = args.hide)

    if args.bias_matrix and args.values: 
        figuresOfBiasMatrixNumbers = []    
        for i in range(num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
            for j in range(num_params): 
                for k in range(num_params):
                    if(j >= k):
                        ax = figure.add_subplot(num_params,num_params, num_params * j + k + 1)
                        drawPlot(ax, BiasM_images[param_names[i],param_names[j],param_names[k]]) 
                        if(k == 0):
                            ax.set_ylabel(param_names[j] )
                        if(j == num_params - 1):
                            ax.set_xlabel(param_names[k])
                        ax.text(-20,0, str(round(BiasM[param_names[i],param_names[j],param_names[k]],5)), fontweight='bold')
            figuresOfBiasMatrixNumbers.append(figure)

        for i, figure in enumerate(figuresOfBiasMatrixNumbers): 
            SaveFigureToPdf(figure, 'figure' + str(10) + '_' + str(i) + '.png', args.wdir, args.plots_dir, hide = args.hide) 

    elif args.bias_matrix:
        figuresOfBiasMatrix = [] 
        for i in range(num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
            for j in range(num_params): 
                for k in range(num_params):
                    if(j >= k):
                        ax = figure.add_subplot(num_params,num_params, num_params * j + k + 1)
                        drawPlot(ax, BiasM_images[param_names[i],param_names[j],param_names[k]]) 
                        if(k == 0):
                            ax.set_ylabel(param_names[j])
                        if(j == num_params - 1):
                            ax.set_xlabel(param_names[k])
            figuresOfBiasMatrix.append(figure)
        for i, figure in enumerate(figuresOfBiasMatrix): 
            SaveFigureToPdf(figure, 'figure' + str(9) + '_' + str(i) + '.png', args.wdir, args.plots_dir, hide = args.hide) 

    if args.biases and args.values:
        figure12 = plt.figure() 
        figure12.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 14)
        for i in range(num_params):
            ax = figure12.add_subplot(1,num_params,i+1)
            drawPlot(ax, bias_images[param_names[i]], title = param_names[i])
            ax.text(-20,0, str(round(biases[param_names[i]],5)), fontweight='bold')

        SaveFigureToPdf(figure12, 'figure12.png', args.wdir, args.plots_dir, hide = args.hide)

    elif args.biases:
        figure11 = plt.figure() 
        figure11.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 15)
        for i in range(num_params):
            ax = figure11.add_subplot(1,num_params,i+1)
            drawPlot(ax, bias_images[param_names[i]], title = param_names[i])

        SaveFigureToPdf(figure11, 'figure11.png', args.wdir, args.plots_dir, hide = args.hide)
     




if __name__ == '__main__':
    main()
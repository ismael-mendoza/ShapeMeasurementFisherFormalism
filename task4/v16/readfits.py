#!/usr/bin/env python
"""This gathers the residuals and computes the bias from multiple .csv files scattered in a directory. 
Each file will not necessarily have only one row of results (although it will most of the times"""

import csv

import numpy as np 

import os

import sys

import functions as fns

import triangle

import fisher

import math

import matplotlib.pyplot as plt

import matplotlib.mlab as mlab

import error_ellipse

import defaults

def main(argv):
    wdir, galaxy_file, rltsdir, snr, plots_dir, verbose, info_file, N = argv[1], argv[2], argv[3], float(argv[4]), argv[5], argv[6], argv[7], int(argv[8])

    names =defaults.names()
    params = fns.getParamsCsv(wdir, galaxy_file)
    #print params

    if not os.path.isdir(os.path.join(wdir, rltsdir)):
        os.mkdir(os.path.join(wdir, rltsdir))

    residuals = {}
    path = os.path.join(wdir,rltsdir)

    for filename in os.listdir(path):
        with open(os.path.join(path,filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for param in row.keys():
                    #store every residual in a list and then can calculate mean. 
                    if param not in residuals.keys(): 
                        residuals[param] = [] 
                    #calculate residual for each row. 
                    residuals[param].append(float(params[param]) - float(row[param]))

    biases = {param:np.mean(residuals[param]) for param in residuals.keys()}
    dummy_image, variance_noise = fns.drawGalaxy(params, snr = snr, noise_seed = 0) #just to get the variance_noise
    
    #fisher object contains results of fisher analysis
    fish = fisher.fisher(params = params, gal_image = fns.drawGalaxy(params), sigma_n = math.sqrt(variance_noise))


    info = defaults.info(params, fish, fits_biases = biases, N = N)
    if verbose: 
        for line in info.galaxy + info.fisher:
            print line

    if info_file:
        with open(os.path.join(wdir, names.info + '.txt'), 'w') as txtfile:
            for line in info.galaxy + info.fisher + info.fits:
                txtfile.write(line + '\n') #last character to skip lines.
    
    #produce triangle plots.


    #draw Gaussians, 
            #     import matplotlib.pyplot as plt
            # import numpy as np
            # import matplotlib.mlab as mlab
            # import math

            # mean = 0
            # variance = 1
            # sigma = math.sqrt(variance)
            # x = np.linspace(-3,3,100)
            # plt.plot(x,mlab.normpdf(x,mean,sigma))

            # plt.show()

    #draw error ellipses, 
        #triangle.error_ellipse()

    #first we want triangle plot of fisher analysis. 
    fish_figure = plt.figure() 
    for i in range(fish.num_params): 
        for j in range(fish.num_params):
            #obtain 2d mean 
            mean = fish.biases[fish.param_names[i]], fish.biases[fish.param_names[j]]
            #obtain 2x2 cov matrix of the corresponding i,j elements.  
            cov = np.array([[fish.covariance_matrix[fish.param_names[a],fish.param_names[b]] for a in [i,j]] for b in [i,j]])
            #print cov
            if(i == j):
                sigma = math.sqrt(cov[0][0])
                #draw a gaussian for the diagonals
                    #in theory the width should some multiple of the variance of a given parameter. maybe 3, just a guess, 
                    #the number of subdivisions does not matter as long as its large enough.

                x = np.linspace(mean[0] - sigma*3, mean[0] + sigma*3, 1000)
                ax = fish_figure.add_subplot(fish.num_params,fish.num_params, fish.num_params * i + j + 1)
                ax.plot(x,mlab.normpdf(x,mean[0],sigma))
            elif(i > j):
                #draw an error_ellipse in off-diagonals.
                ellip = error_ellipse.plot_cov_ellipse(pos = mean, cov = cov, nstd = 1)
                ax = fish_figure.add_subplot(fish.num_params,fish.num_params, fish.num_params * i + j + 1)
                ax.add_artist(ellip)
                ax.set_xlim(mean[0] - 2, mean[0] + 2)
                ax.set_ylim(mean[1] - 2, mean[1] + 2)

            if(j == 0):
                ax.set_ylabel(fish.param_names[i])
            if(i == fish.num_params - 1):
                ax.set_xlabel(fish.param_names[j])

    #now we want to to do a triangle plot from triangle.py and overlap this figure. 

    #have to transform residual.values to a convenient form for the triangle plot. 
    # figure = triangle.corner(np.array(residuals.values()).transpose(), labels=fish.param_names,
    #                      truths= fish.biases.values(),
    #                      show_titles=True, title_args={"fontsize": 12})
    # figure.gca().annotate("Triangle plot", xy=(0.5, 1.0), xycoords="figure fraction",
    #                       xytext=(0, -5), textcoords="offset points",
    #                       ha="center", va="top")
    fish_figure.savefig(os.path.join(wdir,plots_dir,"demo.png"))



if __name__ == "__main__":
    main(sys.argv)
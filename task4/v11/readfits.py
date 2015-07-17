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

def main(argv):
    wdir, galaxy_file, rltsdir, snr = argv[1], argv[2], argv[3], argv[4]

    params = fns.getParamsCsv(wdir, galaxy_file)

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
                    if not residuals[param] == []: 
                        residuals[param] = [] 
                    #calculate residual for each row. 
                    residuals[param].append(float(params[param]) - float(row[param]))


    biases = {param:np.mean(residuals[param]) for param in residuals.keys()}
    fisher_analysis = fisher.fisher(params, fns.drawGalaxy(params), snr)
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
    figure = plt.figure() 
    for i in range(fisher.num_params): 
        for j in range(fisher.num_params):
            #obtain 2d mean 
            mean = fisher.biases[fisher.param_names[i]], fisher.biases[fisher.param_names[j]]
            #obtain 2x2 cov matrix of the corresponding i,j elements.  
            cov = np.array([[fisher.covariance_matrix[fisher.param_names[a],fisher.param_names[b]] for a in [i,j]] for b in [i,j]])
            if(i == j):
                sigma = math.sqrt(cov[i][i])
                #draw a gaussian for the diagonals
                    #in theory the width should some multiple of the variance of a given parameter. maybe 3, just a guess, 
                    #the number of subdivisions does not matter 

                x = np.linspace(mean[0] - sigma*3, mean[0] + sigma*3, 1000)
                ax = figure.add_subplot(fisher.num_params,fisher.num_params, fisher.num_params * i + j + 1)
                ax.plot(x,mlab.normpdf(x,mean[0],sigma))
                # ax.text(0,20, "std:" + str(fisher.second_derivatives_galaxy_images[self.param_names[i],self.param_names[j]].std().round(4)), fontsize = 8, fontweight='bold') 
            elif(i > j):
                #draw an error_ellipse in off-diagonals.

                ellip = error_ellipse.plot_cov_ellipse(mean, cov, nstd = 1)
                ax = figure.add_subplot(fisher.num_params,fisher.num_params, fisher.num_params * i + j + 1)
                ax.add_artist(ellip)

            if(j == 0):
                ax.set_ylabel(fisher.param_names[i])
            if(i == fisher.num_params - 1):
                ax.set_xlabel(fisher.param_names[j])

    #now we want to to do a triangle plot from triangle.py and overlap this figure. 
    figure = triangle.corner(residuals.values(), labels=fisher.param_names,
                         truths= biases.values(),
                         show_titles=True, title_args={"fontsize": 12})
    figure.gca().annotate("A Title", xy=(0.5, 1.0), xycoords="figure fraction",
                          xytext=(0, -5), textcoords="offset points",
                          ha="center", va="top")
    figure.savefig("demo.png")



if __name__ == "__main__":
    main(sys.argv)
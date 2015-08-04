#!/usr/bin/env python
"""This gathers the residuals and computes the bias from multiple .csv files scattered in a directory. 
Each file will not necessarily have only one row of results (although it will most of the times"""

import csv

import numpy as np 

import os

import sys

import triangle

import fisher

import math

import matplotlib.pyplot as plt

import matplotlib.mlab as mlab

import matplotlib.patches as mpatch

import defaults

import galfun

import info

def errorEllipseCor(centroid,cor_xy,alpha = 1.52):
    """Calculation for ellipse according to Dan Coe paper arxiv.0906.4123.
    One change we did was assuming that the units of x,y are not necessarily
    the same so that we instead plot the dimensionless correlation_matrix
    instead of the covariance_matrix.
    Confidence level by default is 68.3% -> alpha = 1.52 according to the 
    same paper. 
    """
    a = 1 + cor_xy
    b = 1 - cor_xy
    theta = 45
    width = 2 * alpha * a
    height = 2 * alpha * b

    ellipse = mpatch.Ellipse(xy=centroid, width=width, height=height, 
                             angle = theta, edgecolor='r', fc='None', lw=2)

    return ellipse

def main(argv):
    (project, snr, verbose, info_file) = (
     argv[1], float(argv[2]), argv[3], argv[4])

    g_parameters = galfun.GParameters(project)
    image = galfun.drawGalaxies(g_parameters=g_parameters, image=True)
    _, variance_noise = galfun.addNoise(image, snr, 0) 
    
    #fisher object contains results of fisher analysis
    fish = fisher.Fisher(g_parameters=g_parameters, snr=snr)

    if not os.path.isdir(os.path.join(project, defaults.RESULTS_DIR)):
        os.mkdir(os.path.join(project, defaults.RESULTS_DIR))

    residuals = {}
    pulls = {}
    rltsdir = os.path.join(project,defaults.RESULTS_DIR)
    number_fits = 0

    #generate residuals by reading every file in rltsdir
    for filename in os.listdir(rltsdir):
        number_fits += 1
        with open(os.path.join(rltsdir,filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for param in row.keys():
                    if param not in residuals.keys(): 
                        residuals[param] = []
                    if param not in pulls.keys():
                        pulls[param] = [] 

                    residuals[param].append(float(g_parameters.params[param])
                                            -float(row[param]))
                    #pull
                    pulls[param].append((float(g_parameters.params[param])
                                            -float(row[param]))/(math.sqrt(fish.covariance_matrix[param,param])))

    biases = {param:np.mean(residuals[param]) for param in residuals.keys()}
    stds = {param:np.std(residuals[param]) for param in residuals.keys()}

    #produce pull plot. 
    norm_area = (defaults.EXTENT_PULL[1]-defaults.EXTENT_PULL[0]) * (number_fits/defaults.BINS_PULL)

    fish_figure = plt.figure() 
    for i in range(fish.num_params): 
        for j in range(fish.num_params):
            sigma_i = math.sqrt(fish.covariance_matrix[fish.param_names[i],
                                                       fish.param_names[i]])
            sigma_j = math.sqrt(fish.covariance_matrix[fish.param_names[j],
                                                       fish.param_names[j]])
            mean = (fish.biases[fish.param_names[i]]/sigma_i,
                    fish.biases[fish.param_names[j]]/sigma_j)
            cor_xy = fish.correlation_matrix[fish.param_names[i],
                                             fish.param_names[j]
                                            ]
            sigma_fits = stds[fish.param_names[i]]
            if(i == j):

                sigma_gauss = 1 
                x = np.linspace(defaults.EXTENT_PULL[0], defaults.EXTENT_PULL[1], 1000)
                ax = fish_figure.add_subplot(fish.num_params,fish.num_params, fish.num_params * i + j + 1)
                ax.plot(x, norm_area * mlab.normpdf(x,mean[0],sigma_gauss))

                ax.text(.5,.4,'sfh:' + str(round(sigma_i,defaults.SIG_DIGITS)),
                        transform=ax.transAxes, ha='center',
                        fontsize=defaults.FONTSIZE_VALUE,fontweight='bold')
                ax.text(.5,.3,'sfs:' + str(round(sigma_fits,defaults.SIG_DIGITS)) + '+-' + str(round((sigma_fits)/(math.sqrt(2*number_fits)),defaults.SIG_DIGITS)),
                    fontsize=defaults.FONTSIZE_VALUE, ha='center',
                        transform=ax.transAxes,fontweight='bold')
                ax.text(.5,.2,'bfh:' + str(round(fish.biases[fish.param_names[i]],defaults.SIG_DIGITS)), fontsize=defaults.FONTSIZE_VALUE, transform=ax.transAxes,ha='center',fontweight='bold')
                ax.text(.5,.1,'bfs:' + str(round(biases[fish.param_names[i]],defaults.SIG_DIGITS)) + '+-' + str(round(sigma_fits/math.sqrt(number_fits),defaults.SIG_DIGITS)), 
                 fontsize = defaults.FONTSIZE_VALUE,ha='center',
                 transform=ax.transAxes,fontweight='bold')
            elif(i > j):
                #draw an error_ellipse in off-diagonals.
                ellip = errorEllipseCor(mean,cor_xy)
                ax = fish_figure.add_subplot(fish.num_params,fish.num_params, 
                     fish.num_params * i + j + 1)
                ax.add_patch(ellip)
            else:
                ax = fish_figure.add_subplot(fish.num_params,fish.num_params, fish.num_params * i + j + 1)

    #now we want to to do a triangle plot from triangle.py and overlap this figure. 

    #have to transform point.values to a convenient form for the triangle plot. 
    extents = [defaults.EXTENT_PULL]*6
    points_plot = []
    truths = []
    #just so the ordering is correct.
    for param in fish.param_names:
        points_plot.append(pulls[param])
        truths.append(biases[param])
    
    figure = triangle.corner(np.array(points_plot).transpose(), 
                             bins = defaults.BINS_PULL, 
                             labels=fish.param_names,
                             extents = extents,
                             truths= truths, plot_contours=True,
                             show_titles=True,
                             title_args={"fontsize": 12}, fig = fish_figure)

    plt.show()
    #strange error.
    #figure.savefig(os.path.join(project,defaults.TRIANGLE_NAME))

    information = info.Info(g_parameters, fish, fits_biases = biases, 
                         number_fits = number_fits)

    if verbose == 'True':
        information.printInfo()

    if info_file == 'True':
        information.writeInfo()

if __name__ == "__main__":
    main(sys.argv)
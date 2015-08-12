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
    project, snr, number_fits = argv[1], float(argv[2]), float(argv[3])

    g_parameters = galfun.GParameters(project)
    orig_image = galfun.drawGalaxies(g_parameters=g_parameters, image=True)
    mins = defaults.getMinimums(g_parameters,orig_image)
    maxs = defaults.getMaximums(g_parameters,orig_image)
    fish = fisher.Fisher(g_parameters=g_parameters, snr=snr)

    if not os.path.isdir(os.path.join(project, defaults.RESULTS_DIR)):
        os.mkdir(os.path.join(project, defaults.RESULTS_DIR))

    residuals = {}
    pulls = {}
    rltsdir = os.path.join(project,defaults.RESULTS_DIR)

    #produce results from rltsdir's files.
    for filename in os.listdir(rltsdir):
        with open(os.path.join(rltsdir,filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for param in row: 
                    if param not in residuals:
                        residuals[param] = []
                    if param not in pulls:
                        pulls[param] = []
                    residual = (float(row[param])-
                               float(g_parameters.params[param]))
                    pull = (residual/ 
                            math.sqrt(fish.covariance_matrix[param,param]))

                    residuals[param].append(residual)
                    pulls[param].append(pull)

    biases = {param:np.mean(residuals[param]) for param in residuals}
    pull_means = {param:np.mean(pulls[param]) for param in residuals}
    res_stds = {param:np.std(residuals[param]) for param in residuals}
    pull_mins = {param:((mins[param]-float(g_parameters.params[param]))/ 
                 math.sqrt(fish.covariance_matrix[param,param])) for
                 param in residuals}
    pull_maxs = {param:((maxs[param]-float(g_parameters.params[param]))/ 
                 math.sqrt(fish.covariance_matrix[param,param])) for
                 param in residuals}
    #produce pull plot. 
    norm_area = (defaults.EXTENT_PULL[1]-defaults.EXTENT_PULL[0]) * (
                 number_fits/defaults.BINS_PULL)

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
            sigma_fits = res_stds[fish.param_names[i]]
            if(i == j):

                sigma_gauss = 1 
                x = np.linspace(defaults.EXTENT_PULL[0], defaults.EXTENT_PULL[1], 1000)
                ax = fish_figure.add_subplot(fish.num_params,fish.num_params, fish.num_params * i + j + 1)
                ax.plot(x, norm_area * mlab.normpdf(x,mean[0],sigma_gauss))

                ax.text(.5,.4,'sfh:' + str(round(sigma_i,defaults.SIG_DIGITS)),
                        transform=ax.transAxes, ha='center',
                        fontsize=defaults.FONTSIZE_VALUE,fontweight='bold')
                ax.text(.5,.3,'sft:' + str(round(sigma_fits,defaults.SIG_DIGITS)) + '+-' + str(round((sigma_fits)/(math.sqrt(2*number_fits)),defaults.SIG_DIGITS)),
                    fontsize=defaults.FONTSIZE_VALUE, ha='center',
                        transform=ax.transAxes,fontweight='bold')
                ax.text(.5,.2,'bfh:' + str(round(fish.biases[fish.param_names[i]],defaults.SIG_DIGITS)), fontsize=defaults.FONTSIZE_VALUE, transform=ax.transAxes,ha='center',fontweight='bold')
                ax.text(.5,.1,'bft:' + str(round(biases[fish.param_names[i]],defaults.SIG_DIGITS)) + '+-' + str(round(sigma_fits/math.sqrt(number_fits),defaults.SIG_DIGITS)), 
                 fontsize = defaults.FONTSIZE_VALUE,ha='center',
                 transform=ax.transAxes,fontweight='bold')
                ax.axvline(x=pull_mins[fish.param_names[j]],color='g')
                ax.axvline(x=pull_maxs[fish.param_names[j]],color='g')
            elif(i > j):
                #draw an error_ellipse in off-diagonals.
                ellip = errorEllipseCor(mean,cor_xy)
                ax = fish_figure.add_subplot(fish.num_params,fish.num_params, 
                     fish.num_params * i + j + 1)
                ax.add_patch(ellip)
                #add bound lines.
                ax.axhline(y=pull_mins[fish.param_names[i]],color='g')
                ax.axhline(y=pull_maxs[fish.param_names[i]],color='g')
                ax.axvline(x=pull_mins[fish.param_names[j]],color='g')
                ax.axvline(x=pull_maxs[fish.param_names[j]],color='g')
            else:
                ax = fish_figure.add_subplot(fish.num_params,
                     fish.num_params, fish.num_params * i + j + 1)



    #have to transform point.values to a convenient form for the triangle plot. 
    extents = [defaults.EXTENT_PULL]*6
    points_plot = []
    truths = []

    #just so the ordering is consistent.
    for param in fish.param_names:
        points_plot.append(pulls[param])
        truths.append(pull_means[param])
    
    figure = triangle.corner(np.array(points_plot).transpose(), 
                             bins = defaults.BINS_PULL, 
                             labels=fish.param_names,
                             extents = extents,
                             truths= truths, plot_contours=True,
                             show_titles=True,
                             title_args={"fontsize": defaults.FONTSIZE_LABEL},
                             fig = fish_figure)


    figure.savefig(os.path.join(project,defaults.TRIANGLE_NAME))

if __name__ == "__main__":
    main(sys.argv)
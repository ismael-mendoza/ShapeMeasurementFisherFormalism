#!/usr/bin/env python
"""This gathers the residuals and computes the bias from multiple .csv files
scattered in a directory.
Each file will not necessarily have only one row of results (although it will
most of the times.
"""

#todo
#looks very messy, maybe pass this to an ipython notebook? maybe only the fancy part... 
#not super urgent though. 
# 

import csv

import numpy as np

import os

import sys

import triangle

import math

import matplotlib.pyplot as plt

import matplotlib.mlab as mlab

import matplotlib.patches as mpatch

import defaults

import galfun

import fisher

import copy

def errorEllipseCor(centroid, cor_xy, alpha=1.52):
    """Calculation for cor. ellipse according to
    One change we did was assuming that the units of x,y are not necessarily
    the same so that we instead plot the dimensionless correlation_matrix
    instead of the covariance_matrix.
    Confidence level by default is 68.3% -> alpha = 1.52 according to the
    same paper.
    """

    """Return a correlation ellipse according to Dan Coe paper arxiv.0906.4123.

    We use correlations instead of covariances for the ellipses because most of
    the parameters in our models have dimensions.
    Confidence level by default is 68.3% -> alpha = 1.52 according to the same
    paper.

    Args:
        centroid(:py:tuple): Position of the center of the ellipse in
                             (x,y) form.

    Returns:
        A :py:dict.
    """
    a = 1 + cor_xy
    b = 1 - cor_xy
    theta = 45
    width = 2 * alpha * a
    height = 2 * alpha * b

    ellipse = mpatch.Ellipse(xy=centroid, width=width, height=height,
                             angle=theta, edgecolor='r', fc='None', lw=2)

    return ellipse


def main(argv):
    project, snr, number_fits = argv[1], float(argv[2]), float(argv[3])

    g_parameters = galfun.GParameters(project)
    orig_image = galfun.drawGalaxies(g_parameters=g_parameters, image=True)
    mins = defaults.getMinimums(g_parameters, orig_image)
    maxs = defaults.getMaximums(g_parameters, orig_image)
    fish = fisher.Fisher(g_parameters=g_parameters, snr=snr)
    param_names = fish.param_names
    num_params = fish.num_params

    if not os.path.isdir(os.path.join(project, defaults.RESULTS_DIR)):
        os.mkdir(os.path.join(project, defaults.RESULTS_DIR))

    residuals = {}
    pulls = {}
    redchis = [] #list containing values of reduced chi2 for each fit.
    rltsdir = os.path.join(project, defaults.RESULTS_DIR)

    # produce results from rltsdir's files.
    for filename in os.listdir(rltsdir):
        with open(os.path.join(rltsdir, filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                redchis.append(float(row['redchi']))
                for param in g_parameters.fit_params:
                    if param not in residuals:
                        residuals[param] = []
                    if param not in pulls:
                        pulls[param] = []
                    residual = (float(row[param]) -
                                float(g_parameters.params[param]))
                    pull = (residual /
                            math.sqrt(fish.covariance_matrix[param, param]))

                    residuals[param].append(residual)
                    pulls[param].append(pull)

    biases = {param: np.mean(residuals[param]) for param in residuals}
    pull_means = {param: np.mean(pulls[param]) for param in residuals}
    res_stds = {param: np.std(residuals[param]) for param in residuals}
    pull_mins = {param: ((mins[param] - float(g_parameters.params[param])) /
                         math.sqrt(fish.covariance_matrix[param, param])) for
                 param in residuals}
    pull_maxs = {param: ((maxs[param] - float(g_parameters.params[param])) /
                         math.sqrt(fish.covariance_matrix[param, param])) for
                 param in residuals}
    for param in pulls:


    # produce pull plot.
    normalized_area = (defaults.EXTENT_PULL[1] - defaults.EXTENT_PULL[0]) * (
        number_fits / defaults.BINS_PULL)

    if defaults.FANCY:
        fish_figure = plt.figure(figsize=(11, 11))
    else:
        fish_figure = plt.figure()

    dim = num_params + 1  # for redchi add extra dimension.
    for i in range(dim):
        for j in range(dim):
            if i < num_params and j < num_params and i >= j:
                param_i = param_names[i]
                param_j = param_names[j]
                sigma_i = math.sqrt(fish.covariance_matrix[param_i,
                                                           param_i])
                sigma_j = math.sqrt(fish.covariance_matrix[param_j,
                                                           param_j])
                mean = (fish.biases[param_i] / sigma_i,
                        fish.biases[param_j] / sigma_j)
                cor_xy = fish.correlation_matrix[param_i, param_j]
                if i == j:
                    sigma_gauss = 1 #normalized normal pdf.
                    x = np.linspace(defaults.EXTENT_PULL[0],
                                    defaults.EXTENT_PULL[1], 1000)
                    ax = fish_figure.add_subplot(dim, dim,
                                                 dim * i + j + 1)
                    ax.plot(x, normalized_area * mlab.normpdf(x, mean[0],
                                                        sigma_gauss))
                    sigma_fisher = round(sigma_i, defaults.SIG_DIGITS)
                    sigma_fits = round(res_stds[param_i], defaults.SIG_DIGITS)

                    #error of bias_fits
                    sigma_1 = round(sigma_fits / math.sqrt(number_fits),
                                    defaults.SIG_DIGITS)
                    #error of sigma_fits
                    sigma_2 = round(sigma_fits / math.sqrt(2*number_fits),
                                    defaults.SIG_DIGITS)
                    bias_fisher = round(fish.biases[param_i],
                                        defaults.SIG_DIGITS)
                    bias_fits = round(biases[param_i],
                                      defaults.SIG_DIGITS)

                    data = (r'{\setlength\arraycolsep{0.1em}' +
                          r'\begin{eqnarray*}' +
                          r'\sigma_{F}' + '& = &' +
                          str(sigma_fisher) +
                          r'\\' +
                          r'\sigma_{fits}' + '& = &' +
                          str(sigma_fits) +
                          r'\pm' +
                          str(sigma_2) +
                          r'\\' +
                          r'b_{F}' + '& = &' +
                          str(bias_fisher) +
                          r'\\' +
                          r'b_{fits}' + '& = &' +
                          str(bias_fits) +
                          r'\pm' +
                          str(sigma_1) +
                          r'\end{eqnarray*}'+
                          r'}')


                    ax.text(1.65, .35, data,
                            transform=ax.transAxes, ha='center',
                            fontsize=10,
                            fontweight='bold')

                    # add green lines indicating bounds.
                    ax.axvline(x=pull_mins[param_j], color='g')
                    ax.axvline(x=pull_maxs[param_j], color='g')

                elif i > j:
                    # draw an error_ellipse in off-diagonals axises.
                    ellip = errorEllipseCor(mean, cor_xy)
                    ax = fish_figure.add_subplot(dim, dim,
                                                 dim * i + j + 1)
                    ax.add_patch(ellip)

                    ax.axhline(y=pull_mins[param_i], color='g')
                    ax.axhline(y=pull_maxs[param_i], color='g')
                    ax.axvline(x=pull_mins[param_j], color='g')
                    ax.axvline(x=pull_maxs[param_j], color='g')

            else:
                ax = fish_figure.add_subplot(dim,
                                             dim,
                                             dim * i + j + 1)

    # have to transform point.values to a convenient form for the plot.
    extents = [defaults.EXTENT_PULL] * num_params + [(min(redchis),
                                                      max(redchis))]
    points_plot = []
    truths = []
    plot_names = []

    #do not include subscript in plots if there is only one galaxy
    if defaults.FANCY:
        for param in param_names:
            if 'e1' in param:
                tex_param = r'$e_{1}$'
            elif 'e2' in param:
                tex_param = r'$e_{2}$'
            elif 'x0' in param:
                tex_param = r'$x_{0}$'
            elif 'y0' in param:
                tex_param = r'$y_{0}$'
            else:
                new_param = param.replace('_1', '')
                tex_param = r'$' + new_param + r'$'
            plot_names.append(tex_param)

        plot_names.append(r'$\chi^{2}/dof$')

    else:
        plot_names = copy.deepcopy(param_names)
        plot_names.append('redchi')

    # just so the ordering is consistent.
    for param in param_names:
        points_plot.append(pulls[param])
        truths.append(pull_means[param])

    #chi2, redchi should be centered around 1
    truths.append(1)
    points_plot.append(redchis)

    figure1 = triangle.corner(np.array(points_plot).transpose(),
                              bins=defaults.BINS_PULL,
                              labels=plot_names,
                              extents=extents,
                              truths=truths, plot_contours=False,
                              show_titles=True,
                              fig=fish_figure)

    if defaults.FANCY:
        figure1.subplots_adjust(hspace=.2, wspace=.2) #adjust to avoid overlap.

    # add redchi individual histogram
    figure2 = plt.figure()
    ax = figure2.add_subplot(111)
    ax.hist(redchis)
    ax.set_title('Reduced chi2 histogram from data.',
                 fontsize=defaults.FONTSIZE_TITLE)

    plt.rc('text', usetex=True)

    figure1.savefig(os.path.join(project, defaults.TRIANGLE_NAME))
    figure2.savefig(os.path.join(project, defaults.REDCHI_HIST_NAME))

if __name__ == "__main__":
    main(sys.argv)

#!/usr/bin/env python
import matplotlib.pyplot as plt

import matplotlib.cm as cm

import defaults

import os

import numpy as np

import copy

import galfun

import fisher

import math

import galsim

def annotateAxisCenter(ax, text):
    ax.text(.5, .5, text,
            ha='center', va='center', transform=ax.transAxes,
            fontweight='bold', fontsize=defaults.FONTSIZE_VALUE)


def saveFigureToPdf(figure, file_name, project, plotsdir='', hide=False):
    # save and preview pdf.
    if not os.path.isdir(project):
        os.mkdir(project)

    if not os.path.isdir(os.path.join(project, plotsdir)):
        os.mkdir(os.path.join(project, plotsdir))

    file_name = os.path.join(project, plotsdir, file_name)
    figure.savefig(file_name, bbox_inches='tight', dpi=defaults.DPI)

    if not hide:
        os.system("open " + file_name)

    # close it so does not consume memory.
    plt.close(figure)


def drawImage(ax, plot, title="", xlabel="", ylabel=""):
    """draws a given plot with default values using imshow in the given
    subplot if they are not subplots should just pass plt. It also adds a
    title with name. -- OWN"""

    # remove numbers and ticks but can use label.
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])

    # labels figure.
    ax.set_title(title, fontsize=defaults.FONTSIZE_LABEL)
    ax.set_xlabel(xlabel, fontsize=defaults.FONTSIZE_LABEL)
    ax.set_ylabel(ylabel, fontsize=defaults.FONTSIZE_LABEL)

    ax.imshow(plot, interpolation='nearest', rasterized=False,
              cmap=cm.RdYlGn, origin='lower',
              extent=[-defaults.NX, defaults.NX, -defaults.NY, defaults.NY],
              vmax=abs(plot).max(), vmin=-abs(plot).max())


class Plots(object):

    """Produce plots for a fish object and displays them, in a given
    plots_dir that is in a given project. Hide if just save images but not
    display them.
    """

    def __init__(self, fish, project, plots_dir, hide, error_bars,
                 bias_sigma):
        self.fish = fish
        self.project = project
        self.plots_dir = plots_dir
        self.num_params = fish.num_params
        self.param_names = fish.param_names
        self.gal_image = fish.image
        self.hide = hide
        self.error_bars = error_bars
        self.bias_sigma = bias_sigma

        # delete plots_dir directory for nameImages to work if already there.
        if os.path.isdir(os.path.join(self.project, self.plots_dir)):
            os.system('rm -r ' + os.path.join(self.project, self.plots_dir))
        else:
            os.mkdir(os.path.join(self.project, self.plots_dir))

    def nameImages(self):
        """Return first base#n.extension string where n is an integer
        that is not a file in project/pltsdir.
        """
        n = 0
        while True:
            n += 1
            filename = ''.join([defaults.FIGURE_BASENAME, str(n),
                                defaults.FIGURE_EXTENSION])
            path = os.path.join(self.project, self.plots_dir, filename)
            figure_name_exists = os.path.isfile(path)
            if not figure_name_exists:
                break

        return filename

    def galaxy(self):
        figure, subplt = plt.subplots(1, 1)
        figure.suptitle('Initial Galaxy', fontsize=defaults.FONTSIZE_TITLE)
        drawImage(subplt, self.gal_image.array)

        saveFigureToPdf(figure, self.nameImages(),
                        self.project, self.plots_dir, hide=self.hide)

    def derivatives(self):
        figure = plt.figure()
        figure.suptitle(
            'Derivatives of model with respect to each parameter',
            fontsize=defaults.FONTSIZE_TITLE)

        for i in range(self.num_params):
            ax = figure.add_subplot(1, self.num_params, i + 1)
            drawImage(ax, self.fish.derivatives_images[self.param_names[i]],
                      title=self.param_names[i])

        saveFigureToPdf(figure, self.nameImages(), self.project,
                        self.plots_dir,
                        hide=self.hide)

    def fisherMatrix(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements',
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params):
            for j in range(self.num_params):
                if i >= j:
                    ax = figure.add_subplot(self.num_params,
                                            self.num_params,
                                            self.num_params * i + j + 1)
                    drawImage(ax,
                              self.fish.fisher_matrix_images[
                                  self.param_names[i],
                                  self.param_names[j]
                              ])

                    if j == 0:
                        ax.set_ylabel(self.param_names[i])
                    if i == self.num_params - 1:
                        ax.set_xlabel(self.param_names[j])

        saveFigureToPdf(figure, self.nameImages(), self.project,
                        self.plots_dir,
                        hide=self.hide)

    def fisherMatrixValues(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values ',
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params):
            for j in range(self.num_params):
                if i >= j:
                    ax = figure.add_subplot(self.num_params,
                                            self.num_params,
                                            self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[
                              self.param_names[i], self.param_names[j]])

                    if j == 0:
                        ax.set_ylabel(self.param_names[i])

                    if i == self.num_params - 1:
                        ax.set_xlabel(self.param_names[j])

                    annotateAxisCenter(ax, str(round(self.fish.fisher_matrix[
                        self.param_names[i], self.param_names[j]], defaults.SIG_DIGITS)))

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir,
                        hide=self.hide)

    def fisherMatrixChi2(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values of chi2 method ',
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params):
            for j in range(self.num_params):
                if i >= j:
                    ax = figure.add_subplot(
                        self.num_params, self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[
                              self.param_names[i], self.param_names[j]])
                    if j == 0:
                        ax.set_ylabel(self.param_names[i])
                    if i == self.num_params - 1:
                        ax.set_xlabel(self.param_names[j])
                    annotateAxisCenter(ax,
                                       str(round(self.fish.fisher_matrix_chi2[
                                           self.param_names[i],
                                           self.param_names[j]],
                                           defaults.SIG_DIGITS)))
        saveFigureToPdf(figure, self.nameImages(), self.project,
                        self.plots_dir, hide=self.hide)

    def secondDerivatives(self):
        figure = plt.figure()
        figure.suptitle('Second derivatives of the galaxies with respect to'
                        'parameters',
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params):
            for j in range(self.num_params):
                if i >= j:
                    ax = figure.add_subplot(
                        self.num_params, self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.second_derivatives_images[
                              self.param_names[i], self.param_names[j]])
                    if j == 0:
                        ax.set_ylabel(self.param_names[i])
                    if i == self.num_params - 1:
                        ax.set_xlabel(self.param_names[j])
        saveFigureToPdf(figure, self.nameImages(), self.project,
                        self.plots_dir, hide=self.hide)

    def biasMatrixValues(self):
        figuresOfBiasMatrixNumbers = []
        for i in range(self.num_params):
            figure = plt.figure()
            param_i = self.param_names[i]
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) +
                            'th partial (D' +
                            param_i + ')',
                            fontsize=defaults.FONTSIZE_TITLE,
                            fontweight='bold')
            for j in range(self.num_params):
                for k in range(self.num_params):
                    param_j = self.param_names[i]
                    param_k = self.param_names[k]
                    if j >= k:
                        ax = figure.add_subplot(
                            self.num_params, self.num_params, self.num_params *
                            j + k + 1)
                        drawImage(ax,
                                  self.fish.bias_matrix_images[
                                    param_i, param_j, param_k])
                        if k == 0:
                            ax.set_ylabel(param_j)
                        if j == self.num_params - 1:
                            ax.set_xlabel(param_k)



                        annotateAxisCenter(ax, str(round(self.fish.bias_matrix[
                                                            param_i, param_j,
                                                            param_k], defaults.SIG_DIGITS)))
            figuresOfBiasMatrixNumbers.append(figure)

        for i, figure in enumerate(figuresOfBiasMatrixNumbers):
            basename = self.nameImages().strip(defaults.FIGURE_EXTENSION)
            filename = ''.join([basename, '_',
                                str(i), defaults.FIGURE_EXTENSION])
            saveFigureToPdf(figure, filename, self.project, self.plots_dir,
                            hide=self.hide)

    def biasMatrix(self):
        figuresOfBiasMatrix = []
        for i in range(self.num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) +
                            'th partial (D' +
                            self.param_names[i] + ')',
                            fontsize=defaults.FONTSIZE_TITLE,
                            fontweight='bold')
            for j in range(self.num_params):
                for k in range(self.num_params):
                    if j >= k:
                        ax = figure.add_subplot(
                            self.num_params, self.num_params, self.num_params * j + k + 1)
                        drawImage(ax, self.fish.bias_matrix_images[self.param_names[
                                  i], self.param_names[j], self.param_names[k]])
                        if k == 0:
                            ax.set_ylabel(self.param_names[j])
                        if j == self.num_params - 1:
                            ax.set_xlabel(self.param_names[k])
            figuresOfBiasMatrix.append(figure)

        for i, figure in enumerate(figuresOfBiasMatrix):
            basename = self.nameImages().replace(defaults.FIGURE_EXTENSION, '')
            filename = ''.join(
                [basename, '_', str(i), defaults.FIGURE_EXTENSION])
            saveFigureToPdf(figure, filename, self.project, self.plots_dir,
                            hide=self.hide)

    def covarianceMatrix(self):
        figure = plt.figure()
        figure.suptitle('Covariance matrix elements',
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]
                if i >= j:
                    ax = figure.add_subplot(
                        self.num_params, self.num_params, self.num_params * i + j + 1)
                    # figure out better way.
                    drawImage(ax, self.fish.fisher_matrix_images[
                              self.param_names[i], self.param_names[j]] * 0)
                    if j == 0:
                        ax.set_ylabel(self.param_names[i])
                    if i == self.num_params - 1:
                        ax.set_xlabel(self.param_names[j])
                    annotateAxisCenter(ax,
                                       str(round(self.fish.covariance_matrix[
                                           param_i, param_j],
                                           defaults.SIG_DIGITS)))

        saveFigureToPdf(figure, self.nameImages(), self.project,
                        self.plots_dir, hide=self.hide)

    def correlationMatrix(self):
        figure = plt.figure()
        figure.suptitle('Correlation matrix elements',
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params):
            for j in range(self.num_params):
                if i >= j:
                    ax = figure.add_subplot(
                        self.num_params, self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[
                              self.param_names[i], self.param_names[j]] * 0)

                    if j == 0:
                        ax.set_ylabel(self.param_names[i])
                    if i == self.num_params - 1:
                        ax.set_xlabel(self.param_names[j])

                    annotateAxisCenter(ax, str(round(self.fish.correlation_matrix[
                                       self.param_names[i], self.param_names[j]], defaults.SIG_DIGITS)))

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir,
                        hide=self.hide)

    def biasValues(self):
        figure = plt.figure()
        figure.suptitle('Contribution to Bias from each pixel by parameter',
                        fontsize=defaults.FONTSIZE_TITLE)
        for i in range(self.num_params):
            ax = figure.add_subplot(1, self.num_params, i + 1)
            param = self.param_names[i]
            drawImage(ax, self.fish.bias_images[param], title=param)
            value = round(self.fish.biases[param], defaults.SIG_DIGITS)

            ax.text(.5, .05, str(value),
                    transform=ax.transAxes, fontsize=defaults.FONTSIZE_VALUE,
                    fontweight='bold')

        saveFigureToPdf(figure, self.nameImages(), self.project,
                        self.plots_dir, hide=self.hide)

    def biases(self):
        figure = plt.figure()
        figure.suptitle('Contribution to Bias from each pixel by parameter',
                        fontsize=defaults.FONTSIZE_TITLE)
        for i in range(self.num_params):
            ax = figure.add_subplot(1, self.num_params, i + 1)
            drawImage(ax, self.fish.bias_images[self.param_names[i]],
                      title=self.param_names[i])

        saveFigureToPdf(figure, self.nameImages(), self.project,
                        self.plots_dir, hide=self.hide)

    def biasPlot1(self):
        """Plot of bias as a function of hlr with fixed snr."""
        steps = 10
        hrange = (.15, .5)
        hlrs = np.linspace(hrange[0], hrange[1], steps)
        figure = plt.figure()
        figure.suptitle('Plot of bias as a function of hlr')
        id_params = copy.deepcopy(self.fish.g_parameters.id_params)
        snr = self.fish.snr
        for i in range(self.num_params):
            param = self.param_names[i]
            biases = []
            sigmas = []
            bias_sigmas = []
            if self.num_params % 2 == 0:
                ax = figure.add_subplot(2, self.num_params / 2, i + 1)
            else:
                ax = figure.add_subplot(2, self.num_params / 2 + 1, i + 1)
            for hlr in hlrs:
                gal_id = id_params.keys()[0]
                id_params[gal_id]['hlr'] = hlr
                g_parameters = galfun.GParameters(id_params=id_params)
                fish = fisher.Fisher(g_parameters, snr)
                bias = fish.biases[param]
                sigma = math.sqrt(fish.covariance_matrix[param, param])
                biases.append(bias)
                sigmas.append(sigma)
                bias_sigma = bias / sigma
                bias_sigmas.append(bias_sigma)

            if self.bias_sigma:
                ax.errorbar(hlrs, bias_sigmas)
            elif self.error_bars:
                ax.errorbar(hlrs, biases, yerr=sigmas)
            else:
                ax.errorbar(hlrs, biases)
            ax.set_title(param, fontsize=defaults.FONTSIZE_LABEL)

        saveFigureToPdf(figure, self.nameImages(),
                        self.project, self.plots_dir, hide=self.hide)

    def biasPlot2(self):
        """Plot of bias*(snr_nrom/snr)**2 as a function of snr with fixed size (hlr)."""
        steps = 10
        snr_range = (20, 60)
        snrs = np.linspace(snr_range[0], snr_range[1], steps)
        figure = plt.figure()
        figure.suptitle('Plot of bias*(snr)**2 as a function of snr')
        id_params = copy.deepcopy(self.fish.g_parameters.id_params)
        for i in range(self.num_params):
            param = self.param_names[i]
            ys = []
            if self.num_params % 2 == 0:
                ax = figure.add_subplot(2, self.num_params / 2, i + 1)
            else:
                ax = figure.add_subplot(2, self.num_params / 2 + 1, i + 1)
            for snr in snrs:
                g_parameters = galfun.GParameters(id_params=id_params)
                fish = fisher.Fisher(g_parameters, snr)
                bias = fish.biases[param]
                y = bias * (snr)**2
                ys.append(y)

            ax.errorbar(snrs, ys)
            ax.set_title(param, fontsize=defaults.FONTSIZE_LABEL)

        saveFigureToPdf(figure, self.nameImages(),
                        self.project, self.plots_dir, hide=self.hide)

    def biasPlot3(self):
        """Plot of bias*(snr_norm/snr)**2 as a function of hlr/psf_fwhm with a fixed snr."""
        fancy_plot = True
        if fancy_plot:
            print 'Fancy plotting something'

        steps = 25
        x_range = (.2, 1.5)  # x = hlr_gal / psf_fwhm
        xs = np.linspace(x_range[0], x_range[1], steps)
        figure = plt.figure(figsize=(25, 25))
        id_params = copy.deepcopy(self.fish.g_parameters.id_params)
        snr = self.fish.snr #normally use 20.
        figure.suptitle('Plot of ' + r'$b(a_{i})' + r'\left(' +
                        str(defaults.SNR_NORM) + r'/snr\right)^{2}$' +
                        'as a function of ' + r'hlr/fwhm$_{psf}$', fontsize=20)
        ys = {} # y= bias*snr2
        for x in xs:
            gal_id = id_params.keys()[0]
            hlr = x * id_params[gal_id]['psf_fwhm']
            id_params[gal_id]['hlr'] = hlr
            g_parameters = galfun.GParameters(id_params=id_params)
            fish = fisher.Fisher(g_parameters, snr)
            biases = fish.biases
            for i in range(self.num_params):
                param = self.param_names[i]
                if param not in ys:
                    ys[param] = []
                bias = biases[param]
                y = bias
                ys[param].append(y)


        if fancy_plot:
            #do not plot x0 and y0
            extra = r'$\left(' + str(defaults.SNR_NORM) + r'/snr\right)^{2}$'
            y_names = [r'$b(x_{0})$', r'$b(y_{0})$', r'$b(flux)$',
                        r'$b(hlr)$', r'$b(e_{1})$',
                        r'$b(e_{2})$']
            y_titles = [elt + extra for elt in y_names]
            x_titles = [r'hlr/fwhm$_{psf}$'] * 6
            for i in range(self.num_params):
                param = self.param_names[i]
                ax = figure.add_subplot(2, self.num_params / 2,
                                        i + 1)
                ax.scatter(xs, ys[param])
                ax.errorbar(xs, ys[param])
                #ax.set_title(titles[i], fontsize=14)
                ax.set_xlabel(x_titles[i], fontsize=22)
                ax.set_ylabel(y_titles[i], fontsize=22)
                ax.tick_params(labelsize=20)
                #scientic notation force.
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        else:
            for i in range(self.num_params):
                param = self.param_names[i]
                ax = figure.add_subplot(2, self.num_params / 2,
                                        i + 1)
                ax.scatter(xs, ys[param])
                ax.errorbar(xs, ys[param])
                ax.set_title(param, fontsize=14)


        saveFigureToPdf(figure, self.nameImages(),
                        self.project, self.plots_dir, hide=self.hide)

    def biasPlot4(self):
        pass

    def biasPlot5(self):
        separation_range = (.0, 3.0)  # arcsecs.
        steps = 50
        figure = plt.figure(figsize=(50, 50))
        figure.suptitle('Bias as a function of distance between two galaxies')
        id_params = copy.deepcopy(self.fish.g_parameters.id_params)
        snr = self.fish.snr
        separations = np.linspace(separation_range[0], separation_range[1],
                                  steps)
        biases = {}
        sigmas = {}
        bias_sigmas = {}
        for separation in separations:
            pos1 = separation / 2
            pos2 = -separation / 2
            ids = id_params.keys()
            id1 = ids[0]
            id2 = ids[1]
            id_params[id1]['x0'] = pos1
            id_params[id2]['x0'] = pos2
            g_parameters = galfun.GParameters(id_params=id_params)
            fish = fisher.Fisher(g_parameters, snr)
            for i in range(self.num_params):
                param = self.param_names[i]
                if param not in biases:
                    biases[param] = []
                if param not in sigmas:
                    sigmas[param] = []
                if param not in bias_sigmas:
                    bias_sigmas[param] = []
                bias = fish.biases[param]
                sigma = math.sqrt(fish.covariance_matrix[param, param])
                biases[param].append(bias)
                sigmas[param].append(sigma)
                bias_sigmas[param].append(bias / sigma)

        if defaults.FANCY:
            y_titles = [r'$b(x^{1}_{0})$', r'$b(y^{1}_{0})$', r'$b(flux^{1})$',
                       r'$b(hlr^{1})$', r'$b(e^{1}_{1})$',
                       r'$b(e^{1}_{2})$',
                       r'$b(x^{2}_{0})$', r'$b(y^{2}_{0})$', r'$b(flux^{2})$',
                       r'$b(hlr^{2})$', r'$b(e^{2}_{1})$',
                       r'$b(e^{2}_{2})$']
            x_titles = [r'$\left|x^{2}_{0} - x^{1}_{0}\right|$'] * 12

            for i in range(self.num_params):
                param = self.param_names[i]
                ax = figure.add_subplot(2, self.num_params / 2,
                                        i + 1)
                ax.scatter(separations, biases[param])
                ax.errorbar(separations, biases[param])
                #ax.set_title(titles[i], fontsize=14)
                ax.set_xlabel(x_titles[i], fontsize=30)
                ax.set_ylabel(y_titles[i], fontsize=30)
                ax.tick_params(labelsize=22)
                #scientic notation force.
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        else:
            for i in range(self.num_params):
                param = self.param_names[i]
                ax = figure.add_subplot(2, self.num_params / 2, i + 1)
                if self.bias_sigma:
                    ax.errorbar(separations, bias_sigmas[param])
                elif self.error_bars:
                    ax.errorbar(separations, biases[param], yerr=sigmas[param])
                else:
                    ax.scatter(separations, biases[param]) #dots
                    ax.errorbar(separations, biases[param]) #lines

                ax.set_title(param, fontsize=defaults.FONTSIZE_LABEL)

        if defaults.FANCY:
            figure.subplots_adjust(wspace=.3)

        saveFigureToPdf(figure, self.nameImages(),
                        self.project, self.plots_dir, hide=self.hide)

    def biasPlot6(self):
        """Must parametrize at least one galaxy with e/beta. Assume also
        that at least one of the galaxies is non-circular (first galaxy
        given has e1,e2 different from zero).
        """
        angle_range = (0, math.pi)
        steps = 30
        figure = plt.figure()
        figure.suptitle('Bias as a function of beta of a galaxy.')
        id_params = copy.deepcopy(self.fish.g_parameters.id_params)
        snr = self.fish.snr
        angles = np.linspace(angle_range[0], angle_range[1], steps)
        biases = {}
        sigmas = {}
        bias_sigmas = {}
        for angle in angles:
            ids = id_params.keys()
            id1 = ids[0]
            id_params[id1]['beta'] = angle
            g_parameters = galfun.GParameters(id_params=id_params)
            fish = fisher.Fisher(g_parameters, snr)
            for i in range(self.num_params):
                param = self.param_names[i]

                if param not in biases:
                    biases[param] = []
                if param not in sigmas:
                    sigmas[param] = []
                if param not in bias_sigmas:
                    bias_sigmas[param] = []

                bias = fish.biases[param]
                sigma = math.sqrt(fish.covariance_matrix[param, param])
                biases[param].append(bias)
                sigmas[param].append(sigma)
                bias_sigmas[param].append(bias / sigma)

        for i in range(self.num_params):
            param = self.param_names[i]
            ax = figure.add_subplot(2, self.num_params / 2, i + 1)
            if self.bias_sigma:
                ax.errorbar(angles, bias_sigmas[param])
            elif self.error_bars:
                ax.errorbar(angles, biases[param], yerr=sigmas[param])
            else:
                ax.scatter(angles, biases[param])
                ax.errorbar(angles, biases[param])

            ax.set_title(param, fontsize=defaults.FONTSIZE_LABEL)

        saveFigureToPdf(figure, self.nameImages(),
                        self.project, self.plots_dir, hide=self.hide)

    def biasPlot7(self):
        """Parametrize both galaxies with e/beta"""
        vsteps = 3
        hsteps = 3
        hrange = (0, math.pi)
        vrange = (0, math.pi)
        figure = plt.figure()
        figure.suptitle('Bias as a function two galaxies angular position.')
        id_params = copy.deepcopy(self.fish.g_parameters.id_params)
        ids = id_params.keys()
        id1 = ids[0]
        id2 = ids[1]
        snr = self.fish.snr
        biases = {}
        hangles = np.linspace(hrange[0], hrange[1], hsteps)
        vangles = np.linspace(vrange[0], vrange[1], vsteps)
        for hangle in hangles:
            for vangle in vangles:
                id_params[id1]['beta'] = hangle
                id_params[id2]['beta'] = vangle
                g_parameters = galfun.GParameters(id_params=id_params)
                fish = fisher.Fisher(g_parameters, snr)
                for k in range(self.num_params):
                    param = self.param_names[k]

                    if param not in biases:
                        biases[param] = {}
                    if hangle not in biases:
                        biases[param][hangle] = []

                    biases[param][hangle].append(fish.biases[param])

        figure = plt.figure()
        #transform each hangle list to a matrix of biases and plot the matrix.
        for k in range(self.num_params):
            param = self.param_names[k]
            matrix = []
            for hangle in biases[param]:
                row = biases[param][hangle]
                matrix.append(row)
            plot = np.array(matrix)
            ax = figure.add_subplot(2, self.num_params / 2, k + 1)
            ax.imshow(plot, interpolation='none', rasterized=False,
                      cmap=cm.RdYlGn, origin='lower',
                      extent=[hrange[0], hrange[1], vrange[0], vrange[1]],
                      vmax=abs(plot).max(), vmin=-abs(plot).max())

        saveFigureToPdf(figure, self.nameImages(),
                                self.project, self.plots_dir, hide=self.hide)

        def plotConditionNumber1(self):
            """Plot condition number of fisher matrix as a function of
            separation between the galaxies.
            """

            separation_range = (.0, 3.0)  # arcsecs.
            steps = 50
            figure = plt.figure(figsize=(50, 50))
            figure.suptitle('Condition Number as a function of distance'
                            'between two galaxies')
            id_params = copy.deepcopy(self.fish.g_parameters.id_params)
            snr = self.fish.snr
            separations = np.linspace(separation_range[0], separation_range[1],
                                      steps)
            condition_numbers = {}
            for separation in separations:
                pos1 = separation / 2
                pos2 = -separation / 2
                ids = id_params.keys()
                id1 = ids[0]
                id2 = ids[1]
                id_params[id1]['x0'] = pos1
                id_params[id2]['x0'] = pos2
                g_parameters = galfun.GParameters(id_params=id_params)
                fish = fisher.Fisher(g_parameters, snr)
                for i in range(self.num_params):
                    param = self.param_names[i]
                    if param not in values:
                        biases[param] = []


                for i in range(self.num_params):
                    param = self.param_names[i]
                    ax = figure.add_subplot(2, self.num_params / 2,
                                            i + 1)
                    ax.scatter(separations, biases[param])
                    ax.errorbar(separations, biases[param])
                    #ax.set_title(titles[i], fontsize=14)
                    ax.set_xlabel(x_titles[i], fontsize=30)
                    ax.set_ylabel(y_titles[i], fontsize=30)
                    ax.tick_params(labelsize=22)
                    #scientic notation force.
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,
                    0))


    def noisy_image(self):
        """Visualization of noise on the given galaxy"""
        figure, subplt = plt.subplots(1, 1)
        figure.suptitle('Noisy Initial Galaxy',
                        fontsize=defaults.FONTSIZE_TITLE)
        noise_image, _ = galfun.addNoise(self.gal_image, self.fish.snr)
        drawImage(subplt, noise_image.array)


        saveFigureToPdf(figure, self.nameImages(),
                        self.project, self.plots_dir, hide=self.hide)

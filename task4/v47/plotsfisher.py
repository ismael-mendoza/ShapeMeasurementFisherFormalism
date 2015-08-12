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

plt.rc('text', usetex=False) #ignore latex commands. 


def annotateAxisCenter(ax, text):
    ax.text(.5,.5, text,
    ha='center', va ='center', transform=ax.transAxes, 
    fontweight='bold', fontsize = defaults.FONTSIZE_VALUE)


def saveFigureToPdf(figure,file_name, project, plotsdir='', hide = False): 

     #save and preview pdf.
    if not os.path.isdir(project):
        os.mkdir(project)

    if not os.path.isdir(os.path.join(project,plotsdir)):
        os.mkdir(os.path.join(project,plotsdir))

    file_name = os.path.join(project, plotsdir, file_name) 
    figure.savefig(file_name, bbox_inches='tight', dpi=defaults.DPI)

    if not hide: 
        os.system("open " + file_name)

    #close it so does not consume memory.
    plt.close(figure)


def drawImage(ax, plot, title = "", xlabel = "", ylabel = ""): 
    """draws a given plot with default values using imshow in the given 
    subplot if they are not subplots should just pass plt. It also adds a 
    title with name. -- OWN"""
    
    #remove numbers and ticks but can use label.
    ax.axes.get_xaxis().set_ticks([]) 
    ax.axes.get_yaxis().set_ticks([])

    #labels figure. 
    ax.set_title(title, fontsize = defaults.FONTSIZE_LABEL)
    ax.set_xlabel(xlabel, fontsize=defaults.FONTSIZE_LABEL)
    ax.set_ylabel(ylabel, fontsize=defaults.FONTSIZE_LABEL)

    ax.imshow(plot, interpolation='nearest', rasterized=False, 
                cmap=cm.RdYlGn, origin='lower', 
                extent=[-defaults.NX,defaults.NX,-defaults.NY,defaults.NY],
                vmax=abs(plot).max(), vmin=-abs(plot).max())


class FisherPlots:
    """Produce plots for a fish object and displays them, in a given 
    plots_dir that is in a given project. Hide if just save images but not 
    display them.
    """
    def __init__(self, fish, project, plots_dir, hide): 
        self.fish = fish
        self.project = project
        self.plots_dir = plots_dir
        self.num_params = fish.num_params
        self.param_names = fish.param_names
        self.gal_image = fish.image.array
        self.hide = hide

        #delete plots_dir directory for nameImages to work if already there.
        if os.path.isdir(os.path.join(self.project,self.plots_dir)): 
            os.system('rm -r ' + os.path.join(self.project,self.plots_dir))
        else: 
            os.mkdir(os.path.join(self.project,self.plots_dir))

    def nameImages(self):
        """Returns first base#n.extension string where n is an integer 
        that is not a file in project/pltsdir.
        """
        n = 0
        while(True):
            n+= 1
            filename = ''.join([defaults.FIGURE_BASENAME,str(n),
                                defaults.FIGURE_EXTENSION])
            path = os.path.join(self.project, self.plots_dir,filename)
            figure_name_exists = os.path.isfile(path)
            if not figure_name_exists:
                break

        return filename

    def galaxy(self):
        figure, subplt= plt.subplots(1,1)
        figure.suptitle('Initial Galaxy', fontsize=defaults.FONTSIZE_TITLE)
        drawImage(subplt, self.gal_image)

        saveFigureToPdf(figure, self.nameImages(), 
            self.project, self.plots_dir, hide=self.hide)

    def derivatives(self):
        figure = plt.figure() 
        figure.suptitle(
                       'Derivatives of model with respect to each parameter', 
                        fontsize = defaults.FONTSIZE_TITLE)

        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawImage(ax, self.fish.derivatives_images[self.param_names[i]],
                      title=self.param_names[i])

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, 
                        hide=self.hide)

    def fisherMatrix(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements', 
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params): 
                for j in range(self.num_params):
                    if(i >= j):
                        ax = figure.add_subplot(self.num_params,
                                                self.num_params, 
                                                self.num_params * i + j + 1)
                        drawImage(ax, 
                                  self.fish.fisher_matrix_images[
                                  self.param_names[i], 
                                  self.param_names[j]
                                  ])

                        if(j == 0):
                            ax.set_ylabel(self.param_names[i])
                        if(i == self.num_params - 1):
                            ax.set_xlabel(self.param_names[j])

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, 
                       hide=self.hide)

    def fisherMatrixValues(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values ', 
                        fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,
                                            self.num_params, 
                                            self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[
                              self.param_names[i],self.param_names[j]])

                    if j == 0:
                        ax.set_ylabel(self.param_names[i])

                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])

                    annotateAxisCenter(ax, str(round(self.fish.fisher_matrix[
                        self.param_names[i],self.param_names[j]],defaults.SIG_DIGITS)))

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, 
                        hide = self.hide)

    def fisherMatrixChi2(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values of chi2 method ', fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[self.param_names[i],self.param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
                    annotateAxisCenter(ax, str(round(self.fish.fisher_matrix_chi2[self.param_names[i],self.param_names[j]],defaults.SIG_DIGITS)))
        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, hide = self.hide)

    def secondDerivatives(self):
        figure = plt.figure()
        figure.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.second_derivatives_images[self.param_names[i],self.param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, hide = self.hide)

    def biasMatrixValues(self):
        figuresOfBiasMatrixNumbers = []    
        for i in range(self.num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + self.param_names[i] + ')' , fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
            for j in range(self.num_params): 
                for k in range(self.num_params):
                    if(j >= k):
                        ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * j + k + 1)
                        drawImage(ax, self.fish.bias_matrix_images[self.param_names[i],self.param_names[j],self.param_names[k]]) 
                        if(k == 0):
                            ax.set_ylabel(self.param_names[j])
                        if(j == self.num_params - 1):
                            ax.set_xlabel(self.param_names[k])

                        annotateAxisCenter(ax, str(round(self.fish.bias_matrix
                                [self.param_names[i],
                                self.param_names[j],
                                self.param_names[k]],defaults.SIG_DIGITS)))
            figuresOfBiasMatrixNumbers.append(figure)

        for i, figure in enumerate(figuresOfBiasMatrixNumbers):
            basename = self.nameImages().strip(defaults.FIGURE_EXTENSION)
            filename = ''.join([basename,'_',
                                str(i),defaults.FIGURE_EXTENSION])
            saveFigureToPdf(figure, filename, self.project, self.plots_dir, 
                            hide = self.hide) 

    def biasMatrix(self):
        figuresOfBiasMatrix = [] 
        for i in range(self.num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + self.param_names[i] + ')' , fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
            for j in range(self.num_params): 
                for k in range(self.num_params):
                    if(j >= k):
                        ax = figure.add_subplot(self.num_params,self.num_params,self.num_params * j + k + 1)
                        drawImage(ax, self.fish.bias_matrix_images[self.param_names[i],self.param_names[j],self.param_names[k]]) 
                        if(k == 0):
                            ax.set_ylabel(self.param_names[j])
                        if(j == self.num_params - 1):
                            ax.set_xlabel(self.param_names[k])
            figuresOfBiasMatrix.append(figure)

        for i, figure in enumerate(figuresOfBiasMatrix): 
            basename = self.nameImages().strip(defaults.FIGURE_EXTENSION)
            filename = ''.join([basename,'_',str(i),defaults.FIGURE_EXTENSION])
            saveFigureToPdf(figure, filename, self.project, self.plots_dir, 
                            hide = self.hide)  

    def covarianceMatrix(self):
        figure = plt.figure()
        figure.suptitle('Covariance matrix elements', fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[self.param_names[i],self.param_names[j]] * 0) #figure out better way.
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
                    annotateAxisCenter(ax, str(round(self.fish.covariance_matrix[self.param_names[i],self.param_names[j]],defaults.SIG_DIGITS)))

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, 
                        hide = self.hide)

    def correlationMatrix(self):
        figure = plt.figure()
        figure.suptitle('Correlation matrix elements', fontsize=defaults.FONTSIZE_TITLE, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[self.param_names[i],self.param_names[j]] * 0)

                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])

                    annotateAxisCenter(ax, str(round(self.fish.correlation_matrix[self.param_names[i],self.param_names[j]],defaults.SIG_DIGITS)))

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, 
                        hide = self.hide)

    def biasValues(self):
        figure = plt.figure() 
        figure.suptitle('Contribution to Bias from each pixel by parameter', 
                        fontsize = defaults.FONTSIZE_TITLE)
        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawImage(ax, self.fish.bias_images[self.param_names[i]], 
                      title = self.param_names[i])

            annotateAxisCenter(ax, str(round(self.fish.biases[self.param_names
                                                              [i]],defaults.SIG_DIGITS)))

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, 
                        hide = self.hide)

    def biases(self):
        figure = plt.figure() 
        figure.suptitle('Contribution to Bias from each pixel by parameter', 
                        fontsize = defaults.FONTSIZE_TITLE)
        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawImage(ax, self.fish.bias_images[self.param_names[i]], 
                      title = self.param_names[i])

        saveFigureToPdf(figure, self.nameImages(), self.project, self.plots_dir, 
                        hide = self.hide)



    #plot of two identical circular galaxies and their bias as a function of 
    #distance. 

    #user writes the 2 galaxies and only two galaxies can be accepted, place check in display.py to avoid errors. 
    #user writes range of separations to observe. 
    #in default we write the step in the range (or also given by the user?)
    #maybe need to write methods for changing g_params because fisher only 
    #nreads g_params, do not know how necessary...
    #plot
    
    # def biasPlot1(self):
    #     extension = (.9,1.2) #arcsecs
    #     steps = 4 #how many step from initial separation to final separation.
    #     figure = plt.figure()
    #     figure.suptitle('Bias as a function of distance between two galaxies')
    #     id_params = copy.deepcopy(self.fish.g_parameters.id_params)
    #     snr = self.fish.snr
    #     for i in range(self.num_params):
    #         param = self.param_names[i]
    #         biases = []
    #         stds = []
    #         ax = figure.add_subplot(1,self.num_params,i+1)
    #         separations = np.linspace(extension[0],extension[1],steps)
    #         for separation in separations:
    #             pos1 = separation/2
    #             pos2 = -separation/2
    #             ids = id_params.keys()
    #             id1 = ids[0]
    #             id2 =ids[1]
    #             id_params[id1]['x0'] = pos1
    #             id_params[id2]['x0'] = pos2
    #             g_parameters = galfun.GParameters(id_params=id_params)
    #             fish = fisher.Fisher(g_parameters, snr)
    #             biases.append(fish.biases[param])
    #             stds.append(math.sqrt(fish.covariance_matrix[param][param]))


    #         ax.errorbar(separations,biases,xerr=stds)
    #         ax.set_title(param, fontsize = defaults.FONTSIZE_LABEL)

    # def biasPlot2(self):
    #     """Must parametrize at least one galaxy with e/beta. Assume also 
    #     that at least one of the galaxies is non-circular (first galaxy 
    #     given). 

    #     """
    #     #distance and snr specified by user.
    #     extension=(0,math.pi)
    #     steps=10         
    #     figure = plt.figure()
    #     figure.suptitle('Bias as a function of beta of a galaxy.')
    #     id_params = copy.deepcopy(self.fish.g_parameters.id_params)
    #     snr = self.fish.snr
    #     for i in range(self.num_params):
    #         param = self.param_names[i]
    #         biases = []
    #         stds = []
    #         ax = figure.add_subplot(1,self.num_params,i+1)
    #         angles = np.linspace(extension[0],extension[1],steps)
    #         for angle in angles:
    #             ids = id_params.keys()
    #             id1 = ids[0]
    #             id_params[id1]['beta'] = angle
    #             g_parameters = galfun.GParameters(id_params=id_params)
    #             fish = fisher.Fisher(g_parameters, snr)
    #             biases.append(fish.biases[param])
    #             stds.append(math.sqrt(fish.covariance_matrix[param][param]))

    #         ax.errorbar(angles,biases,xerr=stds)
    #         ax.set_title(param, fontsize = defaults.FONTSIZE_LABEL)

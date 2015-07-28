#!/usr/bin/env python
import matplotlib.pyplot as plt

import matplotlib.cm as cm

import defaults

import os

plt.rc('text', usetex=False) #ignore latex commands. 
cts = defaults.constants()


def saveFigureToPdf(figure,file_name, wdir, plotsdir, hide = False): 

     #save and preview pdf.
    if not os.path.isdir(wdir):
        os.mkdir(wdir)

    if not os.path.isdir(os.path.join(wdir,plotsdir)):
        os.mkdir(os.path.join(wdir,plotsdir))

    file_name = os.path.join(wdir, plotsdir, file_name) 
    figure.savefig(file_name, bbox_inches='tight', dpi = cts.dpi)

    if not hide: 
        os.system("open " + file_name)


def drawImage(axis, plot, title = "", xlabel = "", ylabel = ""): 
    """draws a given plot with default values using imshow in the given 
    subplot if they are not subplots should just pass plt. It also adds a 
    title with name. -- OWN"""
    
    #remove numbers and ticks but can use label.
    axis.axes.get_xaxis().set_ticks([]) 
    axis.axes.get_yaxis().set_ticks([])

    #labels figure. 
    axis.set_title(title, fontsize = cts.fontsize_label)
    axis.set_xlabel(xlabel, fontsize=cts.fontsize_label)
    axis.set_ylabel(ylabel, fontsize=cts.fontsize_label)

    axis.imshow(plot, interpolation='nearest', rasterized=False, 
                cmap=cm.RdYlGn, origin='lower', 
                extent=[-cts.nx,cts.nx,-cts.ny,cts.ny],
                vmax=abs(plot).max(), vmin=-abs(plot).max())


class FisherPlots:
    """Produce plots for a fish object and displays them, in a given 
    plots_dir that is in a given wdir. Hide if just save images but not 
    display them.
    """
    def __init__(self, fish, wdir, plots_dir, hide): 
        self.fish = fish
        self.wdir = wdir
        self.plots_dir = plots_dir
        self.num_params = fish.num_params
        self.param_names = fish.param_names
        self.gal_image = fish.image.array
        self.hide = hide

        #delete plots_dir directory for nameImages to work if already there.
        if os.path.isdir(os.path.join(self.wdir,self.plots_dir)): 
            os.system('rm -r ' + os.path.join(self.wdir,self.plots_dir))
        else: 
            os.mkdir(os.path.join(self.wdir,self.plots_dir))

    def nameImages(self):
        """Returns first base#n.extension string where n is an integer 
        that is not a file in wdir/pltsdir.
        """
        n = 0
        while(True):
            n+= 1
            filename = ''.join([cts.figure_basename,str(n),
                                cts.figure_extension])
            path = os.path.join(self.wdir, self.plots_dir,filename)
            figure_name_exists = os.path.isfile(path)
            if not figure_name_exists:
                break

        return filename

    def galaxy(self):
        figure, subplt= plt.subplots(1,1)
        figure.suptitle('Initial Galaxy', fontsize=cts.fontsize_titles)
        drawImage(subplt, self.gal_image)

        saveFigureToPdf(figure, self.nameImages(), 
            self.wdir, self.plots_dir, hide=self.hide)

    def derivatives(self):
        figure = plt.figure() 
        figure.suptitle(
                       'Derivatives of model with respect to each parameter', 
                        fontsize = cts.fontsize_titles)

        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawImage(ax, self.fish.derivatives_images[self.param_names[i]],
                      title=self.param_names[i])

        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, 
                        hide=self.hide)

    def fisherMatrix(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements', 
                        fontsize=cts.fontsize_titles, fontweight='bold')
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

        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, 
                       hide=self.hide)

    def fisherMatrixValues(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values ', 
                        fontsize=cts.fontsize_titles, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,
                                            self.num_params, 
                                            self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fish_matrix_images[
                              self.param_names[i],self.param_names[j]])

                    if j == 0:
                        ax.set_ylabel(self.param_names[i])

                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])

                    ax.text(-20,0, str(round(self.fish.fisher_matrix[
                        self.param_names[i],self.param_names[j]],5)), 
                        fontweight='bold', fontsize = cts.fontsize_values)

        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, 
                        hide = self.hide)

    def fisherMatrixChi2(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values of chi2 method ', fontsize=cts.fontsize_titles, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[self.param_names[i],self.param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
                    ax.text(-20,0, str(round(self.fish.fisher_matrix_chi2[self.param_names[i],self.param_names[j]],5)), 
                        fontweight='bold', fontsize = cts.fontsize_values)
        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, hide = self.hide)

    def secondDerivatives(self):
        figure = plt.figure()
        figure.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=cts.fontsize_titles, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.second_derivatives_images[self.param_names[i],self.param_names[j]])
                    # ax.text(0,20, "std:" + str(fish.second_derivatives_galaxy_images[self.param_names[i],self.param_names[j]].std().round(4)), fontsize = 8, fontweight='bold') 
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, hide = self.hide)

    def biasMatrixValues(self):
        figuresOfBiasMatrixNumbers = []    
        for i in range(self.num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + self.param_names[i] + ')' , fontsize=cts.fontsize_titles, fontweight='bold')
            for j in range(self.num_params): 
                for k in range(self.num_params):
                    if(j >= k):
                        ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * j + k + 1)
                        drawImage(ax, self.fish.bias_matrix_images[self.param_names[i],self.param_names[j],self.param_names[k]]) 
                        if(k == 0):
                            ax.set_ylabel(self.param_names[j])
                        if(j == self.num_params - 1):
                            ax.set_xlabel(self.param_names[k])
                        ax.text(-20,0, 
                            str(round(self.fish.bias_matrix
                                [self.param_names[i],
                                self.param_names[j],
                                self.param_names[k]],5)), 
                            fontweight='bold', fontsize = cts.fontsize_values)
            figuresOfBiasMatrixNumbers.append(figure)

        for i, figure in enumerate(figuresOfBiasMatrixNumbers):
            basename = self.nameImages().strip(cts.figure_extension)
            filename = ''.join([basename,'_',str(i),cts.figure_extension])
            saveFigureToPdf(figure, filename, self.wdir, self.plots_dir, 
                            hide = self.hide) 

    def biasMatrix(self):
        figuresOfBiasMatrix = [] 
        for i in range(self.num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + self.param_names[i] + ')' , fontsize=cts.fontsize_titles, fontweight='bold')
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
            basename = self.nameImages().strip(cts.figure_extension)
            filename = ''.join([basename,'_',str(i),cts.figure_extension])
            saveFigureToPdf(figure, filename, self.wdir, self.plots_dir, 
                            hide = self.hide)  

    def covarianceMatrix(self):
        figure = plt.figure()
        figure.suptitle('Covariance matrix elements', fontsize=cts.fontsize_titles, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fish_matrix_images[self.param_names[i],self.param_names[j]] * 0) #figure out better way.
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
                    ax.text(-20,0, str(round(self.fish.covariance_matrix[self.param_names[i],self.param_names[j]],5)), fontweight='bold', fontsize= cts.fontsize_values)

        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, 
                        hide = self.hide)

    def correlationMatrix(self):
        figure = plt.figure()
        figure.suptitle('Correlation matrix elements', fontsize=cts.fontsize_titles, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawImage(ax, self.fish.fisher_matrix_images[self.param_names[i],self.param_names[j]] * 0)

                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])

                    ax.text(-20,0, str(round(self.fish.correlation_matrix[self.param_names[i],self.param_names[j]],5)), fontweight='bold', fontsize = cts.fontsize_values)

        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, 
                        hide = self.hide)

    def biasValues(self):
        figure = plt.figure() 
        figure.suptitle('Contribution to Bias from each pixel by parameter', 
                        fontsize = cts.fontsize_titles)
        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawImage(ax, self.fish.bias_images[self.param_names[i]], 
                      title = self.param_names[i])
            ax.text(-20,0, 
                    str(round(self.fish.biases[self.param_names[i]],5)), 
                    fontweight='bold', fontsize = cts.fontsize_values)

        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, 
                        hide = self.hide)

    def biases(self):
        figure = plt.figure() 
        figure.suptitle('Contribution to Bias from each pixel by parameter', 
                        fontsize = cts.fontsize_titles)
        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawImage(ax, self.fish.bias_images[self.param_names[i]], 
                      title = self.param_names[i])

        saveFigureToPdf(figure, self.nameImages(), self.wdir, self.plots_dir, 
                        hide = self.hide)
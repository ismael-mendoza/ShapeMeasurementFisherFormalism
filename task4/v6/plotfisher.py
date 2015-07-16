import argparse
import defaults
import os
import matplotlib.pyplot as plt
from functions import *
import csv
import math
import numpy as np

plt.rc('text', usetex=False) #ignore latex commands. 
class fisherplots:
    """Produce plots for a fisher object and displays them, in a given plots_dir that is in a given wdir. Hide if just save images but not display them."""

    def __init__(self, fisher, wdir, plots_dir, hide): 
        self.fisher = fisher
        self.wdir = wdir
        self.plot_dir = plots_dir
        self.num_params = fisher.num_params
        self.param_names = fisher.param_names
        self.gal_image = fisher.gal_image
        self.hide = hide

    def plotGalaxy(self):
        figure, subplt= plt.subplots(1,1)
        figure.suptitle('Initial Galaxy', fontsize = 20)
        drawPlot(subplt, self.gal_image.array)
        SaveFigureToPdf(figure, 'figure1.png', self.wdir, self.plots_dir, hide = self.hide)

    def plotDerivatives(self):
        figure = plt.figure() 
        figure.suptitle('Derivatives of model with respect to each parameter', fontsize = 20)
        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawPlot(ax, self.fisher.derivatives_galaxy_images[self.param_names[i]], title = self.param_names[i])
        SaveFigureToPdf(figure, 'figure2.png', self.wdir, self.plots_dir, hide = self.hide)

    def fisherMatrix(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements ', fontsize=14, fontweight='bold')
        for i in range(self.num_params): 
                for j in range(self.num_params):
                    if(i >= j):
                        ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                        drawPlot(ax, self.fisher.fisher_matrix_images[self.param_names[i],self.param_names[j]])
                        if(j == 0):
                            ax.set_ylabel(self.param_names[i] )
                        if(i == num_params - 1):
                            ax.set_xlabel(self.param_names[j])
        SaveFigureToPdf(figure, 'figure3.png', self.wdir, self.plots_dir, hide = self.hide)

    def fisherMatrixValues(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values ', fontsize=14, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawPlot(ax, self.fisher.fisher_matrix_images[self.param_names[i],self.param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
                    ax.text(-20,0, str(round(self.fisher.fisher_matrix[param_names[i],param_names[j]],5)), fontweight='bold')
        SaveFigureToPdf(figure, 'figure4.png', self.wdir, self.plots_dir, hide = self.hide)

    def fisherMatrixChi2(self):
        figure = plt.figure()
        figure.suptitle('Fisher matrix elements with values of chi2 method ', fontsize=14, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawPlot(ax, self.fisher.fisher_matrix_images[selfparam_names[i],param_names[j]])
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
                    ax.text(-20,0, str(round(self.fisher.fisher_matrix_chi2[self.param_names[i],self.param_names[j]],5)), fontweight='bold')
        SaveFigureToPdf(figure, 'figure5.png', self.wdir, self.plots_dir, hide = self.hide)

    def SecondDerivativesGalaxy(self):
        figure = plt.figure()
        figure.suptitle('Second derivatives of the galaxies with respect to parameters', fontsize=14, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawPlot(ax, self.fisher.second_derivatives_galaxy_images[self.param_names[i],self.param_names[j]])
                    # ax.text(0,20, "std:" + str(fisher.second_derivatives_galaxy_images[self.param_names[i],self.param_names[j]].std().round(4)), fontsize = 8, fontweight='bold') 
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
        SaveFigureToPdf(figure, 'figure8.png', self.wdir, self.plots_dir, hide = self.hide)

    def biasMatrixValues(self):
        figuresOfBiasMatrixNumbers = []    
        for i in range(self.num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
            for j in range(self.num_params): 
                for k in range(self.num_params):
                    if(j >= k):
                        ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * j + k + 1)
                        drawPlot(ax, self.fisher.bias_matrix_images[self.param_names[i],self.param_names[j],self.param_names[k]]) 
                        if(k == 0):
                            ax.set_ylabel(self.param_names[j])
                        if(j == self.num_params - 1):
                            ax.set_xlabel(self.param_names[k])
                        ax.text(-20,0, str(round(self.fisher.bias_matrix[self.param_names[i],self.param_names[j],self.param_names[k]],5)), fontweight='bold')
            figuresOfBiasMatrixNumbers.append(figure)

        for i, figure in enumerate(figuresOfBiasMatrixNumbers): 
            SaveFigureToPdf(figure, 'figure' + str(10) + '_' + str(i) + '.png', self.wdir, self.plots_dir, hide = self.hide) 

    def biasMatrix(self):
        figuresOfBiasMatrix = [] 
        for i in range(self.num_params):
            figure = plt.figure()
            figure.suptitle('Bias matrix elements j,k for ' + str(i + 1) + 'th partial (D' + param_names[i] + ')' , fontsize=14, fontweight='bold')
            for j in range(self.num_params): 
                for k in range(self.num_params):
                    if(j >= k):
                        ax = figure.add_subplot(self.num_params,self.num_params,self.num_params * j + k + 1)
                        drawPlot(ax, self.fisher.bias_matrix_images[self.param_names[i],self.param_names[j],self.param_names[k]]) 
                        if(k == 0):
                            ax.set_ylabel(self.param_names[j])
                        if(j == self.num_params - 1):
                            ax.set_xlabel(self.param_names[k])
            figuresOfBiasMatrix.append(figure)
        for i, figure in enumerate(figuresOfBiasMatrix): 
            SaveFigureToPdf(figure, 'figure' + str(9) + '_' + str(i) + '.png', self.wdir, self.plots_dir, hide = self.hide) 

    def covarianceMatrix(self):
        figure = plt.figure()
        figure.suptitle('Covariance matrix elements', fontsize=14, fontweight='bold')
        for i in range(self.num_params): 
            for j in range(self.num_params):
                if(i >= j):
                    ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawPlot(ax, self.fisher.fisher_matrix_images[self.param_names[i],self.param_names[j]] * 0) #figure out better way.
                    if(j == 0):
                        ax.set_ylabel(self.param_names[i])
                    if(i == self.num_params - 1):
                        ax.set_xlabel(self.param_names[j])
                    ax.text(-20,0, str(round(self.fisher.covariance_matrix[self.param_names[i],self.param_names[j]],5)), fontweight='bold')
        SaveFigureToPdf(figure, 'figure6.png', self.wdir, self.plots_dir, hide = self.hide)

    def correlationMatrix(self):
        figure = plt.figure()
        figure.suptitle('Correlation matrix elements', fontsize=14, fontweight='bold')
        for i in range(num_params): 
            for j in range(num_params):
                if(i >= j):
                    ax = figure6.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                    drawPlot(ax, self.fisher.fisher_matrix_images[self.param_names[i],self.param_names[j]] * 0)
                    if(j == 0):
                        ax.set_ylabel(param_names[i])
                    if(i == num_params - 1):
                        ax.set_xlabel(param_names[j])
                    ax.text(-20,0, str(round(self.fisher.correlation_matrix[self.param_names[i],self.param_names[j]],5)), fontweight='bold')
        SaveFigureToPdf(figure, 'figure7.png', self.wdir, self.plots_dir, hide = self.hide)

    def biasesValues(self):
        figure = plt.figure() 
        figure.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 14)
        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawPlot(ax, self.fisher.bias_images[self.param_names[i]], title = self.param_names[i])
            ax.text(-20,0, str(round(self.fisher.biases[self.param_names[i]],5)), fontweight='bold')

        SaveFigureToPdf(figure, 'figure12.png', self.wdir, self.plots_dir, hide = self.hide)


    def biases(self):
        figure = plt.figure() 
        figure.suptitle('Contribution to Bias from each pixel by parameter', fontsize = 15)
        for i in range(self.num_params):
            ax = figure.add_subplot(1,self.num_params,i+1)
            drawPlot(ax, self.fisher.bias_images[self.param_names[i]], title = self.param_names[i])

        SaveFigureToPdf(figure, 'figure11.png', self.wdir, self.plots_dir, hide = self.hide)

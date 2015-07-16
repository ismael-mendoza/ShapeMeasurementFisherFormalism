import argparse
import defaults
import os
import matplotlib.pyplot as plt
from functions import *
import csv
import math
import numpy as np



class fisher:
    """Given a galaxy image and the appropiate parameters that describe it, will produce a fisher analysis of it"""

    def __init__(self, params, gal_image, sigma_n): 
        names = defaults.names()
        self.gal_image = gal_image 
        self.params = params 
        self.sigma_n = sigma_n
        self.param_names = names.galaxy_parameters[params['model']]
        self.num_params = len(self.param_names)

        self.derivatives_galaxy = self.derivativesGalaxy()
        self.fisher_matrix_images = self.fisherMatrixImages()
        self.fisher_matrix = self.fisherMatrix()
        self.second_derivatives_galaxy = self.secondDerivativesGalaxy()
        self.bias_matrix_images = self.biasMatrixImages()
        self.biases = self.biasFunc()

    def derivativesGalaxy(self):
        steps = defaults.steps(params).dct
        return {self.param_names[i]:partialDifferentiate(func = drawGalaxy, parameter = self.param_names[i], step = steps[self.param_names[i]])(params).array for i in range(self.num_params)}

    def fisherMatrixImages(self):
        FisherM_images = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                FisherM_images[self.param_names[i],self.param_names[j]] = (self.derivatives_galaxy[self.param_names[i]] * self.derivatives_galaxy[self.param_names[j]]) /(self.sigma_n**2)
        return FisherM_images

    def fisherMatrix(self):
        FisherM = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                FisherM[self.param_names[i],self.param_names[j]] = fisherMatrixImages(self.gal_image,params)[self.param_names[i],self.param_names[j]].sum() #sum over all pixels. 
        return FisherM

    def fisherMatrixChi2(self):
        FisherM_chi2 = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                FisherM_chi2[self.param_names[i],self.param_names[j]] = .5 * secondPartialDifferentiate(chi2, self.param_names[i], self.param_names[j], steps[self.param_names[i]], steps[self.param_names[j]], sigma_n = self.
                    sigma_n, gal_image = self.gal_image)(params)

        return FisherM_chi2

    def secondDerivativesGalaxy(self):
        secondDs_gal = {}
        for i in range(self.num_params): 
            for j in range(self.num_params):
                SecondDs_gal[self.param_names[i],self.param_names[j]] = (secondPartialDifferentiate(drawGalaxy, self.param_names[i], self.param_names[j], steps[self.param_names[i]], steps[self.param_names[j]])(params).array)

        return secondDs_gal

    def biasMatrixImages(self):
        BiasM_images = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                for k in range(self.num_params): 
                    BiasM_images[self.param_names[i],self.param_names[j],self.param_names[k]] = (self.derivatives_galaxy[self.param_names[i]] * self.second_derivatives_galaxy[self.param_names[j],self.param_names[k]]) / (self.sigma_n**2)

        return BiasM_images

    def biasMatrix(self):
        BiasM = {}
        for i in range(self.num_params):
            for  j in range(self.num_params):
                for k in range(self.num_params):
                    BiasM[self.param_names[i],self.param_names[j],self.param_names[k]] = BiasM_images[self.param_names[i],self.param_names[j],self.param_names[k]].sum()

        return BiasM

    def covarianceMatrix(self):
        CovM = {}
        FisherM_array = np.array([[self.fisher_matrix[self.param_names[i],self.param_names[j]] for i in range(self.num_params)] for j in range(self.num_params)])#convert to numpy array because it can be useful.
        CovM_array = np.linalg.inv(FisherM_array) #need to be an array to inverse.

        for i in range(self.num_params):
            for j in range(self.num_params):
                CovM[self.param_names[i],self.param_names[j]] = CovM_array[i][j]

        return CovM

    def biasImages(self):
        bias_images = {}
        for i in range(self.num_params):
            sumation = 0 
            for j in range(self.num_params):
                for k in range(self.num_params):
                    for l in range(self.num_params):
                        sumation += CovM[self.param_names[i],self.param_names[j]]*CovM[self.param_names[k],self.param_names[l]]*BiasM_images[self.param_names[j],self.param_names[k],self.param_names[l]]
            bias_images[self.param_names[i]] = (-.5) * sumation

        return bias_images


    def biasesFunc(self):
        return {self.param_names[i]:self.bias_images[self.param_names[i]].sum() for i in range(self.num_params)}







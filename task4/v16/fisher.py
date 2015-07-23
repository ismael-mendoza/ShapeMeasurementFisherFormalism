#!/usr/bin/env python
"""This module contains functions necessary to produce statistical results of the fisher formalism from a given galaxy"""

import defaults

import functions as fns

import math

import numpy as np



class fisher:
    """Given a galaxy image and the appropiate parameters that describe it, will produce a fisher objec that contains the analysis of it using the fisher formalism"""

    def __init__(self, gals_params, gals_image, sigma_n): 
        names = defaults.names()
        self.gals_image = gals_image 
        self.gals_params = gals_params 
        self.steps = defaults.steps(self.params).dict
        self.sigma_n = float(sigma_n) #this is the jiggle in one pixel due to the noise, uniform for all pixels for now too.
        self.param_names = [self.gals_params[gal_id].keys() for gal_id in self.gals_params.keys()]
        #names.galaxy_parameters[params['model']]
        self.num_params = len(self.param_names)

        self.derivatives_images = self.derivativesImages()
        self.fisher_matrix_images = self.fisherMatrixImages()
        self.fisher_matrix = self.fisherMatrix()
        self.fisher_matrix_chi2 = self.fisherMatrixChi2()
        self.covariance_matrix = self.covarianceMatrix()
        self.correlation_matrix = self.correlationMatrix()
        self.second_derivatives_images = self.secondDerivativesImages()
        self.bias_matrix_images = self.biasMatrixImages()
        self.bias_matrix = self.biasMatrix()
        self.bias_images = self.biasImages()
        self.biases = self.biases()      

    def derivativesImages(self):
        #create a dictionary with the derivatives of the model with respect to each parameter.
        return {self.param_names[i]:fns.partialDifferentiate(func = fns.drawGalaxies, parameter = self.param_names[i], step = self.steps[self.param_names[i]])(self.params) for i in range(self.num_params)}

    def fisherMatrixImages(self):
        FisherM_images = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                FisherM_images[self.param_names[i],self.param_names[j]] = (self.derivatives_images[self.param_names[i]] * self.derivatives_images[self.param_names[j]]) /(self.sigma_n**2)
        return FisherM_images

    def fisherMatrix(self):
        FisherM = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                FisherM[self.param_names[i],self.param_names[j]] = self.fisher_matrix_images[self.param_names[i],self.param_names[j]].sum() #sum over all pixels. 
        return FisherM

    def fisherMatrixChi2(self):
        FisherM_chi2 = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                FisherM_chi2[self.param_names[i],self.param_names[j]] = .5 * fns.secondPartialDifferentiate(fns.chi2, self.param_names[i], self.param_names[j], self.steps[self.param_names[i]], self.steps[self.param_names[j]], sigma_n = self.
                    sigma_n, gal_image = self.gal_image)(self.params)

        return FisherM_chi2

    def covarianceMatrix(self):
        #Covariance matrix is inverse of Fisher Matrix:
        CovM = {}
        FisherM_array = np.array([[self.fisher_matrix[self.param_names[i],self.param_names[j]] for i in range(self.num_params)] for j in range(self.num_params)])#convert to numpy array because it can be useful.
        CovM_array = np.linalg.inv(FisherM_array) #need to be an array to inverse.

        for i in range(self.num_params):
            for j in range(self.num_params):
                CovM[self.param_names[i],self.param_names[j]] = CovM_array[i][j]

        return CovM
    
    def correlationMatrix(self):
        correlation_matrix = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                correlation_matrix[self.param_names[i],self.param_names[j]] = self.covariance_matrix[self.param_names[i],self.param_names[j]] / math.sqrt(self.covariance_matrix[self.param_names[i],self.param_names[i]] * self.covariance_matrix[self.param_names[j],self.param_names[j]])

        return correlation_matrix

    def secondDerivativesImages(self):
        secondDs_gal = {}
        for i in range(self.num_params): 
            for j in range(self.num_params):
                secondDs_gal[self.param_names[i],self.param_names[j]] = (fns.secondPartialDifferentiate(fns.drawGalaxy, self.param_names[i], self.param_names[j], self.steps[self.param_names[i]], self.steps[self.param_names[j]])(self.params))

        return secondDs_gal

    def biasMatrixImages(self):
        BiasM_images = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                for k in range(self.num_params): 
                    BiasM_images[self.param_names[i],self.param_names[j],self.param_names[k]] = (self.derivatives_images[self.param_names[i]] * self.second_derivatives_images[self.param_names[j],self.param_names[k]]) / (self.sigma_n**2)

        return BiasM_images

    def biasMatrix(self):
        BiasM = {}
        for i in range(self.num_params):
            for  j in range(self.num_params):
                for k in range(self.num_params):
                    BiasM[self.param_names[i],self.param_names[j],self.param_names[k]] = self.bias_matrix_images[self.param_names[i],self.param_names[j],self.param_names[k]].sum()

        return BiasM

    def biasImages(self):
        #now we want bias of each parameter per pixel, so we can see how each parameter contributes.
        bias_images = {}
        for i in range(self.num_params):
            sumation = 0 
            for j in range(self.num_params):
                for k in range(self.num_params):
                    for l in range(self.num_params):
                        sumation += self.covariance_matrix[self.param_names[i],self.param_names[j]]*self.covariance_matrix[self.param_names[k],self.param_names[l]]*self.bias_matrix_images[self.param_names[j],self.param_names[k],self.param_names[l]]
            bias_images[self.param_names[i]] = (-.5) * sumation

        return bias_images


    def biases(self):
        return {self.param_names[i]:self.bias_images[self.param_names[i]].sum() for i in range(self.num_params)}







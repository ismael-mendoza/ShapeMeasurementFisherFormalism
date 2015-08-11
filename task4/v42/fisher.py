#!/usr/bin/env python
"""This module contains functions necessary to produce statistical results of the fisher formalism from a given galaxy"""

import math

import numpy as np

import copy 

import galfun

import defaults


def partialDifferentiate(func, param, steps, **kwargs):
    """Partially derive f with respect to a parameter with a certain step.
    We are assuming that the function has a certain structure, namely, 
    one of its arguments is a dictionary of variables that can be changed and
    other (**kwargs) arguments are requisites or extra variables the
    function needs to be evaluated. Assume steps is a dictionary.
    """

    def Dfunc(params):
        """Evaluate the partial derivative at params."""
        params_up = copy.deepcopy(params) #avoids altering params later.
        params_up[param] += steps[param] #increment the value of the parameter by step. 

        params_down = copy.deepcopy(params)
        params_down[param] -= steps[param]

        return ((func(params_up, **kwargs) - func(params_down, **kwargs)) /
               (2 * steps[param]))

    return Dfunc


def secondPartialDifferentiate(func, param1, param2, steps,**kwargs): 
    Df = partialDifferentiate(func, param1, steps, **kwargs)
    return partialDifferentiate(Df, param2, steps)


def chi2(params, gal_image, var_noise, **kwargs): 
    """Returns chi2 given the modified parameters and the original galaxy,
    assume var_noise is the same for all pixels
    """
    pixel_diff = (gal_image - galfun.drawGalaxies(params, **kwargs)).array
    chi = pixel_diff.sum()/math.sqrt(var_noise)
    return chi**2


class Fisher:
    """Given a galaxy image and the appropiate parameters that describe it, will produce a fisher objec that contains the analysis of it using the fisher formalism"""

    def __init__(self, g_parameters, snr): 
        self.g_parameters = g_parameters
        self.snr = snr
        self.image = galfun.drawGalaxies(g_parameters=self.g_parameters, 
                                         image = True)
        _,self.var_noise = galfun.addNoise(self.image, self.snr, 0)
        self.steps = defaults.getSteps(self.g_parameters)
        self.param_names = g_parameters.ordered_fit_names
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
        partials_images = {}
        for i in range(self.num_params):
            partials_images[self.param_names[i]] = partialDifferentiate(
            func = galfun.drawGalaxies, param = self.param_names[i], 
            steps = self.steps)(params = self.g_parameters.params)
        return partials_images


    def fisherMatrixImages(self):
        FisherM_images = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                FisherM_images[self.param_names[i],self.param_names[j]] = (self.derivatives_images[self.param_names[i]] * self.derivatives_images[self.param_names[j]]) /self.var_noise
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
                FisherM_chi2[self.param_names[i],self.param_names[j]] = .5 * secondPartialDifferentiate(func = chi2, param1 = self.param_names[i], param2 = self.param_names[j], steps = self.steps, var_noise = self.var_noise, gal_image = self.image)(params = self.g_parameters.params)

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
                secondDs_gal[self.param_names[i],self.param_names[j]] = secondPartialDifferentiate(func = galfun.drawGalaxies, param1 = self.param_names[i], param2 =self.param_names[j], steps = self.steps)(params = self.g_parameters.params)

        return secondDs_gal

    def biasMatrixImages(self):
        BiasM_images = {}
        for i in range(self.num_params): 
            for j in range(self.num_params): 
                for k in range(self.num_params): 
                    BiasM_images[self.param_names[i],self.param_names[j],self.param_names[k]] = (self.derivatives_images[self.param_names[i]] * self.second_derivatives_images[self.param_names[j],self.param_names[k]]) / (self.var_noise)

        return BiasM_images

    def biasMatrix(self):
        BiasM = {}
        for i in range(self.num_params):
            for  j in range(self.num_params):
                for k in range(self.num_params):
                    BiasM[self.param_names[i], 
                    self.param_names[j],
                    self.param_names[k]
                    ] = self.bias_matrix_images[self.param_names[i],self.param_names[j],self.param_names[k]].sum()

        return BiasM

    def biasImages(self):
        #now we want bias of each parameter per pixel, so we can see how each parameter contributes.
        bias_images = {}
        for i in range(self.num_params):
            sumation = 0 
            for j in range(self.num_params):
                for k in range(self.num_params):
                    for l in range(self.num_params):
                        sumation += (self.covariance_matrix[
                            self.param_names[i],
                            self.param_names[j]
                            ]*self.covariance_matrix[self.param_names[k],self.param_names[l]]*self.bias_matrix_images[self.param_names[j],self.param_names[k],self.param_names[l]])
            bias_images[self.param_names[i]] = (-.5) * sumation

        return bias_images

    def biases(self):
        return {
            self.param_names[i]:self.bias_images[self.param_names[i]].sum()
            for i in range(self.num_params)
            }
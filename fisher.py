#!/usr/bin/env python
"""This module contains functions necessary to produce statistical results of
the fisher formalism from a given galaxy.
"""

import math

import numpy as np

import copy

import galfun

import defaults

import scipy.optimize as scipyopt

def partialDifferentiate(func, param, steps, **kwargs):
    """Partially derive f with respect to param with a certain step.

    We are assuming that the function has a certain structure, namely,
    one of its arguments is a dictionary of variables that can be changed and
    other (**kwargs) arguments are requisites or extra variables the
    function needs to be evaluated. Assume steps is a dictionary.
    """
    def Dfunc(params):
        """Evaluate the partial derivative at params."""
        params_up = copy.deepcopy(params)  # avoids altering params later.
        # increment the value of the parameter by step.
        params_up[param] += steps[param]

        params_down = copy.deepcopy(params)
        params_down[param] -= steps[param]

        return ((func(params_up, **kwargs) - func(params_down, **kwargs)) /
                (2 * steps[param]))

    return Dfunc


def secondPartialDifferentiate(func, param1, param2, steps, **kwargs):
    """Find second partial derivative of the given function"""
    Df = partialDifferentiate(func, param1, steps, **kwargs)
    return partialDifferentiate(Df, param2, steps)


class Fisher(object):
    """Produce fisher object(containing fisher analysis) for a given set of
    galaxy parameters.

    Given a galaxy image and the appropiate parameters that describe it,
    will produce a fisher object that contains the analysis of it using the
    fisher formalism.
    """

    def __init__(self, g_parameters, snr):
        self.g_parameters = g_parameters
        self.snr = snr
        self.image = galfun.drawGalaxies(g_parameters=self.g_parameters,
                                         image=True)
        _, self.var_noise = galfun.addNoise(self.image, self.snr, 0)
        self.steps = defaults.getSteps(self.g_parameters)
        self.param_names = g_parameters.ordered_fit_names
        self.num_params = len(self.param_names)
        self.num_galaxies = self.g_parameters.num_galaxies

        self.derivatives_images = self.derivativesImages()
        self.fisher_matrix_images = self.fisherMatrixImages()
        self.fisher_matrix = self.fisherMatrix()
        self.covariance_matrix = self.covarianceMatrix()
        self.correlation_matrix = self.correlationMatrix()
        self.second_derivatives_images = self.secondDerivativesImages()
        self.bias_matrix_images = self.biasMatrixImages()
        self.bias_matrix = self.biasMatrix()
        self.bias_images = self.biasImages()
        self.biases = self.getBiases()

        self.fisher_condition_number = self.fisherConditionNumber()

    def derivativesImages(self):
        """Return images of the partial derivatives of the galaxy.

        The partial differentiation includes each of the different parameters
        that describe the galaxy.
        """
        partials_images = {}
        for i in range(self.num_params):
            partials_images[self.param_names[i]] = partialDifferentiate(
                func=galfun.drawGalaxies, param=self.param_names[i],
                steps=self.steps)(params=self.g_parameters.params)
        return partials_images

    def fisherMatrixImages(self):
        """Produce images of fisher matrix)."""
        FisherM_images = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]
                derivative1 = self.derivatives_images[param_i]
                derivative2 = self.derivatives_images[param_j]
                FisherM_images[param_i, param_j] = (
                    derivative1 * derivative2 / self.var_noise)
        return FisherM_images

    def fisherMatrix(self):
        """Calculate the actual values of the fisher matrix."""
        FisherM = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]
                FisherM[param_i, param_j] = (
                    self.fisher_matrix_images[param_i, param_j].sum())
        return FisherM

    def matrixToNumpyArray(self, matrix):
        """Convert matrix dictionary to a numpy array"""
        array = np.zeros([self.num_params, self.num_params])
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]
                element = matrix[param_i, param_j]
                array[i][j] = element
        return array

    def numpyArrayToMatrix(self, array):
        """Convert numpy array to matrix dictionary"""
        matrix = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]
                matrix[param_i, param_j] = array[i][j]
        return matrix

    def covarianceMatrix(self):
        """Calculate the covariance matrix by inverting fisher matrix."""
        fisher_array = self.matrixToNumpyArray(self.fisher_matrix)
        covariance_array = np.linalg.inv(fisher_array)
        return self.numpyArrayToMatrix(covariance_array)

    def correlationMatrix(self):
        """Calculate correlation matrix from the covariance matrix."""
        correlation_matrix = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]
                sigma_ij = self.covariance_matrix[param_i, param_j]
                sigma_i = math.sqrt(self.covariance_matrix[param_i, param_i])
                sigma_j = math.sqrt(self.covariance_matrix[param_j, param_j])
                correlation_matrix[param_i, param_j] = (sigma_ij /
                                                       (sigma_i*sigma_j))

        return correlation_matrix

    def secondDerivativesImages(self):
        """Return the images for the second derivatives of the given galaxy."""
        secondDs_gal = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]
                secondDs_gal[param_i, param_j] = secondPartialDifferentiate(
                    func=galfun.drawGalaxies, param1=param_i, param2=param_j,
                    steps=self.steps)(params=self.g_parameters.params)

        return secondDs_gal

    def biasMatrixImages(self):
        """Produce images of each element of the bias matrix"""
        BiasM_images = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                for k in range(self.num_params):
                    param_i = self.param_names[i]
                    param_j = self.param_names[j]
                    param_k = self.param_names[k]
                    BiasM_images[param_i, param_j, param_k] = (
                        self.derivatives_images[param_i] *
                        self.second_derivatives_images[param_j, param_k] /
                        self.var_noise)

        return BiasM_images

    def biasMatrix(self):
        """Return bias matrix from the images of the bias matrix"""
        BiasM = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                for k in range(self.num_params):
                    param_i = self.param_names[i]
                    param_j = self.param_names[j]
                    param_k = self.param_names[k]
                    BiasM[param_i, param_j, param_k] = self.bias_matrix_images[
                        param_i, param_j, param_k].sum()
        return BiasM

    def biasImages(self):
        """Construct the bias of each parameter per pixel."""
        bias_images = {}
        for i in range(self.num_params):
            sumation = 0
            for j in range(self.num_params):
                for k in range(self.num_params):
                    for l in range(self.num_params):
                        param_i = self.param_names[i]
                        param_j = self.param_names[j]
                        param_k = self.param_names[k]
                        param_l = self.param_names[l]
                        sumation += (self.covariance_matrix[param_i, param_j] *
                                     self.covariance_matrix[param_k, param_l] *
                                     self.bias_matrix_images[param_j, param_k,
                                                             param_l])
            bias_images[self.param_names[i]] = (-.5) * sumation
        return bias_images

    def getBiases(self):
        """Return the value of the bias of each parameter in vector form."""
        return {
            self.param_names[i]: self.bias_images[self.param_names[i]].sum()
            for i in range(self.num_params)
        }

    def fisherConditionNumber(self):
        fisher_array = self.matrixToNumpyArray(self.fisher_matrix)
        return np.linalg.cond(fisher_array)

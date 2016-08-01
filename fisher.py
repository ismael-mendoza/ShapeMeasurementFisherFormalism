"""This module contains functions necessary to produce statistical results of
the fisher formalism from a given galaxy.
"""

import math
import numpy as np
import copy
import galfun
import defaults

class Fisher(object):
    """Produce fisher object (containing fisher analysis) for a given set of
    galaxy parameters.

    Given a galaxy image and the appropiate parameters that describe it,
    will produce a fisher object that contains the analysis of it using the
    Fisher Formalism.
    """

    def __init__(self, g_parameters, image_renderer, snr):
        self.g_parameters = g_parameters
        self.snr = snr
        self.model = galfun.getGalaxiesModels(g_parameters=self.g_parameters)
        self.image_renderer = image_renderer

        #we do not want to mask or crop the images used to obtain the partials. 
        self.image_renderer_partials = galfun.ImageRenderer(stamp=self.image_renderer.stamp)
        self.image = self.image_renderer.getImage(self.model)
        _, self.var_noise = galfun.addNoise(self.image, self.snr, 0)

        self.steps = defaults.getSteps(self.g_parameters, self.image_renderer)
        self.param_names = g_parameters.ordered_fit_names
        self.num_params = len(self.param_names)
        self.num_galaxies = self.g_parameters.num_galaxies

        self.derivatives_images = self.derivativesImages()
        self.second_derivatives_images = self.secondDerivativesImages()
        self.fisher_matrix_images = self.fisherMatrixImages()
        self.fisher_matrix = self.fisherMatrix()
        self.covariance_matrix = self.covarianceMatrix()
        self.correlation_matrix = self.correlationMatrix()
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
            param = self.param_names[i]
            params_up = copy.deepcopy(self.g_parameters.params)
            params_up[param] += self.steps[param]
            params_down = copy.deepcopy(self.g_parameters.params)
            params_down[param] -= self.steps[param]
            gal_up = galfun.getGalaxiesModels(params_up)
            gal_down = galfun.getGalaxiesModels(params_down)
            img_up = self.image_renderer_partials.getImage(gal_up)
            img_down = self.image_renderer_partials.getImage(gal_down)
            partials_images[param] =  ((img_up - img_down)/(2 * self.steps[param])).array
        return partials_images


    def secondDerivativesImages(self):
        """Return the images for the second derivatives of the given galaxy."""
        secondDs_gal = {}
        for i in range(self.num_params):
            for j in range(self.num_params):
                param_i = self.param_names[i]
                param_j = self.param_names[j]

                params_iup_jup = copy.deepcopy(self.g_parameters.params)
                params_iup_jup[param_i] += self.steps[param_i]
                params_iup_jup[param_j] += self.steps[param_j]

                params_idown_jup = copy.deepcopy(self.g_parameters.params)
                params_idown_jup[param_i] -= self.steps[param_i]
                params_idown_jup[param_j] += self.steps[param_j]   

                params_iup_jdown = copy.deepcopy(self.g_parameters.params)
                params_iup_jdown[param_i] += self.steps[param_i]
                params_iup_jdown[param_j] -= self.steps[param_j] 


                params_idown_jdown = copy.deepcopy(self.g_parameters.params)
                params_idown_jdown[param_i] -= self.steps[param_i]
                params_idown_jdown[param_j] -= self.steps[param_j]

                gal_iup_jup = galfun.getGalaxiesModels(params_iup_jup)
                gal_idown_jup = galfun.getGalaxiesModels(params_idown_jup)
                gal_iup_jdown = galfun.getGalaxiesModels(params_iup_jdown)
                gal_idown_jdown = galfun.getGalaxiesModels(params_idown_jdown)

                img_iup_jup = self.image_renderer_partials.getImage(gal_iup_jup)
                img_idown_jup = self.image_renderer_partials.getImage(gal_idown_jup)
                img_iup_jdown = self.image_renderer_partials.getImage(gal_iup_jdown)
                img_idown_jdown = self.image_renderer_partials.getImage(gal_idown_jdown)

                secondDs_gal[param_i, param_j] = ((img_iup_jup + img_idown_jdown - 
                                                   img_idown_jup - img_iup_jdown)/
                                                   (4*self.steps[param_i]*self.steps[param_j])).array

        return secondDs_gal

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
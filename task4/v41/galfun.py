#!/usr/bin/env python

import sys 

import os

import csv

import defaults

import galsim

import plotsfisher

import matplotlib.pyplot as plt

import copy

import names

sys.path.insert(0,'tools')

import prof


def drawGalaxy(params): 
    """Draws an image of a galaxy with noise and convoluted with a given psf if desired. Also sets 
    GaussianNoise to the image if SNR is given.
    """

    #check each model of galaxy
    if(params['model'] == 'gaussian'):
        if('hlr' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], half_light_radius=params['hlr'])
        elif('sigma' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

    elif(params['model'] == 'exponential'):
        gal = galsim.Exponential(flux = params['flux'], half_light_radius = params['hlr'])

    elif(params['model'] == 'bulge+disk'):
        pass


    #set shear
    if('e1' and 'e2' in params.keys()): 
        gal = gal.shear(e1=params['e1'], e2 = params['e2'])
    elif('eta1' and 'eta2' in params.keys()):
        gal = gal.shear(eta1=params['eta1'], eta2=params['eta2'])
    elif('q' and 'beta' in params.keys()):
        gal = gal.shear(q=params['q'], beta = params['beta'] * galsim.radians)

    #shift galaxy.
    gal = gal.shift(params['x0'],params['y0'])

    #generate psf, if psf_flux is not 0 and it was specified by the user.
    if(params.get('psf_flux', 0) != 0):
        #each psf model. 
        if(params['psf_model'] == 'gaussian'):
            if('psf_hlr' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], half_light_radius=params['psf_hlr'])
            elif('psf_sigma' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], 
                                      sigma=params['psf_sigma'])
            elif('psf_fwhm' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], 
                                      fwhm=params['psf_fwhm'])

        psf = psf.shear(e1=params['psf_e1'],e2=params['psf_e2'])

        final = galsim.Convolve([gal, psf])
    else:
        final = gal

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    return final.drawImage(scale=defaults.PIXEL_SCALE, nx=defaults.NX, 
                            ny=defaults.NY)


def drawGalaxies(fit_params=None, id_params=None, image=False, 
                 g_parameters=None, **kwargs):
    """Draw each of the galaxy represented in gals_params, important to note: 
    *Convolution can be done at the end since we will assume that blended galaxies are in the same field of view/stamp 
    so that they are convoluted with the same psf and convolution is an additive operation (integrating is additive)
    *Convolution can also be done for each galaxy separately and then adding them, we will always assume galaxies have same psf. 
    *Noise should only be added once to the image (otherwise two different noises combined), so it should go here and not in the original drawGalaxy function. 
    """

    if id_params is None and g_parameters is None:
        fit_params.update(kwargs)
        id_params = GParameters.convertParams_Id(fit_params)

    if g_parameters is not None:
        id_params = g_parameters.id_params

    gals = []
    
    for gal_id in id_params:
        gals.append(drawGalaxy(id_params[gal_id]))
    final = sum(gals)

    if image is False:
        return final.array
    else:
        return final


def displayGalaxy(image):
    """Displays a galaxy given an image of a galaxy, mostly used just for 
    testing.
    """
    figure, subplt= plt.subplots(1,1)
    figure.suptitle('Galaxy', fontsize = defaults.FONTSIZE_TITLE)
    plotsfisher.drawImage(subplt, image.array)
    plt.show()


def addNoise(image, snr, noise_seed):
    #set gaussian noise with given noise seed.
    noisy_image = copy.deepcopy(image) # do not alter original image.
    bd = galsim.BaseDeviate(noise_seed)
    noise = galsim.GaussianNoise(rng=bd)
    variance_noise = noisy_image.addNoiseSNR(noise, snr,  
                                            preserve_flux =True) 
    return noisy_image, variance_noise


class GParameters:
    """Class that manages given galaxy parameters."""
    def __init__(self, project): 
        #some error handling.
        if not os.path.isdir(project):
            raise OSError('Directory given does not exist.')

        filename = os.path.join(project, defaults.GALAXY_FILE)
        if not os.path.isfile(filename):
            raise OSError('The given file name is not in the directory:')

        #extract params from each of the rows (galaxies) in the given csvfile.
        #also remove empty params. 
        with open(filename, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            id_params = {}
            for row in reader:
                row_to_store = copy.deepcopy(row) 
                gal_id = row.pop('id')
                for param in row:
                    if not row[param]:
                        row_to_store.pop(param)
                id_params[gal_id] = row_to_store

        #convert all appropiate values to floats,
        for gal_id in id_params: 
            for key, value in id_params[gal_id].iteritems():
                try:
                    id_params[gal_id][key] = float(value)
                except ValueError:
                    pass

        self.omit_fit = names.omit_fit
        self.id_params = id_params
        self.params = GParameters.convertId_Params(self.id_params)
        self.fit_params = GParameters.convertId_Params(self.id_params,
                                                       self.omit_fit)
        self.ordered_fit_names = self.sortModelParamsNames() 
        self.num_galaxies = len(self.id_params.keys())
        print self.id_params, self.params, self.fit_params, self.ordered_fit_names,self.num_galaxies 

    # def convertModel_Fit(self):
    #     fit_id_params = GParameters.convertParams_Id(self.model_params)
    #     for gal_id in self.id_params.keys():
    #          e1 = self.id_params[gal_id]['e1']
    #          e2 = self.id_params[gal_id]['e2']
    #          hlr = self.id_params[gal_id]['hlr']
    #          profile = prof.Profile(e1=e1,e2=e2,hlr=hlr)
    #          fit_id_params[gal_id].pop('e1')
    #          fit_id_params[gal_id].pop('e2')
    #          fit_id_params[gal_id]['eta1'] = profile.eta1
    #          fit_id_params[gal_id]['eta2'] = profile.eta2

    #     return GParameters.convertId_Params(fit_id_params)

    # @staticmethod
    # def convertFit_Model(fit_params):
    #     fit_id_params = GParameters.convertParams_Id(fit_params)
    #     for gal_id in fit_id_params:
    #          eta1 = fit_id_params[gal_id]['eta1']
    #          eta2 = fit_id_params[gal_id]['eta2']
    #          hlr = fit_id_params[gal_id]['hlr']
    #          profile = prof.Profile(eta1=eta1,eta2=eta2,hlr=hlr)
    #          fit_id_params[gal_id].pop('eta1')
    #          fit_id_params[gal_id].pop('eta2')
    #          fit_id_params[gal_id]['e1'] = profile.e1
    #          fit_id_params[gal_id]['e2'] = profile.e2

    #     return GParameters.convertId_Params(fit_id_params)


    def sortModelParamsNames(self):
        """Returns param names in the order _1, _2,..., and in each 
        subscript follow 
        order of defaults.py, retuns list of ordered params names. With it we 
        can change order just by changing order in names file gal_parameters.
        Avoids omit. 
        """
        param_names = []
        for gal_id in self.id_params:
            galaxy_model = self.id_params[gal_id]['galaxy_model']
            for name in names.gal_models_parameters[galaxy_model]:
                for param in self.id_params[gal_id]:
                    if param not in self.omit_fit:
                        if param == name:
                            param_names.append(param + '_'+ str(gal_id))
        return param_names

    @staticmethod
    def convertId_Params(id_params, omit_fit = {}): 
        """Uses a dictionary of the form id:params, where params is a dictionary of parameters and converts it to a dictionary of the form param_id:value(param)
        *omit (if given) is a list of strings that are not desired in the final dictionary that are in the original idDict,
        useful because to remove items later have to access it by param_1 for example.
        """

        params = {}
        for gal_id in id_params:
            galaxy_model = id_params[gal_id]['galaxy_model'] 
            for param in id_params[gal_id]:
                if param != 'id':
                    if param not in omit_fit.get(galaxy_model,[]):
                        params[param + '_'+ str(gal_id)] = (
                        id_params[gal_id][param])
        return params

    @staticmethod
    def convertParams_Id(params):
        id_params ={}
        ids = []
        for param in params.keys():
            if param[-1] not in ids:
                ids.append(param[-1]) #appends last character of param

        for gal_id in ids:
            ID_params = {}
            for param in params.keys():
                if param[-1] == gal_id:
                    #slice last 2 characters to avoid '_1'
                    ID_params[param[:-2]] = params[param]
            id_params[gal_id] = ID_params

        return id_params
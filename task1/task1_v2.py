import sys
import os
import math
import logging
import galsim
from lmfit import minimize, Parameters, Parameter, fit_report
import numpy as np

def drawImage(Params, logger, filename = "", AddNoise = False, ProduceFile = False): 
    """Draws image with noise when are simulating the galaxy and without noise when we are creating the model. 
    This is a good way to recycle code."""

    params = Params.valuesdict()

    # Define the galaxy profile
    gal = galsim.Gaussian(flux=params['gal_flux'], sigma=params['gal_sigma'])
    logger.debug('Made galaxy profile')

    #set shear, if desired. 
    gal.shear = galsim.Shear(e1 = params['e1'], e2 = params['e2'])
    logger.debug('set shear for the galaxy')

    #shift galaxy if desired. be careful to do shear and shift in this order, because otherwise might act weirdly (takes center as the original for shear)
    gal = gal.shift(params['x0'],params['y0']) 
    
    # Define the PSF profile
    psf = galsim.Gaussian(flux=params['psf_flux'], sigma=params['psf_sigma']) # PSF flux should always = 1
    logger.debug('Made PSF profile')

    # Final profile is the convolution of these
    # Can include any number of things in the list, all of which are convolved 
    # together to make the final flux profile.
    final = galsim.Convolve([gal, psf])
    logger.debug('Convolved components into final profile')

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    # The returned image has a member, added_flux, which is gives the total flux actually added to 
    # the image.  One could use this value to check if the image is large enough for some desired
    # accuracy level.  Here, we just ignore it.
    image = final.drawImage(scale=params['pixel_scale'], nx = params['nx'], ny = params['ny'])
    logger.debug('Made image of the profile: flux = %f, added_flux = %f',params['gal_flux'],image.added_flux)

    if(AddNoise == True): 

        #dict with variance and images 
        temp_dct = dict()

        #Add Gaussian noise to the image with specified sigma
        variance = image.addNoiseSNR(galsim.GaussianNoise(), params['SNR_gal']) 

        temp_dct['target_image'] = image
        temp_dct['variance_noise'] = variance

        if(ProduceFile == True and filename != ""): 
            outputImageFile(image, filename, logger)

        #have to return variance of noise added as well as image.
        return temp_dct

    if(ProduceFile == True and filename != ""): 
        outputImageFile(image, filename, logger)

    return image

def outputImageFile(image, filename, logger): 
    """Write the image to a file in the output directory (creates this directory if needed"""
    if not os.path.isdir('output'):
        os.mkdir('output')
    file_name = os.path.join('output', filename + '.fits')
    # Note: if the file already exists, this will overwrite it.
    image.write(file_name)
    logger.info('Wrote image to %r' % file_name)  # using %r adds quotes around filename for us


def getFilename(): 
    """Returns file name without any extensions"""
    return str(os.path.basename(__file__))[0:str(os.path.basename(__file__)).find('.py')] 


def mergeParamters(params1, params2): 
    """Merge two parameter dictionaries, only adds the missing entries in the first dictionary"""
    for param in params2:
        if(param not in params1): 
            params1.add( params2[param].name, value = params2[param].value, min = params2[param].min, max = params2[param].max, vary = params2[param].vary, expr = params2[param].expr)

    return params1

#we know want to create the objective function. 
def objfunc(fit_params, data, init_parameters, logger, variance, chi_fit):

    #merge fit_params and init_parameters
    modified_init_parameters = mergeParamters(fit_params, init_parameters)

    #important to set bounds otherwise cannot compare the 2 images
    model = drawImage(modified_init_parameters, logger).array.ravel()

    if(chi_fit == False):
        return model - data

    else: 
        return (model - data) / variance

def main(argv):
    """
    First draw a galaxy with the following characteristics, 
      - Use a circular Gaussian profile for the galaxy.
      - Convolve it by a circular Gaussian PSF.
      - Add Gaussian noise to the image.
      - It is sheared so that e1 != 0 and e2!= 0  and shifted. 

    Then by using lmfit we obtain the parameters solely from the resulting image of the galaxy with the method of least squares. 
    """

    #this way you can get the name of the file you are currently in without any extensions.
    filename = getFilename()


    # In non-script code, use getLogger(__name__) at module scope instead.
    #logging levels describe what messages of logging to output in what situation, this is convenient look at documentation. 
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger(filename)    

    #constants (never going to use in a fit): (do not put in dictionary)
    pixel_scale,= .4


    
    #first we want to set here the true values. Or the ones we want our galaxy to have.    
    #possible parameters for the galaxy formation. 
    init_parameters = Parameters()
    init_parameters.add('gal_sigma', value = 3., min = 0) # arcsec 
    init_parameters.add('gal_flux', value = 100., min = 0)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    init_parameters.add('e1', value = 0.) #ellipticity: e1 
    init_parameters.add('e2', value = 0.)#ellipticity: e2
    init_parameters.add('x0', value = 0.) #shift in x origin. 
    init_parameters.add('y0', value = 0.)     #shift in y

    #list of parameters to fit, defined above for order and to only pass them to fit_params.
    fit_params = init_parameters

    #other initial values given
    init_parameters.add('SNR_gal', value = 20.) # should probably be 20 - 40  since it is what we are interested in.
    init_parameters.add('psf_flux', value = 1.)
    init_parameters.add('psf_sigma', value = 1., min =0)
    init_parameters.add('nx', value = 68)
    init_parameters.add('ny', value = 68)

    #get image and variance from list returned
    temp_dct = drawImage(init_parameters, logger, filename, AddNoise = True, ProduceFile = True)
    target_image = temp_dct['target_image']
    variance_noise = temp_dct['variance_noise']


    results = target_image.FindAdaptiveMom() # (Moments) algorithm that finds moments for the profile displayed.

    logger.info('HSM reports that the image has observed shape and size:')
    logger.info('    e1 = %.3f, e2 = %.3f, sigma = %.3f (pixels)', results.observed_shape.e1,
                results.observed_shape.e2, results.moments_sigma)
    logger.info('Expected values in the limit that pixel response and noise are negligible:')
    logger.info('    e1 = %.3f, e2 = %.3f, sigma = %.3f', 0.0, 0.0, 
                math.sqrt(init_parameters['gal_sigma']**2 + init_parameters['psf_sigma']**2)/init_parameters['pixel_scale']) 
    

    ######## now we do the fitting in the image. 

    #we initialize to random values and set bounds accordingly.
    fit_params.add('gal_sigma', value = 1., min = 0) #the min here is important because the program will crash if you have a negative size.
    fit_params.add('gal_flux', value = 10., min = 0)
    fit_params.add('e1', value = 0, min = -.7, max = .7) #can be both positive or ngative 
    fit_params.add('e2', value = 0, min = -.7, max = .7) #best way to do it satisfactorily, does not necessarily work always but it is enough for the galaxies we are looking at now.
    fit_params.add('x0', value = 0, min = target_image.getXMin(), max= target_image.getXMax())
    fit_params.add('y0', value = 0, min = target_image.getYMin(), max =target_image.getYMax())

    #this is the function that obtains the fit. 
    chi_fit = False #if want to minimize chi (take into account the noise bias)
    minimize(objfunc, fit_params, args=(target_image.array.ravel(),init_parameters, logger, variance_noise, chi_fit))


    print("")
    print("Start fit_rerpot report:")
    print(fit_report(fit_params)) #it is very interesting how it obtains the original sigma from the galaxy and not the convolution, need to look at this in more detail. 


if __name__ == "__main__":
    main(sys.argv)

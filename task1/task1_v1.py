import sys
import os
import math
import logging
import galsim
from lmfit import minimize, Parameters, Parameter, fit_report
import numpy as np

def drawImage(params, filename, AddNoise = False, ProduceFile = False, ): 
    """Draws image with noise when are simulating the galaxy and without noise when we are creating the model. 
    This is a good way to recycle code."""

    # Define the galaxy profile
    gal = galsim.Gaussian(flux=params['gal_flux'], sigma=params['gal_sigma'])
    logger.debug('Made galaxy profile')

    #set shear, 
    gal.shear = galsim.Shear(e1 = params['e1'], e2 = params['e2'])
    logger.debug('set shear for the galaxy')
    
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
    image = final.drawImage(scale=pixel_scale)
    logger.debug('Made image of the profile: flux = %f, added_flux = %f',gal_flux,target_image.added_flux)


    if(AddNoise == True): 
        #Add Gaussian noise to the image with specified sigma
        image.addNoise(galsim.GaussianNoise(sigma=noise))
        logger.debug('Added Gaussian noise')


    if(ProduceFile == True): 
        outputImageFile(filename)

    return image







def outputImageFile(filename): 
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

    #possible parameters for the galaxy formation. 
    init_parameters = dict()

    #first we want to set here the true values. Or the ones we want our galaxy to have.
    #also encodes psf and other but we can get rid of those later, put unwanted parameters at the end
    
    #parameters to fit
    init_parameters['gal_sigma'] = 100. # arcsec 
    init_parameters['gal_flux'] = 3.   # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    init_parameters['e1'] = 0. #ellipticity: e1 
    init_parameters['e2'] = 0. #ellipticity: e2
    init_parameters['x0'] = 0. #shift in x origin. 
    init_parameters['y0'] = 0.     #shift in y

    #list of parameters to fit, have to be defined above this function call
    parameters_to_fit = init_parameters.keys() 

    #other initial values given 
    init_parameters['pixel_scale'] = .4 # arcsec / pixel 
    init_parameters['SNR_gal'] = 20. # should probably be 20 - 40  since it is what we are interested in.
    init_parameters['psf_flux'] = 1.
    init_parameters['psf_simga'] = 1.


    target_image = drawImage(init_parameters, filename, AddNoise = True, ProduceFile = False)
    
    results = target_image.FindAdaptiveMom() # (Moments) algorithm that finds moments for the profile displayed.

    logger.info('HSM reports that the image has observed shape and size:')
    logger.info('    e1 = %.3f, e2 = %.3f, sigma = %.3f (pixels)', results.observed_shape.e1,
                results.observed_shape.e2, results.moments_sigma)
    logger.info('Expected values in the limit that pixel response and noise are negligible:')
    logger.info('    e1 = %.3f, e2 = %.3f, sigma = %.3f', 0.0, 0.0, 
                math.sqrt(gal_sigma**2 + psf_sigma**2)/pixel_scale) 



    ###################fit to obtain parameters, afterwards the only thing I can work with is the image. 
    #want to shear it first, then shift it because otherwise it will shear relative to original origin 
    #and so it will shear it incorrectly. 

    #apparently image(x,y) represents the flux of that given pixel of the image.

    #x = np.linspace(1,21,21) #can also use range as we do not need floating point numbers for pixels, 
    # y = np.linspace(1,68,68)


    #do not need x,y or whatever variables to create a range because the image is an array itself

    #we know want to create the objective function. 
    def objfunc(params,data,image_bounds):
        #the model is a gaussian.  

        #obtain a dictionary with the values 

        values = params.valuesdict()

        gal_sigma = values['gal_sigma']
        gal_flux = values['gal_flux']

        #no way to extract the pixel_scale from the image. 
        new_gal = galsim.Gaussian(flux=gal_flux, sigma=gal_sigma)

        #shear the galaxy

        #shift the galayx. 

        #convolve with psf which we assume it is given for now. 

        image_model = 

        #important to set bounds otherwise cannot compare the 2 images
        model = new_gal.drawImage(scale=pixel_scale, bounds=image_bounds).array.ravel()

        return model - data
    #create set of params 

    params.add('x0', value = 0, min = target_image.getXMin(), max= target_image.getXMax())
    params.add('y0', value = 0, min = target_image.getYMin(), max =target_image.getYMax())


    #do not forget to set min and max of x0 and y0 now that we have the image


    #this is the function that obtains the fit. 
    minimize(objfunc, params, args=(target_image.array.ravel(),target_image.getBounds()))


    print("")
    print("Start fit_rerpot report:")
    print(fit_report(params)) #it is very interesting how it obtains the original sigma from the galaxy and not the convolution, need to look at this in more detail. 


if __name__ == "__main__":
    main(sys.argv)

"""This program fits 6 parameters for a galaxy using minimizer from lmfit,
the image for the galaxy is generated using galsim. Also reports correlations in the fit_reprt"""


import sys
import os
import math
import logging
import galsim
from lmfit import minimize, Parameters, Parameter, fit_report
import numpy as np



###need to find a way to not clutter functions with these global variables
###but not put them as global variables!!

#global constants for image (never going to use in a fit), they are used and will never change but do not want to put them in params: 
pixel_scale = .4
nx = 100
ny = 100

#this way you can get the name of the file you are currently in without any extensions.
filename = str(os.path.basename(__file__))[0:str(os.path.basename(__file__)).find('.py')]

# In non-script code, use getLogger(__name__) at module scope instead.
#logging levels describe what messages of logging to output in what situation, this is convenient look at documentation. 
logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
logger = logging.getLogger(filename)   

def drawGalaxy(Params): 
    """Draws image with noise when are simulating the galaxy and without noise when we are creating the model. 
    This is a good way to recycle code."""

    params = Params.valuesdict()

    # Define the galaxy profile
    gal = galsim.Gaussian(flux=params['gal_flux'], sigma=params['gal_sigma'])
    logger.debug('Made galaxy profile')

    #set shear, if desired. note that shear generates a new sheared galaxy, should not try to change attribute of the galaxy directly
    gal = gal.shear(e1=params['e1'], e2 = params['e2'])
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
    image = final.drawImage(scale=pixel_scale, nx = nx, ny = ny)
    logger.debug('Made image of the profile: flux = %f, added_flux = %f',params['gal_flux'],image.added_flux)

    return image

def outputImageFile(image): 
    """Write the image to a file in the output directory (creates this directory if needed)"""
    if not os.path.isdir('output'):
        os.mkdir('output')
    file_name = os.path.join('output', filename + '.fits')
    # Note: if the file already exists, this will overwrite it.
    image.write(file_name)
    logger.info('Wrote image to %r' % file_name)  # using %r adds quotes around filename for us


def getFilename(): 
    """Returns file name without any extensions"""
    return str(os.path.basename(__file__))[0:str(os.path.basename(__file__)).find('.py')] 


#we know want to create the objective function. 
def objfuncDif(params, data):
    model = drawGalaxy(params).array.ravel()
    return model - data

def objfuncChi(params, data, variance): 
    model = drawGalaxy(params).array.ravel()
    return (model - data) / math.sqrt(variance)

def main(argv):
    """
    First draw a galaxy with the following characteristics, 
      - Use a circular Gaussian profile for the galaxy.
      - Convolve it by a circular Gaussian PSF.
      - Add Gaussian noise to the image.
      - It is sheared so that e1 != 0 and e2!= 0  and shifted. 

    Then by using lmfit we obtain the parameters solely from the resulting image of the galaxy with the method of least squares. 
    """


    #a good idea is to create a module that has all this parameters stored in them.


    #general configuration of the program. 

    AddNoise = True      #want to add noise to a galaxy. 
    ProduceFile = True  #want to output image file, may take longer. 
    chi_fit = False      #if want to minimize chi (take into account the noise bias)

    #first we want to set here the true values. Or the ones we want our galaxy to have.    
    #possible parameters for the galaxy formation. 
    #in the fit


    #more robust if params = gal_params
    params = Parameters()
    params.add('gal_sigma', value = 3.) # arcsec 
    params.add('gal_flux', value = 100.)  # total counts on the image, watch out if this is too large, can cause problems because FT transform on narrow functions. 
    params.add('e1', value = .5) #ellipticity: e1 
    params.add('e2', value = .5)#ellipticity: e2
    params.add('x0', value = 2.) #shift in x origin. 
    params.add('y0', value = 3.)     #shift in y

    #not in the fit
    params.add('SNR_gal', value = 30., vary = False) # should probably be 20 - 40  since it is what we are interested in.
    params.add('psf_flux', value = 1., vary = False)
    params.add('psf_sigma', value = 1., vary = False)

    #get image and variance from list returned
    target_image = drawGalaxy(params)

    #Add Gaussian noise to the image with specified sigma
    if(AddNoise == True):
        variance_noise = target_image.addNoiseSNR(galsim.GaussianNoise(), params['SNR_gal'], True) 

    if(ProduceFile == True and filename != ""): 
        outputImageFile(target_image)



    ######## now we do the fitting in the image. 

    #we initialize to random values and set bounds accordingly if necessary. only for the parameter we want to fit.  
    params.add('gal_sigma', value = 1., min = 0) #the min here is important because the program will crash if you have a negative size.
    params.add('gal_flux', value = 50., min = 0)
    params.add('e1', value = 0.1, min = -.7, max = .7) #can be both positive or ngative 
    params.add('e2', value = 0.1, min = -.7, max = .7) #best way to do it satisfactorily, does not necessarily work always but it is enough for the galaxies we are looking at now.
    params.add('x0', value = 0, min = target_image.getXMin(), max= target_image.getXMax())
    params.add('y0', value = 0, min = target_image.getYMin(), max =target_image.getYMax())

    #this is the function that obtains the fit. 
    if(chi_fit == True):
        minimize(objfuncDif, params, args=(target_image.array.ravel()))

    else:
        minimize(objfuncChi, params, args=(target_image.array.ravel(), variance_noise))


    print("")
    print("Start fit_rerpot report:")
    print(fit_report(params)) #it is very interesting how it obtains the original sigma from the galaxy and not the convolution, need to look at this in more detail. 


if __name__ == "__main__":
    main(sys.argv)

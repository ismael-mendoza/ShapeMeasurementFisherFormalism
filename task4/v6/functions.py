import os
import defaults
import matplotlib.cm as cm
from copy import deepcopy #useful to not change both variables when one is copied.
import galsim

cts = defaults.constants() #some constants used in the functions.

def csvIsEmpty(filename):
    """checks each row and if any is not empty, then the file is not empty"""
    with open(filename, 'r') as f:
        for row in f:
            if len(row) != 0:
                return False
        else: return True

def SaveFigureToPdf(figure,file_name, dir_name, plotdir_name, hide = True): 

     #save and preview pdf.
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)
    if not os.path.isdir(os.path.join(dir_name,plotdir_name)):
        os.mkdir(os.path.join(dir_name,plotdir_name))
    file_name = os.path.join(dir_name, plotdir_name, file_name)  #puts slashes in between things.
    figure.savefig(file_name, bbox_inches='tight')

    if not hide: 
        os.system("open " + file_name)

def drawPlot(figure, plot, title = "", xlabel = "", ylabel = ""): 
    """draws a given plot with default values using imshow in the given subplot
    if they are not subplots should just pass plt. It also adds a title with name. -- OWN"""

    figure.axes.get_xaxis().set_ticks([]) #remove numbers and ticks but can use label.
    figure.axes.get_yaxis().set_ticks([])

    figure.set_title(title)
    figure.set_xlabel(xlabel)
    figure.set_ylabel(ylabel)

    #be careful none must be in quotes. 
    #vmax and vmin useful to center around zero which looks nice.
    figure.imshow(plot, interpolation = 'None', cmap=cm.RdYlGn,
            origin='lower', extent=[-cts.nx,cts.nx,-cts.ny,cts.ny],
            vmax=abs(plot).max(), vmin=-abs(plot).max())

def partialDifferentiate(func, parameter, step, **kwargs):
    """Partially derive f with respect to a parameter with a certain step.
    We are assuming that the function has a certain structure, namely, one of its arguments is a dictionary of 
    variables that can be changed and other (**kwargs) arguments are requisites or extra variables
    the function needs to be evaluated. This is because we are assuming we can add step to params[parameter]."""
    def Dfunc(params):
        """Evaluate the partial derivative at params."""
        params_up = deepcopy(params) #avoids altering params later.
        params_up[parameter] += step #increment the value of the parameter by step. 

        params_down = deepcopy(params)
        params_down[parameter] -= step 

        return (func(params_up, **kwargs) - func(params_down, **kwargs)) / (2 * step)
    return Dfunc

def secondPartialDifferentiate(func, parameter1, parameter2, step1, step2, **kwargs): 
    Df = partialDifferentiate(func, parameter1, step1, **kwargs)
    return partialDifferentiate(Df, parameter2, step2)

def chi2(params, gal_image, sigma_n, **kwargs): 
    """Returns chi2 given the modified parameters and the original galaxy, assume sigma_n is the same for all pixels -- OWN"""
    return ((((gal_image- drawGalaxy(params, **kwargs)).array/ (sigma_n)))**2).sum()

def drawGalaxy(params, SNR = -1, noise_seed = 0, **kwargs): 
    """Draws an image of a galaxy with noise and convoluted with a given psf if desired."""

    params.update(kwargs)

    #have to do this for each type of galaxy
    if(params['model'] == 'gaussian'):
        if('hlr' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], half_light_radius=params['hlr'])
        elif('sigma' in params.keys()):
            gal = galsim.Gaussian(flux=params['flux'], sigma=params['sigma'])

        #set shear, if desired. note that shear generates a new sheared galaxy, should not try to change attribute of the galaxy directly
        if('e1' and 'e2' in params.keys()): 
            gal = gal.shear(e1=params['e1'], e2 = params['e2'])
        elif('q' and 'beta' in params.keys()):
            gal = gal.shear(q=params['q'], beta = params['beta'] * galsim.radians) ##galsim.radians is only useful when you draw it (sheart it) do not need it anywhereelse.

        #shift galaxy if desired. be careful to do shear and shift in this order, because otherwise might act weirdly (takes center as the original for shear)
        gal = gal.shift(params['x0'],params['y0'])

    #add a gaussian psf if its flux it not 0
    if(params['psf_flux'] != 0):
        #have to do this for each type of psf. 
        if(params['psf_model'] == 'gaussian'):
            if('psf_hlr' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], half_light_radius=params['psf_hlr'])
            elif('psf_sigma' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], sigma=params['psf_sigma'])
            elif('psf_fwhm' in params.keys()):
                psf = galsim.Gaussian(flux=params['psf_flux'], fwhm=params['psf_fwhm'])
            final = galsim.Convolve([gal, psf])
    else:
        final = gal

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    image = final.drawImage(scale=cts.pixel_scale, nx = cts.nx, ny = cts.ny)

    #set noise if SNR does not equal the default (probably only true on fits.py)
    if(SNR != -1):
        #set noise seed 
        bd = galsim.BaseDeviate(noise_seed)
        noise = galsim.GaussianNoise(rng=bd)
        variance_noise = image.addNoiseSNR(noise, SNR, preserve_flux =True)
        return  image, variance_noise

    return image
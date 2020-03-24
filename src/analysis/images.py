from copy import deepcopy
import galsim


class ImageRenderer(object):
    """Object used to produce the image of a galaxy.

    Args:
        pixel_scale(float): Pixel_scale to use in the image (ratio of pixels to arcsecs)
        nx(float): Width of the image in pixels.
        ny(float): height of the image in pixels.
        stamp(galsim.Image): optional, galsim.Image of appropiate dimensions to draw the image. does
            not actually use whatever was originally in the stamp.
        bounds(tuple): When drawn, the image will be clipped to these bounds.
        mask(:class:`np.array`): the pixels selected in this mask will be set to 0.

    One of the the following must be specified:
        * stamp
        * nx,ny,pixel_scale

    This object is made so it can be passsed in to a class :class:`analysis.fisher.Fisher` object.
    """

    def __init__(self, pixel_scale=None, nx=None, ny=None, stamp=None, project=None,
                 bounds=None, mask=None):

        self.pixel_scale = pixel_scale
        self.nx = nx
        self.ny = ny
        self.bounds = bounds
        self.mask = mask
        self.stamp = stamp

        if self.stamp is None:
            if self.nx is not None and self.ny is not None and self.pixel_scale is not None:
                self.stamp = galsim.Image(nx=self.nx, ny=self.ny, scale=self.pixel_scale)

        else:
            self.nx = self.stamp.array.shape[0]
            self.ny = self.stamp.array.shape[1]
            self.pixel_scale = self.stamp.scale

        if self.bounds is not None:
            self.stamp = self.stamp[bounds]

    def get_image(self, galaxy):
        img = deepcopy(self.stamp)
        galaxy.drawImage(image=img, use_true_center=False)

        if self.mask is None:
            return img
        else:
            img.array[mask] = 0.
            return img


def add_noise(image, snr, noise_seed=0):
    """Set gaussian noise to the given galsim.Image.

    Args:
        image(:class:`galsim.Image`): Galaxy image that noise is going to be added
                             to.
        snr(float): Signal to noise ratio.
        noise_seed(int): Seed to set to galsim.BaseDeviate which
                         will create the galsim.noise instance.

    Returns:
        A :class:`galsim.Image`, variance_noise tuple. The image is the noisy version
        of the original image and variance_noise is the noise variance on each
        pixel due to the added noise.

    """

    noisy_image = deepcopy(image)  # do not alter original image.
    bd = galsim.BaseDeviate(noise_seed)
    noise = galsim.GaussianNoise(rng=bd)
    variance_noise = noisy_image.addNoiseSNR(noise, snr, preserve_flux=True)
    return noisy_image, variance_noise


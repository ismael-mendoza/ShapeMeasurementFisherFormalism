Models and Parameters
=====================
This document describes the galaxy parameters and models that are implemented in the `generate.py` program. For more detailed information 
take a look into :mod:`analysis.models` and the `galsim <https://github.com/GalSim-developers/GalSim>`_ documentation. 

Galaxy Models
--------------

============= ====================================================================================
Name           Description
============= ====================================================================================
gaussian      Gaussian galaxy profile, :py:class:`galsim.Gaussian`. 
exponential   Exponential galaxy profile, :py:class:`galsim.Exponential`. 
deVaucouleurs DeVaucouleurs galaxy profile, :py:class:`galsim.DeVaucouleurs`.
bulgeDisk     Combined exponential (disk) and deVaucouleurs (bulge) profile. 
============= ====================================================================================

Galaxy Parameters
------------------

======== ======= ====================================================================================
Name     Type    Description
======== ======= ====================================================================================
-------- ------- ------------------------------------------------------------------------------------
**Program Parameters**
-----------------------------------------------------------------------------------------------------
project  string  Pixel offset of left edge of bounding box relative to left edge of survey image.
gal      int32   Pixel offset of right edge of bounding box relative to left edge of survey image.
ymin     int32   Pixel offset of bottom edge of bounding box relative to bottom edge of survey image.
ymax     int32   Pixel offset of top edge of bounding box relative to bottom edge of survey image.
-------- ------- ------------------------------------------------------------------------------------
**Source Properties**
-----------------------------------------------------------------------------------------------------
f_disk   float32 Fraction of total galaxy flux to due a Sersic n=1 disk component.
f_bulge  float32 Fraction of total galaxy flux to due a Sersic n=4 bulge component.
flux     float32 Total detected flux in electrons.
x        float32 Source centroid in x relative to image center in arcseconds.
y        float32 Source centroid in y relative to image center in arcseconds.
e1       float32 Real part (+) of galaxy ellipticity spinor (Q11-Q22+2*Q12)/(Q11+Q22).
e2       float32 Imaginary part (x) of galaxy ellipticity spinor 2*Q12/(Q11+Q22).
g1       float32 Real part (+) of galaxy ellipticity spinor (Q11-Q22)/(Q11+Q22+2\|Q\|**0.5).
g2       float32 Imaginary part (x) of galaxy ellipticity spinor (2*Q12)/(Q11+Q22+2\|Q\|**0.5).
eta1     float32 Conformal shear with a/b = exp(eta).
eta2     float32 Conformal shear with a/b = exp(eta).
sigma    float32 Galaxy size arcseconds calculated as \|Q\|**0.25.
sigma_p  float32 Galaxy size in arcseconds calculated as (0.5*trQ)**0.5.
hlr      float32 The half-light-radius of the profile. 
beta     float32 Position angle of second-moment ellipse in radians, or zero when a = b.
q        float32 b/a ratio of semi-minor axis to semi-major axis. 
fwhm     float32 Full-width-half-max of the profile, defined as radius at half the max peak... 
snr      float32 Signal to noise ration as defined in... 
hlr_b    float32 hlr for the bulge component of the :class:`analysis.models.bulgeDisk` model.
flux_d   float32 flux for the disk component of the :class:`analysis.models.bulgeDisk` model.
flux_b   float32 flux for the bulge component of the :class:`analysis.models.bulgeDisk` model.
R_r      float32 hlr for the bulge component of the :class:`analysis.models.bulgeDisk` model.
======== ======= ====================================================================================

PSF Models
-----------

============= ====================================================================================
Name          Description
============= ====================================================================================
psf_gaussian  Gaussian galaxy profile, :py:class:`galsim.Gaussian`. 
psf_moffat    Moffat galaxy profile, :py:class:`galsim.Moffat`. 
============= ====================================================================================

PSF Parameters 
--------------
For the PSF, the parameters have the same significance as in the table of parameters as above. 
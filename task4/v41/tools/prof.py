# Tools to analytically manipulate ellipses and some astronomical profiles.
#
# Written by Joshua E. Meyers (2014-2015).

import numpy as np
from scipy.optimize import newton
from scipy.special import gammainc, gamma
import math

# Make extensive use of lazy properties, so that we don't calculate anything until
# we absolutely need to, and can access computed attributes through the Python
# `.` syntax.
def lazyprop(fn):
    attr_name = '_lazy_' + fn.__name__
    @property
    def _lazyprop(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazyprop

class Ellipticity(object):
    """Class to represent different parameterizations of ellipticity.  The different known
    variables are:

    q -- the ratio of semiminor to semimajor axes.
    g -- reduced shear |g| = (1-q) / (1+q)   (corresponds to Schneider epsilon)
    e -- distortion |e| = (1-q^2) / (1+q^2)  (corresponds to Schneider chi)
    eta -- conformal shear q = exp(-|eta|)

    The latter three above options may be entered either as components, or with an ellipticity
    phase `theta`.  I.e., the following are the same:

    >>> ellipse = Ellipticity(e1=0.1, e2=0.1)
    >>> ellipse = Ellipticity(e=0.1*(2**0.5), theta=np.pi/4)

    Note that the ellipticity phase is *not* the same as the position angle of the ellipse.  In
    fact, PA = theta/2.0.
    """
    def __init__(self, **kwargs):
        # always immediately set e1 and e2, which are normal class attributes
        # everything else is a lazy property
        if len(kwargs) == 0:
            self.e1 = 0.0
            self.e2 = 0.0
        if 'e' in kwargs:
            self._lazy_e = kwargs.pop('e')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        elif 'g' in kwargs:
            self._lazy_g = kwargs.pop('g')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self._lazy_e = 2.0*self.g / (1 + self.g**2)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        elif 'eta' in kwargs:
            self._lazy_eta = kwargs.pop('eta')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self._lazy_e = math.tanh(self.eta/2.0)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        elif 'q' in kwargs:
            self._lazy_q = kwargs.pop('q')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self._lazy_e = (1.0 - self.q**2)/(1.0 + self.q**2)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        elif 'e1' in kwargs:
            self.e1 = kwargs.pop('e1')
            self.e2 = kwargs.pop('e2', 0.0)
        elif 'e2' in kwargs: # didn't find an e1
            self.e2 = kwargs.pop('e2')
            self.e1 = 0.0
        elif 'g1' in kwargs:
            self._lazy_g1 = kwargs.pop('g1')
            self._lazy_g2 = kwargs.pop('g2', 0.0)
            self._lazy_g = math.sqrt(self.g1**2 + self.g2**2)
            self._lazy_e = 2.0*self.g / (1 + self.g**2)
            self._lazy_theta = math.atan2(self.g2, self.g1)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        elif 'g2' in kwargs: #didn't find a g1
            self._lazy_g2 = kwargs.pop('g2')
            self._lazy_g1 = 0.0
            self._lazy_g = math.sqrt(self.g1**2 + self.g2**2)
            self._lazy_e = 2.0*self.g / (1 + self.g**2)
            self._lazy_theta = math.atan2(self.g2, self.g1)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        elif 'eta1' in kwargs:
            self._lazy_eta1 = kwargs.pop('eta1')
            self._lazy_eta2 = kwargs.pop('eta2', 0.0)
            self._lazy_eta = math.sqrt(self.eta1**2 + self.eta2**2)
            self._lazy_e = math.tanh(self.eta/2.0)
            self._lazy_theta = math.atan2(self.eta2, self.eta1)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        elif 'eta2' in kwargs: #didn't find an eta1
            self._lazy_eta2 = kwargs.pop('eta2')
            self._lazy_eta1 = 0.0
            self._lazy_eta = math.sqrt(self.eta1**2 + self.eta2**2)
            self._lazy_e = math.tanh(self.eta/2.0)
            self._lazy_theta = math.atan2(self.eta2, self.eta1)
            self.e1 = self.e * math.cos(self.theta)
            self.e2 = self.e * math.sin(self.theta)
        if 'theta' in kwargs: # bare theta passed
            self._lazy_theta = kwargs.pop('theta')
            self.e1 = 0.0
            self.e2 = 0.0
        if kwargs:
            raise ValueError('Too many arguments to Ellipticity')

    def __repr__(self):
        return (self.__class__.__name__+"(e1="+str(self.e1)+", e2="+str(self.e2)+")")

    def __str__(self):
        return ("(e1="+str(self.e1)+", e2="+str(self.e2)+")")

    @lazyprop
    def e(self):
        return math.sqrt(self.e1**2 + self.e2**2)

    @lazyprop
    def g(self):
        return self.e / (1.0 + math.sqrt(1.0 - self.e**2))

    @lazyprop
    def eta(self):
        return 2.0 * math.atanh(self.e)

    @lazyprop
    def q(self):
        return math.sqrt((1.0 - self.e)/(1.0 + self.e))

    @lazyprop
    def theta(self):
        return math.atan2(self.e2, self.e1)

    @lazyprop
    def g1(self):
        return self.e1 / (1.0 + math.sqrt(1.0 - self.e**2))

    @lazyprop
    def g2(self):
        return self.e2 / (1.0 + math.sqrt(1.0 - self.e**2))

    @lazyprop
    def eta1(self):
        return self.eta * math.cos(self.theta)

    @lazyprop
    def eta2(self):
        return self.eta * math.sin(self.theta)


# Abstract Base Class for defining arbitrary radial profiles getting elliptified.
# Derived classes must define methods for:
# FWHM
# rsqr
# peak
# flux
class Profile(object):
    def __init__(self, **kwargs):
        # centroid attributes
        self.x0 = kwargs.pop('x0', 0.0)
        self.y0 = kwargs.pop('y0', 0.0)

        # both flux and peak are lazy properties.
        if 'peak' in kwargs:
            self._lazy_peak = kwargs.pop('peak')
        else:
            self._lazy_flux = kwargs.pop('flux', 1.0)

        # Here's the tricky part.  There are many ways to parameterize the ellipticity+size
        # of the profile.  We'll make a, b, phi attributes and everything else lazy properties.
        if 'a' in kwargs and 'b' in kwargs:
            self.a = kwargs.pop('a')*1.0
            self.b = kwargs.pop('b')*1.0
            self.phi = kwargs.pop('phi', 0.0)
        elif 'C' in kwargs: # inverse covariance matrix analog, see Voigt++12.
            self._lazy_C = kwargs.pop('C')
            one_over_a_squared = 0.5 * (self.C[0,0] + self.C[1,1]
                                        + math.sqrt((self.C[0,0] - self.C[1,1])**2
                                                    + 4.0 * self.C[0,1]**2))
            one_over_b_squared = self.C[0,0] + self.C[1,1] - one_over_a_squared
            # there's degeneracy between a, b and phi at this point so enforce a > b
            if one_over_a_squared > one_over_b_squared:
                one_over_a_squared, one_over_b_squared = one_over_b_squared, one_over_a_squared
            self.a = math.sqrt(1.0 / one_over_a_squared)
            self.b = math.sqrt(1.0 / one_over_b_squared)
            if self.a == self.b:
                self.phi = 0.0
            else:
                self.phi = 0.5 * math.atan2(
                    2.0 * self.C[0,1] / (one_over_a_squared - one_over_b_squared),
                    (self.C[0,0] - self.C[1,1]) / (one_over_a_squared - one_over_b_squared))
        elif 'I' in kwargs: # second moment matrix
            self._lazy_I = kwargs.pop('I')
            self.phi = 0.5 * math.atan2(self.I[0,0] - self.I[1,1], 2.0 * self.I[0,1])
            self._lazy_rsqr = 2.0*math.sqrt(np.linalg.det(self.I))
            e = (math.sqrt((self.I[0,0] - self.I[1,1])**2 + 4.0*self.I[0,1]**2)
                 / (self.I[0,0] + self.I[1,1]))
            q = math.sqrt((1.0-e)/(1.0+e))
            self.a = self.half_light_radius / math.sqrt(q)
            self.b = self.half_light_radius * math.sqrt(q)
        else:
            if 'half_light_radius' in kwargs:
                self._lazy_half_light_radius = kwargs.pop('half_light_radius')
            elif 'hlr' in kwargs:
                self._lazy_half_light_radius = kwargs.pop('hlr')
            elif 'HLR' in kwargs:
                self._lazy_half_light_radius = kwargs.pop('HLR')
            elif 'FWHM' in kwargs:
                self._lazy_FWHM = kwargs.pop('FWHM')
            elif 'fwhm' in kwargs:
                self._lazy_FWHM = kwargs.pop('fwhm')
            elif 'rsqr' in kwargs:
                self._lazy_rsqr = kwargs.pop('rsqr')
            else:
                raise ValueError('Profile size not specified')
            if 'phi' in kwargs:
                self.phi = kwargs.pop('phi')
                theta = 2.0*self.phi
                self._lazy_ellip = Ellipticity(theta=theta, **kwargs)
            else:
                self._lazy_ellip = Ellipticity(**kwargs)
                self.phi = 0.5*self.theta
            self.a = self.half_light_radius / math.sqrt(self.q)
            self.b = self.half_light_radius * math.sqrt(self.q)

    @lazyprop
    def Ixx(self):
        return self.I[0, 0]

    @lazyprop
    def Ixy(self):
        return self.I[0, 1]

    @lazyprop
    def Iyy(self):
        return self.I[1, 1]

    @lazyprop
    def half_light_radius(self):
        return math.sqrt(self.a * self.b)

    @lazyprop
    def HLR(self):
        return self.half_light_radius

    @lazyprop
    def hlr(self):
        return self.half_light_radius

    @lazyprop
    def FWHM(self):
        raise NotImplementedError('Derived class must define FWHM property getter')

    @lazyprop
    def fwhm(self):
        return self.FWHM

    @lazyprop
    def rsqr(self):
        raise NotImplementedError('Derived class must define rsqr property getter')

    @lazyprop
    def ellip(self):
        return Ellipticity(q=self.b/self.a, theta=2.0*self.phi)

    @lazyprop
    def peak(self):
        raise NotImplementedError('Derived class must define peak property getter')

    @lazyprop
    def flux(self):
        raise NotImplementedError('Derived class must define flux property getter')

    @lazyprop
    def C(self):
        C = np.matrix(np.identity(2, dtype=float))
        cph = math.cos(self.phi)
        sph = math.sin(self.phi)
        C[0,0] = (cph/self.a)**2 + (sph/self.b)**2
        C[0,1] = C[1,0] = 0.5 * (1.0/self.a**2 - 1.0/self.b**2) * math.sin(2.0 * self.phi)
        C[1,1] = (sph/self.a)**2 + (cph/self.b)**2
        return C

    @lazyprop
    def I(self):
        unrotI = 0.5 * self.rsqr * np.matrix([[1./self.q, 0.0],
                                              [0.0, self.q]], dtype=float)
        cph = math.cos(self.phi)
        sph = math.sin(self.phi)
        R = np.matrix([[cph, -sph],
                       [sph, cph]], dtype=float)
        return R * unrotI * R.T

    @lazyprop
    def e1(self):
        return self.ellip.e1

    @lazyprop
    def e2(self):
        return self.ellip.e2

    @lazyprop
    def e(self):
        return self.ellip.e

    @lazyprop
    def g1(self):
        return self.ellip.g1

    @lazyprop
    def g2(self):
        return self.ellip.g2

    @lazyprop
    def g(self):
        return self.ellip.g

    @lazyprop
    def eta1(self):
        return self.ellip.eta1

    @lazyprop
    def eta2(self):
        return self.ellip.eta2

    @lazyprop
    def eta(self):
        return self.ellip.eta

    @lazyprop
    def q(self):
        return self.ellip.q

    @lazyprop
    def theta(self):
        return self.ellip.theta

class Sersic(Profile):
    """Class to represent a Sersic profile.  Can be initialized many different ways.
    Known variables include:

    profile shape
    -------------
    `n` - the Sersic index.  This is required.

    amplitude
    ---------
    `peak` - the maximum value of the central part of the profile
    `flux` - the value of the spatial integrand of the profile
    These are related; only one should be specified.  Default is flux=1.

    size
    ----
    half_light_radius -- duh.
    `FWHM` -- full width at half maximum
    `rsqr` -- the second-moment squared radius, defined as Ixx + Iyy.
    `a`, `b` -- the semi-major and semi-minor axes.  Note that together these also implicitly
                specify the magnitude of the profile ellipticity.

    ellipticity
    -----------
    You can use any of the parameters used to create an Ellipse object.
    Additionally, you can use:
    `a`, `b` -- as mentioned above, these define the axis ratio `q`.
    `phi` -- the position angle of the major axis of the ellipse, measured CCW from +x.

    """
    def __init__(self, **kwargs):
        # Sersic index
        if 'n' not in kwargs:
            raise ValueError('Sersic initialization requires `n` index')
        self.n = kwargs.pop('n')
        super(Sersic, self).__init__(**kwargs)

    @lazyprop
    def half_light_radius(self):
        if hasattr(self, 'a') and hasattr(self, 'b'):
            return super(Sersic, self).half_light_radius
        elif hasattr(self, '_lazy_rsqr'):
            return math.sqrt(self.rsqr / self.HLR_to_rsqr(1.0, self.n))
        elif hasattr(self, '_lazy_FWHM'):
            return self.FWHM / self.HLR_to_FWHM(1.0, self.n, self.kappa)

    @lazyprop
    def kappa(self):
        return self.n_to_kappa(self.n)

    @lazyprop
    def FWHM(self):
        return self.HLR_to_FWHM(self.half_light_radius, self.n, self.kappa)

    @lazyprop
    def rsqr(self):
        return self.HLR_to_rsqr(self.half_light_radius, self.n)

    @lazyprop
    def peak(self):
        return self.flux_to_peak(self.flux, self.n, self.half_light_radius, self.kappa)

    @lazyprop
    def flux(self):
        return self.peak_to_flux(self.peak, self.n, self.half_light_radius, self.kappa)

    def __call__(self, x, y):
        x = np.array(x)
        y = np.array(y)
        xvec = np.array([x-self.x0, y-self.y0]).T
        if len(xvec.shape) == 1:
            exponent = np.einsum('j,jk,k', xvec, self.C, xvec)
        else:
            exponent = np.einsum('ij,jk,ij->i', xvec, self.C, xvec)
        exponent **= 0.5 / self.n
        exponent *= -self.kappa
        return self.peak * np.exp(exponent)

    @staticmethod
    def n_to_kappa(n):
        '''Compute Sersic exponent factor kappa from the Sersic index'''
        kguess = 1.9992 * n - 0.3271
        return newton(lambda k: gammainc(2.0 * n, k) - 0.5, kguess)

    @classmethod
    def HLR_to_FWHM(cls, half_light_radius, n, kappa=None):
        '''Compute the full-width at half maximum given input parameters.'''
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        return 2.0 * half_light_radius * (math.log(2.0) / kappa)**n

    @classmethod
    def FWHM_to_HLR(cls, FWHM, n, kappa=None):
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        return 0.5 * FWHM * (kappa / math.log(2.0))**n

    @staticmethod
    def HLR_to_rsqr(half_light_radius, n):
        """ Using Mathematica-generated fitting function, compute the circularized second-moment
        squared radius rsqr = Ixx + Iyy; where Ixx = Iyy and Ixy = 0.0
        """
        return (half_light_radius * (0.985444 + n * (0.391016 + n * (0.0739602
                                     + n * (0.00698719 + n * (0.00212432
                                     + n * (-0.000154052 + n * 0.0000219632)))))))**2

    @staticmethod
    def rsqr_to_HLR(cls, rsqr, n):
        return math.sqrt(rsqr / self.HLR_to_rsqr(1.0, n))

    @classmethod
    def flux_to_peak(cls, flux, n, half_light_radius, kappa=None):
        '''Compute peak surface brightness from integrated flux and Sersic params.'''
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        fluxnorm = cls.peak_to_flux(1.0, n, half_light_radius, kappa=kappa)
        return flux/fluxnorm

    @classmethod
    def peak_to_flux(cls, peak, n, half_light_radius, kappa=None):
        '''Compute flux integrated over all space given input parameters.'''
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        return (2.0 * math.pi * peak * n * (kappa**(-2.0 * n)) * (half_light_radius**2.0)
                * gamma(2.0 * n))


class Gaussian(Sersic):
    """A Gaussian is a Sersic with index n=0.5.  See Sersic documentation for class details.
    """
    def __init__(self, sigma=None, **kwargs):
        if sigma is not None:
            kwargs.update({'rsqr':2.0*sigma**2})
            self._lazy_sigma = sigma
        super(Gaussian, self).__init__(n=0.5, **kwargs)

    @lazyprop
    def sigma(self):
        return math.sqrt(0.5 * self.rsqr)


class Exponential(Sersic):
    """An exponential is a Sersic with index n=1.0.  See Sersic documentation for class details.
    """
    def __init__(self, **kwargs):
        super(Exponential, self).__init__(n=1.0, **kwargs)


class DeVaucouleurs(Sersic):
    """A DeVaucouleurs is a Sersic with index n=4.0.  See Sersic documentation for class details.
    """
    def __init__(self, **kwargs):
        super(DeVaucouleurs, self).__init__(n=4.0, **kwargs)


class Moffat(Profile):
    """Class to represent a Moffat profile.  Can be initialized many different ways.
    Known variables include:

    profile shape
    -------------
    `beta` - the Moffat parameter.  This is required.

    amplitude
    ---------
    `peak` - the maximum value of the central part of the profile
    `flux` - the value of the spatial integrand of the profile
    These are related; only one should be specified.  Default is flux=1.

    size
    ----
    half_light_radius -- duh.
    `FWHM` -- full width at half maximum
    `rsqr` -- the second moment radius, defined as Ixx + Iyy.
    `a`, `b` -- the semi-major and semi-minor axes.  Note that together these also implicitly
                specify the magnitude of the profile ellipticity.

    ellipticity
    -----------
    You can use any of the parameters used to create an Ellipse object.
    Additionally, you can use:
    `a`, `b` -- as mentioned above, these define the axis ratio `q`.
    `phi` -- the position angle of the major axis of the ellipse, measured CCW from +x.

    """
    def __init__(self, **kwargs):
        # Moffat beta parameter
        if 'beta' not in kwargs:
            raise ValueError('Moffat initialization requires `beta` index')
        self.beta = kwargs.pop('beta')
        super(Moffat, self).__init__(**kwargs)

    @lazyprop
    def half_light_radius(self):
        if hasattr(self, 'a') and hasattr(self, 'b'):
            return super(Moffat, self).half_light_radius
        elif hasattr(self, '_lazy_rsqr'):
            return math.sqrt(self.rsqr / self.HLR_to_rsqr(1.0, self.beta))
        elif hasattr(self, '_lazy_FWHM'):
            return self.FWHM / self.HLR_to_FWHM(1.0, self.beta)

    @lazyprop
    def FWHM(self):
        return self.HLR_to_FWHM(self.half_light_radius, self.beta)

    @lazyprop
    def rsqr(self):
        return self.HLR_to_rsqr(self.half_light_radius, self.beta)

    @lazyprop
    def peak(self):
        return self.flux_to_peak(self.flux, self.beta, self.C)

    @lazyprop
    def flux(self):
        return self.peak_to_flux(self.peak, self.beta, self.C)

    def __call__(self, x, y):
        x = np.array(x)
        y = np.array(y)
        xvec = np.array([x-self.x0, y-self.y0]).T
        if len(xvec.shape) == 1:
            base = np.einsum('j,jk,k', xvec, self.C, xvec)
        else:
            base = np.einsum('ij,jk,ij->i', xvec, self.C, xvec)
        base *= (2**(1.0 / (self.beta - 1.0)) - 1.0)
        return self.peak * (1+base)**(-self.beta)

    @staticmethod
    def HLR_to_FWHM(half_light_radius, beta):
        return 2.0 * half_light_radius * math.sqrt(
            (2.0**(1.0/beta) - 1.0) * (2.0**(1.0/(beta-1.0)) - 1.0))

    @staticmethod
    def FWHM_to_HLR(FWHM, beta):
        return 0.5 * FWHM / math.sqrt(
            (2.0**(1.0/beta) - 1.0) * (2.0**(1.0/(beta-1.0)) - 1.0))

    @classmethod
    def HLR_to_rsqr(cls, half_light_radius, beta):
        return ((0.5 * cls.HLR_to_FWHM(half_light_radius, beta))**2 /
                (2.0**(1.0/beta) - 1.0) * (beta - 2.0))

    @classmethod
    def rsqr_to_HLR(cls, rsqr, beta):
        return cls.FWHM_to_HLR(2.0 * math.sqrt(rsqr * (2.0**(1.0/beta) - 1.0)*(beta - 2.0)), beta)

    @staticmethod
    def flux_to_peak(flux, beta, C):
        factor = (beta - 1.0) / (math.pi / math.sqrt(abs(np.linalg.det(C))))
        factor *= (2**(1.0 / (beta - 1.0)) - 1.0)
        return flux * factor

    @staticmethod
    def peak_to_flux(peak, beta, C):
        factor = (beta - 1.0) / (math.pi / math.sqrt(abs(np.linalg.det(C))))
        factor *= (2**(1.0 / (beta - 1.0)) - 1.0)
        return peak / factor


class Sum(object):
    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
        self.flux = reduce(lambda x,y:x+y, [a.flux for a in args], 0)
        self.x0 = 1./self.flux * reduce(lambda x,y:x+y, [a.flux * a.x0 for a in args], 0)
        self.y0 = 1./self.flux * reduce(lambda x,y:x+y, [a.flux * a.y0 for a in args], 0)
        self.Ixx = 1./self.flux * reduce(lambda x,y:x+y, [a.flux * (a.Ixx + (a.x0-self.x0)**2)
                                                          for a in args], 0)
        self.Ixy = 1./self.flux * reduce(lambda x,y:x+y, [a.flux * (a.Ixy
                                                                    + (a.x0-self.y0)*(a.y0-self.y0))
                                                          for a in args], 0)
        self.Iyy = 1./self.flux * reduce(lambda x,y:x+y, [a.flux * (a.Iyy + (a.y0-self.y0)**2)
                                                          for a in args], 0)

class Convolution(object):
    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
        self.flux = reduce(lambda x,y:x*y, [a.flux for a in args], 1)
        self.x0 = reduce(lambda x,y:x+y, [a.x0 for a in args], 0)
        self.y0 = reduce(lambda x,y:x+y, [a.y0 for a in args], 0)
        self.Ixx = reduce(lambda x,y:x+y, [a.Ixx for a in args], 0)
        self.Ixy = reduce(lambda x,y:x+y, [a.Ixy for a in args], 0)
        self.Iyy = reduce(lambda x,y:x+y, [a.Iyy for a in args], 0)

def compare_Ellipticities(a, b):
    np.testing.assert_almost_equal(a.e, b.e)
    np.testing.assert_almost_equal(a.e1, b.e1)
    np.testing.assert_almost_equal(a.e2, b.e2)
    np.testing.assert_almost_equal(a.g, b.g)
    np.testing.assert_almost_equal(a.g1, b.g1)
    np.testing.assert_almost_equal(a.g2, b.g2)
    np.testing.assert_almost_equal(a.eta, b.eta)
    np.testing.assert_almost_equal(a.eta1, b.eta1)
    np.testing.assert_almost_equal(a.eta2, b.eta2)
    np.testing.assert_almost_equal(a.theta, b.theta)
    np.testing.assert_almost_equal(a.q, b.q)

def test_Ellipticity():
    a = Ellipticity(e1=0.11, e2=-0.023)

    b = Ellipticity(e=a.e, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(g=a.g, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(eta=a.eta, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(q=a.q, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(e1=a.e1, e2=a.e2)
    compare_Ellipticities(a, b)

    b = Ellipticity(g1=a.g1, g2=a.g2)
    compare_Ellipticities(a, b)

    b = Ellipticity(eta1=a.eta1, eta2=a.eta2)
    compare_Ellipticities(a, b)

def compare_Sersics(a, b):
    np.testing.assert_almost_equal(a.x0, b.x0)
    np.testing.assert_almost_equal(a.y0, b.y0)
    np.testing.assert_almost_equal(a.n, b.n)
    np.testing.assert_almost_equal(a.kappa, b.kappa)
    np.testing.assert_almost_equal(a.a, b.a)
    np.testing.assert_almost_equal(a.b, b.b)
    np.testing.assert_array_almost_equal(a.I, b.I)
    np.testing.assert_array_almost_equal(a.C, b.C)
    np.testing.assert_almost_equal(a.e, b.e)
    np.testing.assert_almost_equal(a.e1, b.e1)
    np.testing.assert_almost_equal(a.e2, b.e2)
    np.testing.assert_almost_equal(a.g, b.g)
    np.testing.assert_almost_equal(a.g1, b.g1)
    np.testing.assert_almost_equal(a.g2, b.g2)
    np.testing.assert_almost_equal(a.eta, b.eta)
    np.testing.assert_almost_equal(a.eta1, b.eta1)
    np.testing.assert_almost_equal(a.eta2, b.eta2)
    np.testing.assert_almost_equal(a.theta, b.theta)
    np.testing.assert_almost_equal(a.q, b.q)
    np.testing.assert_almost_equal(a.FWHM, b.FWHM)
    np.testing.assert_almost_equal(a.rsqr, b.rsqr)
    np.testing.assert_almost_equal(a.half_light_radius, b.half_light_radius)
    np.testing.assert_almost_equal(a.flux, b.flux)
    np.testing.assert_almost_equal(a.peak, b.peak)

def test_Sersic():
    # test radial part
    a = Sersic(n=2.2, FWHM=1.1)

    b = Sersic(n=2.2, rsqr=a.rsqr)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=a.FWHM)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, half_light_radius=a.half_light_radius)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, I=a.I)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, C=a.C)
    compare_Sersics(a, b)

    # test flux part
    a = Sersic(n=2.2, FWHM=1.0, flux=1.1)
    b = Sersic(n=2.2, FWHM=1.0, peak=a.peak)
    compare_Sersics(a, b)

    # test ellipticity part
    a = Sersic(n=2.2, FWHM=1.0, e1=0.21, e2=-0.034)
    b = Sersic(n=2.2, FWHM=1.0, eta1=a.eta1, eta2=a.eta2)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, g1=a.g1, g2=a.g2)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, e1=a.e1, e2=a.e2)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, e=a.e, theta=a.theta)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, g=a.g, theta=a.theta)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, eta=a.eta, theta=a.theta)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, q=a.q, theta=a.theta)
    compare_Sersics(a, b)

    #test half-light-radius
    from scipy.integrate import quad
    for n in np.linspace(0.2, 8.0, 30):
        a = Sersic(n=n, half_light_radius=1.0)
        int_hlr = quad(lambda x:a(0,x)*x*2*math.pi, 0, 1.0)[0]
        np.testing.assert_almost_equal(int_hlr, 0.5)

def compare_Moffats(a, b):
    np.testing.assert_almost_equal(a.x0, b.x0)
    np.testing.assert_almost_equal(a.y0, b.y0)
    np.testing.assert_almost_equal(a.beta, b.beta)
    np.testing.assert_almost_equal(a.a, b.a)
    np.testing.assert_almost_equal(a.b, b.b)
    np.testing.assert_array_almost_equal(a.I, b.I)
    np.testing.assert_array_almost_equal(a.C, b.C)
    np.testing.assert_almost_equal(a.e, b.e)
    np.testing.assert_almost_equal(a.e1, b.e1)
    np.testing.assert_almost_equal(a.e2, b.e2)
    np.testing.assert_almost_equal(a.g, b.g)
    np.testing.assert_almost_equal(a.g1, b.g1)
    np.testing.assert_almost_equal(a.g2, b.g2)
    np.testing.assert_almost_equal(a.eta, b.eta)
    np.testing.assert_almost_equal(a.eta1, b.eta1)
    np.testing.assert_almost_equal(a.eta2, b.eta2)
    np.testing.assert_almost_equal(a.theta, b.theta)
    np.testing.assert_almost_equal(a.q, b.q)
    np.testing.assert_almost_equal(a.FWHM, b.FWHM)
    np.testing.assert_almost_equal(a.rsqr, b.rsqr)
    np.testing.assert_almost_equal(a.half_light_radius, b.half_light_radius)
    np.testing.assert_almost_equal(a.flux, b.flux)
    np.testing.assert_almost_equal(a.peak, b.peak)

def test_Moffat():
    # test radial part
    a = Moffat(beta=3.0, FWHM=1.1)

    b = Moffat(beta=3.0, rsqr=a.rsqr)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=a.FWHM)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, half_light_radius=a.half_light_radius)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, I=a.I)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, C=a.C)
    compare_Moffats(a, b)

    # test flux part
    a = Moffat(beta=3.0, FWHM=1.0, flux=1.1)
    b = Moffat(beta=3.0, FWHM=1.0, peak=a.peak)
    compare_Moffats(a, b)

    # test ellipticity part
    a = Moffat(beta=3.0, FWHM=1.0, e1=0.21, e2=-0.034)
    b = Moffat(beta=3.0, FWHM=1.0, eta1=a.eta1, eta2=a.eta2)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, g1=a.g1, g2=a.g2)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, e1=a.e1, e2=a.e2)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, e=a.e, theta=a.theta)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, g=a.g, theta=a.theta)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, eta=a.eta, theta=a.theta)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, q=a.q, theta=a.theta)
    compare_Moffats(a, b)

    #test half-light-radius
    from scipy.integrate import quad
    for beta in np.linspace(2.5, 8.0, 30):
        a = Moffat(beta=beta, half_light_radius=1.0)
        int_hlr = quad(lambda x:a(0,x)*x*2*math.pi, 0, 1.0)[0]
        np.testing.assert_almost_equal(int_hlr, 0.5)

def test_galsim():
    try:
        import galsim
    except:
        return

    a = Moffat(half_light_radius=1.0, beta=3.0, e1=0.1, e2=-0.2)
    b = galsim.Moffat(half_light_radius=1.0, beta=3.0).shear(e1=0.1, e2=-0.2)
    np.testing.assert_almost_equal(a(0,0), b.xValue(0,0))
    np.testing.assert_almost_equal(a(0,1), b.xValue(0,1))
    np.testing.assert_almost_equal(a(2,1), b.xValue(2,1))

    a = Sersic(half_light_radius=1.0, n=3.0, e1=0.1, e2=-0.2)
    b = galsim.Sersic(half_light_radius=1.0, n=3.0).shear(e1=0.1, e2=-0.2)
    np.testing.assert_almost_equal(a(0,0), b.xValue(0,0))
    np.testing.assert_almost_equal(a(0,1), b.xValue(0,1))
    np.testing.assert_almost_equal(a(2,1), b.xValue(2,1))

if __name__ == '__main__':
    test_galsim()
    test_Ellipticity()
    test_Sersic()
    test_Moffat()

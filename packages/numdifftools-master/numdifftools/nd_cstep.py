"""numerical differentiation functions:
Derivative, Gradient, Jacobian, and Hessian

Author : pbrod, josef-pkt
License : BSD
Notes
-----
These are simple forward differentiation, so that we have them available
without dependencies.
* Jacobian should be faster than numdifftools.core because it doesn't use loop
  over observations.
* numerical precision will vary and depend on the choice of stepsizes
"""

# TODO:
# * some cleanup
# * check numerical accuracy (and bugs) with numdifftools and analytical
#   derivatives
#   - linear least squares case: (hess - 2*X'X) is 1e-8 or so
#   - gradient and Hessian agree with numdifftools when evaluated away from
#     minimum
#   - forward gradient, Jacobian evaluated at minimum is inaccurate, centered
#     (+/- base_step) is ok
# * dot product of Jacobian is different from Hessian, either wrong example or
#   a bug (unlikely), or a real difference
#
#
# What are the conditions that Jacobian dotproduct and Hessian are the same?
#
# See also:
#
# BHHH: Greene p481 17.4.6,  MLE Jacobian = d loglike / d beta , where loglike
# is vector for each observation
#    see also example 17.4 when J'J is very different from Hessian
#    also does it hold only at the minimum, what's relationship to covariance
#    of Jacobian matrix
# http://projects.scipy.org/scipy/ticket/1157
# http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
#    objective: sum((y-f(beta,x)**2),   Jacobian = d f/d beta
#    and not d objective/d beta as in MLE Greene similar:
# http://crsouza.blogspot.com/2009/11/neural-network-learning-by-levenberg_18.html#hessian
#
# in example: if J = d x*beta / d beta then J'J == X'X
#    similar to
#    http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
from __future__ import print_function
import numpy as np
from numdifftools.core import dea3
from collections import namedtuple
from matplotlib import pyplot as plt
from numdifftools.multicomplex import bicomplex
from numdifftools.test_functions import get_test_function, function_names
from numpy import linalg
from scipy import misc
from scipy.ndimage.filters import convolve1d
# NOTE: we only do double precision internally so far
EPS = np.MachAr().eps


def c_atan2(x, y):
    a, b = np.real(x), np.imag(x)
    c, d = np.real(y), np.imag(y)
    return np.arctan2(a, c) + 1j * (c * b - a * d) / (a**2 + c**2)


def c_max(x, y):
    return np.where(x.real < y.real, y, x)


def c_min(x, y):
    return np.where(x.real > y.real, y, x)


def _make_exact(h):
    '''Make sure h is an exact representable number
    This is important when calculating numerical derivatives and is
    accomplished by adding 1 and then subtracting 1..
    '''
    return (h + 1.0) - 1.0


def default_scale(method='forward', n=1):
    return dict(central=3, central2=3, complex=1,
                hybrid=2).get(method, 2) + 1.25*int((n - 1))


def _default_base_step(x, scale, epsilon=None):
    if epsilon is None:
        h = (10 * EPS) ** (1. / scale) * np.maximum(np.log1p(np.abs(x)), 0.1)
    else:
        if np.isscalar(epsilon):
            h = np.ones(x.shape) * epsilon
        else:
            h = np.asarray(epsilon)
            if h.shape != x.shape:
                raise ValueError("If h is not a scalar it must have the same"
                                 " shape as x.")
    return h


_CENTRAL_WEIGHTS_AND_POINTS = {
    (1, 3): (np.array([-1, 0, 1]) / 2.0, np.arange(-1, 2)),
    (1, 5): (np.array([1, -8, 0, 8, -1]) / 12.0, np.arange(-2, 3)),
    (1, 7): (np.array([-1, 9, -45, 0, 45, -9, 1]) / 60.0, np.arange(-3, 4)),
    (1, 9): (np.array([3, -32, 168, -672, 0, 672, -168, 32, -3]) / 840.0,
             np.arange(-4, 5)),
    (2, 3): (np.array([1, -2.0, 1]), np.arange(-1, 2)),
    (2, 5): (np.array([-1, 16, -30, 16, -1]) / 12.0, np.arange(-2, 3)),
    (2, 7): (np.array([2, -27, 270, -490, 270, -27, 2]) / 180.0,
             np.arange(-3, 4)),
    (2, 9): (np.array([-9, 128, -1008, 8064, -14350,
                      8064, -1008, 128, -9]) / 5040.0,
             np.arange(-4, 5))}


def fornberg_weights_all(x, x0, M=1):
    '''
    Return finite difference weights_and_points for derivatives
    of all orders 0, 1, ..., m

    Parameters
    ----------
    x : vector, length n
        x-coordinates for grid points
    x0 : scalar
        location where approximations are to be accurate
    m : scalar integer
        highest derivative that we want to find weights_and_points for

    Returns
    -------
    C :  array, shape n x m+1
        contains coefficients for the j'th derivative in column j (0 <= j <= m)

    See also:
    ---------
    fornberg_weights

    References
    ----------
    B. Fornberg (1998)
    "Calculation of weights_and_points in finite difference formulas",
    SIAM Review 40, pp. 685-691.

    http://www.scholarpedia.org/article/Finite_difference_method
    '''
    N = len(x)
    if M >= N:
        raise ValueError('length(x) must be larger than m')

    c1, c4 = 1, x[0] - x0
    C = np.zeros((N, M + 1))
    C[0, 0] = 1
    for n in range(1, N):
        m = np.arange(0, min(n, M) + 1)
        c2, c5, c4 = 1, c4, x[n] - x0
        for v in range(n):
            c3 = x[n] - x[v]
            c2, c6, c7 = c2 * c3, m * C[v, m-1], C[v, m]
            C[v, m] = (c4 * c7 - c6) / c3
        else:
            C[n, m] = c1 * (c6 - c5 * c7) / c2
        c1 = c2
    return C


def fornberg_weights(x, x0, m=1):
    '''
    Return weights for finite difference approximation of the m'th derivative
    U^m(x0), evaluated at x0, based on n values of U at x[0], x[1],... x[n-1]:

        U^m(x0) = sum weights[i] * U(x[i])

    Parameters
    ----------
    x : vector
        abscissas used for the evaluation for the derivative at x0.
    x0 : scalar
        location where approximations are to be accurate
    m : integer
        order of derivative. Note for m=0 this can be used to evaluate the
        interpolating polynomial itself.

    Notes
    -----
    The x values can be arbitrarily spaced but must be distinct and len(x) > m.

    The Fornberg algorithm is much more stable numerically than regular
    vandermonde systems for large values of n.

    See also
    --------
    fornberg_weights_all
    '''
    return fornberg_weights_all(x, x0, m)[:, -1]


_cmn_doc = """
    Calculate %(derivative)s with finite difference approximation

    Parameters
    ----------
    f : function
       function of one array f(x, `*args`, `**kwargs`)
    steps : float, array-like or StepsGenerator object, optional
       Spacing used, if None, then the spacing is automatically chosen
       according to (10*EPS)**(1/scale)*max(log(1+|x|), 0.1) where scale is
       depending on method and derivative-order (see default_scale).
       A StepsGenerator can be used to extrapolate the results. However,
       the generator must generate minimum 3 steps in order to extrapolate
       the values.
    method : string, optional
        defines method used in the approximation
        'complex': complex-step derivative (scale=%(scale_complex)s)
        'central': central difference derivative (scale=%(scale_central)s)
        'backward': backward difference derivative (scale=%(scale_backward)s)
        'forward': forward difference derivative (scale=%(scale_forward)s)
        %(extra_method)s
    full_output : bool, optional
        If `full_output` is False, only the derivative is returned.
        If `full_output` is True, then (der, r) is returned `der` is the
        derivative, and `r` is a Results object.

    Call Parameters
    ---------------
    x : array_like
       value at which function derivative is evaluated
    args : tuple
        Arguments for function `f`.
    kwds : dict
        Keyword arguments for function `f`.
    %(returns)s
    Notes
    -----
    The complex-step derivative has truncation error O(steps**2), so
    truncation error can be eliminated by choosing steps to be very small.
    The complex-step derivative avoids the problem of round-off error with
    small steps because there is no subtraction. However, the function
    needs to be analytic. This method does not work if f(x) involves non-
    analytic functions such as e.g.: abs, max, min or the derivative-order is 2
    or greater. For this reason the 'central' method is the default method.
    This method is usually very accurate, but sometimes one can only allow
    evaluation in forward or backward direction.

    Be careful in decreasing the step size too small due to round-off errors.

    %(extra_note)s
    References
    ----------
    Ridout, M.S. (2009) Statistical applications of the complex-step method
        of numerical differentiation. The American Statistician, 63, 66-74
    %(example)s
    %(see_also)s
    """


class StepsGenerator(object):
    '''
    Generates a sequence of steps

    where
        steps = base_step * step_ratio ** (np.arange(num_steps) + offset)

    Parameters
    ----------
    base_step : float, array-like, optional
       Defines the base step, if None, then base_step is set to
           (10*EPS)**(1/scale)*max(log(1+|x|), 0.1)
       where x and scale are supplied at runtime through the __call__ method.
    step_ratio : real scalar, optional, default 4
        Ratio between sequential steps generated.
        Note: Ratio > 1
    num_steps : scalar integer, optional, default 10
        defines number of steps generated. It should be larger than
        derivative_order + method_order - 1
    offset : real scalar, optional, default 0
        offset to the base step
    scale : real scalar, optional
        scale used in base step. If set to a value it will override the default
        scale.
    '''

    def __init__(self, base_step=None, step_ratio=4, num_steps=10, offset=0,
                 use_exact_steps=True, scale=None):
        self.base_step = base_step
        self.num_steps = num_steps
        self.step_ratio = step_ratio
        self.offset = offset
        self.scale = scale
        self.use_exact_steps = use_exact_steps

    def __repr__(self):
        class_name = self.__class__.__name__
        kwds = ['%s=%s' % (name, str(getattr(self, name)))
                for name in self.__dict__.keys()]
        return """%s(%s)""" % (class_name, ','.join(kwds))

    def _default_base_step(self, xi, method, n):
        scale = self.scale
        if scale is None:
            scale = default_scale(method, n)

        delta = _default_base_step(xi, scale, self.base_step)
        if self.use_exact_steps:
            return _make_exact(delta)
        return delta

    def __call__(self, x, method='forward', n=1, order=None):
        xi = np.asarray(x)
        delta = self._default_base_step(xi, method, n)

        step_ratio, offset = float(self.step_ratio), self.offset
        for i in range(int(self.num_steps-1), -1, -1):
            h = (delta * step_ratio**(i + offset))
            if (np.abs(h) > 0).all():
                yield h


class StepsGenerator2(object):
    '''
    Generates a sequence of steps

    where
        steps = logspace(log10(step_min), log10(step_max), num_steps)

    Parameters
    ----------
    step_min : float, array-like, optional
       Defines the minimim step. Default value is:
           (10*EPS)**(1/scale)*max(log(1+|x|), 0.1)
       where x and scale are supplied at runtime through the __call__ method.
    step_max : real scalar, optional
        maximum step generated. Default value is:
            exp(log(step_min) * scale / (scale + 1.5))
    num_steps : scalar integer, optional
        defines number of steps generated.
    scale : real scalar, optional
        scale used in base step. If set to a value it will override the scale
        supplied at runtime.
    '''

    def __init__(self, step_min=None, step_max=None, num_steps=10, scale=None):
        self.step_min = step_min
        self.num_steps = num_steps
        self.step_max = step_max
        self.scale = scale

    def __call__(self, x, scale=1.5):
        if self.scale is not None:
            scale = self.scale
        xi = np.asarray(x)
        step_min, step_max = self.step_min, self.step_max
        delta = _default_base_step(xi, scale, step_min)
        if step_min is None:
            step_min = (10 * EPS)**(1. / scale)
        if step_max is None:
            step_max = np.exp(np.log(step_min) * scale / (scale + 1.5))
        steps = np.logspace(0, -np.log10(step_min) + np.log10(step_max),
                            self.num_steps)[::-1]

        for step in steps:
            h = _make_exact(delta * step)
            if (np.abs(h) > 0).all():
                yield h


class _Derivative(object):

    info = namedtuple('info', ['error_estimate', 'final_step', 'index'])

    def __init__(self, f, steps=None, method='central', full_output=False,
                 n=1, order=2):
        self.n = n
        self.order = order
        self.f = f
        self.steps = self._make_callable(steps)
        self.method = method
        self.full_output = full_output

    def _make_callable(self, steps):
        if hasattr(steps, '__call__'):
            return steps
        num_steps = self.n+self.order-1 + 10*int(steps is None)
        return StepsGenerator(base_step=steps, step_ratio=4.0**(1./self.n),
                              num_steps=num_steps)

    def _get_arg_min(self, errors):
        shape = errors.shape
        arg_mins = np.nanargmin(errors, axis=0)
        min_errors = np.nanmin(errors, axis=0)
        for i, min_error in enumerate(min_errors):
            idx = np.flatnonzero(errors[:, i] == min_error)
            arg_mins[i] = idx[idx.size // 2]
        ix = np.ravel_multi_index((arg_mins, np.arange(shape[1])), shape)
        return ix

    def _extrapolate(self, results, steps, shape):
        dont_extrapolate = results.shape[0] < 3
        if dont_extrapolate:
            err = np.empty(shape)
            err.fill(np.NaN)
            der = 0.5 * (results[0] + results[-1])
            final_step = 0.5*(steps[0]+steps[-1])
            return der.reshape(shape), self.info(err, final_step, 0)

        der, errors = dea3(results[0:-2], results[1:-1], results[2:],
                           symmetric=True)
        steps = steps[2:]
        if len(der) > 2:
            der, errors = dea3(der[0:-2], der[1:-1], der[2:], symmetric=True)
            steps = steps[2:]
        ix = self._get_arg_min(errors)
        final_step = steps.flat[ix].reshape(shape)
        err = errors.flat[ix].reshape(shape)
        return der.flat[ix].reshape(shape), self.info(err, final_step, ix)

    def _get_function_name(self):
        return '_%s' % self.method

    def _get_functions(self):
        name = self._get_function_name()
        return getattr(self, name), self.f, self.steps

    def _eval_first(self, f, x, *args, **kwds):
        if self._eval_first_condition():
            return f(x, *args, **kwds)
        return 0

    def _get_finite_difference_rule(self, step_ratio):
        return np.ones((1,))

    def _vstack(self, sequence, steps):
        original_shape = sequence[0].shape
        f_del = np.vstack(list(r.ravel()) for r in sequence)
        h = np.vstack(list((np.ones(original_shape)*step).ravel())
                      for step in steps)
        if f_del.size != h.size:
            raise ValueError('fun did not return data of correct size ' +
                             '(it must be vectorized)')
        return f_del, h, original_shape

    def _apply_fd_rule(self, fd_rule, sequence, steps):
        return self._vstack(sequence, steps)

    def _derivative(self, xi, args, kwds):
        diff, f, step_generator = self._get_functions()
        steps = [step for step in step_generator(xi, self.method, self.n)]
        fxi = self._eval_first(f, xi, *args, **kwds)
        results = [diff(f, fxi, xi, h, *args, **kwds) for h in steps]
        fd_rule = self._get_finite_difference_rule(step_generator.step_ratio)
        return self._apply_fd_rule(fd_rule, results, steps)

    def __call__(self, x, *args, **kwds):
        results = self._derivative(np.asarray(x), args, kwds)
        derivative, info = self._extrapolate(*results)
        if self.full_output:
            return derivative, info
        return derivative


class Derivative(_Derivative):
    __doc__ = _cmn_doc % dict(
        derivative='n-th derivative',
        scale_backward=str(default_scale('backward')),
        scale_central=str(default_scale('central')),
        scale_complex=str(default_scale('complex')),
        scale_forward=str(default_scale('forward')),
        extra_method="",
        extra_note='', returns="""
    Returns
    -------
    der : ndarray
       array of derivatives
    """, example="""
    Examples
    --------
    >>> import numpy as np
    >>> import numdifftools.nd_cstep as ndc

    # 1'st derivative of exp(x), at x == 1

    >>> fd = ndc.Derivative(np.exp)       # 1'st derivative
    >>> np.allclose(fd(1), 2.71828183)
    True

    >>> d2 = fd([1, 2])
    >>> d2
    array([ 2.71828183,  7.3890561 ])""", see_also="""
    See also
    --------
    Gradient,
    Hessian
    """)
    """
    Find the n-th derivative of a function at a point.

    Given a function, use a difference formula with spacing `dx` to
    compute the `n`-th derivative at `x0`.

    Parameters
    ----------
    f : function
        Input function.
    x0 : float
        The point at which `n`-th derivative is found.
    dx : float, optional
        Spacing.
    method : Method of estimation.  Valid options are:
        'central', 'forward' or 'backward'.          (Default 'central')
    n : int, optional (Default 1)
        Order of the derivative.
    order : int, optional       (Default 2)
        defining order of basic method used.
        For 'central' methods, it must be an even number eg. [2,4].

    Notes
    -----
    Decreasing the step size too small can result in round-off error.

    Note on order: higher order methods will generally be more accurate,
             but may also suffer more from numerical problems. First order
             methods would usually not be recommended.
    Note on method: Central difference methods are usually the most accurate,
            but sometimes one can only allow evaluation in forward or backward
            direction.

    Examples
    --------
    >>> def f(x):
    ...     return x**3 + x**2

    >>> df = NDerivative(f)
    >>> np.allclose(df(1), 5)
    True
    >>> ddf = NDerivative(f, n=2)
    >>> np.allclose(ddf(1), 8)
    True
    """

    def _fd_mat(self, step_ratio, parity, nterms):
        ''' Return matrix for finite difference derivation.

        Parameters
        ----------
        step_ratio : real scalar
            ratio between steps in unequally spaced difference rule.
        parity : scalar, integer
            0 (one sided, all terms included but zeroth order)
            1 (only odd terms included)
            2 (only even terms included)
        nterms : scalar, integer
            number of terms
        '''
        srinv = 1.0 / step_ratio
        [i, j] = np.ogrid[0:nterms, 0:nterms]

        try:
            fact = [1, 2, 2][parity]
        except Exception as msg:
            raise ValueError('%s. Parity must be 0, 1 or 2! (%d)' % (str(msg),
                                                                     parity))
        offset = max(1, parity)
        c = 1.0 / misc.factorial(np.arange(offset, fact * nterms + 1, fact))
        mat = c[j] * srinv ** (i * (fact * j + offset))
        return np.atleast_2d(mat)

    def _get_finite_difference_rule(self, step_ratio):
        '''
        Generate finite differencing rule in advance.

        The rule is for a nominal unit step size, and will
        be scaled later to reflect the local step size.

        Member methods used
        -------------------
        _fd_mat

        Member variables used
        ---------------------
        n
        order
        method
        '''
        if self.method == 'complex':
            return np.ones((1,))
        order, method_order = self.n - 1, self.order
        parity = 0
        num_terms = order + method_order
        if self.method.startswith('central'):
            parity = (order % 2) + 1
            num_terms, order = num_terms // 2, order // 2
        fd_mat = self._fd_mat(step_ratio, parity, num_terms)
        fd_rule = linalg.pinv(fd_mat)[order]
        if self.method == 'backward' and self.n % 2 == 0:
            fd_rule *= -1
        return fd_rule

    def _eval_first_condition(self):
        even_order = (self.n % 2 == 0)
        return ((even_order and self.method == 'central') or
                self.method not in ['central', 'complex'])

    @staticmethod
    def _central_even(fun, f_x0i, x0i, h, *args, **kwds):
        return (fun(x0i + h, *args, **kwds) +
                fun(x0i - h, *args, **kwds)) / 2.0 - f_x0i

    @staticmethod
    def _central(fun, f_x0i, x0i, h, *args, **kwds):
        return (fun(x0i + h, *args, **kwds) -
                fun(x0i - h, *args, **kwds)) / 2.0

    @staticmethod
    def _forward(fun, f_x0i, x0i, h, *args, **kwds):
        return (fun(x0i + h, *args, **kwds) - f_x0i)

    @staticmethod
    def _backward(fun, f_x0i, x0i, h, *args, **kwds):
        return (f_x0i - fun(x0i - h))

    @staticmethod
    def _complex(f, fx, x, h, *args, **kwds):
        return f(x + 1j * h, *args, **kwds).imag

    @staticmethod
    def _complex2(f, fx, x, h, *args, **kwds):
        z = bicomplex(x + 1j * h, h)
        return f(z, *args, **kwds).imag12

    def _get_function_name(self):
        name = '_%s' % self.method
        if self.method == 'central' and (self.n % 2) == 0:
            name = name + '_even'
        elif self.method == 'complex' and self.n > 1:
            if self.n == 2:
                name = name + '2'
            else:
                raise ValueError('Complex method only support first and'
                                 'second order derivatives.')
        return name

    def _apply_fd_rule(self, fd_rule, sequence, steps):
        '''
        Return derivative estimates of f at x0 for a sequence of stepsizes h

        Member variables used
        ---------------------
        n
        '''
        f_del, h, original_shape = self._vstack(sequence, steps)

        ne = h.shape[0]
        if ne < fd_rule.size and self.method not in ['complex']:
            raise ValueError('num_steps (%d) must  be larger than '
                             'n + order - 1 (%d)'
                             ' (n=%d, order=%d)' % (ne, fd_rule.size,
                                                    self.n, self.order)
                             )
        f_diff = convolve1d(f_del, fd_rule[::-1], axis=0,
                            origin=(fd_rule.size-1)//2)

        der_init = f_diff / (h ** self.n)
        ne = max(ne - fd_rule.size + 1, 1)
        return der_init[:ne], h[:ne], original_shape


class Gradient(Derivative):
    __doc__ = _cmn_doc % dict(
        derivative='Gradient',
        scale_backward=str(default_scale('backward')),
        scale_central=str(default_scale('central')),
        scale_complex=str(default_scale('complex')),
        scale_forward=str(default_scale('forward')),
        extra_method="",
        returns="""
    Returns
    -------
    grad : array
        gradient
    """, extra_note="", example="""
    Examples
    --------
    >>> import numpy as np
    >>> import numdifftools.nd_cstep as ndc
    >>> fun = lambda x: np.sum(x**2)
    >>> dfun = ndc.Gradient(fun)
    >>> dfun([1,2,3])
    array([ 2.,  4.,  6.])

    # At [x,y] = [1,1], compute the numerical gradient
    # of the function sin(x-y) + y*exp(x)

    >>> sin = np.sin; exp = np.exp
    >>> z = lambda xy: sin(xy[0]-xy[1]) + xy[1]*exp(xy[0])
    >>> dz = ndc.Gradient(z)
    >>> grad2 = dz([1, 1])
    >>> grad2
    array([ 3.71828183,  1.71828183])

    # At the global minimizer (1,1) of the Rosenbrock function,
    # compute the gradient. It should be essentially zero.

    >>> rosen = lambda x : (1-x[0])**2 + 105.*(x[1]-x[0]**2)**2
    >>> rd = ndc.Gradient(rosen)
    >>> grad3 = rd([1,1])
    >>> np.allclose(grad3,[0, 0])
    True""", see_also="""
    See also
    --------
    Derivative, Hessian, Jacobian
    """)

    @staticmethod
    def _central(f, fx, x, h, *args, **kwds):
        n = len(x)
        increments = np.identity(n) * h
        partials = [(f(x + hi, *args, **kwds) - f(x - hi, *args, **kwds)) / 2.0
                    for hi in increments]
        return np.array(partials).T

    @staticmethod
    def _backward(f, fx, x, h, *args, **kwds):
        n = len(x)
        increments = np.identity(n) * h
        partials = [(fx - f(x - hi, *args, **kwds)) for hi in increments]
        return np.array(partials).T

    @staticmethod
    def _forward(f, fx, x, h, *args, **kwds):
        n = len(x)
        increments = np.identity(n) * h
        partials = [(f(x + hi, *args, **kwds) - fx) for hi in increments]
        return np.array(partials).T

    @staticmethod
    def _complex(f, fx, x, h, *args, **kwds):
        # From Guilherme P. de Freitas, numpy mailing list
        # http://mail.scipy.org/pipermail/numpy-discussion/2010-May/050250.html
        n = len(x)
        increments = np.identity(n) * 1j * h
        partials = [f(x + ih, *args, **kwds).imag for ih in increments]
        return np.array(partials).T


class Jacobian(Gradient):
    __doc__ = _cmn_doc % dict(
        derivative='Jacobian',
        scale_backward=str(default_scale('backward')),
        scale_central=str(default_scale('central')),
        scale_complex=str(default_scale('complex')),
        scale_forward=str(default_scale('forward')),
        extra_method="",
        returns="""
    Returns
    -------
    jacob : array
        Jacobian
    """, extra_note="""
    If f returns a 1d array, it returns a Jacobian. If a 2d array is returned
    by f (e.g., with a value for each observation), it returns a 3d array
    with the Jacobian of each observation with shape xk x nobs x xk. I.e.,
    the Jacobian of the first observation would be [:, 0, :]
    """, example='''
     Examples
    --------
    >>> import numdifftools.nd_cstep as ndc

    #(nonlinear least squares)

    >>> xdata = np.reshape(np.arange(0,1,0.1),(-1,1))
    >>> ydata = 1+2*np.exp(0.75*xdata)
    >>> fun = lambda c: (c[0]+c[1]*np.exp(c[2]*xdata) - ydata)**2

    >>> Jfun = ndc.Jacobian(fun)
    >>> Jfun([1,2,0.75]).reshape((3,-1)) # should be numerically zero
    array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])

    >>> fun2 = lambda x : x[0]*x[1]*x[2] + np.exp(x[0])*x[1]
    >>> Jfun3 = ndc.Jacobian(fun2)
    >>> Jfun3([3.,5.,7.])
    array([ 135.42768462,   41.08553692,   15.        ])
    ''', see_also="""
    See also
    --------
    Derivative, Hessian, Gradient
    """)


# class _Hessian(_Derivative):
#
#     @staticmethod
#     def default_scale(method, n=2):
#         return dict(central=8, central2=8, complex_=2.5,
#                     hybrid=6).get(method, 4)


class Hessian(_Derivative):
    def __init__(self, f, steps=None, method='central', full_output=False):
        super(Hessian, self).__init__(f, steps, method, full_output,
                                      n=2, order=1)

    __doc__ = _cmn_doc % dict(
        derivative='Hessian',
        scale_backward=str(default_scale('backward', 2)),
        scale_central=str(default_scale('central', 2)),
        scale_complex=str(default_scale('complex', 2)),
        scale_forward=str(default_scale('forward', 2)),
        extra_method="""'central2' : central difference derivative (scale=%s)
             hybrid : finite difference and complex-step (scale=%s)
        """ % (default_scale('central2', 2), default_scale('hybrid', 2)),
        returns="""
    Returns
    -------
    hess : ndarray
       array of partial second derivatives, Hessian
    """, extra_note="""Computes the Hessian according to method as:
    'forward', Eq. (7):
        1/(d_j*d_k) * ((f(x + d[j]*e[j] + d[k]*e[k]) - f(x + d[j]*e[j])))
    'central2', Eq. (8):
        1/(2*d_j*d_k) * ((f(x + d[j]*e[j] + d[k]*e[k]) - f(x + d[j]*e[j])) -
                         (f(x + d[k]*e[k]) - f(x)) +
                         (f(x - d[j]*e[j] - d[k]*e[k]) - f(x + d[j]*e[j])) -
                         (f(x - d[k]*e[k]) - f(x)))
    'central', Eq. (9):
        1/(4*d_j*d_k) * ((f(x + d[j]*e[j] + d[k]*e[k]) -
                          f(x + d[j]*e[j] - d[k]*e[k])) -
                         (f(x - d[j]*e[j] + d[k]*e[k]) -
                          f(x - d[j]*e[j] - d[k]*e[k]))
    'hybrid', Eq. (10):
        1/(2*d_j*d_k) * imag(f(x + i*d[j]*e[j] + d[k]*e[k]) -
                            f(x + i*d[j]*e[j] - d[k]*e[k]))
    where e[j] is a vector with element j == 1 and the rest are zero and
    d[i] is steps[i].
    """, example="""
    Examples
    --------
    >>> import numpy as np
    >>> import numdifftools.nd_cstep as ndc

    # Rosenbrock function, minimized at [1,1]

    >>> rosen = lambda x : (1.-x[0])**2 + 105*(x[1]-x[0]**2)**2
    >>> Hfun = ndc.Hessian(rosen)
    >>> h = Hfun([1, 1])
    >>> h
    array([[ 842., -420.],
           [-420.,  210.]])

    # cos(x-y), at (0,0)

    >>> cos = np.cos
    >>> fun = lambda xy : cos(xy[0]-xy[1])
    >>> Hfun2 = ndc.Hessian(fun)
    >>> h2 = Hfun2([0, 0])
    >>> h2
    array([[-1.,  1.],
           [ 1., -1.]])""", see_also="""
    See also
    --------
    Derivative, Hessian
    """)

    def _eval_first_condition(self):
        return self.method not in ['complex', 'central']

    @staticmethod
    def _hybrid(f, fx, x, h, *args, **kwargs):
        '''Calculate Hessian with hybrid finite difference and complex-step
        derivative approximation
        The stepsize is the same for the complex and the finite difference part
        '''
        # TODO: might want to consider lowering the step for pure derivatives
        n = len(x)
        # h = _default_base_step(x, 3, base_step, n)
        ee = np.diag(h)
        hess = 2. * np.outer(h, h)

        for i in range(n):
            for j in range(i, n):
                hess[i, j] = (f(x + 1j * ee[i, :] + ee[j, :], *args,
                                **kwargs) -
                              f(*((x + 1j * ee[i, :] - ee[j, :],) + args),
                                  **kwargs)).imag / hess[j, i]
                hess[j, i] = hess[i, j]
        return hess

    @staticmethod
    def _complex(f, fx, x, h, *args, **kwargs):
        '''Calculate Hessian with complex-step derivative approximation
        '''
        n = len(x)
        ee = np.diag(h)
        hess = np.outer(h, h)
        for i in range(n):
            for j in range(i, n):
                zph = bicomplex(x + 1j * ee[i, :], ee[j, :])
                hess[i, j] = (f(zph, *args, **kwargs)).imag12 / hess[j, i]
                hess[j, i] = hess[i, j]
        return hess

    @staticmethod
    def _central(f, fx, x, h, *args, **kwargs):
        '''Eq 9.'''
        n = len(x)
        # h = _default_base_step(x, 4, base_step, n)
        ee = np.diag(h)
        hess = np.outer(h, h)

        for i in range(n):
            for j in range(i, n):
                hess[i, j] = (f(x + ee[i, :] + ee[j, :], *args, **kwargs) -
                              f(x + ee[i, :] - ee[j, :], *args, **kwargs) -
                              f(x - ee[i, :] + ee[j, :], *args, **kwargs) +
                              f(x - ee[i, :] - ee[j, :], *args, **kwargs)
                              ) / (4. * hess[j, i])
                hess[j, i] = hess[i, j]
        return hess

    @staticmethod
    def _central2(f, fx, x, h, *args, **kwargs):
        '''Eq. 8'''
        n = len(x)
        # NOTE: ridout suggesting using eps**(1/4)*theta
        # h = _default_base_step(x, 3, base_step, n)
        ee = np.diag(h)
        f0 = fx  # f(x, *args, **kwargs)
        dtype = np.result_type(f0)
        g = np.empty(n, dtype=dtype)
        gg = np.empty(n, dtype=dtype)
        for i in range(n):
            g[i] = f(x + ee[i, :], *args, **kwargs)
            gg[i] = f(x - ee[i, :], *args, **kwargs)

        hess = np.empty((n, n), dtype=dtype)
        np.outer(h, h, out=hess)
        for i in range(n):
            for j in range(i, n):
                hess[i, j] = (f(x + ee[i, :] + ee[j, :], *args, **kwargs) -
                              g[i] - g[j] + f0 +
                              f(x - ee[i, :] - ee[j, :], *args, **kwargs) -
                              gg[i] - gg[j] + f0) / (2 * hess[j, i])
                hess[j, i] = hess[i, j]

        return hess

    @staticmethod
    def _forward(f, fx, x, h, *args, **kwargs):
        '''Eq. 7'''
        n = len(x)
        ee = np.diag(h)

        f0 = fx  # f(x, *args, **kwargs)
        dtype = np.result_type(f0)
        g = np.empty(n, dtype=dtype)
        for i in range(n):
            g[i] = f(x + ee[i, :], *args, **kwargs)

        hess = np.empty((n, n), dtype=dtype)
        np.outer(h, h, out=hess)
        for i in range(n):
            for j in range(i, n):
                hess[i, j] = (f(x + ee[i, :] + ee[j, :], *args, **kwargs) -
                              g[i] - g[j] + f0) / hess[j, i]
                hess[j, i] = hess[i, j]
        return hess

    def _backward(self, f, fx, x, h, *args, **kwargs):
        return self._forward(f, fx, x, -h, *args, **kwargs)


def main():
    import statsmodels.api as sm

    data = sm.datasets.spector.load()
    data.exog = sm.add_constant(data.exog, prepend=False)
    mod = sm.Probit(data.endog, data.exog)
    _res = mod.fit(method="newton")
    _test_params = [1, 0.25, 1.4, -7]
    _llf = mod.loglike
    _score = mod.score
    _hess = mod.hessian

    def fun(beta, x):
        return np.dot(x, beta).sum(0)

    def fun1(beta, y, x):
        # print(beta.shape, x.shape)
        xb = np.dot(x, beta)
        return (y - xb) ** 2  # (xb-xb.mean(0))**2

    def fun2(beta, y, x):
        # print(beta.shape, x.shape)
        return fun1(beta, y, x).sum(0)

    nobs = 200
    x = np.random.randn(nobs, 3)

    # xk = np.array([1, 2, 3])
    xk = np.array([1., 1., 1.])
    # xk = np.zeros(3)
    beta = xk
    y = np.dot(x, beta) + 0.1 * np.random.randn(nobs)
    xk = np.dot(np.linalg.pinv(x), y)

    epsilon = 1e-6
    args = (y, x)
    from scipy import optimize
    _xfmin = optimize.fmin(fun2, (0, 0, 0), args)  # @UndefinedVariable
    # print(approx_fprime((1, 2, 3), fun, steps, x))
    jac = Gradient(fun1, epsilon, method='forward')(xk, *args)
    jacmin = Gradient(fun1, -epsilon, method='forward')(xk, *args)
    # print(jac)
    print(jac.sum(0))
    print('\nnp.dot(jac.T, jac)')
    print(np.dot(jac.T, jac))
    print('\n2*np.dot(x.T, x)')
    print(2 * np.dot(x.T, x))
    jac2 = (jac + jacmin) / 2.
    print(np.dot(jac2.T, jac2))

    # he = approx_hess(xk,fun2,steps,*args)
    print(Hessian(fun2, 1e-3, method='central2')(xk, *args))
    he = Hessian(fun2, method='central2')(xk, *args)
    print('hessfd')
    print(he)
    print('base_step =', None)
    print(he - 2 * np.dot(x.T, x))

    for eps in [1e-3, 1e-4, 1e-5, 1e-6]:
        print('eps =', eps)
        print(Hessian(fun2, eps, method='central2')(xk, *args) -
              2 * np.dot(x.T, x))

    hcs2 = Hessian(fun2, method='hybrid')(xk, *args)
    print('hcs2')
    print(hcs2 - 2 * np.dot(x.T, x))

    hfd3 = Hessian(fun2, method='central')(xk, *args)
    print('hfd3')
    print(hfd3 - 2 * np.dot(x.T, x))

    hfi = []
    epsi = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]) * 10.
    for eps in epsi:
        h = eps * np.maximum(np.log1p(np.abs(xk)), 0.1)
        hfi.append(Hessian(fun2, h, method='hybrid')(xk, *args))
        print('hfi, eps =', eps)
        print(hfi[-1] - 2 * np.dot(x.T, x))

    import numdifftools as nd
    print('Dea3')
    err = 1000 * np.ones(hfi[0].shape)
    val = np.zeros(err.shape)
    errt = []
    for i in range(len(hfi) - 2):
        tval, terr = nd.dea3(hfi[i], hfi[i + 1], hfi[i + 2])
        errt.append(terr)
        k = np.flatnonzero(terr < err)
        if k.size > 0:
            np.put(val, k, tval.flat[k])
            np.put(err, k, terr.flat[k])
    print(val - 2 * np.dot(x.T, x))
    print(err)
    erri = [v.max() for v in errt]

    plt.loglog(epsi[1:-1], erri)
    plt.show('hold')
    hnd = nd.Hessian(lambda a: fun2(a, y, x))
    hessnd = hnd(xk)
    print('numdiff')
    print(hessnd - 2 * np.dot(x.T, x))
    # assert_almost_equal(hessnd, he[0])
    gnd = nd.Gradient(lambda a: fun2(a, y, x))
    _gradnd = gnd(xk)

    print(Derivative(np.cosh)(0))
    print(nd.Derivative(np.cosh)(0))


def _example3(x=0.0001, fun_name='inv', epsilon=None, method='central',
              scale=None, n=1):
    fun0, dfun = get_test_function(fun_name, n)
    if dfun is None:
        return
    fd = Derivative(fun0, steps=epsilon, method=method, n=n)
    t = []
    scales = np.arange(1, 13, 0.5)
    for scale in scales:
        fd.steps.scale = scale
        try:
            val = fd(x)
        except Exception:
            val = np.nan
        t.append(val)
    t = np.array(t)
    tt = dfun(x)
    relativ_error = np.abs(t - tt) / (np.abs(tt) + 1e-17) + 1e-17
    if not (relativ_error>0).all():
        return
    plt.semilogy(scales, relativ_error)
    plt.vlines(default_scale(fd.method, n), min(relativ_error), 1)
    plt.xlabel('scales')
    plt.ylabel('Relative error')
    txt = ['', "1'st", "2'nd", "3'rd", "4'th", "5'th", "6'th"]
    plt.title("The %s derivative of %s using %s" % (txt[n], fun_name, method))

    plt.axis([min(scales), max(scales), relativ_error.min(), 1])
    plt.figure()
    #plt.show('hold')


def _example2(x=0.0001, fun_name='inv', epsilon=None, method='central',
              scale=None, n=1):
    fun0, dfun = get_test_function(fun_name, n)

    fd = Derivative(fun0, steps=epsilon, method=method, n=n)
    t = []
    orders = n + (n % 2) + np.arange(0, 12, 2)

    for order in orders:
        fd.order = order
        fd.steps.num_steps = n + order - 1
        t.append(fd(x))
    t = np.array(t)
    tt = dfun(x)
    plt.semilogy(orders, np.abs(t - tt) / (np.abs(tt) + 1e-17) + 1e-17)

    plt.show('hold')


def _example(x=0.0001, fun_name='inv', epsilon=None, method='central',
             scale=None):
    '''
    '''
    fun0, dfun = get_test_function(fun_name)

    h = _default_base_step(x, scale=2, epsilon=None)  # 1e-4

    fd = Derivative(fun0, steps=epsilon, method=method, scale=scale,
                    full_output=True)

    t, res = fd(x)

    txt = (' (f(x+h)-f(x))/h = %g\n' %
           ((fun0(x + h) - fun0(x)) / h))
    deltas = np.array([h for h in epsilon(x, fd.scale)])

    print((txt +
           '      true df(x) = %20.15g\n' +
           ' estimated df(x) = %20.15g\n' +
           ' true err = %g\n err estimate = %g\n relative err = %g\n'
           ' delta = %g\n') % (dfun(x), t, dfun(x) - t,
                               res.error_estimate,
                               res.error_estimate / t,
                               deltas.flat[res.index]))
    # plt.show('hold')


def test_docstrings():
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)


if __name__ == '__main__':  # pragma : no cover
    # test_docstrings()
    # main()
    num_extrap = 5
    method = 'complex'
    for name in ['log10'] : #function_names:
        for n in range(1, 3):
            if method != 'complex':
                num_steps = n + 1 + num_extrap
                if method == 'central':
                    num_steps = (n+1) // 2 + num_extrap
            else:
                num_steps = 1 + num_extrap
            step_ratio = 4**(1./n)
            epsilon = StepsGenerator(num_steps=num_steps,
                                     step_ratio=step_ratio,
                                     offset=1, use_exact_steps=True)
            _example3(x=0.5, fun_name=name, epsilon=epsilon, method=method,
                      scale=None, n=n)
    plt.show('hold')
    #     import nxs
#     steps = StepsGenerator(num_steps=7)
#     d = Derivative(np.cos, method='central', steps=steps,
#                    full_output=True)
#     print(d([0, 1e5*np.pi*2]))
#     print(d(1e10*np.pi*2))

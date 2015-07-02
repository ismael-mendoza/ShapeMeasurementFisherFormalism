Introduction - Derivative Estimation
====================================

The general problem of differentiation of a function typically pops up
in three ways in Python.

-  The symbolic derivative of a function.

-  Compute numerical derivatives of a function defined only by a
   sequence of data points.

-  Compute numerical derivatives of a analytically supplied function.

Clearly the first member of this list is the domain of the symbolic
toolbox SymPy, or some set of symbolic tools. Numerical differentiation of a
function defined by data points can be achieved with the function
gradient, or perhaps by differentiation of a curve fit to the data,
perhaps to an interpolating spline or a least squares spline fit.

The third class of differentiation problems is where Numdifftools is valuable.
This document will describe the methods used in Numdifftools and in particular
the Derivative class.


Numerical differentiation of a general function of one variable
===============================================================

Surely you recall the traditional definition of a derivative, in terms
of a limit.

.. math::
    f'(x) = \lim_{\delta \to 0}{\frac{f(x+\delta) - f(x)}{\delta}}
    :label: 1

For small :math:`\delta`, the limit approaches :math:`f'(x)`. This is a
one-sided approximation for the derivative. For a fixed value of
:math:`\delta`, this is also known as a finite difference approximation
(a forward difference.) Other approximations for the derivative are also
available. We will see the origin of these approximations in the Taylor
series expansion of a function :math:`f(x)` around some point
:math:`x_0`.

.. math::
    f(x) &= f(x_0) + (x - x_0)f'(x_0) + \frac{(x - x_0)^2}{2} f''(x_0) +
    :label: 2

    & \frac{(x - x_0)^3}{6} f^{(3)}(x_0) + \frac{(x - x_0)^4}{24} f^{(4)}(x_0) + \\


    & \frac{(x - x_0)^5}{120} f^{(5)}(x_0) + \frac{(x - x_0)^6}{720} f^{(6)}(x_0) +...\\

Truncate the series in :eq:`2` to the first three terms, then form the
forward difference approximation :eq:`1`, where :math:`x = x_0 + \delta`.

.. math::
    f'(x_0) = \frac{f(x_0+\delta) - f(x_0)}{\delta} - \frac{\delta}{2} f''(x_0) - \frac{\delta^2}{6} f'''(x_0) + ...
    :label: 3

When :math:`\delta` is small, :math:`\delta^2` and any higher powers are
vanishingly small. So we tend to ignore those higher powers, and
describe the approximation in :eq:`3` as a first order approximation since
the error in this approximation approaches zero at the same rate as the
first power of :math:`\delta`.  [1]_ The values of :math:`f''(x_0)` and
:math:`f'''(x_0)`, while unknown to us, are fixed constants as
:math:`\delta` varies.

Higher order approximations arise in the same fashion. The central
difference :eq:`4` is a second order approximation.

.. math::
    f'(x_0) = \frac{f(x_0+\delta) - f(x_0-\delta)}{2\delta} - \frac{\delta^2}{3} f'''(x_0) + ...
    :label: 4


Unequally spaced finite difference rules
========================================

While most finite difference rules used to differentiate a function will
use equally spaced points, this fails to be appropriate when one does
not know the final spacing. Adaptive quadrature rules can succeed by
subdividing each sub-interval as necessary. But an adaptive
differentiation scheme must work differently, since differentiation is a
point estimate. Derivative generates a sequence of sample points that follow a log
spacing away from the point in question, then it uses a single rule
(generated on the fly) to estimate the desired derivative. Because the
points are log spaced, the same rule applies at any scale, with only a
scale factor applied.


Odd and even transformations of a function
==========================================

Returning to the Taylor series expansion of :math:`f(x)` around some
point :math:`x_0`, an even function  [2]_ around :math:`x_0` must have
all the odd order derivatives vanish at :math:`x_0`. An odd function has
all its even derivatives vanish from its expansion. Consider the derived
functions :math:`f_{odd}(x)` and :math:`f_{even}(x)`.

.. math::
    f_{odd}(x) = \frac{f(x - x_0) - f(-x - x_0)}{2}
    :label: 5

The Taylor series expansion of :math:`f_{odd}(x)` has the useful
property that we have killed off any even order terms, but the odd order
terms are identical to :math:`f(x)`, as expanded around :math:`x_0`.

.. math::
    f_{odd}(x) &= (x - x_0)f'(x_0) + \frac{(x - x_0)^3}{6} f^{(3)}(x_0) +
    :label: 6

    & \frac{(x - x_0)^5}{120} f^{(5)}(x_0) + \frac{(x - x_0)^7}{5040} f^{(7)}(x_0) +...\\


Likewise, :math:`f_{even}(x)` has no odd order terms or a constant term,
but other even order terms that are identical to :math:`f(x)`.

.. math::
    f_{even}(x) = \frac{f(-x-x_0) - 2f(x_0) + f(x-x_0)}{2}
    :label: 7

.. math::
    f_{even}(x) &= \frac{(x - x_0)^2}{2} f^{(2)}(x_0) + \frac{(x - x_0)^4}{24} f^{(4)}(x_0) +
    :label: 8

    & \frac{(x - x_0)^6}{720} f^{(6)}(x_0) + \frac{(x - x_0)^8}{40320} f^{(8)}(x_0) + ...\\


The point of these transformations is we can rather simply generate a
higher order approximation for any odd order derivatives of :math:`f(x)`
by working with :math:`f_{odd}(x)`. Even order derivatives of
:math:`f(x)` are similarly generated from :math:`f_{even}(x)`. For
example, a second order approximation for :math:`f'(x_0)` is trivially
written in :eq:`9` as a function of :math:`\delta`.

.. math::
    f'(x_0; \delta) = \frac{f_{odd}(x_0 + \delta)}{\delta} - \frac{\delta^2}{6} f^{(3)}(x_0)
    :label: 9

We can do better rather simply, so why not? :eq:`10` shows a fourth order
approximation for :math:`f'(x_0)`.

.. math::
    f'(x_0; \delta) = \frac{8 f_{odd}(x_0+\delta)-f_{odd}(x_0+2\delta)}{6\delta} + \frac{\delta^4}{30} f^{(5)}(x_0)
    :label: 10

Again, the next non-zero term :eq:`11` in that expansion has a higher power
of :math:`\delta` on it, so we would normally ignore it since the lowest
order neglected term should dominate the behavior for small
:math:`\delta`.

.. math::
    \frac{\delta^6}{252} f^{(7)}(x_0)
    :label: 11

Derivative uses similar approximations for all derivatives of :math:`f` up to the
fourth order. Of course, its not always possible for evaluation of a
function on both sides of a point, as central difference rules will
require. In these cases, you can specify forward or backward difference
rules as appropriate.


Romberg extrapolation methodology applied to derivative estimation
==================================================================

Some individuals might suggest that the above set of approximations are
entirely adequate for any sane person. Can we do better?

Suppose we were to generate several different estimates of the
approximation in :eq:`3` for different values of :math:`\delta` at a fixed
:math:`x_0`. Thus, choose a single :math:`\delta`, estimate a
corresponding resulting approximation to :math:`f'(x_0)`, then do the
same for :math:`\delta/2`. If we assume that the error drops off
linearly as :math:`\delta \to 0`, then it is a simple matter to
extrapolate this process to a zero step size. Our lack of knowledge of
:math:`f''(x_0)` is irrelevant. All that matters is :math:`\delta` is
small enough that the linear term dominates so we can ignore the
quadratic term, therefore the error is purely linear.

.. math::
    f'(x_0) = \frac{f(x_0+\delta) - f(x_0)}{\delta} - \frac{\delta}{2} f''(x_0)
    :label: 12

The linear extrapolant for this interval halving scheme as
:math:`\delta \to 0` is given by:

.. math::
    f^{'}_{0} = 2 f^{'}_{\delta/2} - f^{'}_{\delta}
    :label: 13

Since I've always been a big fan of convincing myself that something
will work before I proceed too far, lets try this out in Python.
Consider the function :math:`e^x`. Generate a pair of approximations to
:math:`f'(0)`, once at :math:`\delta` of 0.1, and the second
approximation at :math:`1/2` that value. Recall that
:math:`\frac{d(e^x)}{dx} = e^x`, so at x = 0, the derivative should be
exactly 1. How well will we do?

   >>> from numpy import exp, allclose
   >>> f = exp
   >>> dx = 0.1
   >>> df1 = (f(dx) - f(0))/dx
   >>> allclose(df1, 1.05170918075648)
   True

   >>> df2 = (f(dx/2) - f(0))/(dx/2)
   >>> allclose(df2, 1.02542192752048)
   True

   >>> allclose(2*df2 - df1, 0.999134674284488)
   True


In fact, this worked very nicely, reducing the error to roughly 1
percent of our initial estimates. Should we be surprised at this
reduction? Not if we recall that last term in :eq:`3`. We saw there that the
next term in the expansion was :math:`O(\delta^2)`. Since :math:`\delta`
was 0.1 in our experiment, that 1 percent number makes perfect sense.

The Romberg extrapolant in :eq:`13` assumed a linear process, with a
specific reduction in :math:`\delta` by a factor of 2. Assume the two
term (linear + quadratic) residual term in :eq:`3`, evaluating our
approximation there with a third value of :math:`\delta`. Again, assume
the step size is cut in half again. The three term Romberg extrapolant
is given by:

.. math::
    f'_0 = \frac{1}{3}f'_\delta - 2f'_{\delta/2} + \frac{8}{3}f'_{\delta/4}
    :label: 14

A quick test in Python yields much better results yet.

    >>> from numpy import exp, allclose
    >>> f = exp
    >>> dx = 0.1

    >>> df1 = (f(dx) - f(0))/dx
    >>> allclose(df1,  1.05170918075648)
    True

    >>> df2 = (f(dx/2) - f(0))/(dx/2)
    >>> allclose(df2, 1.02542192752048)
    True

    >>> df3 = (f(dx/4) - f(0))/(dx/4)
    >>> allclose(df3, 1.01260482097715)
    True

    >>> allclose(1./3*df1 - 2*df2 + 8./3*df3, 1.00000539448361)
    True

Again, Derivative uses the appropriate multiple term Romberg extrapolants for all
derivatives of :math:`f` up to the fourth order. This, combined with the
use of high order approximations for the derivatives, allows the use of
quite large step sizes. See [LynessMoler1966]_ and [LynessMoler1969]_.


Uncertainty estimates for Derivative
====================================

We can view the Romberg extrapolation step as a polynomial curve fit in
the step size parameter :math:`\delta`. Our desired extrapolated value
is seen as simply the constant term coefficient in that polynomial
model. Remember though, this polynomial model (see :eq:`10` and :eq:`11`) has
only a few terms in it with known non-zero coefficients. That is, we
will expect a constant term :math:`a_0`, a term of the form
:math:`a_1 \delta^4`, and a third term :math:`a_2 \delta^6`.

A neat trick to compute the statistical uncertainty in the estimate of
our desired derivative is to use statistical methodology for that error
estimate. While I do appreciate that there is nothing truly statistical
or stochastic in this estimate, the approach still works nicely,
providing a very reasonable estimate in practice. A three term
Romberg-like extrapolant, then evaluated at four distinct values for
:math:`\delta`, will yield an estimate of the standard error of the
constant term, with one spare degree of freedom. The uncertainty is then
derived by multiplying that standard error by the appropriate percentile
from the Students-t distribution.

   >>> import scipy.stats as ss
   >>> allclose(ss.t.cdf(12.7062047361747, 1), 0.975)
   True

This critical level will yield a two-sided confidence interval of 95
percent.

These error estimates are also of value in a difference sense. Since
they are efficiently generated at all the different scales, the
particular spacing which yields the minimum predicted error is chosen as
the best derivative estimate. This has been shown to work consistently
well. A spacing too large tends to have large errors of approximation
due to the finite difference schemes used. But a too small spacing is
bad also, in that we see a significant amplification of least
significant fit errors in the approximation. A middle value generally
seems to yield quite good results. For example, Derivative will estimate the
derivative of :math:`e^x` automatically. As we see, the final overall
spacing used was 0.02166085.

    >>> import numdifftools as nd
    >>> from numpy import exp, allclose
    >>> f = nd.Derivative(exp)
    >>> allclose(f(1), 2.71828183)
    True
    >>> allclose(f.error_estimate, 5.21804822e-14)
    True
    >>> allclose(f.final_delta, 0.02166085)
    True


However, if we force the step size to be artificially large, then
approximation error takes over.

    >>> f = nd.Derivative(exp, delta=1, step_nom=1)
    >>> allclose(f(1), 3.19452805)
    True
    >>> allclose(f(1)-exp(1), 0.47624622)
    True
    >>> f.final_delta
    array([ 1.])

And if the step size is forced to be too small, then we see noise
dominate the problem.

   >>> f = nd.Derivative(exp, delta=.0000000001, step_nom=1)
   >>> allclose(f(1), 2.71828093)
   True
   >>> allclose(f(1) - exp(1), -8.97648138e-07)
   True
   >>> allclose(f.final_delta, 1.00000008e-10)
   True


Numdifftools, like Goldilocks in the fairy tale bearing her name, stays comfortably
in the middle ground.

Derivative in action
====================

How does numdifftools.Derivative work in action? A simple nonlinear
function with a well known derivative is :math:`e^x`. At :math:`x = 0`,
the derivative should be 1.

   >>> f = nd.Derivative(exp)
   >>> f(0)
   array([ 1.])

   >>> allclose(f.error_estimate, 5.28466160e-14)
   True

A second simple example comes from trig functions. The first four
derivatives of the sine function, evaluated at :math:`x = 0`, should be
respectively :math:`[cos(0), -sin(0), -cos(0), sin(0)]`, or
:math:`[1,0,-1,0]`.

    >>> from numpy import sin, allclose
    >>> import numdifftools as nd
    >>> df = nd.Derivative(sin, 1)
    >>> allclose(df(0), 1.)
    True

    >>> ddf = nd.Derivative(sin, 2)
    >>> allclose(ddf(0), 0.)
    True

    >>> dddf = nd.Derivative(sin, 3)
    >>> allclose(dddf(0), -1.)
    True

    >>> ddddf = nd.Derivative(sin, 4)
    >>> allclose(ddddf(0), 0.)
    True


Gradient and Hessian  estimation
================================

Estimation of the gradient vector (numdifftools.Gradient) of a function of multiple variables
is a simple task, requiring merely repeated calls to numdifftools.Derivative. Likewise, the
diagonal elements of the hessian matrix are merely pure second partial
derivatives of a function. numdifftools.Hessdiag accomplishes this task, again calling numdifftools.Derivative multiple
times. Efficient computation of the off-diagonal (mixed partial derivative) elements of the
Hessian matrix uses a scheme much like that of numdifftools.Derivative, wherein
numdifftools.Derivative is called to determine an initial step size, then Romberg extrapolation
is used to improve a set of second order finite difference estimates of those mixed partials.

Conclusion
==========

numdifftools.Derivative is an a adaptive scheme that can compute the derivative of
arbitrary (well behaved) functions. It is reasonably fast as an adaptive method.
Many options have been provided for the user who wishes the ultimate amount of
control over the estimation.



Acknowledgments
===============


My thanks are due to Shaun Simmons for convincing me to learn enough
LaTeX to write this document.

References
==========
.. [LynessMoler1966] Lyness, J. M., Moler, C. B. (1966). Vandermonde Systems and Numerical
                     Differentiation. *Numerische Mathematik*.

.. [LynessMoler1969] Lyness, J. M., Moler, C. B. (1969). Generalized Romberg Methods for
                     Integrals of Derivatives. *Numerische Mathematik*.

.. [NAG] *NAG Library*. NAG Fortran Library Document: D04AAF


.. rubric:: Footnotes

.. [1]
   We would normally write these additional terms using O() notation,
   where all that matters is that the error term is :math:`O(\delta)` or
   perhaps :math:`O(\delta^2)`, but explicit understanding of these
   error terms will be useful in the Romberg extrapolation step later
   on.

.. [2]
   An even function is one which expresses an even symmetry around a
   given point. An even symmetry has the property that
   :math:`f(x) = f(-x)`. Likewise, an odd function expresses an odd
   symmetry, wherein :math:`f(x) = -f(-x)`.

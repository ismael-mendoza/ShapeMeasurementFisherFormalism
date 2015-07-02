""" Test functions for numdifftools module

"""
import unittest
import numdifftools.nd_cstep as nd
import numpy as np
from numpy.testing import assert_array_almost_equal


class TestStepGenerator(unittest.TestCase):

    def test_default_generator(self):
        step_gen = nd.StepsGenerator(base_step=None, num_steps=10,
                                     step_ratio=4, offset=-1)
        h = np.array([h for h in step_gen(0)])
        desired = np.array([3.088162e-04, 7.720404e-05,
                            1.930101e-05, 4.825253e-06, 1.206313e-06,
                            3.015783e-07, 7.539457e-08,
                            1.884864e-08, 4.712161e-09, 1.178040e-09])
        assert_array_almost_equal((h - desired) / desired, 0)

    def test_default_base_step(self):
        step_gen = nd.StepsGenerator(num_steps=1, offset=0)
        h = [h for h in step_gen(0)]
        desired = (10 * nd.EPS) ** (1. / 2) * 0.1
        assert_array_almost_equal((h[0] - desired) / desired, 0)

    def test_fixed_base_step(self):
        desired = 0.1
        step_gen = nd.StepsGenerator(base_step=desired, num_steps=1, scale=2,
                                     offset=0)
        h = [h for h in step_gen(0)]
        assert_array_almost_equal((h[0] - desired) / desired, 0)


class TestFornbergWeights(unittest.TestCase):
    def test_weights(self):
        x = np.r_[-1, 0, 1]
        xbar = 0
        k = 1
        weights = nd.fornberg_weights(x, xbar, k)
        np.testing.assert_allclose(weights, [-.5, 0, .5])


class TestDerivative(unittest.TestCase):
    def test_high_order_derivative_cos(self):
        true_vals = (-1.0, 0.0, 1.0, 0.0, -1.0, 0.0)
        for method in ['complex', 'central', 'forward', 'backward']:
            n_max = dict(complex=2, central=6).get(method, 5)
            for n in range(1, n_max+1):
                true_val = true_vals[n-1]
                for order in range(2, 9, 2):
                    d3cos = nd.Derivative(np.cos, n=n, order=order,
                                          method=method, full_output=True)
                    y, info = d3cos(np.pi / 2.0)
                    error = np.abs(y - true_val)
                    small = error < max(info.error_estimate*10, 10**n*1e-13)
                    self.assertTrue(small)

    def test_derivative_of_cos_x(self):
        x = np.r_[0, np.pi / 6.0, np.pi / 2.0]
        true_vals = (-np.sin(x), -np.cos(x), np.sin(x), np.cos(x), -np.sin(x),
                     -np.cos(x))
        for method in ['complex', 'central', 'forward', 'backward']:
            n_max = dict(complex=2, central=6).get(method, 5)
            for n in range(1, n_max+1):
                true_val = true_vals[n-1]
                start, stop, step = dict(central=(2, 7, 2)).get(method,
                                                                (1, 5, 1))
                for order in range(start, stop, step):
                    d3cos = nd.Derivative(np.cos, n=n, order=order,
                                          method=method, full_output=True)
                    y, info = d3cos(x)
                    error = np.abs(y - true_val)
                    small = error < np.maximum(info.error_estimate*15,
                                               10**n*1e-12)
                    self.assertTrue(small.all())

    def test_default_scale(self):
        for method, scale in zip(['complex', 'central', 'forward', 'backward'],
                                 [1, 3, 2, 2]):
            np.testing.assert_allclose(scale,
                                       nd.default_scale(method, n=1))

    def test_derivative_cube(self):
        '''Test for Issue 7'''
        def cube(x):
            return x * x * x
        dcube = nd.Derivative(cube)
        shape = (3, 2)
        x = np.ones(shape) * 2
        dx = dcube(x)
        assert_array_almost_equal(list(dx.shape), list(shape),
                                  decimal=8,
                                  err_msg='Shape mismatch')
        txt = 'First differing element %d\n value = %g,\n true value = %g'
        for i, (val, tval) in enumerate(zip(dx.ravel(), (3 * x**2).ravel())):
            assert_array_almost_equal(val, tval, decimal=8,
                                      err_msg=txt % (i, val, tval))

    def test_derivative_exp(self):
        # derivative of exp(x), at x == 0
        dexp = nd.Derivative(np.exp)
        assert_array_almost_equal(dexp(0), np.exp(0), decimal=8)

    def test_derivative_sin(self):
        # Evaluate the indicated (default = first)
        # derivative at multiple points
        dsin = nd.Derivative(np.sin)
        x = np.linspace(0, 2. * np.pi, 13)
        y = dsin(x)
        np.testing.assert_almost_equal(y, np.cos(x), decimal=8)

    def test_backward_derivative_on_sinh(self):
        # Compute the derivative of a function using a backward difference
        # scheme.  A backward scheme will only look below x0.
        dsinh = nd.Derivative(np.sinh, method='backward')
        self.assertAlmostEqual(dsinh(0.0), np.cosh(0.0))

    def test_central_and_forward_derivative_on_log(self):
        # Although a central rule may put some samples in the wrong places, it
        # may still succeed
        epsilon = nd.StepsGenerator(num_steps=10, offset=-1, step_ratio=2)
        dlog = nd.Derivative(np.log, method='central', steps=epsilon)
        x = 0.001
        self.assertAlmostEqual(dlog(x), 1.0 / x)

        # But forcing the use of a one-sided rule may be smart anyway
        dlog = nd.Derivative(np.log, method='forward', steps=epsilon)
        self.assertAlmostEqual(dlog(x), 1 / x)


class TestJacobian(unittest.TestCase):

    def testjacobian(self):
        xdata = np.reshape(np.arange(0, 1, 0.1), (-1, 1))
        ydata = 1 + 2 * np.exp(0.75 * xdata)

        def fun(c):
            return (c[0] + c[1] * np.exp(c[2] * xdata) - ydata) ** 2

        for method in ['complex', 'central', 'forward', 'backward']:
            for order in [2, 4]:
                Jfun = nd.Jacobian(fun, method=method, order=order)
                J = Jfun([1, 2, 0.75])  # should be numerically zero
                assert_array_almost_equal(J, np.zeros(J.shape))


class TestGradient(unittest.TestCase):

    def testgradient(self):
        def fun(x):
            return np.sum(x ** 2)
        dtrue = [2., 4., 6.]

        for method in ['complex', 'central', 'backward', 'forward']:
            for order in [2, 4]:
                dfun = nd.Gradient(fun, method=method, order=order)
                d = dfun([1, 2, 3])
                assert_array_almost_equal(d, dtrue)


class TestHessian(unittest.TestCase):

    def test_hessian_cosIx_yI_at_I0_0I(self):
        # cos(x-y), at (0,0)
        epsilon = nd.StepsGenerator(num_steps=10)
        cos = np.cos

        def fun(xy):
            return cos(xy[0] - xy[1])
        htrue = [[-1., 1.], [1., -1.]]
        methods = ['hybrid', 'complex', 'central', 'central2', 'forward',
                   'backward']
        for method in methods:
            Hfun2 = nd.Hessian(fun, method=method, steps=epsilon)
            h2 = Hfun2([0, 0])
            assert_array_almost_equal(h2, htrue)


if __name__ == '__main__':
    unittest.main()

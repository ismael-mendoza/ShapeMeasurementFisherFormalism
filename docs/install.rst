Installation
============

Install with GIT
----------------

The code is hosted on `github <https://github.com/ismael2395/ShapeMeasurementFisherFormalism>`_ so the easiest method to perform an initial installation is with `git <http://git-scm.com>`_::

	git clone https://github.com/ismael2395/ShapeMeasurementFisherFormalism.git

This will create a new subdirectory `ShapeMeasurementFisherFormalism` containing the latest stable version of the complete package.

Experts who already have a `correctly configured github account <https://help.github.com/articles/which-remote-url-should-i-use/#cloning-with-ssh>`_ might prefer this alternative::

	git clone git@github.com:ismael2395/ShapeMeasurementFisherFormalism.git

Update with GIT
---------------

You can update your local copy of the package at any time using::

	cd ShapeMeasurementFisherFormalism
	git update

Getting Started
---------------

Programs can be run directly from the top-level directory as long as you have the required packages already installed, e.g.::

	cd ShapeMeasurementFisherFormalism
	python generate.py --help

For an introduction to the available programs, see :doc:`here </programs>` and for examples of running these programs see :doc:`here </examples>`.

Required Packages
-----------------

The following python packages are required by this package:

* `numpy <http://www.numpy.org>`_ (version >= 1.9)
* `galsim <https://github.com/GalSim-developers/GalSim>`_ (version >= 1.2)
* `lmfit <http://cars9.uchicago.edu/software/python/lmfit/>`_ (version >= 0.8.3)

Note that `numpy <http://www.numpy.org>`_ is available in recent `anaconda <https://store.continuum.io/cshop/anaconda/>`_ or `enthought canopy <https://www.enthought.com/products/canopy/>`_ distributions. Installing GalSim is a more involved process, but well worth the effort. The `lmfit` package is only required if you will be running your own fits using the :ref:`prog-fitting`.
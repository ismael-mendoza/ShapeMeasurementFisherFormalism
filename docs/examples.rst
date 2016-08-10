Examples
========

Command-line Options
--------------------

Print out usage info for command-line options::

	./generate.py --help
	./fitting.py --help

Quick Simulation Demo
---------------------

The first thing to do is to generate a desired galaxy model (with an optional psf) using the generate.py file:: 

	python generate.py -p project -gal 1 --galaxy-model gaussian --psf_model psf_gaussian  --g1 0.2 --g2 0.2 --y0 0. --x0 0. --flux 1. --psf_flux 1. --hlr 0.5 --psf_fwhm 0.7 --snr 20.0

Display partial derivatives and save it to the project folder just created with generate.py::

	python display.py -p project --partials --snr 20. 

Display a fisher matrix elements without showing it when the command is executed::

	python display.py -p project --fisher --snr 20. --hide

It is important to always specificy the project that is being run as well as the signal to noise ratio to be used. . 
Finally one can also save all the possible outputs to the current project folder with::

	python display.py -p project --all --snr 20.

This command hides the output by default and saves all files in a pdf format. 

More complicated examples
------------------

Please refer to the :doc:`tutorial notebooks </notebooks>`.


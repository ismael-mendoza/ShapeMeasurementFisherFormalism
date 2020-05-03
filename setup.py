#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

requirements = [
    'numpy',
    'galsim',
    'lmfit',
]

setup(
    name='smff',
    version='1.0',
    description='Fisher Formalism for Weak Lensing using 1 or 2 Galsim galaxies.',
    long_description='Fisher Formalism for Weak Lensing using 1 or 2 Galsim galaxies.',
    author='Ismael Mendoza',
    author_email='imendoza@umich.edu',
    url='https://github.com/ismael2395/ShapeMeasurementFisherFormalism/',
    packages=[
        'smff',
    ],
    package_dir={'smff': 'smff'},
    scripts=[],
    install_requires=requirements,  # requirements,
    license='MIT',
)

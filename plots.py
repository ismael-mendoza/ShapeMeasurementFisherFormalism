#!/usr/bin/env python
"""Interface for generating plots on the spot.
"""

import argparse

import fisher

import draw

import defaults

import galfun

import info


def main():
    parser = argparse.ArgumentParser(description=(
        'Display different plots from the data produced by the package'), formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--project', default=defaults.PROJECT,
                        type=str,
                        help=('Specify a directory name where the project will'
                              'be saved. In this fashion each individual'
                              'project represents one analysis.'))

    parser.add_argument('--hide', action='store_true',
                        help='Do not show plots produced')

    parser.add_argument('--error_bars', action='store_true',
                        help='Show error_bars on bplots.')

    parser.add_argument('--bias_sigma', action='store_true',
                        help='Show bias over sigma instead of bias in bplots')

    parser.add_argument('-d', '--draw-galaxy', action='store_true',
                        help='Show original galaxies')

    parser.add_argument('-fp', '--partials', action='store_true',
                        help='Show partial derivative images.')

    parser.add_argument('-sp', '--second-partials', action='store_true',
                        help='Show second partial derivative images.')

    parser.add_argument('-f', '--fisher', action='store_true',
                        help='Show fisher matrix images.')

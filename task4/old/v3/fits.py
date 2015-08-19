import argparse
import defaults

def main():

    #have to add a require argument: filename to write and extract galaxies from.

    parser = argparse.ArgumentParser(description = 'Displays a triangle plot comparing the expected noise bias calculation with the fisher matrix formalism for a given file of galaxies.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--wdir', default = 'output', metavar = 'DIRECTORY', type = str,
    help = 'Specify a directory name where the output will be inputted and files read from.')
    parser.add_argument('-f', '--filename', default = 'galaxies', metavar = 'FILENAME', type = str,
    help = 'Specify a file where galaxies will be registered or read from.')
    parser.add_argument('-rltsdir', default = 'results', metavar = 'RESULTS', type = str,
    help = 'Specify a directory where output data from fits will be produced')
    parser.add_argument('-number-fits', metavar = 'N', required = True, type = int,
    help = 'Specify how many fits are run.')
    parser.add_argument('--run-fits',
    help = 'Run N different instantiations of noise of the added galaxies to produce a triangle plot. Fits all galaxies in given file N times. (for loop??? only if useful)')
    parser.add_argument('--run-fits-slac', metavar = 'SLAC_COMPUTER',
    help = 'Same as above but have to be in SLAC computer.')


    args = parser.parse_args()


if __name__ == '__main__':
    main()
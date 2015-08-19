import argparse
import defaults

def main():

    #have to add a require argument: filename to write and extract galaxies from.

    parser = argparse.ArgumentParser(description = 'Displays a triangle plot comparing the expected noise bias calculation with the fisher matrix formalism for a given file of galaxies.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)






    fits_group = parser.add_argument_group('Fitting options')
    fits_group.add_argument('--run-fits', metavar = 'N',
        help = 'Run N different instantiations of noise of the added galaxies to produce a triangle plot.')

    args = parser.parse_args()


if __name__ == '__main__':
    main()
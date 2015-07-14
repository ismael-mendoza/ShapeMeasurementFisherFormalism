import argparse
import defaults
import os

def csvIsEmpty(filename):
    #check if header is there so we do not write it again. File is not there if the file is empty. 
    with open(filename, 'r') as f:
        if(f.read() == ''):
            return True
        else:
            return False

def main():

    #initialize names used. 
    names = defaults.names()
    #initialize dictionary with default parameters
    dflt_params = defaults.parameters().dflt_params

    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(description = 'Generate galaxies specified by the user that are added to a file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--defaults', action = 'store_true',
        help = 'Print defaults use for generating the different types of galaxies.')
    parser.add_argument('--no-display', action = 'store_true',
        help = 'Do not display the image on screen (but will be stored as pdf files).')


    #add require argument: working directory
    parser.add_argument('--wd', default = 'output', metavar = 'DIRECTORY', type = str,
        help = 'Specify a directory name where the output will be inputted and files read from.')
    #have to add a require argument: filename to write and extract galaxies from.
    parser.add_argument('-f', '--filename', default = 'galaxies', metavar = 'FILENAME', type = str,
        help = 'Specify a file where galaxies will be registered or read from.')

    
    
    parser.add_argument('--id', required = True, default = None, metavar = 'GALAXYID', type = int,
        help = 'Add galaxy with id GALAXYID, if GALAXYID already exists in the table, you will change the parameters of the galaxy with that GALAXYID but not create a new entry.')
    parser.add_argument('--type', default = None, metavar = 'TYPE', type = str, choices = names.galaxy_choices,
        help = 'Change the galaxy type to TYPE.')

    #add all parameter argument to the parser. this can change 
    #what about if the galaxy is convolved psf? should this go specified in the same row of the table (yes, for now psf gaussian so two parameters only can expand later.)
    for name in names.parameters:
        parser.add_argument('--' + name, default = dflt_params[name], metavar = name.upper(), type = float, help = 'Change the parameter ' + name + ' from its default value.') 
    #for now only default parameters, no sigma or FWHM, will add sersic and exponential parameters later easily.

    #works if adding multiple galaxies?? 


    args = parser.parse_args()

    #check if directory exits. 
    if not os.path.isdir(args.wd):
        os.mkdir(args.wd)


    ##############better way to update???? 
    #check if filename exists, not necessary just open it and append things 
    with open(os.path.join(args.wd, args.filename + '.csv'), 'a') as csvfile:
        writer = csvfile.DictWriter(csvfile, fieldnames=names.header)
        reader = csvfile.DictReader(csvfile)
        row_to_write = dict(
            id = args.id, 
            x= args.x0,
            y0 = args.y0,
            flux = args.y0,
            hlr = args.y0,
            e1 = args.y0,
            e2 = args.y0,
            psf_flux = args.y0,
            psf_fwhm = args.y0
        )

    ##############better way to update???? 

        #if file is empty just write the new row
        if(csvIsEmpty): 
            writer.writeheader()
            writer.writerow(row_to_write)


        #else, check if the row is their and overwrite that row.
        else:
            for row in reader:
                if(row['id'] == args.id):




        #add header if the csv file is empty







    if args.defaults:
        #print defaults from default.py in some format. 
        pass




if __name__ == '__main__':
    main()
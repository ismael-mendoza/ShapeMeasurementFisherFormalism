import argparse
import defaults
import os
import csv
import shutil

def csvIsEmpty(filename):
    """checks each row and if any is not empty, then the file is not empty"""
    with open(filename, 'r') as f:
        for row in f:
            if len(row) != 0:
                return False
        else: return True


def main():

    #initialize names used. 
    names = defaults.names()
    #initialize dictionary with default parameters
    dflt_params = defaults.parameters().dflt_params

    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(description = 'Generate galaxies specified by the user that are added to a file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--defaults', action = 'store_true',
        help = 'Print defaults use for generating the different types of galaxies.')


    #add require argument: working directory
    parser.add_argument('--wd', default = 'output', metavar = 'DIRECTORY', type = str,
        help = 'Specify a directory name where the output will be inputted and files read from.')
    #have to add a require argument: filename to write and extract galaxies from.
    parser.add_argument('-f', '--filename', default = 'galaxies', metavar = 'FILENAME', type = str,
        help = 'Specify a file where galaxies will be registered or read from.')

    parser.add_argument('--id', required = True, default = None, metavar = 'GALAXYID', type = int,
        help = 'Add galaxy with id GALAXYID, if GALAXYID already exists in the table, you will change the parameters of the galaxy with that GALAXYID but not create a new entry.')
    parser.add_argument('--model', default = names.galaxy_models[0], metavar = 'MODEL', type = str, choices = names.galaxy_models,
        help = 'Change the galaxy type to MODEL.')

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

    #initialize name of files.
    filename = os.path.join(args.wd, args.filename + '.csv')
    if not os.path.isfile(filename):
        f = open(filename, 'w+')
        f.close()
    tempname = os.path.join(args.wd, 'temp' + '.csv')
    shutil.copyfile(filename, tempname) #creat a copy to read from and compare.
 
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=names.header)
        row_to_write = dict(
            id = args.id,
            model = args.model,
            x0 = args.x0,
            y0 = args.y0,
            flux = args.flux,
            hlr = args.hlr,
            e1 = args.e1,
            e2 = args.e2,
            psf_flux = args.psf_flux,
            psf_fwhm = args.psf_fwhm
        )
        #if file is empty just write the new row. 
        if csvIsEmpty(tempname):      
            writer.writeheader()
            writer.writerow(row_to_write)
        else:
            with open(tempname, 'r') as tempfile:
                reader = csv.DictReader(tempfile)
                writer.writeheader()        
                for row in reader:
                    print row
                    print row['id']
                    print args.id
                    if(int(row['id']) != args.id):
                        writer.writerow(row)
                writer.writerow(row_to_write)
                        

    os.remove(tempname)
        
    #print defaults from default.py in some format. 
    if args.defaults:
        pass

if __name__ == '__main__':
    main()
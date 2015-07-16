import argparse
import defaults
import os
import csv
import shutil
import functions as fns


def main():
    
    names = defaults.names() #initialize names used. 
    dflt_params = defaults.parameters().dict #initialize dictionary with default parameters 

    # Initialize and parse command-line arguments.
    parser = argparse.ArgumentParser(description = 'Generate galaxies specified by the user that are added to a file.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--defaults', action = 'store_true',
    help = 'Print defaults use for generating the different types of galaxies.')
    parser.add_argument('--wd', default = 'output', metavar = 'DIRECTORY', type = str,
    help = 'Specify a directory name where the output will be inputted and files read from.')
    parser.add_argument('--galaxy-file', default = 'galaxies', metavar = 'FILENAME', type = str,
    help = 'Specify a file where galaxies will be registered or read from, do not specify .csv at the end.')
    parser.add_argument('--id', required = True, default = None, metavar = 'GALAXYID', type = int,
    help = 'Add galaxy with id GALAXYID, if GALAXYID already exists in the table, you will change the parameters of the galaxy with that GALAXYID but not create a new entry.')
    parser.add_argument('--model', default = names.galaxy_models[0], metavar = 'MODEL', type = str, choices = names.galaxy_models,
    help = 'Change the galaxy model to MODEL.')
    parser.add_argument('--psf_model', default = names.psf_models[0], metavar = 'PSF_MODEL', type = str, choices = names.psf_models,
    help = 'Change the psf model to PSF_MODEL.')

    #add all parameter argument to the parser. this can change as time goes on.
    for name in names.parameters:
        parser.add_argument('--' + name, default = dflt_params[name], metavar = name.upper(), type = float, help = 'Change the parameter ' + name + ' from its default value.') 

    args = parser.parse_args()

    #check if directory exits. 
    if not os.path.isdir(args.wd):
        os.mkdir(args.wd)

    #initialize name of files that are required for storing galaxy data. 
    filename = os.path.join(args.wd, args.galaxy_file + '.csv')
    if not os.path.isfile(filename):
        f = open(filename, 'w+')
        f.close()
    tempname = os.path.join(args.wd, 'temp' + '.csv')
    shutil.copyfile(filename, tempname) #creat a copy to read from and compare.
    
    #write galaxy data to a filename.
    with open(filename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=names.fieldnames)
        row_to_write = dict(
            id = args.id,
            model = args.model,
            x0 = args.x0,
            y0 = args.y0,
            flux = args.flux,
            hlr = args.hlr,
            e1 = args.e1,
            e2 = args.e2,
            psf_model = args.psf_model,
            psf_flux = args.psf_flux,
            psf_fwhm = args.psf_fwhm
        )
        #if file is empty just write the new row. 
        if fns.csvIsEmpty(tempname):      
            writer.writeheader()
            writer.writerow(row_to_write)
        else:
            with open(tempname, 'r') as tempfile:
                reader = csv.DictReader(tempfile)
                writer.writeheader()        
                for row in reader:
                    if(int(row['id']) != args.id):
                        writer.writerow(row)
                writer.writerow(row_to_write)
                        

    os.remove(tempname)
        
    #print appropiate defaults for creating galaxies from default.py in some format. 
    if args.defaults:
        pass

if __name__ == '__main__':
    main()
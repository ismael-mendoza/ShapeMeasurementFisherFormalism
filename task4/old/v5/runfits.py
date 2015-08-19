import os
import csv
import defaults
import sys
from functions import *
from lmfit import Parameters, minimize, fit_report
import math

def objfuncChi(fit_params, data, variance_noise, **kwargs): 
    model = drawGalaxy(fit_params.valuesdict(), **kwargs).array.ravel()
    return (model - data) / math.sqrt(variance_noise)

def main(argv):

    names = defaults.names()
    dflt_params = defaults.parameters().dflt_params
    noise_seed, snr, wdir, galaxy_file, rltsdir = int(argv[1]), int(argv[2]), argv[3], argv[4], argv[5]

    if not os.path.isdir(wdir):
        print ('Directory does not exists')
        return -1

    filename = os.path.join(wdir, galaxy_file + '.csv')
    if not os.path.isfile(filename):
        print('Galaxies file does not exist')
        return -1

    if not os.path.isdir(os.path.join(wdir, rltsdir)):
        os.mkdir(os.path.join(wdir, rltsdir))

    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        lst_params = [] #possible here create dictionary with total number of params of both galaxies.
        for row in reader: 
            lst_params.append(row)

    params = lst_params[0] 
    for key, value in params.iteritems():
        try:
            params[key] = float(value)
        except ValueError:
            pass

    gal_image, variance_noise = drawGalaxy(params, SNR = snr, noise_seed = noise_seed)

    mins = defaults.min(gal_image).min_values
    maxs = defaults.max(gal_image).max_values

    fit_params = Parameters()
    for param in names.galaxy_parameters[params['model']]:
        fit_params.add(param, value = dflt_params[param], min = mins[param], max = maxs[param])
    
    nfit_params = dict()
    for param in params.keys(): 
        if param not in names.galaxy_parameters[params['model']]:
            nfit_params[param] = params[param]

    minimize(objfuncChi, fit_params, kws=dict(data = gal_image.array.ravel(), variance_noise = variance_noise, **nfit_params))
    print params
    #printing fit.
    print("")
    print("Start fit_report:")
    print(fit_report(fit_params))

    # result_filename = os.path.join(wdir,rltsdir,'results' +str(noise_seed) + '.cnv')
    # with open(result_filename, 'w') as csvfile: #'a' for appending and not 'w' for overwriting.
    #     row_to_write = fit_params.valuesdict()
    #     # for param in fit_params.keys():
    #     #     row_to_write[param] = fit_params[param].value
    #         #row_to_write[parameter_name + '_stderr'] = fit_params[parameter_name].stderr #maybe not needed, it depends apparently, leave it out for now.
    #     writer = csv.DictWriter(csvfile, fieldnames=row_to_write.keys())
    #     writer.writeheader()
    #     writer.writerow(row_to_write)

if __name__ == '__main__':
    main(sys.argv)
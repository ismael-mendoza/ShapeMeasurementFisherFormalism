#run various tests to validate results from fisher formalism.
#uses default wdir and file_name for simplicity.
#single galaxy. 
#the sigma used and parameters are all defined in refregier paper.
#also id of the galaxy must be 1. 

import math

import sys 

import defaults

import galfun

import fisher
#parameters will always be (for simplicity):['hlr','flux','x0','y0','e1','e2']

def bias(paramFunc, params, param_names, biases_of_params, steps):
    bias = 0 
    for k in range(len(param_names)):
            bias += fisher.partialDifferentiate(paramFunc, param_names[k], steps)(params) * biases_of_params[param_names[k]]
    return bias

def variance(paramFunc1, paramFunc2, params, param_names, 
            CovM_of_params, steps):
    """Uses the existing covariance matrix of the initial parameters to 
    compute the variance of other parameters that were not included initially
    in the model
    """
    var = 0 
    for k in range(len(param_names)):
        for l in range(len(param_names)):
            var += fisher.partialDifferentiate(paramFunc1, param_names[k], steps)(params) * fisher.partialDifferentiate(paramFunc2, param_names[l], steps)(params) * CovM_of_params[param_names[k],param_names[l]]
    return var


def radiusMean_func(params):
    return math.sqrt(a1_func(params)**2 + a2_func(params)**2)


def q_func(params): 
    if('q_1' in params.keys()): 
        return params['q_1']

    elif('e1_1' and 'e2_1' in params.keys()): 
        return math.sqrt((1-math.sqrt(params['e1_1']**2 + params['e2_1']**2))/(1 + math.sqrt(params['e1_1']**2 + params['e2_1']**2)))


def beta_func(params):
    """Calculates the galaxy parameter beta in terms of the other 6 initial known parameters"""

    if('beta_1' in params.keys()):
        return params['beta_1']

    elif('e1_1' and 'e2_1' in params.keys()):
        return .5 * math.atan(params['e2_1'] /params['e1_1']) #should be in radians. 


def a1_func(params):
    """Calculates the galaxy parameter a1 in terms of the other 6 initial known parameters"""
    return a2_func(params) / q_func(params)


def a2_func(params):
    """Calculates the galaxy parameter a2 in terms of the other 6 initial known parameters"""

    if('hlr_1' in params.keys()):
        return sigma_func(params) * math.sqrt(q_func(params))
    else:
        raise ValueError

def sigma_func(params):
    #the log 4 comes from galsim and tis definition of hlr and sigma in 
    #Gaussians.
    return params['hlr_1'] / math.sqrt(math.log(4))


def amplitude_func(params): 
    """Calculates the galaxy parameter amplitude in terms of the other 6 initial known parameters"""

    return params['flux_1']/sigma_func(params)

def flux_func(params):
    return params['flux_1']

def main(argv):

    names = defaults.names()
    g_parameters = galfun.GParameters(wdir = names.wdir, 
                                      galaxy_file= names.galaxy_file)
    fish = fisher.fisher(g_parameters = g_parameters, snr = 60)

    params = g_parameters.model_params
    #print params

########test against refregier, this is a single 6 parameter shaped galaxy. 

    #sigma[a]/a *p[A]

    rho_A = (math.sqrt(variance(amplitude_func,amplitude_func,params,params.keys(), fish.covariance_matrix,fish.steps)) /amplitude_func(params))**(-1)

    #a1 = sqrt(2)
    sigma_a1 = math.sqrt(variance(a1_func,a1_func,params,params.keys(), fish.covariance_matrix,fish.steps))

    print 'sqrt(2)'
    print (sigma_a1 /a1_func(params)) * rho_A


    #flux = sqrt(2)
    sigma_flux = math.sqrt(fish.covariance_matrix['flux_1','flux_1'])
    print 'sqrt(2)'
    print (sigma_flux /params['flux_1']) * rho_A

    #b[a]/a *p[A]**2

    #A = 5/2
    print '5/2'
    bias_A = bias(amplitude_func, params, params.keys(), fish.biases, fish.steps)
    print (bias_A/amplitude_func(params)) * rho_A**2

    #flux = 5/2
    print '5/2'
    bias_flux = fish.biases['flux_1']
    print (bias_flux/params['flux_1']) * rho_A**2

    #radius mean = 1
    print '1'
    bias_radius_mean = bias(radiusMean_func, params, params.keys(), fish.biases, fish.steps)
    print (bias_radius_mean/radiusMean_func(params))*rho_A**2

####end test of refregier. 

if __name__ == '__main__':
    main(sys.argv)
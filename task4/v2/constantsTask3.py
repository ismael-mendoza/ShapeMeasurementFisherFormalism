

#global constants for image (never going to use in a fit): (do not put in dictionary)
pixel_scale = .2
nx = 40
ny = 40
num_params = 6 


#variables for adjusting space between plots. 

left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = .4   # the amount of width reserved for blank space between subplots
hspace = 1.0   # the amount of height reserved for white space between subplots

# #parameter names used throughout the task.
# gal_param_names = ['x0', 'y0', 'hlr', 'flux', 'e1', 'e2']
# psf_param_names =['flux', 'HLR']
#it is much better if we just pass the parameter names around from the first file and get the names from it. 
#much more robust and less confusing. 


# initial_data_filename = 'initial_data.csv'
# results_folders = 'results' #name of results folders
# trial_folders = 'trial' #name of trial folders. 
#same with this folder names. 



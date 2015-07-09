"""This gathers the residuals from the fit.csv file (where all the runs of the fits inputted their data), considering that all the data is in one file."""
import csv
import numpy as np 

galaxy_params_names = ['gal_flux', 'gal_sigma', 'e1', 'e2', 'x0', 'y0']
# galaxy_params_values = [parameter + '_value' for parameter in galaxy_params_names]
# galaxy_params_stderr = [parameter + '_stderr' for parameter in galaxy_params_names]

#get true values
with open('initial_data.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	init_params = dict() 
	row = reader.next()
	for parameter_name in row.keys():
		init_params[parameter_name] = row[parameter_name]

with open('fit.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	residuals = dict()
	#create a list in each residual parameter_name entry to store residuals in each and then obtain mean. 
	for parameter_name in galaxy_params_names:
		residuals[parameter_name] = []
	for row in reader:
		for parameter_name in galaxy_params_names:
			#calculate residula for such a row. 
			residuals[parameter_name].append(float(init_params[parameter_name]) - float(row[parameter_name + '_value']))

biases = {parameter_name:np.mean(residuals[parameter_name]) for parameter_name in residuals.keys()}
print biases
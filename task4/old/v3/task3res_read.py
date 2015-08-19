"""This gathers the residuals from the single fit.csv file (where all the runs of the fits inputted their data), considering that all the data is in one file."""
import csv
import numpy as np 
from constantsTask3 import *

#get true values
with open('initial_data.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	init_params = dict() 
	row = reader.next()
	for param_name in row.keys():
		init_params[param_name] = row[param_name]

with open('fit.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	residuals = dict()
	#create a list in each residual param_name entry to store residuals in each and then obtain mean. 
	for param_name in params_names:
		residuals[param_name] = []
	for row in reader:
		for param_name in params_names:
			#calculate residula for such a row. 
			residuals[param_name].append(float(init_params[param_name]) - float(row[param_name + '_value']))

biases = {param_name:np.mean(residuals[param_name]) for param_name in residuals.keys()}
print biases
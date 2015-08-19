"""This gathers the residuals and computes the bias from multiple .csv files scattered in a directory. 
Each file will not necessarily have only one row of results (although it will most of the times"""

import csv
import numpy as np 
import os
import sys
from constantsTask3 import *


def main(argv):
    if(len(argv) == 3):
		results_number = argv[1] 
    	trial_number = argv[2]
    else:
    	print('need exactly 2 arguments to run program')
    	print('usage: ./task3res_readMultiple.py,result_number, trial_number')
    	return -1



	#path to all the files to read. 
	path = '/Users/Ismael/code/research/task3/v24'
	results_dir = results_folders + results_number
	trial_dir = trial_folders + trial_number

	#get true values
	initial_data_filename = os.path.join(result_dir, initial_data_filename)
	with open(initial_data_filename) as csvfile:
		reader = csv.DictReader(csvfile)
		init_params = dict() 
		row = reader.next()
		for param_name in row.keys():
			init_params[param_name] = row[param_name]

	residuals = dict()
	#create a list in each residual param_name entry to store residuals in each and then obtain mean. 
	for param_name in param_names:
		residuals[param_name] = []

	path = os.path.join(path,results_dir,trial_dir)
	for filename in os.listdir(path):
		with open(os.path.join(path,filename)) as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				for param_name in param_names:
					#calculate residual for such a row. 
					residuals[param_name].append(float(init_params[param_name]) - float(row[param_name + '_value']))

	biases = {param_name:np.mean(residuals[param_name]) for param_name in residuals.keys()}
	print biases

if __name__ == "__main__":
    main(sys.argv)
"""This gathers the residuals from the fits.csv file (where all the runs of the fits inputted their data)"""
import csv

#galaxy_params = ['gal_flux', 'gal_sigma', 'e1', 'e2', 'x0', 'y0']
#get true values
with open('initial_data.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	init_params = dict() 
	row = reader.next()
	for parameter_name in row.keys():
		init_params[parameter_name] = row[parameter_name]

with open('fits.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	residuals = dict()
	for row in reader:
		for parameter_name in galaxy_params:
		residuals[parameter_name] = [row[parameter_name]]
		
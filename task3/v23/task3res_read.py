"""This gathers the residuals from the fits.csv file (where all the runs of the fits inputted their data)"""
import csv

with open('fits.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	residuals = [] 
	for row in reader:
		residuals.append(object)
		
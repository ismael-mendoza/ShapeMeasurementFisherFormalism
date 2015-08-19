"""This gathers the residuals and computes the bias from multiple .csv files scattered in a directory. 
Each file will not necessarily have only one row of results (although it will most of the times"""

import csv

import numpy as np 

import os

import sys

import functions as fns

import triangle

import fisher

def main(argv):
    wdir, galaxy_file, rltsdir, snr = argv[1], argv[2], argv[3], argv[4]

    params = fns.getParamsCsv(wdir, galaxy_file)

    if not os.path.isdir(os.path.join(wdir, rltsdir)):
        os.mkdir(os.path.join(wdir, rltsdir))

    residuals = {}
    path = os.path.join(wdir,rltsdir)

    for filename in os.listdir(path):
        with open(os.path.join(path,filename)) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for param in row.keys():
                    #store every residual in a list and then can calculate mean. 
                    if not residuals[param] == []: 
                        residuals[param] = [] 
                    #calculate residual for each row. 
                    residuals[param].append(float(params[param]) - float(row[param]))

    biases = {param:np.mean(residuals[param]) for param in residuals.keys()}
    fisher_analysis = fisher.fisher(params, fns.drawGalaxy(params), snr)
    #produce triangle plots.


    #draw Gaussians, 
            #     import matplotlib.pyplot as plt
            # import numpy as np
            # import matplotlib.mlab as mlab
            # import math

            # mean = 0
            # variance = 1
            # sigma = math.sqrt(variance)
            # x = np.linspace(-3,3,100)
            # plt.plot(x,mlab.normpdf(x,mean,sigma))

            # plt.show()

    #draw error ellipses, 
        #triangle.error_ellipse()


    for i in range(fisher.num_params): 
        for j in range(fisher.num_params):
            if(i >= j):
                ax = figure.add_subplot(self.num_params,self.num_params, self.num_params * i + j + 1)
                fns.drawImage(ax, self.fisher.second_derivatives_images[self.param_names[i],self.param_names[j]])
                # ax.text(0,20, "std:" + str(fisher.second_derivatives_galaxy_images[self.param_names[i],self.param_names[j]].std().round(4)), fontsize = 8, fontweight='bold') 
                if(j == 0):
                    ax.set_ylabel(self.param_names[i])
                if(i == self.num_params - 1):
                    ax.set_xlabel(self.param_names[j])

if __name__ == "__main__":
    main(sys.argv)
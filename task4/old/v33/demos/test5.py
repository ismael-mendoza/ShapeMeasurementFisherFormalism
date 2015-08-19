import numpy as np
import triangle
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

mean = 1
variance = .00005
sigma = math.sqrt(variance)
x = np.linspace(-3,3,1000)
figure = plt.figure()
for i in range(9):    
    ax= figure.add_subplot(3,3,i+1)
    ax.plot(x,mlab.normpdf(x,mean,sigma))


ndim, nsamples = 3, 10000
samples = np.random.randn(ndim * nsamples).reshape([nsamples, ndim])
figurefinal = triangle.corner(samples, fig = figure)

figurefinal.show()
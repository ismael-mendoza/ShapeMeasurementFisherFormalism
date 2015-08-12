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

# Set up the parameters of the problem.
ndim, nsamples = 3, 50000

# Generate some fake data.
data1 = np.random.randn(ndim * 4 * nsamples / 5).reshape([4 * nsamples / 5,
                                                          ndim])
data2 = (5 * np.random.rand(ndim)[None, :]
         + np.random.randn(ndim * nsamples / 5).reshape([nsamples / 5, ndim]))
data = np.vstack([data1, data2])

# Plot it.
figurefinal = triangle.corner(data, labels=[r"$x$", r"$y$", r"$\log \alpha$",
                                       r"$\Gamma \, [\mathrm{parsec}]$"],
                         truths=[0.0, 0.0, 0.0],
                         quantiles=[0.16, 0.5, 0.84],
                         show_titles=True, title_args={"fontsize": 12}, fig = figure)
figurefinal.gca().annotate("A Title", xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")

figurefinal.show()
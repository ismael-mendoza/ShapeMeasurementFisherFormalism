import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import math
from matplotlib.patches import Ellipse

ellipse = Ellipse(xy=(0, 0), width=0.036, height=0.012, 
                        edgecolor='r', fc='None', lw=2)
mean = 0
variance = 1
sigma = math.sqrt(variance)
x = np.linspace(-3,3,1000)
figure = plt.figure()
ax= figure.add_subplot(111)
ax.plot(x,mlab.normpdf(x,mean,sigma))
ax.add_patch(ellipse)
ax.autoscale(True)

plt.show()

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import math

mean = 0
variance = 1
sigma = math.sqrt(variance)
x = np.linspace(-3,3,1000)
figure = plt.figure()
ax= figure.add_subplot(111)
ax.plot(x,mlab.normpdf(x,mean,sigma))

figure.savefig('1.png')

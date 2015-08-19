from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

figure = plt.figure()
ax = figure.add_subplot(111)

ellipse = Ellipse(xy=(0, 0), width=10, height=1, 
                        edgecolor='r', fc='None', lw=2)
ax.add_patch(ellipse)
# ax.set_xlim(0 - 2, 0 + 2)
# ax.set_ylim(0 - 2, 0 + 2)
ax.autoscale(True)
plt.show()
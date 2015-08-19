from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

plt.figure()
ax = plt.gca()

ellipse = Ellipse(xy=(157.18, 68.4705), width=0.036, height=0.012, 
                        edgecolor='r', fc='None', lw=2)
ax.add_patch(ellipse)
plt.show()
#!/usr/bin/env python

import numpy as np
import triangle
import matplotlib.pyplot as plt
from pylab import figure, show, rand
from matplotlib.patches import Ellipse



fig = figure()
ax = fig.add_subplot(111, aspect='equal')
x = [1, -1,  4.3]
y = [3,  1.1,  0.12]
mu = np.mean(x), np.mean(y)
X = np.vstack((x,y))
cov = np.cov(X)
ellipse = triangle.error_ellipse(mu = mu, cov = cov)
ax.add_artist(ellipse)
ax.set_xlim(np.mean(x) - 2, np.mean(x) + 2)
ax.set_ylim(np.mean(y) - 2, np.mean(y) + 2)

show()





# ell = Ellipse(xy=rand(2)*10, width=rand(), height=rand(), angle=rand()*360)

# fig = figure()
# ax = fig.add_subplot(111, aspect='equal')
# ax.add_artist(ell)
# ell.set_clip_box(ax.bbox)
# ell.set_alpha(rand())
# ell.set_facecolor(rand(3))

# ax.set_xlim(0, 10)
# ax.set_ylim(0, 10)

# show()
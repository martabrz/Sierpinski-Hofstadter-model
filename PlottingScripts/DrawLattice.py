import numpy as np
import matplotlib.pyplot as plt
from beautiful_latex import latex_plot

latex_plot(scale = 0.9)
filename = 'xy_positions.dat'
data = np.genfromtxt(filename)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_rasterization_zorder(1);
cax = ax.scatter(data[:,0], data[:,1], s = 0.5, marker='s', edgecolors='none', color = 'black')
ax.set_xticks([])
ax.set_yticks([])
fig.savefig(filename[0:len(filename)-4] + '.png', dpi = 800)

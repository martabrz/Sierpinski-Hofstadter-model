import numpy as np
import matplotlib.pyplot as plt
from beautiful_latex import latex_plot

latex_plot(scale = 0.9)

filename = 'IPR_000.dat'
filename2 = 'ipr_zeroflux.dat'
filename3 = 'ipr_037flux.dat'
filename4 = 'ipr_01flux.dat'
data_01flux = np.genfromtxt(filename)

print(np.min(data_01flux[:,0]), np.max(data_01flux[:,0]))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_rasterization_zorder(1);
binwidth = 0.0005
# ax.plot([0, 5000], [1/(4*81), 1/(4*81)], '-', color = 'xkcd:gold', lw = 1, label = 'l = 4')
ax.plot([0, 5000], [1/(112), 1/(112)], '-', color = 'xkcd:darkgreen', lw = 1, label = 'l = 3')
ax.plot([0, 5000], [1/(40*8), 1/(40*8)], '-', color = 'red', lw = 1, label = 'l = 2')
ax.plot([0, 5000], [1/(1024), 1/(1024)], '-', color = 'xkcd:orange', lw = 1, label = 'l = 1')
ax.plot([0, 5000], [1/(2465), 1/(2465)], '--', color = 'xkcd:crimson', lw = 1, label = '1/sites')

ax.hist(data_01flux[:,0], bins=np.arange(min(data_01flux[:,0]), max(data_01flux[:,0]) + binwidth, binwidth),\
orientation= 'horizontal', facecolor = 'xkcd:blue', edgecolor = 'black',  linewidth = 0.5);
# ax.set_ylim(0, 0.01)
ax.legend()
ax.set_xscale('log')
# ax.set_yticks(np.linspace(-4, 5, 10))
ax.set_ylim(0, 0.04)
ax.set_xlim(1-0.1, 5000)
# ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_title(r'$\alpha = 0.00$')
ax.set_ylabel('$I_{\psi}$')
ax.set_xlabel('\# states')
fig.savefig('IPR_histogramflux000.svg', dpi = 800)

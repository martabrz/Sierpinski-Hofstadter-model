import numpy as np
import matplotlib as mpl
from matplotlib.colors import LogNorm
import itertools
import sys
from beautiful_latex import latex_plot
latex_plot(scale = 0.9)

my_cmap = mpl.cm.get_cmap('Reds')
my_cmap.set_under('k')
norma = mpl.colors.Normalize(vmin=1e-10, vmax=1)


filename = 'eval.dat'
data = np.loadtxt(filename)
flux = np.loadtxt('flux.dat')
alpha = flux/(2*np.pi)
datar = data.reshape(len(alpha), -1)
eta = 1e-5
E_a = -6
E_b = 6
asp = (max(alpha) - min(alpha))/(0.5*(abs(E_a - E_b)))
E = np.linspace(E_a, E_b, len(alpha))
grid = np.zeros((len(E), len(alpha)))
l = datar.shape[1]

# DOS with Gaussian smearing
print(np.shape(datar))
for i, j in itertools.product(range(len(E)), range(len(alpha))):
    temp = 0
    for k in range(l):
        temp += np.exp((- (E[i] - datar[j, k])**2)/eta)
    grid[i,j] = temp

# data normalization
# normalization = []
# for k in grid:
#     for i in k:
#         normalization.append(i)
# n = np.max(normalization)
# grid /= n
# print(grid)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_rasterization_zorder(1);
cax = ax.imshow(grid, interpolation = 'nearest', \
extent = [min(alpha), max(alpha), E_a, E_b], cmap = my_cmap, \
aspect = 'auto', origin = 'lower', zorder = 0, norm = norma)
# cbar = fig.colorbar(cax)
# cbar.set_ticks([])
ax.set_ylabel('$E/t$')
ax.set_ylim(-4,1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax.set_yticks([-4, -3, -2, -1, 0, 1])
ax.set_xlabel('$p/q$')
fig.savefig(filename[0:len(filename)-4] + 'TEST.svg', \
transparent = True, dpi = 800, bbox_inches = 'tight')

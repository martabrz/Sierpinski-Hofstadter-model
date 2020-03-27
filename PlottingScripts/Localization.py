import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from beautiful_latex import latex_plot
latex_plot(scale = 0.9)

filename = 'localization.dat'
data = np.genfromtxt(filename)
x = data[:,0]
y = data[:,1]
siz = data[:,2]/max(data[:,2])
norma = mpl.colors.Normalize(vmin=0, vmax=1)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_rasterization_zorder(1);

cdict = cm.get_cmap('Reds')._segmentdata
cdict['red'][0] = (0, 0.3, 0.5) # x=0 for bottom color in colormap
cdict['blue'][0] = (0, 0.5, 0.5) # y=0.5 gray
cdict['green'][0] = (0, 0.5, 0.5) # y1=y for simple interpolation
cdict['red'][-1] = (1, 1, 1) # x=1 for top color in colormap
cdict['blue'][-1] = (1, 0.5, 0.5)
cdict['green'][-1] = (1, 0.5, 0.5)
my_cmap = LinearSegmentedColormap('name', cdict)
cax = ax.scatter(x, -y, marker = 's', edgecolors = 'none', s = 1.5+siz, \
c = siz, norm = norma, cmap = my_cmap, alpha = 1, rasterized = True)

# cbar = fig.colorbar(cax)
# cbar.ax.set_title('$|\psi|^2$')
# ax.patch.set_facecolor('black')
# ax.patch.set_alpha(0.9)
ax.set_xticks([])
ax.set_yticks([])
fig.savefig('2140.svg', \
transparent = True, dpi = 1000, bbox_inches = 'tight')

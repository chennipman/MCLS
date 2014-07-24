# compare the results of the numerical output and the analytical results of the Taylor Vortex

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np


# Analytical Result




# make these smaller to increase the resolution
dx, dy = 0.1, 0.1

# generate 2 2d grids for the x & y bounds
y, x = np.mgrid[slice(0, 2 + dy, dy),
                slice(0, 2 + dx, dx)]

# u,v,pres
expo = exp(-2*pi*pi*0.01)
u = (-np.sin(x*math.pi)*np.cos(y*math.pi))*expo

print u

# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
z = z[:-1, :-1]
levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

plt.subplot(2, 1, 1)
im = plt.pcolormesh(x, y, z, cmap=cmap, norm=norm)
plt.colorbar()
# set the limits of the plot to the limits of the data
plt.axis([x.min(), x.max(), y.min(), y.max()])
plt.title('pcolormesh with levels')



plt.subplot(2, 1, 2)
# contours are *point* based plots, so convert our bound into point
# centers
plt.contourf(x[:-1, :-1] + dx / 2.,
             y[:-1, :-1] + dy / 2., z, levels=levels,
             cmap=cmap)
plt.colorbar()
plt.title('contourf with levels')


plt.show()

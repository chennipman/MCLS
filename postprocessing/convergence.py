# makes plot of convergence

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from readvtk import readvtkfunction

x_coordinates = []
y_coordinates = []
z_coordinates = []
velocity_u = []
velocity_v = []
velocity_w = []
vof = []
level_set = []
curvature = []
unsmoothed_curvature = []
pressure = []

filename = '../executable/case_TV/07_13_16_33_22_d94ee5afac97d6e35ceaced0814681ea33012b2d/solution_file.1.pure'

(x_coordinates, y_coordinates, z_coordinates, velocity_u, velocity_v, velocity_w, vof, level_set, curvature, unsmoothed_curvature, pressure) = readvtkfunction(filename)

x_center = x_coordinates[0:20]+[0.05]
y_center = y_coordinates[0:20]+[0.05] 

velocity_u_part = velocity_u[:,:,1]
velocity_v_part = velocity_v[:,:,1]


plt.title("velocity_u ")

levels = MaxNLocator(nbins=200).tick_values(-1, 1)


cmap = plt.get_cmap('PiYG')
cmap = plt.get_cmap('seismic')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

plt.contourf(x_coordinates,
             y_center, velocity_u_part[:,1:21].T, levels=levels,cmap=cmap)

plt.colorbar()
plt.show()




# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 20:39:09 2014
@author: coen
Creates contour plots of the results

"""


import sys
sys.path.append('general_functions')

#import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from hex2colormap import hex2colormap
from readvtk import readvtkfunction


# TU Delft basic colors
glaucous  = '#66BCAA'
violet    = '#0F1150' 
green     = '#007A85' 
turquoise = '#0093AB'
darkblue  = '#002B60' 
skyblue   = '#77C0D7' 

# TU Delft accent colors
lavender  = '#7BA0C9' 
fuchsia   = '#A10058'
orange    = '#EC7F2C'
yelgreen  = '#ADC610'
hotpurple = '#83267C'
yellow    = '#F7EB90' 


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
pressure        = pressure[:,:,0]


print 'something' 

# plot the variables
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

colormap = hex2colormap(yellow,hotpurple)
cmap = plt.get_cmap(colormap)
levels = MaxNLocator(nbins=50).tick_values(-1, 1)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

plt.ion()
plt.show(block=False)

plt.figure()
plt.title("velocity $u$ ")
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
plt.contourf(x_coordinates, y_center, velocity_u_part[:,1:21].T, levels=levels,cmap=cmap)
plt.colorbar()
print 'start saving figure' 
plt.savefig('../figures/TV/u-velocity.pdf')
print 'complete saving figure' 
plt.close

plt.figure()
plt.title("velocity $v$ ")
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
plt.contourf(x_center, y_coordinates, velocity_v_part[1:21,:].T, levels=levels,cmap=cmap)
plt.colorbar()
print 'start saving figure 2' 
plt.savefig('../figures/TV/v-velocity.pdf')
print 'complete saving figure 2' 
plt.close

plt.figure()
plt.title("pressure")
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
levels = MaxNLocator(nbins=50).tick_values(0, 1)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.contourf(x_center, y_center, pressure, levels=levels,cmap=cmap)
plt.colorbar()
print 'start saving figure 3' 
plt.savefig('../figures/TV/pressure.pdf')
print 'complete saving figure 3' 
plt.close
plt.show()

print 'ready' 







# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 22:02:41 2014

@author: coen
"""

import numpy as np
import matplotlib.pyplot as plt
from hex2colormap import hex2colormap
from hex2colormap import hex3colormap


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

x = np.arange(0, np.pi, 0.1)
y = np.arange(0, 2*np.pi, 0.1)
X, Y = np.meshgrid(x,y)
Z = np.cos(X) * np.sin(Y) * 10

# Make the figure:

plt.figure()

colormapname = hex2colormap(yellow,skyblue)
plt.subplot(1,2,1)
plt.imshow(Z, interpolation='nearest', cmap=colormapname)
plt.colorbar()

# the three colors are not working well
colormapname = hex3colormap(hotpurple,yelgreen,skyblue)
plt.subplot(1,2,2)
plt.imshow(Z, interpolation='nearest', cmap=colormapname)
plt.colorbar()

plt.show()


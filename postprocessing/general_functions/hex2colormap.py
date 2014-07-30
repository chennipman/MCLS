# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 21:48:45 2014

@author: coen


Input:   2 Hex
Output: colormap

"""

def hex2colormap(hex1, hex2):

    from matplotlib.colors import hex2color
    from matplotlib.colors import LinearSegmentedColormap
    
    some_name = 'some_name' 
    
    [hex1r, hex1g, hex1b] = hex2color(hex1)
    [hex2r, hex2g, hex2b] = hex2color(hex2)    

    
    cdict = {'red':   ((0.0, hex1r, hex1r),
                       (1.0, hex2r, hex2r)),
    
             'green': ((0.0, hex1g, hex1g),
                       (1.0, hex2g, hex2g)),
    
             'blue':  ((0.0, hex1b, hex1b),
                       (1.0, hex2b, hex2b))
            }
    
    colormapname2  =   LinearSegmentedColormap(some_name, cdict)
    return colormapname2
    
def hex3colormap(hex1, hex2, hex3):

    from matplotlib.colors import hex2color
    from matplotlib.colors import LinearSegmentedColormap
    
    some_name = 'some_name' 
    
    [hex1r, hex1g, hex1b] = hex2color(hex1)
    [hex2r, hex2g, hex2b] = hex2color(hex2)    
    [hex3r, hex3g, hex3b] = hex2color(hex3)    

    
    cdict = {'red':   ((0.0, hex1r, hex1r),
                       (0.5, hex2r, hex2r),
                       (1.0, hex3r, hex3r)),
    
             'green': ((0.0, hex1g, hex1g),
                       (0.5, hex2r, hex2r),
                       (1.0, hex3g, hex3g)),
    
             'blue':  ((0.0, hex1b, hex1b),
                       (0.5, hex2r, hex2r),
                       (1.0, hex3b, hex3b))
            }
    
    colormapname2  =   LinearSegmentedColormap(some_name, cdict)
    return colormapname2
    
    
    
    

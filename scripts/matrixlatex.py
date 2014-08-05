# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 14:39:44 2014

@author: coen
"""
def matrixlatexwrite(x_top_row, y_first_column, matrix, filename, size_x, size_y):
    f = open(filename ,'w')
    f.write(r"y\textbackslash x &")
    for i in range(0,size_x):
        val = x_top_row[i]
        formated = '%.2f' % val
        f.write('%s' % formated)
        if i != size_x-1:                # not last column
            f.write(" & ")
        else:                       # last column
            f.write(r"\\ \hline \hline")
            f.write("\n")        
                
    for j in range(0,size_y):
        for i in range(0,size_x):
            if i == 0:
                val = y_first_column[j]
                formated = '%.2f' % val
                f.write('%s' % formated)
                f.write(" & ")
            e = matrix[i][j]
            formated = '%.5f' % e
            f.write('%s' % formated)
            if i != size_x-1:           # not last column
                f.write(" & ")
            elif j != size_y-1:         # last column, not last row
                f.write(r" \\ \hline")
                f.write("\n")     
            else:                       # last column, last row
                f.write(" ")
                       
                
                
                
import sys
sys.path.append('general_functions')
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

filename = '../executable/case_TV/08_04_19_46_25_8dd559b2a42be3c35e914d28215372a7723fe07c/solution_file.1.pure'

(x_coordinates, y_coordinates, z_coordinates, velocity_u, velocity_v, velocity_w, vof, level_set, curvature, unsmoothed_curvature, pressure) = readvtkfunction(filename)

x_center = x_coordinates[0:20]+[0.05]
y_center = y_coordinates[0:20]+[0.05] 

velocity_u_part = velocity_u[:,:,1]
velocity_v_part = velocity_v[:,:,1]
pressure        = pressure[:,:,0]

# make tables
foldername = '../tables/TV/'

# write u to table
filenametable = foldername + 'u-velocity-first-part.tex' 
size_x = 10
size_y = 20
matrixlatexwrite(x_coordinates[0:size_x], y_center, velocity_u_part[0:size_x,:], filenametable, size_x, size_y)

size_x_2 = 11
filenametable = foldername + 'u-velocity-second-part.tex'
matrixlatexwrite(x_coordinates[size_x:size_x_2+size_x], y_center, velocity_u_part[size_x:size_x_2+size_x ,:], filenametable, size_x_2, size_y)

# write v to table
filenametable = foldername + 'v-velocity-first-part.tex' 
size_x = 10
size_y = 21
matrixlatexwrite(x_center[0:size_x], y_coordinates, velocity_v_part[0:size_x,:], filenametable, size_x, size_y)

size_x_2 = 10
filenametable = foldername + 'v-velocity-second-part.tex'
matrixlatexwrite(x_center[size_x:size_x_2+size_x], y_coordinates, velocity_v_part[size_x:size_x_2+size_x ,:], filenametable, size_x_2, size_y)

# write v to table
filenametable = foldername + 'pressure-first-part.tex' 
size_x = 10
size_y = 20
matrixlatexwrite(x_center[0:size_x], y_center, pressure, filenametable, size_x, size_y)

size_x_2 = 10
filenametable = foldername + 'pressure-second-part.tex'
matrixlatexwrite(x_center[size_x:size_x_2+size_x], y_center, pressure, filenametable, size_x_2, size_y)



print 'finished writing table'



   
   







    
    
    
#    f.write(r"\end{tabular}")


"""
    f.write(r"\begin{tabular}{l")
    for i in range(0,size_x):
        f.write("c")
    f.write("} \n")
"""

""" 
        # Values
    for i in range(0, size_x):
        for j in range(0, size_y):
            if j == 0: # first column
                e = y_first_column[j]
                formated = fcj % e
                f.write('%s' % formated) 
                f.write(" & ")
            e = matrix[i][j]
            fcj = '%.5g'
            formated = fcj % e
            formated = fix(formated, table=True) # fix 1e+2
            f.write('%s' % formated)
            if j != size_y-1:                # not last row
                f.write(" & ")
            else:                       # last row
                f.write(r"\\")
                f.write("\n")
                
  
     
    

    for i in range(0, size_x):
        e = x_top_row[i]
        formated = fcj % e
        formated = fix(formated, table=True) # fix 1e+2
        f.write('%s' % formated)
        f.write(" & ")
        if i == size_x:
            f.write("\\")         
            
"""
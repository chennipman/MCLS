# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 03:11:35 2014

Analyse the runs for different directions. The output in all directions should be the same.

Output: plot with rising velocity

"""



import numpy as np
import sys
sys.path.append('general_functions')
from readvtk import readvtkfunction
import os
#import matplotlib.pyplot as plt
from matplotlib.pyplot import *

lavender  = '#7BA0C9' 
fuchsia   = '#A10058'
orange    = '#EC7F2C'
yelgreen  = '#ADC610'
hotpurple = '#83267C'
yellow    = '#F7EB90' 

    
#go to basis directory a
os.chdir('../')

# choose files based on logfiles
folder = 'logfiles/'
logfile = '07-31-04-20-logfile'
namelogfile = folder + logfile + '.txt'
f = open(namelogfile, 'r')
lines = f.readlines() 


main_directions = []
middle_directions = []

line = lines[2]
line = line[1:-2]
main_directions.extend(np.array([str(val) for val in line.rstrip('\n').split(',') if val != '']))

line = lines[4]
line = line[1:-2]
middle_directions.extend(np.array([str(val) for val in line.rstrip('\n').split(',') if val != '']))

f.close()


main_directions = ['x','y','z']
middle_directions = ['x','y','z']


counter = 0
overview = np.array(range(12), dtype=str).reshape(6,2)

for main_direc in main_directions: 
    for middle_direc in middle_directions: 
        not_execute = 0;


   # do the the parameter file
    # the main directions is set
        if main_direc == 'x':
            
            if middle_direc == 'y':
                overview[counter,:] = [main_direc, middle_direc]
                counter+=1
                #last direction = z
                
            elif middle_direc == 'z':
                overview[counter,:] = [main_direc, middle_direc]
                counter+=1                
                #last direction = y
            else:
                # do nothing
                print 'restricted situation x'  
                not_execute = 1
                
        elif main_direc == 'y':
            
            if middle_direc == 'x':
                overview[counter,:] = [main_direc, middle_direc]
                counter+=1
                #last direction = z
                 
            elif middle_direc == 'z':
                overview[counter,:] = [main_direc, middle_direc]
                counter+=1                
                #last direction = x 
            else:
                # do nothing
                print 'restricted situation y'  
                not_execute = 1
                
        elif main_direc == 'z':

            if middle_direc == 'x':
                overview[counter,:] = [main_direc, middle_direc]
                counter+=1
                #last direction = y

            elif middle_direc == 'y':
                overview[counter,:] = [main_direc, middle_direc]
                counter+=1 
                #last direction = x 
            else:
                # do nothing
                print 'restricted situation z'  
                not_execute = 1

print overview

print overview[1,0]
print overview[0,0]

# visualize
rc('text', usetex=True)
figure()

for i in range(counter):
    print i
    filename = lines[i+5].rstrip(' \n') +'/interface_details.csv' 
    # time,enclosed_volume,centroid_x1,centroid_x2,centroid_x3,centroid_u1,centroid_u2,centroid_u3,centroid_vof_u1,centroid_vof_u2,centroid_vof_u3
    print filename
    my_data = np.genfromtxt(filename, delimiter=',')
    print overview[i,0]
    if overview[i,0] == 'x':
        column = 5
    elif overview[i,0] == 'y':
        column = 6              
    elif overview[i,0] == 'z':
        column = 7
            
    print my_data[1:,0]
    
    print my_data[1:,column]
    
    plot(my_data[1:,0],my_data[1:,column], label=overview[i,:])

for i in range(4,8):
    lastname = 'c1g1l' + str(i)
    filename = 'postprocessing/' + lastname + '.txt'     
    ref_data = np.genfromtxt(filename, delimiter='  ')        
    plot(ref_data[1:,0],ref_data[1:,4], label=lastname)    
    
xlabel('time $t$')
ylabel('rising velocity')    
legend()
show()         




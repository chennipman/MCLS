# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 13:06:23 2014

@author: coen

The input is a logfile with multiple runs. Also the \dt for multiple timesteps should be specified. 
Output: several plots with the convergence


"""
#import time
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



def rms(x, axis=None):
    return np.sqrt(np.mean(x**2, axis=axis))
    
#go to basis directory a
os.chdir('../')


time_stepping_methods = []
time_steps = []


# choose files based on logfiles
folder = 'logfiles/'
logfile = '07-30-23-42-logfile'
namelogfile = folder + logfile + '.txt'
f = open(namelogfile, 'r')
lines = f.readlines() 

line = lines[2]
line = line[1:-2]
time_stepping_methods.extend(np.array([float(val) for val in line.rstrip('\n').split(',') if val != '']))

line = lines[4]
line = line[1:-2]
time_steps.extend(np.array([float(val) for val in line.rstrip('\n').split(',') if val != '']))

f.close()

print time_steps
print type(time_steps)
print len(time_steps)
#


uLinf = np.zeros((len(time_stepping_methods),len(time_steps)-1))
vLinf = np.zeros((len(time_stepping_methods),len(time_steps)-1))
pLinf = np.zeros((len(time_stepping_methods),len(time_steps)-1))

uL2 = np.zeros((len(time_stepping_methods),len(time_steps)-1))
vL2 = np.zeros((len(time_stepping_methods),len(time_steps)-1))
pL2 = np.zeros((len(time_stepping_methods),len(time_steps)-1))

print uLinf


for (counter,time_stepping_method) in enumerate(time_stepping_methods):
    
    filerunnumber = (counter+1)*len(time_steps)+4 # the magical 4 comes from the 5 lines on top of the file and -1 for starting at 0    
    print 'counter: '+str(counter)
    print 'filerun: '+str(filerunnumber)
    
    filename = lines[filerunnumber].rstrip(' \n') +'/solution_file.1.pure' 
    print filename
    # define the benchmark solution
    (x_coordinates, y_coordinates, z_coordinates, velocity_u_bench, velocity_v_bench, velocity_w, vof, level_set, curvature, unsmoothed_curvature, pressure_bench) = readvtkfunction(filename)
     
    for (counter2,time_step_restriction_global) in enumerate(time_steps[:-1]):
        # get all the other solutions
        filerunnumber = counter2+counter*len(time_steps)+5 # the magical 5 comes from the 5 lines on top of the file   
        print 'counter2: '+str(counter2)
        print 'filerun: '+str(filerunnumber)
        
        filename = lines[filerunnumber].rstrip(' \n') +'/solution_file.1.pure' 
        print filename
        (x_coordinates, y_coordinates, z_coordinates, velocity_u, velocity_v, velocity_w, vof, level_set, curvature, unsmoothed_curvature, pressure) = readvtkfunction(filename)
        print velocity_u[3][4]
        print np.amax(velocity_u_bench-velocity_u)

        uLinf[counter,counter2] = np.amax(abs(velocity_u_bench-velocity_u))
        vLinf[counter,counter2] = np.amax(abs(velocity_v_bench-velocity_v))
        pLinf[counter,counter2] = np.amax(abs(pressure_bench-pressure))
        
        uL2[counter,counter2] = rms(velocity_u_bench-velocity_u)
        vL2[counter,counter2] = rms(velocity_v_bench-velocity_v)
        pL2[counter,counter2] = rms(pressure_bench-pressure)



# visualize
rc('text', usetex=True)

figure()       
u25 = plot(time_steps[:-1], uL2[0,:].T, ls='-', color=lavender, label='u25')
u24 = plot(time_steps[:-1], uL2[1,:].T, ls='-.', color=fuchsia, label='u24')
u23 = plot(time_steps[:-1], uL2[2,:].T, ls='--', color=orange, label='u23')
u21 = plot(time_steps[:-1], uL2[3,:].T, ls=':', color=yelgreen, label='u21')
yscale('log')        
xscale('log')
xlabel('timestep $\Delta t$')
ylabel('error $L_2$-norm')
legend( (time_stepping_methods[0], time_stepping_methods[1], time_stepping_methods[2], time_stepping_methods[3]), loc='upper left')
title('$u$-velocity')
savefig('figures/TV/u-velocity_L2.pdf')


figure()
v25 = plot(time_steps[:-1], vL2[0,:].T, ls='-', color=lavender, label='v25')
v24 = plot(time_steps[:-1], vL2[1,:].T, ls='-.', color=fuchsia, label='v24')
v23 = plot(time_steps[:-1], vL2[2,:].T, ls='--', color=orange, label='v23')
v21 = plot(time_steps[:-1], vL2[3,:].T, ls=':', color=yelgreen, label='v21')
yscale('log')        
xscale('log')
xlabel('timestep $\Delta t$')
ylabel('error $L_2$-norm')
legend( (time_stepping_methods[0], time_stepping_methods[1], time_stepping_methods[2], time_stepping_methods[3]), loc='upper left')
title('$v$-velocity')
savefig('figures/TV/v-velocity_L2.pdf')

figure()
p25 = plot(time_steps[:-1], pL2[0,:].T, ls='-', color=lavender, label='p25')
p24 = plot(time_steps[:-1], pL2[1,:].T, ls='-.', color=fuchsia, label='p24')
p23 = plot(time_steps[:-1], pL2[2,:].T, ls='--', color=orange, label='p23')
p21 = plot(time_steps[:-1], pL2[3,:].T, ls=':', color=yelgreen, label='p21')
yscale('log')        
xscale('log')
xlabel('timestep $\Delta t$')
ylabel('error $L_2$-norm')
legend( (time_stepping_methods[0], time_stepping_methods[1], time_stepping_methods[2], time_stepping_methods[3]), loc='upper left')
title('pressure')
savefig('figures/TV/pressure_L2.pdf')

show()         
         
         










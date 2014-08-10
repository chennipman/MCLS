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
import re
from triangle import *


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

# define the testcase
testcase = 'TV'


# choose files based on logfiles
logfilesdir = 'logfiles/' + testcase 

#logfile = '/08-06-13-40-logfile'
#extracomment = 'Run_1' 

#logfile = '/08-07-14-33-logfile'
#extracomment = 'Run_2' 

#logfile = '/08-07-15-16-logfile'
#extracomment = 'Run_3' 

#logfile = '/08-07-15-38-logfile'
#extracomment = 'Run_3_without_pressure' 

logfile = '/08-07-16-32-logfile'
extracomment = 'Run_1_without_pressure' 

namelogfile = logfilesdir + logfile + '.txt'
f = open(namelogfile, 'r')
lines = f.readlines() 

line = lines[2]
line = line[1:-2]
time_stepping_methods.extend(np.array([str(val) for val in line.rstrip('\n').split(',') if val != '']))

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
mLinf = np.zeros((len(time_stepping_methods),len(time_steps)-1))
pLinf = np.zeros((len(time_stepping_methods),len(time_steps)-1))

uL2 = np.zeros((len(time_stepping_methods),len(time_steps)-1))
vL2 = np.zeros((len(time_stepping_methods),len(time_steps)-1))
mL2 = np.zeros((len(time_stepping_methods),len(time_steps)-1))
pL2 = np.zeros((len(time_stepping_methods),len(time_steps)-1))

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
        mLinf[counter,counter2] = max(uLinf[counter,counter2],vLinf[counter,counter2])
        pLinf[counter,counter2] = np.amax(abs(pressure_bench-pressure))
        
        uL2[counter,counter2] = rms(velocity_u_bench-velocity_u)
        vL2[counter,counter2] = rms(velocity_v_bench-velocity_v)
        mL2[counter,counter2] = rms(np.append((velocity_v_bench-velocity_v),(velocity_v_bench-velocity_v)))
        pL2[counter,counter2] = rms(pressure_bench-pressure)


# visualize

rc('text', usetex=True)
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : '11'}
params = {'legend.fontsize': 10}
rcParams.update(params)

rc('font', **font)  # pass in the font dict as kwargs
time_stepping_methods =[w.replace('_', '-') for w in time_stepping_methods]
time_stepping_methods =[w.replace('\'', '') for w in time_stepping_methods]


figure(1, figsize=(10,4))
figuretitle = extracomment.replace('_without_pressure', ', without pressure in the predictor').replace('_', ' ') 

suptitle(figuretitle)
subplot(121)
m25 = plot(time_steps[:-1], mL2[0,:].T, ls='-',  marker ='v', color=lavender, label='m25')
m24 = plot(time_steps[:-1], mL2[1,:].T, ls='-.', marker ='^', color=fuchsia,  label='m24')
m23 = plot(time_steps[:-1], mL2[2,:].T, ls='--', marker ='<', color=orange,   label='m23')
m21 = plot(time_steps[:-1], mL2[3,:].T, ls=':',  marker ='>', color=yelgreen, label='m21')
yscale('log')        
xscale('log')
ylim( (1e-6, 1e-1))
slope_marker((2e-3, 2e-6), (2), size_frac=0.2)
slope_marker((3e-2, 1e-4), (1), size_frac=0.2)
xlabel('timestep $\Delta t$')
ylabel('error $L_2$-norm')
#legend( (time_stepping_methods[0], time_stepping_methods[1], time_stepping_methods[2], time_stepping_methods[3]), loc='best')
title('velocity')


subplot(122)
p25 = plot(time_steps[:-1], pL2[0,:].T, ls='-',  marker ='v', color=lavender, label='p25')
p24 = plot(time_steps[:-1], pL2[1,:].T, ls='-.', marker ='^', color=fuchsia,  label='p24')
p23 = plot(time_steps[:-1], pL2[2,:].T, ls='--', marker ='<', color=orange,   label='p23')
p21 = plot(time_steps[:-1], pL2[3,:].T, ls=':',  marker ='>', color=yelgreen, label='p21')
yscale('log')        
xscale('log')
ylim( (2e-6, 2e-1))
slope_marker((2e-3, 5e-6), (2), size_frac=0.2)
slope_marker((3e-2, 5e-4), (1), size_frac=0.2)
xlabel('timestep $\Delta t$')
ylabel('error $L_2$-norm')
legend( (time_stepping_methods[0], time_stepping_methods[1], time_stepping_methods[2], time_stepping_methods[3]), loc='best')
title('pressure')

filefolder = 'figures/' + testcase
filename = filefolder + logfile +  extracomment + '_L2.pdf'
savefig(filename, bbox_inches='tight')
show()         
         
         










# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 13:06:23 2014

@author: coen

Make and execute several runs for different timesteps needed for a Richardson extrapolation

Output: logfiles are changed in logfiles_for_richardson


"""

import os
import time
import datetime

# define the testcase
testcase = 'TV'

#go to basis directory and make clean (clean only needed if other testcases were run)
os.chdir('../')
#os.system('make clean')

# make timestamp for logfile anc create logfile
st = int(time.time())
ts = datetime.datetime.fromtimestamp(st).strftime('%m-%d-%H-%M')
tsstr = str(ts)
namelogfile = []
namelogfile.append(tsstr)
namelogfile.append("-logfile.txt")
namelogfile = ''.join(namelogfile)
namelogfile = 'logfiles_for_richardson/' + namelogfile

time_stepping_methods = [5,4,3,1]
#time_stepping_methods = [5,4,3]
time_steps = [1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 1e-5]
#time_steps = [1e-1, 5e-2, 1e-4]
#time_steps = [1e-1, 5e-2, 1e-2, 1e-3]

with open(namelogfile, 'w') as file:
    file.write('logfile for richardson extrapolation. The testcase is: ' + testcase +  '. Time of run: ' + tsstr + ' \n')
    file.write('time_stepping_methods: \n')
    file.write(str(time_stepping_methods))
    file.write('\n' + 'time_steps: \n')
    file.write(str(time_steps))
    file.write('\n')
    file.close()



for time_stepping_method in time_stepping_methods:
	for time_step_restriction_global in time_steps:

		print(' time step = ') 
		print( time_step_restriction_global)
		name_parameter_file = 'testcases/' + testcase + '/set_parameters.cpp'           
              
		with open(name_parameter_file, 'r') as file:
		    data = file.readlines()

		#now change the set_parameter file
		string_time = '      time_step_restriction_global			= ' + str(time_step_restriction_global) + ';   \n'
		data[144] = string_time
		string_step_method = '      time_stepping_method 				= ' + str(time_stepping_method) + '; 	// time scheme 1:explicit euler 2: imex, 3: runge-kutta   \n'
		data[153] = string_step_method

		# and write everything back
		with open(name_parameter_file, 'w') as file:
		    file.writelines(data)
		    file.close()

		#run the make file
		cwd = os.getcwd() # get current directory
		makefile_command = 'make case_' + testcase + ' -j2' 
		os.system(makefile_command)
		os.chdir(cwd)

		#execute the executable
		path_executables = 'executable/case_' + testcase
		print cwd
		os.chdir(path_executables)            
		case_direc = os.popen("ls -d */ | tail -n 1").readlines() # most recent made directory
		case_direc = "".join(case_direc)
		case_direc = case_direc.rstrip('\n')
		os.chdir(case_direc)
		os.system("./MCLS </dev/null 1> output.log 2> error.log &")

		# back to the orignal directory
		os.chdir(cwd)
		print 'start tail \n'
		os.system('tail ' + namelogfile)
		print '\n end tail \n' 

		with open(namelogfile, 'a') as file:
		    case_direc = case_direc.rstrip('/')
		    case_direc = path_executables + '/' + case_direc + ' \n'  
		    file.write(case_direc)

		    file.close()


		time.sleep(1)






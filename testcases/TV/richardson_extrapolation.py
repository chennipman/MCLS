import os
import time

#time_stepping_methods = [1,3]
#for time_stepping_method in time_stepping_methods:

#time_steps = [1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3,  5e-4, 2e-4, 1e-4]

time_steps = [1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3]
for time_step_restriction_global in time_steps:
#	print(' time_stepping_method = ')
#	print(time_stepping_method)
	print(' time step = ') 
	print( time_step_restriction_global)
	with open('set_parameters.cpp', 'r') as file:
	#open the file
	    data = file.readlines()

	#now change the set_parameter file
	string_time = '      time_step_restriction_global			= ' + str(time_step_restriction_global) + ';  \n '
	data[147] = string_time
	time_stepping_method = 3
	string_step_method = '      time_stepping_method 				= ' + str(time_stepping_method) + '; 	// time scheme 1:explicit euler 2: imex, 3: runge-kutta   \n'
	data[156] = string_step_method

	# and write everything back
	with open('set_parameters.cpp', 'w') as file:
	    file.writelines(data)
	    file.closed

	#running the make file
	build_dir = "../../"
	cwd = os.getcwd() # get current directory
	os.chdir(build_dir)
	os.system("make")
	os.chdir(cwd)

	case_direc = os.popen("ls -d */ | tail -n 1").readlines()
	case_direc = "".join(case_direc)
	case_direc = case_direc.rstrip('\n')
	os.chdir(case_direc)
	os.system("./MCLS")
	os.chdir(cwd)
	time.sleep(10)
	
	

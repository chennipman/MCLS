# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 18:02:43 2014

@author: coen

Rename a variable in all the parameter files
"""

import os

#go to basis directory 
os.chdir('../')

testcases = ['BB', 'KH', 'LAP', 'LDC', 'RB', 'RT', 'TV', 'undefined_case']

for testcase in testcases:
	name_parameter_file = 'testcases/' + testcase + '/set_parameters_default.cpp'           
	with open(name_parameter_file, 'r') as file:
	    data = file.readlines()
	string_step_method = '      time_stepping_method 				= explict_euler; 	//{none, explicit_euler, imex, runge_kutta, two_pres_solve, two_pres_solve_output};   \n'
	data[153] = string_step_method
		# and write everything back
	with open(name_parameter_file, 'w') as file:
	    file.writelines(data)
	    file.close()
	print 'renamed ' + name_parameter_file


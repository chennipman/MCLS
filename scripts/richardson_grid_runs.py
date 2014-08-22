# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 21:01:47 2014

@author: coen

Richardson Grid Convergence

Several grid sizes are defined. 
The parameter file is changed to the grid size in hand.
The Makefile is run and an executable is created.
The executable is run. 
The location of the run is saved in a logfile

"""


import os
import time
import datetime

# define the case
case = 'BB'

if case == 'RB':
    domain_size_x1 = 1.0
    domain_size_x3 = 2.0
elif case == 'RT':
    domain_size_x1 = 1.0
    domain_size_x3 = 4.0
elif case == 'BB':
    domain_size_x1 = 6.0
    domain_size_x3 = 0.1
elif case == 'KH':
    domain_size_x1 = 0.1
    domain_size_x3 = 6.0
else:
    print 'undefined testcase'
    exit


#go to basis directory and make clean 
#(clean only needed if other testcases were run)
os.chdir('../')
cwd = os.getcwd()
print cwd
#os.system('make clean')

# make timestamp for logfile and create logfile
st = int(time.time())
ts = datetime.datetime.fromtimestamp(st).strftime('%m-%d-%H-%M')
tsstr = str(ts)
namelogfile = []
namelogfile.append(tsstr)
namelogfile.append("-logfile.txt")
logfilesdir = 'logfiles/'
command = "mkdir -p " + logfilesdir + case
os.system(command)
namelogfile =case + '/'+ ''.join(namelogfile)
namelogfile = 'logfiles/' + namelogfile

number_primary_cells_i_all = [10] #,100]

if case == 'BB':
	number_primary_cells_i_all = [10, 20, 40] #,100]


# write the log file
with open(namelogfile, 'w') as file:
    file.write('logfile for richardson extrapolation. The case is: ' + case +  '. Time of run: ' + tsstr + ' \n')
    file.write('number_primary_cells_i_all: \n')
    file.write(str(number_primary_cells_i_all) + '\n')
    file.close()

# loop over all the different gridsizes
for number_primary_cells_i in number_primary_cells_i_all:

    # calculate all the parameters depending on the grid size
    number_primary_cells_k = int(number_primary_cells_i*domain_size_x3/domain_size_x1)
    number_primary_cells_j = 1;
    domain_size_x2 = domain_size_x1/number_primary_cells_i

    # update the user about the progress
    print('number_primary_cells_i = ')
    print( number_primary_cells_i)
    name_parameter_file = 'testcases/' + case + '/set_parameters.cpp'


    # open the parameter file and obtain all the data in it
    with open(name_parameter_file, 'r') as file:
        data = file.readlines()

    # change the parameter file to the paraameter set before
    # it is important that the structure of the parameter file is not changed
    # because this script depends on the line numbers in the parameter file
    string_domain_size_x2 = '      domain_size_x2='+str(domain_size_x2) +';  \n'
    data[141] = string_domain_size_x2
    string_number_primary_cells_i = '      number_primary_cells_i='+str(number_primary_cells_i)+';  \n'
    data[171] = string_number_primary_cells_i
    string_number_primary_cells_j = '      number_primary_cells_j='+str(number_primary_cells_j)+';  \n'
    data[172] = string_number_primary_cells_j
    string_number_primary_cells_k = '      number_primary_cells_k='+str(number_primary_cells_k)+';  \n'
    data[173] = string_number_primary_cells_k

    # and write everything back to the parameter file
    with open(name_parameter_file, 'w') as file:
        file.writelines(data)
        file.close()
        print 'finished writing parameter file'

    #run the makefile
    cwd = os.getcwd() # get current directory
    makefile_command = 'make case_' + case + ' -j2'
    os.system(makefile_command)

    #execute the executable
    path_executables = 'executable/case_' + case
    os.chdir(path_executables)
    case_direc = os.popen("ls -d */ | tail -n 1").readlines() # most recent made directory
    case_direc = "".join(case_direc)
    case_direc = case_direc.rstrip('\n')
    os.chdir(case_direc)
    os.system("./MCLS </dev/null 1> output.log 2> error.log &")

    # back to the orignal working directory
    os.chdir(cwd)
    print 'start tail \n'
    os.system('tail ' + namelogfile)
    print '\n end tail \n'

    # add the run to the logfile
    with open(namelogfile, 'a') as file:
        case_direc = case_direc.rstrip('/')
        case_direc = path_executables + '/' + case_direc + ' \n'
        file.write(case_direc)
        file.close()

    # extra second added to prevent creating runs in the same second
    # and consequently the same directory
    time.sleep(1)



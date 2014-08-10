# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 21:01:47 2014

@author: coen

Run the grid convergence
"""


import os
import time
import datetime

# define the testcase
testcase = 'RB'

if testcase == 'RB':
    domain_size_x1 = 1.0
    domain_size_x3 = 2.0
elif testcase == 'RT':
    domain_size_x1 = 1.0
    domain_size_x3 = 4.0
elif testcase == 'BB':
    domain_size_x1 = 0.1
    domain_size_x3 = 6.0
elif testcase == 'KH':
    domain_size_x1 = 0.1
    domain_size_x3 = 6.0
else:
    print 'undefined testcase'
    exit


#go to basis directory and make clean 
#(clean only needed if other testcases were run)
os.chdir('../')
#os.system('make clean')

# make timestamp for logfile anc create logfile
st = int(time.time())
ts = datetime.datetime.fromtimestamp(st).strftime('%m-%d-%H-%M')
tsstr = str(ts)
namelogfile = []
namelogfile.append(tsstr)
namelogfile.append("-logfile.txt")
logfilesdir = 'logfiles/'
command = "mkdir -p " + logfilesdir + testcase
os.system(command)
namelogfile =testcase + '/'+ ''.join(namelogfile)
namelogfile = 'logfiles/' + namelogfile


number_primary_cells_i_all = [50] #,100]


with open(namelogfile, 'w') as file:
    file.write('logfile for richardson extrapolation. The testcase is: ' + testcase +  '. Time of run: ' + tsstr + ' \n')
    file.write('number_primary_cells_i_all: \n')
    file.write(str(number_primary_cells_i_all) + '\n')
#    file.write('\n' + 'time_steps: \n')
    file.close()


for number_primary_cells_i in number_primary_cells_i_all:
    number_primary_cells_k = int(number_primary_cells_i*domain_size_x3/domain_size_x1)
    number_primary_cells_j = 1;
    domain_size_x2 = domain_size_x1/number_primary_cells_i


    print('number_primary_cells_i = ')
    print( number_primary_cells_i)
    name_parameter_file = 'testcases/' + testcase + '/set_parameters.cpp'



    with open(name_parameter_file, 'r') as file:
        data = file.readlines()

    #now change the set_parameter file
    string_domain_size_x2 = '      domain_size_x2='+str(domain_size_x2) +';  \n'
    data[141] = string_domain_size_x2
    string_number_primary_cells_i = '      number_primary_cells_i='+str(number_primary_cells_i)+';  \n'
    data[171] = string_number_primary_cells_i
    string_number_primary_cells_j = '      number_primary_cells_j='+str(number_primary_cells_j)+';  \n'
    data[172] = string_number_primary_cells_j
    string_number_primary_cells_k = '      number_primary_cells_k='+str(number_primary_cells_k)+';  \n'
    data[173] = string_number_primary_cells_k


    # and write everything back
    with open(name_parameter_file, 'w') as file:
        file.writelines(data)
        file.close()
        print 'finished writing parameter file'

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



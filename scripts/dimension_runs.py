# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 01:50:31 2014

Make and execute several runs for different directions. 

Output: logfiles are changed in logfiles


"""

import os
import time
import datetime

# define the testcase
testcase = 'RB'

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
namelogfile = 'logfiles/' + namelogfile

main_directions = ['x','y','z']
middle_directions = ['x','y','z']

main_directions = ['y']
middle_directions = ['z']

with open(namelogfile, 'w') as file:
    file.write('logfile for dimension check. The testcase is: ' + testcase +  '. Time of run: ' + tsstr + ' \n')
    file.write('main_directions: \n')
    file.write(str(main_directions))
    file.write('\n' + 'middle_directions: \n')
    file.write(str(middle_directions))
    file.write('\n')
    file.close()

name_parameter_file = 'testcases/' + testcase + '/set_parameters.cpp'  
with open(name_parameter_file, 'r') as file:
    data = file.readlines()

for main_direc in main_directions: 
    # select the proper boundary conditions
    command = 'cp ' + 'testcases/' + testcase + '/set_boundary_conditions_' + main_direc + '.cpp '  + 'testcases/' + testcase + '/set_boundary_conditions.cpp'    
    print command
    print os.getcwd()
    os.system(command)
    
    for middle_direc in middle_directions: 
        not_execute = 0;
    # do the the parameter file
    # the main directions is set
        if main_direc == 'x':
            data[135] = '      domain_size_x1=2.0;   \n'
            data[163] = '      number_primary_cells_i =100;   \n'
            data[217] = '      gravity.u1 = -0.98;   \n'
            
            if middle_direc == 'y':
                data[136] = '      domain_size_x2=1.0;   \n'
                data[164] = '      number_primary_cells_j =50;   \n'
                data[218] = '      gravity.u2 = 0.0;   \n'
                
                #last direction = z
                data[137] = '      domain_size_x3=1.0/50.0;   \n'
                data[165] = '      number_primary_cells_k =1;   \n'
                data[219] = '      gravity.u3 = 0.0;   \n'

                data[194] ='      the_bubbles[0].center_location.x1  				= 0.5;   \n'               
                data[195] ='      the_bubbles[0].center_location.x2  				= 0.5;   \n'               
                data[196] ='      the_bubbles[0].center_location.x3  				= 0.5*domain_size_x3;   \n' 
                
            elif middle_direc == 'z':
                data[137] = '      domain_size_x3=1.0;   \n'
                data[165] = '      number_primary_cells_k =50;   \n'
                data[219] = '      gravity.u3 = 0.0;   \n'
                
                #last direction = y
                data[136] = '      domain_size_x2=1.0/50.0;   \n'
                data[164] = '      number_primary_cells_j =1;   \n'
                data[218] = '      gravity.u2 = 0.0;   \n'

                data[194] ='      the_bubbles[0].center_location.x1  				= 0.5;   \n'               
                data[195] ='      the_bubbles[0].center_location.x2  				= 0.5*domain_size_x2;   \n'               
                data[196] ='      the_bubbles[0].center_location.x3  				= 0.5;   \n'  
            else:
                # do nothing
                print 'restricted situation x'  
                not_execute = 1
                
        elif main_direc == 'y':
            data[136] = '      domain_size_x2=2.0;   \n'
            data[164] = '      number_primary_cells_j =100;   \n'
            data[218] = '      gravity.u2 = -0.98;   \n'
            
            if middle_direc == 'x':
                data[135] = '      domain_size_x1=1.0;   \n'
                data[163] = '      number_primary_cells_i =50;   \n'
                data[217] = '      gravity.u1 = 0.0;   \n'

                #last direction = z
                data[137] = '      domain_size_x3=1.0/50.0;   \n'
                data[165] = '      number_primary_cells_k =1;   \n'
                data[219] = '      gravity.u3 = 0.0;   \n'

                data[194] ='      the_bubbles[0].center_location.x1  				= 0.5;   \n'               
                data[195] ='      the_bubbles[0].center_location.x2  				= 0.5;   \n'               
                data[196] ='      the_bubbles[0].center_location.x3  				= 0.5*domain_size_x3;   \n' 
                 
            elif middle_direc == 'z':
                data[137] = '      domain_size_x3=1.0;   \n'
                data[165] = '      number_primary_cells_k =50;   \n'
                data[219] = '      gravity.u3 = 0.0;   \n'
                
                #last direction = x 
                data[135] = '      domain_size_x1=1.0/50;   \n'
                data[163] = '      number_primary_cells_i =1;   \n'
                data[217] = '      gravity.u1 = 0.0;   \n'

                data[194] ='      the_bubbles[0].center_location.x1  				= 0.5*domain_size_x1;   \n'               
                data[195] ='      the_bubbles[0].center_location.x2  				= 0.5;   \n'               
                data[196] ='      the_bubbles[0].center_location.x3  				= 0.5;   \n'  
            else:
                # do nothing
                print 'restricted situation y'  
                not_execute = 1
                
        elif main_direc == 'z':
            data[137] = '      domain_size_x3=2.0;   \n'
            data[165] = '      number_primary_cells_k =100;   \n'
            data[219] = '      gravity.u3 = -0.98;   \n'

            if middle_direc == 'x':
                data[135] = '      domain_size_x1=1.0;   \n'
                data[163] = '      number_primary_cells_i =50;   \n'
                data[217] = '      gravity.u1 = 0.0;   \n'

                #last direction = y
                data[136] = '      domain_size_x2=1.0/50.0;   \n'
                data[164] = '      number_primary_cells_j =1;   \n'
                data[218] = '      gravity.u2 = 0.0;   \n'

                data[194] ='      the_bubbles[0].center_location.x1  				= 0.5;   \n'               
                data[195] ='      the_bubbles[0].center_location.x2  				= 0.5*domain_size_x2;   \n'               
                data[196] ='      the_bubbles[0].center_location.x3  				= 0.5;   \n'  

            elif middle_direc == 'y':
                data[136] = '      domain_size_x2=1.0;   \n'
                data[164] = '      number_primary_cells_j =50;   \n'
                data[218] = '      gravity.u2 = 0.0;   \n'
 
                #last direction = x 
                data[135] = '      domain_size_x1=1.0/50;   \n'
                data[163] = '      number_primary_cells_i =1;   \n'
                data[217] = '      gravity.u1 = 0.0;   \n'

                data[194] ='      the_bubbles[0].center_location.x1  				= 0.5*domain_size_x1;   \n'               
                data[195] ='      the_bubbles[0].center_location.x2  				= 0.5;   \n'               
                data[196] ='      the_bubbles[0].center_location.x3  				= 0.5;   \n'  
            else:
                # do nothing
                print 'restricted situation z'  
                not_execute = 1

        if not_execute != 1:
            
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
            
            # write to the logfile
            with open(namelogfile, 'a') as file:
                case_direc = case_direc.rstrip('/')
                case_direc = path_executables + '/' + case_direc + ' \n'  
                file.write(case_direc)
                file.close()
                
            time.sleep(1)


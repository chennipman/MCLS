#include "../headers/array.h"
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/********************************************************************************/
/*  Function to write a restart file containing the main unknowns               */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The restart file is a binary file with the complete solution. 		*/
/* The solution is written in the following order				*/
/* 										*/
/*	  level_set_new								*/	
/*	  volume_of_fluid 							*/		
/*	  u_1_velocity_new 							*/		
/*	  u_2_velocity_new 							*/		
/*	  u_3_velocity_new							*/		
/*	  pressure								*/		
/********************************************************************************/
EXPORT void write_restart_file(
	  double time_of_restart,			// time for which restart file is written
	  Array3<double> level_set, 			// level set field 
							// mass conserving
	  Array3<double> volume_of_fluid, 		// volume of fluid field
	  Array3<double> u_1_velocity, 			// velocity field  x1 direction
	  Array3<double> u_2_velocity, 			// velocity field  x2 direction
	  Array3<double> u_3_velocity,			// velocity field  x3 direction
	  Array3<double> pressure,			// pressure field
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  restart_parameters my_restart_parameters	// all parameters for reading/writing restart files
	   )
     {     
     
     
     
	  ofstream OutputFile;			        // the stream is set to the output file
	  string OutputFile_name;			// the name of the output file

	  /* initialize the output file for writing */
	  
	  OutputFile_name=my_restart_parameters.name_restart_file_to_write;
	  OutputFile.open(OutputFile_name.c_str(), ios::out | ios:: binary);
	  
	  
	  /* write the time at which the restart file is written */
	  cerr<<"time of restart is "<< time_of_restart<<" \n";
	  OutputFile.write((char *) &time_of_restart, sizeof(double));

    /* write the level-set, volume of fluid, velocity and pressure */

          level_set.write( OutputFile );
          volume_of_fluid.write( OutputFile );
          u_1_velocity.write( OutputFile );
          u_2_velocity.write( OutputFile );
          u_3_velocity.write( OutputFile );
          pressure.write( OutputFile );

          OutputFile.close();
}

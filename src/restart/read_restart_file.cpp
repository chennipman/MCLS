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
/*  Function to read a restart file containing the main unknowns                */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The restart file is a binary file with the complete solution. 		*/
/* The solution is read in the following order					*/
/* 										*/
/*	  level_set								*/	
/*	  volume_of_fluid 							*/		
/*	  u_1_velocity 								*/		
/*	  u_2_velocity 								*/		
/*	  u_3_velocity								*/		
/*	  pressure								*/		
/********************************************************************************/
EXPORT void read_restart_file(
	  double &time_of_restart,			// time for which restart file is available
	  Array3<double> level_set, 			// level set field at new time level
							// mass conserving
	  Array3<double> volume_of_fluid, 		// volume of fluid field
	  Array3<double> u_1_velocity, 			// velocity field x1 direction
	  Array3<double> u_2_velocity, 			// velocity field x2 direction
	  Array3<double> u_3_velocity,			// velocity field x3 direction
	  Array3<double> pressure,			// pressure field
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  restart_parameters my_restart_parameters	// all parameters for reading/writing restart files
	   )
     {     
     
     
     
	    ifstream InputFile;				// the stream is set to the output file
	    string InputFile_name;			// the name of the output file
// 	    int final_index_i;				// final index first dimension 
// 	    int final_index_j;				// final index second dimension 
// 	    int final_index_k;				// final index third dimension 
// 	    int i_index, j_index;  			// local variables for loop indexing
// 	    double gnoef;
     

	  
	  /* initialize the output file for writing */
	  
	  InputFile_name=my_restart_parameters.name_restart_file_to_read;
	  InputFile.open(InputFile_name.c_str(), ios::in | ios:: binary);

	  /* read the time at which the restart file is available */
	  
	  InputFile.read((char *) &time_of_restart, sizeof(time_of_restart));
	  
	  cerr<<"time of restart is "<< time_of_restart<<" \n";

    /* read the level-set, volume of fluid, velocity and pressure */

          level_set.read( InputFile );
          volume_of_fluid.read( InputFile );
          u_1_velocity.read( InputFile );
          u_2_velocity.read( InputFile );
          u_3_velocity.read( InputFile );
          pressure.read( InputFile );

          InputFile.close();
}

#include "../headers/array.h"
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
/********************************************************************************/
/********************************************************************************/
/*  Function to interpolate velocity u1 to the cell center of primary cell      */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes							http://www.kustzeilers.nl/cms/driehoek/homehttp://www.kustzeilers.nl/cms/driehoek/home		*/
/*  										*/
/*  										*/
/*  										*/
/*  										*/
/********************************************************************************/
EXPORT void interpolate_velocity_u1_vertex(
	  Array3<double> u_1_velocity_new, 			// velocity field at new time level x1 direction
	  Array3<double> u_1_velocity_vertex,		// velocity in cell vertex, x1 component
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k			// number of primary (pressure) cells in x3 direction
			    )
  {
	int i_index, j_index, k_index; 			// local variables for loop indexing
	    
	    for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
	    {
		for(j_index=1;j_index<number_primary_cells_j+2;j_index++)
		{
		    for(k_index=1;k_index<number_primary_cells_k+2;k_index++)
		    {
		      u_1_velocity_vertex[i_index][j_index-1][k_index-1]=0.25*(
			      u_1_velocity_new[i_index][j_index-1][k_index  ] +
			      u_1_velocity_new[i_index][j_index-1][k_index-1] +
			      u_1_velocity_new[i_index][j_index  ][k_index  ] +
			      u_1_velocity_new[i_index][j_index  ][k_index-1] 
								       );
		    }
		}
	    }
  }

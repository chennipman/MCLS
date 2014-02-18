#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
/********************************************************************************/
/********************************************************************************/
/*  Function to write the cell coordinates to file in vtk format                */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/*  										*/
/*  										*/
/*  										*/
/*  										*/
/********************************************************************************/
EXPORT void  write_coordinates_vtk( 
	  ofstream& output_stream, 		// stream connected to output file
	  int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
	  double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	  double mesh_width_x3			// grid spacing in x3 direction (uniform)
       )
  {
	double x1_coordinate;			// x1 coordinate of cell center
	double x2_coordinate;			// x2 coordinate of cell center
	double x3_coordinate;			// x3 coordinate of cell center
	
	int i_index, j_index, k_index; 		// local variables for loop indexing
	int full_row=7;				// a full row of 8 numbers has been written to file

	/* write x1 coordinates */
	
	output_stream << " X_COORDINATES "<< number_primary_cells_i+1 << " float \n";
	for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
	{
	    full_row--;
	    x1_coordinate=mesh_width_x1*i_index;
	    output_stream << x1_coordinate << " ";
	    if(!full_row)
	    {
		 output_stream <<"\n";
		 full_row=7;
	    }
	}
	if(full_row!=7)
	{
	    full_row=7;
	    output_stream << "\n";
	}

	/* write x2 coordinates */
	
	output_stream << " Y_COORDINATES "<< number_primary_cells_j+1 << " float \n";
	for(j_index=0;j_index<number_primary_cells_j+1;j_index++)
	{
	    full_row--;
	    x2_coordinate=mesh_width_x2*j_index;
	    output_stream << x2_coordinate << " ";
	    if(!full_row)
	    {
		 output_stream <<"\n";
		 full_row=7;
	    }
	}
	if(full_row!=7)
	{
	    full_row=7;
	    output_stream << "\n";
	}

	/* write x3 coordinates */
	
	output_stream << " Z_COORDINATES "<< number_primary_cells_k+1 << " float \n";
	for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
	{
	    full_row--;
	    x3_coordinate=mesh_width_x3*k_index;
	    output_stream << x3_coordinate << " ";
	    if(!full_row)
	    {
		 output_stream <<"\n";
		 full_row=7;
	    }
	}
	if(full_row!=7)
	{
	    full_row=7;
	    output_stream << "\n";
	}
		  
  }

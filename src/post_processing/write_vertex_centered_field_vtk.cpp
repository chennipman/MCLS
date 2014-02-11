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
/*  Function to write cell centered scalar field to file                        */
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

      void  write_vertex_centered_field_vtk( 
	ofstream& output_stream, 			// stream connected to output file
	string scalar_name,				// name of the scalar field to write to file 
	string look_up_table_name,			// name of the look-up table
	double ***cell_centered_field, 			// cell centered scalar field
	int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	int number_primary_cells_k			// number of primary (pressure) cells in x3 direction
	    )
  {
	
      int i_index, j_index, k_index; 		// local variables for loop indexing
      int full_row=7;				// a full row of 8 numbers has been written to file

      /* write short header */
      
      output_stream << "SCALARS " << scalar_name << " float \n";
      output_stream << "LOOKUP_TABLE " << look_up_table_name << "\n";
	    
      for(k_index=1;k_index<number_primary_cells_k+2;k_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+2;j_index++)
	  {
	      for(i_index=1;i_index<number_primary_cells_i+2;i_index++)
	      {
		   full_row--;
		   output_stream << 0.125*(
		     cell_centered_field[i_index  ][j_index ][k_index   ]+
		     cell_centered_field[i_index-1][j_index ][k_index   ]+
		     cell_centered_field[i_index  ][j_index-1][k_index  ]+
		     cell_centered_field[i_index-1][j_index-1][k_index  ]+
		     cell_centered_field[i_index  ][j_index  ][k_index-1]+
		     cell_centered_field[i_index-1][j_index  ][k_index-1]+
		     cell_centered_field[i_index  ][j_index-1][k_index-1]+
		     cell_centered_field[i_index-1][j_index-1][k_index-1])
		     << " ";
		   if(!full_row)
		   {
			output_stream<<"\n";
			full_row=7;
		   }
	      }
	  }
      }
      
      if(full_row!=7)
      {
	  output_stream<<"\n";
      }
  }
      

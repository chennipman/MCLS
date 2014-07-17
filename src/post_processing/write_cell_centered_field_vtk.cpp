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

EXPORT void  write_cell_centered_field_vtk( 
	std::ofstream& output_stream, 			// stream connected to output file
	std::string scalar_name,			// name of the scalar field to write to file 
	std::string look_up_table_name,			// name of the look-up table
	Array3<double> cell_centered_field, 		// cell centered scalar field
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
	    
      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	      {
		   full_row--;
		   double value=0.0;
		   if(fabs(cell_centered_field[i_index][j_index][k_index])>1e-18)
		   {
		      value=	  cell_centered_field[i_index][j_index][k_index];
		   }
		   // if(cell_centered_field[i_index][j_index][k_index]<0)
		   // {
		   // 	   output_stream << cell_centered_field[i_index][j_index][k_index] << "     ";
		   // }
		   // else
		   // {
		   // 	   output_stream << " "<< cell_centered_field[i_index][j_index][k_index] << "     ";
		   // }
		   if(value<0)
		   {
		   	   output_stream << value << "     ";
		   }
		   else
		   {
		   	   output_stream << " "<< value << "     ";
		   }
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
      

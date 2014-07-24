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
/*  Programmer	: Coen Hennipman						*/
/*  Date	: 11-07-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/*  										*/
/*  										*/
/*  										*/
/*  										*/
/********************************************************************************/

EXPORT void  write_face_field_pure( 
	std::ofstream& output_stream, 			// stream connected to output file
	Array3<double> face_field, 			// face field, 1 component
	int first_dimension,				// number of elements in first dimension
	int second_dimension,				// number of elements in second dimension
	int third_dimension				// number of elements in third dimension
	    )
{
	
      int i_index, j_index, k_index; 		// local variables for loop indexing
      int full_row=0;				// a full row of 3*3=9 numbers has been written to file
  for(i_index=0;i_index<first_dimension;i_index++)
  {
     for(j_index=0;j_index<second_dimension;j_index++)
      {
	  for(k_index=0;k_index<third_dimension;k_index++)
	  {
		   full_row++;
		   output_stream << face_field[i_index][j_index][k_index] << " ";
 		   if(full_row>=7) //full row 
 		   {
			output_stream<<"\n";
			full_row=0;
		   }
	      }
      }
  }
  if (full_row != 0 ) // if the last line was not ended with a newline, do so now
    {
	output_stream<<"\n";
    }
}


#include "../headers/array.h"
/********************************************************************************/
/********************************************************************************/
/*  Function to copy a cell centered field including virtual cells              */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* when we switch to more efficient definition for the multi-dimensional arrays */
/* this function can be simplified of even discarded				*/
/********************************************************************************/
      void copy_general_field( 
	    Array3<double> source_field, 		// original field
	    Array3<double> target_field,		// copy of the original field
	    int start_index_i,			// start index first dimension 
	    int final_index_i,			// final index first dimension 
	    int start_index_j,			// start index second dimension 
	    int final_index_j,			// final index second dimension 
	    int start_index_k,			// start index third dimension 
	    int final_index_k			// final index third dimension 
	   )
     {
     int i_index, j_index, k_index;  		// local variables for loop indexing
     
     
          for(i_index=start_index_i;i_index<final_index_i+1;i_index++)
          {
		for(j_index=start_index_j;j_index<final_index_j+1;j_index++)
           	{
           	      for(k_index=start_index_k;k_index<final_index_k+1;k_index++)
           	      {
			    target_field[i_index][j_index][k_index]=
				  source_field[i_index][j_index][k_index];
           	      }
           	}
          }
      }		  
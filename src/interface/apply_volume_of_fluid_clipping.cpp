#include "../headers/array.h"

#include<cstdlib>
#include<iostream>
#include<algorithm>
/********************************************************************************/
/********************************************************************************/
/*  Function to clip the value of the volume of fluid function     		*/
/*  method. 										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* Due to different reasons the volume of fluid field can reach unphysical      */
/* values, outside the allowed interval [0,1]. This can be corrected without    */
/* jeopardizing the mass conservation or by simply clipping the values off      */
/* This functions simply clips off the excess values and monitors the number    */
/* of cells that need correcting.                                               */
/********************************************************************************/
      int apply_volume_of_fluid_clipping(				
		 Array3<double> volume_of_fluid, 		// volume of fluid field		
		 int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
		 int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
		 int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
		 double volume_of_fluid_tolerance	// tolerance for volume of fluid value   
		
      )
      {
      int i_index, j_index, k_index;  		// local variables for loop indexing

      int number_of_clipped_cells=0;		// the number of cells where clipping was
						// necessary to bring the volume of fluid
						// field in the interval [0,1]
     
       
      for(i_index=0;i_index<number_primary_cells_i+2;i_index++)
      {
	  for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
	  {
	      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
	      {
		    if(volume_of_fluid[i_index][j_index][k_index]>
			( 1.0+volume_of_fluid_tolerance )||
			 volume_of_fluid[i_index][j_index][k_index]< 
			 (-1.0*volume_of_fluid_tolerance))
		    {
			volume_of_fluid[i_index][j_index][k_index]=std::max(
			  std::min(volume_of_fluid[i_index][j_index][k_index], 1.0), 0.0);
			number_of_clipped_cells++;
		    }
	      }
	  }
      }
      return number_of_clipped_cells;
}      

#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to check if the volume of fluid field is completely within the     */
/*  proper bounds								*/
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/********************************************************************************/

      int check_volume_of_fluid(
		 double ***volume_of_fluid,		// volume of fluid field field
		 int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
		 int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
		 int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
		 double volume_of_fluid_tolerance	// tolerance for volume of fluid value   
	  
	  )
      
{
     int i_index, j_index, k_index;  // local variables for loop indexing
     int number_of_error_cells=0;    // Number of cells where the volume of fluid fiels is outside the
				     // allowed interval [0,1] (with proper tolerances)
				     
     double minimum_value_vof=10.0;	// mimimum value of the volume of fluid field
     double maximum_value_vof=-10.0;	// maximum value of the volume of fluid field
     
     for(i_index=0;i_index<number_primary_cells_i+2;i_index++)
      {
	  for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
	  {
	      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
	      {
		    maximum_value_vof=
			  std::max(maximum_value_vof,
				    volume_of_fluid[i_index][j_index][k_index]);
		    minimum_value_vof=
			  std::min(minimum_value_vof,
				    volume_of_fluid[i_index][j_index][k_index]);
		    if(volume_of_fluid[i_index][j_index][k_index]>
			( 1.0+volume_of_fluid_tolerance )||
			 volume_of_fluid[i_index][j_index][k_index]< 
			 (-1.0*volume_of_fluid_tolerance))
		    {
			number_of_error_cells++;
		    }
	      }
	  }
      }
      
      if(number_of_error_cells)
      {
      /* cells have been detected that have incorrect values for the volume of fluid */
	    std::cerr<<"ERROR:some of the volume of fluid values are out of bounds \n";
	    std::cerr<<"minimum value "<<minimum_value_vof ;
	    std::cerr<<" maximum value "<<maximum_value_vof <<"\n";
	    std::cerr<<" in function check_volume_of_fluid line 59 \n";
      
      }
      return number_of_error_cells;
}      

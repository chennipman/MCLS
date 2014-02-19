#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the volume of fluid field                               */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/********************************************************************************/
EXPORT int compute_volume_of_fluid(					
	Array3<double> level_set, 			// level set field field
	Array3<double> d_level_set_d_x1, 		// first partial derivative of level-set with 
						// respect to x1, central approximation
	Array3<double> d_level_set_d_x2, 		// first partial derivative of level-set with 
						// respect to x2, central approximation 
	Array3<double> d_level_set_d_x3, 		// first partial derivative of level-set with 
						// respect to x3, central approximation
	Array3<double> volume_of_fluid,		// volume of fluid field
	int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
	int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
	int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
	double lower_bound_derivatives    	// lower bound for the first partial derivatives
						// to consider it a limiting case of vanishing
						// partial derivatives
	  )
     {
      int i_index, j_index, k_index;  		// local variables for loop indexing
      int no_conversion=0;			// =0, successful conversion of level-set field
						// to volume of fluid field. Otherwise, the 
						// conversion failed.


		for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			 for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			 {
			   
			   /* function level_set_2_vof shoud return '0' */
			   /* for successful conversion */

			    if(level_set_2_vof(level_set[i_index][j_index][k_index],
					d_level_set_d_x1[i_index][j_index][k_index],
					  d_level_set_d_x2[i_index][j_index][k_index],
					    d_level_set_d_x3[i_index][j_index][k_index],
					      volume_of_fluid[i_index][j_index][k_index],
					      lower_bound_derivatives))
			    {
			      /* something went wrong in the conversion from level-set to */
			      /* volume of fluid */
			      
				no_conversion=1;
				std::cerr << "**************************************************** \n";
				std::cerr << "ERROR \n";
				std::cerr << "level-set to vof conversion failed in call from compute_volume_of_fluid \n";
				std::cerr << "for cell index "<< i_index<<" "<<j_index<<" "<<k_index<< "\n";
				std::cerr << " in function compute_volume_of_fluid line 66 \n";
				std::cerr << "**************************************************** \n";
				break;
			    }
			 }
			 if(no_conversion)
			 {
			      break;
			 }
		    }
		    if(no_conversion)
		    {
			break;
		    }	      
		}
	if(!no_conversion)
	{
	  /* if the conversion was successful then apply the boundary conditions */
	  /* to fill the virtual cells 						 */
	  
	  field_neumann_boundary(volume_of_fluid, 
				 number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k);
	}
	return (no_conversion);
	
     }					    

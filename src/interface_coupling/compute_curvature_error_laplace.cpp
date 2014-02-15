#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the error in the curvature for the laplace test case    */
/*  proper bounds								*/
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/********************************************************************************/
double compute_curvature_error_laplace(
      Array3<double> curvature,				// interface curvature 
      Array3<double> volume_of_fluid,			// volume of fluid field
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3				// grid spacing in x3 direction (uniform)
       )
      {
      double exact_curvature;				// exact interface curvature, for Laplace testcase
      double distance_to_bubble_center;			// distance of cell center to the center
							// of the bubble
      double maximum_relative_curvature_error=0;	// maximum relative error in the curvature:
							//    (K_h-K)/K
      double x1_coordinate_cell_center;			// x1 coordinate of the cell center
      double x2_coordinate_cell_center;			// x2 coordinate of the cell center
      double x3_coordinate_cell_center;			// x3 coordinate of the cell center

		/* for the moment this is set to 0 */
		/* should ideally use input value  */
		
      double x1_coordinate_bubble_center=0;		// x1 coordinate of the bubble center
      double x2_coordinate_bubble_center=0;		// x2 coordinate of the bubble center
      double x3_coordinate_bubble_center=0;		// x3 coordinate of the bubble center
      double volume_of_fluid_lower_bound=0.025;		// lower bound for volume of fluid for
							// which the curvature is of any interest
      double volume_of_fluid_upper_bound=0.025;		// lower bound for volume of fluid for
							// which the curvature is of any interest
      int i_index, j_index, k_index;  // local variables for loop indexing
      
      
      
     for(i_index=2;i_index<number_primary_cells_i;i_index++)
      {
	
	  x1_coordinate_cell_center= (i_index-0.5)*mesh_width_x1;
	  
	  for(j_index=2;j_index<number_primary_cells_j;j_index++)
	  {
	      
	      x2_coordinate_cell_center= (j_index-0.5)*mesh_width_x1;
	  
	      for(k_index=2;k_index<number_primary_cells_k;k_index++)
	      {
		  x3_coordinate_cell_center= (k_index-0.5)*mesh_width_x1;
		  distance_to_bubble_center=sqrt(
		    (x1_coordinate_cell_center-x1_coordinate_bubble_center)*
		      (x1_coordinate_cell_center-x1_coordinate_bubble_center)+
			(x2_coordinate_cell_center-x2_coordinate_bubble_center)*
			  (x2_coordinate_cell_center-x2_coordinate_bubble_center)+
			    (x3_coordinate_cell_center-x3_coordinate_bubble_center)*
			      (x3_coordinate_cell_center-x3_coordinate_bubble_center)
					    );
		  /* in 3D the curvature is 2/R */
		   exact_curvature=2.0/distance_to_bubble_center;
		   if(volume_of_fluid[i_index][j_index][k_index]> volume_of_fluid_lower_bound &&
			  volume_of_fluid[i_index][j_index][k_index]< volume_of_fluid_upper_bound)
		   {
			maximum_relative_curvature_error=std::max(maximum_relative_curvature_error, 
					       fabs(-1.0*curvature[i_index][j_index][k_index]-
						    exact_curvature)/exact_curvature);
		   }
									
		      
	      }
	  }
      }
      
      return maximum_relative_curvature_error;
 }
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to  update the volume of fluid field due to fluxes in x1          	*/
/*  direction										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
//
      void update_volume_of_fluid_x1_flux(				
	 double ***level_set, 		//level-set field			
	 double ***d_level_set_d_x1,		// first partial derivative of
						// the level-set field wrt x1
						// second order central approximation				
	 double ***d_level_set_d_x2,		// first partial derivative of
						// the level-set field wrt x2
						// second order central approximation 				
	 double ***d_level_set_d_x3,		// first partial derivative of
						// the level-set field wrt x3
						// second order central approximation
	 double ***volume_of_fluid, 		// volume of fluid field
	 double ***vof_after_x1_update,	// intermediate volume of fluid field, 
						// after update in x1 direction
	 double ***u_1_velocity_new, 		// velocity field at new time level x1 direction
	 int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction 
	 int number_primary_cells_j,		// number of primary (pressure) cells in x1 direction 
	 int number_primary_cells_k,		// number of primary (pressure) cells in x1 direction
	 double actual_time_step_level_set,	// time step used for level-set advection
						// computed from all stability restrictions and 
						// possibly subscycling
	 double mesh_width_x1,		// grid spacing in x1 direction (uniform)
	 double lower_bound_derivatives	// lower bound for the first partial derivatives
						// to consider it a limiting case of vanishing
						// partial derivatives
      
  )
     {
      /* function definitions */
      
     void compute_vof_flux_x1(		// compute volume of fluid	
	double ***level_set, 			// fluxes in x1 direction
        double ***u_1_velocity_new, 			
	double ***d_level_set_d_x1,			
	double ***d_level_set_d_x2,			
	double ***d_level_set_d_x3,			
	double ***flux_x1,				
	int number_primary_cells_i,			      
	int number_primary_cells_j,			      
	int number_primary_cells_k,			
	double actual_time_step_level_set,		
	double mesh_width_x1,				        
	double lower_bound_derivatives   					
	);
      double ***double_Matrix2(             // allocate memory for a three-
	  int number_primary_cells_i,		// dimensional array of doubles 
	  int number_primary_cells_j, 		
	  int number_primary_cells_k
	);
      void free_double_Matrix2( 		// deallocate memory for a three			
	  double ***doubleMatrix2,       	// dimensional array of doubles
	  int number_primary_cells_i,	
	  int number_primary_cells_j
	);
      void field_neumann_boundary(          // apply neumann boundary
	  double ***field, 			// condition to cel centered field	
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
	);
	double ***flux_x1;			// volume of fluid flux in x1 direction
        int i_index, j_index, k_index;  	// local variables for loop indexing
	double one_over_dx1	=    		// 1/(grid spacing in x1 direction)
		    1.0/( mesh_width_x1);

     

     /* allocate temporary storage for the fluxes */
      
      flux_x1=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
				  number_primary_cells_k+2);
      
      /* compute fluxes in x1 direction */
      
      compute_vof_flux_x1(level_set, u_1_velocity_new, 			
			    d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,			
			      flux_x1,				
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				 actual_time_step_level_set, mesh_width_x1, lower_bound_derivatives);
      
      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
		  vof_after_x1_update[i_index][j_index][k_index]=
			( volume_of_fluid[i_index][j_index][k_index] -
			    ( flux_x1[i_index][j_index][k_index]-
				  flux_x1[i_index-1][j_index][k_index]))/
		    (1.0-actual_time_step_level_set*one_over_dx1*
			  (u_1_velocity_new[i_index][j_index][k_index]-
			      u_1_velocity_new[i_index-1][j_index][k_index]));
				    
		}
	  }
      }

 /* apply the homogeneous neumann boundary conditions to the level-set field */
 
      field_neumann_boundary(vof_after_x1_update, number_primary_cells_i, 
			     number_primary_cells_j, number_primary_cells_k);

/* de-allocate the temporary storage for the fluxes */

      free_double_Matrix2(flux_x1, number_primary_cells_i+1, number_primary_cells_j+2);
      }

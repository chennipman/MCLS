#include "../headers/array.h"

/********************************************************************************/
/********************************************************************************/
/*  Function to  update the volume of fluid field due to fluxes in x2          	*/
/*  direction										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
//
      void update_volume_of_fluid_x2_flux(				
	 Array3<double> level_set, 		//level-set field			
	 Array3<double> d_level_set_d_x1,		// first partial derivative of
						// the level-set field wrt x1
						// second order central approximation				
	 Array3<double> d_level_set_d_x2,		// first partial derivative of
						// the level-set field wrt x2
						// second order central approximation 				
	 Array3<double> d_level_set_d_x3,		// first partial derivative of
						// the level-set field wrt x3
						// second order central approximation
	 Array3<double> volume_of_fluid, 		// volume of fluid field
	 Array3<double> vof_after_x2_update,	// intermediate volume of fluid field, 
						// after update in x1 direction
	 Array3<double> u_2_velocity_new, 		// velocity field at new time level x1 direction
	 int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction 
	 int number_primary_cells_j,		// number of primary (pressure) cells in x1 direction 
	 int number_primary_cells_k,		// number of primary (pressure) cells in x1 direction
	 double actual_time_step_level_set,	// time step used for level-set advection
						// computed from all stability restrictions and 
						// possibly subscycling
	 double mesh_width_x2,		// grid spacing in x1 direction (uniform)
	 double lower_bound_derivatives	// lower bound for the first partial derivatives
						// to consider it a limiting case of vanishing
						// partial derivatives			
      
  )
     {
      /* function definitions */
      
     void compute_vof_flux_x2(		// compute volume of fluid	
	Array3<double> level_set, 			// fluxes in x2 direction
        Array3<double> u_2_velocity_new, 			
	Array3<double> d_level_set_d_x1,			
	Array3<double> d_level_set_d_x2,			
	Array3<double> d_level_set_d_x3,			
	Array3<double> flux_x2,				
	int number_primary_cells_i,			
	int number_primary_cells_j,			
	int number_primary_cells_k,			
	double actual_time_step_level_set,		
	double mesh_width_x2,				
	double lower_bound_derivatives    		
	);
      void field_neumann_boundary(		// apply neumann boundary
	  Array3<double> field,			// condition to cel centered field	
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
      );
      
      Array3<double> flux_x2;			// volume of fluid flux in x1 direction
      int i_index, j_index, k_index;  	// local variables for loop indexing
      double one_over_dx2	=    		// 1/(grid spacing in x2 direction)
		    1.0/(mesh_width_x2);

     

     /* allocate temporary storage for the fluxes */
      
      flux_x2.create(number_primary_cells_i+2, number_primary_cells_j+1,
				  number_primary_cells_k+2);
      
      /* compute fluxes in x2 direction */
      
      compute_vof_flux_x2(level_set, u_2_velocity_new, 			
			    d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,			
			      flux_x2,				
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				 actual_time_step_level_set, mesh_width_x2, lower_bound_derivatives);
      
      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
		  vof_after_x2_update[i_index][j_index][k_index]=
			( volume_of_fluid[i_index][j_index][k_index] -
			    ( flux_x2[i_index][j_index][k_index]-
				  flux_x2[i_index][j_index-1][k_index]))/
		    (1.0-actual_time_step_level_set*one_over_dx2*
			  (u_2_velocity_new[i_index][j_index][k_index]-
			      u_2_velocity_new[i_index][j_index-1][k_index]));
				    
		}
	  }
      }

 /* apply the homogeneous neumann boundary conditions to the level-set field */
 
      field_neumann_boundary(vof_after_x2_update, number_primary_cells_i, 
			     number_primary_cells_j, number_primary_cells_k);

/* de-allocate the temporary storage for the fluxes */

      flux_x2.destroy();
      }

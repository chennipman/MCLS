#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the body force in the x1 direction                      */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
	
      void compute_body_force_x1( 
      double ***curvature,				// interface curvature 
      double ***level_set, 				// level set field
      double ***surface_tension_body_force_x1,	// x1 component of the body force due to
							// CSF formulation of surface tension model
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double rho_plus_over_rho_minus, 		// rho_plus / rho_minus indicated 
      double sigma_over_rho_minus,			// sigma / rho_minus (scaled sigma)
      double maximum_body_force_x1,			// maximum body force in x1 direction			
      double total_body_force_x1,			// total body force in x1 direction 
      double maximum_weighted_curvature,		// maximum 'active' value of the curvature 
							// used to evaluate the capillary time step
							// restriction 
      double smoothing_distance_factor		// the smoothing distance is smoothing_distance_factor
       )						
       
      /* function definitions */
 {
     void set_constant_matrix2(                    	// set three-dimensional array to constant value
	    int first_dimension,	
	    int second_dimension,	
	    int third_dimension,	
	    double ***matrix2_to_set,	
	    double constant_value	
	  );
      double computed_derivative_heaviside_function(	// use numerical approximation to the 	 
	      double level_set_left,				// derivative of the smoothed heaviside function
	      double level_set_right,	  			// (H(phi_left)-H(phi_right))/(phi_left-phi_right)
	      double mesh_width_x1,			
	      double mesh_width_x2,			
	      double mesh_width_x3,			
	      double smoothing_distance_factor		
	  );
      
      double rho_minus_over_rho_mean=			//  rho_minus/( 0.5* (rho_plus + rho_minus))
		2.0/(1.0+rho_plus_over_rho_minus); 
      double weighted_curvature=0;				// curvature weighted by the derivative of
								// the derivative of the heaviside of the level-set
      double one_over_dx1	=    				// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      int i_index, j_index, k_index;  			// local variables for loop indexing
      
      /* reset total and maximum body force */
      
     maximum_body_force_x1=0;
     total_body_force_x1=0;
	
     /* set body force in x1 direction to zero */
     
     set_constant_matrix2(number_primary_cells_i+1,number_primary_cells_j+2, number_primary_cells_k+2,
				surface_tension_body_force_x1, 0.0);


	
     for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
      {
	
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      
	      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	      {
		/* weight the curvature with the derivative of the heaviside function */
		/* this basically makes it only 'active' near the interface */
		
		weighted_curvature=0.5*( curvature[i_index][j_index][k_index]+
					      curvature[i_index+1][j_index][k_index])*
						computed_derivative_heaviside_function(
						    level_set[i_index][j_index][k_index],
							 level_set[i_index+1][j_index][k_index],	  		
							    mesh_width_x1, mesh_width_x2, mesh_width_x3,			
							      smoothing_distance_factor);
		
		
		surface_tension_body_force_x1[i_index][j_index][k_index]= 
					     -1.0*rho_minus_over_rho_mean*sigma_over_rho_minus*
						one_over_dx1*(level_set[i_index+1][j_index][k_index]-
								level_set[i_index][j_index][k_index])*
								    weighted_curvature;
								    
		/* update the maximum body force in the x1 direction */
		
		maximum_body_force_x1=std::max(maximum_body_force_x1, 
					       fabs(surface_tension_body_force_x1[i_index][j_index][k_index]));

		/* update the maximum weighted curvature in the x1 direction */
		
		maximum_weighted_curvature=std::max(maximum_weighted_curvature,
						      fabs(weighted_curvature));
		total_body_force_x1+=surface_tension_body_force_x1[i_index][j_index][k_index];
		      
	      }
	  }
      }
	  	std::cerr<< " maximum body force x1 direction "<< maximum_body_force_x1 <<"\n";
 	  	std::cerr<< " maximum weighted curvature x1 "<< maximum_weighted_curvature <<"\n";
 }
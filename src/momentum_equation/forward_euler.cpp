#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
      
/********************************************************************************/
/* Function for one forward euler step						*/
/*  										*/
/*  Programmer	: Coen Hennipman       						*/
/*  Date	: 06-03-2014    						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This functions computes one explicit euler step.				*/
/* 									*/
/********************************************************************************/

 EXPORT void forward_euler(
      Array3<double> u_1_velocity_star, 		// velocity field at star time level x1 direction
      Array3<double> u_2_velocity_star, 		// velocity field at star time level x2 direction
      Array3<double> u_3_velocity_star,			// velocity field at star time level x3 direction
      Array3<double> u_1_new_con_diff,                  // convection and diffusion in x1 direction   
      Array3<double> u_2_new_con_diff,                  // convection and diffusion in x2 direction     
      Array3<double> u_3_new_con_diff,                  // convection and diffusion in x3 direction      
      Array3<double> u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction
      Array3<double> u_1_old_con_diff,                  // convection and diffusion in x1 direction   
      Array3<double> u_2_old_con_diff,                  // convection and diffusion in x2 direction     
      Array3<double> u_3_old_con_diff,                  // convection and diffusion in x3 direction     
      Array3<double> scaled_density_u1,			// scaled density for the controlvolumes
      Array3<double> scaled_density_u2,			// scaled density for the controlvolumes
      Array3<double> scaled_density_u3,			// scaled density for the controlvolumes      
      
      Array3<double> momentum_source_term_u_1,		// source term of the momentum equation in x1 direction
					       		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2,		// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3,		// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> level_set, 			// level-set field
      double actual_time_step_navier_stokes,	        // time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      double sigma,					// constant for the convection and diffusion terms of the previous stage
      double zeta,					// constant for the convection and diffusion terms of the actual stage
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      double smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
      double rho_plus_over_rho_minus,
      double rho_minus_over_mu_minus,			// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      boundary_face boundary_faces[6],		        // array with all the information
							// for the boundary conditions 
      int source_terms_in_momentum_predictor,    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation  
      double actual_time				// actual time 
       )

  {
       // calculate the convection and diffusion(without momentum_source_terms)
       // the momentum source terms are included below or in the pressure correction
       // depending on 'source _terms_in_momentum_predictor' 
	convection_diffussion_source_terms(
       u_1_new_con_diff,u_2_new_con_diff,u_3_new_con_diff,
       u_1_velocity_old,u_2_velocity_old,u_3_velocity_old,
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       0);

	// calculate u* in all 3 directions
  	if(source_terms_in_momentum_predictor) // the source terms are in the predictor
      	{
       add_4_arrays(u_1_velocity_star,	1.0, u_1_velocity_old,	
       actual_time_step_navier_stokes*sigma, u_1_new_con_diff, 	
       actual_time_step_navier_stokes*zeta, u_1_old_con_diff, 
       actual_time_step_navier_stokes*(sigma+zeta), momentum_source_term_u_1, 
       number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       
       add_4_arrays(u_2_velocity_star,	1.0, u_2_velocity_old,	
       actual_time_step_navier_stokes*sigma, u_2_new_con_diff,
       actual_time_step_navier_stokes*zeta, u_2_old_con_diff,
       actual_time_step_navier_stokes*(sigma+zeta), momentum_source_term_u_2,  	
       number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       
       add_4_arrays(u_3_velocity_star,	1.0,u_3_velocity_old,	
       actual_time_step_navier_stokes*sigma, u_3_new_con_diff,
       actual_time_step_navier_stokes*zeta, u_3_old_con_diff,
       actual_time_step_navier_stokes*(sigma+zeta), momentum_source_term_u_3, 	
       number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
      	}
      	else 	// the source terms are in the pressure correction		
      	{
       add_3_arrays(u_1_velocity_star,	1.0,u_1_velocity_old,	
       actual_time_step_navier_stokes*sigma, u_1_new_con_diff, 	
       actual_time_step_navier_stokes*zeta, u_1_old_con_diff, 	
       number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       
       add_3_arrays(u_2_velocity_star,	1.0,u_2_velocity_old,	
       actual_time_step_navier_stokes*sigma, u_2_new_con_diff, 
       actual_time_step_navier_stokes*zeta, u_2_old_con_diff,	
       number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       
       add_3_arrays(u_3_velocity_star,	1.0,u_3_velocity_old,	
       actual_time_step_navier_stokes*sigma, u_3_new_con_diff, 
       actual_time_step_navier_stokes*zeta, u_3_old_con_diff, 	
       number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
       }
       
	// shift the convection and diffusion term
      copy_general_field(u_1_new_con_diff, u_1_old_con_diff,
                       0, number_primary_cells_i,
                         0, number_primary_cells_j+1,
                           0, number_primary_cells_k+1);
      
      copy_general_field(u_2_new_con_diff, u_2_old_con_diff,
                       0, number_primary_cells_i+1,
                         0, number_primary_cells_j,
                           0, number_primary_cells_k+1);
  
      copy_general_field(u_3_new_con_diff, u_3_old_con_diff,
                       0, number_primary_cells_i+1,
                         0, number_primary_cells_j+1,
                           0, number_primary_cells_k);
                           
	// apply the boundary conditions 
       apply_boundary_conditions_velocity(boundary_faces,		
	  u_1_velocity_star, u_2_velocity_star, u_3_velocity_star, 			 
	    mesh_width_x1, mesh_width_x2, mesh_width_x3,				 
	      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
	        actual_time+(sigma+zeta)*actual_time_step_navier_stokes);	 


  }      










      	

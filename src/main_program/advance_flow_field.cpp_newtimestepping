#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to advance the flow field to the next time                   t     */
/*                     								*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* Advancing the flow field consists of three steps:                            */
/* -First a tentative velocity field is computed at the new time level          */
/* -Next this tentative velocity field is projected on the space of divergence  */
/* free velocity field                                                          */
/* -Finally, the velocity fields are shifted: old becomes new                   */
/********************************************************************************/
EXPORT void advance_flow_field(
      Array3<double> level_set,				// level-set field
      Array3<double> pressure,				// pressure field
      Array3<double> u_1_velocity_old,			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old,			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction
      Array3<double> u_1_velocity_new,			// velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new,			// velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,			// velocity field at new time level x3 direction
      Array3<double> momentum_source_term_u_1,		// source term of the momentum equation in x1 direction
					       		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2,		// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3,		// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x1,	// source term of the momentum equation in x1 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x2,	// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x3,	// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> scaled_density_u1,			// scaled density for the controlvolumes
							// of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,			// scaled density for the controlvolumes
							// of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,			// scaled density for the controlvolumes
							// of the momentum equation in x3 direction
      boundary_face boundary_faces[6],		// array with all the information
							// for the boundary conditions 
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      int number_matrix_connections,			// number of connections in momentum matrix				       
      vector gravity,					// gravitational acceleration vector 
      double actual_time_step_navier_stokes,	// time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      double rho_plus_over_rho_minus,		// ratio of the densities of the two phases
      double smoothing_distance_factor,		// the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
      double rho_minus_over_mu_minus,		// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      double tolerance_pressure,	  		// the tolerance with which the system for the pressure is solved	
      int maximum_iterations_allowed_pressure,	// maximum number of iterations allowed for the
							// conjugate gradient method
      double tolerance_velocity,	  		// the tolerance with which the system for the momentum predictor is solved	
      int maximum_iterations_allowed_velocity,	// maximum number of iterations allowed for the
							// conjugate gradient method
      int continuous_surface_force_model,       	// =1, the continuous surface force model is applied
					        	// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation
	 )
        {
       Array3<double> u_1_velocity_star; 		     // velocity field at star time level x1 direction   
       Array3<double> u_2_velocity_star; 		     // velocity field at star time level x2 direction     
       Array3<double> u_3_velocity_star;		     // velocity field at star time level x3 direction
       
       Array3<double> u_1_con_diff; 		     	// convection and diffusion in  x1 direction   
       Array3<double> u_2_con_diff; 		     	// convection and diffusion in x2 direction     
       Array3<double> u_3_con_diff;		     	// convection and diffusion in x3 direction


    	/* allocate memory for tentative velocity field u_star */
      u_1_velocity_star.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_star.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_star.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

	// these terms are the output of the momentum_predictor
      u_1_con_diff.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
       

	/* Explicit Euler  */           								
    	// calculate momentum predictor 
	momentum_predictor(
       u_1_con_diff,u_2_con_diff,u_3_con_diff,
       u_1_velocity_old,u_2_velocity_old,u_3_velocity_old,
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       source_terms_in_momentum_predictor);

       add_arrays(u_1_velocity_star,	1.0,u_1_velocity_old,	actual_time_step_navier_stokes,u_1_con_diff, 	number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       add_arrays(u_2_velocity_star,	1.0,u_2_velocity_old,	actual_time_step_navier_stokes,u_2_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       add_arrays(u_3_velocity_star,	1.0,u_3_velocity_old,	actual_time_step_navier_stokes,u_3_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
	// end Explicit Euler //
	

	/* Runge-Kutta not working yet           								
	// the parameters are based on Wesseling's Principles of Computational Fluid Dynamics, p199
	// u_star is used as an intermediate velocity in one of the steps and as the final output for the rest of the calculation 
	// allocate memory for tentative velocity field u_star_star
       Array3<double> u_1_velocity_star_star; 		     // velocity field at star_star time level x1 direction   
       Array3<double> u_2_velocity_star_star; 		     // velocity field at star_star time level x2 direction     
       Array3<double> u_3_velocity_star_star;		     // velocity field at star_star time level x3 direction
      u_1_velocity_star_star.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_star_star.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_star_star.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

    	// allocate memory for tentative velocity field u_input
       Array3<double> u_1_velocity_input; 		     // velocity field at input level x1 direction   
       Array3<double> u_2_velocity_input; 		     // velocity field at input level x2 direction     
       Array3<double> u_3_velocity_input;		     // velocity field at input level x3 direction
      u_1_velocity_input.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_input.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_input.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      
      // first step Runge-Kutta  
      // calculate input 
       add_arrays(u_1_velocity_input,	8.0/15.0,u_1_velocity_old,	0.0,u_1_velocity_old, 	number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       add_arrays(u_2_velocity_input,	8.0/15.0,u_2_velocity_old,	0.0,u_2_velocity_old, 	number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       add_arrays(u_3_velocity_input,	8.0/15.0,u_3_velocity_old,	0.0,u_3_velocity_old, 	number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
    
	momentum_predictor(
       u_1_con_diff,u_2_con_diff,u_3_con_diff,
       u_1_velocity_input,u_2_velocity_input,u_3_velocity_input,
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       source_terms_in_momentum_predictor);

       // calculate u_star
       add_arrays(u_1_velocity_star,	1.0,u_1_velocity_old,	actual_time_step_navier_stokes,u_1_con_diff, 	number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       add_arrays(u_2_velocity_star,	1.0,u_2_velocity_old,	actual_time_step_navier_stokes,u_2_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       add_arrays(u_3_velocity_star,	1.0,u_3_velocity_old,	actual_time_step_navier_stokes,u_3_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);

	// second step Runge-Kutta
	// calculate input for second step Runge-Kutta
       add_arrays(u_1_velocity_input,	17.0/60.0,u_1_velocity_old,	5.0/12.0,u_1_velocity_star, 	number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       add_arrays(u_2_velocity_input,	17.0/60.0,u_2_velocity_old,	5.0/12.0,u_2_velocity_star, 	number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       add_arrays(u_3_velocity_input,	17.0/60.0,u_3_velocity_old,	5.0/12.0,u_3_velocity_star, 	number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
             
     	momentum_predictor(
       u_1_con_diff,u_2_con_diff,u_3_con_diff,
       u_1_velocity_input,u_2_velocity_input,u_3_velocity_input,
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       source_terms_in_momentum_predictor); 
      
      // calculate u_star_star
       add_arrays(u_1_velocity_star_star,	1.0,u_1_velocity_star,	actual_time_step_navier_stokes,u_1_con_diff, 	number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       add_arrays(u_2_velocity_star_star,	1.0,u_2_velocity_star,	actual_time_step_navier_stokes,u_2_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       add_arrays(u_3_velocity_star_star,	1.0,u_3_velocity_star,	actual_time_step_navier_stokes,u_3_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
      
      // calculate input for third step Runge-Kutta
       add_arrays(u_1_velocity_input,	17.0/60.0,u_1_velocity_old,	5.0/12.0,u_1_velocity_star, 	number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       add_arrays(u_2_velocity_input,	17.0/60.0,u_2_velocity_old,	5.0/12.0,u_2_velocity_star, 	number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       add_arrays(u_3_velocity_input,	17.0/60.0,u_3_velocity_old,	5.0/12.0,u_3_velocity_star, 	number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
        
     	momentum_predictor(
       u_1_con_diff,u_2_con_diff,u_3_con_diff,
       u_1_velocity_input,u_2_velocity_input,u_3_velocity_input,
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       source_terms_in_momentum_predictor); 
       
       // calculate output hereby reusing u_star     
       add_arrays(u_1_velocity_star,	1.0,u_1_velocity_star_star,	actual_time_step_navier_stokes,u_1_con_diff, 	number_primary_cells_i+1,number_primary_cells_j+2,number_primary_cells_k+2);
       add_arrays(u_2_velocity_star,	1.0,u_2_velocity_star_star,	actual_time_step_navier_stokes,u_2_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+1,number_primary_cells_k+2);
       add_arrays(u_3_velocity_star,	1.0,u_3_velocity_star_star,	actual_time_step_navier_stokes,u_3_con_diff, 	number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+1);
	// end Runge-Kutta */
	

      /* apply the boundary conditions to u star, and fill the virtual cells */
      apply_boundary_conditions_velocity(boundary_faces,		
					  u_1_velocity_star, u_2_velocity_star, u_3_velocity_star, 			 
					    mesh_width_x1, mesh_width_x2, mesh_width_x3,				 
					      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);			 
	 
     /* solve momentum corrector equation */
     /* compute a correction to the velocity field u star, to make it divergence free */
     /* and apply this correction to u star */
     
      solve_momentum_corrector(	level_set, pressure,			
				momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,	
				  surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
				    scaled_density_u1, scaled_density_u2, scaled_density_u3,
				    u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,	        
				      mesh_width_x1, mesh_width_x2, mesh_width_x3,		        
					number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	        
					  gravity, tolerance_pressure, actual_time_step_navier_stokes,    
					    rho_plus_over_rho_minus, continuous_surface_force_model,       
					      source_terms_in_momentum_predictor, maximum_iterations_allowed_pressure,	 	
						boundary_faces);
  
      /* shift the velocity field */
      
      // shift the velocity fields
      // u new -> u old
      // u star -> u new
      
       shift_velocity_field( u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,	     	
			      u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,	     	
				u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,	     	
				  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
       
       
       
      /* deallocate  memory for tentative velocity field u star */
       
	u_1_velocity_star.destroy();
	u_2_velocity_star.destroy();
	u_3_velocity_star.destroy();

      u_1_con_diff.destroy();
      u_2_con_diff.destroy();
      u_3_con_diff.destroy();


	
     }
	
	
	
	
	
	
	
	
	

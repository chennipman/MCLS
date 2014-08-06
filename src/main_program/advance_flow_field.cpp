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
      boundary_face boundary_faces[6],		        // array with all the information
							// for the boundary conditions 
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      int number_matrix_connections,			// number of connections in momentum matrix				       
      vector gravity,					// gravitational acceleration vector 
      double actual_time_step_navier_stokes,	        // time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      time_stepping_methods time_stepping_method, 	// time scheme 1:explicit euler 2: imex, 3: runge-kutta 
      double rho_plus_over_rho_minus,		        // ratio of the densities of the two phases
      double smoothing_distance_factor,		        // the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
      double rho_minus_over_mu_minus,		        // this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      double tolerance_pressure,	  		// the tolerance with which the system for the pressure is solved	
      int maximum_iterations_allowed_pressure,	        // maximum number of iterations allowed for the
							// conjugate gradient method
      double tolerance_velocity,	  		// the tolerance with which the system for the momentum predictor is solved	
      int maximum_iterations_allowed_velocity,	        // maximum number of iterations allowed for the
							// conjugate gradient method
      int continuous_surface_force_model,       	// =1, the continuous surface force model is applied
					        	// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor,    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation
      double actual_time
	 )
        {
       Array3<double> u_1_velocity_star; 		// velocity field at star time level x1 direction   
       Array3<double> u_2_velocity_star; 		// velocity field at star time level x2 direction     
       Array3<double> u_3_velocity_star;		// velocity field at star time level x3 direction

    /* allocate memory for tentative velocity field u star */
    
      u_1_velocity_star.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_star.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_star.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      
    /* solve momentum predictor equation */
    /* compute a new velocity field u star, that is not divergence free */
    
    if (time_stepping_method==explicit_euler) // Explicit Euler
    {
      printf("explicit euler \n"); 
      
	/* Explicit Euler  */           								
    	// calculate momentum predictor 
      solve_momentum_predictor_explicit(
       u_1_velocity_star,u_2_velocity_star,u_3_velocity_star,
       u_1_velocity_old,u_2_velocity_old,u_3_velocity_old,
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set, actual_time_step_navier_stokes,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       boundary_faces, source_terms_in_momentum_predictor, actual_time);      
    }
    
    else if (time_stepping_method==imex)
    {
      printf("imex \n");
      solve_momentum_predictor_imex(	level_set, u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,			
				  u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,			
				    momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3, 	
				     scaled_density_u1, scaled_density_u2, scaled_density_u3,
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
					number_matrix_connections, actual_time_step_navier_stokes,		
					  mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					    rho_plus_over_rho_minus, smoothing_distance_factor, rho_minus_over_mu_minus,			
					      mu_plus_over_mu_minus, boundary_faces,
					      tolerance_velocity, maximum_iterations_allowed_velocity,actual_time); 
    }

    else if ( (time_stepping_method==runge_kutta) || (time_stepping_method==two_pres_solve) || (time_stepping_method==two_pres_solve_output) )
    {
      printf("runge kutta \n"); 
      solve_momentum_predictor_rk(
       u_1_velocity_star,u_2_velocity_star,u_3_velocity_star,
       u_1_velocity_old,u_2_velocity_old,u_3_velocity_old,
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set, actual_time_step_navier_stokes,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       boundary_faces,source_terms_in_momentum_predictor,actual_time);           
    }
    
    else 
    {
      printf("unkown time stepping method \n"); 
    }

  
     /* solve momentum corrector equation 						*/
     /* compute the new pressure that makes the velocity field u star divergence free 	*/
     /* apply the new computed pressure to u star 					*/

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
						boundary_faces, actual_time);


    if (time_stepping_method==two_pres_solve) // extra pressure correction step every time step
    {
       printf("second pressure correction step \n"); 
     /* the fuction below gives a pressure field of the same order as the velocity field 	*/
     /* it is based on eq 44 of 								*/
     /* New explicit Runge-Kutta methods for the incompressible Navier_Stokes equations		*/
     /* by: Bejamin Sanderse and B. Koren Bibref:SanRKPCM1 					*/
      solve_final_pressure_corrector(	level_set, pressure,			
				u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,	        
				momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
				scaled_density_u1, scaled_density_u2, scaled_density_u3,
 				mesh_width_x1, mesh_width_x2, mesh_width_x3,		        
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	        
				gravity, tolerance_pressure,smoothing_distance_factor,
				rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
				maximum_iterations_allowed_pressure);
    }
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
	
     }
	
	
	
	
	
	
	
	
	

#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
      
/********************************************************************************/
/* Function for momentum prediction in all directions				*/
/*  										*/
/*  Programmer	: Coen Hennipman       						*/
/*  Date	: 06-03-2014    						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* u_star is reused for u***							*/
/* Both this Runge Kutta procdeure and the forward_euler are designed for 	*/
/* low storage schemes, therefore only the previous terms from the con_diff	*/
/* can be incorporated. 							*/
/********************************************************************************/

 EXPORT void solve_momentum_predictor_rk(
      Array3<double> u_1_velocity_star, 		// velocity field at star time level x1 direction
      Array3<double> u_2_velocity_star, 		// velocity field at star time level x2 direction
      Array3<double> u_3_velocity_star,			// velocity field at star time level x3 direction
      Array3<double> u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction
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
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      double smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
      double rho_plus_over_rho_minus,			// ratio of the densities of the two phases
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
	double sigma, zeta; 
 
       Array3<double> u_1_old_con_diff; 	// contains the convection and diffusion terms of a previous RK-stage in the x1 direction
       Array3<double> u_2_old_con_diff; 	// contains the convection and diffusion terms of a previous RK-stage in the x2 direction
       Array3<double> u_3_old_con_diff;		// contains the convection and diffusion terms of a previous RK-stage in the x3 direction

      u_1_old_con_diff.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_old_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_old_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);	

       Array3<double> u_1_new_con_diff; 	// contains the convection and diffusion terms of a RK-stage in the x1 direction
       Array3<double> u_2_new_con_diff; 	// contains the convection and diffusion terms of a RK-stage in the x2 direction	
       Array3<double> u_3_new_con_diff;		// contains the convection and diffusion terms of a RK-stage in the x3 direction	

      u_1_new_con_diff.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_new_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_new_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);	


	// first step Runge-Kutta
	sigma = 8.0/15.0; // parameter for new stage
	zeta  = 0.0; // parameter for previous stage

      forward_euler(
	u_1_velocity_star,u_2_velocity_star,u_3_velocity_star,			
	u_1_new_con_diff,u_2_new_con_diff,u_3_new_con_diff,                      
	u_1_velocity_old,u_2_velocity_old,u_3_velocity_old,		
	u_1_old_con_diff,u_2_old_con_diff,u_3_old_con_diff,               
 	scaled_density_u1,scaled_density_u2,scaled_density_u3,			   
 	momentum_source_term_u_1,momentum_source_term_u_2,momentum_source_term_u_3,	
	level_set, actual_time_step_navier_stokes,	       
  	sigma,zeta,					
  	number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,			
	mesh_width_x1,mesh_width_x2,mesh_width_x3,				
	smoothing_distance_factor,
 	rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,			
 	boundary_faces, source_terms_in_momentum_predictor, actual_time    	
       );
       
       actual_time = actual_time+ (sigma+zeta)*actual_time_step_navier_stokes; 

	// second step Runge-Kutta
       Array3<double> u_1_velocity_star_star; 		// velocity field at star_star time level x1 direction   
       Array3<double> u_2_velocity_star_star; 		// velocity field at star_star time level x2 direction     
       Array3<double> u_3_velocity_star_star;		// velocity field at star_star time level x3 direction

    	/* allocate memory for tentative velocity field u star_star */
      u_1_velocity_star_star.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_star_star.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_star_star.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

	sigma =  5.0/12.0; // parameter for new stage
	zeta  = -17.0/60.0; // parameter for previous stage
		
      forward_euler(
	u_1_velocity_star_star,u_2_velocity_star_star,u_3_velocity_star_star,			
	u_1_new_con_diff,u_2_new_con_diff,u_3_new_con_diff,                      
	u_1_velocity_star,u_2_velocity_star,u_3_velocity_star,		
	u_1_old_con_diff,u_2_old_con_diff,u_3_old_con_diff,               
 	scaled_density_u1,scaled_density_u2,scaled_density_u3,			   
 	momentum_source_term_u_1,momentum_source_term_u_2,momentum_source_term_u_3,	
	level_set, actual_time_step_navier_stokes,	       
  	sigma,zeta,					
  	number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,			
	mesh_width_x1,mesh_width_x2,mesh_width_x3,				
	smoothing_distance_factor,
 	rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,			
 	boundary_faces, source_terms_in_momentum_predictor, actual_time
       );
       
       actual_time = actual_time+ (sigma+zeta)*actual_time_step_navier_stokes; 

	// third step Runge-Kutta
	sigma = 3.0/4.0; // parameter for new stage
	zeta  = -5.0/12.0; // parameter for previous stage
		
      forward_euler(
	u_1_velocity_star,u_2_velocity_star,u_3_velocity_star,		// u_star is reused to reduce the number of allocations it can also be seen as u_star_star_star	
	u_1_new_con_diff,u_2_new_con_diff,u_3_new_con_diff,                      
	u_1_velocity_star_star,u_2_velocity_star_star,u_3_velocity_star_star,		
	u_1_old_con_diff,u_2_old_con_diff,u_3_old_con_diff,               
 	scaled_density_u1,scaled_density_u2,scaled_density_u3,			   
 	momentum_source_term_u_1,momentum_source_term_u_2,momentum_source_term_u_3,	
	level_set, actual_time_step_navier_stokes,	       
  	sigma,zeta,					
  	number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,			
	mesh_width_x1,mesh_width_x2,mesh_width_x3,				
	smoothing_distance_factor,
 	rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,			
 	boundary_faces, source_terms_in_momentum_predictor, actual_time    	
       );


	u_1_old_con_diff.destroy();
	u_2_old_con_diff.destroy();
	u_3_old_con_diff.destroy();
	
	u_1_new_con_diff.destroy();
	u_2_new_con_diff.destroy();
	u_3_new_con_diff.destroy();

	u_1_velocity_star_star.destroy();
	u_2_velocity_star_star.destroy();
	u_3_velocity_star_star.destroy();

  }      










      	

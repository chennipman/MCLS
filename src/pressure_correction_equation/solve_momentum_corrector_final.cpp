#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to solve for the correction to the pressure field			*/
/*  This function is partly a copy of solve_momentum_corrector.			*/
/*  										*/
/*  Programmer	: Coen Hennipman						*/
/*  Date	: 20-05-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* Note that in all cases the PRESSURE is computed in the pressure correction   */
/* equation, not the pressure correction.                                       */
/* Currently we assume a Dirichlet boundary condition for all normal velocities.*/
/* The fuction below gives a pressure field of the same order as the velocity  	*/
/* field, it is based on eq 44 of 						*/
/* 'New explicit Runge-Kutta methods for the incompressible Navier_Stokes	*/
/* equations' written by: Bejamin Sanderse and B. Koren Bibref:SanRKPCM1	*/
/* This function is almost similar to solve_momentum_corrector,			*/
/* the difference is the build of the rhs. 					*/
/* eq 44: Lp_{n+1}=MF_{n+1}-r_1^.(t_{n+1})					*/
/* where:									*/
/* L is the laplaciaan								*/
/* p_{n+1} is the pressure at the new time_stepping_method			*/
/* M is the discrete convergence 						*/
/* F_{n+1} are the convection_diffussion terms at the new time 			*/
/* r_1^.(t_{n+1}) is the time derivative of the divergence which is asssumed to */
/* be zero in our case. 							*/
/* By applying this method the order of the pressure field is of the same order	*/
/* as the velocity field. 							*/
/********************************************************************************/

EXPORT void solve_momentum_corrector_final(
      Array3<double> level_set,					// level-set field
      Array3<double> pressure,					// pressure field
      Array3<double> u_1_velocity_star, 			// velocity field at star time level x1 direction
      Array3<double> u_2_velocity_star, 			// velocity field at star time level x2 direction
      Array3<double> u_3_velocity_star,				// velocity field at star time level x3 direction
      Array3<double> momentum_source_term_u_1,		// source term of the momentum equation in x1 direction
					       		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2,		// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3,		// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> scaled_density_u1,				// scaled density for the controlvolumes
								// of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,				// scaled density for the controlvolumes
								// of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,				// scaled density for the controlvolumes
								// of the momentum equation in x3 direction
      double mesh_width_x1,		        		// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		        		// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		        		// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	        		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	        		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	        		// number of primary (pressure) cells in x3 direction
      vector gravity,						// gravitational acceleration vector 
      double tolerance,	  					// the tolerance with which the system is solved	
      double smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
      double rho_plus_over_rho_minus,			// ratio of the densities of the two phases
      double rho_minus_over_mu_minus,			// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      int source_terms_in_momentum_predictor,    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation  
      int maximum_iterations_allowed_pressure,	 			// maximum number of iterations allowed for the
								// conjugate gradient method
      boundary_face boundary_faces[6],				// array with all the information
								// for the boundary conditions 
      double actual_time					// actual time
	
      )
   
   
   {
	// initialize the convection and diffusion terms based on the u_star velocity field  
       Array3<double> u_1_new_con_diff; 	
       Array3<double> u_2_new_con_diff; 		
       Array3<double> u_3_new_con_diff;			

      u_1_new_con_diff.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_new_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_new_con_diff.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      
	set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2, 
			    number_primary_cells_k+2, u_1_new_con_diff, 0.0); 
	set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, u_2_new_con_diff, 0.0); 
	set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, u_3_new_con_diff, 0.0); 

	// compute the convection and diffusion terms for in the pressure corrector
	convection_diffussion_source_terms(
       u_1_new_con_diff,u_2_new_con_diff,u_3_new_con_diff,
       u_1_velocity_star,u_2_velocity_star,u_3_velocity_star, // the terms are calculated based on u_star, normally based on u_n
       scaled_density_u1,scaled_density_u2,scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set,
       number_primary_cells_i,number_primary_cells_j,number_primary_cells_k,
       mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,
       rho_plus_over_rho_minus,rho_minus_over_mu_minus,mu_plus_over_mu_minus,
       0); // this function is used here only for calculating the convection_diffussion terms, the momentum_source_terms are not incorporated here

      Array2<double> pressure_matrix; 					// pressure matrix
      Array1<double> pressure_rhside;					// pressure right hand side
      Array3<double> pressure_boundary_condition_x1;
      Array3<double> pressure_boundary_condition_x2;
      Array3<double> pressure_boundary_condition_x3;
      int total_number_pressure_points;				// total number of points with pressure

    /* allocate memory for the pressure correction matrix and right hand side */
   
      total_number_pressure_points=number_primary_cells_i*number_primary_cells_j*number_primary_cells_k;
    
      pressure_matrix.create(4,total_number_pressure_points);
      pressure_rhside.create(total_number_pressure_points);

    /* allocate memory for the pressure boundary condition */
   
      pressure_boundary_condition_x1.create(2,
                                       number_primary_cells_j+2,
                                         number_primary_cells_k+2);
      pressure_boundary_condition_x2.create(2,
                                       number_primary_cells_i+2,
                                         number_primary_cells_k+2);
      pressure_boundary_condition_x3.create(2,
                                       number_primary_cells_i+2,
                                         number_primary_cells_j+2);
   
   /* build the system of equations for the pressure correction equation */
  
      build_pressure_system_final(pressure_matrix,   
					pressure_rhside, level_set,			
					    scaled_density_u1, scaled_density_u2, scaled_density_u3,
					     u_1_new_con_diff, u_2_new_con_diff, u_3_new_con_diff,
                                              pressure_boundary_condition_x1, pressure_boundary_condition_x2,
                                                pressure_boundary_condition_x3,
						  mesh_width_x1, mesh_width_x2, mesh_width_x3,		        
						   number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	        
						     rho_plus_over_rho_minus,	        
						        gravity);

   
      /* solve the system of equations for the pressure correction equation */

      solve_pressure_correction_system(pressure_matrix, pressure_rhside, pressure,			  
					number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		  
					  tolerance, maximum_iterations_allowed_pressure);

     /* apply boundary conditions to pressure field */
   
      apply_boundary_conditions_pressure( pressure,
                                          pressure_boundary_condition_x1,
                                             pressure_boundary_condition_x2,
                                                pressure_boundary_condition_x3,
                                                  scaled_density_u1, scaled_density_u2, scaled_density_u3,
                                                    mesh_width_x1, mesh_width_x2, mesh_width_x3,
                                                      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
   
    /* deallocate memory for the pressure correction matrix and right hand side */
   
      pressure_matrix.destroy();
      pressure_rhside.destroy();

    /* deallocate the memory for the pressure boundary conditions 		*/

      pressure_boundary_condition_x1.destroy();
      pressure_boundary_condition_x2.destroy();
      pressure_boundary_condition_x3.destroy();
      
    /* deallocate the memory for the convection and diffusion terms 		*/

	u_1_new_con_diff.destroy();
	u_2_new_con_diff.destroy();
	u_3_new_con_diff.destroy();
      
  }

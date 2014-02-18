#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to initialize the pressure field in the initial condition. It uses */
/*  the same set of functions as the solution of the 'normal' pressure equation */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* Note that in all cases the PRESSURE is computed in the pressure correction   */
/* equation, not the pressure correction.                                       */
/* Currently we assume a Dirichlet boundary condition for all normal velocities.*/
/********************************************************************************/

EXPORT void initialize_pressure(
      Array3<double> level_set,					// level-set field
      Array3<double> pressure,					// pressure field
      Array3<double> momentum_source_term_u_1,			// source term of the momentum equation in x1 direction
					        		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2,  			// source term of the momentum equation in x2 direction
					        		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3,			// source term of the momentum equation in x3 direction
					        		// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x1,  		// source term of the momentum equation in x1 direction
					        		// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x2,  		// source term of the momentum equation in x2 direction
					        		// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x3,		// source term of the momentum equation in x3 direction
					        		// defined on all u1 points (including boundaries)
      Array3<double> scaled_density_u1,				// scaled density for the controlvolumes
								// of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,				// scaled density for the controlvolumes
								// of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,				// scaled density for the controlvolumes
								// of the momentum equation in x3 direction
      Array3<double> u_1_velocity, 	        		// velocity field x1 direction
      Array3<double> u_2_velocity, 	       		 	// velocity field x2 direction
      Array3<double> u_3_velocity,	        		       // velocity field x3 direction
      double mesh_width_x1,		        		// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		        		// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		        		// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	        		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	        		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	        		// number of primary (pressure) cells in x3 direction
      vector gravity,						// gravitational acceleration vector 
      double tolerance,	  				// the tolerance with which the system is solved	
      double actual_time_step_navier_stokes,    		// actual time step for Navier-Stokes solution algorithm 
      double rho_plus_over_rho_minus,	        	// ratio of density where (level set >0) and 
					        		// density where (level set < 0)
      int continuous_surface_force_model,       		// =1, the continuous surface force model is applied
					        		// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor,   		// =1, the source terms are applied in the momentum predictor
					        		// equation
					        		// =0, the source terms are applied in the momentum corrector
					        		// equation
      int maximum_iterations_allowed,	 		// maximum number of iterations allowed for the
								// conjugate gradient method
      boundary_face boundary_faces[6]			// array with all the information
								// for the boundary conditions 
	
      )
   
   
   {
      Array2<double> pressure_matrix; 				   // pressure matrix
      Array1<double> pressure_rhside;				   // pressure right hand side
      int total_number_pressure_points;		          // total number of points with pressure

    /* allocate memory for the pressure correction matrix and right hand side */
   
      total_number_pressure_points=number_primary_cells_i*number_primary_cells_j*number_primary_cells_k;
    
      pressure_matrix.create(4,total_number_pressure_points);
      pressure_rhside.create(total_number_pressure_points);
   
   
   /* build the system of equations for the pressure correction equation */
  
      build_pressure_system_initialization(pressure_matrix,
					pressure_rhside, level_set,			
					  momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,  		        
					    surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
					     scaled_density_u1, scaled_density_u2, scaled_density_u3,
					      u_1_velocity, u_2_velocity, u_3_velocity,
						mesh_width_x1, mesh_width_x2, mesh_width_x3,		        
						 number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	        
						  actual_time_step_navier_stokes, rho_plus_over_rho_minus,	        
						    continuous_surface_force_model, source_terms_in_momentum_predictor,   
						     gravity);

   
   /* solve the system of equations for the pressure correction equation */

      solve_pressure_correction_system(pressure_matrix, pressure_rhside, pressure,			  
					number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		  
					  tolerance, maximum_iterations_allowed);
   
     field_extrapolate_boundary(pressure, number_primary_cells_i,
                        number_primary_cells_j,number_primary_cells_k);
   
     
    /* allocate memory for the pressure correction matrix and right hand side */
   
      pressure_matrix.destroy();
      pressure_rhside.destroy();

  }

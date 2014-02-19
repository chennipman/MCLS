#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
/********************************************************************************/
/*  Function to solve the tentative velocity field u star, x3 components        */
/*                              							*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/*  											*/
/*  											*/
/*  											*/
/*  											*/
/********************************************************************************/
  
EXPORT void solve_momentum_predictor_u3(
      Array3<double> level_set, 				// level-set field
      Array3<double> momentum_source_term_u_3, 		// complete source term for the momentum equation
							// in x3 direction=(-p,3+ g_3 +F3)
      Array3<double> scaled_density_u3,			// scaled density for the controlvolumes
							// of the momentum equation in x3 direction
      Array3<double> u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction
      Array3<double> u_3_velocity_star,			// velocity field at star time level x1 direction

      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      int number_matrix_connections,			// number of connections in momentum matrix				       
      double actual_time_step_navier_stokes,	// time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      double rho_plus_over_rho_minus,		// ratio of the densities of the two phases
      double smoothing_distance_factor,		// the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
      double rho_minus_over_mu_minus,		// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      boundary_face boundary_faces[6],		// array with all the information
							// for the boundary conditions 
      double tolerance_velocity,	  		// the tolerance with which the system for the momentum predictor is solved	
      int maximum_iterations_allowed_velocity	// maximum number of iterations allowed for the
							// conjugate gradient method
							// for the boundary conditions 
	   )
      {
      Array2<double> momentum_matrix_u3;			// momentum matrix velocity x3 direction
      Array1<double> momentum_rhside_u3;			// momentum rhside velocity x3 direction
      Array1<double> preconditioner_matrix_M;		// preconditioner matrix (diagonal)
      Array1<double> compressed_velocity_u3;		// compressed solution vector
      double relative_L2_norm_residual; 	  	// the L2 norm of the residual, scaled with
							// the L2 norm of the right hand side
      double relative_Linfinity_norm_residual;  	// the L infinity norm of the residual
							// scaled with the maximum difference 
							// between two components of the residual
      int iteration_number;				// the number of iterations where the iterative
							// process was terminated
      int total_number_u_3_points;			// total number of grid points with u_3
  
      
      
      
    
    /* allocate memory for momentum matrix, preconditioner and rhside */

      total_number_u_3_points=number_primary_cells_i
				    *number_primary_cells_j*(number_primary_cells_k+1);

      momentum_rhside_u3.create(total_number_u_3_points);
      momentum_matrix_u3.create(number_matrix_connections,total_number_u_3_points); 
      preconditioner_matrix_M.create(total_number_u_3_points); 
      compressed_velocity_u3.create(total_number_u_3_points);
      
    /* build the system of linear equations */

      build_momentum_predictor_u3( momentum_matrix_u3, momentum_rhside_u3, level_set, 				
				    momentum_source_term_u_3, scaled_density_u3,		
				      u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,			
					number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k,			
					  actual_time_step_navier_stokes,		
					    mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					      rho_plus_over_rho_minus,			
						smoothing_distance_factor,			
						  rho_minus_over_mu_minus,			
						    mu_plus_over_mu_minus,			
						      boundary_faces,
							number_matrix_connections);
    
      /* build the preconditioner matrix M */
      
      build_preconditioner( number_primary_cells_i, number_primary_cells_j, number_primary_cells_k+1,
			    momentum_matrix_u3, preconditioner_matrix_M);
      
      /* compress the solution vector */
      
      compress_solution_velocity_u3(u_3_velocity_star, compressed_velocity_u3, 
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);

      /* solve the linear system for the new pressure */

      if(conjugate_gradient_method( number_primary_cells_i, number_primary_cells_j, number_primary_cells_k+1,
				  momentum_matrix_u3, preconditioner_matrix_M, momentum_rhside_u3, compressed_velocity_u3,
				    tolerance_velocity, iteration_number, relative_L2_norm_residual,
				      relative_Linfinity_norm_residual, maximum_iterations_allowed_velocity))
      {
	std::cout << " No convergence in linear solver for momentum equation u3.\n";
	exit(1);
      }
      else
      {
	std::cout << " The momentum equation u3 converged in " << iteration_number <<" iterations,\n";
	std::cout << " with relative L2 norm of residual " << relative_L2_norm_residual <<" \n";
	
	/* decompress the solution vector */
      
	decompress_solution_velocity_u3(u_3_velocity_star, compressed_velocity_u3, 
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
      }
    
     /* export the matrix to matlab for inspection if requested */

//       export_matrix_matlab(number_primary_cells_i, number_primary_cells_j, number_primary_cells_k+1,
// 			    4, momentum_matrix_u3, momentum_rhside_u3, compressed_velocity_u3, "velocity_u3");

    
    /* deallocate memory for momentum matrix and rhside */
    
      momentum_rhside_u3.destroy();
      momentum_matrix_u3.destroy();
      preconditioner_matrix_M.destroy();
      compressed_velocity_u3.destroy();
      
      }

#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
enum variable{velocity_u1, velocity_u2, velocity_u3, level_set, pressure};
enum boundary_conditions_type{dirichlet, neumann, periodic};
enum boundary_conditions_rule{constant, function};
enum cell_centerings{cell_centered, vertex_centered};


class boundary_variable
{
public:
  variable variable_name;
  boundary_conditions_type boundary_condition_type;
  boundary_conditions_rule boundary_condition_rule;
  cell_centerings cell_centering;
  double boundary_condition_value;
  boundary_variable(variable varname, boundary_conditions_type bound_type,
				     boundary_conditions_rule bound_rule,
				     cell_centerings  cell_cent,
					double bound_value );
  boundary_variable(variable varname);
};

class boundary_face
{
public:
    boundary_variable boundary_variables[5];
    boundary_face(void);
   
};

#include <string>
#include <sstream>
#include <fstream>
using namespace std;
/********************************************************************************/
/*  Function to solve the tentative velocity field u star, x2 components        */
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
  
      void solve_momentum_predictor_u2(
      double ***level_set, 				// level-set field
      double ***momentum_source_term_u_2, 		// complete source term for the momentum equation
							// in x3 direction=(-p,2+ g_2 +F2)
      double ***scaled_density_u2,			// scaled density for the controlvolumes
							// of the momentum equation in x2 direction
      double ***u_1_velocity_old, 			// velocity field at old time level x1 direction
      double ***u_2_velocity_old, 			// velocity field at old time level x2 direction
      double ***u_3_velocity_old,			// velocity field at old time level x3 direction
      double ***u_2_velocity_star,			// velocity field at star time level x1 direction
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
      double tolerance_velocity,	  		// the tolerance with which the system for the momentum predictor is solved	
      int maximum_iterations_allowed_velocity	// maximum number of iterations allowed for the
							// conjugate gradient method
							// for the boundary conditions 
	   )
      {


      double **double_Matrix(				// allocate memory for a two
		int number_primary_cells_i,		// dimensional array
		int number_primary_cells_j
	      );
      void   free_double_Matrix( 			// deallocate memory for a two
		double **doubleMatrix, 			// dimensional array
		int number_primary_cells_i
	      );
      double *double_Vector(				// allocate memory for a one
	      int number_primary_cells_i		// dimensional array
	      );
      void free_double_Vector(				// deallocate memory for a one
	      double *double_Vector			// dimensional array
	      );
				    
				    

                                         
    void build_momentum_predictor_u2(			// build system of linear
      double **momentum_matrix_u2,			// equations for momentum equation velocity
      double *momentum_rhside_u2,			// in x2 direction
      double ***level_set, 				
      double ***momentum_source_term_u_2, 		
      double ***scaled_density_u2,
      double ***u_1_velocity_old, 			
      double ***u_2_velocity_old, 			
      double ***u_3_velocity_old,			
      int number_primary_cells_i,			
      int number_primary_cells_j,			
      int number_primary_cells_k,			
      double actual_time_step_navier_stokes,		
      double mesh_width_x1,				
      double mesh_width_x2,				
      double mesh_width_x3,				
      double rho_plus_over_rho_minus,			
      double smoothing_distance_factor,			
      double rho_minus_over_mu_minus,			
      double mu_plus_over_mu_minus,			
      boundary_face boundary_faces[6],
      int number_matrix_connections
	 );
   void compress_solution_velocity_u2(			// compress solution vector
      double ***full_solution,		     
      double   *compressed_solution_vector,  
      int number_primary_cells_i,	     
      int number_primary_cells_j,	     
      int number_primary_cells_k	     
     );
   void decompress_solution_velocity_u2(		// decompress solution vector
      double ***full_solution,		     
      double   *compressed_solution_vector,  
      int number_primary_cells_i,	     
      int number_primary_cells_j,	     
      int number_primary_cells_k	     
     );
   int export_matrix_matlab(           			// export the matrix to matlab file
      int i_dimension,    		
      int j_dimension, 		
      int k_dimension, 		
      int number_matrix_connections,  
      double **matrix_A,		
      double *rhside_vector,		
      double *solution_vector,	
      string variable_name
	);
   void build_preconditioner(				// build incomplete choleski preconditioner
      int i_dimension,    	
      int j_dimension, 	
      int k_dimension, 	
      double **matrix_A,	
      double *preconditioner_matrix_M 	
	);
   int conjugate_gradient_method(			// solve linear system with conjugate gradient
      int i_dimension,   				// method (only SPD systems)
      int j_dimension,   		
      int k_dimension,   		
      double  **matrix_A,   		
      double  *preconditioner_matrix_M,  	
      double   *rhside_vector_b,   		
      double   *solution_vector_x, 		
      double   tolerance,	  		
      int     &iteration_number,  	
      double  &relative_L2_norm_residual, 
      double  &relative_Linfinity_norm_residual,
      int maximum_iterations_allowed	
        );


      double **momentum_matrix_u2;			// momentum matrix velocity x2 direction
      double *momentum_rhside_u2;			// momentum rhside velocity x2 direction
      double *preconditioner_matrix_M;		// preconditioner matrix (diagonal)
      double *compressed_velocity_u2;		// compressed solution vector
      double relative_L2_norm_residual; 	  	// the L2 norm of the residual, scaled with
							// the L2 norm of the right hand side
      double relative_Linfinity_norm_residual;  	// the L infinity norm of the residual
							// scaled with the maximum difference 
							// between two components of the residual
      int iteration_number;				// the number of iterations where the iterative
							// process was terminated
      int total_number_u_2_points;			// total number of grid points with u_2
  
      
      
      
    
    /* allocate memory for momentum matrix, preconditioner, rhside */
    /* and compressed solution vector */

      total_number_u_2_points=number_primary_cells_i
				    *(number_primary_cells_j+1)*number_primary_cells_k;

      momentum_rhside_u2=double_Vector(total_number_u_2_points);
      momentum_matrix_u2=double_Matrix(number_matrix_connections, total_number_u_2_points); 
      preconditioner_matrix_M=double_Vector(total_number_u_2_points); 
      compressed_velocity_u2=double_Vector(total_number_u_2_points);
      
    /* build the system of linear equations */

      build_momentum_predictor_u2( momentum_matrix_u2, momentum_rhside_u2, level_set, 				
				    momentum_source_term_u_2, scaled_density_u2,		
				      u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,			
					number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
					  actual_time_step_navier_stokes,		
					    mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					      rho_plus_over_rho_minus,			
						smoothing_distance_factor,			
						  rho_minus_over_mu_minus,			
						    mu_plus_over_mu_minus,			
						      boundary_faces,
							number_matrix_connections);
    
      /* build the preconditioner matrix M */
      
      build_preconditioner( number_primary_cells_i, number_primary_cells_j+1, number_primary_cells_k,
			    momentum_matrix_u2, preconditioner_matrix_M);
      
      /* compress the solution vector */
      
      compress_solution_velocity_u2(u_2_velocity_star, compressed_velocity_u2, 
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);

      /* solve the linear system for the new pressure */

      if(conjugate_gradient_method( number_primary_cells_i, number_primary_cells_j+1, number_primary_cells_k,
				  momentum_matrix_u2, preconditioner_matrix_M, momentum_rhside_u2, compressed_velocity_u2,
				    tolerance_velocity, iteration_number, relative_L2_norm_residual,
				      relative_Linfinity_norm_residual, maximum_iterations_allowed_velocity))
      {
	std::cout << " No convergence in linear solver for momentum equation u2.\n";
	exit(1);
      }
      else
      {
	std::cout << " The momentum equation u2 converged in " << iteration_number <<" iterations,\n";
	std::cout << " with relative L2 norm of residual " << relative_L2_norm_residual <<" \n";
	
	/* decompress the solution vector */
      
	decompress_solution_velocity_u2(u_2_velocity_star, compressed_velocity_u2, 
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
      }
     /* export the matrix to matlab for inspection if requested */

//       export_matrix_matlab(number_primary_cells_i, number_primary_cells_j+1, number_primary_cells_k,
// 			    4, momentum_matrix_u2, momentum_rhside_u2, compressed_velocity_u2, "velocity_u2");

    
    /* deallocate memory for momentum matrix and rhside */
    
      free_double_Vector(momentum_rhside_u2);
      free_double_Matrix(momentum_matrix_u2, number_matrix_connections);
      free_double_Vector(preconditioner_matrix_M);
      free_double_Vector(compressed_velocity_u2);
      }
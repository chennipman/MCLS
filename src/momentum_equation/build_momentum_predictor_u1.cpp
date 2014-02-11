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
         
      
/********************************************************************************/
/*  Function to build the whole system of linear equations for the velocity     */
/*  u1                                                                          */
/*                                                                              */
/*  Programmer: Duncan van der Heul                                             */
/*  Date      : 10-03-2013                                                      */
/*  Update    :                                                                 */
/********************************************************************************/
/* Notes                                                                        */
/* In the current implementation the generation of the matrix, right-hand-side  */
/* and the application of the boundary conditions are completely separated.     */
/* In this function all the required tasks are collected:                       */
/* - Build matrix without boundary conditions                                   */
/* - Build rhside without boundary conditions                                   */
/* - Apply boundary conditions to the right hand side, using the matrixi        */
/* - Apply boundary conditions to the matrix (matrix folding)                   */
/********************************************************************************/
    void build_momentum_predictor_u1(
      double **momentum_matrix_u1,                      // momentum matrix velocity x1 direction
      double *momentum_rhside_u1,                       // momentum rhside velocity x1 direction
      double ***level_set,                              // level-set field
      double ***momentum_source_term_u_1,               // complete source term for the momentum
                                                        // in x1 direction=(-p,1+ g_1 +F1)
      double ***scaled_density_u1,                      // scaled density for the controlvolumes
                                                        // of the momentum equation in x1 direction
      double ***u_1_velocity_old,                       // velocity field at old time level x1 direction
      double ***u_2_velocity_old,                       // velocity field at old time level x2 direction
      double ***u_3_velocity_old,                       // velocity field at old time level x3 direction

      int number_primary_cells_i,                       // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,                       // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,                       // number of primary (pressure) cells in x3 direction
      double actual_time_step_navier_stokes,            // time step used for level-set advection
                                                        // computed from all stability restrictions and
                                                        // possibly subscycling
      double mesh_width_x1,                             // grid spacing in x1 direction (uniform)
      double mesh_width_x2,                             // grid spacing in x2 direction (uniform)
      double mesh_width_x3,                             // grid spacing in x3 direction (uniform)
      double rho_plus_over_rho_minus,                   // ratio of the densities of the two phases
      double smoothing_distance_factor,                 // the smoothing distance is smoothing_distance_factor
                                                        // times the smallest mesh width
      double rho_minus_over_mu_minus,                   // this was the 'Reynolds' number
                                                        // in the original implementation of Sander
      double mu_plus_over_mu_minus,                     // ratio of the viscosities of both phases
      boundary_face boundary_faces[6],                  // array with all the information
                                                        // for the boundary conditions
      int number_matrix_connections                     // number of connections in momentum matrix
      
	 )


      {
    void build_momentum_matrix_u1(
      double **momentum_matrix_u1,			// build matrix for the momentum equation
      double ***level_set, 				// for the velocity in x1 direction
      double ***scaled_density_u1,
      int number_primary_cells_i,			
      int number_primary_cells_j,			
      int number_primary_cells_k,			
      double actual_time_step_navier_stokes,		
      double mesh_width_x1,				
      double mesh_width_x2,				
      double mesh_width_x3,				
      double rho_plus_over_rho_minus,			
      double rho_minus_over_mu_minus,			
      double mu_plus_over_mu_minus,			
      double smoothing_distance_factor			
       );
    void build_momentum_rhs_u1(
      double *momentum_rhside_u1,			// build rhside for the momentum equation
      double ***level_set, 				// for the velocity in x1 direction
      double ***scaled_density_u1,
      double ***momentum_source_term_u_1, 		 
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
      double mu_plus_over_mu_minus			 
       );
  void fold_momentum_rhside_u1(
      boundary_face boundary_faces[6],			// apply boundary conditions for the 
      double *momentum_rhside_u1,			// rhside of the momentum equation 
      double **momentum_matrix_u1,			// velocity in x1 direction
      int number_primary_cells_i,			 
      int number_primary_cells_j,			 
      int number_primary_cells_k,			 
      double mesh_width_x1,				 
      double mesh_width_x2,				 
      double mesh_width_x3				 
	   );
  void fold_momentum_matrix_u1(
      boundary_face boundary_faces[6],			// apply boundary conditions for the
      double **momentum_matrix_u1,			// matrix of the momentum equation
      int number_primary_cells_i,			// velocity in x1 direction
      int number_primary_cells_j,			 
      int number_primary_cells_k,
      int number_matrix_connections
	   );
      

    /* build momentum equation matrix for velocity u1 */
    /* without considering the boundary conditions */

      build_momentum_matrix_u1( momentum_matrix_u1, level_set, scaled_density_u1,
				 number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				  actual_time_step_navier_stokes,		
				    mesh_width_x1, mesh_width_x2, mesh_width_x3,				
				      rho_plus_over_rho_minus, rho_minus_over_mu_minus, mu_plus_over_mu_minus,			
					smoothing_distance_factor);

    /* build momentum equation rhside for velocity u1 */
    /* without considering the boundary conditions */

      build_momentum_rhs_u1( momentum_rhside_u1, level_set, scaled_density_u1,				
			      momentum_source_term_u_1, 		 
				u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,			 
				  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				    actual_time_step_navier_stokes,		 
				      mesh_width_x1, mesh_width_x2,mesh_width_x3,				 
					rho_plus_over_rho_minus,			 
					  smoothing_distance_factor,			 
					    rho_minus_over_mu_minus,			 
					      mu_plus_over_mu_minus);

    /* apply boundary conditions to the momentum equation rhside for velocity u1  */
    /* this has to be done first, because in that case the matrix can be utilised */

      fold_momentum_rhside_u1(boundary_faces,	
				momentum_rhside_u1, momentum_matrix_u1,			
				  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			 
				    mesh_width_x1, mesh_width_x2, mesh_width_x3);


    /* apply boundary conditions to the momentum equation matrix for velocity u1 */
 

      fold_momentum_matrix_u1(boundary_faces,			
				momentum_matrix_u1,		
				  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				    number_matrix_connections);




      }




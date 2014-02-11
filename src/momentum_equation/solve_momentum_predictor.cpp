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
/*  Function to solve the tentative velocity field u star, all components       */
/*                              						*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/*  										*/
/*  										*/
/*  										*/
/*  										*/
/********************************************************************************/
 void solve_momentum_predictor(
      double ***level_set, 				// level-set field
      double ***u_1_velocity_old, 			// velocity field at old time level x1 direction
      double ***u_2_velocity_old, 			// velocity field at old time level x2 direction
      double ***u_3_velocity_old,			// velocity field at old time level x3 direction
      double ***u_1_velocity_star, 			// velocity field at star time level x1 direction
      double ***u_2_velocity_star, 			// velocity field at star time level x2 direction
      double ***u_3_velocity_star,			// velocity field at star time level x3 direction
      double ***momentum_source_term_u_1, 		// complete source term for the momentum equation
							// in x1 direction=(-p,1+ g_1 +F1)
      double ***momentum_source_term_u_2, 		// complete source term for the momentum equation
							// in x1 direction=(-p,2+ g_2 +F2)
      double ***momentum_source_term_u_3, 		// complete source term for the momentum equation
							// in x1 direction=(-p,3+ g_3 +F3)
      double ***scaled_density_u1,			// scaled density for the controlvolumes
							// of the momentum equation in x1 direction
      double ***scaled_density_u2,			// scaled density for the controlvolumes
							// of the momentum equation in x2 direction
      double ***scaled_density_u3,			// scaled density for the controlvolumes
							// of the momentum equation in x3 direction
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      int number_matrix_connections,			// number of connections in momentum matrix				       
      double actual_time_step_navier_stokes,		// time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      double rho_plus_over_rho_minus,			// ratio of the densities of the two phases
      double smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
      double rho_minus_over_mu_minus,			// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      boundary_face boundary_faces[6],			// array with all the information
							// for the boundary conditions
      double tolerance_velocity,	  		// the tolerance with which the system for the momentum predictor is solved	
      int maximum_iterations_allowed_velocity		// maximum number of iterations allowed for the
							// conjugate gradient method
	
      )
      {
      void solve_momentum_predictor_u1(			// solve predictor equation for star velocity
      	    double ***level_set, 				// x1 component
      	    double ***momentum_source_term_u_1, 	
	    double ***scaled_density_u1,		
      	    double ***u_1_velocity_old, 		
      	    double ***u_2_velocity_old, 		
      	    double ***u_3_velocity_old,	
	    double ***u_1_velocity_star,
      	    int number_primary_cells_i,		
      	    int number_primary_cells_j,		
      	    int number_primary_cells_k,		
      	    int number_matrix_connections,		      
      	    double actual_time_step_navier_stokes,	
      	    double mesh_width_x1,			
      	    double mesh_width_x2,			
      	    double mesh_width_x3,			
      	    double rho_plus_over_rho_minus,		
      	    double smoothing_distance_factor,		
      	    double rho_minus_over_mu_minus,		
      	    double mu_plus_over_mu_minus,		
      	    boundary_face boundary_faces[6],
	    double tolerance_velocity,
	    int maximum_iterations_allowed_velocity
	   );
      void solve_momentum_predictor_u2(			// solve predictor equation for star velocity
      	    double ***level_set, 			// x1 component
      	    double ***momentum_source_term_u_2, 	
	    double ***scaled_density_u2,		
      	    double ***u_1_velocity_old, 		
      	    double ***u_2_velocity_old, 		
      	    double ***u_3_velocity_old,		
	    double ***u_2_velocity_star,
      	    int number_primary_cells_i,		
      	    int number_primary_cells_j,		
      	    int number_primary_cells_k,		
      	    int number_matrix_connections,		      
      	    double actual_time_step_navier_stokes,	
      	    double mesh_width_x1,			
      	    double mesh_width_x2,			
      	    double mesh_width_x3,			
      	    double rho_plus_over_rho_minus,		
      	    double smoothing_distance_factor,		
      	    double rho_minus_over_mu_minus,		
      	    double mu_plus_over_mu_minus,		
      	    boundary_face boundary_faces[6],
	    double tolerance_velocity,
	    int maximum_iterations_allowed_velocity
	   );
      void solve_momentum_predictor_u3(			// solve predictor equation for star velocity
      	    double ***level_set, 			// x1 component
      	    double ***momentum_source_term_u_3, 	
	    double ***scaled_density_u3,		
      	    double ***u_1_velocity_old, 		
      	    double ***u_2_velocity_old, 		
      	    double ***u_3_velocity_old,		
	    double ***u_3_velocity_star,
       	    int number_primary_cells_i,		
      	    int number_primary_cells_j,		
      	    int number_primary_cells_k,		
      	    int number_matrix_connections,		      
      	    double actual_time_step_navier_stokes,	
      	    double mesh_width_x1,			
      	    double mesh_width_x2,			
      	    double mesh_width_x3,			
      	    double rho_plus_over_rho_minus,		
      	    double smoothing_distance_factor,		
      	    double rho_minus_over_mu_minus,		
      	    double mu_plus_over_mu_minus,		
      	    boundary_face boundary_faces[6],
	    double tolerance_velocity,
	    int maximum_iterations_allowed_velocity
	   );
      void apply_boundary_conditions_velocity(		// apply the boundary conditions to the 
	   boundary_face boundary_faces[6],		// velocity field
	   double ***u_1_velocity, 			 
	   double ***u_2_velocity, 			 
	   double ***u_3_velocity, 			 
	   double mesh_width_x1,				 
	   double mesh_width_x2,				 
	   double mesh_width_x3,				 
	   int number_primary_cells_i,			 
	   int number_primary_cells_j,			 
	   int number_primary_cells_k			 
           );
      void output_predictor_velocityfields(		// visualize the predictor velocity field
       	   double ***u_1_velocity_new, 		
       	   double ***u_2_velocity_new, 		
       	   double ***u_3_velocity_new,		
       	   int number_primary_cells_i,		
       	   int number_primary_cells_j,		
       	   int number_primary_cells_k,		
       	   double mesh_width_x1,			
       	   double mesh_width_x2,			
       	   double mesh_width_x3			
           );
       
      
      /* solve the momentum predictor equation for the velocity u star, x1 component */
      		std::cout<<"solve momentum predictor u1 "<< " \n";

      solve_momentum_predictor_u1( level_set, momentum_source_term_u_1, 	scaled_density_u1,
				  u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,	u_1_velocity_star,		
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
				      number_matrix_connections, actual_time_step_navier_stokes,	
					mesh_width_x1, mesh_width_x2, mesh_width_x3,			
					  rho_plus_over_rho_minus, smoothing_distance_factor,		
					    rho_minus_over_mu_minus, mu_plus_over_mu_minus,		
					      boundary_faces,
					      tolerance_velocity, maximum_iterations_allowed_velocity);
      
      /* solve the momentum predictor equation for the velocity u star, x2 component */
      		std::cout<<"solve momentum predictor u2 "<< " \n";

      solve_momentum_predictor_u2( level_set, momentum_source_term_u_2, 	scaled_density_u2,
				  u_1_velocity_old, u_2_velocity_old, u_3_velocity_old, u_2_velocity_star,		
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
				      number_matrix_connections, actual_time_step_navier_stokes,	
					mesh_width_x1, mesh_width_x2, mesh_width_x3,			
					  rho_plus_over_rho_minus, smoothing_distance_factor,		
					    rho_minus_over_mu_minus, mu_plus_over_mu_minus,		
					      boundary_faces,
					      tolerance_velocity, maximum_iterations_allowed_velocity);
      
      /* solve the momentum predictor equation for the velocity u star, x3 component */
            		std::cout<<"solve momentum predictor u3 "<< " \n";

      solve_momentum_predictor_u3( level_set, momentum_source_term_u_3, 	scaled_density_u3,
				  u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,	u_3_velocity_star,		
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
				      number_matrix_connections, actual_time_step_navier_stokes,	
					mesh_width_x1, mesh_width_x2, mesh_width_x3,			
					  rho_plus_over_rho_minus, smoothing_distance_factor,		
					    rho_minus_over_mu_minus, mu_plus_over_mu_minus,		
					      boundary_faces,
					      tolerance_velocity, maximum_iterations_allowed_velocity);
   
      /* apply the boundary conditions to u star, and fill the virtual cells */

      apply_boundary_conditions_velocity(boundary_faces,		
					  u_1_velocity_star, u_2_velocity_star, u_3_velocity_star, 			 
					    mesh_width_x1, mesh_width_x2, mesh_width_x3,				 
					      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);			 
   
      /* visualize the predictor velocity field */
   
//       output_predictor_velocityfields(u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,		
// 					   number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
// 					     mesh_width_x1, mesh_width_x2, mesh_width_x3);
      
      
      }      
      
      

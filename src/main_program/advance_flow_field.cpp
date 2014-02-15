#include "../headers/array.h"
class vector
{
public:
  double u1,u2,u3;
};
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
    void advance_flow_field(
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
   
      void solve_momentum_predictor(			// solve momentum predictor equation
      	 Array3<double> level_set, 			
      	 Array3<double> u_1_velocity_old, 	
      	 Array3<double> u_2_velocity_old, 	
      	 Array3<double> u_3_velocity_old,	
      	 Array3<double> u_1_velocity_star, 	
      	 Array3<double> u_2_velocity_star, 	
      	 Array3<double> u_3_velocity_star,	
      	 Array3<double> momentum_source_term_u_1, 
      	 Array3<double> momentum_source_term_u_2, 
      	 Array3<double> momentum_source_term_u_3, 
      	 Array3<double> scaled_density_u1,			
	 Array3<double> scaled_density_u2,			
	 Array3<double> scaled_density_u3,			
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
     
      void solve_momentum_corrector(				// solve momentum corrector equation
      	   Array3<double> level_set,		
	   Array3<double> pressure,		
	   Array3<double> momentum_source_term_u_1,
	   Array3<double> momentum_source_term_u_2,  
	   Array3<double> momentum_source_term_u_3,
	   Array3<double> surface_tension_body_force_x1, 
	   Array3<double> surface_tension_body_force_x2, 
	   Array3<double> surface_tension_body_force_x3,
	   Array3<double> scaled_density_u1,
	   Array3<double> scaled_density_u2,
	   Array3<double> scaled_density_u3,
	   Array3<double> u_1_velocity_star, 	       
	   Array3<double> u_2_velocity_star, 	       
	   Array3<double> u_3_velocity_star,	       
	   double mesh_width_x1,		       
	   double mesh_width_x2,		       
	   double mesh_width_x3,		       
          int number_primary_cells_i,	       
          int number_primary_cells_j,	       
          int number_primary_cells_k,	       
	   vector gravity,			
	   double tolerance_pressure,	  
          double actual_time_step_navier_stokes,   
          double rho_plus_over_rho_minus,	       
          int continuous_surface_force_model,      
          int source_terms_in_momentum_predictor,  
	   int maximum_iterations_allowed_pressure,
          boundary_face boundary_faces[6]
      );
 
      void shift_velocity_field(			     // shift the velocity fields
          Array3<double> u_1_velocity_old, 		             // u new -> u old
          Array3<double> u_2_velocity_old, 		             // u star -> u new
          Array3<double> u_3_velocity_old,		
          Array3<double> u_1_velocity_new, 		
          Array3<double> u_2_velocity_new, 		
          Array3<double> u_3_velocity_new,		
          Array3<double> u_1_velocity_star,		
          Array3<double> u_2_velocity_star,		
          Array3<double> u_3_velocity_star,		
      	  int number_primary_cells_i,	
      	  int number_primary_cells_j,	
      	  int number_primary_cells_k	
	 );	
    
      
       Array3<double> u_1_velocity_star; 		     // velocity field at star time level x1 direction   
       Array3<double> u_2_velocity_star; 		     // velocity field at star time level x2 direction     
       Array3<double> u_3_velocity_star;		     // velocity field at star time level x3 direction
       
       
       
    /* allocate memory for tentative velocity field u star */
    
      u_1_velocity_star.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_star.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_star.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
  
    /* solve momentum predictor equation */
    /* compute a new velocity field u star, that is not divergence free */
    
      solve_momentum_predictor(	level_set, u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,			
				  u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,			
				    momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3, 	
				     scaled_density_u1, scaled_density_u2, scaled_density_u3,
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
					number_matrix_connections, actual_time_step_navier_stokes,		
					  mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					    rho_plus_over_rho_minus, smoothing_distance_factor, rho_minus_over_mu_minus,			
					      mu_plus_over_mu_minus, boundary_faces,
					      tolerance_velocity, maximum_iterations_allowed_velocity);

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
	
     }
	
	
	
	
	
	
	
	
	

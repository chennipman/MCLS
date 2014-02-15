#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
enum variable{velocity_u1, velocity_u2, velocity_u3, level_set, pressure};
enum boundary_conditions_type{dirichlet, neumann, periodic};
enum boundary_conditions_rule{constant, function};
enum cell_centerings{cell_centered, vertex_centered};
enum geometry{bubbly_flow, wavy_flow};


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
class coordinate
{
public:
  double x1,x2,x3;
  coordinate(double xx1=0, double xx2=0, double xx3=0){x1=xx1;x2=xx2;x3=xx3;}
};
class bubble
{
public:
  double principle_axis_x1;
  double principle_axis_x2;
  double principle_axis_x3;
  int label;
  coordinate center_location;
  bubble(int number, coordinate bubble_center, double bubble_radius);
};
class surface
{
public:
  int active;
  int orientation;
  double height;
};

class vector
{
public:
  double u1,u2,u3;
};
class restart_parameters
{
public:
      int start_from_restart_file;		
      int write_solution_to_restart_file;
      string name_restart_file_to_write;
      string name_restart_file_to_read;
      restart_parameters(void);
};

/********************************************************************************/
/********************************************************************************/
/*  initialize_all_variables						     	       */
/*  Description : perform initialization of all variables	      		       */
/*  Programmer	: Duncan van der Heul       					       */
/*  Date	: 03-01       							       */
/*  Update	:        							       */
/********************************************************************************/

void initialize_all_variables(
      Array3<double> u_1_velocity_new, 			       // velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new, 			       // velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,			       // velocity field at new time level x3 direction
      Array3<double> u_1_velocity_old, 			       // velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			       // velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			       // velocity field at old time level x3 direction
      Array3<double> pressure,				       // pressure
      Array3<double> level_set_old,			       // level-set field at old time level
      Array3<double> volume_of_fluid,			       // volume of fluid field
      Array3<double> curvature,				       // interface curvature
      Array3<double> smoothed_curvature,			       // smoothed interface curvature
      Array3<double> momentum_source_term_u_1,                 // source term of the momentum equation in x1 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2,                 // source term of the momentum equation in x2 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3,                 // source term of the momentum equation in x3 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x1,            // source term of the momentum equation in x1 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x2,            // source term of the momentum equation in x2 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x3,            // source term of the momentum equation in x3 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> scaled_density_u1,                        // scaled density for the controlvolumes
                                                          // of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,                        // scaled density for the controlvolumes
                                                          // of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,                        // scaled density for the controlvolumes
                                                          // of the momentum equation in x3 direction
      boundary_face boundary_faces[6],		       // array with all the information
							       // for the boundary conditions 
      geometry flow_type,				       // the kind of initial condition that has to be applied
      bubble *the_bubbles,				       // array with definitions of the bubbles
      int number_of_bubbles,				       // number of bubbles in the initial condition
      surface *the_free_surfaces,			       // array with definitions of the free surfaces
      int number_of_free_surfaces, 			       // number of bubbles in the domain (<10)
      double lower_bound_derivatives,			// lower bound for the first partial derivatives
							       // to consider it a limiting case of vanishing
							       // partial derivatives
      double mesh_width_x1,				       // grid spacing in x1 direction (uniform)
      double mesh_width_x2,				       // grid spacing in x2 direction (uniform)
      double mesh_width_x3,				       // grid spacing in x3 direction (uniform)
      int number_primary_cells_i,			       // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			       // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			       // number of primary (pressure) cells in x3 direction
      vector initial_velocity,			       // initial velocity field at t=0
      restart_parameters my_restart_parameters, 	       // all parameters for reading/writing restart files
      double &start_time_simulation,			       // starting time for the simulation
      int apply_curvature_smoothing,			       // =1, apply curvature smoothing
							       // =0, use unsmoothed curvature
      int number_curvature_smoothing_steps,		       // number of iterations applied in the
							       // curvature smoothing algorithm
      int apply_curvature_smoothing_filter, 	       // =1, apply curvature smoothing filter
	    						       // =0, do not apply curvature smoothing filter		
      double smoothing_distance_factor,		       // the heaviside function is smoothed over
							       // an interval of width 2*smoothing_distance
      int debugging_mode,				       // =1, program is run in 
      							       // debugging mode with additional output and checks
      							       // =0, progam is run in production mode
      vector gravity,                                     // gravitational acceleration vector       double tolerance,                            // the tolerance with which the system is solved
      double tolerance,                                   // the tolerance with which the system is solved
      double actual_time_step_navier_stokes,              // actual time step for Navier-Stokes solution algorithm
      double rho_plus_over_rho_minus,                     // ratio of density where (level set >0) and
                                                          // density where (level set < 0)
      double sigma_over_rho_minus,                        // sigma / rho_minus (scaled sigma)
      int continuous_surface_force_model,                 // =1, the continuous surface force model is applied
                                                          // =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor,             // =1, the source terms are applied in the momentum predictor
                                                          // equation
                                                          // =0, the source terms are applied in the momentum corrector
                                                          // equation
      int maximum_iterations_allowed                      // maximum number of iterations allowed for the
                                                          // conjugate gradient method
      )
      {
     void initialize_flow_field(			// initialize the flow field: velocity and pressure
	  Array3<double> u_1_velocity_new, 			
	  Array3<double> u_2_velocity_new, 			
	  Array3<double> u_3_velocity_new,			
	  Array3<double> u_1_velocity_old, 			
	  Array3<double> u_2_velocity_old, 			
	  Array3<double> u_3_velocity_old,			
	  Array3<double> pressure,				
         Array3<double> level_set,                         
         Array3<double> momentum_source_term_u_1,
         Array3<double> momentum_source_term_u_2,
         Array3<double> momentum_source_term_u_3,
         Array3<double> surface_tension_body_force_x1,
         Array3<double> surface_tension_body_force_x2,
         Array3<double> surface_tension_body_force_x3,
         Array3<double> scaled_density_u1,
         Array3<double> scaled_density_u2,
         Array3<double> scaled_density_u3,
	  boundary_face boundary_faces[6],			
	  double mesh_width_x1,				
	  double mesh_width_x2,				
	  double mesh_width_x3,				
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k,
	  vector initial_velocity,	
         vector gravity,                                                       
         double tolerance,
         double actual_time_step_navier_stokes,
         double rho_plus_over_rho_minus,
         double sigma_over_rho_minus,
         int continuous_surface_force_model,
         int source_terms_in_momentum_predictor,
         int maximum_iterations_allowed
                                                   
     );

     void initialize_interface(				// intialize the interface: level-set and volume of fluid
	  geometry flow_type,			
	  bubble *the_bubbles,			
	  int number_of_bubbles,		
	  surface *the_free_surfaces,		
	  int number_of_free_surfaces, 		
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k,		
	  double mesh_width_x1,			
	  double mesh_width_x2,			
	  double mesh_width_x3,			
	  double lower_bound_derivatives,	
	  Array3<double> level_set,			
	  Array3<double> volume_of_fluid,	
	  Array3<double> curvature,
	  Array3<double> smoothed_curvature,
      	  int apply_curvature_smoothing,		
         int number_curvature_smoothing_steps,	
         int apply_curvature_smoothing_filter,
	  double smoothing_distance_factor,
	  int debugging_mode
	  );
     void read_restart_file(            		       // read restart file
	  double &time_of_restart,
	  Array3<double> level_set_new, 			
							
	  Array3<double> volume_of_fluid, 			
	  Array3<double> u_1_velocity_new, 			
	  Array3<double> u_2_velocity_new, 			
	  Array3<double> u_3_velocity_new,			
	  Array3<double> pressure,				
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k,			
	  restart_parameters my_restart_parameters
	  );
     void initialize_coupling(                          // initialize all variables
           Array3<double> level_set,                         // related to the coupling between
           Array3<double> pressure,                          // flow and the interface
           Array3<double> curvature,
           Array3<double> unsmoothed_curvature,
           Array3<double> surface_tension_body_force_x1,
           Array3<double> surface_tension_body_force_x2,
           Array3<double> surface_tension_body_force_x3,
           Array3<double> momentum_source_term_u_1,
           Array3<double> momentum_source_term_u_2,
           Array3<double> momentum_source_term_u_3,
           Array3<double> scaled_density_u1,
           Array3<double> scaled_density_u2,
           Array3<double> scaled_density_u3,
           int number_primary_cells_i,
           int number_primary_cells_j,
           int number_primary_cells_k,
           double mesh_width_x1,
           double mesh_width_x2,
           double mesh_width_x3,
           double rho_plus_over_rho_minus,
           double actual_time_step_navier_stokes,
           vector gravity,
           double sigma_over_rho_minus,
           int apply_curvature_smoothing,
           int number_curvature_smoothing_steps,
           int apply_curvature_smoothing_filter,
           double smoothing_distance_factor,
           double lower_bound_derivatives
      );
 
      if(my_restart_parameters.start_from_restart_file)
      {
	
	/* a restart file is available, the computation starts from an old solution */
	/* and the start time is set to the time of the restart file */
	
	  read_restart_file(start_time_simulation, level_set_old, volume_of_fluid, 			
			      u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,			
				pressure,				
				  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				    my_restart_parameters);
	   
      }
      else
      {
      

      /* initialize the interface */
      

      
	  initialize_interface( flow_type, the_bubbles, number_of_bubbles, the_free_surfaces, number_of_free_surfaces, 		
			      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
				mesh_width_x1, mesh_width_x2, mesh_width_x3, lower_bound_derivatives,	
				  level_set_old, volume_of_fluid, curvature, smoothed_curvature, apply_curvature_smoothing,
			    		number_curvature_smoothing_steps, apply_curvature_smoothing_filter, smoothing_distance_factor,
				  debugging_mode
 				);

     /* initialize the coupling between interface and flow field */


         initialize_coupling( level_set_old, pressure, smoothed_curvature, curvature,
                     surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                        scaled_density_u1, scaled_density_u2, scaled_density_u3,
                         number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                           mesh_width_x1, mesh_width_x2, mesh_width_x3,
                            rho_plus_over_rho_minus, actual_time_step_navier_stokes,
                              gravity, sigma_over_rho_minus,
                                  apply_curvature_smoothing, number_curvature_smoothing_steps,
                                   apply_curvature_smoothing_filter,
                                     smoothing_distance_factor, lower_bound_derivatives);

      
      /* initialize the flow field */

         initialize_flow_field( u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,
                              u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,
                                pressure, level_set_old,
                                  momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                                   surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                                     scaled_density_u1, scaled_density_u2, scaled_density_u3,
                                      boundary_faces, mesh_width_x1, mesh_width_x2, mesh_width_x3,
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                                          initial_velocity, gravity, tolerance, actual_time_step_navier_stokes,
                                            rho_plus_over_rho_minus, sigma_over_rho_minus, continuous_surface_force_model,
                                             source_terms_in_momentum_predictor, maximum_iterations_allowed);


      }
}  

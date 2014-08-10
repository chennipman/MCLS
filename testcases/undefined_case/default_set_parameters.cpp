#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/********************************************************************************/
/*  Function to set the computation parameters   (undefined_case)                               */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function should be replaced by a function that reads all parameters     */
/* from an input file.                                                           */
/********************************************************************************/
EXPORT void set_parameters(
      int &number_primary_cells_i,				// number of primary (pressure) cells in x1 direction
      int &number_primary_cells_j,				// number of primary (pressure) cells in x2 direction
      int &number_primary_cells_k,				// number of primary (pressure) cells in x3 direction
      double &actual_time_step_level_set,			// time step used for level-set advection
								// computed from all stability restrictions and 
								// possibly subscycling
      double &mesh_width_x1,					// grid spacing in x1 direction (uniform)
      double &mesh_width_x2,					// grid spacing in x2 direction (uniform)
      double &mesh_width_x3,					// grid spacing in x3 direction (uniform)
      int &apply_mass_distribution_algorithm,   		// =1:apply the mass redistribution algorithm
								// to avoid numerical vapour
      int &apply_mass_conservation_correction, 		        // =1:apply the mass conservation correction
								// to the level-set field
      double &volume_of_fluid_tolerance,			// tolerance in volume of fluid field
								// for admissable values
      double &lower_bound_derivatives,			        // lower bound for the first partial derivatives
								// to consider it a limiting case of vanishing
								// partial derivatives
      int &number_vof_2_level_set_iterations,		        // number of OUTER iterations in the conversion 
								// from volume of fluid to level-set
      int &number_iterations_ridder,				// maximum number of iterations allowed in the
								// nonlinear root finding algorithm
								// this is the number of INNER iterations
      double &vof_2_level_set_tolerance,			// tolerance in the conversion from volume
								// of fluid value to level-set value
      double &cfl_number_reinitialization,			// courant-friedrichs-lewy number for the reinitialization
								// equation
      double &cfl_number_navier_stokes,			        // courant-friedrichs-lewy number for the navier-stokes
								// equations
      int &maximum_reinitialization_steps,			// maximum number of time steps in the reinitialization
								// algorithm
      double &tolerance_reinitialization,			// stop the reinitialization when the infinite norm of
								// the pseudo time derivative has fallen below this
								// tolerance value
      int &number_matrix_connections,			        // number of connections in momentum matrix				       
      vector &gravity,					        // gravitational acceleration vector 
      double &actual_time_step_navier_stokes,		        // time step used for level-set advection
								// computed from all stability restrictions and 
								// possibly subscycling
      double &rho_plus_over_rho_minus,			        // ratio of the densities of the two phases
      double &smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
								// times the smallest mesh width
      double &rho_minus_over_mu_minus,			        // this was the 'Reynolds' number
								// in the original implementation of Sander
      double &mu_plus_over_mu_minus,				// ratio of the viscosities of both phases
      double &tolerance_pressure,	  			// the tolerance with which the system for the pressure is solved	
      int &maximum_iterations_allowed_pressure,		        // maximum number of iterations allowed for the
								// conjugate gradient method
      double &tolerance_velocity,	  			// the tolerance with which the system for the momentum predictor is solved	
      int &maximum_iterations_allowed_velocity,		        // maximum number of iterations allowed for the
								// conjugate gradient method
      int &continuous_surface_force_model,     	 	        // =1, the continuous surface force model is applied
					        		// =0, the exact interface boundary conditions are applied
      int &source_terms_in_momentum_predictor,  		// =1, the source terms are applied in the momentum predictor
					        		// equation
					        		// =0, the source terms are applied in the momentum corrector
      double &sigma_over_rho_minus,				// sigma / rho_minus (scaled sigma)
      double &time_step_restriction_global,    		        // upper bound on the time step, set by user
      int &fixed_time_step,					// =1 the time step is not restricted by 
								// physical constraints 
								// =0 the time step is adapted to the physical
								// constraints
      int &apply_curvature_smoothing,			        // =1, apply curvature smoothing
								// =0, use unsmoothed curvature
      int &number_curvature_smoothing_steps,		        // number of iterations applied in the
								// curvature smoothing algorithm
      int &apply_curvature_smoothing_filter,		        // =1, apply curvature smoothing filter
								// =0, do not apply curvature smoothing filter
      int &number_of_subcycles,				        // for every time step for the Navier-Stokes
								// equations number_of_subcycles steps are taken
								// with the interface evolution equations
      int &vtk_output,					        // =1, write output in vtk format
								// =0, skip output in vtk format
      int &tecplot_output,					// =1, write output in tecplot format
								// =0, skip output in tecplot format
      int &pure_output,						// =1, write output in pure format
								// =0, skip output in pure format
      double &time_interval_for_output,			        // interval in time for which the solution is written to file
      double &domain_size_x1,					// length of edges of domain in x1 direction 
      double &domain_size_x2,					// length of edges of domain in x1 direction 
      double &domain_size_x3,					// length of edges of domain in x1 direction 
      geometry &flow_type,					// the kind of initial condition that has to be applied
      bubble *the_bubbles,					// array with definitions of the bubbles
      int &number_of_bubbles,					// number of bubbles in the initial condition
      surface *the_free_surfaces,				// array with definitions of the free surfaces
      int &number_of_free_surfaces, 				// number of bubbles in the domain (<10)
      double &start_time_simulation,				// starting time for the simulation
      double &end_time_simulation,				// end time for the simulation
      time_stepping_methods &time_stepping_method, 		// time scheme {none, explicit_euler, imex, runge_kutta, two_pres_solve, two_pres_solve_output};   
      vector &initial_velocity,				        // initial velocity at t=0
      restart_parameters &my_restart_parameters,		// all parameters for reading/writing restart files
      int &maximum_number_mass_redistribution_iterations,       // number of iterations allowed to make
								// the volume of fluid field valid
								// these are the sweeps on the vof error
      double &time_step_mass_redistribution,		        // time step for the mass redistribution
								// algorithm
      double &redistribution_vof_tolerance, 		        // threshold value of time-derivative 
      int &number_of_phases,					// number of phases in the domain:
								// 1: single phase flow
								// 2: two phase flow
      int &debugging_mode,					// =1, program is run in 
      								// debugging mode with additional output and checks
      								// =0, progam is run in production mode
      double &time_interval_for_reinitialization                // time interval between reinitialization
                                                                // of level-set field
                                        
					
      )
      {
      
/*-------------------------------------------------------------------------------------------------*/


      /* general settings */
      debugging_mode					= 0;


      /* computational domain geometry information */
      /* the origin is always chosen in the lower, left corner */
      domain_size_x1                                            = 2.0;
      domain_size_x2                                            = 2.0;
      domain_size_x3                                            = 1.0/20;


      /* settings for modelling */
      number_of_phases				= 2;


      /* settings for time-stepping */
      cfl_number_navier_stokes			        = 1.0;
      time_step_restriction_global			= 0.1;      
      actual_time_step_navier_stokes    		= time_step_restriction_global;
      time_interval_for_output			        = time_step_restriction_global*1;
      time_interval_for_reinitialization                = time_step_restriction_global*1e10; 
      number_of_subcycles				= 1;
      actual_time_step_level_set			= actual_time_step_navier_stokes/number_of_subcycles;
      fixed_time_step					= 1;	
      start_time_simulation				= 0.0;
      end_time_simulation				= 1.0;
      time_stepping_method 				= explicit_euler;


      /* settings for restart from solution file and solution file writing */
      my_restart_parameters.start_from_restart_file		        = 0;		
      my_restart_parameters.write_solution_to_restart_file		= 0;
      my_restart_parameters.name_restart_file_to_write		        = "restart_file_mcls_out";
      my_restart_parameters.name_restart_file_to_read		        = "restart_file_mcls_in";


      /* grid parameters */
      number_primary_cells_i                                    = 20;
      number_primary_cells_j                                    = 20;
      number_primary_cells_k                                    = 1;	
      mesh_width_x1=domain_size_x1/number_primary_cells_i;		
      mesh_width_x2=domain_size_x2/number_primary_cells_j;			
      mesh_width_x3=domain_size_x3/number_primary_cells_k;		


      /* interface handling parameters */
      apply_mass_distribution_algorithm  				= 1;    
      apply_mass_conservation_correction 				= 1;    
      volume_of_fluid_tolerance		 			        = 1e-6;
      lower_bound_derivatives		 				= 0.00000001;
      number_vof_2_level_set_iterations  				= 25;	
      number_iterations_ridder		 			        = 100;	
      vof_2_level_set_tolerance		 			        = 1e-6;	
      cfl_number_reinitialization	 				= 0.5;
      maximum_reinitialization_steps	 				= 400;
      tolerance_reinitialization	 				= 0.0001;
      apply_curvature_smoothing	        			        =  0;
      number_curvature_smoothing_steps   				= 4;
      apply_curvature_smoothing_filter	 			        =  0;
      maximum_number_mass_redistribution_iterations		        = 1;
      time_step_mass_redistribution					= 0.01;	
      redistribution_vof_tolerance 					= 0.0001;		


      /* initial condition */
      flow_type				 			        = bubbly_flow;
      the_bubbles[0].principle_axis_x1	 			        = 0.01;
      the_bubbles[0].principle_axis_x2	 			        = 0.01;
      the_bubbles[0].principle_axis_x3	 			        = 0.01;
      the_bubbles[0].center_location.x1  				= 0.5*domain_size_x1;
      the_bubbles[0].center_location.x2  				= 0.5*domain_size_x2;
      the_bubbles[0].center_location.x3  				= 0.025;
      number_of_bubbles			 			        = 0;
      number_of_free_surfaces		 				= 0;	
      the_free_surfaces[0].active	 				= 1;
      the_free_surfaces[0].orientation   				= 3;
      the_free_surfaces[0].height	 				= 0.49;
      initial_velocity.u1		 				= 0.0;
      initial_velocity.u2		 				= 0.0;
      initial_velocity.u3		 				= 0.0;
     

      /* material properties */
      rho_minus_over_mu_minus		 	=    100.0;
      mu_plus_over_mu_minus		 	=    1.0;
      rho_plus_over_rho_minus		 	=    1.0;
      sigma_over_rho_minus                      =     0.0;
      smoothing_distance_factor		        =1.5;	


      /* physics */
      gravity.u1			 =0.0;				
      gravity.u2			 =0.0;				
      gravity.u3			 =0.0;				
      
      /* linear solvers and matrices*/

      tolerance_pressure		 		= 1e-6;	 
      tolerance_velocity		 		= 1e-5;
      maximum_iterations_allowed_pressure		=      400;	
      maximum_iterations_allowed_velocity		=      100;
      number_matrix_connections		 	        =        7;		       

      /* coupling between interface and flow model */

      continuous_surface_force_model	 	 = 1;      
      source_terms_in_momentum_predictor	 = 1;


      /* output */
      vtk_output			 = 1;
      tecplot_output			 = 0;
      pure_output			 = 1;

      }

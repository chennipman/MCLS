#include "../headers/array.h"
/*********************************************/
/* include the necessary standard functions  */
/*********************************************/
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/*********************************************/
/* include the necessary class  definitions  */
/*********************************************/
#include "header_clsses.h"


/*++++++++++++++++++++++++++++++++++++*/ 
/* Primary variables for Navier-Stokes*/
/*++++++++++++++++++++++++++++++++++++*/ 

Array3<double> u_1_velocity_old; 		// 3-array with x1-velocity at old time level at u1 points, old time level
Array3<double> u_2_velocity_old; 		// 3-array with x2-velocity at old time level at u2 points, new time level
Array3<double> u_3_velocity_old; 		// 3-array with x3-velocity at old time level at u3 points, old time level
Array3<double> u_1_velocity_new; 		// 3-array with x1-velocity at new time level at u1 points, new time level
Array3<double> u_2_velocity_new; 		// 3-array with x2-velocity at new time level at u2 points, old time level
Array3<double> u_3_velocity_new; 		// 3-array with x3-velocity at new time level at u3 points, new time level
Array3<double> pressure;	    		// 3-array with pressure field at cell centers, old time level 


/*++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for interface capturing  */
/*++++++++++++++++++++++++++++++++++++++++++++*/  
Array3<double> volume_of_fluid;		// 3-array with volume_of_fluid at cell centers
Array3<double> level_set_old;        	// 3-array with level set field at cell centers, old time level
Array3<double> level_set_new;        	// 3-array with level set field at cell centers, new time level
Array3<double> curvature; 			// 3-array with interface curvature
Array3<double> unsmoothed_curvature;	// 3-array with interface curvature without smoothing

/*++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Primary variables for interface to flow coupling */
/*++++++++++++++++++++++++++++++++++++++++++++++++++*/

Array3<double> momentum_source_term_u_1;       // 3-array with momentum equation source term x1 component
Array3<double> momentum_source_term_u_2;       // 3-array with momentum equation source term x2 component
Array3<double> momentum_source_term_u_3;       // 3-array with momentum equation source term x3 component
Array3<double> surface_tension_body_force_x1;  // 3-array with momentum equation body force x1 component
Array3<double> surface_tension_body_force_x2;  // 3-array with momentum equation body force x2 component
Array3<double> surface_tension_body_force_x3;  // 3-array with momentum equation body force x3 component
Array3<double> scaled_density_u1;              // 3-array with scaled density for u1 control volume
Array3<double> scaled_density_u2;              // 3-array with scaled density for u2 control volume
Array3<double> scaled_density_u3;              // 3-array with scaled density for u3 control volume



/*++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for boundary conditions  */
/*++++++++++++++++++++++++++++++++++++++++++++*/  

boundary_face boundary_faces[6];		// array with all the information
						// for the boundary conditions 
						
/*++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for geometry and grid    */
/*++++++++++++++++++++++++++++++++++++++++++++*/  
						
int number_primary_cells_i;			// number of primary (pressure) cells in x1 direction
int number_primary_cells_j;			// number of primary (pressure) cells in x2 direction
int number_primary_cells_k;			// number of primary (pressure) cells in x3 direction
double mesh_width_x1;				// grid spacing in x1 direction (uniform)
double mesh_width_x2;				// grid spacing in x2 direction (uniform)
double mesh_width_x3;				// grid spacing in x3 direction (uniform)
double domain_size_x1;			// length of edges of domain in x1 direction 
double domain_size_x2;			// length of edges of domain in x1 direction 
double domain_size_x3;			// length of edges of domain in x1 direction 






/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* General variables 						    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 

int debugging_mode; 			// =1, program is run in debugging mode with additional output
					// and checks
					// =0, program is run in production mode

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for restarting from existing solution    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 

restart_parameters my_restart_parameters;

/*++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for temporal discretisation    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++*/  

double vof_2_level_set_tolerance;		// tolerance in the conversion from volume
						// of fluid value to level-set value
double cfl_number_reinitialization;		// courant-friedrichs-lewy number for the reinitialization
						// equation
double cfl_number_navier_stokes;		// courant-friedrichs-lewy number for the navier-stokes
						// equations
double actual_time_step_level_set;		// time step used for level-set advection
						// computed from all stability restrictions and 
						// possibly subscycling
double actual_time_step_navier_stokes;	// time step used for level-set advection
						// computed from all stability restrictions and 
						// possibly subscycling
double time_step_restriction_global;    	// upper bound on the time step, set by user
int fixed_time_step;				// =1 the time step is not restricted by 
						// physical constraints 
						// =0 the time step is adapted to the physical
						// constraints
int number_of_subcycles;			// for every time step for the Navier-Stokes
						// equations number_of_subcycles steps are taken
						// with the interface evolution equations
double start_time_simulation;			// starting time for the simulation
double end_time_simulation;			// end time for the simulation
double time_interval_for_reinitialization; 	// time interval between reinitialization
						// of level-set field
int time_stepping_method;			// select time schmeme 1:explicit euler
						// 2: imex, 3: runge-kutta
					
					
/*++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for interface handling	    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++*/  
				
int number_of_phases;					// number of phases in the domain:
							// 1: single phase flow
							// 2: two phase flow
int apply_mass_distribution_algorithm;  		// =1:apply the mass redistribution algorithm
							// to avoid numerical vapour
int apply_mass_conservation_correction; 		// =1:apply the mass conservation correction
							// to the level-set field
double volume_of_fluid_tolerance;			// tolerance in volume of fluid field
							// for admissable values
double lower_bound_derivatives;			        // lower bound for the first partial derivatives
							// to consider it a limiting case of vanishing
							// partial derivatives
int number_vof_2_level_set_iterations;		        // number of OUTER iterations in the conversion 
							// from volume of fluid to level-set
int number_iterations_ridder;				// maximum number of iterations allowed in the
							// nonlinear root finding algorithm
							// this is the number of INNER iterations
int maximum_reinitialization_steps;			// maximum number of time steps in the reinitialization
							// algorithm
double tolerance_reinitialization;			// stop the reinitialization when the infinite norm of
							// the pseudo time derivative has fallen below this
							// tolerance value
int apply_curvature_smoothing;			        // =1, apply curvature smoothing
							// =0, use unsmoothed curvature
int number_curvature_smoothing_steps;		        // number of iterations applied in the
							// curvature smoothing algorithm
int apply_curvature_smoothing_filter;		        // =1, apply curvature smoothing filter
							// =0, do not apply curvature smoothing filter
int number_redistribution_iterations; 		        // number of iterations allowed to make
							// the volume of fluid field valid
							// these are the sweeps on the vof error
int maximum_number_mass_redistribution_iterations;      // number of iterations allowed to make
							// the volume of fluid field valid
							// these are the sweeps on the vof error
double time_step_mass_redistribution;		        // time step for the mass redistribution
							// algorithm
double redistribution_vof_tolerance; 		        // threshold value of time-derivative 
							// in volume of fluid redistribution equation
double mass_redistribution_diffusion_coefficient;       // diffusion coefficient for mass redistribution equation
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for linear solvers handling	    	  */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
					
double tolerance_pressure;	  		// the tolerance with which the system for the pressure is solved	
int maximum_iterations_allowed_pressure;	// maximum number of iterations allowed for the
						// conjugate gradient method
					
double tolerance_velocity;	  		// the tolerance with which the system for the momentum predictors is solved	
int maximum_iterations_allowed_velocity;	// maximum number of iterations allowed for the
						// conjugate gradient method
					
int number_matrix_connections;		// number of connections in momentum matrix	


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* gravitational acceleration				    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  

vector gravity(0.0,0.0,0.0);			// gravitational acceleration vector 

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* material properties				      	    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  

double rho_plus_over_rho_minus;		// ratio of the densities of the two phases
double smoothing_distance_factor;		// the smoothing distance is smoothing_distance_factor
						// times the smallest mesh width
double rho_minus_over_mu_minus;		// this was the 'Reynolds' number
						// in the original implementation of Sander
double mu_plus_over_mu_minus;			// ratio of the viscosities of both phases
double sigma_over_rho_minus;			// sigma / rho_minus (scaled sigma)


					
					
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for output handling        	    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
					
int vtk_output;				// =1, write output in vtk format
						// =0, skip output in vtk format
int tecplot_output;				// =1, write output in tecplot format
						// =0, skip output in tecplot format
double time_interval_for_output;		// interval in time for which the solution is written to file


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for flow/interface coupling       	    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  

int source_terms_in_momentum_predictor; 	// =1, the source terms are applied in the momentum predictor
						// equation
						// =0, the source terms are applied in the momentum corrector

int continuous_surface_force_model;     	// =1, the continuous surface force model is applied
						// =0, the exact interface boundary conditions are applied

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for initial condition      	    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/  

geometry flow_type;				// the kind of initial condition that has to be applied
bubble the_bubbles[10];			// array with definitions of the bubbles
int number_of_bubbles;			// number of bubbles in the initial condition
surface the_free_surfaces[10];		// array with definitions of the free surfaces
int number_of_free_surfaces; 			// number of bubbles in the domain (<10)
vector initial_velocity(0.0,0.0,0.0);	// the initial velocity at t=0

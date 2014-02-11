#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
enum variable{velocity_u1, velocity_u2, velocity_u3, level_set, pressure_field};
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
/*  Function to do the time-stepping sequence from start to end                 */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
void time_stepping_sequence(
      double ***level_set_new, 				// level set field at new time level
								// mass conserving
      double ***level_set_old, 				// level set field at old time level
								// mass conserving
      double ***volume_of_fluid, 				// volume of fluid field
      double ***curvature,					// interface curvature
      double ***unsmoothed_curvature,			// interface curvature without smoothing
      double ***u_1_velocity_old,				// velocity field at old time level x1 direction
      double ***u_2_velocity_old,				// velocity field at old time level x2 direction
      double ***u_3_velocity_old,				// velocity field at old time level x3 direction
      double ***u_1_velocity_new, 				// velocity field at new time level x1 direction
      double ***u_2_velocity_new, 				// velocity field at new time level x2 direction
      double ***u_3_velocity_new,				// velocity field at new time level x3 direction
      double ***pressure,					// pressure field
      double ***surface_tension_body_force_x1,            // x1 component of the body force due to
                                                          // CSF formulation of surface tension model
      double ***surface_tension_body_force_x2,            // x2 component of the body force due to
                                                          // CSF formulation of surface tension model
      double ***surface_tension_body_force_x3,            // x3 component of the body force due to
                                                          // CSF formulation of surface tension model
      double ***momentum_source_term_u_1,                 // complete source term for the momentum equation
                                                          // in x1 direction=(-p,1+ g_1 +F1), new time level
      double ***momentum_source_term_u_2,                 // complete source term for the momentum equation
                                                          // in x2 direction=(-p,2+ g_2 +F2), new time level
      double ***momentum_source_term_u_3,                 // complete source term for the momentum equation
      double ***scaled_density_u1,                        // scaled density value for the controlvolumes
                                                          // of the momentum equation in x1 direction
      double ***scaled_density_u2,                        // scaled density value for the controlvolumes
                                                          // of the momentum equation in x2 direction
      double ***scaled_density_u3,                        // scaled density value for the controlvolumes
                                                          // of the momentum equation in x3 direction
      boundary_face boundary_faces[6],			// array with all the information
								// for the boundary conditions 
      int number_primary_cells_i,				// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,				// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,				// number of primary (pressure) cells in x3 direction
      double actual_time_step_level_set,			// time step used for level-set advection
								// computed from all stability restrictions and 
								// possibly subscycling
      double mesh_width_x1,					// grid spacing in x1 direction (uniform)
      double mesh_width_x2,					// grid spacing in x2 direction (uniform)
      double mesh_width_x3,					// grid spacing in x3 direction (uniform)
      int apply_mass_distribution_algorithm,    		// =1:apply the mass redistribution algorithm
								// to avoid numerical vapour
      int apply_mass_conservation_correction,    		// =1:apply the mass conservation correction
								// to the level-set field
      double volume_of_fluid_tolerance,			// tolerance in volume of fluid field
								// for admissable values
      double lower_bound_derivatives,			// lower bound for the first partial derivatives
								// to consider it a limiting case of vanishing
								// partial derivatives
      int number_vof_2_level_set_iterations,		// number of OUTER iterations in the conversion 
								// from volume of fluid to level-set
      int number_iterations_ridder,				// maximum number of iterations allowed in the
								// nonlinear root finding algorithm
								// this is the number of INNER iterations
      double vof_2_level_set_tolerance,			// tolerance in the conversion from volume
								// of fluid value to level-set value
      double cfl_number_reinitialization,			// courant-friedrichs-lewy number for the reinitialization
								// equation
      double cfl_number_navier_stokes,			// courant-friedrichs-lewy number for the navier-stokes
								// equations
      int maximum_reinitialization_steps,			// maximum number of time steps in the reinitialization
								// algorithm
      double tolerance_reinitialization,			// stop the reinitialization when the infinite norm of
								// the pseudo time derivative has fallen below this
								// tolerance value
      int number_matrix_connections,				// number of connections in momentum matrix				       
      vector gravity,						// gravitational acceleration vector 
      double actual_time_step_navier_stokes,		// time step used for level-set advection
								// computed from all stability restrictions and 
								// possibly subscycling
      double rho_plus_over_rho_minus,			// ratio of the densities of the two phases
      double smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
								// times the smallest mesh width
      double rho_minus_over_mu_minus,			// this was the 'Reynolds' number
								// in the original implementation of Sander
      double mu_plus_over_mu_minus,				// ratio of the viscosities of both phases
      double tolerance_pressure,	  			// the tolerance with which the system for the pressure is solved	
      int maximum_iterations_allowed_pressure,		// maximum number of iterations allowed for the
								// conjugate gradient method
      double tolerance_velocity,	  			// the tolerance with which the system for the momentum predictor is solved	
      int maximum_iterations_allowed_velocity,		// maximum number of iterations allowed for the
								// conjugate gradient method
      int continuous_surface_force_model,       		// =1, the continuous surface force model is applied
					        		// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor,   		// =1, the source terms are applied in the momentum predictor
					        		// equation
					        		// =0, the source terms are applied in the momentum corrector
      double sigma_over_rho_minus,				// sigma / rho_minus (scaled sigma)
      double time_step_restriction_global,      		// upper bound on the time step, set by user
      int fixed_time_step,					// =1 the time step is not restricted by 
								// physical constraints 
								// =0 the time step is adapted to the physical
								// constraints
      int apply_curvature_smoothing,				// =1, apply curvature smoothing
								// =0, use unsmoothed curvature
      int number_curvature_smoothing_steps,			// number of iterations applied in the
								// curvature smoothing algorithm
      int apply_curvature_smoothing_filter,			// =1, apply curvature smoothing filter
								// =0, do not apply curvature smoothing filter
      int number_of_subcycles,				// for every time step for the Navier-Stokes
								// equations number_of_subcycles steps are taken
								// with the interface evolution equations
      int vtk_output,						// =1, write output in vtk format
								// =0, skip output in vtk format
      int tecplot_output,					// =1, write output in tecplot format
								// =0, skip output in tecplot format
      double time_interval_for_output,			// interval in time for which the solution is written to file
      double start_time_simulation,				// starting time for the simulation
      double end_time_simulation,				// end time for the simulation
      restart_parameters my_restart_parameters,  		// all parameters for reading/writing restart files
      int maximum_number_mass_redistribution_iterations, 	// number of iterations allowed to make
								// the volume of fluid field valid
								// these are the sweeps on the vof error
      double time_step_mass_redistribution,			// time step for the mass redistribution
								// algorithm
      double redistribution_vof_tolerance, 			// threshold value of time-derivative 
								// in volume of fluid redistribution equation
      int number_of_phases,					// number of phases in the domain:
								// 1: single phase flow
								// 2: two phase flow
      int debugging_mode,					// =1, program is run in 
      								// debugging mode with additional output and checks
      								// =0, progam is run in production mode
      double time_interval_for_reinitialization 		// time interval between reinitialization
								// of level-set field

						
  
	 )

{
	/* function definitions */
       void compute_time_step_size(				// compute time step size
          double ***u_1_velocity_new, 		
          double ***u_2_velocity_new, 		
          double ***u_3_velocity_new,		
          double rho_plus_over_rho_minus,		
          double maximum_weighted_curvature,	
          double sigma_over_rho_minus,		
          double cfl_number_navier_stokes,		
          int number_primary_cells_i,		
          int number_primary_cells_j,		
          int number_primary_cells_k,		
          double mesh_width_x1,			
          double mesh_width_x2,			
          double mesh_width_x3,			
          double time_step_restriction_global,	
          int fixed_time_step,			
          double actual_time_step_navier_stokes	
	   );
       void advance_interface(				// advance the interface to the new time level
          double ***level_set_new, 			
          double ***level_set_old, 			
          double ***volume_of_fluid, 		
          double ***u_1_velocity_new, 		
          double ***u_2_velocity_new, 		
          double ***u_3_velocity_new,		
          int number_primary_cells_i,		
          int number_primary_cells_j,		
          int number_primary_cells_k,		
          double actual_time_step_level_set,	
          double mesh_width_x1,			
          double mesh_width_x2,			
          double mesh_width_x3,			
          int apply_mass_distribution_algorithm,    
          int apply_mass_conservation_correction,   
          double volume_of_fluid_tolerance,		
          double lower_bound_derivatives,		
          int number_vof_2_level_set_iterations,	
          int number_iterations_ridder,		
          double vof_2_level_set_tolerance,		
          double cfl_number_reinitialization,	
          int maximum_reinitialization_steps,	
          double tolerance_reinitialization,		
      	   int maximum_number_mass_redistribution_iterations, 
      	   double time_step_mass_redistribution,		
	   double redistribution_vof_tolerance,
	   double &time_to_reinitialize,
	   double time_interval_reinitialization,
	   double actual_time
      );
      void advance_coupling_part1(		             // update the coupling between interface and flow
 	   double ***level_set, 			      // density and surface tension
	   double ***pressure,				
	   double ***curvature, 				
	   double ***smoothed_curvature,			
	   double ***surface_tension_body_force_x1,	
	   double ***surface_tension_body_force_x2,	
	   double ***surface_tension_body_force_x3,	
	   double ***momentum_source_term_u_1, 	
	   double ***momentum_source_term_u_2, 	
	   double ***momentum_source_term_u_3, 
	   double ***scaled_density_u1,
	   double ***scaled_density_u2,
	   double ***scaled_density_u3,
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
	   double maximum_weighted_curvature,		
	   int apply_curvature_smoothing,		
	   int number_curvature_smoothing_steps,	
	   int apply_curvature_smoothing_filter,
	   double smoothing_distance_factor,
	   double lower_bound_derivatives
      );
      void advance_coupling_part2(                     // update the coupling between interface and flow
          double ***level_set,                         // full momentum source term
          double ***pressure,
          double ***curvature,
          double ***smoothed_curvature,
          double ***surface_tension_body_force_x1,
          double ***surface_tension_body_force_x2,
          double ***surface_tension_body_force_x3,
          double ***momentum_source_term_u_1,
          double ***momentum_source_term_u_2,
          double ***momentum_source_term_u_3,
          double ***scaled_density_u1,
          double ***scaled_density_u2,
          double ***scaled_density_u3,
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
          double maximum_weighted_curvature,
          int apply_curvature_smoothing,
          int number_curvature_smoothing_steps,
          int apply_curvature_smoothing_filter,
          double smoothing_distance_factor,
          double lower_bound_derivatives
      );
      void advance_flow_field(				// advance the flow field to the new time level
          double ***level_set,			
          double ***pressure,			
          double ***u_1_velocity_old,		
          double ***u_2_velocity_old,		
          double ***u_3_velocity_old,		
          double ***u_1_velocity_new,		
          double ***u_2_velocity_new,		
          double ***u_3_velocity_new,		
          double ***momentum_source_term_u_1,	
          double ***momentum_source_term_u_2,	
          double ***momentum_source_term_u_3,	
          double ***surface_tension_body_force_x1,	
          double ***surface_tension_body_force_x2,	
          double ***surface_tension_body_force_x3,	
	   double ***scaled_density_u1,
	   double ***scaled_density_u2,
	   double ***scaled_density_u3,
          boundary_face boundary_faces[6],		
          double mesh_width_x1,			
          double mesh_width_x2,			
          double mesh_width_x3,			
          int number_primary_cells_i,		
          int number_primary_cells_j,		
          int number_primary_cells_k,		
          int number_matrix_connections,		 
          vector gravity,				
          double actual_time_step_navier_stokes,	
          double rho_plus_over_rho_minus,		
          double smoothing_distance_factor,		
          double rho_minus_over_mu_minus,		
          double mu_plus_over_mu_minus,		
          double tolerance_pressure,	  	
          int maximum_iterations_allowed_pressure,	
          double tolerance_velocity,	  	
          int maximum_iterations_allowed_velocity,	
          int continuous_surface_force_model,       
          int source_terms_in_momentum_predictor    
	   );
    void output_solution(					// output the solution to file
	  double ***level_set_new, 		
	  double ***volume_of_fluid, 	
	  double ***curvature,
	  double ***smoothed_curvature,
	  double ***u_1_velocity_new, 		
	  double ***u_2_velocity_new, 		
	  double ***u_3_velocity_new,		
	  double ***pressure,			
	  int vtk_output,			
	  int tecplot_output,			
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k,		
	  double mesh_width_x1,			
	  double mesh_width_x2,			
	  double mesh_width_x3,			
	  int index_of_output_file
	  );
      double ***double_Matrix2(				// allocate memory for a three-dimensional array of doubles
	  int number_primary_cells_i,		
	  int number_primary_cells_j, 		
	  int number_primary_cells_k
	  );
      void   free_double_Matrix2( 				// deallocate memory for a three-dimensional array of doubles
	  double ***doubleMatrix2, 
	  int number_primary_cells_i,	
	  int number_primary_cells_j
	  );
      void match_level_set_to_volume_of_fluid(		// compute a new level-set field
	  double ***level_set_not_mass_conserving, 		// that corresponds to given
	  double ***volume_of_fluid,				// volume of fluid field
	  double ***level_set_mass_conserving,		// using an existing level-set
	  int number_primary_cells_i, 			// field as a starting value
	  int number_primary_cells_j, 
	  int number_primary_cells_k,			
	  double volume_of_fluid_tolerance,		
	  double lower_bound_derivatives,		
	  int number_vof_2_level_set_iterations,		
	  int number_iterations_ridder,			
	  double vof_2_level_set_tolerance		
	  );
      void write_restart_file(				// write solution to restart file
	  double time_of_restart,
	  double ***level_set_new, 			
	  double ***volume_of_fluid, 			
	  double ***u_1_velocity_new, 			
	  double ***u_2_velocity_new, 			
	  double ***u_3_velocity_new,			
	  double ***pressure,				
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k,			
	  restart_parameters my_restart_parameters	
	  );
      void analyse_interface_properties(
      	  double ***volume_of_fluid, 	
	  double ***level_set_new,
      	  double ***u_1_velocity_new,
      	  double ***u_2_velocity_new,
      	  double ***u_3_velocity_new,
      	  int number_primary_cells_i,
      	  int number_primary_cells_j,
      	  int number_primary_cells_k,
      	  double mesh_width_x1,	
      	  double mesh_width_x2,		
      	  double mesh_width_x3,		
	  double actual_time,
	  double lower_bound_derivatives
         );

    int subcycle_index;					// index of the subcycle for the interface evolution  
    int index_of_output_file=1;				// index of the output file
    double actual_time;					// actual time in the simulation for which the solution
    double time_to_write_output;				// time at which new output has to be written
    double time_to_reinitialize;				// time at which reinitialization of the 
    								// level-set field will be performed
    double maximum_weighted_curvature;			// maximum 'active' value of the curvature 
								// used to evaluate the capillary time step
								// restriction 
    

      /* define the start settings for time integration */
      
      actual_time=start_time_simulation;
      time_to_write_output=start_time_simulation+time_interval_for_output;
      time_to_reinitialize=start_time_simulation+time_interval_for_reinitialization;
    

    
    while(actual_time<end_time_simulation)
    {

       /* compute the time step size for advancing the flow field */
       /* and advancing the interface */
	  
	  compute_time_step_size(u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,		
				  rho_plus_over_rho_minus, maximum_weighted_curvature,	
				    sigma_over_rho_minus, cfl_number_navier_stokes,		
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
					mesh_width_x1, mesh_width_x2, mesh_width_x3,			
					  time_step_restriction_global,	fixed_time_step,			
					    actual_time_step_navier_stokes);

	  
	  actual_time+=actual_time_step_navier_stokes;
	  actual_time_step_level_set=actual_time_step_navier_stokes/number_of_subcycles;
	  std::cerr<<"----------------------------------------------------------\n";
	  std::cerr<<"actual time : "<< actual_time <<" \n";
	  std::cerr<<"actual time step size Navier-Stokes: "<<actual_time_step_navier_stokes<<"\n";
	  std::cerr<<"actual time step size interface: "<<actual_time_step_level_set<< " \n";
	  std::cerr<<"end_time_simulation: "<<end_time_simulation<< " \n";
	  std::cerr<<"----------------------------------------------------------\n";
	  std::cerr<<"start computing solution for t= "<< actual_time <<" \n";
	 
	  if(number_of_phases>1)
	  {
       /* advance interface over one time step       */
       /* but this may involve a number of subcycles */
       
	  	for(subcycle_index=0;subcycle_index<number_of_subcycles;subcycle_index++)
	  	{
			advance_interface( level_set_new, level_set_old, volume_of_fluid,
			  u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,
			    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
			      actual_time_step_level_set,
				mesh_width_x1, mesh_width_x2, mesh_width_x3,
				  apply_mass_distribution_algorithm,
				    apply_mass_conservation_correction,
				      volume_of_fluid_tolerance,
					lower_bound_derivatives,
					  number_vof_2_level_set_iterations,
					    number_iterations_ridder,
					      vof_2_level_set_tolerance,
						cfl_number_reinitialization,
						  maximum_reinitialization_steps,
						    tolerance_reinitialization,
						      maximum_number_mass_redistribution_iterations,
      	   					 time_step_mass_redistribution, redistribution_vof_tolerance,
						time_to_reinitialize, time_interval_for_reinitialization,
			 			actual_time
   					);
	  	}
	  
	  	std::cerr<<"interface updated \n";
		
		/* analyse the properties of the volume enclosed by the interface */
		/* and write those to a file */
		
		analyse_interface_properties( volume_of_fluid, level_set_new,	
      						 u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,	
      						   number_primary_cells_i, number_primary_cells_j,	number_primary_cells_k,	
      							mesh_width_x1,	mesh_width_x2,	mesh_width_x3,	actual_time,
				    			   lower_bound_derivatives
     					);

	  }

       /* compute interface coupling */
       /* the coupling between interface and the flow field  is through the body force  */
       /* and the density                                                               */
	
         advance_coupling_part1( level_set_new, pressure, curvature, unsmoothed_curvature,
                     surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                        scaled_density_u1, scaled_density_u2, scaled_density_u3,
                         number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                           mesh_width_x1, mesh_width_x2, mesh_width_x3,
                            rho_plus_over_rho_minus, actual_time_step_navier_stokes,
                              gravity, sigma_over_rho_minus,
                                maximum_weighted_curvature,
                                  apply_curvature_smoothing, number_curvature_smoothing_steps,
                                   apply_curvature_smoothing_filter,
                                     smoothing_distance_factor, lower_bound_derivatives);
		
	  	std::cerr<<"coupling updated part 1: surface tension and density \n";
      
       /* advance flow field over one time step */
       
         advance_flow_field( level_set_new, pressure,
                       u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,
                         u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,
                           momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                            surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                             scaled_density_u1, scaled_density_u2, scaled_density_u3,
                              boundary_faces,
                                mesh_width_x1, mesh_width_x2, mesh_width_x3,
                                  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                                   number_matrix_connections, gravity,
                                     actual_time_step_navier_stokes,
                                       rho_plus_over_rho_minus,
                                         smoothing_distance_factor,
                                          rho_minus_over_mu_minus, mu_plus_over_mu_minus,
                                           tolerance_pressure, maximum_iterations_allowed_pressure,
                                             tolerance_velocity, maximum_iterations_allowed_velocity,
                                                continuous_surface_force_model,
                                                 source_terms_in_momentum_predictor);
	  std::cerr<<"flow field updated \n";

         advance_coupling_part2( level_set_new, pressure, curvature, unsmoothed_curvature,
                     surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                        scaled_density_u1, scaled_density_u2, scaled_density_u3,
                         number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                           mesh_width_x1, mesh_width_x2, mesh_width_x3,
                            rho_plus_over_rho_minus, actual_time_step_navier_stokes,
                              gravity, sigma_over_rho_minus,
                                maximum_weighted_curvature,
                                  apply_curvature_smoothing, number_curvature_smoothing_steps,
                                   apply_curvature_smoothing_filter,
                                     smoothing_distance_factor, lower_bound_derivatives);

         std::cerr<<"coupling updated part 2: complete momentum source term \n";
	  
      /* output the solution to file for postprocessing */
       
	if(actual_time>=time_to_write_output)
	{
	   /* if the amount of time between two subsequent outputs has passed 	*/
	   /* output the solution and set the time for when the next  output 	*/
	   /* should occur							*/


	   /* build the  output file name for this instant of outputting */
	   /* the index of the solution file is included in the name    */
	   /* of the output file					*/
	 
	      output_solution( level_set_new, volume_of_fluid, curvature, unsmoothed_curvature,		
			  u_1_velocity_new, u_2_velocity_new, u_3_velocity_new, pressure,			
			    vtk_output, tecplot_output,			
			      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
				mesh_width_x1, mesh_width_x2, mesh_width_x3, index_of_output_file);
	      index_of_output_file++;
	      time_to_write_output+=time_interval_for_output;
	 }
	  std::cerr<<"completed solution for t= "<< actual_time <<" \n";
	  std::cerr<<"----------------------------------------------------------\n \n";
     
    }
    
    /* if a restart file should be written, it is written to binary file now */
    
    if(my_restart_parameters.write_solution_to_restart_file)
    {
	      write_restart_file(actual_time, level_set_new, volume_of_fluid, 			
				    u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,			
				      pressure,				
					 number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
					  my_restart_parameters);
      
    }
    
    /* deallocate memory for the coupling terms between flow field and interface */
    
    free_double_Matrix2(scaled_density_u1, number_primary_cells_i+1, 
		   				    number_primary_cells_j+2);
    free_double_Matrix2(scaled_density_u2, number_primary_cells_i+2, 
						    number_primary_cells_j+1);
    free_double_Matrix2(scaled_density_u3, number_primary_cells_i+2, 
						    number_primary_cells_j+2);
    

    free_double_Matrix2(surface_tension_body_force_x1, number_primary_cells_i+1, 
		    number_primary_cells_j+2);
    free_double_Matrix2(surface_tension_body_force_x2, number_primary_cells_i+2, 
		    number_primary_cells_j+1);
    free_double_Matrix2(surface_tension_body_force_x3, number_primary_cells_i+2, 
		    number_primary_cells_j+2);
    
    free_double_Matrix2(momentum_source_term_u_1, number_primary_cells_i+1,
                         number_primary_cells_j+2);
    free_double_Matrix2(momentum_source_term_u_2, number_primary_cells_i+2,
                         number_primary_cells_j+1);
    free_double_Matrix2(momentum_source_term_u_3, number_primary_cells_i+2,
                         number_primary_cells_j+2);
   
    /* deallocate storage for the interface curvature, both the directly computed and smoothed one */
	  
     free_double_Matrix2(curvature, number_primary_cells_i+2,number_primary_cells_j+2); 
     free_double_Matrix2(unsmoothed_curvature, number_primary_cells_i+2,number_primary_cells_j+2);
  
}
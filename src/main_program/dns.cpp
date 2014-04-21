#include "../headers/array.h"
/********************************************************************************/
/********************************************************************************/
/*  Main program for the Cartesian Mass Conserving Level Set Method         	*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 03-01       								*/
/*  Update	:        								*/
/********************************************************************************/
//
//
/********************************************************************************/
#include "../headers/header_main.h" 
/********************************************************************************/

int main ()
{
 

/* initialization */
   
      initialize_computation( u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,		
				  u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,		
				   pressure, level_set_old, level_set_new, volume_of_fluid,		
				    curvature, unsmoothed_curvature,
                                   momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                                    surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                                    scaled_density_u1, scaled_density_u2, scaled_density_u3, boundary_faces,
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
					actual_time_step_level_set, mesh_width_x1, mesh_width_x2, mesh_width_x3,			
					  apply_mass_distribution_algorithm, apply_mass_conservation_correction,   
					    volume_of_fluid_tolerance, lower_bound_derivatives, number_vof_2_level_set_iterations,	
					      number_iterations_ridder, vof_2_level_set_tolerance, cfl_number_reinitialization,	
						cfl_number_navier_stokes, maximum_reinitialization_steps, tolerance_reinitialization,	
						  number_matrix_connections, gravity, actual_time_step_navier_stokes,	
						    rho_plus_over_rho_minus, smoothing_distance_factor, rho_minus_over_mu_minus,		
						      mu_plus_over_mu_minus, tolerance_pressure, maximum_iterations_allowed_pressure,
							tolerance_velocity, maximum_iterations_allowed_velocity,
							  continuous_surface_force_model, source_terms_in_momentum_predictor,   
							    sigma_over_rho_minus, time_step_restriction_global,      
							      fixed_time_step, apply_curvature_smoothing, number_curvature_smoothing_steps,	
								apply_curvature_smoothing_filter, number_of_subcycles,			
								 vtk_output, tecplot_output, time_interval_for_output,		
								   domain_size_x1, domain_size_x2, domain_size_x3, flow_type,			
								  the_bubbles, number_of_bubbles, the_free_surfaces,		
								 number_of_free_surfaces, 
								start_time_simulation, end_time_simulation, time_stepping_method, initial_velocity,
							      my_restart_parameters, maximum_number_mass_redistribution_iterations, 
  							    time_step_mass_redistribution, redistribution_vof_tolerance,
							  number_of_phases, debugging_mode, time_interval_for_reinitialization
				);
   
    std::cerr << "Initialization completed \n";

/* time stepping sequence */


      time_stepping_sequence( level_set_new, level_set_old, volume_of_fluid, 	
			     curvature, unsmoothed_curvature,
			      u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,		
				u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,		
				  pressure, momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                                surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                                  scaled_density_u1, scaled_density_u2, scaled_density_u3, boundary_faces,
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
					actual_time_step_level_set, mesh_width_x1, mesh_width_x2, mesh_width_x3,			
					  apply_mass_distribution_algorithm, apply_mass_conservation_correction,   
					    volume_of_fluid_tolerance, lower_bound_derivatives,		
					      number_vof_2_level_set_iterations, number_iterations_ridder,		
						vof_2_level_set_tolerance, cfl_number_reinitialization,	cfl_number_navier_stokes,		
						  maximum_reinitialization_steps, tolerance_reinitialization, number_matrix_connections,				       
						    gravity, actual_time_step_navier_stokes, rho_plus_over_rho_minus,		
						      smoothing_distance_factor, rho_minus_over_mu_minus, mu_plus_over_mu_minus,		
							tolerance_pressure, maximum_iterations_allowed_pressure, 
							  tolerance_velocity, maximum_iterations_allowed_velocity, 
							    continuous_surface_force_model,       
							      source_terms_in_momentum_predictor, sigma_over_rho_minus,		
								time_step_restriction_global, fixed_time_step, apply_curvature_smoothing,		
								  number_curvature_smoothing_steps,	apply_curvature_smoothing_filter,	
								    number_of_subcycles, vtk_output, tecplot_output, time_interval_for_output,
								      start_time_simulation, end_time_simulation, time_stepping_method,
								    my_restart_parameters, maximum_number_mass_redistribution_iterations,
								time_step_mass_redistribution, redistribution_vof_tolerance, number_of_phases,
								debugging_mode, time_interval_for_reinitialization, mass_redistribution_diffusion_coefficient
				);

    std::cout << "Time stepping completed \n";

    return 0;
}













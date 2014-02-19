CXX=g++
RM=rm -f
CFLAGS= -c -O3 -include ../objects/funcdefs.h -Wall
#############################################
OBJS_DIR:=../objects/
EXEC_DIR:=../executable/
HEADER_DIR:=./headers/
MAIN_DIR:=./main_program/
INITIAL_DIR:=./initialization/
UTILS_DIR:=./utils/
MOMENTUM_DIR:=./momentum_equation/
INTERFACE_DIR:=./interface/
LINEAR_DIR:=./linear_solver/
PRESS_DIR:=./pressure_correction_equation/
COUPL_DIR:=./interface_coupling/
OUTPUT_DIR:=./post_processing/
RESTART_DIR:=./restart/
#############################################

OBJS=$(addprefix $(OBJS_DIR),dns.o time_stepping_sequence.o advance_flow_field.o \
				compute_time_step_size.o \
				initialize_computation.o set_parameters.o \
				set_boundary_conditions.o allocate_main_variables.o \
				initialize_all_variables.o \
				initialize_flow_field.o utilities.o apply_boundary_conditions_velocity.o  \
				apply_boundary_conditions_velocity_u1.o apply_boundary_conditions_velocity_u2.o \
				apply_boundary_conditions_velocity_u3.o initialize_interface.o \
				initialize_level_set.o initialize_volume_of_fluid.o \
				initialize_bubbles.o initialize_free_surface.o compute_volume_of_fluid.o \
				compute_level_set_gradient.o level_set_2_vof.o \
				advance_interface.o advect_level_set.o \
				apply_ridders_method.o apply_volume_of_fluid_clipping.o \
				check_volume_of_fluid.o compute_derivatives_level_set.o compute_redistribution_velocity.o \
				compute_vof_flux_x1.o compute_vof_flux_x2.o compute_vof_flux_x3.o copy_cell_centered_field.o \
				field_neumann_boundary.o level_set_2_vof_phi_negative.o \
				make_level_set_mass_conserving.o match_level_set_to_volume_of_fluid.o \
				modify_volume_of_fluid_values.o reinitialize_level_set.o \
				update_volume_of_fluid_flux_x1.o update_volume_of_fluid_flux_x2.o update_volume_of_fluid_flux_x3.o \
				upwind_flux_mass_redistribution.o vof_2_level_set.o \
				compute_level_set_flux.o build_momentum_matrix_u1.o \
				build_momentum_matrix_u2.o build_momentum_matrix_u3.o \
				build_momentum_predictor_u1.o build_momentum_predictor_u2.o \
				build_momentum_predictor_u3.o build_momentum_rhs_u1.o \
				build_momentum_rhs_u2.o build_momentum_rhs_u3.o \
				compute_convection_term.o  \
				compute_scaled_density.o compute_scaled_viscosity.o \
				compute_weighted_average.o  \
				fold_momentum_matrix_u1.o fold_momentum_matrix_u2.o \
				fold_momentum_matrix_u3.o fold_momentum_rhside_u1.o \
				fold_momentum_rhside_u2.o fold_momentum_rhside_u3.o \
				shift_velocity_field.o solve_momentum_predictor.o \
				solve_momentum_predictor_u1.o solve_momentum_predictor_u2.o \
				solve_momentum_predictor_u3.o apply_preconditioner.o \
				build_preconditioner.o conjugate_gradient_method.o \
				apply_pressure_correction.o build_pressure_matrix.o \
				build_pressure_rhs_boundary.o build_pressure_rhs_momentum_part.o \
				build_pressure_system.o set_pressure_boundary_condition.o \
				solve_momentum_corrector.o solve_pressure_correction_system.o \
				advance_coupling.o compute_body_force_x1.o compute_body_force_x2.o \
				compute_body_force_x3.o compute_curvature.o compute_curvature_error_laplace.o \
				computed_derivative_heaviside_function.o compute_momentum_source_terms.o \
				compute_surface_tension_body_force.o curvature_filter.o \
				delta_function.o smooth_curvature.o \
				interpolate_velocity_u1_center.o interpolate_velocity_u2_center.o \
				interpolate_velocity_u3_center.o output_solution.o \
				interpolate_velocity_u1_vertex.o interpolate_velocity_u2_vertex.o \
				interpolate_velocity_u3_vertex.o \
				write_cell_centered_field_tecplot.o write_cell_centered_field_vtk.o \
				write_vertex_centered_vector_field_vtk.o write_coordinates_tecplot.o \
				write_coordinates_vtk.o  write_vertex_centered_field_vtk.o \
				apply_bisection.o map_index_u1.o map_index_u2.o map_index_u3.o \
				map_index_pressure.o  \
				decompress_solution_pressure.o copy_general_field.o \
				export_matrix_matlab.o compress_solution_velocity_u1.o \
				compress_solution_velocity_u2.o compress_solution_velocity_u3.o \
				decompress_solution_velocity_u1.o decompress_solution_velocity_u2.o \
				decompress_solution_velocity_u3.o field_extrapolate_boundary.o \
				project_pressure_rhside.o shift_interface.o \
				write_restart_file.o read_restart_file.o \
				compute_normal_derivative_at_faces.o \
				redistribute_volume_of_fluid_error.o \
				apply_volume_of_fluid_redistribution.o initialize_curvature.o \
				dump_solution_for_debugging.o dump_curvature_for_debugging.o \
				check_symmetry_scalars.o compute_weighted_curvature.o check_symmetry_velocities.o \
				write_interface_solution.o dump_adapted_vof_for_debugging.o \
				compute_mass_in_domain.o analyse_interface_properties.o \
				compute_vof_at_u1_points.o compute_vof_at_u2_points.o \
				compute_vof_at_u3_points.o compute_scaled_density_u1.o \
				compute_scaled_density_u2.o compute_scaled_density_u3.o \
				compute_density_u_controlvolumes.o shift_pressure_solution.o \
		  		initialize_coupling.o initialize_pressure.o \
		  		build_pressure_rhs_initialization.o build_pressure_system_initialization.o \
		  		advance_coupling_part1.o advance_coupling_part2.o \
		  		dump_divergence_for_debugging.o apply_boundary_conditions_pressure.o \
				output_predictor_velocityfield.o compress_solution_pressure.o \
				dump_reinitialization_for_debugging.o compute_first_stage_RK.o \
				compute_second_stage_RK.o compute_new_time_RK.o \
				advect_level_set_higher_order.o weno_flux_computation.o \
				evaluate_convection_operator_weno.o)
#                               apply_redistribution_volume_of_fluid_correction.o 
all:$(EXEC_DIR)MCLS

$(EXEC_DIR)MCLS:../objects/funcdefs.h $(OBJS)
	$(CXX) $(OBJS) -O3 -o $(EXEC_DIR)MCLS

../objects/funcdefs.h:
	./gen_funcdefs.h

#---------------------------#
#  main program functions   #
#---------------------------#

$(OBJS_DIR)dns.o:$(MAIN_DIR)dns.cpp $(HEADER_DIR)header_main.h
	$(CXX) $(CFLAGS) $(MAIN_DIR)dns.cpp -o $(OBJS_DIR)dns.o
$(OBJS_DIR)time_stepping_sequence.o:$(MAIN_DIR)time_stepping_sequence.cpp 
	$(CXX) $(CFLAGS) $(MAIN_DIR)time_stepping_sequence.cpp -o $(OBJS_DIR)time_stepping_sequence.o
$(OBJS_DIR)advance_flow_field.o:$(MAIN_DIR)advance_flow_field.cpp 
	$(CXX) $(CFLAGS) $(MAIN_DIR)advance_flow_field.cpp -o $(OBJS_DIR)advance_flow_field.o
$(OBJS_DIR)compute_time_step_size.o:$(MAIN_DIR)compute_time_step_size.cpp 
	$(CXX) $(CFLAGS) $(MAIN_DIR)compute_time_step_size.cpp -o $(OBJS_DIR)compute_time_step_size.o

#-----------------------------#
#  initialization functions   #
#-----------------------------#

$(OBJS_DIR)initialize_computation.o:$(INITIAL_DIR)initialize_computation.cpp \
				$(INITIAL_DIR)set_parameters.cpp \
				$(INITIAL_DIR)set_boundary_conditions.cpp \
				$(INITIAL_DIR)allocate_main_variables.cpp \
				$(INITIAL_DIR)initialize_flow_field.cpp \
				$(INITIAL_DIR)initialize_all_variables.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_computation.cpp -o $(OBJS_DIR)initialize_computation.o
$(OBJS_DIR)set_parameters.o:$(INITIAL_DIR)set_parameters.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)set_parameters.cpp -o $(OBJS_DIR)set_parameters.o
$(OBJS_DIR)set_boundary_conditions.o:$(INITIAL_DIR)set_boundary_conditions.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)set_boundary_conditions.cpp -o $(OBJS_DIR)set_boundary_conditions.o
$(OBJS_DIR)allocate_main_variables.o:$(INITIAL_DIR)allocate_main_variables.cpp
	$(CXX) $(CFLAGS) $(INITIAL_DIR)allocate_main_variables.cpp -o $(OBJS_DIR)allocate_main_variables.o
$(OBJS_DIR)initialize_all_variables.o:$(INITIAL_DIR)initialize_all_variables.cpp \
 				$(INITIAL_DIR)initialize_flow_field.cpp \
 				$(INITIAL_DIR)initialize_interface.cpp
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_all_variables.cpp -o $(OBJS_DIR)initialize_all_variables.o
$(OBJS_DIR)initialize_flow_field.o:$(INITIAL_DIR)initialize_flow_field.cpp $(UTILS_DIR)utilities.cpp	
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_flow_field.cpp -o $(OBJS_DIR)initialize_flow_field.o
$(OBJS_DIR)initialize_interface.o:$(INITIAL_DIR)initialize_interface.cpp $(UTILS_DIR)utilities.cpp\
 				$(INITIAL_DIR)initialize_level_set.cpp \
 				$(INITIAL_DIR)initialize_volume_of_fluid.cpp
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_interface.cpp -o $(OBJS_DIR)initialize_interface.o
$(OBJS_DIR)initialize_level_set.o:$(INITIAL_DIR)initialize_level_set.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_level_set.cpp -o $(OBJS_DIR)initialize_level_set.o
$(OBJS_DIR)initialize_volume_of_fluid.o:$(INITIAL_DIR)initialize_volume_of_fluid.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_volume_of_fluid.cpp -o $(OBJS_DIR)initialize_volume_of_fluid.o
$(OBJS_DIR)initialize_bubbles.o:$(INITIAL_DIR)initialize_bubbles.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_bubbles.cpp -o $(OBJS_DIR)initialize_bubbles.o
$(OBJS_DIR)initialize_free_surface.o:$(INITIAL_DIR)initialize_free_surface.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_free_surface.cpp -o $(OBJS_DIR)initialize_free_surface.o
$(OBJS_DIR)initialize_curvature.o:$(INITIAL_DIR)initialize_curvature.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_curvature.cpp -o $(OBJS_DIR)initialize_curvature.o
$(OBJS_DIR)write_interface_solution.o:$(INITIAL_DIR)write_interface_solution.cpp 
	$(CXX) $(CFLAGS) $(INITIAL_DIR)write_interface_solution.cpp -o $(OBJS_DIR)write_interface_solution.o
$(OBJS_DIR)initialize_coupling.o:$(INITIAL_DIR)initialize_coupling.cpp
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_coupling.cpp -o $(OBJS_DIR)initialize_coupling.o
$(OBJS_DIR)initialize_pressure.o:$(INITIAL_DIR)initialize_pressure.cpp
	$(CXX) $(CFLAGS) $(INITIAL_DIR)initialize_pressure.cpp -o $(OBJS_DIR)initialize_pressure.o
$(OBJS_DIR)build_pressure_rhs_initialization.o:$(INITIAL_DIR)build_pressure_rhs_initialization.cpp
	$(CXX) $(CFLAGS) $(INITIAL_DIR)build_pressure_rhs_initialization.cpp -o $(OBJS_DIR)build_pressure_rhs_initialization.o
$(OBJS_DIR)build_pressure_system_initialization.o:$(INITIAL_DIR)build_pressure_system_initialization.cpp
	$(CXX) $(CFLAGS) $(INITIAL_DIR)build_pressure_system_initialization.cpp -o $(OBJS_DIR)build_pressure_system_initialization.o

#-----------------------#
#  coupling functions   #
#-----------------------#

$(OBJS_DIR)advance_coupling.o:$(COUPL_DIR)advance_coupling.cpp
	$(CXX) $(CFLAGS) $(COUPL_DIR)advance_coupling.cpp -o $(OBJS_DIR)advance_coupling.o
$(OBJS_DIR)advance_coupling_part1.o:$(COUPL_DIR)advance_coupling_part1.cpp
	$(CXX) $(CFLAGS) $(COUPL_DIR)advance_coupling_part1.cpp -o $(OBJS_DIR)advance_coupling_part1.o
$(OBJS_DIR)advance_coupling_part2.o:$(COUPL_DIR)advance_coupling_part2.cpp
	$(CXX) $(CFLAGS) $(COUPL_DIR)advance_coupling_part2.cpp -o $(OBJS_DIR)advance_coupling_part2.o
$(OBJS_DIR)compute_body_force_x1.o:$(COUPL_DIR)compute_body_force_x1.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_body_force_x1.cpp -o $(OBJS_DIR)compute_body_force_x1.o
$(OBJS_DIR)compute_body_force_x2.o:$(COUPL_DIR)compute_body_force_x2.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_body_force_x2.cpp -o $(OBJS_DIR)compute_body_force_x2.o
$(OBJS_DIR)compute_body_force_x3.o:$(COUPL_DIR)compute_body_force_x3.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_body_force_x3.cpp -o $(OBJS_DIR)compute_body_force_x3.o
$(OBJS_DIR)compute_curvature.o:$(COUPL_DIR)compute_curvature.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_curvature.cpp -o $(OBJS_DIR)compute_curvature.o
$(OBJS_DIR)compute_curvature_error_laplace.o:$(COUPL_DIR)compute_curvature_error_laplace.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_curvature_error_laplace.cpp -o $(OBJS_DIR)compute_curvature_error_laplace.o
$(OBJS_DIR)computed_derivative_heaviside_function.o:$(COUPL_DIR)computed_derivative_heaviside_function.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)computed_derivative_heaviside_function.cpp -o $(OBJS_DIR)computed_derivative_heaviside_function.o
$(OBJS_DIR)compute_momentum_source_terms.o:$(COUPL_DIR)compute_momentum_source_terms.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_momentum_source_terms.cpp -o $(OBJS_DIR)compute_momentum_source_terms.o
$(OBJS_DIR)compute_surface_tension_body_force.o:$(COUPL_DIR)compute_surface_tension_body_force.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_surface_tension_body_force.cpp -o $(OBJS_DIR)compute_surface_tension_body_force.o
$(OBJS_DIR)curvature_filter.o:$(COUPL_DIR)curvature_filter.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)curvature_filter.cpp -o $(OBJS_DIR)curvature_filter.o
$(OBJS_DIR)delta_function.o:$(COUPL_DIR)delta_function.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)delta_function.cpp -o $(OBJS_DIR)delta_function.o
$(OBJS_DIR)smooth_curvature.o:$(COUPL_DIR)smooth_curvature.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)smooth_curvature.cpp -o $(OBJS_DIR)smooth_curvature.o
$(OBJS_DIR)dump_curvature_for_debugging.o:$(COUPL_DIR)dump_curvature_for_debugging.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)dump_curvature_for_debugging.cpp -o $(OBJS_DIR)dump_curvature_for_debugging.o
$(OBJS_DIR)compute_weighted_curvature.o:$(COUPL_DIR)compute_weighted_curvature.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_weighted_curvature.cpp -o $(OBJS_DIR)compute_weighted_curvature.o
$(OBJS_DIR)compute_scaled_density_u1.o:$(COUPL_DIR)compute_scaled_density_u1.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_scaled_density_u1.cpp -o $(OBJS_DIR)compute_scaled_density_u1.o
$(OBJS_DIR)compute_scaled_density_u2.o:$(COUPL_DIR)compute_scaled_density_u2.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_scaled_density_u2.cpp -o $(OBJS_DIR)compute_scaled_density_u2.o
$(OBJS_DIR)compute_scaled_density_u3.o:$(COUPL_DIR)compute_scaled_density_u3.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_scaled_density_u3.cpp -o $(OBJS_DIR)compute_scaled_density_u3.o
$(OBJS_DIR)compute_density_u_controlvolumes.o:$(COUPL_DIR)compute_density_u_controlvolumes.cpp 
	$(CXX) $(CFLAGS) $(COUPL_DIR)compute_density_u_controlvolumes.cpp -o $(OBJS_DIR)compute_density_u_controlvolumes.o

#----------------------#
#  utility functions   #
#----------------------#

$(OBJS_DIR)utilities.o:$(UTILS_DIR)utilities.cpp 
	$(CXX) $(CFLAGS) $(UTILS_DIR)utilities.cpp -o $(OBJS_DIR)utilities.o

#----------------------------------#
#  momentum equation functions     #
#----------------------------------#
$(OBJS_DIR)apply_boundary_conditions_velocity.o:$(MOMENTUM_DIR)apply_boundary_conditions_velocity.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)apply_boundary_conditions_velocity.cpp -o $(OBJS_DIR)apply_boundary_conditions_velocity.o
$(OBJS_DIR)apply_boundary_conditions_velocity_u1.o:$(MOMENTUM_DIR)apply_boundary_conditions_velocity_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)apply_boundary_conditions_velocity_u1.cpp -o $(OBJS_DIR)apply_boundary_conditions_velocity_u1.o
$(OBJS_DIR)apply_boundary_conditions_velocity_u2.o:$(MOMENTUM_DIR)apply_boundary_conditions_velocity_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)apply_boundary_conditions_velocity_u2.cpp -o $(OBJS_DIR)apply_boundary_conditions_velocity_u2.o
$(OBJS_DIR)apply_boundary_conditions_velocity_u3.o:$(MOMENTUM_DIR)apply_boundary_conditions_velocity_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)apply_boundary_conditions_velocity_u3.cpp -o $(OBJS_DIR)apply_boundary_conditions_velocity_u3.o
$(OBJS_DIR)build_momentum_matrix_u1.o:$(MOMENTUM_DIR)build_momentum_matrix_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_matrix_u1.cpp -o $(OBJS_DIR)build_momentum_matrix_u1.o
$(OBJS_DIR)build_momentum_matrix_u2.o:$(MOMENTUM_DIR)build_momentum_matrix_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_matrix_u2.cpp -o $(OBJS_DIR)build_momentum_matrix_u2.o
$(OBJS_DIR)build_momentum_matrix_u3.o:$(MOMENTUM_DIR)build_momentum_matrix_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_matrix_u3.cpp -o $(OBJS_DIR)build_momentum_matrix_u3.o
$(OBJS_DIR)build_momentum_predictor_u1.o:$(MOMENTUM_DIR)build_momentum_predictor_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_predictor_u1.cpp -o $(OBJS_DIR)build_momentum_predictor_u1.o
$(OBJS_DIR)build_momentum_predictor_u2.o:$(MOMENTUM_DIR)build_momentum_predictor_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_predictor_u2.cpp -o $(OBJS_DIR)build_momentum_predictor_u2.o
$(OBJS_DIR)build_momentum_predictor_u3.o:$(MOMENTUM_DIR)build_momentum_predictor_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_predictor_u3.cpp -o $(OBJS_DIR)build_momentum_predictor_u3.o
$(OBJS_DIR)build_momentum_rhs_u1.o:$(MOMENTUM_DIR)build_momentum_rhs_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_rhs_u1.cpp -o $(OBJS_DIR)build_momentum_rhs_u1.o
$(OBJS_DIR)build_momentum_rhs_u2.o:$(MOMENTUM_DIR)build_momentum_rhs_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_rhs_u2.cpp -o $(OBJS_DIR)build_momentum_rhs_u2.o
$(OBJS_DIR)build_momentum_rhs_u3.o:$(MOMENTUM_DIR)build_momentum_rhs_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)build_momentum_rhs_u3.cpp -o $(OBJS_DIR)build_momentum_rhs_u3.o
$(OBJS_DIR)compute_convection_term.o:$(MOMENTUM_DIR)compute_convection_term.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)compute_convection_term.cpp -o $(OBJS_DIR)compute_convection_term.o
$(OBJS_DIR)compute_scaled_density.o:$(MOMENTUM_DIR)compute_scaled_density.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)compute_scaled_density.cpp -o $(OBJS_DIR)compute_scaled_density.o
$(OBJS_DIR)compute_scaled_viscosity.o:$(MOMENTUM_DIR)compute_scaled_viscosity.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)compute_scaled_viscosity.cpp -o $(OBJS_DIR)compute_scaled_viscosity.o
$(OBJS_DIR)compute_weighted_average.o:$(MOMENTUM_DIR)compute_weighted_average.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)compute_weighted_average.cpp -o $(OBJS_DIR)compute_weighted_average.o
$(OBJS_DIR)fold_momentum_matrix_u1.o:$(MOMENTUM_DIR)fold_momentum_matrix_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)fold_momentum_matrix_u1.cpp -o $(OBJS_DIR)fold_momentum_matrix_u1.o
$(OBJS_DIR)fold_momentum_matrix_u2.o:$(MOMENTUM_DIR)fold_momentum_matrix_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)fold_momentum_matrix_u2.cpp -o $(OBJS_DIR)fold_momentum_matrix_u2.o
$(OBJS_DIR)fold_momentum_matrix_u3.o:$(MOMENTUM_DIR)fold_momentum_matrix_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)fold_momentum_matrix_u3.cpp -o $(OBJS_DIR)fold_momentum_matrix_u3.o
$(OBJS_DIR)fold_momentum_rhside_u1.o:$(MOMENTUM_DIR)fold_momentum_rhside_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)fold_momentum_rhside_u1.cpp -o $(OBJS_DIR)fold_momentum_rhside_u1.o
$(OBJS_DIR)fold_momentum_rhside_u2.o:$(MOMENTUM_DIR)fold_momentum_rhside_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)fold_momentum_rhside_u2.cpp -o $(OBJS_DIR)fold_momentum_rhside_u2.o
$(OBJS_DIR)fold_momentum_rhside_u3.o:$(MOMENTUM_DIR)fold_momentum_rhside_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)fold_momentum_rhside_u3.cpp -o $(OBJS_DIR)fold_momentum_rhside_u3.o
$(OBJS_DIR)shift_velocity_field.o:$(MOMENTUM_DIR)shift_velocity_field.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)shift_velocity_field.cpp -o $(OBJS_DIR)shift_velocity_field.o
$(OBJS_DIR)solve_momentum_predictor.o:$(MOMENTUM_DIR)solve_momentum_predictor.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)solve_momentum_predictor.cpp -o $(OBJS_DIR)solve_momentum_predictor.o
$(OBJS_DIR)solve_momentum_predictor_u1.o:$(MOMENTUM_DIR)solve_momentum_predictor_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)solve_momentum_predictor_u1.cpp -o $(OBJS_DIR)solve_momentum_predictor_u1.o
$(OBJS_DIR)solve_momentum_predictor_u2.o:$(MOMENTUM_DIR)solve_momentum_predictor_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)solve_momentum_predictor_u2.cpp -o $(OBJS_DIR)solve_momentum_predictor_u2.o
$(OBJS_DIR)solve_momentum_predictor_u3.o:$(MOMENTUM_DIR)solve_momentum_predictor_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)solve_momentum_predictor_u3.cpp -o $(OBJS_DIR)solve_momentum_predictor_u3.o
$(OBJS_DIR)map_index_u1.o:$(MOMENTUM_DIR)map_index_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)map_index_u1.cpp -o $(OBJS_DIR)map_index_u1.o
$(OBJS_DIR)map_index_u2.o:$(MOMENTUM_DIR)map_index_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)map_index_u2.cpp -o $(OBJS_DIR)map_index_u2.o
$(OBJS_DIR)map_index_u3.o:$(MOMENTUM_DIR)map_index_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)map_index_u3.cpp -o $(OBJS_DIR)map_index_u3.o
$(OBJS_DIR)copy_general_field.o:$(MOMENTUM_DIR)copy_general_field.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)copy_general_field.cpp -o $(OBJS_DIR)copy_general_field.o
$(OBJS_DIR)decompress_solution_velocity_u1.o:$(MOMENTUM_DIR)decompress_solution_velocity_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)decompress_solution_velocity_u1.cpp -o $(OBJS_DIR)decompress_solution_velocity_u1.o
$(OBJS_DIR)decompress_solution_velocity_u2.o:$(MOMENTUM_DIR)decompress_solution_velocity_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)decompress_solution_velocity_u2.cpp -o $(OBJS_DIR)decompress_solution_velocity_u2.o
$(OBJS_DIR)decompress_solution_velocity_u3.o:$(MOMENTUM_DIR)decompress_solution_velocity_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)decompress_solution_velocity_u3.cpp -o $(OBJS_DIR)decompress_solution_velocity_u3.o
$(OBJS_DIR)compress_solution_velocity_u1.o:$(MOMENTUM_DIR)compress_solution_velocity_u1.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)compress_solution_velocity_u1.cpp -o $(OBJS_DIR)compress_solution_velocity_u1.o
$(OBJS_DIR)compress_solution_velocity_u2.o:$(MOMENTUM_DIR)compress_solution_velocity_u2.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)compress_solution_velocity_u2.cpp -o $(OBJS_DIR)compress_solution_velocity_u2.o
$(OBJS_DIR)compress_solution_velocity_u3.o:$(MOMENTUM_DIR)compress_solution_velocity_u3.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)compress_solution_velocity_u3.cpp -o $(OBJS_DIR)compress_solution_velocity_u3.o
$(OBJS_DIR)output_predictor_velocityfield.o:$(MOMENTUM_DIR)output_predictor_velocityfield.cpp 
	$(CXX) $(CFLAGS) $(MOMENTUM_DIR)output_predictor_velocityfield.cpp -o $(OBJS_DIR)output_predictor_velocityfield.o


#---------------------------------#
#  pressure_correction_equation   #
#---------------------------------#

$(OBJS_DIR)apply_pressure_correction.o:$(PRESS_DIR)apply_pressure_correction.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)apply_pressure_correction.cpp -o $(OBJS_DIR)apply_pressure_correction.o
$(OBJS_DIR)build_pressure_matrix.o:$(PRESS_DIR)build_pressure_matrix.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)build_pressure_matrix.cpp -o $(OBJS_DIR)build_pressure_matrix.o
$(OBJS_DIR)build_pressure_rhs_boundary.o:$(PRESS_DIR)build_pressure_rhs_boundary.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)build_pressure_rhs_boundary.cpp -o $(OBJS_DIR)build_pressure_rhs_boundary.o
$(OBJS_DIR)build_pressure_rhs_momentum_part.o:$(PRESS_DIR)build_pressure_rhs_momentum_part.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)build_pressure_rhs_momentum_part.cpp -o $(OBJS_DIR)build_pressure_rhs_momentum_part.o
$(OBJS_DIR)build_pressure_system.o:$(PRESS_DIR)build_pressure_system.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)build_pressure_system.cpp -o $(OBJS_DIR)build_pressure_system.o
$(OBJS_DIR)set_pressure_boundary_condition.o:$(PRESS_DIR)set_pressure_boundary_condition.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)set_pressure_boundary_condition.cpp -o $(OBJS_DIR)set_pressure_boundary_condition.o
$(OBJS_DIR)solve_momentum_corrector.o:$(PRESS_DIR)solve_momentum_corrector.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)solve_momentum_corrector.cpp -o $(OBJS_DIR)solve_momentum_corrector.o
$(OBJS_DIR)solve_pressure_correction_system.o:$(PRESS_DIR)solve_pressure_correction_system.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)solve_pressure_correction_system.cpp -o $(OBJS_DIR)solve_pressure_correction_system.o
$(OBJS_DIR)map_index_pressure.o:$(PRESS_DIR)map_index_pressure.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)map_index_pressure.cpp -o $(OBJS_DIR)map_index_pressure.o
$(OBJS_DIR)decompress_solution_pressure.o:$(PRESS_DIR)decompress_solution_pressure.cpp
	$(CXX) $(CFLAGS) $(PRESS_DIR)decompress_solution_pressure.cpp -o $(OBJS_DIR)decompress_solution_pressure.o
$(OBJS_DIR)compress_solution_pressure.o:$(PRESS_DIR)compress_solution_pressure.cpp
	$(CXX) $(CFLAGS) $(PRESS_DIR)compress_solution_pressure.cpp -o $(OBJS_DIR)compress_solution_pressure.o
$(OBJS_DIR)field_extrapolate_boundary.o:$(PRESS_DIR)field_extrapolate_boundary.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)field_extrapolate_boundary.cpp -o $(OBJS_DIR)field_extrapolate_boundary.o
$(OBJS_DIR)project_pressure_rhside.o:$(PRESS_DIR)project_pressure_rhside.cpp 
	$(CXX) $(CFLAGS) $(PRESS_DIR)project_pressure_rhside.cpp -o $(OBJS_DIR)project_pressure_rhside.o
$(OBJS_DIR)shift_pressure_solution.o:$(PRESS_DIR)shift_pressure_solution.cpp
	$(CXX) $(CFLAGS) $(PRESS_DIR)shift_pressure_solution.cpp -o $(OBJS_DIR)shift_pressure_solution.o
$(OBJS_DIR)apply_boundary_conditions_pressure.o:$(PRESS_DIR)apply_boundary_conditions_pressure.cpp
	$(CXX) $(CFLAGS) $(PRESS_DIR)apply_boundary_conditions_pressure.cpp -o $(OBJS_DIR)apply_boundary_conditions_pressure.o
$(OBJS_DIR)dump_divergence_for_debugging.o:$(PRESS_DIR)dump_divergence_for_debugging.cpp
	$(CXX) $(CFLAGS) $(PRESS_DIR)dump_divergence_for_debugging.cpp -o $(OBJS_DIR)dump_divergence_for_debugging.o

#----------------------------#
#  linear solver functions   #
#----------------------------#

$(OBJS_DIR)apply_preconditioner.o:$(LINEAR_DIR)apply_preconditioner.cpp 
	$(CXX) $(CFLAGS) $(LINEAR_DIR)apply_preconditioner.cpp -o $(OBJS_DIR)apply_preconditioner.o
$(OBJS_DIR)build_preconditioner.o:$(LINEAR_DIR)build_preconditioner.cpp 
	$(CXX) $(CFLAGS) $(LINEAR_DIR)build_preconditioner.cpp -o $(OBJS_DIR)build_preconditioner.o
$(OBJS_DIR)conjugate_gradient_method.o:$(LINEAR_DIR)conjugate_gradient_method.cpp 
	$(CXX) $(CFLAGS) $(LINEAR_DIR)conjugate_gradient_method.cpp -o $(OBJS_DIR)conjugate_gradient_method.o
$(OBJS_DIR)export_matrix_matlab.o:$(LINEAR_DIR)export_matrix_matlab.cpp 
	$(CXX) $(CFLAGS) $(LINEAR_DIR)export_matrix_matlab.cpp -o $(OBJS_DIR)export_matrix_matlab.o

#------------------------#
#  interface functions   #
#------------------------#

$(OBJS_DIR)compute_volume_of_fluid.o:$(INTERFACE_DIR)compute_volume_of_fluid.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_volume_of_fluid.cpp -o $(OBJS_DIR)compute_volume_of_fluid.o
$(OBJS_DIR)compute_level_set_gradient.o:$(INTERFACE_DIR)compute_level_set_gradient.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_level_set_gradient.cpp -o $(OBJS_DIR)compute_level_set_gradient.o
$(OBJS_DIR)level_set_2_vof.o:$(INTERFACE_DIR)level_set_2_vof.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)level_set_2_vof.cpp -o $(OBJS_DIR)level_set_2_vof.o
$(OBJS_DIR)advance_interface.o:$(INTERFACE_DIR)advance_interface.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)advance_interface.cpp -o $(OBJS_DIR)advance_interface.o
$(OBJS_DIR)advect_level_set.o:$(INTERFACE_DIR)advect_level_set.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)advect_level_set.cpp -o $(OBJS_DIR)advect_level_set.o
$(OBJS_DIR)apply_redistribution_volume_of_fluid_correction.o:$(INTERFACE_DIR)apply_redistribution_volume_of_fluid_correction.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)apply_redistribution_volume_of_fluid_correction.cpp -o $(OBJS_DIR)apply_redistribution_volume_of_fluid_correction.o
$(OBJS_DIR)apply_ridders_method.o:$(INTERFACE_DIR)apply_ridders_method.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)apply_ridders_method.cpp -o $(OBJS_DIR)apply_ridders_method.o
$(OBJS_DIR)apply_bisection.o:$(INTERFACE_DIR)apply_bisection.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)apply_bisection.cpp -o $(OBJS_DIR)apply_bisection.o
$(OBJS_DIR)apply_volume_of_fluid_clipping.o:$(INTERFACE_DIR)apply_volume_of_fluid_clipping.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)apply_volume_of_fluid_clipping.cpp -o $(OBJS_DIR)apply_volume_of_fluid_clipping.o
$(OBJS_DIR)apply_volume_of_fluid_redistribution.o:$(INTERFACE_DIR)apply_volume_of_fluid_redistribution.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)apply_volume_of_fluid_redistribution.cpp -o $(OBJS_DIR)apply_volume_of_fluid_redistribution.o
$(OBJS_DIR)check_volume_of_fluid.o:$(INTERFACE_DIR)check_volume_of_fluid.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)check_volume_of_fluid.cpp -o $(OBJS_DIR)check_volume_of_fluid.o
$(OBJS_DIR)compute_derivatives_level_set.o:$(INTERFACE_DIR)compute_derivatives_level_set.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_derivatives_level_set.cpp -o $(OBJS_DIR)compute_derivatives_level_set.o
$(OBJS_DIR)compute_redistribution_velocity.o:$(INTERFACE_DIR)compute_redistribution_velocity.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_redistribution_velocity.cpp -o $(OBJS_DIR)compute_redistribution_velocity.o
$(OBJS_DIR)compute_vof_flux_x1.o:$(INTERFACE_DIR)compute_vof_flux_x1.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_vof_flux_x1.cpp -o $(OBJS_DIR)compute_vof_flux_x1.o
$(OBJS_DIR)compute_vof_flux_x2.o:$(INTERFACE_DIR)compute_vof_flux_x2.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_vof_flux_x2.cpp -o $(OBJS_DIR)compute_vof_flux_x2.o
$(OBJS_DIR)compute_vof_flux_x3.o:$(INTERFACE_DIR)compute_vof_flux_x3.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_vof_flux_x3.cpp -o $(OBJS_DIR)compute_vof_flux_x3.o
$(OBJS_DIR)copy_cell_centered_field.o:$(INTERFACE_DIR)copy_cell_centered_field.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)copy_cell_centered_field.cpp -o $(OBJS_DIR)copy_cell_centered_field.o
$(OBJS_DIR)field_neumann_boundary.o:$(INTERFACE_DIR)field_neumann_boundary.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)field_neumann_boundary.cpp -o $(OBJS_DIR)field_neumann_boundary.o
$(OBJS_DIR)level_set_2_vof_phi_negative.o:$(INTERFACE_DIR)level_set_2_vof_phi_negative.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)level_set_2_vof_phi_negative.cpp -o $(OBJS_DIR)level_set_2_vof_phi_negative.o
$(OBJS_DIR)make_level_set_mass_conserving.o:$(INTERFACE_DIR)make_level_set_mass_conserving.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)make_level_set_mass_conserving.cpp -o $(OBJS_DIR)make_level_set_mass_conserving.o
$(OBJS_DIR)match_level_set_to_volume_of_fluid.o:$(INTERFACE_DIR)match_level_set_to_volume_of_fluid.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)match_level_set_to_volume_of_fluid.cpp -o $(OBJS_DIR)match_level_set_to_volume_of_fluid.o
$(OBJS_DIR)modify_volume_of_fluid_values.o:$(INTERFACE_DIR)modify_volume_of_fluid_values.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)modify_volume_of_fluid_values.cpp -o $(OBJS_DIR)modify_volume_of_fluid_values.o
$(OBJS_DIR)reinitialize_level_set.o:$(INTERFACE_DIR)reinitialize_level_set.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)reinitialize_level_set.cpp -o $(OBJS_DIR)reinitialize_level_set.o
$(OBJS_DIR)update_volume_of_fluid_flux_x1.o:$(INTERFACE_DIR)update_volume_of_fluid_flux_x1.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)update_volume_of_fluid_flux_x1.cpp -o $(OBJS_DIR)update_volume_of_fluid_flux_x1.o
$(OBJS_DIR)update_volume_of_fluid_flux_x2.o:$(INTERFACE_DIR)update_volume_of_fluid_flux_x2.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)update_volume_of_fluid_flux_x2.cpp -o $(OBJS_DIR)update_volume_of_fluid_flux_x2.o
$(OBJS_DIR)update_volume_of_fluid_flux_x3.o:$(INTERFACE_DIR)update_volume_of_fluid_flux_x3.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)update_volume_of_fluid_flux_x3.cpp -o $(OBJS_DIR)update_volume_of_fluid_flux_x3.o
$(OBJS_DIR)upwind_flux_mass_redistribution.o:$(INTERFACE_DIR)upwind_flux_mass_redistribution.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)upwind_flux_mass_redistribution.cpp -o $(OBJS_DIR)upwind_flux_mass_redistribution.o
$(OBJS_DIR)vof_2_level_set.o:$(INTERFACE_DIR)vof_2_level_set.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)vof_2_level_set.cpp -o $(OBJS_DIR)vof_2_level_set.o
$(OBJS_DIR)compute_level_set_flux.o:$(INTERFACE_DIR)compute_level_set_flux.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_level_set_flux.cpp -o $(OBJS_DIR)compute_level_set_flux.o
$(OBJS_DIR)shift_interface.o:$(INTERFACE_DIR)shift_interface.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)shift_interface.cpp -o $(OBJS_DIR)shift_interface.o
$(OBJS_DIR)compute_normal_derivative_at_faces.o:$(INTERFACE_DIR)compute_normal_derivative_at_faces.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_normal_derivative_at_faces.cpp -o $(OBJS_DIR)compute_normal_derivative_at_faces.o
$(OBJS_DIR)redistribute_volume_of_fluid_error.o:$(INTERFACE_DIR)redistribute_volume_of_fluid_error.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)redistribute_volume_of_fluid_error.cpp -o $(OBJS_DIR)redistribute_volume_of_fluid_error.o
$(OBJS_DIR)dump_solution_for_debugging.o:$(INTERFACE_DIR)dump_solution_for_debugging.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)dump_solution_for_debugging.cpp -o $(OBJS_DIR)dump_solution_for_debugging.o
$(OBJS_DIR)dump_adapted_vof_for_debugging.o:$(INTERFACE_DIR)dump_adapted_vof_for_debugging.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)dump_adapted_vof_for_debugging.cpp -o $(OBJS_DIR)dump_adapted_vof_for_debugging.o
$(OBJS_DIR)compute_mass_in_domain.o:$(INTERFACE_DIR)compute_mass_in_domain.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_mass_in_domain.cpp -o $(OBJS_DIR)compute_mass_in_domain.o
$(OBJS_DIR)analyse_interface_properties.o:$(INTERFACE_DIR)analyse_interface_properties.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)analyse_interface_properties.cpp -o $(OBJS_DIR)analyse_interface_properties.o
$(OBJS_DIR)compute_vof_at_u1_points.o:$(INTERFACE_DIR)compute_vof_at_u1_points.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_vof_at_u1_points.cpp -o $(OBJS_DIR)compute_vof_at_u1_points.o
$(OBJS_DIR)compute_vof_at_u2_points.o:$(INTERFACE_DIR)compute_vof_at_u2_points.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_vof_at_u2_points.cpp -o $(OBJS_DIR)compute_vof_at_u2_points.o
$(OBJS_DIR)compute_vof_at_u3_points.o:$(INTERFACE_DIR)compute_vof_at_u3_points.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_vof_at_u3_points.cpp -o $(OBJS_DIR)compute_vof_at_u3_points.o
$(OBJS_DIR)dump_reinitialization_for_debugging.o:$(INTERFACE_DIR)dump_reinitialization_for_debugging.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)dump_reinitialization_for_debugging.cpp -o $(OBJS_DIR)dump_reinitialization_for_debugging.o
$(OBJS_DIR)advect_level_set_higher_order.o:$(INTERFACE_DIR)advect_level_set_higher_order.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)advect_level_set_higher_order.cpp -o $(OBJS_DIR)advect_level_set_higher_order.o
$(OBJS_DIR)compute_first_stage_RK.o:$(INTERFACE_DIR)compute_first_stage_RK.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_first_stage_RK.cpp -o $(OBJS_DIR)compute_first_stage_RK.o
$(OBJS_DIR)compute_second_stage_RK.o:$(INTERFACE_DIR)compute_second_stage_RK.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_second_stage_RK.cpp -o $(OBJS_DIR)compute_second_stage_RK.o
$(OBJS_DIR)compute_new_time_RK.o:$(INTERFACE_DIR)compute_new_time_RK.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)compute_new_time_RK.cpp -o $(OBJS_DIR)compute_new_time_RK.o
$(OBJS_DIR)evaluate_convection_operator_weno.o:$(INTERFACE_DIR)evaluate_convection_operator_weno.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)evaluate_convection_operator_weno.cpp -o $(OBJS_DIR)evaluate_convection_operator_weno.o
$(OBJS_DIR)weno_flux_computation.o:$(INTERFACE_DIR)weno_flux_computation.cpp 
	$(CXX) $(CFLAGS) $(INTERFACE_DIR)weno_flux_computation.cpp -o $(OBJS_DIR)weno_flux_computation.o

#------------------------------------#
#  output/postprocessing functions   #
#------------------------------------#


$(OBJS_DIR)interpolate_velocity_u1_center.o:$(OUTPUT_DIR)interpolate_velocity_u1_center.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)interpolate_velocity_u1_center.cpp -o $(OBJS_DIR)interpolate_velocity_u1_center.o
$(OBJS_DIR)interpolate_velocity_u2_center.o:$(OUTPUT_DIR)interpolate_velocity_u2_center.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)interpolate_velocity_u2_center.cpp -o $(OBJS_DIR)interpolate_velocity_u2_center.o
$(OBJS_DIR)interpolate_velocity_u3_center.o:$(OUTPUT_DIR)interpolate_velocity_u3_center.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)interpolate_velocity_u3_center.cpp -o $(OBJS_DIR)interpolate_velocity_u3_center.o
$(OBJS_DIR)interpolate_velocity_u1_vertex.o:$(OUTPUT_DIR)interpolate_velocity_u1_vertex.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)interpolate_velocity_u1_vertex.cpp -o $(OBJS_DIR)interpolate_velocity_u1_vertex.o
$(OBJS_DIR)interpolate_velocity_u2_vertex.o:$(OUTPUT_DIR)interpolate_velocity_u2_vertex.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)interpolate_velocity_u2_vertex.cpp -o $(OBJS_DIR)interpolate_velocity_u2_vertex.o
$(OBJS_DIR)interpolate_velocity_u3_vertex.o:$(OUTPUT_DIR)interpolate_velocity_u3_vertex.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)interpolate_velocity_u3_vertex.cpp -o $(OBJS_DIR)interpolate_velocity_u3_vertex.o
$(OBJS_DIR)output_solution.o:$(OUTPUT_DIR)output_solution.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)output_solution.cpp -o $(OBJS_DIR)output_solution.o
$(OBJS_DIR)write_cell_centered_field_tecplot.o:$(OUTPUT_DIR)write_cell_centered_field_tecplot.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)write_cell_centered_field_tecplot.cpp -o $(OBJS_DIR)write_cell_centered_field_tecplot.o
$(OBJS_DIR)write_cell_centered_field_vtk.o:$(OUTPUT_DIR)write_cell_centered_field_vtk.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)write_cell_centered_field_vtk.cpp -o $(OBJS_DIR)write_cell_centered_field_vtk.o
$(OBJS_DIR)write_vertex_centered_vector_field_vtk.o:$(OUTPUT_DIR)write_vertex_centered_vector_field_vtk.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)write_vertex_centered_vector_field_vtk.cpp -o $(OBJS_DIR)write_vertex_centered_vector_field_vtk.o
$(OBJS_DIR)write_vertex_centered_field_vtk.o:$(OUTPUT_DIR)write_vertex_centered_field_vtk.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)write_vertex_centered_field_vtk.cpp -o $(OBJS_DIR)write_vertex_centered_field_vtk.o
$(OBJS_DIR)write_coordinates_tecplot.o:$(OUTPUT_DIR)write_coordinates_tecplot.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)write_coordinates_tecplot.cpp -o $(OBJS_DIR)write_coordinates_tecplot.o
$(OBJS_DIR)write_coordinates_vtk.o:$(OUTPUT_DIR)write_coordinates_vtk.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)write_coordinates_vtk.cpp -o $(OBJS_DIR)write_coordinates_vtk.o
$(OBJS_DIR)check_symmetry_scalars.o:$(OUTPUT_DIR)check_symmetry_scalars.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)check_symmetry_scalars.cpp -o $(OBJS_DIR)check_symmetry_scalars.o
$(OBJS_DIR)check_symmetry_velocities.o:$(OUTPUT_DIR)check_symmetry_velocities.cpp 
	$(CXX) $(CFLAGS) $(OUTPUT_DIR)check_symmetry_velocities.cpp -o $(OBJS_DIR)check_symmetry_velocities.o


#------------------------------------#
#  restart functions   		     #
#------------------------------------#

$(OBJS_DIR)write_restart_file.o:$(RESTART_DIR)write_restart_file.cpp 
	$(CXX) $(CFLAGS) $(RESTART_DIR)write_restart_file.cpp -o $(OBJS_DIR)write_restart_file.o
$(OBJS_DIR)read_restart_file.o:$(RESTART_DIR)read_restart_file.cpp 
	$(CXX) $(CFLAGS) $(RESTART_DIR)read_restart_file.cpp -o $(OBJS_DIR)read_restart_file.o


clean:
	$(RM) $(OBJS) ../objects/funcdefs.h
# backup
# classes
# headers
# initialization
# interface :
# interface_coupling :
# linear_solver :
# main_program
# momentum_equation
# navier_stokes
# post_processing
# pressure_correction_equation
# utils

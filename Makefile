CFLAGS=-O3 -Wall

BUILD_ROOT=objects
EXECUTABLE_DIR=executable

COMMON_OBJS=\
	$(addprefix $(BUILD_ROOT)/initialization/,\
		allocate_main_variables.o\
		build_pressure_rhs_initialization.o\
		build_pressure_system_initialization.o\
		initialize_all_variables.o\
		initialize_bubbles.o\
		initialize_computation.o\
		initialize_coupling.o\
		initialize_curvature.o\
		initialize_flow_field.o\
		initialize_free_surface.o\
		initialize_interface.o\
		initialize_level_set.o\
		initialize_pressure.o\
		initialize_volume_of_fluid.o\
		set_boundary_conditions.o\
		set_parameters.o\
		write_interface_solution.o\
	)\
	$(addprefix $(BUILD_ROOT)/interface/,\
		advance_interface.o\
		advect_level_set.o\
		advect_level_set_higher_order.o\
		analyse_interface_properties.o\
		apply_bisection.o\
		apply_ridders_method.o\
		apply_volume_of_fluid_clipping.o\
		apply_volume_of_fluid_redistribution.o\
		check_volume_of_fluid.o\
		compute_derivatives_level_set.o\
		compute_first_stage_RK.o\
		compute_level_set_flux.o\
		compute_level_set_gradient.o\
		compute_mass_in_domain.o\
		compute_new_time_RK.o\
		compute_normal_derivative_at_faces.o\
		compute_redistribution_velocity.o\
		compute_second_stage_RK.o\
		compute_vof_at_u1_points.o\
		compute_vof_at_u2_points.o\
		compute_vof_at_u3_points.o\
		compute_vof_flux_x1.o\
		compute_vof_flux_x2.o\
		compute_vof_flux_x3.o\
		compute_volume_of_fluid.o\
		copy_cell_centered_field.o\
		dump_adapted_vof_for_debugging.o\
		dump_reinitialization_for_debugging.o\
		dump_solution_for_debugging.o\
		evaluate_convection_operator_weno.o\
		field_neumann_boundary.o\
		level_set_2_vof.o\
		level_set_2_vof_phi_negative.o\
		make_level_set_mass_conserving.o\
		match_level_set_to_volume_of_fluid.o\
		modify_volume_of_fluid_values.o\
		redistribute_volume_of_fluid_error.o\
		reinitialize_level_set.o\
		shift_interface.o\
		update_volume_of_fluid_flux_x1.o\
		update_volume_of_fluid_flux_x2.o\
		update_volume_of_fluid_flux_x3.o\
		upwind_flux_mass_redistribution.o\
		vof_2_level_set.o\
		weno_flux_computation.o\
	)\
	$(addprefix $(BUILD_ROOT)/interface_coupling/,\
		advance_coupling.o\
		advance_coupling_part1.o\
		advance_coupling_part2.o\
		compute_body_force_x1.o\
		compute_body_force_x2.o\
		compute_body_force_x3.o\
		compute_curvature.o\
		compute_curvature_error_laplace.o\
		compute_density_u_controlvolumes.o\
		compute_momentum_source_terms.o\
		compute_scaled_density_u1.o\
		compute_scaled_density_u2.o\
		compute_scaled_density_u3.o\
		compute_surface_tension_body_force.o\
		compute_weighted_curvature.o\
		computed_derivative_heaviside_function.o\
		curvature_filter.o\
		delta_function.o\
		dump_curvature_for_debugging.o\
		smooth_curvature.o\
	)\
	$(addprefix $(BUILD_ROOT)/linear_solver/,\
		apply_preconditioner.o\
		build_preconditioner.o\
		conjugate_gradient_method.o\
		export_matrix_matlab.o\
	)\
	$(addprefix $(BUILD_ROOT)/main_program/,\
		advance_flow_field.o\
		compute_time_step_size.o\
		time_stepping_sequence.o\
	)\
	$(addprefix $(BUILD_ROOT)/momentum_equation/,\
		apply_boundary_conditions_velocity.o\
		apply_boundary_conditions_velocity_u1.o\
		apply_boundary_conditions_velocity_u2.o\
		apply_boundary_conditions_velocity_u3.o\
		compute_scaled_density.o\
		compute_scaled_viscosity.o\
		compute_weighted_average.o\
		copy_general_field.o\
		momentum_predictor.o\
		output_predictor_velocityfield.o\
		shift_velocity_field.o\
	)\
	$(addprefix $(BUILD_ROOT)/post_processing/,\
		check_symmetry_scalars.o\
		check_symmetry_velocities.o\
		interpolate_velocity_u1_center.o\
		interpolate_velocity_u1_vertex.o\
		interpolate_velocity_u2_center.o\
		interpolate_velocity_u2_vertex.o\
		interpolate_velocity_u3_center.o\
		interpolate_velocity_u3_vertex.o\
		output_solution.o\
		write_cell_centered_field_tecplot.o\
		write_cell_centered_field_vtk.o\
		write_coordinates_tecplot.o\
		write_coordinates_vtk.o\
		write_vertex_centered_field_vtk.o\
		write_vertex_centered_vector_field_vtk.o\
	)\
	$(addprefix $(BUILD_ROOT)/pressure_correction_equation/,\
		apply_boundary_conditions_pressure.o\
		apply_pressure_correction.o\
		build_pressure_matrix.o\
		build_pressure_rhs_boundary.o\
		build_pressure_rhs_momentum_part.o\
		build_pressure_system.o\
		compress_solution_pressure.o\
		decompress_solution_pressure.o\
		dump_divergence_for_debugging.o\
		field_extrapolate_boundary.o\
		map_index_pressure.o\
		project_pressure_rhside.o\
		set_pressure_boundary_condition.o\
		shift_pressure_solution.o\
		solve_momentum_corrector.o\
		solve_pressure_correction_system.o\
	)\
	$(addprefix $(BUILD_ROOT)/restart/,\
		read_restart_file.o\
		write_restart_file.o\
	)\
	$(addprefix $(BUILD_ROOT)/utils/,\
		utilities.o\
	)

ALL_TARGETS:=$(COMMON_OBJS)


# MCLS

MCLS_OBJS=$(BUILD_ROOT)/main_program/dns.o
ALL_TARGETS:=$(ALL_TARGETS) $(MCLS_OBJS) $(EXECUTABLE_DIR)/MCLS

$(EXECUTABLE_DIR)/MCLS: $(COMMON_OBJS) $(MCLS_OBJS) | target_dirs
	$(CXX) $^ -o $@

.PHONY: MCLS
MCLS: $(EXECUTABLE_DIR)/MCLS


# TODO: unit tests
# momentum_predictor_unit_test.o
# advance_flow_field_unit_test.o
# add_arrays_unit_test.o


# common rules

ALL_TARGETS:=$(ALL_TARGETS) $(BUILD_ROOT)/funcdefs.h
$(BUILD_ROOT)/funcdefs.h: | target_dirs
	./gen_funcdefs.h --output $@ -- $(patsubst $(BUILD_ROOT)/%.o,src/%.cpp,$(COMMON_OBJS))

.PHONY: target_dirs
target_dirs:
	mkdir -p $(sort $(dir $(ALL_TARGETS)))

.PHONY: clean
clean:
	$(RM) $(ALL_TARGETS)

$(BUILD_ROOT)/%.o: src/%.c $(BUILD_ROOT)/funcdefs.h | target_dirs
	$(CC) -c -include $(BUILD_ROOT)/funcdefs.h $(CFLAGS) $< -o $@

$(BUILD_ROOT)/%.o: src/%.cpp $(BUILD_ROOT)/funcdefs.h | target_dirs
	$(CXX) -c -include $(BUILD_ROOT)/funcdefs.h $(CFLAGS) $< -o $@

# vim: noet

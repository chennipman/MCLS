 CFLAGS=-O3 -Wall 
 CXXFLAGS=-O3 -Wall 
# CC=gcc
# CXX=g++
# CFLAGS=-g -static -Wall
# CXXFLAGS=-g -static -Wall

BUILD_DIR=objects

COMMIT_NUMBER:=$(shell git log -n 1 | head -1 | sed 's/commit //')
NOW:=$(shell date +%m_%d_%H_%M_%S)
RUN_DIR:=$(NOW)_$(COMMIT_NUMBER)

COMMON_SRCS= \
	$(addprefix src/initialization/, \
		allocate_main_variables.cpp \
		build_pressure_rhs_initialization.cpp \
		build_pressure_system_initialization.cpp \
		initialize_all_variables.cpp \
		initialize_bubbles.cpp \
		initialize_computation.cpp \
		initialize_coupling.cpp \
		initialize_curvature.cpp \
		initialize_free_surface.cpp \
		initialize_interface.cpp \
		initialize_level_set.cpp \
		initialize_pressure.cpp \
		initialize_rayleigh_taylor.cpp \
		initialize_volume_of_fluid.cpp \
		write_interface_solution.cpp \
	) \
	$(addprefix src/interface/, \
		advance_interface.cpp \
		advect_level_set.cpp \
		advect_level_set_higher_order.cpp \
		analyse_interface_properties.cpp \
		analyse_validity_vof_correction.cpp \
		analyse_validity_vof_field.cpp \
		apply_bisection.cpp \
		apply_ridders_method.cpp \
		apply_volume_of_fluid_clipping.cpp \
		apply_volume_of_fluid_redistribution.cpp \
		check_volume_of_fluid.cpp \
		compute_derivatives_level_set.cpp \
		compute_first_stage_RK.cpp \
		compute_level_set_flux.cpp \
		compute_level_set_gradient.cpp \
		compute_mass_in_domain.cpp \
		compute_new_time_RK.cpp \
		compute_normal_derivative_at_faces.cpp \
		compute_redistribution_time_derivative.cpp \
		compute_redistribution_velocity.cpp \
		compute_redistribution_velocity_field.cpp \
		compute_second_stage_RK.cpp \
		compute_vof_at_u1_points.cpp \
		compute_vof_at_u2_points.cpp \
		compute_vof_at_u3_points.cpp \
		compute_vof_flux_x1.cpp \
		compute_vof_flux_x2.cpp \
		compute_vof_flux_x3.cpp \
		compute_volume_of_fluid.cpp \
		copy_cell_centered_field.cpp \
		determine_cell_is_pure.cpp \
		dump_adapted_vof_for_debugging.cpp \
		dump_reinitialization_for_debugging.cpp \
		dump_solution_for_debugging.cpp \
		dump_redistribution_for_debugging.cpp \
		evaluate_convection_operator_weno.cpp \
		field_neumann_boundary.cpp \
		level_set_2_vof.cpp \
		level_set_2_vof_phi_negative.cpp \
		make_level_set_mass_conserving.cpp \
		make_vof_field_valid.cpp \
		match_level_set_to_volume_of_fluid.cpp \
		modify_volume_of_fluid_values.cpp \
		redistribute_volume_of_fluid_error.cpp \
		reinitialize_level_set.cpp \
		shift_interface.cpp \
		update_volume_of_fluid_flux_x1.cpp \
		update_volume_of_fluid_flux_x2.cpp \
		update_volume_of_fluid_flux_x3.cpp \
		upwind_flux_mass_redistribution.cpp \
		vof_2_level_set.cpp \
		weno_flux_computation.cpp \
	) \
	$(addprefix src/interface_coupling/, \
		advance_coupling.cpp \
		advance_coupling_part1.cpp \
		advance_coupling_part2.cpp \
		compute_body_force_x1.cpp \
		compute_body_force_x2.cpp \
		compute_body_force_x3.cpp \
		compute_curvature.cpp \
		compute_curvature_error_laplace.cpp \
		compute_density_u_controlvolumes.cpp \
		compute_momentum_source_terms.cpp \
		compute_scaled_density_u1.cpp \
		compute_scaled_density_u2.cpp \
		compute_scaled_density_u3.cpp \
		compute_surface_tension_body_force.cpp \
		compute_weighted_curvature.cpp \
		computed_derivative_heaviside_function.cpp \
		curvature_filter.cpp \
		delta_function.cpp \
		dump_curvature_for_debugging.cpp \
		smooth_curvature.cpp \
	) \
	$(addprefix src/linear_solver/, \
		apply_preconditioner.cpp \
		build_preconditioner.cpp \
		conjugate_gradient_method.cpp \
		export_matrix_matlab.cpp \
	) \
	$(addprefix src/main_program/, \
		advance_flow_field.cpp \
		compute_time_step_size.cpp \
		time_stepping_sequence.cpp \
	) \
	$(addprefix src/momentum_equation/, \
		apply_boundary_conditions_velocity.cpp \
		apply_boundary_conditions_velocity_u1.cpp \
		apply_boundary_conditions_velocity_u2.cpp \
		apply_boundary_conditions_velocity_u3.cpp \
		compute_scaled_density.cpp \
		compute_scaled_viscosity.cpp \
		compute_weighted_average.cpp \
		copy_general_field.cpp \
		output_predictor_velocityfield.cpp \
		shift_velocity_field.cpp \
		build_momentum_matrix_u1.cpp \
		build_momentum_matrix_u2.cpp \
		build_momentum_matrix_u3.cpp \
		build_momentum_predictor_u1.cpp \
		build_momentum_predictor_u2.cpp \
		build_momentum_predictor_u3.cpp \
		build_momentum_rhs_u1.cpp \
		build_momentum_rhs_u2.cpp \
		build_momentum_rhs_u3.cpp \
		compress_solution_velocity_u1.cpp \
		compress_solution_velocity_u2.cpp \
		compress_solution_velocity_u3.cpp \
		compute_convection_term.cpp \
		copy_general_field.cpp \
		decompress_solution_velocity_u1.cpp \
		decompress_solution_velocity_u2.cpp \
		decompress_solution_velocity_u3.cpp \
		fold_momentum_matrix_u1.cpp \
		fold_momentum_matrix_u2.cpp \
		fold_momentum_matrix_u3.cpp \
		fold_momentum_rhside_u1.cpp \
		fold_momentum_rhside_u2.cpp \
		fold_momentum_rhside_u3.cpp \
		map_index_u1.cpp \
		map_index_u2.cpp \
		map_index_u3.cpp \
		solve_momentum_predictor_imex.cpp \
		solve_momentum_predictor_rk.cpp \
		solve_momentum_predictor_explicit.cpp \
		forward_euler.cpp \
		convection_diffussion_source_terms.cpp \
		solve_momentum_predictor_u1.cpp \
		solve_momentum_predictor_u2.cpp \
		solve_momentum_predictor_u3.cpp \
		dump_to_check_pressure.cpp \
	) \
	$(addprefix src/post_processing/, \
		check_symmetry_scalars.cpp \
		check_symmetry_velocities.cpp \
		interpolate_velocity_u1_center.cpp \
		interpolate_velocity_u1_vertex.cpp \
		interpolate_velocity_u2_center.cpp \
		interpolate_velocity_u2_vertex.cpp \
		interpolate_velocity_u3_center.cpp \
		interpolate_velocity_u3_vertex.cpp \
		output_solution.cpp \
		write_cell_centered_field_tecplot.cpp \
		write_cell_centered_field_vtk.cpp \
		write_coordinates_tecplot.cpp \
		write_coordinates_vtk.cpp \
		write_vertex_centered_field_vtk.cpp \
		write_vertex_centered_vector_field_vtk.cpp \
		write_face_field_pure.cpp \
	) \
	$(addprefix src/pressure_correction_equation/, \
		apply_boundary_conditions_pressure.cpp \
		apply_pressure_correction.cpp \
		build_pressure_matrix.cpp \
		build_pressure_rhs_boundary.cpp \
		build_pressure_rhs_momentum_part.cpp \
		build_pressure_system.cpp \
		build_pressure_system_final.cpp \
		compress_solution_pressure.cpp \
		decompress_solution_pressure.cpp \
		dump_divergence_for_debugging.cpp \
		field_extrapolate_boundary.cpp \
		map_index_pressure.cpp \
		project_pressure_rhside.cpp \
		set_pressure_boundary_condition.cpp \
		shift_pressure_solution.cpp \
		solve_momentum_corrector.cpp \
		solve_pressure_correction_system.cpp \
		solve_final_pressure_corrector.cpp \
	) \
	$(addprefix src/restart/, \
		read_restart_file.cpp \
		write_restart_file.cpp \
	) \
	$(addprefix src/utils/, \
		utilities.cpp \
	) 	
COMMON_OBJS=$(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(basename $(COMMON_SRCS))))
ALL_TARGETS:=$(COMMON_OBJS)
EXECUTABLE_DIR=executable

# MCLS(default)

CASE_DEFAULT=undefined_case
MCLS_OBJS=$(BUILD_DIR)/src/main_program/dns.o
ALL_TARGETS:=$(ALL_TARGETS) $(MCLS_OBJS) $(EXECUTABLE_DIR)/$(CASE_DEFAULT)/$(RUN_DIR)/MCLS

.PHONY: make

make: copy_case_$(CASE_DEFAULT) $(EXECUTABLE_DIR)/$(CASE_DEFAULT)/$(RUN_DIR)/MCLS

$(EXECUTABLE_DIR)/$(CASE_DEFAULT)/$(RUN_DIR)/MCLS: $(BUILD_DIR)/testcases/$(CASE_DEFAULT)/set_boundary_conditions.o $(BUILD_DIR)/testcases/$(CASE_DEFAULT)/set_parameters.o $(BUILD_DIR)/testcases/$(CASE_DEFAULT)/initialize_flow_field.o $(COMMON_OBJS) $(MCLS_OBJS) 
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) $^ -o $@  
	
	cp testcases/$(CASE_DEFAULT)/set_boundary_conditions.cpp          $(EXECUTABLE_DIR)/$(CASE_DEFAULT)/$(RUN_DIR)
	cp testcases/$(CASE_DEFAULT)/set_parameters.cpp                   $(EXECUTABLE_DIR)/$(CASE_DEFAULT)/$(RUN_DIR)
	cp testcases/$(CASE_DEFAULT)/initialize_flow_field.cpp            $(EXECUTABLE_DIR)/$(CASE_DEFAULT)/$(RUN_DIR)
	
	@test -z "`git status --porcelain`" || echo WARNING: there are uncommited changes

.PHONY: MCLS
MCLS:
	$(EXECUTABLE_DIR)/MCLS

# UNIT TESTS

UNIT_TESTS=$(addprefix test_, $(notdir $(basename $(wildcard unittest/*.cpp))))

ALL_TARGETS:= \
	$(ALL_TARGETS) \
	$(addprefix $(EXECUTABLE_DIR)/, $(UNIT_TESTS)) \
	$(addprefix $(BUILD_DIR)/unittest/, $(addsuffix .o, $(UNIT_TESTS)))

$(EXECUTABLE_DIR)/test_%: $(BUILD_DIR)/unittest/%.o $(COMMON_OBJS) 
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) $^ -o $@

.PHONY: test $(UNIT_TESTS)

test: $(UNIT_TESTS)

$(UNIT_TESTS): $(addprefix $(EXECUTABLE_DIR)/,$(UNIT_TESTS))
	$(EXECUTABLE_DIR)/$@


# MAKE TEST CASES

define CASE_RULES

ALL_TARGETS:= \
	$(ALL_TARGETS) \
	$(BUILD_DIR)/testcases/$(1)/set_boundary_conditions.o \
	$(BUILD_DIR)/testcases/$(1)/set_parameters.o \
	$(BUILD_DIR)/testcases/$(1)/initialize_flow_field.o
	


$(EXECUTABLE_DIR)/case_$(1)/$$(RUN_DIR)/MCLS: $$(BUILD_DIR)/testcases/$(1)/set_boundary_conditions.o $$(BUILD_DIR)/testcases/$(1)/set_parameters.o $$(BUILD_DIR)/testcases/$(1)/initialize_flow_field.o $$(COMMON_OBJS) $$(BUILD_DIR)/src/main_program/dns.o
	@mkdir -p $$(dir $$@)
	$$(CXX) $$(LDFLAGS) $$^ -o $$@
	cp testcases/$(1)/set_boundary_conditions.cpp $(EXECUTABLE_DIR)/case_$(1)/$$(RUN_DIR)/
	cp testcases/$(1)/set_parameters.cpp          $(EXECUTABLE_DIR)/case_$(1)/$$(RUN_DIR)/
	cp testcases/$(1)/initialize_flow_field.cpp   $(EXECUTABLE_DIR)/case_$(1)/$$(RUN_DIR)/
	
	@test -z "`git status --porcelain`" || echo WARNING: there are uncommited changes

.PHONY: case_$(1)

case_$(1): copy_case_$(1) $(EXECUTABLE_DIR)/case_$(1)/$$(RUN_DIR)/MCLS 
#	$$^  #can be used to directly run the executable

endef

TEST_CASES=$(notdir $(wildcard testcases/*))
$(foreach case,$(TEST_CASES),$(eval $(call CASE_RULES,$(case))))


# COPY TEST CASES FROM DEFAULT TO NORMAL

define COPY_RULES

.PHONY: copy_case_$(1)

copy_case_$(1): 
	cp -n testcases/$(1)/default_set_boundary_conditions.cpp testcases/$(1)/set_boundary_conditions.cpp
	cp -n testcases/$(1)/default_set_parameters.cpp          testcases/$(1)/set_parameters.cpp 
	cp -n testcases/$(1)/default_initialize_flow_field.cpp   testcases/$(1)/initialize_flow_field.cpp 


#copy_case_$(1): case_$(1)_copy
#	$$^  #can be used to directly run the executable

endef
$(foreach case,$(TEST_CASES),$(eval $(call COPY_RULES,$(case))))


# COMMON RULES

ALL_TARGETS:=$(ALL_TARGETS) $(BUILD_DIR)/funcdefs.h

$(BUILD_DIR)/funcdefs.h:
	@mkdir -p $(dir $@)
#	./gen_funcdefs.h --output $@ -- $(COMMON_SRCS) testcases/TV/initialize_flow_field.cpp testcases/TV/set_boundary_conditions.cpp testcases/TV/set_parameters.cpp 
	./gen_funcdefs.h --output $@ -- $(COMMON_SRCS) testcases/TV/default_initialize_flow_field.cpp testcases/TV/default_set_boundary_conditions.cpp testcases/TV/default_set_parameters.cpp 
	#the TV case is used, because this is the only case with additional function definitions in the set_boundary_conditions 
	#if other function definitions are used the 'TV', should be replace by the then used function definition
	#the default case set is used because then this part is directly working after pulling the repository 

.PHONY: clean
clean:
	$(RM) $(ALL_TARGETS)

$(BUILD_DIR)/%.o: %.c $(BUILD_DIR)/funcdefs.h
	@mkdir -p $(dir $@)
	$(CC) -c -include $(BUILD_DIR)/funcdefs.h $(CPPFLAGS) $(CFLAGS) $< -o $@

$(BUILD_DIR)/%.o: %.cpp $(BUILD_DIR)/funcdefs.h
	@mkdir -p $(dir $@)
	$(CXX) -c -include $(BUILD_DIR)/funcdefs.h $(CPPFLAGS) $(CXXFLAGS) $< -o $@

# vim: noet

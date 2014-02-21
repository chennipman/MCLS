#include "../headers/array.h"

/***********************************************************************************/
/*  Function to initialize the coupling between the interface and the flow field.  */
/*            									          */
/*  										          */
/*  Programmer	: Duncan van der Heul       					          */
/*  Date	: 10-03-2013       						          */
/*  Update	:        							          */
/***********************************************************************************/
/* Notes									          */
/* The coupling is through the surface tension forces that appear in the           */
/* momentum equation. In the momentum equation the option exists to include        */
/* the sources in the predictor or in the corrector equation.                      */
/* When the sources are included in the predictor equation the discretisation      */
/* is far less dissipative.                                                        */
/*                                                                                 */
/*                                                                                 */
/***********************************************************************************/
EXPORT void initialize_coupling(
 	    Array3<double> level_set, 				// level-set field
	    Array3<double> pressure,				// pressure field
	    Array3<double> curvature, 				// interface curvature
	    Array3<double> unsmoothed_curvature,			// interface curvature without smoothing
	    Array3<double> surface_tension_body_force_x1,	// x1 component of the body force due to
								// CSF formulation of surface tension model
	    Array3<double> surface_tension_body_force_x2,	// x2 component of the body force due to
								// CSF formulation of surface tension model
	    Array3<double> surface_tension_body_force_x3,	// x3 component of the body force due to
								// CSF formulation of surface tension model
	    Array3<double> momentum_source_term_u_1, 		// complete source term for the momentum equation
								// in x1 direction=(-p,1+ g_1 +F1)
	    Array3<double> momentum_source_term_u_2, 		// complete source term for the momentum equation
								// in x2 direction=(-p,2+ g_2 +F2)
	    Array3<double> momentum_source_term_u_3, 		// complete source term for the momentum equation
								// in x3 direction=(-p,3+ g_3 +F3)
	    Array3<double> scaled_density_u1,			// scaled density for the controlvolumes
								// of the momentum equation in x1 direction
	    Array3<double> scaled_density_u2,			// scaled density for the controlvolumes
								// of the momentum equation in x2 direction
	    Array3<double> scaled_density_u3,			// scaled density for the controlvolumes
								// of the momentum equation in x3 direction
	    int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	    int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	    int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	    double mesh_width_x1,				// grid spacing in x1 direction (uniform)
	    double mesh_width_x2,				// grid spacing in x2 direction (uniform)
	    double mesh_width_x3,				// grid spacing in x3 direction (uniform)
	    double rho_plus_over_rho_minus,			// ratio of the densities of the two phases
	    double actual_time_step_navier_stokes, 		// actual time step for Navier-Stokes solution algorithm 
	    vector gravity,					// gravitational acceleration vector 
	    double sigma_over_rho_minus,			// sigma / rho_minus (scaled sigma)
	    int apply_curvature_smoothing,			// =1, apply curvature smoothing
								// =0, use unsmoothed curvature
	    int number_curvature_smoothing_steps,		// number of iterations applied in the
								// curvature smoothing algorithm
	    int apply_curvature_smoothing_filter,		// =1, apply curvature smoothing filter
								// =0, do not apply curvature smoothing filter
	    double smoothing_distance_factor,
           double lower_bound_derivatives			// lower bound for the first partial derivatives
								// to consider it a limiting case of vanishing
								// partial derivatives
      )
	{
        double maximum_weighted_curvature=0;                    // maximum 'active' value of the curvature
                                                                // used to evaluate the capillary time step
                                                                // restriction

	/* compute the density field at the velocity points */

	compute_density_u_controlvolumes( level_set, scaled_density_u1, scaled_density_u2, scaled_density_u3,		
						     number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
       							rho_plus_over_rho_minus, lower_bound_derivatives);

 
	/* compute surface tension body forces in the CSF model */


	compute_surface_tension_body_force(level_set, 				
						surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,	
					 	 curvature, unsmoothed_curvature,
						  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,				
						    mesh_width_x1, mesh_width_x2, mesh_width_x3,				
						      rho_plus_over_rho_minus, sigma_over_rho_minus,			
							maximum_weighted_curvature,			
							  apply_curvature_smoothing,			
							    number_curvature_smoothing_steps, 
							      apply_curvature_smoothing_filter,
							      smoothing_distance_factor);								
    
	/* compute the complete source terms for the momentum equation */

       /* intialize the pressure with constant value 0 */

       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2,
                         number_primary_cells_k+2, pressure, 0.0);

	compute_momentum_source_terms( level_set, pressure,					
				      surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,		
					momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                                     scaled_density_u1, scaled_density_u2, scaled_density_u3,
					  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,				
					     mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					      rho_plus_over_rho_minus, actual_time_step_navier_stokes, gravity);


  

	}

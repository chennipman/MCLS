#include "../headers/array.h"
/********************************************************************************/
/*  Function to allocate the main variables for the computation                 */
/*  										       */
/*  Programmer	: Duncan van der Heul       					       */
/*  Date	: 10-03-2013       						       */
/*  Update	:        							       */
/********************************************************************************/
/* Notes									       */
/********************************************************************************/
EXPORT void allocate_main_variables(
      Array3<double> & u_1_velocity_new, 			// velocity field at new time level x1 direction
      Array3<double> & u_2_velocity_new, 			// velocity field at new time level x2 direction
      Array3<double> & u_3_velocity_new,				// velocity field at new time level x3 direction
      Array3<double> & u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> & u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> & u_3_velocity_old,				// velocity field at old time level x3 direction
      Array3<double> & pressure,					// pressure
      Array3<double> & level_set_old,				// level-set field at old time level
      Array3<double> & level_set_new,				// level-set field at new time level
      Array3<double> & volume_of_fluid,				// volume of fluid field
      Array3<double> & curvature,					// interface curvature
      Array3<double> & unsmoothed_curvature,			// interface curvature without smoothing
      Array3<double> & momentum_source_term_u_1,              // source term of the momentum equation in x1 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> & momentum_source_term_u_2,              // source term of the momentum equation in x2 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> & momentum_source_term_u_3,              // source term of the momentum equation in x3 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> & surface_tension_body_force_x1,         // source term of the momentum equation in x1 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> & surface_tension_body_force_x2,         // source term of the momentum equation in x2 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> & surface_tension_body_force_x3,         // source term of the momentum equation in x3 direction
                                                          // defined on all u1 points (including boundaries)
      Array3<double> & scaled_density_u1,                     // scaled density for the controlvolumes
                                                          // of the momentum equation in x1 direction
      Array3<double> & scaled_density_u2,                     // scaled density for the controlvolumes
                                                          // of the momentum equation in x2 direction
      Array3<double> & scaled_density_u3,                     // scaled density for the controlvolumes
                                                          // of the momentum equation in x3 direction
      int number_primary_cells_i,				// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,				// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k				// number of primary (pressure) cells in x3 direction
	 )
    {
      
   

      /* allocate all velocity fields */
      
      u_1_velocity_old.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_old.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_old.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      u_1_velocity_new.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_new.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_new.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
     
      /* allocate pressure field and the fields associated with the interface */
      
      pressure.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      volume_of_fluid.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      level_set_old.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      level_set_new.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      curvature.create(number_primary_cells_i+2,number_primary_cells_j+2, number_primary_cells_k+2);
      unsmoothed_curvature.create(number_primary_cells_i+2,number_primary_cells_j+2, number_primary_cells_k+2);


      /* allocate all fields related to the interface to flow coupling */
      /* allocate memory for the coupling terms between flow field and interface */

      scaled_density_u1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      scaled_density_u2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      scaled_density_u3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

      surface_tension_body_force_x1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      surface_tension_body_force_x2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      surface_tension_body_force_x3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

      momentum_source_term_u_1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      momentum_source_term_u_2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      momentum_source_term_u_3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

    }
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

/********************************************************************************/
/*  Function to allocate the main variables for the computation                 */
/*  										       */
/*  Programmer	: Duncan van der Heul       					       */
/*  Date	: 10-03-2013       						       */
/*  Update	:        							       */
/********************************************************************************/
/* Notes									       */
/********************************************************************************/
    void allocate_main_variables(
      double *** & u_1_velocity_new, 			// velocity field at new time level x1 direction
      double *** & u_2_velocity_new, 			// velocity field at new time level x2 direction
      double *** & u_3_velocity_new,				// velocity field at new time level x3 direction
      double *** & u_1_velocity_old, 			// velocity field at old time level x1 direction
      double *** & u_2_velocity_old, 			// velocity field at old time level x2 direction
      double *** & u_3_velocity_old,				// velocity field at old time level x3 direction
      double *** & pressure,					// pressure
      double *** & level_set_old,				// level-set field at old time level
      double *** & level_set_new,				// level-set field at new time level
      double *** & volume_of_fluid,				// volume of fluid field
      double *** & curvature,					// interface curvature
      double *** & unsmoothed_curvature,			// interface curvature without smoothing
      double *** & momentum_source_term_u_1,              // source term of the momentum equation in x1 direction
                                                          // defined on all u1 points (including boundaries)
      double *** & momentum_source_term_u_2,              // source term of the momentum equation in x2 direction
                                                          // defined on all u1 points (including boundaries)
      double *** & momentum_source_term_u_3,              // source term of the momentum equation in x3 direction
                                                          // defined on all u1 points (including boundaries)
      double *** & surface_tension_body_force_x1,         // source term of the momentum equation in x1 direction
                                                          // defined on all u1 points (including boundaries)
      double *** & surface_tension_body_force_x2,         // source term of the momentum equation in x2 direction
                                                          // defined on all u1 points (including boundaries)
      double *** & surface_tension_body_force_x3,         // source term of the momentum equation in x3 direction
                                                          // defined on all u1 points (including boundaries)
      double *** & scaled_density_u1,                     // scaled density for the controlvolumes
                                                          // of the momentum equation in x1 direction
      double *** & scaled_density_u2,                     // scaled density for the controlvolumes
                                                          // of the momentum equation in x2 direction
      double *** & scaled_density_u3,                     // scaled density for the controlvolumes
                                                          // of the momentum equation in x3 direction
      int number_primary_cells_i,				// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,				// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k				// number of primary (pressure) cells in x3 direction
	 )
    {
      
   
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

      /* allocate all velocity fields */
      
      u_1_velocity_old	=	double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_old	=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_old	=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      u_1_velocity_new	=	double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_new	=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_new	=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
     
      /* allocate pressure field and the fields associated with the interface */
      
      pressure			=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      volume_of_fluid		=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      level_set_old		=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      level_set_new		=	double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      curvature		=	double_Matrix2(number_primary_cells_i+2,number_primary_cells_j+2, number_primary_cells_k+2);
      unsmoothed_curvature	=	double_Matrix2(number_primary_cells_i+2,number_primary_cells_j+2, number_primary_cells_k+2);


      /* allocate all fields related to the interface to flow coupling */
      /* allocate memory for the coupling terms between flow field and interface */

      scaled_density_u1=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      scaled_density_u2=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      scaled_density_u3=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

      surface_tension_body_force_x1=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      surface_tension_body_force_x2=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      surface_tension_body_force_x3=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

      momentum_source_term_u_1=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      momentum_source_term_u_2=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      momentum_source_term_u_3=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);

    }
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
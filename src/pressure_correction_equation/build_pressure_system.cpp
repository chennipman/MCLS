class vector
{
public:
  double u1,u2,u3;
  vector(double u1, double u2, double u3);
  vector( void);
};

/********************************************************************************/
/********************************************************************************/
/*  Function to construct the whole system of linear equations for the          */
/*  pressure correction equation equation                                       */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* For the pressure, an inhomogeneous neumann boundary condition hold on the    */
/* boundaries where the velocity is zero. The normal derivative of the pressure */
/* balances with the gravitational acceleration.				*/
/********************************************************************************/

  void build_pressure_system(
      double **pressure_correction_matrix,    		// pressure correction matrix
      double *final_pressure_rhs,	     		// right hand side of pressure correction equation
							// including contributions 
							// inhomogeneous boundary conditions
      double ***level_set_new,				// level-set field
      double ***sources_x1,		        	// source term of the momentum equation in x1 direction
					        	// defined on all u1 points (including boundaries)
      double ***sources_x2,  		        	// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      double ***sources_x3,  		        	// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      double ***surface_tension_body_force_x1,		// source term of the momentum equation in x1 direction
					        	// defined on all u1 points (including boundaries)
      double ***surface_tension_body_force_x2,		// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      double ***surface_tension_body_force_x3,		// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      double ***scaled_density_u1,			// scaled density for the controlvolumes
							// of the momentum equation in x1 direction
      double ***scaled_density_u2,			// scaled density for the controlvolumes
							// of the momentum equation in x2 direction
      double ***scaled_density_u3,			// scaled density for the controlvolumes
							// of the momentum equation in x3 direction
      double ***u_1_velocity_star, 	        	// velocity field at star time level x1 direction
      double ***u_2_velocity_star, 	        	// velocity field at star time level x2 direction
      double ***u_3_velocity_star,	        	// velocity field at star time level x3 direction
      double ***pressure_boundary_condition_x1,    	//inhomogeneous boundary condition for
                                                   	// the pressure planes with normal in x1 direction
      double ***pressure_boundary_condition_x2,    	//inhomogeneous boundary condition for the pressure
                                                   	// the pressure planes with normal in x1 direction
      double ***pressure_boundary_condition_x3,    	//inhomogeneous boundary condition for the pressure
                                                   	// the pressure planes with normal in x1 direction
      double mesh_width_x1,		        	// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		        	// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		        	// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	        	// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	        	// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	        	// number of primary (pressure) cells in x3 direction
      double actual_time_step_navier_stokes,    	// actual time step for Navier-Stokes solution algorithm 
      double rho_plus_over_rho_minus,	       		// ratio of density where (level set >0) and 
					        	// density where (level set < 0)
      int continuous_surface_force_model,       	// =1, the continuous surface force model is applied
					        	// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor,    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation
      vector gravity					// gravitational acceleration vector 

       )

  {
    double ***double_Matrix2(                      	// allocate memory for a three-
        int number_primary_cells_i,                	// dimensional array of doubles
	 int number_primary_cells_j, 			
	 int number_primary_cells_k
        );
    
    void free_double_Matrix2(                      	// deallocate memory for a three
        double ***doubleMatrix2,                   	// dimensional array of doubles
        int number_primary_cells_i,	
	 int number_primary_cells_j
        );

    void build_pressure_matrix (		       // build pressure correction matrix
        double **pressure_correction_matrix, 
        double mesh_width_x1,		     
        double mesh_width_x2,		     
        double mesh_width_x3,		     
        double ***level_set_new, 	
	double ***scaled_density_u1,
	double ***scaled_density_u2,
	double ***scaled_density_u3,
        int number_primary_cells_i,	     
        int number_primary_cells_j,	     
        int number_primary_cells_k,	     
        double rho_plus_over_rho_minus
        );
    
    void build_pressure_rhs_momentum_part_CSF( 				// build pressure correction right hand side
        double ***initial_pressure_rhs,	       				// with the source terms and the divergence
        double ***sources_x1,		       				// of u star
        double ***sources_x2,  		  
        double ***sources_x3,  		  
        double ***surface_tension_body_force_x1,	 	  
        double ***surface_tension_body_force_x2,	          
        double ***surface_tension_body_force_x3,		  
        double ***u_1_velocity_star, 	  
        double ***u_2_velocity_star, 	  
        double ***u_3_velocity_star,	  
        double mesh_width_x1,		  
        double mesh_width_x2,		  
        double mesh_width_x3,		  
        int number_primary_cells_i,	  
        int number_primary_cells_j,	  
        int number_primary_cells_k,	    
        double actual_time_step_navier_stokes, 
        int continuous_surface_force_model,    
        int source_terms_in_momentum_predictor);
    void build_pressure_rhs_boundary(					// apply boundary conditions to the pressure
	double ***initial_pressure_rhs,	     				// correction equation right hand side
	double *final_pressure_rhs,	     
	double mesh_width_x1,		     
	double mesh_width_x2,		     
	double mesh_width_x3,		     
	int number_primary_cells_i,	     
	int number_primary_cells_j,	     
	int number_primary_cells_k,	     
	double ***pressure_boundary_condition_x1, 
	double ***pressure_boundary_condition_x2,	
	double ***pressure_boundary_condition_x3);
   void set_pressure_boundary_condition(				// set inhomogeneous boundary
        vector gravity,				    			// conditions for the pressure
	double ***pressure_boundary_condition_x1,	
	double ***pressure_boundary_condition_x2,
	double ***pressure_boundary_condition_x3,
	int number_primary_cells_i,
	int number_primary_cells_j,
	int number_primary_cells_k);
     
       double ***initial_pressure_rhs;	     	    			// right hand side of pressure correction equation
						    			// excluding contributions 
						    			// inhomogeneous boundary conditions
	
      /* allocate memory for the initial pressure right hand side */
       
      initial_pressure_rhs=double_Matrix2(number_primary_cells_i+2, 
					    number_primary_cells_j+2,
					      number_primary_cells_k+2);
      
      set_pressure_boundary_condition(gravity,
				       pressure_boundary_condition_x1,
					pressure_boundary_condition_x2,
					  pressure_boundary_condition_x3,
					      number_primary_cells_i,
					       number_primary_cells_j,
						number_primary_cells_k);
      
      /* build pressure correction matrix */
      
      build_pressure_matrix(pressure_correction_matrix, 
			      mesh_width_x1, mesh_width_x2, mesh_width_x3,		     
				level_set_new, 		  
				 scaled_density_u1, scaled_density_u2, scaled_density_u3,
				  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	     
				    rho_plus_over_rho_minus);

      /* build pressure correction right hand side, the part that originates from the momentum equation */

      build_pressure_rhs_momentum_part_CSF(initial_pressure_rhs,	       
					    sources_x1, sources_x2, sources_x3,	  		  
					      surface_tension_body_force_x1, 
						surface_tension_body_force_x2, 
						  surface_tension_body_force_x3,		  
						    u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,	  
						      mesh_width_x1, mesh_width_x2, mesh_width_x3,		  
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	    
							  actual_time_step_navier_stokes, 
							    continuous_surface_force_model,    
							      source_terms_in_momentum_predictor);

      /* apply the boundary conditions to the pressure correction equation */

      build_pressure_rhs_boundary(initial_pressure_rhs,	final_pressure_rhs,	
				    mesh_width_x1, mesh_width_x2, mesh_width_x3,		
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
					pressure_boundary_condition_x1, 
					  pressure_boundary_condition_x2,	
					    pressure_boundary_condition_x3);
      
      
      
      /* deallocate the memory for the initial pressure right hand side 	*/
      
      free_double_Matrix2(initial_pressure_rhs,number_primary_cells_i+2, 
					    number_primary_cells_j+2);

  }
      
      
 
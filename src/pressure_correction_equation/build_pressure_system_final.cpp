#include "../headers/array.h"

/********************************************************************************/
/********************************************************************************/
/*  Function to construct the whole system of linear equations for the          */
/*  pressure correction equation equation                                       */
/*  										*/
/*  Programmer	: Coen Hennipman						*/
/*  Date	: 20-05-2014       						*/
/*  Update	:         							*/
/********************************************************************************/
/* Notes									*/
/* For the pressure, an inhomogeneous neumann boundary condition hold on the    */
/* boundaries where the velocity is zero. The normal derivative of the pressure */
/* balances with the gravitational acceleration.				*/
/* This function is almost similar to build_pressure_system,			*/
/* the difference is the build of the rhs. 					*/
/********************************************************************************/

EXPORT void build_pressure_system_final(
      Array2<double> pressure_correction_matrix,    	// pressure correction matrix
      Array1<double> final_pressure_rhs,	     	// right hand side of pressure correction equation
							// including contributions 
							// inhomogeneous boundary conditions
      Array3<double> level_set_new,			// level-set field
      Array3<double> scaled_density_u1,			// scaled density for the controlvolumes
							// of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,			// scaled density for the controlvolumes
							// of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,			// scaled density for the controlvolumes
							// of the momentum equation in x3 direction
      Array3<double> u_1_new_con_diff, 			// contains the convection and diffusion terms 
      Array3<double> u_2_new_con_diff, 			// contains the convection and diffusion terms
      Array3<double> u_3_new_con_diff,			// contains the convection and diffusion terms 
      Array3<double> pressure_boundary_condition_x1,	//inhomogeneous boundary condition for
                                                   	// the pressure planes with normal in x1 direction
      Array3<double> pressure_boundary_condition_x2, 	//inhomogeneous boundary condition for the pressure
                                                   	// the pressure planes with normal in x1 direction
      Array3<double> pressure_boundary_condition_x3,   	//inhomogeneous boundary condition for the pressure
                                                   	// the pressure planes with normal in x1 direction
      double mesh_width_x1,		        	// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		        	// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		        	// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	        	// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	        	// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	        	// number of primary (pressure) cells in x3 direction
      double rho_plus_over_rho_minus,	       		// ratio of density where (level set >0) and 
					        	// density where (level set < 0)
      vector gravity					// gravitational acceleration vector 

       )

  {
       Array3<double> initial_pressure_rhs;	     	    		// right hand side of pressure correction equation
						    			// excluding contributions 
						    			// inhomogeneous boundary conditions
	
      /* allocate memory for the initial pressure right hand side */
       
      initial_pressure_rhs.create(number_primary_cells_i+2, 
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

      /* build pressure correction right hand side, 			*/
      /* compute the divergence of the convection and diffusion terms	*/
      /* the divergence of the convection and diffusion terms is the 	*/
      /* inital pressure in the rhs 					*/
      
      divergence_of_vector_field(initial_pressure_rhs, 
      				u_1_new_con_diff, u_2_new_con_diff, u_3_new_con_diff, 
      				mesh_width_x1, mesh_width_x2, mesh_width_x3,
				number_primary_cells_i,number_primary_cells_j,number_primary_cells_k);


      /* apply the boundary conditions to the pressure correction equation */

      build_pressure_rhs_boundary(initial_pressure_rhs,	final_pressure_rhs,	
				    mesh_width_x1, mesh_width_x2, mesh_width_x3,		
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
					pressure_boundary_condition_x1, 
					  pressure_boundary_condition_x2,	
					    pressure_boundary_condition_x3);
      
      
      
      /* deallocate the memory for the initial pressure right hand side 	*/
      
      initial_pressure_rhs.destroy();

  }
      
      
 

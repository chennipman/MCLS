#include<cstdlib>
#include<iostream>

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the advection of the level_set function          	*/
/*  with a 3rd order weno scheme                                                */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The level-set field is advected using a 3rd order weno scheme in space       */
/* and an explicit TVD Runge-Kutta time integration method.                     */
/* We can write the level-set equation as                               	*/
/*      d phi                                                                   */
/*      ------=  Operator(phi)                                                  */
/* 	d  t                                                                    */
/*   where Operator(phi) is the nonlinear discretisation of the convective      */
/* terms of the level-set equation                                              */
/*  A third order explicit TVD Runge-Kutta method can then be defined as:       */
/*   phi_stage_1= phi_old + dt* Operator(phi_old)                               */
/*   phi_stage_2= 3/4 phi_old + 1/4 phi_stage_1 + 1/4 dt Operator(phi_stage_1)  */
/*   phi_new    = 1/3 phi_old + 2/3 phi_stage_1 + 2/3 dt Operator(phi_stage_2)  */
/********************************************************************************/
//
void advect_level_set_higher_order(
      double ***level_set_old, 		        // level set field at old time level
      double ***level_set_star, 		// level set field at star time level
      double ***u_1_velocity_new, 	        // velocity field at new time level x1 direction
      double ***u_2_velocity_new, 	        // velocity field at new time level x2 direction
      double ***u_3_velocity_new,	        // velocity field at new time level x3 direction
      int number_primary_cells_i,	        // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	        // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	        // number of primary (pressure) cells in x3 direction
      double actual_time_step_level_set,        // time step used for level-set advection
					        // computed from all stability restrictions and 
					        // possibly subscycling
      double mesh_width_x1,		        // grid spacing in x1 direction (uniform)
      double mesh_width_x2,		        // grid spacing in x2 direction (uniform)
      double mesh_width_x3		        // grid spacing in x3 direction (uniform)
      
)
{
/*******************************************************************************************/      
/* 				function definitions 					   */
/*******************************************************************************************/      
   double ***double_Matrix2(                     // allocate memory for a three-
          int number_primary_cells_i,            // dimensional array of doubles
          int number_primary_cells_j,             
          int number_primary_cells_k
          );
   void free_double_Matrix2(                     // deallocate memory for a three
          double ***doubleMatrix2,               // dimensional array of doubles
          int number_primary_cells_i,     
          int number_primary_cells_j
          );
   void field_neumann_boundary(			// apply neumann boundary condition to
    	double ***field, 			// cell centered field
    	int number_primary_cells_i,	
    	int number_primary_cells_j,	
    	int number_primary_cells_k	
	  );
   void  field_extrapolate_boundary(      	// extrapolate field to virtual cells
        double ***field, 			
        int number_primary_cells_i,	
        int number_primary_cells_j,	
        int number_primary_cells_k	
	);
  double compute_level_set_flux( 		// function to interpolate the 
    	double cell_face_velocity,		// level-set to the cell faces
    	double left_phi, 			// and evaluate the level-set flux
    	double right_phi
	  );
  void evaluate_convection_operator_weno(       // evaluate the convection operator
      double ***level_set,                      // discretised with weno scheme
      double ***convection_operator,            
      double ***u_1_velocity_new,               
      double ***u_2_velocity_new,               
      double ***u_3_velocity_new,               
      int number_primary_cells_i,               
      int number_primary_cells_j,               
      int number_primary_cells_k,               
      double mesh_width_x1,                     
      double mesh_width_x2,                     
      double mesh_width_x3                      
       );
  void compute_first_stage_RK(                  // compute the first stage in RK-method
      double ***level_set,                        
      double ***level_set_stage_1,                
      double ***convection_operator,              
      int number_primary_cells_i,                 
      int number_primary_cells_j,                 
      int number_primary_cells_k,                 
      double actual_time_step_level_set           
      );
  void compute_second_stage_RK(                 // compute the second stage in RK-method
      double ***level_set,                        
      double ***level_set_stage_1,                
      double ***level_set_stage_2,                
      double ***convection_operator,              
      int number_primary_cells_i,                 
      int number_primary_cells_j,                 
      int number_primary_cells_k,                 
      double actual_time_step_level_set           
      );
  void compute_new_time_RK(                     // compute the new time level in RK-method
      double ***level_set_star,
      double ***level_set_old,
      double ***level_set_stage_1,                
      double ***level_set_stage_2,                
      double ***convection_operator,              
      int number_primary_cells_i,                 
      int number_primary_cells_j,                 
      int number_primary_cells_k,                 
      double actual_time_step_level_set          
       );
/*******************************************************************************************/      
/* 				                     					   */
/*******************************************************************************************/      

    double one_over_dx1	=    		        // 1/(grid spacing in x1 direction)
	1.0/(mesh_width_x1);
    double one_over_dx2	=    		        // 1/(grid spacing in x2 direction)
	1.0/(mesh_width_x1);
    double one_over_dx3	=    		        // 1/(grid spacing in x3 direction)
	1.0/(mesh_width_x3);
    double cell_face_velocity; 		        // velocity at the cell face (no interpolation required)
					        // because of staggered mesh
    double flux_level_set; 		        // advective flux
    int i_index, j_index, k_index;              // local variables for loop indexing
    
    
    double ***level_set_stage_1;                // stage 1 for RK time integration
    double ***level_set_stage_2;                // stage 2 for RK time integration
    double ***convection_operator;              // right hand side of the system of ode's
    
    
    /* allocate memory for the stages and the operator */
    
      level_set_stage_1=double_Matrix2(number_primary_cells_i+2, 
                                                    number_primary_cells_j+2,
                                                      number_primary_cells_k+2);
      level_set_stage_2=double_Matrix2(number_primary_cells_i+2, 
                                                    number_primary_cells_j+2,
                                                      number_primary_cells_k+2);
      convection_operator=double_Matrix2(number_primary_cells_i+2, 
                                                    number_primary_cells_j+2,
                                                      number_primary_cells_k+2);

    /* evaluate the operator for the first stage, using the solution of the previous time */
    
      evaluate_convection_operator_weno( level_set_old, convection_operator,            
                                         u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,               
                                           number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,               
                                             mesh_width_x1, mesh_width_x2, mesh_width_x3);
   
    /* compute the first stage */
    
      compute_first_stage_RK(level_set_old, level_set_stage_1, convection_operator,              
                                number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                 
                                   actual_time_step_level_set);

      field_extrapolate_boundary(level_set_stage_1, number_primary_cells_i, 
                           number_primary_cells_j,number_primary_cells_k);      
      
    /* evaluate the operator for the second stage, using the first stage solution */

      evaluate_convection_operator_weno( level_set_stage_1, convection_operator,            
                                         u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,               
                                           number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,               
                                             mesh_width_x1, mesh_width_x2, mesh_width_x3);
    
    /* compute the second stage */

      compute_second_stage_RK(level_set_old, level_set_stage_1, level_set_stage_2, convection_operator,              
                                number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                 
                                  actual_time_step_level_set);

      field_extrapolate_boundary(level_set_stage_2, number_primary_cells_i, 
                           number_primary_cells_j,number_primary_cells_k);      
    
 
    /* evaluate the operator for the third stage, using the second stage solution */

      evaluate_convection_operator_weno( level_set_stage_2, convection_operator,            
                                         u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,               
                                           number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,               
                                             mesh_width_x1, mesh_width_x2, mesh_width_x3);
    
    /* compute the solution at the new time */
    
      compute_new_time_RK(level_set_star, level_set_old,
                            level_set_stage_1, level_set_stage_2, convection_operator,              
                             number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                 
                                actual_time_step_level_set);

    /* extend the level-set field to the virtual cells */
    
      field_extrapolate_boundary(level_set_star, number_primary_cells_i, 
			   number_primary_cells_j,number_primary_cells_k);	
     
    /* deallocate temporary storage */
          
      free_double_Matrix2(convection_operator, number_primary_cells_i+2, 
                                  number_primary_cells_j+2);
      free_double_Matrix2(level_set_stage_1, number_primary_cells_i+2,
                                  number_primary_cells_j+2);
      free_double_Matrix2(level_set_stage_2, number_primary_cells_i+2,
                                  number_primary_cells_j+2);
   
}




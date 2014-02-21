#include "../headers/array.h"
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
EXPORT void advect_level_set_higher_order(
      Array3<double> level_set_old, 		// level set field at old time level
      Array3<double> level_set_star, 		// level set field at star time level
      Array3<double> u_1_velocity_new, 	        // velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new, 	        // velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,	        // velocity field at new time level x3 direction
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
    
    
    Array3<double> level_set_stage_1;                // stage 1 for RK time integration
    Array3<double> level_set_stage_2;                // stage 2 for RK time integration
    Array3<double> convection_operator;              // right hand side of the system of ode's
    
    
    /* allocate memory for the stages and the operator */
    
      level_set_stage_1.create(number_primary_cells_i+2, 
                                                    number_primary_cells_j+2,
                                                      number_primary_cells_k+2);
      level_set_stage_2.create(number_primary_cells_i+2, 
                                                    number_primary_cells_j+2,
                                                      number_primary_cells_k+2);
      convection_operator.create(number_primary_cells_i+2, 
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
          
      convection_operator.destroy();
      level_set_stage_1.destroy();
      level_set_stage_2.destroy();
   
}




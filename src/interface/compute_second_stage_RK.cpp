 
 
 
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the second stage of three stage Runge-Kutta method      */
/*  for advection equation of level set.                                        */
/*                                                                              */
/*  Programmer  : Duncan van der Heul                                           */
/*  Date        : 10-03-2013                                                    */
/*  Update      :                                                               */
/********************************************************************************/
/* Notes                                                                        */
/* The level-set field is advected using a 3rd order weno scheme in space       */
/* and an explicit TVD Runge-Kutta time integration method.                     */
/* We can write the level-set equation as                                       */
/*      d phi                                                                   */
/*      ------=  Operator(phi)                                                  */
/*      d  t                                                                    */
/*   where Operator(phi) is the nonlinear discretisation of the convective      */
/* terms of the level-set equation                                              */
/*  A third order explicit TVD Runge-Kutta method can then be defined as:       */
/*   phi_stage_1= phi_old + dt* Operator(phi_old)                               */
/*   phi_stage_2= 3/4 phi_old + 1/4 phi_stage_1 + 1/4 dt Operator(phi_stage_1)  */
/*   phi_new    = 1/3 phi_old + 2/3 phi_stage_1 + 2/3 dt Operator(phi_stage_2)  */
/********************************************************************************/
//
 void compute_second_stage_RK(
    double ***level_set,                        // level set at previous time level
    double ***level_set_stage_1,                // stage 1 for RK time integration
    double ***level_set_stage_2,                // stage 2 for RK time integration
    double ***convection_operator,              // right hand side of the system of ode's
    int number_primary_cells_i,                 // number of primary (pressure) cells in x1 direction
    int number_primary_cells_j,                 // number of primary (pressure) cells in x2 direction
    int number_primary_cells_k,                 // number of primary (pressure) cells in x3 direction
    double actual_time_step_level_set           // time step used for level-set advection
                                                // computed from all stability restrictions and 
                                                // possibly subscycling
    )
 { 
    int i_index, j_index, k_index;              // local variables for loop indexing
 
      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
            for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
            {
                   for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
                   {
     
                        level_set_stage_2[i_index  ][j_index][k_index]=
                                               0.75*level_set[i_index  ][j_index][k_index]+
                                                 0.25*level_set_stage_1[i_index  ][j_index][k_index]+
                                                  0.25*actual_time_step_level_set*
                                                        convection_operator[i_index  ][j_index][k_index]; 
  
                   }  
  
            }  
     
      } 
 
 
 }   
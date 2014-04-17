#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the redistribution velocity field                       */
/*  method.                                                                     */
/*                                                                              */
/*  Programmer  : Duncan van der Heul                                           */
/*  Date        : 10-03-2013                                                    */
/*  Update      :                                                               */
/********************************************************************************/
/* Notes                                                                        */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/
EXPORT void compute_redistribution_velocity_field(
        Array3<double> level_set,                               // level set field 
                                                                // mass conserving
        Array3<double> redistribution_velocity_x1,              // artificial redistribution velocity x1 direction
        Array3<double> redistribution_velocity_x2,              // artificial redistribution velocity x2 direction
        Array3<double> redistribution_velocity_x3,              // artificial redistribution velocity x3 direction
        double mesh_width_x1,                                   // grid spacing in x1 direction (uniform)
        double mesh_width_x2,                                   // grid spacing in x2 direction (uniform)
        double mesh_width_x3,                                   // grid spacing in x3 direction (uniform)
        int number_primary_cells_i,                             // number of primary (pressure) cells in x1 direction
        int number_primary_cells_j,                             // number of primary (pressure) cells in x2 direction
        int number_primary_cells_k                              // number of primary (pressure) cells in x3 direction
)
{
        int i_index, j_index, k_index;                          // local variables for loop indexing
      
         set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2, 
                            number_primary_cells_k+2, redistribution_velocity_x1, 0.0);
         set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
                            number_primary_cells_k+2, redistribution_velocity_x2, 0.0);
         set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
                            number_primary_cells_k+1, redistribution_velocity_x3, 0.0);
        
        /* first compute the error distribution velocities */
        /* these remain fixed during the error distribution */
        /* sweeps, because they depend on the level-set field */
        /* and that does not change during the whole process */

        /* compute advection velocity in 1-direction */
        
        for( i_index=0;i_index<number_primary_cells_i+1;i_index++){
            for(j_index=0;j_index<number_primary_cells_j+2;j_index++){
                for(k_index=0;k_index<number_primary_cells_k+2;k_index++){
                  
                  redistribution_velocity_x1[i_index][j_index][k_index]=
                                              compute_redistribution_velocity(level_set[i_index][j_index][k_index],
                                                level_set[i_index+1][j_index][k_index],
                                                  mesh_width_x1);
                }
            }
        }
        
        /* compute advection velocity in 2-direction */
        
        for( i_index=0;i_index<number_primary_cells_i+2;i_index++){
            for(j_index=0;j_index<number_primary_cells_j+1;j_index++){
                for(k_index=0;k_index<number_primary_cells_k+2;k_index++){
                  
                  redistribution_velocity_x2[i_index][j_index][k_index]=
                                              compute_redistribution_velocity(level_set[i_index][j_index][k_index],
                                                level_set[i_index][j_index+1][k_index],
                                                  mesh_width_x2);
                }
            }
        }
        
        
        /* compute advection velocity in 3-direction */
        
        
        for( i_index=0;i_index<number_primary_cells_i+2;i_index++){
            for(j_index=0;j_index<number_primary_cells_j+2;j_index++){
                for(k_index=0;k_index<number_primary_cells_k+1;k_index++){
                  
                  redistribution_velocity_x3[i_index][j_index][k_index]=
                                              compute_redistribution_velocity(level_set[i_index][j_index][k_index],
                                                level_set[i_index][j_index][k_index+1],
                                                  mesh_width_x3);
                }
            }
        }
        
//         output_predictor_velocityfields(redistribution_velocity_x1, redistribution_velocity_x2, redistribution_velocity_x3,               
//                                         number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,              
//                                           mesh_width_x1, mesh_width_x2, mesh_width_x3);
}
#include "../headers/array.h"
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the advection operator working on a cell                */
/*  centered field, using a third/fifth order weno scheme.                      */
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
/* this function evaluates the convection operator.                             */
/********************************************************************************/
//
EXPORT void evaluate_convection_operator_weno(
      Array3<double> level_set,                      // level set field
      Array3<double> convection_operator,            // the convection operator
      Array3<double> u_1_velocity_new,               // velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new,               // velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,               // velocity field at new time level x3 direction
      int number_primary_cells_i,               // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,               // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,               // number of primary (pressure) cells in x3 direction
      double mesh_width_x1,                     // grid spacing in x1 direction (uniform)
      double mesh_width_x2,                     // grid spacing in x2 direction (uniform)
      double mesh_width_x3                      // grid spacing in x3 direction (uniform)
       )
{
      double cell_face_velocity;                // velocity at the cell face (no interpolation required)
                                                // because of staggered mesh
      double flux_level_set;                    // advective flux
      double one_over_dx1 =                     // 1/(grid spacing in x1 direction)
        1.0/(mesh_width_x1);
      double one_over_dx2 =                     // 1/(grid spacing in x2 direction)
        1.0/(mesh_width_x1);
      double one_over_dx3 =                     // 1/(grid spacing in x3 direction)
        1.0/(mesh_width_x3);
      int i_index, j_index, k_index;            // local variables for loop indexing
    
    /* intialize the convection_operator with constant value 0 */
    
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
                            number_primary_cells_k+2, convection_operator, 0.0);

   /* handle all the convective transport for internal faces in x1 direction */
   /* Note: first real cell has index 1, last real cell index number_primary_cells_n */

      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
            for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
            {
                   for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
                   {
     
                        cell_face_velocity = 0.5*(u_1_velocity_new[i_index][j_index][k_index]+
                                                      u_1_velocity_new[i_index-1][j_index][k_index]);
                        if(cell_face_velocity>=0)
                        {
//                              if(i_index>2 && i_index<number_primary_cells_i-1)
                             if(i_index>2000 && i_index<number_primary_cells_i-1)
                            {
                                  flux_level_set =  weno_flux_computation( 
                                                level_set[i_index-3][j_index][k_index],           
                                                level_set[i_index-2][j_index][k_index],             
                                                level_set[i_index-1][j_index][k_index],             
                                                level_set[i_index][j_index][k_index],                   
                                                level_set[i_index+1][j_index][k_index],              
                                                level_set[i_index+2][j_index][k_index],              
                                                mesh_width_x1);
                            }
                            else
                            {
                                   
                                   flux_level_set =  (level_set[i_index][j_index  ][k_index]-
                                                        level_set[i_index-1][j_index][k_index])*
                                                                one_over_dx1;
                                  
                            }
                        }
                        else
                        {
//                             if(i_index>2 && i_index<number_primary_cells_i-1)
                            if(i_index>2000 && i_index<number_primary_cells_i-1)
                            {
                                  flux_level_set =  weno_flux_computation( 
                                                -level_set[i_index+3][j_index][k_index],           
                                                -level_set[i_index+2][j_index][k_index],             
                                                -level_set[i_index+1][j_index][k_index],             
                                                -level_set[i_index][j_index][k_index],                   
                                                -level_set[i_index-1][j_index][k_index],              
                                                -level_set[i_index-2][j_index][k_index],              
                                                mesh_width_x1) ;
                            }
                            else
                            {
                                   flux_level_set =  (level_set[i_index+1][j_index  ][k_index]-
                                                        level_set[i_index][j_index][k_index])*
                                                                one_over_dx1;
                            }
                        }
                        flux_level_set *= cell_face_velocity;
                        convection_operator[i_index  ][j_index][k_index]-=  flux_level_set; 
  
                   }  
  
            }  
     
      } 
  
/* handle all the convective transport for internal faces in x2 direction */
/* Note: first real cell has index 1, last real cell index number_primary_cells_n */

      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
            for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
            {
                   for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
                   {
      
                        cell_face_velocity = 0.5*(u_2_velocity_new[i_index][j_index][k_index]+
                                                       u_2_velocity_new[i_index][j_index-1][k_index]);
                        if(cell_face_velocity>=0)
                        {
//                             if(j_index>2 && j_index<number_primary_cells_j-1)
                            if(j_index>2000 && j_index<number_primary_cells_j-1)
                            {
                                  flux_level_set =  weno_flux_computation( 
                                                level_set[i_index][j_index-3][k_index],           
                                                level_set[i_index][j_index-2][k_index],             
                                                level_set[i_index][j_index-1][k_index],             
                                                level_set[i_index][j_index  ][k_index],                   
                                                level_set[i_index][j_index+1][k_index],              
                                                level_set[i_index][j_index+2][k_index],              
                                                mesh_width_x2);
                            }
                            else
                            {
                                   
                                   flux_level_set =  (level_set[i_index][j_index  ][k_index]-
                                                        level_set[i_index][j_index-1][k_index])*
                                                                one_over_dx2;
                                  
                            }
                        }
                        else
                        {
//                             if(j_index>2 && j_index<number_primary_cells_j-1)
                            if(j_index>2000 && j_index<number_primary_cells_j-1)
                            {
                                  flux_level_set =   weno_flux_computation( 
                                                -level_set[i_index][j_index+3][k_index],           
                                                -level_set[i_index][j_index+2][k_index],             
                                                -level_set[i_index][j_index+1][k_index],             
                                                -level_set[i_index][j_index  ][k_index],                   
                                                -level_set[i_index][j_index-1][k_index],              
                                                -level_set[i_index][j_index-2][k_index],              
                                                mesh_width_x2);
                             }
                            else
                            {
                                   
                                   flux_level_set =  (level_set[i_index][j_index+1][k_index]-
                                                        level_set[i_index][j_index][k_index])*
                                                                one_over_dx2;
                                  
                            }
                       }
                        flux_level_set *= cell_face_velocity;
                        convection_operator[i_index][j_index][k_index]  -=  flux_level_set;
 
                    }  
   
             }  
    
      }  

/* handle all the convective transport for internal faces in x3 direction */
/* Note: first real cell has index 1, last real cell index number_primary_cells_n */

      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
             for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
             {
                   for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
                   {
      
                        cell_face_velocity = 0.5*(u_3_velocity_new[i_index][j_index][k_index]+
                                                     u_3_velocity_new[i_index][j_index][k_index-1]);
                        if(cell_face_velocity>=0)
                        {
                            if(k_index>2 && k_index<number_primary_cells_k-1)
//                             if(k_index>2000 && k_index<number_primary_cells_k-1)
                            {
                                  flux_level_set =  weno_flux_computation( 
                                                level_set[i_index][j_index][k_index-3],           
                                                level_set[i_index][j_index][k_index-2],             
                                                level_set[i_index][j_index][k_index-1],             
                                                level_set[i_index][j_index][k_index  ],                   
                                                level_set[i_index][j_index][k_index+1],              
                                                level_set[i_index][j_index][k_index+2],              
                                                mesh_width_x3);
                            }
                            else
                            {
                                    flux_level_set =  (level_set[i_index][j_index  ][k_index]-
                                                        level_set[i_index][j_index][k_index-1])*
                                                                one_over_dx3;
                           }
                        }
                        else
                        {
                            if(k_index>2 && k_index<number_primary_cells_k-1)
//                             if(k_index>2000 && k_index<number_primary_cells_k-1)
                            {
                                  flux_level_set =    weno_flux_computation( 
                                                -level_set[i_index][j_index][k_index+3],           
                                                -level_set[i_index][j_index][k_index+2],             
                                                -level_set[i_index][j_index][k_index+1],             
                                                -level_set[i_index][j_index][k_index  ],                   
                                                -level_set[i_index][j_index][k_index-1],              
                                                -level_set[i_index][j_index][k_index-2],              
                                                mesh_width_x3);
                            }
                            else
                            {
                                  flux_level_set =  (level_set[i_index][j_index][k_index+1]-
                                                        level_set[i_index][j_index][k_index])*
                                                                one_over_dx3;
                            }
                        }
                        flux_level_set *= cell_face_velocity;
                        convection_operator[i_index][j_index][k_index]  -=  flux_level_set;

  
                   }  
   
             }  
      }  
}   

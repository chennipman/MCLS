#include "../headers/array.h"
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the source terms for the momentum equation              */
/*                                                                              */
/*  Programmer: Duncan van der Heul                                             */
/*  Date      : 10-03-2013                                                      */
/*  Update    :                                                                 */
/********************************************************************************/
/* Notes                                                                        */
/* In the original implementation of Sander the source terms of the momentum    */
/* equation were computed in the function that applied the pressure correction  */
/* to the velocity field. Now it is contained in this separate function.        */
/********************************************************************************/
EXPORT void compute_momentum_source_terms(
       Array3<double> level_set,                             // level-set field
       Array3<double> pressure,                              // pressure field
       Array3<double> u_1_velocity_star,                     // velocity field at star time level x1 direction
       Array3<double> u_2_velocity_star,                     // velocity field at star time level x2 direction
       Array3<double> u_3_velocity_star,                     // velocity field at star time level x3 direction
       Array3<double> surface_tension_body_force_x1,         // x1 component of the body force due to
                                                        // CSF formulation of surface tension model
       Array3<double> surface_tension_body_force_x2,         // x2 component of the body force due to
                                                        // CSF formulation of surface tension model
       Array3<double> surface_tension_body_force_x3,         // x3 component of the body force due to
                                                        // CSF formulation of surface tension model
       Array3<double> momentum_source_term_u_1,              // complete source term for the momentum equation
                                                        // in x1 direction=(-p,1+ g_1 +F1)
       Array3<double> momentum_source_term_u_2,              // complete source term for the momentum equation
                                                        // in x2 direction=(-p,2+ g_2 +F2)
       Array3<double> momentum_source_term_u_3,              // complete source term for the momentum equation
                                                        // in x3 direction=(-p,3+ g_3 +F3)
       int number_primary_cells_i,                      // number of primary (pressure) cells in x1 direction
       int number_primary_cells_j,                      // number of primary (pressure) cells in x2 direction
       int number_primary_cells_k,                      // number of primary (pressure) cells in x3 direction
       double mesh_width_x1,                            // grid spacing in x1 direction (uniform)
       double mesh_width_x2,                            // grid spacing in x2 direction (uniform)
       double mesh_width_x3,                            // grid spacing in x3 direction (uniform)
       double rho_plus_over_rho_minus,                  // ratio of the densities of the two phases
       double actual_time_step_navier_stokes,           // actual time step for Navier-Stokes solution algorithm
       vector gravity                                   // gravitational acceleration vector
              )
       {
       double density_cell_face_x1;                     // density at cell face, normal in x1 direction
                                                        // (location of u1 velocity)
       double density_cell_face_x2;                     // density at cell face, normal in x2 direction
                                                        // (location of u2 velocity)
       double density_cell_face_x3;                     // density at cell face, normal in x3 direction
                                                        // (location of u3 velocity)
       double one_over_dx1=                             // 1/(grid spacing in x1 direction)
              1.0/(mesh_width_x1);
       double one_over_dx2=                             // 1/(grid spacing in x2 direction)
              1.0/(mesh_width_x2);
       double one_over_dx3=                             // 1/(grid spacing in x3 direction)^2
              1.0/(mesh_width_x3);
       double level_set_left;                           // level-set value in left hand neighbouring cell
       double level_set_right;                          // level-set value in right hand neighbouring cell
       double vector_length_momentum_source_u1;         // number unknowns in array momentum_source_term_u_1
       double vector_length_momentum_source_u2;         // number unknowns in array momentum_source_term_u_2
       double vector_length_momentum_source_u3;         // number unknowns in array momentum_source_term_u_3
       int i_index, j_index, k_index;                   // local variables for loop indexing

       vector_length_momentum_source_u1=
       ( number_primary_cells_i+1 )*number_primary_cells_j * number_primary_cells_k;
       vector_length_momentum_source_u2=
              number_primary_cells_i * ( number_primary_cells_j+1) * number_primary_cells_k;
       vector_length_momentum_source_u3=
              number_primary_cells_i *   number_primary_cells_j * (number_primary_cells_k+1);
       set_constant_vector(vector_length_momentum_source_u1, **momentum_source_term_u_1, 0.0);
       set_constant_vector(vector_length_momentum_source_u2, **momentum_source_term_u_2, 0.0);
       set_constant_vector(vector_length_momentum_source_u3, **momentum_source_term_u_3, 0.0);

       /* compute the source term for the momentum equation, x1 component */

       for(i_index=1;i_index<number_primary_cells_i;i_index++)
       {
              for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
              {
              for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
              {
                     level_set_left =level_set[i_index  ][j_index][k_index];
                     level_set_right=level_set[i_index+1][j_index][k_index];
                     density_cell_face_x1=compute_scaled_density(
                     level_set_left, level_set_right, rho_plus_over_rho_minus);
                     momentum_source_term_u_1[i_index][j_index][k_index]=
                            -1.0* one_over_dx1*(pressure[i_index+1][j_index][k_index]-
                            pressure[i_index][j_index][k_index])/density_cell_face_x1+
                            gravity.u1 + surface_tension_body_force_x1[i_index][j_index][k_index];
              }
              }
       }

       /* compute the source term for the momentum equation, x1 component */

       for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
       {
              for(j_index=1;j_index<number_primary_cells_j;j_index++)
              {
              for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
              {
                     level_set_left =level_set[i_index][j_index  ][k_index];
                     level_set_right=level_set[i_index][j_index+1][k_index];
                     density_cell_face_x2=compute_scaled_density(
                     level_set_left, level_set_right, rho_plus_over_rho_minus);
                     momentum_source_term_u_2[i_index][j_index][k_index]=
                            -1.0* one_over_dx2*(pressure[i_index][j_index+1][k_index]-
                            pressure[i_index][j_index][k_index])/density_cell_face_x2+
                            gravity.u2 + surface_tension_body_force_x2[i_index][j_index][k_index];
              }
              }
       }

       /* compute the source term for the momentum equation, x3 component */

       for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
       {
              for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
              {
              for(k_index=1;k_index<number_primary_cells_k;k_index++)
              {
                     level_set_left =level_set[i_index][j_index][k_index  ];
                     level_set_right=level_set[i_index][j_index][k_index+1];
                     density_cell_face_x3=compute_scaled_density(
                     level_set_left, level_set_right, rho_plus_over_rho_minus);
                     momentum_source_term_u_3[i_index][j_index][k_index]=
                            -1.0* one_over_dx3*(pressure[i_index][j_index][k_index+1]-
                            pressure[i_index][j_index][k_index])/density_cell_face_x3+
                            gravity.u3 + surface_tension_body_force_x3[i_index][j_index][k_index];
              }
              }
       }

       }

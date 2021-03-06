#include "../../src/headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/********************************************************************************/
/*  Function to initialize the flow field                                       */
/*  										       */
/*  Programmer	: Duncan van der Heul       					       */
/*  Date	: 10-03-2013       						       */
/*  Update	:        							       */
/********************************************************************************/
/* Notes									       */
/* For the moment the whole solution is initially set to zero, this should be   */
/* extended to more advanced cases, e.g. a parabolic profile etc.		       */
/********************************************************************************/
EXPORT void initialize_flow_field(
      Array3<double> u_1_velocity_new, 			// velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new, 			// velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,			// velocity field at new time level x3 direction
      Array3<double> u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction
      Array3<double> pressure,				// pressure
      Array3<double> level_set,                         // level-set field
      Array3<double> momentum_source_term_u_1,          // source term of the momentum equation in x1 direction
                                                   // defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2,          // source term of the momentum equation in x2 direction
                                                   // defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3,          // source term of the momentum equation in x3 direction
                                                   // defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x1,     // source term of the momentum equation in x1 direction
                                                   // defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x2,     // source term of the momentum equation in x2 direction
                                                   // defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x3,     // source term of the momentum equation in x3 direction
                                                   // defined on all u1 points (including boundaries)
      Array3<double> scaled_density_u1,                 // scaled density for the controlvolumes
                                                   // of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,                 // scaled density for the controlvolumes
                                                   // of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,                 // scaled density for the controlvolumes
                                                   // of the momentum equation in x3 direction
      boundary_face boundary_faces[6],		// array with all the information
							// for the boundary conditions 
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      vector initial_velocity	,			// initial velocity field at t=0
      vector gravity,                              // gravitational acceleration vector       double tolerance,                            // the tolerance with which the system is solved
      double tolerance,                            // the tolerance with which the system is solved
      double actual_time_step_navier_stokes,       // actual time step for Navier-Stokes solution algorithm
      double rho_plus_over_rho_minus,              // ratio of density where (level set >0) and
                                                   // density where (level set < 0)
      double sigma_over_rho_minus,                 // sigma / rho_minus (scaled sigma)
      int continuous_surface_force_model,          // =1, the continuous surface force model is applied
                                                   // =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor,      // =1, the source terms are applied in the momentum predictor
                                                   // equation
                                                   // =0, the source terms are applied in the momentum corrector
                                                   // equation
      int maximum_iterations_allowed               // maximum number of iterations allowed for the
                                                   // conjugate gradient method
	)
   /* function definitions */
    {
    /* initialize the velocity field */
    /* if no restart file is available, simply set solution to constant value */
    /* as specified in the parameter listing */
    
      /* old time level */
      
        double x,y,z;
        int i,j,k;

  for(i=0;i<number_primary_cells_i+2;i++)
  {
     for(j=0;j<number_primary_cells_j+2;j++)
      {
          for(k=0;k<number_primary_cells_k+2;k++)
          {

          x = mesh_width_x1*(i-0.5);
          y = mesh_width_x2*(j-0.5);
          z = mesh_width_x3*(k-0.5);

          if(i!=number_primary_cells_i+1)
                {x += 0.5*mesh_width_x1;
                u_1_velocity_old[i][j][k]= initial_velocity.u1;
                u_1_velocity_new[i][j][k]= initial_velocity.u1;
                x -= 0.5*mesh_width_x1;}
          if(j!=number_primary_cells_j+1)
                {y += 0.5*mesh_width_x2;
                u_2_velocity_old[i][j][k]= initial_velocity.u2;
                u_2_velocity_new[i][j][k]= initial_velocity.u2;
                y -= 0.5*mesh_width_x2;}
          if(k!=number_primary_cells_k+1)
                {z += 0.5*mesh_width_x3;
                u_3_velocity_old[i][j][k]= initial_velocity.u3;
                u_3_velocity_new[i][j][k]= initial_velocity.u3;
                z -= 0.5*mesh_width_x3;}
          }
      }
  }


     /* compute the pressure corresponding to the initial velocity field and the initial velocity field */
     /* because the initial velocity field is chosen constant, it can be ignored in the computation     */
     /* of the initial pressure. We do have to take into account the initial body force field distribution */


     /* initialize the pressure */

     initialize_pressure( level_set, pressure,
                            momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                              surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                                scaled_density_u1, scaled_density_u2, scaled_density_u3,
                                   u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,
                                     mesh_width_x1, mesh_width_x2, mesh_width_x3,
                                      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                                        gravity, tolerance, actual_time_step_navier_stokes,
                                         rho_plus_over_rho_minus, continuous_surface_force_model,
                                           source_terms_in_momentum_predictor, maximum_iterations_allowed,
                                            boundary_faces);



     
     
    }
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

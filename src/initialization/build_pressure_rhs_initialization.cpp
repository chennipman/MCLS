#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to construct the rhside for the pressure correction                */
/*  specifically for the computation of the initial pressure field.             */
/*  The rhside will only contain the contribution of the divergence 		*/
/*  of the body force field due to the continuous surface force model.       	*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/  
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/********************************************************************************/
void build_pressure_rhs_initialization(
      Array3<double> initial_pressure_rhs,	     		// right hand side of pressure correction equation
					     		// excluding contributions 
					     		// inhomogeneous boundary conditions
      Array3<double> momentum_source_term_u_1,          	// complete source term for the momentum equation
                                                   	// in x1 direction=(-p,1+ g_1 +F1)
      Array3<double> momentum_source_term_u_2,          	// complete source term for the momentum equation
                                                   	// in x2 direction=(-p,2+ g_2 +F2)
      Array3<double> momentum_source_term_u_3,          	// complete source term for the momentum equation
                                                   	// in x3 direction=(-p,3+ g_3 +F3)
      Array3<double> csf_force_x1,	 	     		// source term of the momentum equation in x1 direction
					     		// defined on all u1 points (including boundaries)
      Array3<double> csf_force_x2,	             		// source term of the momentum equation in x2 direction
					     		// defined on all u1 points (including boundaries)
      Array3<double> csf_force_x3,		     		// source term of the momentum equation in x3 direction
					    		// defined on all u1 points (including boundaries)
      Array3<double> u_1_velocity_star, 	     		// velocity field at star time level x1 direction
      Array3<double> u_2_velocity_star, 	     		// velocity field at star time level x2 direction
      Array3<double> u_3_velocity_star,	     		// velocity field at star time level x3 direction
			      
      double mesh_width_x1,		     		// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		     		// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		     		// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	     		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	     		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	     		// number of primary (pressure) cells in x3 direction
      double actual_time_step_navier_stokes, 		// actual time step for Navier-Stokes solution algorithm 
      int continuous_surface_force_model,    		// =1, the continuous surface force model is applied
					     		// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor 		// =1, the source terms are applied in the momentum predictor
					     		// equation
					     		// =0, the source terms are applied in the momentum corrector
					     		// equation
      )
{
					  
      int i_index, j_index, k_index;  			// local variables for loop indexing
      double one_over_dx1	=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      double one_over_dx2	=    			// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
      double one_over_dx3	=    			// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);
      double one_over_actual_time_step_navier_stokes=	// 1/( actual time step used 		    
	    1.0/actual_time_step_navier_stokes;		// in navier stokes solution algorithm
					  

/* compute the divergence of the body force field                  */
/* in all the primary cells of the grid 			     	 */


      if(source_terms_in_momentum_predictor)
      {
	    for( i_index=1; i_index< number_primary_cells_i+1;i_index++)
	    {
		for( j_index=1; j_index< number_primary_cells_j+1 ; j_index++)
		{
		    for( k_index=1; k_index< number_primary_cells_k+1; k_index++)
		    {
			initial_pressure_rhs[i_index][j_index][k_index]=
			  one_over_dx1*(momentum_source_term_u_1[i_index][j_index][k_index] -
			                  momentum_source_term_u_1[i_index-1][j_index  ][k_index  ])+
			  one_over_dx2*(momentum_source_term_u_2[i_index][j_index][k_index] -
			                  momentum_source_term_u_2[i_index  ][j_index-1][k_index  ])+
			  one_over_dx3*(momentum_source_term_u_3[i_index][j_index][k_index] -
			                  momentum_source_term_u_3[i_index  ][j_index  ][k_index-1]);

		    }
		}
	    }
      }
     
 }


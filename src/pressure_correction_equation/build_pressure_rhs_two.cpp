#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to construct the pressure correction rhside                        */
/*  This function computes the contribution of the divergence of the velocity   */
/*  field and the source terms in the momentum predictor equation   		*/
/*  Use is made of the continuous surface model (CSF) for the surface tension	*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The pressure correction equation is a discretisation of                      */
/* -div(  rho_minus/rho grad p) = -f						*/
/* where the minus sign is introduced to make the 				*/
/* main diagonal element positive						*/
/* The boundary conditions are all of inhomogeneous Neumann type		*/
/* The initial pressure right hand side contains the divergence of the 	        */
/* velocity field u* etc., but it is here extended with the contributions 	*/
/* from the inhomogeneous boundary conditions that are applied			*/
/* The divergence of the velocity field is computed using a central 		*/
/* approximation of the derivatives:						*/
/*										*/
/*     										*/
/*   u3(i,j,k)-u3(i,j,k-1)  u2(i,j,k)-u2(i,j-1,k)  u1(i,j,k)-u1(i-1,j,k)	*/
/*   ------------------- + ------------------- + -------------------  = div u	*/
/*            d x3                  d x2                  d x1 			*/
/*										*/
/* Note that the boundary conditions are incorporated in the velocity field     */
/*			       	   n						*/
/* sxa= rho_minus/rho *( fxa -p,a) + ga						*/
/* where:									*/
/* rho_minus/rho = dimensionless density 					*/
/* sxa = source term of momentum equation in direction a			*/
/* fxa = continuous surface force in direction a				*/
/* ga  = component of the gravitational acceleration in direction a		*/		
/* this term is used because it should be 'balanced' ideally 			*/
/* There are currently two versions implemented					*/
/* Version 1:									*/			
/* set all sxa to zero, then 							*/
/* 	n+1	   * 	                      n+1		 		*/
/* div u    = div u  + div (1/rho *(f - grad p   )				*/
/* (Note that the constant gravitational accceleration does not contribute      */
/* to the divergence )								*/
/* this means that only the pressure gradient term is added in the correction   */
/* Version 2:									*/
/* compute all sxa terms							*/
/* 	n+1	   *                      n+1  n				*/
/* div u    = div u  + div(1/rho * grad (p  - p )				*/
/* this means the old pressure is included in the momentum predictor equation   */ 
/********************************************************************************/
EXPORT void build_pressure_rhs_two(
      Array3<double> initial_pressure_rhs,	     		// right hand side of pressure correction equation
					     		// excluding contributions 
					     		// inhomogeneous boundary conditions
      Array3<double> u_1_new_con_diff, 			// contains the convection and diffusion terms 
      Array3<double> u_2_new_con_diff, 			// contains the convection and diffusion terms
      Array3<double> u_3_new_con_diff,			// contains the convection and diffusion terms 
      double mesh_width_x1,		     		// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		     		// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		     		// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	     		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	     		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	     		// number of primary (pressure) cells in x3 direction
      double actual_time_step_navier_stokes, 		// actual time step for Navier-Stokes solution algorithm 
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
					  

/* compute the divergence of the convection and diffusion terms 	 */
/* in all the primary cells of the grid 			     	 */
      for( i_index=1; i_index< number_primary_cells_i+1;i_index++)
      {
	    for( j_index=1; j_index< number_primary_cells_j+1 ; j_index++)
	    {
		  for( k_index=1; k_index< number_primary_cells_k+1; k_index++)
		  {
		    initial_pressure_rhs[i_index][j_index][k_index]=
			  one_over_actual_time_step_navier_stokes*(
		    ( u_3_new_con_diff[i_index][j_index][k_index]-
			u_3_new_con_diff[i_index  ][j_index  ][k_index-1] )*one_over_dx3+
		    ( u_2_new_con_diff[i_index][j_index][k_index]-
			u_2_new_con_diff[i_index  ][j_index-1][k_index  ] )*one_over_dx2+
		    ( u_1_new_con_diff[i_index][j_index][k_index]-
			u_1_new_con_diff[i_index-1][j_index  ][k_index  ] )*one_over_dx1); 
			  
		  }  
	    }
      }

 }


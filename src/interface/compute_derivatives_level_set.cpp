#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the derivatives of the level-set field at the cell      */
/* faces.
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The level-set field is reinitialized by solving a diffusion type equation    */
/* in the vicinity of the interface.                                            */
/* the reinitialization equation is given in Sander's thesis as equations 6.4   */
/* which should be studied with equation 5.6                                    */
/* The standard Sussman term                                                    */
/* N_h(Phi^k,Phi^0) is weighted with (1-q(Phi^0)) and a term is added           */
/*     Phi^0- Phi^k   								*/
/*    --------------  * q(Phi^0)						*/
/*        Delta t'								*/
/* this function provides the first derivatives of the level-set field at the   */
/* cell faces for the diffusive fluxes.						*/
/********************************************************************************/
     void  compute_derivatives_level_set(				
		double ***level_set, 			// level set field at new time level
							// after convection and reinitialization
							// not mass conserving
		double ***d_level_set_d_x1,		// first partial derivative of
							// the level-set field wrt x1
							// second order central approximation
							// at cell faces 
		double ***d_level_set_d_x2,		// first partial derivative of 
							// the level-set field wrt x2
							// second order central approximation
							// at cell faces 
		double ***d_level_set_d_x3,		// first partial derivative of
 							// the level-set field wrt x3
							// second order central approximation
							// at cell faces 
		 double mesh_width_x1,			// grid spacing in x1 direction (uniform)
		 double mesh_width_x2,			// grid spacing in x2 direction (uniform)
		 double mesh_width_x3,			// grid spacing in x3 direction (uniform)
		 int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
		 int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
		 int number_primary_cells_k		// number of primary (pressure) cells in x3 direction
	      )

     {
      
		double one_over_dx1	=    		// 1/(grid spacing in x1 direction)
		    1.0/(mesh_width_x1);
		double one_over_dx2	=    		// 1/(grid spacing in x2 direction)
		    1.0/(mesh_width_x2);
		double one_over_dx3	=    		// 1/(grid spacing in x3 direction)
		    1.0/(mesh_width_x3);
		int i_index, j_index, k_index;  // local variables for loop indexing
     
      
      /* derivative in x1-direction */
      /* note the i-index runs from 0 to number_primary_cells_i so all the fluxes of the */
      /* non virtual cells are computed */
      
 	      for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
	      {
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			  {
			    d_level_set_d_x1[i_index][j_index][k_index]=
				(level_set[i_index+1][j_index  ][k_index  ]-
				    level_set[i_index][j_index][k_index])*one_over_dx1;
			  }
		    }
	      }

      /* derivative in x1-direction */
      /* note the i-index runs from 0 to number_primary_cells_i so all the fluxes of the */
      /* non virtual cells are computed */

	      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	      {
		    for(j_index=0;j_index<number_primary_cells_j+1;j_index++)
		    {
			  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			  {
			    d_level_set_d_x1[i_index][j_index][k_index]=
				(level_set[i_index  ][j_index+1][k_index  ]-
				    level_set[i_index][j_index][k_index])*one_over_dx1;
			  }
		    }
	      }

      /* derivative in x1-direction */
      /* note the i-index runs from 0 to number_primary_cells_i so all the fluxes of the */
      /* non virtual cells are computed */

 	      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	      {
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			  for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
			  {
			    d_level_set_d_x1[i_index][j_index][k_index]=
				(level_set[i_index][j_index][k_index]-
				    level_set[i_index][j_index][k_index])*one_over_dx1;
			  }
		    }
	      }
    }
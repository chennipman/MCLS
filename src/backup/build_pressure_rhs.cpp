#include "../headers/array.h"
/********************************************************************************/
/********************************************************************************/
/*  Function to construct the pressure correction rhside, for the CSF approach  */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The pressure correction equation is a discretisation of                      */
/* -div(  rho_minus/rho grad p) = -f						*/
/* where the minus sign is introduced to make the main diagonal element positive*/
/* The boundary conditions are all of inhomogeneous Neumann type		*/
/* The initial pressure right hand side contains the divergence of the 	        */
/* velocity field u* etc., but it is here extended with the contributions 	*/
/* from the inhomogeneous boundary conditions that are applied			*/
/********************************************************************************/
void build_pressure_rhs_boundary(
      Array3<double> initial_pressure_rhs,	     // right hand side of pressure correction equation
					     // excluding contributions 
					     // inhomogeneous boundary conditions
      Array3<double> final_pressure_rhs,	     // right hand side of pressure correction equation
					     // including contributions 
					     // inhomogeneous boundary conditions
      double mesh_width_x1,		     // grid spacing in x1 direction (uniform)
      double mesh_width_x2,		     // grid spacing in x2 direction (uniform)
      double mesh_width_x3,		     // grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	     // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	     // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	     // number of primary (pressure) cells in x3 direction
      Array3<double> pressure_boundary_condition_x1, // contribution to inhomogeneous boundary condition
						// pressure correction matrix for back/forward face
      Array3<double> pressure_boundary_condition_x2,	// contribution to inhomogeneous boundary condition
						// pressure correction matrix for left/right face
      Array3<double> pressure_boundary_condition_x3  // contribution to inhomogeneous boundary condition
						// pressure correction matrix for top/bottom face
     )

{
      int i_index, j_index, k_index;  // local variables for loop indexing
      int stencil_index;	      // index of coefficient in the stencil
      double one_over_dx1	=    // 1/(grid spacing in x1 direction)
	    1.0/(number_primary_cells_i);
						      
      double one_over_dx2	=    // 1/(grid spacing in x2 direction)
	    1.0/(number_primary_cells_j);
						      ; // 1/(grid spacing in x2 direction)^2 
      double one_over_dx3	=    // 1/(grid spacing in x3 direction)^2 
	    1.0/(number_primary_cells_k);

/*  initialize the final right hand side with minus the initial right hand side  */
/*  of the pressure correction equation. Contributions of inhomogeneous    	 */
/*  boundary condition are added later.					         */
	
    for(i_index=0;i_index<number_primary_cells_i;i_index++)
    {
	for(j_index=0;j_index<number_primary_cells_j;j_index++)
	{
	    for(k_index=0;k_index<number_primary_cells_k;k_index++)
	    {
		final_pressure_rhs[i_index][j_index][k_index]=
		    -1.0*initial_pressure_rhs[i_index][j_index][k_index];
		
	    }
	}
    }

    
/* add the contribution of the inhomogeneous neumann boundary condition to the right hand side */
/* this is for the front face: i_index=number_primary_cells_i-1*/

    for(j_index=0;j_index<number_primary_cells_j;j_index++)
    {
	for(k_index=0;k_index<number_primary_cells_k;k_index++)
	{
	  final_pressure_rhs[number_primary_cells_i-1][j_index][k_index]+=
	      one_over_dx1*pressure_boundary_condition_x1[1][j_index][k_index];
	  
	}
    }
    
/* add the contribution of the inhomogeneous neumann boundary condition to the right hand side */
/* this is for the back face i_index=0 */
    
    for(j_index=0;j_index<number_primary_cells_j;j_index++)
    {
	for(k_index=0;k_index<number_primary_cells_k;k_index++)
	{
	  final_pressure_rhs[0][j_index][k_index]+=
	      one_over_dx1*pressure_boundary_condition_x1[0][j_index][k_index];
	  
	}
    }
   
/* add the contribution of the inhomogeneous neumann boundary condition to the right hand side */
/* this is for the right face: j_index=number_primary_cells_j-1 */

    for(i_index=0;i_index<number_primary_cells_i;i_index++)
    {
	for(k_index=0;k_index<number_primary_cells_k;k_index++)
	{
	  final_pressure_rhs[i_index][number_primary_cells_j-1][k_index]+=
	      one_over_dx1*pressure_boundary_condition_x2[1][i_index][k_index];
	  
	}
    }
    
/* add the contribution of the inhomogeneous neumann boundary condition to the right hand side */
/* this is for the left face: j_index=0 */
    
    for(i_index=0;i_index<number_primary_cells_i;i_index++)
    {
	for(k_index=0;k_index<number_primary_cells_k;k_index++)
	{
	  final_pressure_rhs[i_index][0][k_index]+=
	      one_over_dx1*pressure_boundary_condition_x2[0][i_index][k_index];
	  
	}
    }
/* add the contribution of the inhomogeneous neumann boundary condition to the right hand side */
/* this is for the top face: k_index=number_primary_cells_k-1 */

    for(i_index=0;i_index<number_primary_cells_i;i_index++)
    {
	for(j_index=0;j_index<number_primary_cells_j;j_index++)
	{
	  final_pressure_rhs[i_index][number_primary_cells_j][number_primary_cells_k-1]+=
	      one_over_dx1*pressure_boundary_condition_x3[1][i_index][j_index];
	  
	}
    }
    
/* add the contribution of the inhomogeneous neumann boundary condition to the right hand side */
/* this is for the left face: j_index=0 */
    
    for(i_index=0;i_index<number_primary_cells_i;i_index++)
    {
	for(j_index=0;j_index<number_primary_cells_j;j_index++)
	{
	  final_pressure_rhs[i_index][j_index][0]+=
	      one_over_dx1*pressure_boundary_condition_x3[0][i_index][j_index];
	  
	}
    }
 
}
#include "../headers/array.h"


/********************************************************************************/
/********************************************************************************/
/*  Function to construct the pressure correction matrix, for the CSF approach  */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* The pressure correction equation is a discretisation of                      */
/* -div(  rho_minus/rho grad p) = -f							*/
/* where the minus sign is introduced to make the main diagonal element positive*/
/* The boundary conditions are all of inhomogeneous Neumann type			*/
/* at the moment, corresponding to Dirichlet for the velocity			*/
/* Therefore, the pressure matrix is singular with Dim(Ker(A))=1.			*/
/* The matrix is symmetric, and only the upper triangular part is computed.     */
/* In the Cartesian case the divgrad operator does not have any mixed 		*/
/* derivatives and the stencil contains six points, so the upper triangular     */
/* part of the stencil will only contain 4 points:                          	*/
/* - the central point		index 0	in the matrix				*/
/* - the forward neighbour	index 1	in the matrix					*/
/* - the righthand neighbour	index 2	in the matrix					*/
/* - the upper neighbour	index 3	in the matrix  				*/
/********************************************************************************/
void build_pressure_matrix (
      Array2<double> pressure_correction_matrix,   	// pressure correction matrix
      double mesh_width_x1,		     		// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		     		// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		     		// grid spacing in x3 direction (uniform)
      Array3<double> level_set_new, 		     	// level set field at new time level
      Array3<double> scaled_density_u1,			// scaled density for the controlvolumes
							// of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,			// scaled density for the controlvolumes
							// of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,			// scaled density for the controlvolumes
							// of the momentum equation in x3 direction
      int number_primary_cells_i,	     		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	     		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	     		// number of primary (pressure) cells in x3 direction
      double rho_plus_over_rho_minus	     		// ratio of density where (level set >0) and 
					     		// density where (level set < 0)
    
    )
{
      double compute_scaled_density(         	// compute the scaled density 		
	  double level_set_left, 			
	  double level_set_right,			
	  double rho_plus_over_rho_minus		
		);
      int map_index_pressure(                	// map 3-D array index to 1-D array
	  int i_index,			     		// index
	  int j_index, 				
	  int k_index,  			
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k		
      );
      double rho;                     	     	// density value
      double one_over_dx1_squared=           	// 1/(grid spacing in x1 direction)^2 
	  1.0/(mesh_width_x1*
		    mesh_width_x1);
      double one_over_dx2_squared=           	// 1/(grid spacing in x2 direction)^2 
	    1.0/(mesh_width_x2*
		    mesh_width_x2);
      double one_over_dx3_squared=           	// 1/(grid spacing in x3 direction)^2 
	    1.0/(mesh_width_x3*
		    mesh_width_x3);

      int i_index, j_index, k_index;         	// local variables for loop indexing
      int stencil_index;	      	     		// index of coefficient in the stencil
      int total_number_pressure_points;	     	// total number of points with pressure
      int one_dimensional_index;	     		// index of point in 1-D array
      int one_dimensional_index_source;	     	// index of point in 1-D array, when the
					     		// symmetry of the matrix is used to find
					     		// the value of the correct center coefficient
      int one_dimensional_index_target;      	// index of point in 1-D array, when the
					     		// symmetry of the matrix is used to find
					     		// the value of the correct center coefficient

						      
      
/* Step one: initialize the matrix to zero */
    
    total_number_pressure_points=number_primary_cells_i*
				  number_primary_cells_j*
				    number_primary_cells_k;

    for(stencil_index=0;stencil_index<4;stencil_index++)
    {
	for(one_dimensional_index=0;
		    one_dimensional_index<total_number_pressure_points;
							  one_dimensional_index++)
	{
	      pressure_correction_matrix[stencil_index][one_dimensional_index]=0.0;
	}
    }
	 

     
/* Step two: build off diagonal part of the upper triangular part of the matrix */

/* for all cells, but the front slice, set the i+ coefficient */

    for(i_index=1;i_index<number_primary_cells_i;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	{
	    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	    {
              rho=scaled_density_u1[i_index][j_index][k_index];
// 		rho=compute_scaled_density(level_set_new[i_index][j_index][k_index],
// 					  level_set_new[i_index+1][j_index][k_index], 
// 						rho_plus_over_rho_minus);
		one_dimensional_index=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[1][one_dimensional_index]=
			    -1.0*one_over_dx1_squared/rho;
	    }
	}
    }
    
/* for all cells, but the righthand slice, set the j+ coefficient */

    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j;j_index++)
	{
	    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	    {
               rho=scaled_density_u2[i_index][j_index][k_index];
// 		rho=compute_scaled_density(level_set_new[i_index][j_index][k_index],
// 					  level_set_new[i_index][j_index+1][k_index], 
// 						rho_plus_over_rho_minus);
		one_dimensional_index=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[2][one_dimensional_index]=
			    -1.0*one_over_dx2_squared/rho;
	    }
	}
    }

/* for all cells, but the top slice, set the k+ coefficient */

    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	{
	    for(k_index=1;k_index<number_primary_cells_k;k_index++)
	    {
               rho=scaled_density_u3[i_index][j_index][k_index];
// 		rho=compute_scaled_density(level_set_new[i_index][j_index][k_index],
// 					  level_set_new[i_index][j_index][k_index+1], 
// 						rho_plus_over_rho_minus);
		one_dimensional_index=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[3][one_dimensional_index]=
			    -1.0*one_over_dx3_squared/rho;
	    }
	}
    }
  


/********************************************************/  
/* matrix contributions from fluxes in the x1 direction */
/********************************************************/  

      
/* the center coefficient is adjusted for all cells that have a right hand neighbour */
/* all the real cells in one line, but the last one */
      
    for(i_index=1;i_index<number_primary_cells_i;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	{
	    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	    {
		one_dimensional_index=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[0][one_dimensional_index]-=
			   pressure_correction_matrix[1][one_dimensional_index];
	    }
	}
    }
    
/* the center coefficient is adjusted for all cells that have a left hand neighbour */
/* all the real cells in one line, but the first one */
/* the left hand coefficient is part of the LOWER triangular matrix, so use is made */
/* of the symmetry property */

    for(i_index=2;i_index<number_primary_cells_i+1;i_index++) // this index set is special
    {
	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	{
	    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	    {
		one_dimensional_index_target=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		one_dimensional_index_source=map_index_pressure(i_index-1,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[0][one_dimensional_index_target]-=
			   pressure_correction_matrix[1][one_dimensional_index_source];
	    }
	}
    }
    
/********************************************************/  
/* matrix contributions from fluxes in the x2 direction */
/********************************************************/  

 /* the center coefficient is adjusted for all cells that have now forward neighbour */
/* all the real cells in one line, but the last one */
      
    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j;j_index++)// this index set is special
	{
	    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	    {
		one_dimensional_index=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[0][one_dimensional_index]-=
			   pressure_correction_matrix[2][one_dimensional_index];
	    }
	}
    }
    
/* the center coefficient is adjusted for all cells that have a backward neighbour */
/* all the real cells in one line, but the first one */
/* the backward coefficient is part of the LOWER triangular matrix, so use is made */
/* of the symmetry property */

    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
    {
	for(j_index=2;j_index<number_primary_cells_j+1;j_index++)// this index set is special
	{
	    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	    {
		one_dimensional_index_target=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		one_dimensional_index_source=map_index_pressure(i_index,j_index-1,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[0][one_dimensional_index_target]-=
			   pressure_correction_matrix[2][one_dimensional_index_source];
		
	    }
	}
    }
       
 
    
/********************************************************/  
/* matrix contributions from fluxes in the x3 direction */
/********************************************************/  

 /* the center coefficient is adjusted for all cells that have a top neighbour */
  /* all the real cells in one line, but the last one */
      
    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	{
	    for(k_index=1;k_index<number_primary_cells_k;k_index++)// this index set is special
	    {
		one_dimensional_index=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[0][one_dimensional_index]-=
			   pressure_correction_matrix[3][one_dimensional_index];
	    }
	}
    }
    
/* the center coefficient is adjusted for all cells that have a bottom neighbour */
/* all the real cells in one line, but the first one */
/* the bottom coefficient is part of the LOWER triangular matrix, so use is made */
/* of the symmetry property */

    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	{
	    for(k_index=2;k_index<number_primary_cells_k+1;k_index++)// this index set is special
	    {
		one_dimensional_index_target=map_index_pressure(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		one_dimensional_index_source=map_index_pressure(i_index,j_index,k_index-1,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		
		pressure_correction_matrix[0][one_dimensional_index_target]-=
			   pressure_correction_matrix[3][one_dimensional_index_source];
	    }
	}
    }
       
}
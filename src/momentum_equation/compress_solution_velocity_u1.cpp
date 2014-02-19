#include "../headers/array.h"
/********************************************************************************/
/********************************************************************************/
/*  Function to construct the compressed solution vector for the velocity u1    */
/*  from the full solution vector						*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/********************************************************************************/
EXPORT void compress_solution_velocity_u1(
      Array3<double> full_solution,		     // 3-D array with solution, virtual cells included
      Array1<double> compressed_solution_vector,  // 1-D array with solution, virtual cells excluded
      int number_primary_cells_i,	     // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	     // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k	     // number of primary (pressure) cells in x3 direction
     )
 {
      int one_dimensional_index;	     	// index of point in 1-D array
      int i_index, j_index, k_index;  		// local variables for loop indexing
   
	
    for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
    {
	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	{
	    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	    {
		
		one_dimensional_index=map_index_u1(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j,  
								  number_primary_cells_k);
		compressed_solution_vector[one_dimensional_index]=
		    full_solution[i_index][j_index][k_index];
		
	    }
	}
    }

   
 }

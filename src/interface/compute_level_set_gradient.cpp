#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the gradient of the level-set field 	                */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The gradient of the level-set field is computed using a second order central */
/* approximation. At the boundary, the gradient is set to the nul vector        */
/* Note that the gradient is scaled by the vector with the mesh width           */
/* So it is an UNDIVIDED difference approximation                               */
/* and NOT the real gradient 							*/
/********************************************************************************/
EXPORT void   compute_level_set_gradient(				
	Array3<double> level_set, 			// level set field at new time level
						// after convection and reinitialization
						// not mass conserving
	Array3<double> d_level_set_d_x1,			// first partial derivative of
							// the level-set field wrt x1
							// second order central approximation
	Array3<double> d_level_set_d_x2,			// first partial derivative of 
							// the level-set field wrt x2
							// second order central approximation
	Array3<double> d_level_set_d_x3,			// first partial derivative of
 							// the level-set field wrt x3
							// second order central approximation
        int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
	int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
	int number_primary_cells_k		// number of primary (pressure) cells in x3 direction
	  )
   {
	int i_index, j_index, k_index;  			// local variables for loop indexing
// 	int vector_length= 					// length of the one dimensional array
// 	    (number_primary_cells_i+2)*				// that results from reshaping the 3D
// 		(number_primary_cells_j+2)*			// array of unknowns, with 1 virtual cell
// 		    (number_primary_cells_k+2);			// on all sides

		    
	/* set all first derivatives to zero, this takes care of all the virtual cells */

        set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, d_level_set_d_x1, 0.0);
        set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, d_level_set_d_x2, 0.0);
        set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, d_level_set_d_x3, 0.0);

	/* compute the gradient in all nonvirtual cells, using a central approximation        */
	/* to the first derivative 						      	      */
	/* note the nonvirtual cells start at index 1 and end at index number_primary_cells_* */
	
	
	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	{
	    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	    {
		    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
		  
		    d_level_set_d_x1[i_index][j_index][k_index]=
				     0.5*(level_set[i_index+1][j_index  ][k_index  ]-
					    level_set[i_index-1][j_index  ][k_index  ]);
		    d_level_set_d_x2[i_index][j_index][k_index]=
				     0.5*(level_set[i_index  ][j_index+1][k_index  ]-
					    level_set[i_index  ][j_index-1][k_index  ]);
		    d_level_set_d_x3[i_index][j_index][k_index]=
				     0.5*(level_set[i_index  ][j_index  ][k_index+1]-
					    level_set[i_index  ][j_index  ][k_index-1]);
		}

	      
	    }
	}
     
   }

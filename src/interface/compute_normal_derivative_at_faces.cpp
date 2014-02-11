#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the normal derivative of a cell centered field at       */
/*  the cell faces	                                                        */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The normal derivative of a scalar field defined at the cell centers is       */
/* approximated with central differences. Note that in this case the            */
/* actual derivative is computed, so a divided difference.          h           */
/********************************************************************************/
   void   compute_normal_derivative_at_faces(				
	double ***scalar_field, 			// scalar field defined at cell centers
	double ***d_field_d_x1_face,			// first normal derivative of
							// scalar field at cell face orthogonal
							// to x1 direction                       
	double ***d_field_d_x2_face,			// first normal derivative of
							// scalar field at cell face orthogonal
							// to x2 direction
	double ***d_field_d_x3_face,			// first normal derivative of
							// scalar field at cell face orthogonal
							// to x3 direction
        int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	double mesh_width_x1,				// grid spacing in x1 direction (uniform)
	double mesh_width_x2,				// grid spacing in x2 direction (uniform)
	double mesh_width_x3				// grid spacing in x3 direction (uniform)
	  )
   {
     
	/* function definitions */ 
	
       void set_constant_matrix2(
	    int first_dimension,			// number of elements in first dimension
	    int second_dimension,			// number of elements in second dimension
	    int third_dimension,			// number of elements in third dimension
	    double ***matrix2_to_set,			// the name of the array that has to be set
	    double constant_value			// the constant value the vector has to be set to
	    );
	
	
	int i_index, j_index, k_index;  		// local variables for loop indexing
	int vector_length= 				// length of the one dimensional array
	    (number_primary_cells_i+2)*			// that results from reshaping the 3D
		(number_primary_cells_j+2)*		// array of unknowns, with 1 virtual cell
		    (number_primary_cells_k+2);		// on all sides

	double one_over_dx1	=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
	double one_over_dx2	=    			// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
	double one_over_dx3	=    			// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);

	
	
	for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
	{
	    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	    {
		    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
		  
		    d_field_d_x1_face[i_index][j_index][k_index]=
				     (scalar_field[i_index+1][j_index  ][k_index  ]-
					    scalar_field[i_index][j_index  ][k_index  ])*one_over_dx1;
		}

	      
	    }
	}
	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	{
	    for(j_index=0;j_index<number_primary_cells_j+1;j_index++)
	    {
		    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
		  
		    d_field_d_x2_face[i_index][j_index][k_index]=
				     (scalar_field[i_index  ][j_index+1][k_index  ]-
					    scalar_field[i_index  ][j_index][k_index  ])*one_over_dx2;
		}

	      
	    }
	}
     
	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	{
	    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	    {
		    for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
		{
		  
		    d_field_d_x3_face[i_index][j_index][k_index]=
				     (scalar_field[i_index  ][j_index  ][k_index+1]-
					    scalar_field[i_index  ][j_index  ][k_index])*one_over_dx3;
		}

	      
	    }
	}
   }
#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to apply homogeneous neumann boundary condition			*/
/*			to a cell centered field       				*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function will extrapolate centered field values from the real cells     */
/* to the virtual cells. The field can be any quantity that is defined in the   */
/* centers of the primary cells 						*/
/* We need to consider the 6 faces as well as the 12 edges and the 8 vertices   */
/********************************************************************************/
//
void  field_neumann_boundary(
      Array3<double> field, 			// cell centered field
      int number_primary_cells_i,	// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k	// number of primary (pressure) cells in x3 direction
      
)

{
     /* Ranges for cells to extrapolate */
     /* edge 1 i_index=0, j_index=1..number_primary_cells_j, k_index=0 */
     /* edge 2 i_index=number_primary_cells_i+1, j_index=1..number_primary_cells_j, k_index=0 */
     /* edge 3 i_index=1..number_primary_cells_i, j_index=0, k_index=0 */
     /* edge 4 i_index=1..number_primary_cells_i, j_index=number_primary_cells_j+1, k_index=0 */
     /* edge 5 i_index=0, j_index=1..number_primary_cells_j, k_index=number_primary_cells_k+1 */
     /* edge 6 i_index=number_primary_cells_i+1, j_index=1..number_primary_cells_j, k_index=number_primary_cells_k+1 */
     /* edge 7 i_index=1..number_primary_cells_i, j_index=0, k_index=number_primary_cells_k+1 */
     /* edge 8 i_index=1..number_primary_cells_i, j_index=number_primary_cells_j+1, k_index=number_primary_cells_k+1 */
     /* edge 9 i_index=number_primary_cells_i+1, j_index=0, k_index=1..number_primary_cells_k */
     /* edge 10 i_index=number_primary_cells_i+1, j_index=number_primary_cells_j+1, k_index=1..number_primary_cells_k */
     /* edge 11 i_index=0, j_index=0, k_index=1..number_primary_cells_k */
     /* edge 12 i_index=0, j_index=number_primary_cells_j+1, k_index=1..number_primary_cells_k */

     double table_edge_extrapolation_range[12][6]=
     { 
	{			1,			 1,			  1,  number_primary_cells_j, 			    1, 			     1},
	{  number_primary_cells_i,  number_primary_cells_i,			  1,  number_primary_cells_j,			    1,			     1},
	{			1,  number_primary_cells_i,			  1,			   1,                       1,  		     1},
	{			1,  number_primary_cells_i,  number_primary_cells_j,  number_primary_cells_j,			    1,			     1},
	{			1,			 1,			  1,  number_primary_cells_j,  number_primary_cells_k,  number_primary_cells_k},
	{  number_primary_cells_i,  number_primary_cells_i,			  1,  number_primary_cells_j,  number_primary_cells_k,  number_primary_cells_k},
	{			1,  number_primary_cells_i,			  1,			   1,  number_primary_cells_k,  number_primary_cells_k},
	{			1,  number_primary_cells_i,  number_primary_cells_j,  number_primary_cells_j,  number_primary_cells_k,  number_primary_cells_k},
	{  number_primary_cells_i,  number_primary_cells_i,			  1,                       1,                       1,  number_primary_cells_k},
	{  number_primary_cells_i,  number_primary_cells_i,  number_primary_cells_j,  number_primary_cells_j,                       1,  number_primary_cells_k},
	{			1,                       1,			  1,			   1,                       1,  number_primary_cells_k},
	{			1,                       1,  number_primary_cells_j,  number_primary_cells_j,                       1,  number_primary_cells_k}
	
     };

     /* Shifts for cells to extrapolate along edges */
     
     /* edge 1   i_shift=-1, j_shift= 0, k_shift=-1 */
     /* edge 2   i_shift=+1, j_shift= 0, k_shift=-1 */
     /* edge 3   i_shift= 0, j_shift=-1, k_shift=-1 */
     /* edge 4   i_shift= 0, j_shift=+1, k_shift=-1 */
     /* edge 5   i_shift=-1, j_shift= 0, k_shift=+1 */
     /* edge 6   i_shift=+1, j_shift= 0, k_shift=+1 */
     /* edge 7   i_shift= 0, j_shift=-1, k_shift=+1 */
     /* edge 8   i_shift= 0, j_shift=+1, k_shift=+1 */
     /* edge 9   i_shift=+1, j_shift=-1, k_shift= 0 */
     /* edge 10  i_shift=+1, j_shift=+1, k_shift= 0 */
     /* edge 11  i_shift=-1, j_shift=-1, k_shift= 0 */
     /* edge 12  i_shift=-1, j_shift=+1, k_shift= 0 */
     
     
     int table_edge_extrapolation_shift[12][3]=
     { 
	{-1, 0,-1 },
	{ 1, 0,-1 },
	{ 0,-1,-1 },
	{ 0, 1,-1 },
	{-1, 0, 1 },
	{ 1, 0, 1 },
	{ 0,-1, 1 },
	{ 0, 1, 1 },
	{ 1,-1, 0 },
	{ 1, 1, 0 },
	{-1,-1, 0 },
	{-1, 1, 0 }
	
     };
     
     double table_face_extrapolation_range[6][6]=
     { 
	{			1,			 1,			  1,  number_primary_cells_j, 			    1,  number_primary_cells_k},
	{  number_primary_cells_i,  number_primary_cells_i,			  1,  number_primary_cells_j,			    1,	number_primary_cells_k},
	{			1,  number_primary_cells_i,			  1,			   1,                       1,  number_primary_cells_k},
	{			1,  number_primary_cells_i,  number_primary_cells_j,  number_primary_cells_j,			    1,	number_primary_cells_k},
	{			1,  number_primary_cells_i,			  1,  number_primary_cells_j,			    1,			     1},
	{			1,  number_primary_cells_i,			  1,  number_primary_cells_j,  number_primary_cells_k,  number_primary_cells_k}
	
     };

     /* Shifts for cells to extrapolate along faces */
     
     /* edge 1   i_shift=-1, j_shift= 0, k_shift= 0  */
     /* edge 2   i_shift= 1, j_shift= 0, k_shift= 0  */
     /* edge 3   i_shift= 0, j_shift=-1, k_shift= 0  */
     /* edge 4   i_shift= 0, j_shift= 1, k_shift= 0  */
     /* edge 5   i_shift= 0, j_shift= 0, k_shift= -1 */
     /* edge 6   i_shift= 0, j_shift= 0, k_shift=  1 */

     int table_face_extrapolation_shift[6][3]=
     { 
	{-1, 0, 0 },
	{ 1, 0, 0 },
	{ 0,-1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, -1},
	{ 0, 0,  1}
     };
     
   int i_index, j_index, k_index;  	// local variables for loop indexing
   int face_index, edge_index;		// local variables for loop indexing
   int i_shift, j_shift, k_shift;  	// shifts underlying the interpolation
   int i_index_start; 	 		// first i_index for extrapolation
   int i_index_end; 	 		// last i_index for extrapolation
   int j_index_start; 	 		// first j_index for extrapolation
   int j_index_end; 	 		// last j_index for extrapolation
   int k_index_start; 	 		// first k_index for extrapolation
   int k_index_end; 	 		// last k_index for extrapolation
    
  
  
  
  
  

/* Apply a homogeneous Neumann boundary condition to the field */
/* This means the lower index virtual cell will get the value of the first real cell */
/* and the higher index virtual cell has the value of the last real cell */


      /* handle the six faces */
    
      for(face_index=0;face_index<6;face_index++){

	  i_index_start	= table_face_extrapolation_range[face_index][0];
	  i_index_end		= table_face_extrapolation_range[face_index][1];
	  j_index_start	= table_face_extrapolation_range[face_index][2];
	  j_index_end		= table_face_extrapolation_range[face_index][3];
	  k_index_start	= table_face_extrapolation_range[face_index][4];
	  k_index_end		= table_face_extrapolation_range[face_index][5];
	  i_shift		= table_face_extrapolation_shift[face_index][0];
	  j_shift		= table_face_extrapolation_shift[face_index][1];
	  k_shift		= table_face_extrapolation_shift[face_index][2];

	  for(i_index=i_index_start;i_index<i_index_end+1;i_index++)
	  {
	      for(j_index=j_index_start;j_index<j_index_end+1;j_index++)
	      {
		  for(k_index=k_index_start;k_index<k_index_end+1;k_index++)
		  {
		      field[i_index+i_shift][j_index+j_shift][k_index+k_shift] =  
				      field[i_index][j_index][k_index];

		  }  
	      }  
	  }
      }
      
      /* handle the 12 edges */
    
      for(edge_index=0;edge_index<12;edge_index++){

	  i_index_start	= table_edge_extrapolation_range[edge_index][0];
	  i_index_end		= table_edge_extrapolation_range[edge_index][1];
	  j_index_start	= table_edge_extrapolation_range[edge_index][2];
	  j_index_end		= table_edge_extrapolation_range[edge_index][3];
	  k_index_start	= table_edge_extrapolation_range[edge_index][4];
	  k_index_end		= table_edge_extrapolation_range[edge_index][5];
	  i_shift		= table_edge_extrapolation_shift[edge_index][0];
	  j_shift		= table_edge_extrapolation_shift[edge_index][1];
	  k_shift		= table_edge_extrapolation_shift[edge_index][2];

	  for(i_index=i_index_start;i_index<i_index_end+1;i_index++)
	  {
	      for(j_index=j_index_start;j_index<j_index_end+1;j_index++)
	      {
		  for(k_index=k_index_start;k_index<k_index_end+1;k_index++)
		  {
		      field[i_index+i_shift][j_index+j_shift][k_index+k_shift] =  
				      field[i_index][j_index][k_index];

		  }  
	      }  
	  }
      }
      
      
      /* handle the eight corners */
      
      
      /* corner (0,0,0)  from 
       *		(1,1,1) */
      /* corner (number_primary_cells_i+1,0,0) from 
       *		(number_primary_cells_i,1,1) */
      /* corner (number_primary_cells_i+1, number_primary_cells_j+1,0) from 
       *		(number_primary_cells_i, number_primary_cells_j,1) */
      /* corner (0, number_primary_cells_j+1,0) from 
       *		(1, number_primary_cells_j,1) */
      /* corner (0,0,number_primary_cells_k+1) from 
       *		(1,1,number_primary_cells_k) */
      /* corner (number_primary_cells_i+1,0,number_primary_cells_k+1) from 
       *		(number_primary_cells_i,1,number_primary_cells_k) */
      /* corner (number_primary_cells_i+1,number_primary_cells_j+1, number_primary_cells_k+1) from 
         		(number_primary_cells_i, number_primary_cells_j, number_primary_cells_k) */
      /* corner (0, number_primary_cells_j+1,number_primary_cells_k+1) from
        		 (1, number_primary_cells_j, number_primary_cells_k) */
      field[0][0][0]							=
			    field[1][1][1];
      field[number_primary_cells_i+1][0][0]				=
			    field[number_primary_cells_i][1][1];
      field[number_primary_cells_i+1][number_primary_cells_j+1][0]	=
			    field[number_primary_cells_i][number_primary_cells_j][0];
      field[0][number_primary_cells_j+1][0]				=
			    field[1][number_primary_cells_j][1];
      field[0][0][number_primary_cells_k+1]				=
			    field[1][1][number_primary_cells_k];
      field[number_primary_cells_i+1][0][number_primary_cells_k+1]	=
			    field[number_primary_cells_i][1][number_primary_cells_k];
      field[number_primary_cells_i+1][number_primary_cells_j+1][number_primary_cells_k+1]	=
			    field[number_primary_cells_i][number_primary_cells_j][number_primary_cells_k];
      field[0][number_primary_cells_j+1][number_primary_cells_k+1]				=
			    field[1][number_primary_cells_j][number_primary_cells_k];
     
}
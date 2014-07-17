#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to establish the validity of the volume of fluid field             */
/*  method.                                                                     */
/*                                                                              */
/*  Programmer  : Duncan van der Heul                                           */
/*  Date        : 10-03-2013                                                    */
/*  Update      :                                                               */
/********************************************************************************/
/* Notes                                                                        */
/*    Check for cells with volume of fluid outside the                          */
/*    interval [-eps,1+eps] , with eps the volume of fluid tolerance            */
/*    and 'vapor' cells.                                                        */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

      EXPORT int determine_cell_is_pure(
        Array3<double> level_set,                               // level set field 
                                                                // mass conserving
        int i_index,						// i-index of cell to be analyzed
        int j_index,						// j-index of cell to be analyzed
        int k_index,						// k-index of cell to be analyzed
        int number_primary_cells_i,                             // number of primary (pressure) cells in x1 direction
        int number_primary_cells_j,                             // number of primary (pressure) cells in x2 direction
        int number_primary_cells_k                              // number of primary (pressure) cells in x3 direction
     	)                                                          
{                                                                  
                                                                   
      const int number_indices_neighbours_2D=4;			// number of neigbouring cells
      								// to consider, two dimensional cases
      const int number_indices_neighbours_3D=14;			// number of neighbouring cells
      								// to consider, three dimensional case
                                                                
      int number_indices_neighbours;				// actual number of neigbouring
      								// cells to consider, after the dimension
      								// of the problem has been determined
      int sign_of_centerpoint;					// sign of level set field at the point
      								// under consideration
      int sign_of_this_point;					// sign of level set field at the neighbouring point
      								// under consideration
      int neighbour_index;
      
      int index_set_neighbours[number_indices_neighbours_3D][3];// set of all index shifts to neighbouring
      								// points to consider
      int pure_cell=1;						//=1, cell is a pure cell,=0 cell is a vapor cell								
      
  /* compute the index-set of neighbours for a vapor cell */
  
  
  /* compute how many cells are in the index set */
  
  if( number_primary_cells_i==1 || number_primary_cells_j==1 || number_primary_cells_k==1)
  {  
  	/* 2-D case  */


    	number_indices_neighbours=number_indices_neighbours_2D;
    	if(number_primary_cells_i==1)
    	{
    		/* 2_D case, constant i_index */
    		/* set up the default table and copy */
    		
    		int table[number_indices_neighbours_3D][3]=
    		{ 
    			{ 0,-1, -1}, // 1
    			{ 0,-1,  1}, // 2
    			{ 0, 1, -1}, // 3
    			{ 0, 1,  1}, // 4
   			{ 0, 0,  0}, // 5
   			{ 0, 0,  0}, // 6
   			{ 0, 0,  0}, // 7
   			{ 0, 0,  0}, // 8
   			{ 0, 0,  0}, // 9
   			{ 0, 0,  0}, //10
   			{ 0, 0,  0}, //11
    			{ 0, 0,  0}, //12
  			{ 0, 0,  0}, //13
   			{ 0, 0,  0}  //14
 		};
      		std::copy(table[0], table[0]+3*number_indices_neighbours_3D, index_set_neighbours[0]);
   	}
    	else 
    	{
    		if(number_primary_cells_j==1)
    		{
    		
    		/* 2_D case, constant j_index */
    		/* set up the default table and copy */

     			int table[number_indices_neighbours_3D][3]=
     			{ 
     				{ -1, 0,-1},// 1
     				{ -1, 0, 1},// 2
     				{  1, 0,-1},// 3
     				{  1, 0, 1},// 4
     				{ 0, 0,  0},// 5
     				{ 0, 0,  0},// 6
     				{ 0, 0,  0},// 7
     				{ 0, 0,  0},// 8
      				{ 0, 0,  0},// 9
     				{ 0, 0,  0},//10
     				{ 0, 0,  0},//11
     				{ 0, 0,  0},//12
     				{ 0, 0,  0},//13
     				{ 0, 0,  0} //14
    			};
      		        std::copy(table[0], table[0]+3*number_indices_neighbours_3D, index_set_neighbours[0]);
   			
 		}
    		else
    		{
    		
     		/* 2_D case, constant k_index */
     		/* set up the default table and copy */

     			int table[number_indices_neighbours_3D][3]=
     			{ 
     				{-1, -1, 0},// 1
     				{-1,  1, 0},// 2 
     				{ 1, -1, 0},// 3
     				{ 1,  1, 0},// 4
     				{ 0, 0,  0},// 5
     				{ 0, 0,  0},// 6
     				{ 0, 0,  0},// 7
     				{ 0, 0,  0},// 8
      				{ 0, 0,  0},// 9
     				{ 0, 0,  0},//10
     				{ 0, 0,  0},//11
     				{ 0, 0,  0},//12
     				{ 0, 0,  0},//13
     				{ 0, 0,  0} //14
  			};
     		        std::copy(table[0], table[0]+3*number_indices_neighbours_3D, index_set_neighbours[0]);
  		}
    	}
  }
  else
  {  
  	/* 3-D case  */
        /* set up the default table and copy */

  	number_indices_neighbours=number_indices_neighbours_3D;
    		int table[number_indices_neighbours_3D][3]=
    		{ 
 			{1,   1,  1},// 1
    			{1,  -1,  1},// 2
    			{-1,  1,  1},// 3
    			{-1, -1,  1},// 4
  			{ 1,  1, -1},// 5
    			{ 1, -1, -1},// 6
    			{-1,  1, -1},// 7
    			{-1, -1, -1},// 8
    		        {  1, 0,  0},// 9
     		        { -1, 0,  0},//10
     		        { 0,  1,  0},//11
     		        { 0, -1,  0},//12
       		        {  0, 0,  1},//13
     		        {  0, 0, -1} //14
  		};
   		
    		std::copy(table[0], table[0]+3*number_indices_neighbours_3D, index_set_neighbours[0]);

  }
  
   	/* determine the sign of the level set field at the center point */
  	
  	sign_of_centerpoint=(int) round(sign(1.0,level_set[i_index][j_index][k_index]));
 
  	/* loop over all neighbours in the index set, and compare signs */
  	
  	for(neighbour_index=0;neighbour_index<number_indices_neighbours;neighbour_index++)
  	{
  		int neighbour_index_i=i_index+index_set_neighbours[neighbour_index][0];
  		int neighbour_index_j=j_index+index_set_neighbours[neighbour_index][1];
  		int neighbour_index_k=k_index+index_set_neighbours[neighbour_index][2];

  		sign_of_this_point=
        		(int) round(sign(1.0,
        		  	  level_set[neighbour_index_i][neighbour_index_j][neighbour_index_k]));
        	if(sign_of_centerpoint!=sign_of_this_point)
        	{
        		pure_cell=0;
        		break;
        		
        	}
        }
        
        return pure_cell;
 }

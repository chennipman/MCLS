#include "../headers/array.h"

/********************************************************************************/
/********************************************************************************/
/*  Function to apply neumann boundary condition to the level set field         */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* All the level set values for the internal cells are updated in               */
/* advect_level_set. To be able to compute the curvature and normal vector near */
/* the boundary, every boundary has one layer of virtual cells.                 */
/* The value of the level set field in these virtual cells is determined by     */
/* the application of a homogeneous neumann boundary condition                	*/
/* This corresponds to constant extrapolation of the last real layers of cells  */
/* to the virtual cell layers                                                   */
/********************************************************************************/
//
void  level_set_neumann_boundary(
      Array3<double> level_set_new, 		// level set field at new time level
      int number_primary_cells_i,	// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k	// number of primary (pressure) cells in x3 direction
      
)

{

  int i_index, j_index, k_index;  // local variables for loop indexing

/* Apply a homogeneous Neumann boundary condition to the level-set field */
/* This means the lower index virtual cell will get the value of the first real cell */
/* and the higher index virtual cell has the value of the last real cell */
/* The indication of the face refers to a standard right handed frame of reference */
  
      for(j_index=1;j_index<number_primary_cells_j;j_index++)
      {
	  for(k_index=1;k_index<number_primary_cells_k;k_index++)
	  {

	    /* slice with i_index = 0 'back face' */  
     
	      level_set_new[0][j_index][k_index]  =  
				level_set_new[1][j_index][k_index];

  
	  }  
  
      }  
     
      for(j_index=1;j_index<number_primary_cells_j;j_index++)
      {
	  for(k_index=1;k_index<number_primary_cells_k;k_index++)
	  {

	    /* slice with i_index = number_primary_cells_i+1 'front face' */  

	      level_set_new[number_primary_cells_i+1][j_index][k_index]  =  
				level_set_new[number_primary_cells_i][j_index][k_index];

  
	  }  
  
      }  
     
 
      for(i_index=1;i_index<number_primary_cells_i;i_index++)
      {
	  for(k_index=1;k_index<number_primary_cells_k;k_index++)
	  {
	    /* slice with j_index = 0 'left face' */  
    
	      level_set_new[i_index][0][k_index]  =  
				level_set_new[i_index][1][k_index];

  
	  }  
   
  
      }  
      for(i_index=1;i_index<number_primary_cells_i;i_index++)
      {
	  for(k_index=1;k_index<number_primary_cells_k;k_index++)
	  {
	    /* slice with j_index = 0 'left face' */  
    
	      level_set_new[i_index][number_primary_cells_j+1][k_index]  =  
				level_set_new[i_index][number_primary_cells_j][k_index];

  
	  }  
   
  
      }  
     
      for(i_index=1;i_index<number_primary_cells_i;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j;j_index++)
	  {
	    /* slice with k_index = 0 'bottom face' */  
	      level_set_new[i_index][j_index][0] =  
		    level_set_new[i_index][j_index][1];

 	  }  
   
      }  

      for(i_index=1;i_index<number_primary_cells_i;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j;j_index++)
	  {
	    /* slice with k_index = number_primary_cells_k+1 'top face' */
	      level_set_new[i_index][j_index][0] =  
		    level_set_new[i_index][j_index][1];

 	  }  
   
      }  
     
}
/********************************************************************************/
/*  Function to map the index of a point in the grid to 1-dimensional array     */
/*               						*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The mapping is different for the three components of the velocity and for    */
/* the pressure, because of the staggering of the unknowns.              .      */
/* This function should be in-lined.                                            */
/********************************************************************************/
  int map_index_u1(
      int i_index,				// i_index of the point in the 3-d array
      int j_index, 				// j_index of the point in the 3-d array
      int k_index,  				// k_index of the point in the 3-d array
      int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k		// number of primary (pressure) cells in x3 direction
      )
  {
    int one_dimensional_index;			// index of the point in the 1-d array
    
    /* in construction of the matrix and rhside for the u_1 matrix the */
    /* j and k index run from 1 to number_primary_cells_j/k to have a  */
    /* nice correspondence with the fields that do include the virtual */
    /* cells							       */
    
    one_dimensional_index=i_index+(j_index-1)*(number_primary_cells_i+1)+
		(k_index-1)*((number_primary_cells_i+1)*number_primary_cells_j);
    
  return  one_dimensional_index; 
  }
  
  
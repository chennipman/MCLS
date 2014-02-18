#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to apply inhomogeneous neumann boundary condition			*/
/*			to a cell centered field       				*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* This function uses linear extrapolation to set the virtual cells. This is    */
/* necessary for the pressure, in the case of gravity.                          */
/* We need to consider the 6 faces as well as the 12 edges and the 8 vertices   */
/* to minimize the coding, all the necessary shifts are tabulated for the       */
/* different faces, edges and corners						*/
/********************************************************************************/
//
EXPORT void  field_extrapolate_boundary(
      Array3<double> field, 			// cell centered field
      int number_primary_cells_i,	// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k	// number of primary (pressure) cells in x3 direction
      
)
{
	field_neumann_boundary(field, number_primary_cells_i, number_primary_cells_j, 
				  	number_primary_cells_k);

     
}

#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to shift the level-set field from the new to the old time-level    */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/*  											*/
/*  											*/
/********************************************************************************/

EXPORT void shift_interface(
      Array3<double> level_set_new, 		// level set field at new time level
						// mass conserving
      Array3<double> level_set_old, 		// level set field at old time level
						// mass conserving
      int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k		// number of primary (pressure) cells in x3 direction
	 )		      
      {
      /* shift the level-set field from the new to old time level*/
      

        copy_cell_centered_field(level_set_new, level_set_old, 
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      

      }

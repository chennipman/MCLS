#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to check the symmetry of all quantities in the case of a           */
/*  symmetric problem 								*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
    void check_symmetry_scalars(
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      Array3<double> field					// field of which symmetry must be checked
	  )
    {
    
      int i_index, j_index, k_index;  		// local variables for loop indexing
      int symmetry_plane_index;			// index of the symmetry plane 
     
      double max_error_ij=0;				// max difference between i and j constant planes
      double max_error_ik=0;				// max difference between i and k constant planes
      double max_error_jk=0;				// max difference between j and j constant planes

      
      	   for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
          {
		for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
           	{
			max_error_ij=std::max(fabs(
			field[symmetry_plane_index][i_index][j_index]-
				field[i_index][symmetry_plane_index][j_index]),
			    		max_error_ij);
			max_error_ik=std::max(fabs(
			field[symmetry_plane_index][i_index][j_index]-
				field[i_index][j_index][symmetry_plane_index]),
			    		max_error_ik);
			max_error_jk=std::max(fabs(
			field[i_index][symmetry_plane_index][j_index]-
				field[i_index][j_index][symmetry_plane_index]),
			    		max_error_jk);
			
            	}
          }

          std::cerr<<"max error ij "<< max_error_ij <<"\n";
          std::cerr<<"max error ik "<< max_error_ik <<"\n";
          std::cerr<<"max error jk "<< max_error_jk <<"\n";
    }     
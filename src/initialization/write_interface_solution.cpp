
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/********************************************************************************/
/********************************************************************************/
/*  Function to dump the solution quantities related to the interface to file   */
/*  for direct comparison with the fortran code.					*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
 void write_interface_solution(
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	  double mesh_width_x3,			// grid spacing in x3 direction (uniform)
	  double ***level_set,			// level-set field
	  double ***volume_of_fluid,			// volume of fluid field
	  double ***curvature,			// interface curvature
	  double ***unsmoothed_curvature		// interface curvature without smoothing
      )
 {
      double cell_center_x1;				// cell center x1 coordinate
      double cell_center_x2;				// cell center x2 coordinate
      double cell_center_x3;				// cell center x3 coordinate
      int i_index, j_index, k_index;  		// local variables for loop indexing
     
		  ofstream interface_solution( "interface_solution_c++.txt");
		  if(!interface_solution)
		  {
		      /* the contructor returned a 0-pointer :-( */
		      cout << "Cannot open file.\n";
		      exit(1);
		  }

	  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	  {
	  	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  	{
      			for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      			{ 
	    		cell_center_x1=(i_index-0.5)*mesh_width_x1;
	    		cell_center_x2=(j_index-0.5)*mesh_width_x2;
	    		cell_center_x3=(k_index-0.5)*mesh_width_x3;
			interface_solution.setf(ios::fixed, ios::floatfield);
			interface_solution.width(4);
			interface_solution<<i_index<<"    "<<j_index<<"    "<<k_index<<"    ";
			interface_solution.precision(8);
			interface_solution<<cell_center_x1<<"    ";
			interface_solution<<cell_center_x2<<"    ";
			interface_solution<<cell_center_x3<<"    ";
			interface_solution<<curvature[i_index][j_index][k_index]<<"    ";
			interface_solution<<level_set[i_index][j_index][k_index]<<"    ";
			interface_solution<<volume_of_fluid[i_index][j_index][k_index]<<" \n";
			}
		}
	  }
 }
#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to initialize the level-set field for the rayleigh taylor case  	*/
/*  										*/
/*  										*/
/*  Programmer	: Coen Hennipman						*/
/*  Date	: 19-06-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* It is assumed the free surfaces are aligned with the coordinate planes.	*/
/* The applied cosine disturbance is hardcoded and 0.05.			*/
/********************************************************************************/

EXPORT void initialize_rayleigh_taylor(
      surface *free_surfaces, 				// array with the definition of the free surfaces
      int number_of_free_surfaces, 			// number of bubbles in the domain (<10)
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      Array3<double> level_set				// level-set field
	    
    )
    {
    double x1_coordinate_cell_center;			// x1 coordinate cell center
    double x2_coordinate_cell_center;			// x2 coordinate cell center
    double x3_coordinate_cell_center;			// x3 coordinate cell center
    double surface_height;				// distance between origin and free surface
    int i_index, j_index, k_index;  			// local variables for loop indexing
    int surface_is_active;				//(=1) the presence of this free surface has to be
							// taken into account, (=0), ignore this surface
    int surface_orientation;				// the nonzero component of the normal vector on this plane
    double PI=atan(1)*4;

       if (number_of_free_surfaces > 1)
       {
	std::cerr << "**************************************************** \n";
	std::cerr << "ERROR \n";
	std::cerr << "more than one free surface in the rayleigh taylor case, terminating...";
	std::cerr << "in function initialize_rayleigh_taylor, line 49 \n";
	std::cerr << "**************************************************** \n";
	exit(1);
	}

       /* for all points determine the distance to all free surfaces */
	    
      for(i_index=0;i_index<number_primary_cells_i+2;i_index++)
      {
	  x1_coordinate_cell_center=(i_index-0.5)*mesh_width_x1;
	  
	  for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
	  {
	      x2_coordinate_cell_center=(j_index-0.5)*mesh_width_x2;
	      
	      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
	      {
		x3_coordinate_cell_center=(k_index-0.5)*mesh_width_x3;
		
		      surface_is_active   	=free_surfaces[0].active;
		      surface_orientation 	=free_surfaces[0].orientation;
		      surface_height		=free_surfaces[0].height;

		      if(surface_is_active)
		      {
			  switch(abs(surface_orientation))
			  {
			      case 1:
				  level_set[i_index][j_index][k_index]= 
						sign(1.0, (double) surface_orientation)* 
						  (x1_coordinate_cell_center-surface_height)
						  +0.05*cos(2*PI*x3_coordinate_cell_center);
				  break;
				  
			      case 2:
				  level_set[i_index][j_index][k_index]= 
						sign( 1.0, (double) surface_orientation)*
						  (x2_coordinate_cell_center-surface_height)
						  +0.05*cos(2*PI*x1_coordinate_cell_center);
				  break;
			    
			      case 3:
				  level_set[i_index][j_index][k_index]= 
						sign(1.0, (double) surface_orientation)*
						  (x3_coordinate_cell_center-surface_height)
						  +0.05*cos(2*PI*x1_coordinate_cell_center);
				  break;
			    
			      default:
		      
                                std::cerr << "**************************************************** \n";
				std::cerr << "ERROR \n";
				std::cerr << "invalid orientation of free surface, terminating...";
				std::cerr << "in function initialize_rayleigh_taylor, line 100 \n";
				std::cerr << "**************************************************** \n";
				exit(1);
			  }
		      }
	      }
	  }
      }
    }

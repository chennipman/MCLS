#include "../headers/array.h"
     double compute_mass_in_domain(
    	   Array3<double> volume_of_fluid,			// volume of fluid field
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	  double mesh_width_x3			// grid spacing in x3 direction (uniform)
	)
     {   
     	int i_index, j_index, k_index;  		// local variables for loop indexing
       double mass_in_domain=0;			// total mass in the domain
     	
      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	      {
		      mass_in_domain+=1.0-volume_of_fluid[i_index][j_index][k_index];
		      
	      }
	  }
      }
      mass_in_domain=mass_in_domain*mesh_width_x1*mesh_width_x2*mesh_width_x3;
      return mass_in_domain;
     }
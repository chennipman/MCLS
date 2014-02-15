#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/*  Function to initialize the volume of fluid field                            */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
      void initialize_volume_of_fluid(
	  Array3<double> level_set, 			// level set field 
	  Array3<double> volume_of_fluid, 			// volume of fluid field
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  double lower_bound_derivatives 		// lower bound for the first partial derivatives
							// to consider it a limiting case of vanishing
							// partial derivatives
	    )
      {
      /* function definitions */
      
      int    compute_volume_of_fluid(		// compute volume of fluid field
		Array3<double> level_set, 		// corresponding to a given 
		Array3<double> d_level_set_d_x1, 		// level-set field
		Array3<double> d_level_set_d_x2, 
		Array3<double> d_level_set_d_x3,
		Array3<double> volume_of_fluid,
		int number_primary_cells_i, 
		int number_primary_cells_j, 
		int number_primary_cells_k,
		double lower_bound_derivatives
		);
      void 	compute_level_set_gradient(		// compute gradient of level-set field
		Array3<double> level_set_star, 
		Array3<double> d_level_set_d_x1, 
		Array3<double> d_level_set_d_x2, 
		Array3<double> d_level_set_d_x3,
		int number_primary_cells_i, 
		int number_primary_cells_j, 
		int number_primary_cells_k
		);
      void field_neumann_boundary(			// apply neumann boundary condition to
	  	Array3<double> field, 			// cell centered field
	  	int number_primary_cells_i,	
	  	int number_primary_cells_j,	
	  	int number_primary_cells_k	
	  	);
      	void  field_extrapolate_boundary(      	// extrapolate field to virtual cells
            Array3<double> field, 			
            int number_primary_cells_i,	
            int number_primary_cells_j,	
            int number_primary_cells_k	
	    );

      Array3<double> d_level_set_d_x1;			// first partial derivative of
							// the level-set field wrt x1
							// second order central approximation
      Array3<double> d_level_set_d_x2;			// first partial derivative of 
							// the level-set field wrt x2
							// second order central approximation
      Array3<double> d_level_set_d_x3;			// first partial derivative of
 							// the level-set field wrt x3
							// second order central approximation



	/* allocate  memory for the gradient of the level-set field */
	
	d_level_set_d_x1.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x2.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	 
	
	/* compute the gradient of the level-set field, necessary for the level-set/vof conversion */
	
	compute_level_set_gradient(level_set, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );


	/* compute the volume of fluid field */
	
	if(compute_volume_of_fluid(level_set, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				volume_of_fluid,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				lower_bound_derivatives
				  ))
	{
	  std::cerr << "**************************************************** \n";
	  std::cerr<<"ERROR \n";
	  std::cerr<<"Something went wrong in the level-set to volume of fluid conversion \n";
	  std::cerr<<"Aborting...\n";
	  std::cerr<<"in function initialize_volume_of_fluid line 97 \n";
	  std::cerr << "**************************************************** \n";
	};
	
    /* extend the volume of fluid field to the virtual cells */
    
//     	field_neumann_boundary(volume_of_fluid, number_primary_cells_i, 
// 			   number_primary_cells_j,number_primary_cells_k);	
     	field_extrapolate_boundary(volume_of_fluid, number_primary_cells_i, 
			   number_primary_cells_j,number_primary_cells_k);	
	
	/* de-allocate  memory for the gradient of the level-set field */
	
	d_level_set_d_x1.destroy();
	d_level_set_d_x2.destroy();
	d_level_set_d_x3.destroy();
      }

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
EXPORT void initialize_volume_of_fluid(
	  Array3<double> level_set, 			// level set field 
	  Array3<double> volume_of_fluid, 		// volume of fluid field
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  double lower_bound_derivatives, 		// lower bound for the first partial derivatives
							// to consider it a limiting case of vanishing
							// partial derivatives
	  double volume_of_fluid_tolerance		// tolerance for volume of fluid value
	    )
      {
        Array3<double> d_level_set_d_x1;		// first partial derivative of
							// the level-set field wrt x1
							// second order central approximation
        Array3<double> d_level_set_d_x2;		// first partial derivative of 
							// the level-set field wrt x2
							// second order central approximation
        Array3<double> d_level_set_d_x3;		// first partial derivative of
 							// the level-set field wrt x3
							// second order central approximation
        Array3<double> invalid_vof_cells;               // indication field to show what is wrong
                                                        // indicator field showing cells that are either
                                                        // within bounds =0
                                                        // underfilled   =-1
                                                        // overfilled    =+1
                                                        // vapour cells  = 5
        int number_cells_vof_out_of_bounds;             // number of control volumes where the volume of fluid
                                                        // function is OUTSIDE the interval [0,1]
        int number_cells_numerical_vapor;               // number of control volumes where the volume of fluid
                                                        // function is INSIDE the interval [0,1]
                                                        // while the cell has 6 neighbours with the same sign:
                                                        // the cell is NOT an interface cell
        int number_cells_invalid_volume_of_fluid;       // sum of number of vapour cells and number of cells
                                                        // with the volume of fluid outside [0,1];



	/* allocate  memory for the gradient of the level-set field */
	/* and the visualization of invalid volume of fluid cells */
	
	d_level_set_d_x1.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x2.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	invalid_vof_cells.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);


	
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
	  std::cerr<<"in function initialize_volume_of_fluid line 86 \n";
	  std::cerr << "**************************************************** \n";
	  exit(1);
	};

	
	/* analyse the validity of the volume of fluid field */
     
        analyse_validity_vof_field( level_set, volume_of_fluid, invalid_vof_cells ,                               
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                    
                                          volume_of_fluid_tolerance, number_cells_vof_out_of_bounds, number_cells_numerical_vapor,                             
                                                number_cells_invalid_volume_of_fluid);
        if( number_cells_invalid_volume_of_fluid)
        {
 	        std::cerr << "**************************************************** \n";
	        std::cerr<<"ERROR \n";
	        std::cerr<<" The volume of fluid field is INVALID.\n";
        	std::cerr<<" total number of invalid cells: "<< number_cells_invalid_volume_of_fluid<<" \n";
        	std::cerr<<" number of vapor cells: "<< number_cells_numerical_vapor<<" \n";
        	std::cerr<<" number of cells with volume of fluid out of bounds: "<< number_cells_vof_out_of_bounds<<" \n";
        	std::cerr<<"in function initialize_volume_of_fluid line 106 \n";
        	std::cerr << "**************************************************** \n";
        	exit(1);
        }
 
    /* extend the volume of fluid field to the virtual cells */
    
     	field_extrapolate_boundary(volume_of_fluid, number_primary_cells_i, 
			   number_primary_cells_j,number_primary_cells_k);	
	
	/* de-allocate memory for the gradient of the level-set field */
	/* and the visualisation of the invalid volume of fluid cells */
	
	d_level_set_d_x1.destroy();
	d_level_set_d_x2.destroy();
	d_level_set_d_x3.destroy();
      	invalid_vof_cells.destroy();
      }

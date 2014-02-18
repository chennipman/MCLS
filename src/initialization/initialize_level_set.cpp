#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to initialize the level-set field                                  */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* A number of cases are predefined for the exact computation of the level-set  */
/* field.                           							*/
/* The cases include:                                                           */
/* - an ellipsoid, defined by its center and the three principle axes           */
/* - a symmetric capillary wave with elevation in the x3 direction              */
/*  in the initialization, a spherical bubble is recognized as one that has 3   */
/* identical principal axes								*/
/********************************************************************************/
EXPORT void initialize_level_set(
	  geometry flow_type,				// the kind of initial condition that has to be applied
	  bubble *the_bubbles,			// array with definitions of the bubbles
	  int number_of_bubbles,			// number of bubbles in the initial condition
	  surface *the_free_surfaces,			// array with definitions of the free surfaces
	  int number_of_free_surfaces, 		// number of bubbles in the domain (<10)
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	  double mesh_width_x3,			// grid spacing in x3 direction (uniform)
	  Array3<double> level_set				// level-set field
	  )
    {
   /* First set the level-set field to a large constant value, so we can minimize */
   /* with respect to this value						  */
   
     set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, level_set, 100000.0); 
      
      
    switch (flow_type) 
    {
      
	case bubbly_flow:
	  
	  /* initialize a  bubble or multiple bubbles if there are nonzero number */
	  /* of bubbles */
	  
	    if(number_of_bubbles)
	    {
		/* there is at least one bubble present */
		
		initialize_bubbles( the_bubbles, number_of_bubbles, 	                      
				number_primary_cells_i,	number_primary_cells_j,	 number_primary_cells_k,
				  mesh_width_x1, mesh_width_x2, mesh_width_x3, level_set);
	    }
	    else
	    {
		std::cerr << "**************************************************** \n";
		std::cerr << "WARNING \n";
		std::cerr << "You chose a flow with bubbles, but there are no bubbles defined. \n";
		std::cerr << "this means you are either using a single phase flow model, \n";
		std::cerr << "or considering surfaces only. \n";
		std::cerr << "In function initialize_level_set, line  101 \n";
		std::cerr << "**************************************************** \n";
	      
	    }
	    if(number_of_free_surfaces)
	    {
		/* there is at least one free surface present */
		
		initialize_free_surface( the_free_surfaces, number_of_free_surfaces,
				number_primary_cells_i,	number_primary_cells_j,	 number_primary_cells_k,
				  mesh_width_x1, mesh_width_x2, mesh_width_x3, level_set);
	    }
	    break;
	case wavy_flow:
	
	  /* initialize a capillary wave */
	  
	    std::cerr << "**************************************************** \n";
	    std::cerr << "ERROR \n";
	    std::cout << "capillary wave is not implemented in this version yet, terminating...\n";
	    std::cerr << "**************************************************** \n";
	   
	    exit(1);
	    break;
	  
	default:
	  
	    std::cerr << "**************************************************** \n";
	    std::cerr << "ERROR \n";
	    std::cout << "invalid geometry for initial condition, terminating...\n";
	    std::cerr << "**************************************************** \n";
	    exit(1);
	    break;
    }
    
    /* extend the level-set field to the virtual cells */
    
      
    field_extrapolate_boundary(level_set, number_primary_cells_i, 
			   number_primary_cells_j,number_primary_cells_k);	

}

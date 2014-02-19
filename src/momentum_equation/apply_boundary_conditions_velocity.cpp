#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to apply the boundary conditions to a given velocity field         */
/*               						                */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function replaces the function 'bound' from Sander's original           */
/* implementation. It sets the value of the velocity components on the          */
/* boundary and the value of the velocity components in the virtual cells.      */
/*                                                                              */
/********************************************************************************/
    
EXPORT void apply_boundary_conditions_velocity(
	  boundary_face boundary_faces[6],		// array with all the information
							// for the boundary conditions 
	  Array3<double> u_1_velocity, 			// velocity field x1 direction
	  Array3<double> u_2_velocity, 			// velocity field x2 direction
	  Array3<double> u_3_velocity, 			// velocity field x3 direction
	  double mesh_width_x1,				// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,				// grid spacing in x2 direction (uniform)
	  double mesh_width_x3,				// grid spacing in x3 direction (uniform)
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k			// number of primary (pressure) cells in x3 direction
     )
  {
     apply_boundary_conditions_velocity_u1(boundary_faces,		
					      u_1_velocity, 			
						  mesh_width_x1, mesh_width_x2, mesh_width_x3, 
						    number_primary_cells_i, 
						      number_primary_cells_j,			
							number_primary_cells_k);

     
     apply_boundary_conditions_velocity_u2(boundary_faces,		
					      u_2_velocity, 			
						  mesh_width_x1, mesh_width_x2, mesh_width_x3, 
						    number_primary_cells_i, 
						      number_primary_cells_j,			
							number_primary_cells_k);

     apply_boundary_conditions_velocity_u3(boundary_faces,		
					      u_3_velocity, 			
						  mesh_width_x1, mesh_width_x2, mesh_width_x3, 
						    number_primary_cells_i, 
						      number_primary_cells_j,			
							number_primary_cells_k);

 
  }







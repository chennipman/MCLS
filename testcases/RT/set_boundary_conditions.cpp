#include<math.h>
#include "../../src/headers/header_constants.h"
#include<iostream>

/********************************************************************************/
/********************************************************************************/
/*  Function to set the boundary conditions                                     */
/*  method. 										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* This is one of the major changes with respect to the code of Sander:         */
/* an attempt to more easily control the boundary conditions. 			*/
/* All the information is contained in an array boundary_faces.			*/
/* The idea is to use a default definition of the boundary conditions and only  */
/* set the only ones that are different from the default.				*/
/*  											*/
/*  The numbering of the faces is as follows:                                   */
/*  Face 0: n=( 1, 0, 0)                             				*/
/*  Face 1: n=(-1, 0, 0)                             				*/
/*  Face 2: n=( 0, 1, 0)                              				*/
/*  Face 3: n=( 0,-1, 0)                              				*/
/*  Face 4: n=( 0, 0, 1)                             				*/
/*  Face 5: n=( 0, 0,-1)                                              		*/
/*  											*/
/*  unknowns are                                                                */
/*  boundary variable 0: u1                                                    	*/
/*  boundary variable 1: u2                                                     */
/*  boundary variable 2: u3                                                     */
/*  boundary variable 3: level-set                                              */
/*  boundary variable 4: pressure                                               */
/*  											*/
/*  Default values are:                                                         */
/*  boundary_condition_type=neumann                                             */
/*  cell_centering=cell_centered                                               	*/
/*  boundary_condition_rule=constant                                           	*/
/*  boundary_condition_value=0                                                 	*/
/*  											*/
/*  											*/
/*  Note the idea of cell_centering is necessary because due to the staggering  */
/*  of the unknowns depending on which boundary is considered and which         */
/*  unknowns, the latter can be either in a cell centered or vertex centered    */
/*  position.										*/
/*  											*/
/*  											*/
/*  											*/
/*  											*/
/********************************************************************************/
EXPORT void set_boundary_conditions( 
	  boundary_face boundary_faces[6],			// array with all the information
								// for the boundary conditions 
	  vector initial_velocity
	  
	 )
    {
// Face 0: n=( 1, 0, 0) */
// Freeslip boundary condition:
// set the boundary condition for u_1 to homogeneous dirichlet= only non default */
//-----
// cell-centered: u2,u3, level-set, pressure
// boundary condition type: homogeneous neumann
//-----
// vertex-centered: u1
// boundary condition type: homogeneous dirichlet

//      boundary_faces[0].boundary_variables[0].boundary_condition_type=dirichlet;
    // u1
    boundary_faces[0].boundary_variables[0].boundary_condition_value=0.0;
    boundary_faces[0].boundary_variables[0].boundary_condition_type=dirichlet;
    // u2
//     boundary_faces[0].boundary_variables[1].boundary_condition_value=0.0;
//     boundary_faces[0].boundary_variables[1].boundary_condition_type=dirichlet;
//     boundary_faces[0].boundary_variables[1].cell_centering=cell_centered;
    // u3
//     boundary_faces[0].boundary_variables[2].boundary_condition_value=0.0;
//     boundary_faces[0].boundary_variables[2].boundary_condition_type=dirichlet;
//     boundary_faces[0].boundary_variables[2].cell_centering=cell_centered;

// Face 1: n=( -1, 0, 0) */
// Freeslip boundary condition:
// set the boundary condition for u_1 to homogeneous dirichlet= only non default */
//-----
// cell-centered: u2,u3, level-set, pressure
// boundary condition type: homogeneous neumann
//-----
// vertex-centered: u1
// boundary condition type: homogeneous dirichlet
    // u1
    boundary_faces[1].boundary_variables[0].boundary_condition_value=0.0;
    boundary_faces[1].boundary_variables[0].boundary_condition_type=dirichlet;
    // u2
//     boundary_faces[1].boundary_variables[1].boundary_condition_value=0.0;
//     boundary_faces[1].boundary_variables[1].boundary_condition_type=dirichlet;
//    boundary_faces[1].boundary_variables[1].cell_centering=cell_centered;
    // u3
//     boundary_faces[1].boundary_variables[2].boundary_condition_value=0.0;
//     boundary_faces[1].boundary_variables[2].boundary_condition_type=dirichlet;
//     boundary_faces[1].boundary_variables[2].cell_centering=cell_centered;


// Face 2: n=( 0, 1, 0) */
// Freeslip boundary condition:
// set the boundary condition for u_2 to homogeneous dirichlet= only non default */
//-----
// cell-centered: u1,u3, level-set, pressure
// boundary condition type: homogeneous neumann
//-----
// vertex-centered: u2
// boundary condition type: homogeneous dirichlet
    // u1
//     boundary_faces[2].boundary_variables[0].boundary_condition_value=0.0;
//     boundary_faces[2].boundary_variables[0].boundary_condition_type=dirichlet;
//     boundary_faces[2].boundary_variables[0].cell_centering=cell_centered;
    // u2
    boundary_faces[2].boundary_variables[1].boundary_condition_value=0.0;
    boundary_faces[2].boundary_variables[1].boundary_condition_type=dirichlet;
    // u3
//     boundary_faces[2].boundary_variables[2].boundary_condition_value=1.0;
//     boundary_faces[2].boundary_variables[2].boundary_condition_type=dirichlet;
//     boundary_faces[2].boundary_variables[2].cell_centering=cell_centered;

// Face 3: n=( 0, -1, 0) */
// Freeslip boundary condition:
// set the boundary condition for u_2 to homogeneous dirichlet= only non default */
//-----
// cell-centered: u1,u3, level-set, pressure
// boundary condition type: homogeneous neumann
//-----
// vertex-centered: u2
// boundary condition type: homogeneous dirichlet

    // u1
//     boundary_faces[3].boundary_variables[0].boundary_condition_value=0.0;
//     boundary_faces[3].boundary_variables[0].boundary_condition_type=dirichlet;
//     boundary_faces[3].boundary_variables[0].cell_centering=cell_centered;
    // u2
    boundary_faces[3].boundary_variables[1].boundary_condition_value=0.0;
    boundary_faces[3].boundary_variables[1].boundary_condition_type=dirichlet;
    // u3
//     boundary_faces[3].boundary_variables[2].boundary_condition_value=1.0;
//     boundary_faces[3].boundary_variables[2].boundary_condition_type=dirichlet;
//     boundary_faces[3].boundary_variables[2].cell_centering=cell_centered;


// Face 4: n=( 0, 0, 1) */
// Freeslip boundary condition:
// set the boundary condition for u_3 to homogeneous dirichlet= only non default */
//-----
// cell-centered: u1,u2, level-set, pressure
// boundary condition type: homogeneous neumann
//-----
// vertex-centered: u3
// boundary condition type: homogeneous dirichlet
    // u1
 //    boundary_faces[4].boundary_variables[0].boundary_condition_value=0.0;
 //    boundary_faces[4].boundary_variables[0].boundary_condition_type=dirichlet;
//     boundary_faces[4].boundary_variables[0].cell_centering=cell_centered;
    // u2
//     boundary_faces[4].boundary_variables[1].boundary_condition_value=0.0;
//     boundary_faces[4].boundary_variables[1].boundary_condition_type=dirichlet;
//     boundary_faces[4].boundary_variables[1].cell_centering=cell_centered;
    // u3
    boundary_faces[4].boundary_variables[2].boundary_condition_value=0.0;
    boundary_faces[4].boundary_variables[2].boundary_condition_type=dirichlet;
 

// Face 5: n=( 0, 0, -1) */
// Freeslip boundary condition:
// set the boundary condition for u_3 to homogeneous dirichlet= only non default */
//-----
// cell-centered: u1,u2, level-set, pressure
// boundary condition type: homogeneous neumann
//-----
// vertex-centered: u3
// boundary condition type: homogeneous dirichlet

     // u1
 //    boundary_faces[5].boundary_variables[0].boundary_condition_value=0.0;
 //    boundary_faces[5].boundary_variables[0].boundary_condition_type=dirichlet;
 //    boundary_faces[5].boundary_variables[0].cell_centering=cell_centered;
    // u2
//     boundary_faces[5].boundary_variables[1].boundary_condition_value=0.0;
//    boundary_faces[5].boundary_variables[1].boundary_condition_type=dirichlet;
//     boundary_faces[5].boundary_variables[1].cell_centering=cell_centered;
    // u3
    boundary_faces[5].boundary_variables[2].boundary_condition_value=0.0;
    boundary_faces[5].boundary_variables[2].boundary_condition_type=dirichlet;
   
}

#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to apply the boundary conditions to a given velocity field in the  */
/*  x1 direction.						                */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function replaces part of the function bound in the implementation of   */
/* Sander. Now the boundary conditions are applied on all components            */
/* separately. This function handles the boundary conditions for the x1         */
/* components.									*/
/********************************************************************************/
EXPORT void apply_boundary_conditions_velocity_u1(
	  boundary_face boundary_faces[6],		// array with all the information
							// for the boundary conditions 
	  Array3<double> u_1_velocity, 			// velocity field x1 direction
	  double mesh_width_x1,				// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,				// grid spacing in x2 direction (uniform)
	  double mesh_width_x3,				// grid spacing in x3 direction (uniform)
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  double actual_time				// actual time 
     )
     {
       double boundary_value;				// value of the boundary condition prescibed:
							// dirichlet-> function value
							// neumann -> normal derivative
       int face_index;					// index of the face under consideration
       int cell_label_boundary;				// index of the slice of cells adjacent to the
							// boundary
       int i_index, j_index, k_index;  			// local variables for loop indexing
       int increment_label_adjacent;			// increment of index for slice of cells on the
							// 'virtual' side of the face, with respect to
							// index of the slice of cells adjacent to the
							// boundary on the 'real' side.
 	double x,y,z; 					// spatial coordinates
	int boundary_var = 0; 
       /******************************************************************/
       /*   +/- I-index faces						 */
       /******************************************************************/
      
       for(face_index=0;face_index<=1; face_index++)
       {
	  if(face_index==0)
	  {
	      cell_label_boundary=number_primary_cells_i;
	      x = number_primary_cells_i*mesh_width_x1;
	  }
	  else
	  {
	      cell_label_boundary=0;
	      x = 0.0; 
	  }
	  if(boundary_faces[face_index].boundary_variables[boundary_var].boundary_condition_type==dirichlet)
	  {
		 for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
		 {
		      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
		      {
			y = mesh_width_x2*(j_index-0.5);
			z = mesh_width_x3*(k_index-0.5);			
			u_1_velocity[cell_label_boundary][j_index][k_index]=boundary_faces[face_index].boundary_variables[boundary_var].get_boundary_condition_value(actual_time,x,y,z);
		      }	  
  
		 }  
	  }  
	  else
	  {
	   std::cout<< "this is not implemented yet, check your input, aborting... \n";
	   std::cout<< "in apply_boundary_conditions_velocity_u1 \n";
	   exit(1);
	  }
  
       }  
       /******************************************************************/
       /*   +/- J-index faces						 */
       /******************************************************************/

       for(face_index=2;face_index<=3; face_index++)
       {
	  if(face_index==2)
	  {
	      cell_label_boundary=number_primary_cells_j+1;
	      y = number_primary_cells_j*mesh_width_x2;
	  }
	  else
	  {
	      cell_label_boundary=0;
	      y = 0; 
	  }
	  
	      /* increment_label_adjacent is the difference wrt to the */
	      /* label of the cells in the adjacent slice   */
	      /* =-1 for face_index=2, =+1 for face_index=3 */
	      
	  increment_label_adjacent= 2*face_index-5;
	  
	  if(boundary_faces[face_index].boundary_variables[0].boundary_condition_type==dirichlet)
	  {
		  /* DIRICHLET BOUNDARY CONDITION */
		  /* boundary_value is the prescribed DIRICHLET value */
		 for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
		 {
		      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
		      {
			x = mesh_width_x1*(i_index-0.5);
			z = mesh_width_x3*(k_index-0.5); 
			/* the boundary condition is applied using linear interpolation */
			/* U_virtual=2*U_boundary - U_real */
			boundary_value = boundary_faces[face_index].boundary_variables[boundary_var].get_boundary_condition_value(actual_time,x,y,z);
			u_1_velocity[i_index][cell_label_boundary][k_index]=
			    2.0*boundary_value
			    -u_1_velocity[i_index][cell_label_boundary+increment_label_adjacent][k_index];
		      }	  
  
		 }  
	  }  
	  else
	  {
	      if(boundary_faces[face_index].boundary_variables[0].boundary_condition_type==neumann)
	      {
		  /* NEUMANN BOUNDARY CONDITION */
		  /* boundary_value is the prescribed NEUMANN value */
		  /* the sign of increment_label_adjacent is used to specify the sign of the */
		  /* outward pointing normal vector */
		  
		 boundary_value= (double) (increment_label_adjacent)*
			boundary_faces[face_index].boundary_variables[0].boundary_condition_value;
			
		 for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
		 {
		      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
		      {
			
			/* the boundary condition is applied using second order approximation */
			/* to the normal derivative */
			/* U_virtual=U_real + dx* dU/dn */
			
			u_1_velocity[i_index][cell_label_boundary][k_index]=
			    mesh_width_x2*boundary_value
			    +u_1_velocity[i_index][cell_label_boundary+increment_label_adjacent][k_index];
		      }	  
  
		 }  
	      }
	      else
	      {
	    
		  std::cout<< "this is not implemented yet, check your input, aborting... \n";
		  std::cout<< "in apply_boundary_conditions_velocity_u1 \n";
		  exit(1);
	      }
	  }
  
       }  

      
       /******************************************************************/
       /*   +/- K-index faces						 */
       /******************************************************************/

       for(face_index=4;face_index<=5; face_index++)
       {
	  if(face_index==4)
	  {
	      cell_label_boundary=number_primary_cells_k+1;
	      z = mesh_width_x3*number_primary_cells_k;
	  }
	  else
	  {
	      cell_label_boundary=0;
	      z = 0.0; 
	  }
	  
	      /* increment_label_adjacent is the difference wrt to the */
	      /* label of the cells in the adjacent slice   */
	      /* =-1 for face_index=4, =+1 for face_index=5 */
	      
	      increment_label_adjacent= 2*face_index-9;
	      
	  if(boundary_faces[face_index].boundary_variables[0].boundary_condition_type==dirichlet)
	  {
		/* DIRICHLET BOUNDARY CONDITION */
		/* boundary_value is the prescribed DIRICHLET value */
		/* U_virtual=2*U_boundary - U_real */
		 for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
		 {
		      for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
		      {
			x = mesh_width_x1*(i_index-0.5);
			y = mesh_width_x2*(j_index-0.5); 
			boundary_value=boundary_faces[face_index].boundary_variables[boundary_var].get_boundary_condition_value(actual_time,x,y,z);
			
			/* the boundary condition is applied using linear interpolation */
			
			u_1_velocity[i_index][j_index][cell_label_boundary]=
			    2.0*boundary_value
			    -u_1_velocity[i_index][j_index][cell_label_boundary+increment_label_adjacent];
		      }	  
  
		 }  
	  }  
	  else
	  {
	      if(boundary_faces[face_index].boundary_variables[0].boundary_condition_type==neumann)
	      {
		  /* NEUMANN BOUNDARY CONDITION */
		  /* boundary_value is the prescribed NEUMANN value */
		  /* the sign of increment_label_adjacent is used to specify the sign of the */
		  /* outward pointing normal vector */
		  
		 boundary_value=(double) (increment_label_adjacent)*
			boundary_faces[face_index].boundary_variables[0].boundary_condition_value;
			
		 for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
		 {
		      for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
		      {
			
			/* the boundary condition is applied using second order approximation */
			/* to the normal derivative */
			/* U_virtual=U_real + dx* dU/dn */
			
			u_1_velocity[i_index][j_index][cell_label_boundary]=
			    mesh_width_x2*boundary_value
			    +u_1_velocity[i_index][j_index][cell_label_boundary+increment_label_adjacent];
		      }	  
  
		 }  
	      }  
	      else
	      {
	    
		  std::cout<< "this is not implemented yet, check your input, aborting... \n";
		  std::cout<< "in apply_boundary_conditions_velocity_u1 \n";
		  exit(1);
	      }
	  }
  
       }  
   }

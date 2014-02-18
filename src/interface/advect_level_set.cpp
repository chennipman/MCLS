#include "../headers/array.h"
#include<cstdlib>
#include<iostream>

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the advection of the level_set function          	*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The level-set field is advected using a first-order upwind scheme in space   */
/* and an explicit Euler time discretisation 					*/
/* For each cell face the computed flux is used for both neigbouring cells 	*/
/* Because the level_set field and the velocity field are staggered in time     */
/* we can use the velocity field at the new time level to define the advection	*/
/********************************************************************************/
//
EXPORT void advect_level_set(
      Array3<double> level_set_old, 		        // level set field at old time level
      Array3<double> level_set_star, 		// level set field at star time level
      Array3<double> u_1_velocity_new, 	        // velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new, 	        // velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,	        // velocity field at new time level x3 direction
      int number_primary_cells_i,	        // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	        // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	        // number of primary (pressure) cells in x3 direction
      double actual_time_step_level_set,        // time step used for level-set advection
					        // computed from all stability restrictions and 
					        // possibly subscycling
      double mesh_width_x1,		        // grid spacing in x1 direction (uniform)
      double mesh_width_x2,		        // grid spacing in x2 direction (uniform)
      double mesh_width_x3		        // grid spacing in x3 direction (uniform)
      
)
{
    double one_over_dx1	=    		        // 1/(grid spacing in x1 direction)
	1.0/(mesh_width_x1);
    double one_over_dx2	=    		        // 1/(grid spacing in x2 direction)
	1.0/(mesh_width_x1);
    double one_over_dx3	=    		        // 1/(grid spacing in x3 direction)
	1.0/(mesh_width_x3);
    double cell_face_velocity; 		        // velocity at the cell face (no interpolation required)
					        // because of staggered mesh
    double flux_level_set; 		        // advective flux
    int i_index, j_index, k_index;              // local variables for loop indexing
    

/* handle all the convective transport for internal faces in x1 direction */
/* Note: first real cell has index 1, last real cell index number_primary_cells_n */

  for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
  {
      for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
      {
	  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	  {
     
	      cell_face_velocity = u_1_velocity_new[i_index][j_index][k_index];
	      flux_level_set = compute_level_set_flux(cell_face_velocity,level_set_old[i_index][j_index][k_index],
					    level_set_old[i_index+1][j_index][k_index]);
	      flux_level_set *= actual_time_step_level_set*one_over_dx1;
	      level_set_star[i_index  ][j_index][k_index]  -=  flux_level_set; 
	      level_set_star[i_index+1][j_index][k_index]  +=  flux_level_set;
  
	  }  
  
      }  
     
  } 
  
/* handle all the convective transport for internal faces in x2 direction */
/* Note: first real cell has index 1, last real cell index number_primary_cells_n */

  for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
  {
      for(j_index=0;j_index<number_primary_cells_j+1;j_index++)
      {
	  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	  {
      
 	      cell_face_velocity = u_2_velocity_new[i_index][j_index][k_index];
	      flux_level_set = compute_level_set_flux(cell_face_velocity,level_set_old[i_index][j_index][k_index],
					    level_set_old[i_index][j_index+1][k_index]);
	      flux_level_set *= actual_time_step_level_set*one_over_dx2;
	      level_set_star[i_index][j_index  ][k_index]  -=  flux_level_set; 
	      level_set_star[i_index][j_index+1][k_index]  +=  flux_level_set;
 
	  }  
   
      }  
    
  }  

/* handle all the convective transport for internal faces in x3 direction */
/* Note: first real cell has index 1, last real cell index number_primary_cells_n */

  for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
  {
      for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
      {
	  for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
	  {
      
  	      cell_face_velocity = u_3_velocity_new[i_index][j_index][k_index];
	      flux_level_set = compute_level_set_flux(cell_face_velocity,level_set_old[i_index][j_index][k_index],
					    level_set_old[i_index][j_index][k_index+1]);
	      flux_level_set *= actual_time_step_level_set*one_over_dx3;
	      level_set_star[i_index][j_index][k_index  ]  -=  flux_level_set; 
	      level_set_star[i_index][j_index][k_index+1]  +=  flux_level_set;

  
	  }  
   
       }  
  }  
  
    /* extend the level-set field to the virtual cells */
    
     field_extrapolate_boundary(level_set_star, number_primary_cells_i, 
			   number_primary_cells_j,number_primary_cells_k);	
   
}




#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to apply inhomogeneous neumann boundary condition			*/
/*			for the pressure       d       				*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function will apply an inhomogeneous boundary condition on the pressure */
/* to the virtual cells that are direct neighbours of the nonvirtual cells      */
/* this is necessary because the pressure in these virtual cells is used to     */
/* compute the source term in the momentum equation for all velocity points     */
/* even those that are prescribed.                                              */
/********************************************************************************/
//
EXPORT void apply_boundary_conditions_pressure(
      Array3<double> pressure,                            // pressure field
      Array3<double> pressure_boundary_condition_x1,      // inhomogeneous boundary condition for
                                                          // the pressure planes with normal in x1 direction
      Array3<double> pressure_boundary_condition_x2,      // inhomogeneous boundary condition for the pressure
                                                          // the pressure planes with normal in x1 direction
      Array3<double> pressure_boundary_condition_x3,      // inhomogeneous boundary condition for the pressure
                                                          // the pressure planes with normal in x1 direction
      Array3<double> scaled_density_u1,                   // scaled density for the controlvolumes
                                                          // of the momentum equation in x1 direction
      Array3<double> scaled_density_u2,                   // scaled density for the controlvolumes
                                                          // of the momentum equation in x2 direction
      Array3<double> scaled_density_u3,                   // scaled density for the controlvolumes
                                                          // of the momentum equation in x3 direction
      double mesh_width_x1,                               // grid spacing in x1 direction (uniform)
      double mesh_width_x2,                               // grid spacing in x2 direction (uniform)
      double mesh_width_x3,                               // grid spacing in x3 direction (uniform)
      int number_primary_cells_i,                         // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,                         // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k                          // number of primary (pressure) cells in x3 direction
      
)

{
     int table_face_extrapolation_range[6][6]=
     { 
	{			1,			 1,			  1,  number_primary_cells_j, 			    1,  number_primary_cells_k},
	{  number_primary_cells_i,  number_primary_cells_i,			  1,  number_primary_cells_j,			    1,	number_primary_cells_k},
	{			1,  number_primary_cells_i,			  1,			   1,                       1,  number_primary_cells_k},
	{			1,  number_primary_cells_i,  number_primary_cells_j,  number_primary_cells_j,			    1,	number_primary_cells_k},
	{			1,  number_primary_cells_i,			  1,  number_primary_cells_j,			    1,			     1},
	{			1,  number_primary_cells_i,			  1,  number_primary_cells_j,  number_primary_cells_k,  number_primary_cells_k}
	
     };

     /* Shifts for cells to extrapolate along faces */
     
     /* edge 1   i_shift=-1, j_shift= 0, k_shift= 0  */
     /* edge 2   i_shift= 1, j_shift= 0, k_shift= 0  */
     /* edge 3   i_shift= 0, j_shift=-1, k_shift= 0  */
     /* edge 4   i_shift= 0, j_shift= 1, k_shift= 0  */
     /* edge 5   i_shift= 0, j_shift= 0, k_shift= -1 */
     /* edge 6   i_shift= 0, j_shift= 0, k_shift=  1 */

     int table_face_extrapolation_shift[6][3]=
     { 
	{-1, 0, 0 },
	{ 1, 0, 0 },
	{ 0,-1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, -1},
	{ 0, 0,  1}
     };
     
   int i_index, j_index, k_index;  	// local variables for loop indexing
   int face_index;	                // local variables for loop indexing
   int i_shift, j_shift, k_shift;  	// shifts underlying the interpolation
   int i_index_start; 	 	        // first i_index for extrapolation
   int i_index_end; 	 		// last i_index for extrapolation
   int j_index_start; 	 	        // first j_index for extrapolation
   int j_index_end; 	 		// last j_index for extrapolation
   int k_index_start; 	 	        // first k_index for extrapolation
   int k_index_end; 	 		// last k_index for extrapolation
   int  boundary_face_index;            // index of the face on the boundary
    
  
  
  
  
  

/* Apply an inhomogeneous Neumann boundary condition to the pressure. */
/* This means the lower index virtual cell will get the value of the first real cell */
/* and the higher index virtual cell has the value of the last real cell */


      /* handle the six faces */

    
      for(face_index=0;face_index<1+1;face_index++){

	  i_index_start	        = table_face_extrapolation_range[face_index][0];
	  i_index_end		= table_face_extrapolation_range[face_index][1];
	  j_index_start	        = table_face_extrapolation_range[face_index][2];
	  j_index_end		= table_face_extrapolation_range[face_index][3];
	  k_index_start	        = table_face_extrapolation_range[face_index][4];
	  k_index_end		= table_face_extrapolation_range[face_index][5];
	  i_shift		= table_face_extrapolation_shift[face_index][0];
	  j_shift		= table_face_extrapolation_shift[face_index][1];
	  k_shift		= table_face_extrapolation_shift[face_index][2];
          
          if(face_index==0)
          {
                   boundary_face_index=0;  
          }
          else
          {
                   boundary_face_index=number_primary_cells_i;  
          }
                 

	  for(i_index=i_index_start;i_index<i_index_end+1;i_index++)
	  {
	      for(j_index=j_index_start;j_index<j_index_end+1;j_index++)
	      {
		  for(k_index=k_index_start;k_index<k_index_end+1;k_index++)
		  {

                       pressure[i_index+i_shift][j_index+j_shift][k_index+k_shift] =
				      pressure[i_index][j_index][k_index]+mesh_width_x1*
				         scaled_density_u1[boundary_face_index][j_index][k_index]*
				              pressure_boundary_condition_x1[face_index][j_index][k_index];

		  }  
	      }  
	  }
      }


      for(face_index=2;face_index<3+1;face_index++){

          i_index_start      = table_face_extrapolation_range[face_index][0];
          i_index_end        = table_face_extrapolation_range[face_index][1];
          j_index_start      = table_face_extrapolation_range[face_index][2];
          j_index_end        = table_face_extrapolation_range[face_index][3];
          k_index_start      = table_face_extrapolation_range[face_index][4];
          k_index_end        = table_face_extrapolation_range[face_index][5];
          i_shift            = table_face_extrapolation_shift[face_index][0];
          j_shift            = table_face_extrapolation_shift[face_index][1];
          k_shift            = table_face_extrapolation_shift[face_index][2];
          
          if(face_index==2)
          {
                   boundary_face_index=0;  
          }
          else
          {
                   boundary_face_index=number_primary_cells_j;  
          }
         
         for(i_index=i_index_start;i_index<i_index_end+1;i_index++)
         {
             for(j_index=j_index_start;j_index<j_index_end+1;j_index++)
             {
                for(k_index=k_index_start;k_index<k_index_end+1;k_index++)
                {
                       pressure[i_index+i_shift][j_index+j_shift][k_index+k_shift] =
                                  pressure[i_index][j_index][k_index]+mesh_width_x2*
                                         scaled_density_u2[i_index][boundary_face_index][k_index]*
                                                 pressure_boundary_condition_x2[face_index-2][i_index][k_index];

                }
             }
         }
      }
      for(face_index=4;face_index<5+1;face_index++){

          i_index_start      = table_face_extrapolation_range[face_index][0];
          i_index_end        = table_face_extrapolation_range[face_index][1];
          j_index_start      = table_face_extrapolation_range[face_index][2];
          j_index_end        = table_face_extrapolation_range[face_index][3];
          k_index_start      = table_face_extrapolation_range[face_index][4];
          k_index_end        = table_face_extrapolation_range[face_index][5];
          i_shift            = table_face_extrapolation_shift[face_index][0];
          j_shift            = table_face_extrapolation_shift[face_index][1];
          k_shift            = table_face_extrapolation_shift[face_index][2];
          
          if(face_index==4)
          {
                   boundary_face_index=0;  
          }
          else
          {
                   boundary_face_index=number_primary_cells_k;  
          }


         for(i_index=i_index_start;i_index<i_index_end+1;i_index++)
         {
             for(j_index=j_index_start;j_index<j_index_end+1;j_index++)
             {
                for(k_index=k_index_start;k_index<k_index_end+1;k_index++)
                {

                    pressure[i_index+i_shift][j_index+j_shift][k_index+k_shift] =
                                  pressure[i_index][j_index][k_index]+mesh_width_x3*
                                         scaled_density_u3[i_index][j_index][boundary_face_index]*
                                                pressure_boundary_condition_x3[face_index-4][i_index][j_index];

                }
             }
         }
      }  
      
}

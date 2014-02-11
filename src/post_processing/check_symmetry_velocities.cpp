#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to check the symmetry of all quantities in the case of a           */
/*  symmetric problem 								*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
    void check_symmetry_velocities(
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double ***field_x1,				// field of which symmetry must be checked
      double ***field_x2,				// field of which symmetry must be checked
      double ***field_x3				// field of which symmetry must be checked
	  )
    {
    
      int i_index, j_index, k_index;  		// local variables for loop indexing
      int symmetry_plane_index=20;			// index of the symmetry plane 
     
      double max_error_jk_11=0;				// max difference between i and j constant planes
      double max_error_ik_22=0;				// max difference between i and k constant planes
      double max_error_ij_33=0;				// max difference between j and j constant planes
      double max_error_ij_21=0;
      double max_error_ik_31=0;
      double max_error_jk_32=0;
      
      

      
      	   for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
          {
		for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
           	{
			max_error_ij_21=std::max(fabs(
			field_x2[symmetry_plane_index][i_index][j_index]-
				field_x1[i_index][symmetry_plane_index][j_index]),
			    		max_error_ij_21);
			max_error_ij_33=std::max(fabs(
			field_x3[j_index][symmetry_plane_index][i_index]-
				field_x3[symmetry_plane_index][j_index][i_index]),
			    		max_error_ij_33);
			
			max_error_ik_31=std::max(fabs(
			field_x3[symmetry_plane_index][j_index][i_index]-
				field_x1[i_index][j_index][symmetry_plane_index]),
			    		max_error_ik_31);
			max_error_ik_22=std::max(fabs(
			field_x2[symmetry_plane_index][i_index][j_index]-
				field_x2[j_index][i_index][symmetry_plane_index]),
			    		max_error_ik_22);
			
			
			max_error_jk_32=std::max(fabs(
			field_x3[j_index][symmetry_plane_index][i_index]-
				field_x2[j_index][i_index][symmetry_plane_index]),
			    		max_error_jk_32);
			max_error_jk_11=std::max(fabs(
			field_x1[i_index][symmetry_plane_index][j_index]-
				field_x1[i_index][j_index][symmetry_plane_index]),
			    		max_error_jk_11);
			
            	}
          }

          std::cerr<<"Assymmetries in the weighted curvature \n";
          
          std::cerr<<"max_error_jk_11 "<< max_error_jk_11 <<"\n";
          std::cerr<<"max_error_ik_22 "<< max_error_ik_22 <<"\n";
          std::cerr<<"max_error_ij_33 "<< max_error_ij_33 <<"\n";
          std::cerr<<"max_error_ij_21 "<< max_error_ij_21<<"\n";
          std::cerr<<"max_error_ik_31 "<< max_error_ik_31 <<"\n";
          std::cerr<<"max_error_jk_32 "<< max_error_jk_32 <<"\n";
    }     
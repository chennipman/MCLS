#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include<string>
#include<sstream>
#include<fstream>
using namespace std;

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
   void compute_weighted_curvature(
	double ***level_set,				// level set field
	double ***curvature,				// interface curvature
       double mesh_width_x1,				// grid spacing in x1 direction (uniform)
       double mesh_width_x2,				// grid spacing in x2 direction (uniform)
       double mesh_width_x3,				// grid spacing in x3 direction (uniform)
       int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
       int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
       int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	double smoothing_distance_factor		// the heaviside function is smoothed over
							// an interval of width 2*smoothing_distance
	   
)
{		 
		 
      double computed_derivative_heaviside_function(// evaluate approximated derivative of heavi-side		 
	      double level_set_left,		
	      double level_set_right,	  	
	      double mesh_width_x1,			
	      double mesh_width_x2,			
	      double mesh_width_x3,			
	      double smoothing_distance_factor		
	  );
      void  write_coordinates_vtk( 			// write coordinates in vtk format 
	  ofstream& output_stream, 		
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k,		
	  double mesh_width_x1,			
	  double mesh_width_x2,			
	  double mesh_width_x3			
       );
      void  write_vertex_centered_vector_field_vtk(      // write vertex centered vector field in vtk format
	ofstream& output_stream, 			
	string vector_name,				
	string look_up_table_name,			
	double ***cell_centered_vector_field_1, 	
	double ***cell_centered_vector_field_2, 	
	double ***cell_centered_vector_field_3, 	
	int number_primary_cells_i,			
	int number_primary_cells_j,			
	int number_primary_cells_k			
	    );
		 
      double ***double_Matrix2(			// allocate memory for a three-dimensional array of doubles
		int number_primary_cells_i,				
		int number_primary_cells_j, 		
		int number_primary_cells_k
		);
      void   free_double_Matrix2( 			// deallocate memory for a three
		double ***doubleMatrix2, 		// dimensional array of doubles
		int number_primary_cells_i,	
		int number_primary_cells_j
		);
      void check_symmetry_velocities(		// check symmetry of vector field
      		int number_primary_cells_i,			
      		int number_primary_cells_j,			
      		int number_primary_cells_k,			
      		double ***field_x1,				
      		double ***field_x2,				
      		double ***field_x3				
	  );
      void interpolate_velocity_u1_vertex(		// interpolate velocity u1 to cell vertices
	  double ***u_1_velocity_new, 			
	  double ***u_1_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void interpolate_velocity_u2_vertex(		// interpolate velocity u2 to cell vertices
	  double ***u_2_velocity_new, 			
	  double ***u_2_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void interpolate_velocity_u3_vertex(		// interpolate velocity u3 to cell vertices
	  double ***u_3_velocity_new, 			
	  double ***u_3_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
	double ***weighted_curvature_x1;		// weighted curvature at u1 points
	double ***weighted_curvature_x2;		// weighted curvature at u2 points
	double ***weighted_curvature_x3;		// weighted curvature at u3 points
	double ***weighted_curvature_x1_vertex;    // weighted curvature x1 interpolated to vertices
	double ***weighted_curvature_x2_vertex;    // weighted curvature x2 interpolated to vertices
	double ***weighted_curvature_x3_vertex;    // weighted curvature x3 interpolated to vertices
	
	int i_index, j_index, k_index;  		// local variables for loop indexing
       int total_number_primary_cells=		// total number of primary cells
	  number_primary_cells_i*
	    number_primary_cells_j*
	      number_primary_cells_k;
       int total_number_vertices=			// total number of vertices
	  (number_primary_cells_i+1)*
	    (number_primary_cells_j+1)*
	      (number_primary_cells_k+1);
      string vector_name;				// name of the vector field to be written 
      string look_up_table_name;			// name of the look-up table to be used

	
	/* allocate memory for the weighted curvature fields */
	
	weighted_curvature_x1=double_Matrix2(number_primary_cells_i+1,number_primary_cells_j+2, number_primary_cells_k+2);
	weighted_curvature_x2=double_Matrix2(number_primary_cells_i+2,number_primary_cells_j+1, number_primary_cells_k+2);
	weighted_curvature_x3=double_Matrix2(number_primary_cells_i+2,number_primary_cells_j+2, number_primary_cells_k+1);
      
 	    /* allocate memory for the velocity at the cell vertices */
      
	weighted_curvature_x1_vertex=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
	weighted_curvature_x2_vertex=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
	weighted_curvature_x3_vertex=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
     
	
     for(i_index=0;i_index<number_primary_cells_i+1;i_index++)
      {
	
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      
	      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		      
	      {
		/* weight the curvature with the derivative of the heaviside function */
		/* this basically makes it only 'active' near the interface */
		
		
		weighted_curvature_x1[i_index][j_index][k_index]=0.5*( curvature[i_index][j_index][k_index]+
					      curvature[i_index+1][j_index][k_index])*
						computed_derivative_heaviside_function(
						    level_set[i_index][j_index][k_index],
							 level_set[i_index+1][j_index][k_index],	  		
							    mesh_width_x1, mesh_width_x2, mesh_width_x3,			
							      smoothing_distance_factor);
		      
	      }
	  }
      }
      
     for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	
	  for(j_index=0;j_index<number_primary_cells_j+1;j_index++)
	  {
	      
	      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	      {
		/* weight the curvature with the derivative of the heaviside function */
		/* this basically makes it only 'active' near the interface */
		
		weighted_curvature_x2[i_index][j_index][k_index]=0.5*( curvature[i_index][j_index][k_index]+
					      curvature[i_index][j_index+1][k_index])*
						computed_derivative_heaviside_function(
						    level_set[i_index][j_index][k_index],
							 level_set[i_index][j_index+1][k_index],	  		
							    mesh_width_x1, mesh_width_x2, mesh_width_x3,			
							      smoothing_distance_factor);
						
		      
	      }
	  }
      }
      
     for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      
	      for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
	      {
		/* weight the curvature with the derivative of the heaviside function */
		/* this basically makes it only 'active' near the interface */
		
		weighted_curvature_x3[i_index][j_index][k_index]=0.5*( curvature[i_index][j_index][k_index]+
					      curvature[i_index][j_index][k_index+1])*
						computed_derivative_heaviside_function(
						    level_set[i_index][j_index][k_index],
							 level_set[i_index][j_index][k_index+1],	  		
							    mesh_width_x1, mesh_width_x2, mesh_width_x3,			
							      smoothing_distance_factor);
						
		      
	      }
	  }
      }
      
      /* check the symmetry of the weighted curvature */
	 
//       check_symmetry_velocities( number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
//       					weighted_curvature_x1, weighted_curvature_x2, weighted_curvature_x3);

	ofstream interface_surfaceforce( "interface_surfaceforce_c++.txt");
	if(!interface_surfaceforce)
	{
	   /* the contructor returned a 0-pointer :-( */
		cout << "Cannot open file.\n";
		exit(1);
	}

	  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	  {
	  	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  	{
      			for(i_index=1;i_index<number_primary_cells_i;i_index++)
      			{ 
			interface_surfaceforce.setf(ios::fixed, ios::floatfield);
			interface_surfaceforce.width(4);
			interface_surfaceforce<<i_index<<"    "<<j_index<<"    "<<k_index<<"    ";
			interface_surfaceforce.precision(8);
// 			interface_surfaceforce<<weighted_curvature_x1[i_index][j_index][k_index]<<"    ";
// 			interface_surfaceforce<<0.5*( curvature[i_index][j_index][k_index]+
// 					      curvature[i_index+1][j_index][k_index])<<"    ";
// 			interface_surfaceforce<<computed_derivative_heaviside_function(
// 						    level_set[i_index][j_index][k_index],
// 							 level_set[i_index+1][j_index][k_index],	  		
// 							    mesh_width_x1, mesh_width_x2, mesh_width_x3,			
// 							      smoothing_distance_factor)<<"    ";
// 			interface_surfaceforce<<level_set[i_index][j_index][k_index]<<"    ";
 			interface_surfaceforce<<curvature[i_index][j_index][k_index]<<"    ";
 			interface_surfaceforce<<curvature[i_index+1][j_index][k_index]<<"    ";
 			interface_surfaceforce<<0.5*(curvature[i_index][j_index][k_index]+
 							curvature[i_index+1][j_index][k_index])<<"    \n";
			}
		}
	  }
      
      /* interpolate the weighted curvature to the vertices */
      
      /* interpolate the velocity to the vertices of the primary cells */
    
       interpolate_velocity_u1_vertex( weighted_curvature_x1, weighted_curvature_x1_vertex,		
 			number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
       interpolate_velocity_u2_vertex( weighted_curvature_x2, weighted_curvature_x2_vertex,		
 			number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
       interpolate_velocity_u3_vertex( weighted_curvature_x3, weighted_curvature_x3_vertex,		
 			number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);

      /* dump the weighted curvature to a vtk file for visualisation */
      
	ofstream output_vtk ( "debug_weighted_curvature.vtk");
	if(!output_vtk)
	{
	    /* the contructor returned a 0-pointer :-( */
	    cout << "Cannot open file.\n";
	    exit(1);
	}
	
	/* write the header for the tecplot file */
	/* write the header for the tecplot file */
	
	output_vtk << "# vtk DataFile Version 3.0 \n";
	output_vtk << "Solution file MCLS \n";
	output_vtk << "ASCII \n";
	output_vtk << "DATASET RECTILINEAR_GRID\n";
	output_vtk << "DIMENSIONS "<< number_primary_cells_i+1 << " "<< number_primary_cells_j+1;
	output_vtk << " "<< number_primary_cells_k+1 << "\n";
	output_vtk.setf(ios::scientific);

	write_coordinates_vtk( output_vtk, 		
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
				  mesh_width_x1, mesh_width_x2, mesh_width_x3);
	
	output_vtk << "POINT_DATA " << total_number_vertices << "\n";
	
	/* write velocity field */
	
	vector_name="weighted_curvature";
	look_up_table_name="vel_tbl";
	
	write_vertex_centered_vector_field_vtk( output_vtk, vector_name, look_up_table_name, 
						  weighted_curvature_x1_vertex, weighted_curvature_x2_vertex, weighted_curvature_x3_vertex, 		
						    number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
	
	
	free_double_Matrix2(weighted_curvature_x1, number_primary_cells_i+1, number_primary_cells_j+2);
	free_double_Matrix2(weighted_curvature_x2, number_primary_cells_i+2, number_primary_cells_j+1);
	free_double_Matrix2(weighted_curvature_x3, number_primary_cells_i+2, number_primary_cells_j+2);
       free_double_Matrix2(weighted_curvature_x1_vertex, number_primary_cells_i+1, number_primary_cells_j+1);
       free_double_Matrix2(weighted_curvature_x2_vertex, number_primary_cells_i+1, number_primary_cells_j+1);
       free_double_Matrix2(weighted_curvature_x3_vertex, number_primary_cells_i+1, number_primary_cells_j+1);
}			      
			      
			      
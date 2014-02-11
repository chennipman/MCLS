#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

    void output_solution(
	  double ***level_set_new, 		// level set field at new time level
						// mass conserving
	  double ***volume_of_fluid, 		// volume of fluid field
	  double ***curvature, 		// interface curvature
	  double ***unsmoothed_curvature,	// interface curvature without smoothing
	  double ***u_1_velocity_new, 	// velocity field at new time level x1 direction
	  double ***u_2_velocity_new, 	// velocity field at new time level x2 direction
	  double ***u_3_velocity_new,		// velocity field at new time level x3 direction
	  double ***pressure,			// pressure field
	  int vtk_output,			// =1, write output in vtk format
						// =0, skip output in vtk format
	  int tecplot_output,			// =1, write output in tecplot format
						// =0, skip output in tecplot format
	  int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
	  double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	  double mesh_width_x3,			// grid spacing in x3 direction (uniform)
	  int index_of_output_file		// index of the output file
		    )
	      {  
      void  write_coordinates_tecplot(         // write coordinates of the cell vertices to file in
	  ofstream& output_stream, 	       // tecplot format 			
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k,		
	  double mesh_width_x1,			
	  double mesh_width_x2,			
	  double mesh_width_x3			
			    );
      void interpolate_velocity_u1_center(	// interpolate velocity u1 to cell centers
	  double ***u_1_velocity_new, 			
	  double ***u_1_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void interpolate_velocity_u2_center(	// interpolate velocity u2 to cell centers
	  double ***u_2_velocity_new, 			
	  double ***u_2_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void interpolate_velocity_u3_center(	// interpolate velocity u3 to cell centers
	  double ***u_3_velocity_new, 			
	  double ***u_3_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void interpolate_velocity_u1_vertex(	// interpolate velocity u1 to cell vertices
	  double ***u_1_velocity_new, 			
	  double ***u_1_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void interpolate_velocity_u2_vertex(	// interpolate velocity u2 to cell vertices
	  double ***u_2_velocity_new, 			
	  double ***u_2_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void interpolate_velocity_u3_vertex(	// interpolate velocity u3 to cell vertices
	  double ***u_3_velocity_new, 			
	  double ***u_3_velocity_center,		
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
			    );
      void  write_cell_centered_field_tecplot( 		// write cell centered field to file
	ofstream& output_tecplot, 			// in tecplot format
	double ***pressure, 
	int number_primary_cells_i,
	int number_primary_cells_j, 
	int number_primary_cells_k
			    );
      void  write_coordinates_vtk( 		// write coordinates in vtk format 
	  ofstream& output_stream, 		
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k,		
	  double mesh_width_x1,			
	  double mesh_width_x2,			
	  double mesh_width_x3			
       );
      void  write_cell_centered_field_vtk(     // write cell centered field in vtk format
	ofstream& output_stream, 			
	string scalar_name,				
	string look_up_table_name,				
	double ***cell_centered_field, 			
	int number_primary_cells_i,			
	int number_primary_cells_j,			
	int number_primary_cells_k			
	    );
      void  write_vertex_centered_field_vtk(       // write vertex centered field in vtk format
	ofstream& output_stream, 			
	string scalar_name,				
	string look_up_table_name,			
	double ***cell_centered_field, 			
	int number_primary_cells_i,			
	int number_primary_cells_j,			
	int number_primary_cells_k			
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
      void   free_double_Matrix2( 		// deallocate memory for a three-dimensional array of doubles
		double ***doubleMatrix2, 
		int number_primary_cells_i,	
		int number_primary_cells_j
	    );
      
      double ***u_1_velocity_center;		// velocity in cell center, x1 component
      double ***u_2_velocity_center;		// velocity in cell center, x2 component
      double ***u_3_velocity_center;		// velocity in cell center, x3 component
      double ***u_1_velocity_vertex;		// velocity in cell vertex, x1 component
      double ***u_2_velocity_vertex;		// velocity in cell vertex, x2 component
      double ***u_3_velocity_vertex;		// velocity in cell vertex, x3 component
      int i_index, j_index, k_index; 		// local variables for loop indexing
      int full_row;				// a full row of 8 numbers has been written to file
      int total_number_primary_cells=		// total number of primary cells
	  number_primary_cells_i*
	    number_primary_cells_j*
	      number_primary_cells_k;
     int total_number_vertices=			// total number of vertices
	  (number_primary_cells_i+1)*
	    (number_primary_cells_j+1)*
	      (number_primary_cells_k+1);
     
      string basic_filename="solution_file.";	// the common part of the names of all output files
      string filename_tecplot;			// the filename of the current output file for tecplot format
      string filename_vtk;			// the filename of the current output file for vtk format
      string scalar_name;			// name of the scalar field to be written 
      string vector_name;			// name of the vector field to be written 
      string look_up_table_name;		// name of the look-up table to be used
      
      

	    string solution_file_index=convertInt(index_of_output_file); 
    
	    filename_tecplot = basic_filename + solution_file_index + ".dat";
	    filename_vtk     = basic_filename + solution_file_index + ".vtk";
      
	    /* allocate memory for the velocity at the cell centers */
      
	    u_1_velocity_center=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
						  number_primary_cells_k+2);
	    u_2_velocity_center=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
						  number_primary_cells_k+2);
	    u_3_velocity_center=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
						  number_primary_cells_k+2);
      
 	    /* allocate memory for the velocity at the cell vertices */
      
	    u_1_velocity_vertex=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
	    u_2_velocity_vertex=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
	    u_3_velocity_vertex=double_Matrix2(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
     
	    /* interpolate the velocity to the cell centers of the primary cells */
	    
	    interpolate_velocity_u1_center( u_1_velocity_new, u_1_velocity_center,		
				number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
	    interpolate_velocity_u2_center( u_2_velocity_new, u_2_velocity_center,		
				number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
	    interpolate_velocity_u3_center( u_3_velocity_new, u_3_velocity_center,		
				number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
  
	    /* interpolate the velocity to the vertices of the primary cells */
	    
 	    interpolate_velocity_u1_vertex( u_1_velocity_new, u_1_velocity_vertex,		
 				number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
 	    interpolate_velocity_u2_vertex( u_2_velocity_new, u_2_velocity_vertex,		
 				number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
 	    interpolate_velocity_u3_vertex( u_3_velocity_new, u_3_velocity_vertex,		
 				number_primary_cells_i,	number_primary_cells_j,	number_primary_cells_k);
  

     

      
	    if(tecplot_output)
	    {
	
	    /* generate an output file in tecplot format */
	
	    /* open the output file */
      
		  ofstream output_tecplot( filename_tecplot.c_str());
		  if(!output_tecplot)
		  {
		      /* the contructor returned a 0-pointer :-( */
		      cout << "Cannot open file.\n";
		      exit(1);
		  }
	    
		  /* write the header for the tecplot file */
	    
		  output_tecplot << "VARIABLES = x1, x2, x3, u1_velocity_center, ";
		  output_tecplot << "u2_velocity_center, u3_velocity_center, volume_of_fluid,  level_set, pressure \n";
		  output_tecplot << "ZONE I=" << number_primary_cells_i+1 << ", J=" << number_primary_cells_j+1;
		  output_tecplot << ", K=" << number_primary_cells_k+1 << ", ";
		  output_tecplot << "DATAPACKING=BLOCK, VARLOCATION = ([4,5,6,7,8,9]=CELLCENTERED) \n";
		  /* for the floating point data, use scientific notation */
		  output_tecplot.setf(ios::scientific);

		  /* write coordinates of cell vertices to file */
	    
		  write_coordinates_tecplot(output_tecplot, number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
				      mesh_width_x1, mesh_width_x2, mesh_width_x3);
	    
		  /* write velocity x1 component */

		  write_cell_centered_field_tecplot(output_tecplot, u_1_velocity_center, number_primary_cells_i,
				       number_primary_cells_j, number_primary_cells_k);
 	      
		  /* write velocity x2 component */

		  write_cell_centered_field_tecplot(output_tecplot, u_2_velocity_center, number_primary_cells_i,
				       number_primary_cells_j, number_primary_cells_k);
 	      
		  /* write velocity x3 component */

		  write_cell_centered_field_tecplot(output_tecplot, u_3_velocity_center, number_primary_cells_i,
				       number_primary_cells_j, number_primary_cells_k);
 	      
		  /* write volume of fluid */

		  write_cell_centered_field_tecplot(output_tecplot, volume_of_fluid, number_primary_cells_i,
				       number_primary_cells_j, number_primary_cells_k);
 	      
		  /* write level-set field */
	    
		  write_cell_centered_field_tecplot(output_tecplot, level_set_new, number_primary_cells_i,
				       number_primary_cells_j, number_primary_cells_k);
 	      
		  /* write pressure field */

		  write_cell_centered_field_tecplot(output_tecplot, pressure, number_primary_cells_i,
				       number_primary_cells_j, number_primary_cells_k);
 	      
	    }
 	    if(vtk_output)
	    {
		  /* generate an output file in vtk format */

		  /* generate the header for the vtk file */
	
	    /* open the output file */
      
		  ofstream output_vtk ( filename_vtk.c_str());
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
		  
		  /* write all cell centered data to file */

		  output_vtk << "CELL_DATA " << total_number_primary_cells << "\n";
		  
		  
		  /* write the volume of fluid */
		  
		  scalar_name="volume_of_fluid";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, volume_of_fluid, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the level-set field */
		  
		  scalar_name="level_set";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, level_set_new, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the curvature field */
		  
		  scalar_name="curvature";
		  look_up_table_name="curv_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, curvature, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the smoothed curvature field */
		  
		  scalar_name="unsmoothed_curvature";
		  look_up_table_name="smcurv_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, unsmoothed_curvature, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the pressure field */
  
		  scalar_name="pressure";
		  look_up_table_name="pres_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, pressure, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write all vertex centered data to file */

		  output_vtk << "POINT_DATA " << total_number_vertices << "\n";
		  
		  /* write velocity field */
		  
		  vector_name="velocity";
		  look_up_table_name="vel_tbl";
	 
		  write_vertex_centered_vector_field_vtk( output_vtk, vector_name, look_up_table_name, 
							  u_1_velocity_vertex, u_2_velocity_vertex, u_3_velocity_vertex, 		
							    number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
		  
		  /* write the volume of fluid */
		  
		  scalar_name="volume_of_fluid_vertices";
		  look_up_table_name="vof_tbl";
		 
		  write_vertex_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, volume_of_fluid, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the level-set field */
		  
		  scalar_name="level_set_vertices";
		  look_up_table_name="lvst_tbl";
		  
		  write_vertex_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, level_set_new, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the pressure field */
  
		  scalar_name="pressure_vertices";
		  look_up_table_name="pres_tbl";
		  
		  write_vertex_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, pressure, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

	
      
      
	    free_double_Matrix2(u_1_velocity_center, number_primary_cells_i+2, number_primary_cells_j+2);
	    free_double_Matrix2(u_2_velocity_center, number_primary_cells_i+2, number_primary_cells_j+2);
	    free_double_Matrix2(u_3_velocity_center, number_primary_cells_i+2, number_primary_cells_j+2);
	    free_double_Matrix2(u_1_velocity_vertex, number_primary_cells_i+1, number_primary_cells_j+1);
	    free_double_Matrix2(u_2_velocity_vertex, number_primary_cells_i+1, number_primary_cells_j+1);
	    free_double_Matrix2(u_3_velocity_vertex, number_primary_cells_i+1, number_primary_cells_j+1);
	    
	    
     }

  }
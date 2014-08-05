#include "../headers/array.h"
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;
// EXPORT string convertInt(int number)
// {
//    stringstream ss;//create a stringstream
//    ss << number;//add number to the stream
//    return ss.str();//return a string with the contents of the stream
// }

EXPORT void output_solution(
	  Array3<double> level_set_new, 	// level set field at new time level
						// mass conserving
	  Array3<double> volume_of_fluid, 	// volume of fluid field
	  Array3<double> curvature, 		// interface curvature
	  Array3<double> unsmoothed_curvature,	// interface curvature without smoothing
	  Array3<double> u_1_velocity_new, 	// velocity field at new time level x1 direction
	  Array3<double> u_2_velocity_new, 	// velocity field at new time level x2 direction
	  Array3<double> u_3_velocity_new,	// velocity field at new time level x3 direction
	  Array3<double> pressure,		// pressure field
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
      Array3<double> u_1_velocity_center;	// velocity in cell center, x1 component
      Array3<double> u_2_velocity_center;	// velocity in cell center, x2 component
      Array3<double> u_3_velocity_center;	// velocity in cell center, x3 component
      Array3<double> u_1_velocity_vertex;	// velocity in cell vertex, x1 component
      Array3<double> u_2_velocity_vertex;	// velocity in cell vertex, x2 component
      Array3<double> u_3_velocity_vertex;	// velocity in cell vertex, x3 component
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
      string filename_pure;			// the filename of the current output file for pure format
      string scalar_name;			// name of the scalar field to be written 
      string vector_name;			// name of the vector field to be written 
      string look_up_table_name;		// name of the look-up table to be used
      double pure_output = 1;
      

	    string solution_file_index=convertInt(index_of_output_file); 
    
	    filename_tecplot = basic_filename + solution_file_index + ".dat";
	    filename_vtk     = basic_filename + solution_file_index + ".vtk";
 	    filename_pure    = basic_filename + solution_file_index + ".pure";

     
	    /* allocate memory for the velocity at the cell centers */
      
	    u_1_velocity_center.create(number_primary_cells_i+2, number_primary_cells_j+2, 
						  number_primary_cells_k+2);
	    u_2_velocity_center.create(number_primary_cells_i+2, number_primary_cells_j+2, 
						  number_primary_cells_k+2);
	    u_3_velocity_center.create(number_primary_cells_i+2, number_primary_cells_j+2, 
						  number_primary_cells_k+2);
      
 	    /* allocate memory for the velocity at the cell vertices */
      
	    u_1_velocity_vertex.create(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
	    u_2_velocity_vertex.create(number_primary_cells_i+1, number_primary_cells_j+1, 
						  number_primary_cells_k+1);
	    u_3_velocity_vertex.create(number_primary_cells_i+1, number_primary_cells_j+1, 
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
				       
 	      	  std::cout<< "Finished Tecplot output \n";
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
		  output_vtk.precision(16);
		  
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

	
      }
      
	    u_1_velocity_center.destroy();
	    u_2_velocity_center.destroy();
	    u_3_velocity_center.destroy();
	    u_1_velocity_vertex.destroy();
	    u_2_velocity_vertex.destroy();
	    u_3_velocity_vertex.destroy();
	    
 	    std::cout<< "Finished vtk_plot output \n";
     

 	    if(pure_output)
	    {
		  /* generate an output file in pure format */

	
	    /* open the output file */
      
		  ofstream output_pure ( filename_pure.c_str());
		  if(!output_pure)
		  {
		      /* the contructor returned a 0-pointer :-( */
		      cout << "Cannot open file.\n";
		      exit(1);
		  }
	    
		  /* write the header for the tecplot file */
		  /* write the header for the tecplot file */
	    

		  if(!output_pure)
		  {
		      /* the contructor returned a 0-pointer :-( */
		      cout << "Cannot open file.\n";
		      exit(1);
		  }
	    
		  /* write the header for the tecplot file */
		  /* write the header for the tecplot file */
	    
		  output_pure << "# Unaveraged Datafile \n";
		  output_pure << "Solution file MCLS \n";
		  output_pure << "ASCII \n";
		  output_pure << "DATASET RECTILINEAR_GRID\n";
		  output_pure << "DIMENSIONS " << "\n" << number_primary_cells_i+1;
		  output_pure  << " "<< number_primary_cells_j+1 << " "<< number_primary_cells_k+1 << "\n";
		  output_pure.setf(ios::scientific);
		  output_pure.precision(16);
		  
		  write_coordinates_vtk( output_pure, 		
					number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
					  mesh_width_x1, mesh_width_x2, mesh_width_x3); // the same function as for the vtk_output is used
		  
		  
		  /* write all face data to file */
		  /* the ghost velocity cell are also written */
		  
		  
		  output_pure << "FACE_DATA " "\n";
		  output_pure << "VELOCITY X " << (number_primary_cells_i+1)*(number_primary_cells_j+2)*(number_primary_cells_k+2)<<"\n";
		  write_face_field_pure( output_pure,u_1_velocity_new,
					 number_primary_cells_i+1, number_primary_cells_j+2,number_primary_cells_k+2);
		  output_pure << "END VELOCITY X " << (number_primary_cells_i+1)*(number_primary_cells_j+2)*(number_primary_cells_k+2)<<"\n";
		  
		  output_pure << "VELOCITY Y " << (number_primary_cells_i+2)*(number_primary_cells_j+1)*(number_primary_cells_k+2)<<"\n";
		  write_face_field_pure( output_pure,u_2_velocity_new,
					 number_primary_cells_i+2, number_primary_cells_j+1,number_primary_cells_k+2);
		  output_pure << "END VELOCITY Y " << (number_primary_cells_i+1)*(number_primary_cells_j+2)*(number_primary_cells_k+2)<<"\n";

		  output_pure << "VELOCITY Z " << (number_primary_cells_i+2)*(number_primary_cells_j+2)*(number_primary_cells_k+1)<<"\n";
		  write_face_field_pure( output_pure,u_3_velocity_new,
					 number_primary_cells_i+2, number_primary_cells_j+2,number_primary_cells_k+1);
		  output_pure << "END VELOCITY Z " << (number_primary_cells_i+1)*(number_primary_cells_j+2)*(number_primary_cells_k+2)<<"\n";

		  /* write all cell centered data to file */

		  output_pure << "CELL_DATA " << total_number_primary_cells << "\n";
		  
		  
		  /* write the volume of fluid */
		  
		  scalar_name="volume_of_fluid";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_pure, scalar_name, look_up_table_name, volume_of_fluid, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
		  output_pure << "END " << scalar_name << "\n";

		  /* write the level-set field */
		  
		  scalar_name="level_set";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_pure, scalar_name, look_up_table_name, level_set_new, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
		  output_pure << "END " << scalar_name << "\n";


		  /* write the curvature field */
		  
		  scalar_name="curvature";
		  look_up_table_name="curv_tbl";
		  
		  write_cell_centered_field_vtk( output_pure, scalar_name, look_up_table_name, curvature, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
		  output_pure << "END " << scalar_name << "\n";

		  /* write the smoothed curvature field */
		  
		  scalar_name="unsmoothed_curvature";
		  look_up_table_name="smcurv_tbl";
		  
		  write_cell_centered_field_vtk( output_pure, scalar_name, look_up_table_name, unsmoothed_curvature, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
		  output_pure << "END " << scalar_name << "\n";

		  /* write the pressure field */
  
		  scalar_name="pressure";
		  look_up_table_name="pres_tbl";
		  
		  write_cell_centered_field_vtk( output_pure, scalar_name, look_up_table_name, pressure, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
		  output_pure << "END " << scalar_name << "\n";


     }
  }

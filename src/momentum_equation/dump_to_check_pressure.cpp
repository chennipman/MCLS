#include "../headers/array.h"
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/*
string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}
*/

EXPORT void dump_to_check_pressure(
	  Array3<double> u_1_velocity_new, 	// velocity field at new time level x1 direction
	  Array3<double> u_2_velocity_new, 	// velocity field at new time level x2 direction
	  Array3<double> u_3_velocity_new,	// velocity field at new time level x3 direction
	  Array3<double> pressure,		// pressure field
	  int vtk_output,			// =1, write output in vtk format
						// =0, skip output in vtk format
	  int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
	  double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	  double mesh_width_x3,			// grid spacing in x3 direction (uniform)
	  int index_of_output_file,		// index of the output file
	  int first_or_second_call		// index of the output file
	  
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
     
      string basic_filename="wrong_filename.";	// the common part of the names of all output files

     
         if (first_or_second_call == 1) // before momentum correction
    {
      printf("Before momentum correction vtk file \n"); 
            basic_filename="pressure_file_before.";	// the common part of the names of all output files

    }
    
    else if (first_or_second_call == 2) // after momentum correction
    {
      printf("After momentum correction vtk file \n");
      basic_filename="pressure_file_after.";	// the common part of the names of all output files
    }
    
    else if (first_or_second_call == 3) // after momentum correction
    {
      printf("Velocity_field_u_star \n");
      basic_filename="velocity_field_u_star.";	// the common part of the names of all output files
    }

    else 
    {
      basic_filename="wrong_filename_second.";	// the common part of the names of all output files
    }


      
      string filename_tecplot;			// the filename of the current output file for tecplot format
      string filename_vtk;			// the filename of the current output file for vtk format
      string scalar_name;			// name of the scalar field to be written 
      string vector_name;			// name of the vector field to be written 
      string look_up_table_name;		// name of the look-up table to be used
      
 
   
	    stringstream ss;//create a stringstream
	    ss << index_of_output_file;//add number to the stream
	    string solution_file_index= ss.str(); // convertInt(index_of_output_file); 
	    ss.str( std::string() );
	    ss.clear();
	    filename_vtk     = basic_filename + solution_file_index + ".vtk";
      
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
		  
		  
		  /* write the pressure field */
  
		  scalar_name="pressure";
		  look_up_table_name="pres_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, pressure, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write all vertex centered data to file */

		  output_vtk << "POINT_DATA " << total_number_vertices << "\n";
		  
		  /* write velocity field */
		  
		  vector_name="Convection_Diffusion";
		  look_up_table_name="vel_tbl";
	 
		  write_vertex_centered_vector_field_vtk( output_vtk, vector_name, look_up_table_name, 
							  u_1_velocity_vertex, u_2_velocity_vertex, u_3_velocity_vertex, 		
							    number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
		  
		  /* write the pressure field */
  
		  scalar_name="pressure_vertices";
		  look_up_table_name="pres_tbl";
		  
		  write_vertex_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, pressure, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

	
      
      
	    u_1_velocity_center.destroy();
	    u_2_velocity_center.destroy();
	    u_3_velocity_center.destroy();
	    u_1_velocity_vertex.destroy();
	    u_2_velocity_vertex.destroy();
	    u_3_velocity_vertex.destroy();
	    
 	    std::cout<< "Finished vtk_plot output \n";
     }

  }

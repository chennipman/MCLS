
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/********************************************************************************/
/********************************************************************************/
/*  Function to the solution of the interface quantities to file for inspection */
/*  method. 										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* The level set field is advected, but does not remain mass conserving even    */
/* when it is advected in a conservative way. In this function the level-set    */
/* field is adapted to the volume of fluid field in an iterative way.           */
/* When the maximum number of correction steps has been applied, and the        */
/* level-set field is still not completely complying with the volume of fluid   */
/* field, the computation is terminated.	                                	*/
/* For debugging purposes all variables involved in the corrective process are  */
/* written to file for visual inspection.						*/
/********************************************************************************/
   void dump_curvature_for_debugging(
       	double ***d_level_set_d_x1,			// first partial derivative wrt x1 of level-set
       	double ***d_level_set_d_x2,			// first partial derivative wrt x2 of level-set
       	double ***d_level_set_d_x3,			// first partial derivative wrt x3 of level-set
       	double ***d_2_level_set_d_x1_2,		// pure second partial derivative wrt x1 of level-set
       	double ***d_2_level_set_d_x2_2,		// pure second partial derivative wrt x2 of level-set
       	double ***d_2_level_set_d_x3_2,		// pure second partial derivative wrt x3 of level-set
       	double ***d_2_level_set_d_x1_d_x2,		// mixed second partial derivative wrt x1 and x2 of level-set
       	double ***d_2_level_set_d_x1_d_x3,		// mixed second partial derivative wrt x1 and x3 of level-set
       	double ***d_2_level_set_d_x2_d_x3,		// mixed second partial derivative wrt x2 and x3 of level-set
 	double ***length_gradient,			// length of gradient vector at center of cell
	double ***curvature,				// interface curvature
	double ***curvature_error,			// error in the curvature for laplace testcase
	int number_primary_cells_i,			// number of primary (pressure) 
							// cells in x1 direction
	int number_primary_cells_j,			// number of primary (pressure) 
							// cells in x2 direction
	int number_primary_cells_k,			// number of primary (pressure) 
							// cells in x3 direction
	  double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	  double mesh_width_x3			// grid spacing in x3 direction (uniform)
	
	)
	{
      void  write_coordinates_vtk( 			// write coordinates in vtk format 
	  ofstream& output_stream, 		
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k,		
	  double mesh_width_x1,			
	  double mesh_width_x2,			
	  double mesh_width_x3			
       );
      void  write_cell_centered_field_vtk(     	// write cell centered field in vtk format
	  ofstream& output_stream, 			
	  string scalar_name,				
	  string look_up_table_name,				
	  double ***cell_centered_field, 			
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
	    );
      int total_number_primary_cells=		// total number of primary cells
	  number_primary_cells_i*
	    number_primary_cells_j*
	      number_primary_cells_k;
      int total_number_vertices=			// total number of vertices
	  (number_primary_cells_i+1)*
	    (number_primary_cells_j+1)*
	      (number_primary_cells_k+1);
	
       int i_index, j_index, k_index;  		// local variables for loop indexing
       	string scalar_name;				// name of the scalar field to be written 
      	string look_up_table_name;			// name of the look-up table to be used
	
       /* dump the solution to a text file for inspection */
     
		  ofstream interface_curvature( "interface_curvature_c++.txt");
		  if(!interface_curvature)
		  {
		      /* the contructor returned a 0-pointer :-( */
		      cout << "Cannot open file.\n";
		      exit(1);
		  }

	  for(k_index=2;k_index<number_primary_cells_k;k_index++)
	  {
	  	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  	{
      			for(i_index=2;i_index<number_primary_cells_i;i_index++)
      			{ 
			interface_curvature.setf(ios::fixed, ios::floatfield);
			interface_curvature.width(4);
			interface_curvature<<i_index<<"    "<<j_index<<"    "<<k_index<<"    ";
			interface_curvature.precision(8);
			interface_curvature<<d_level_set_d_x1[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_level_set_d_x2[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_level_set_d_x3[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_2_level_set_d_x1_2[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_2_level_set_d_x2_2[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_2_level_set_d_x3_2[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_2_level_set_d_x1_d_x2[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_2_level_set_d_x1_d_x3[i_index][j_index][k_index]<<"    ";
			interface_curvature<<d_2_level_set_d_x2_d_x3[i_index][j_index][k_index]<<"    ";
			interface_curvature<<length_gradient[i_index][j_index][k_index]<<"    ";
			interface_curvature<<curvature[i_index][j_index][k_index]<<"     \n";
			}
		}
	  }
      
      /* dump the solution to a vtk file for visualisation */
      
		  ofstream output_vtk ( "debug_curvature.vtk");
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
		  
		  
		  /* write the curvature*/
		  
		  scalar_name="curvature";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, curvature, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  
		  /* write the first partial derivative wrt x1*/
		  
		  scalar_name="curvature_error";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, curvature_error, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  
		  /* write the first partial derivative wrt x1*/
		  
		  scalar_name="d_level_set_d_x1";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_level_set_d_x1, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  
		  /* write the first partial derivative wrt x2*/
		  
		  scalar_name="d_level_set_d_x2";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_level_set_d_x2, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the first partial derivative wrt x3*/
		  
		  scalar_name="d_level_set_d_x3";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_level_set_d_x3, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the second pure partial derivative wrt x1*/
		  
		  scalar_name="d_2_level_set_d_x1_2";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_2_level_set_d_x1_2, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);


		  /* write the second pure partial derivative wrt x2*/
		  
		  scalar_name="d_2_level_set_d_x2_2";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_2_level_set_d_x2_2, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the second pure partial derivative wrt x3*/
		  
		  scalar_name="d_2_level_set_d_x3_2";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_2_level_set_d_x3_2, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the length of the gradient */
		  
		  scalar_name="length_gradient";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, length_gradient, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the level-set field correction*/
		  
		  scalar_name="d_2_level_set_d_x1_d_x3";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_2_level_set_d_x1_d_x3, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the level-set field correction*/
		  
		  scalar_name="d_2_level_set_d_x2_d_x3";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, d_2_level_set_d_x2_d_x3, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);


		}
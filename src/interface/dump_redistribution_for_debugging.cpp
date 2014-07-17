#include "../headers/array.h"
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
EXPORT std::string convertInt(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}
using namespace std;

/********************************************************************************/
/********************************************************************************/
/*  Function to dump the solution of the interface quantities to file for       */
/*  inspection. 								*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The level set field is advected, but does not remain mass conserving even    */
/* when it is advected in a conservative way. In this function the level-set    */
/* field is adapted to the volume of fluid field in an iterative way.           */
/* When the maximum number of correction steps has been applied, and the        */
/* level-set field is still not completely complying with the volume of fluid   */
/* field, the computation is terminated.	                                */
/* For debugging purposes all variables involved in the corrective process are  */
/* written to file for visual inspection.					*/
/********************************************************************************/
EXPORT void dump_redistribution_for_debugging(
	Array3<double> level_set_original,	        // level set field before the redistribution
	Array3<double> volume_of_fluid_original,	// volume of fluid field before the redistribution
	Array3<double> level_set_updated,	        // level set field after 1 redistribution sweep
	Array3<double> invalid_vof_cells,		// indicator field showing cells that are either
                                                        // within bounds =0
                                                        // underfilled   =-1
                                                        // overfilled    =+1
                                                        // vapour cells  = 5
	Array3<double> volume_of_fluid_error,           // deviation of the current vof value from
							// an acceptable value (vapour cell, underfilled,
                                                        // overfilled)
	int number_primary_cells_i,			// number of primary (pressure) 
							// cells in x1 direction
	int number_primary_cells_j,			// number of primary (pressure) 
							// cells in x2 direction
	int number_primary_cells_k,			// number of primary (pressure) 
							// cells in x3 direction
	double mesh_width_x1,			        // grid spacing in x1 direction (uniform)
	double mesh_width_x2,			        // grid spacing in x2 direction (uniform)
	double mesh_width_x3,			        // grid spacing in x3 direction (uniform)
        int index_redistribution_attempt	        // number of attempts to achieve a
						        // valid volume of fluid field
						        // through the redistribution algorithm
	
	)
	{
        int total_number_primary_cells=		        // total number of primary cells
	  number_primary_cells_i*
	    number_primary_cells_j*
	      number_primary_cells_k;
	
       	string scalar_name;						// name of the scalar field to be written 
      	string look_up_table_name;					// name of the look-up table to be used
        string basic_filename="debug_interface_redistribution.";	// the common part of the names of all output files
        string filename_vtk;						// the filename of the current output file for vtk format
	
      
 	    string solution_file_index=convertInt(index_redistribution_attempt); 
    
	    filename_vtk     = basic_filename + solution_file_index + ".vtk";
     
      /* dump the solution to a vtk file for visualisation */
      
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
		  
		  scalar_name="volume_of_fluid_original";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, volume_of_fluid_original, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

                  /* write the volume of fluid error: the deviation from valid values in the interval [0,1] */
		  
		  scalar_name="volume_of_fluid_error";
		  look_up_table_name="vof_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, volume_of_fluid_error, 		
						number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the original level-set field */
		  
		  scalar_name="level_set_original";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, level_set_original, 		
			 			number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		  /* write the level-set field obtained after the first update to the volume of fluid field */
		  
		  scalar_name="level_set_updated";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, level_set_updated, 		
						 number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);


		  /* write the field that shows the cells with an invalid volume of fluid value */
		  
		  scalar_name="invalid_vof_cells";
		  look_up_table_name="lvst_tbl";
		  
		  write_cell_centered_field_vtk( output_vtk, scalar_name, look_up_table_name, invalid_vof_cells, 		
			 			number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

		}

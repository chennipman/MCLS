
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
/*  Function to export a matrix to matlab for inspection                        */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The function creates an m-file that can be run in matlab and generates       */
/* the matrix from the internal sparse format.                                  */
/********************************************************************************/

  int export_matrix_matlab(
	int i_dimension,    		// number of unknowns in the system in i-direction
	int j_dimension, 		// number of unknowns in the system in i-direction
	int k_dimension, 		// number of unknowns in the system in i-direction
	int number_matrix_connections,  // number of connections in the matrix
	double **matrix_A,		// matrix under consideration
	double *rhside_vector,		// rhside vector under consideration	
	double *solution_vector,	// solution vector under consbgideration
	string variable_name		// variable name to which matrix corresponds
    
    
      )
  {
      string filename_end=".m";			// the common part of the names of all output files
      string filename_matlab;			// the filename of the current output file for tecplot format
      int one_dimensional_index;
      int dof_in_row;				// number of unknowns in the first dimension
						// of three dimensional grid
      int dof_in_slice;				// number of unknowns in the first two dimensions
						// of three dimensional grid "a slice"
      int dof_total;				// total number of unknonws in the 
						// three dimensional grid
      int connection_index;
    
      filename_matlab=variable_name+filename_end;
	
	ofstream output_matlab( filename_matlab.c_str());
	if(!output_matlab)
	{
	    /* the contructor returned a 0-pointer :-( */
	    std::cout << "Cannot open file....\n";
	    exit(1);
	}
	
	/* write the initial part for the m file */
	/* this information is necessary to check the matrix */
	/* because the sparsity pattern will change in the vicinity */
	/* of those dofs that are near the end of a row/slice */
	
		 
	dof_in_row   = i_dimension;
	dof_in_slice = i_dimension*j_dimension;
	dof_total    = i_dimension*j_dimension*k_dimension;
	
	output_matlab << "clear all ;\n";
	output_matlab << "dof_in_row		= "<< dof_in_row	<< ";\n";
	output_matlab << "dof_in_slice		= "<< dof_in_slice	<< ";\n";
	output_matlab << "dof_total		= "<< dof_total		<< ";\n";
	
	/* write all the coefficients of the matrix to the file */
	/* only one triangular part is  exported, because the   */
	/* matrix is symmetric					*/
       
	for(connection_index=0;connection_index<4;connection_index++)
	{
	    output_matlab << "d"<<connection_index+1<<" = [" <<"\n";
	    output_matlab.setf(ios::scientific);
	    for(one_dimensional_index=0;one_dimensional_index<dof_total;one_dimensional_index++)
	    {
	      output_matlab<<matrix_A[connection_index][one_dimensional_index]<<"\n";
	    }
	    output_matlab <<"]; \n";
	}
	
	/* write the right hand side vector to file */
	
	output_matlab << "rhside=[ \n";
	for(one_dimensional_index=0;one_dimensional_index<dof_total;one_dimensional_index++)
	{
	   output_matlab<<rhside_vector[one_dimensional_index]<<"\n";
	}
	output_matlab <<"]; \n";

	/* write the solution vector to file */
	
	output_matlab << "solution=[ \n";
	for(one_dimensional_index=0;one_dimensional_index<dof_total;one_dimensional_index++)
	{
	   output_matlab<<solution_vector[one_dimensional_index]<<"\n";
	}
	output_matlab <<"]; \n";

        output_matlab << "A=spdiags([d1, d2, d3, d4],";
	output_matlab << "[0,-1,"<< -dof_in_row<<","<<-dof_in_slice<<"],";
	output_matlab << dof_total<<","<<dof_total<<"); \n";
	output_matlab << "A=A+(A'-spdiags(diag(A),0,"<<dof_total<<","<<dof_total<< "));\n";
	
	output_matlab.close();
	return 0;
  }
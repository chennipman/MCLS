#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

#include <string>
#include <sstream>
#include <fstream>
using namespace std;
/********************************************************************************/
/********************************************************************************/
/*  Function to solve the system of linear equations for the pressure equation  */
/*  using preconditioned conjugate gradient method                              */
/*  										       */
/*  Programmer	: Duncan van der Heul       					       */
/*  Date	: 10-03-2013       						       */
/*  Update	:        							       */
/********************************************************************************/
/* Notes									       */
/* For the pressure, an inhomogeneous neumann boundary condition hold on the    */
/* boundaries where the velocity is zero. The normal derivative of the pressure */
/* balances with the gravitational acceleration.				       */
/********************************************************************************/
 void solve_pressure_correction_system(
      Array2<double> pressure_matrix, 		          // pressure matrix
      Array1<double> pressure_rhside,		 	   // pressure rhside
      Array3<double> pressure,			          // pressure field
      int number_primary_cells_i,	 	          // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	  	          // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,		          // number of primary (pressure) cells in x3 direction
      double   tolerance_pressure,	  	          // the tolerance with which the system is solved	
      int maximum_iterations_allowed_pressure         // maximum number of iterations allowed for the
						          //conjugate gradient method
      )
   {
     void build_preconditioner(			   // build incomplete choleski preconditioner
	int i_dimension,    	
	int j_dimension, 	
	int k_dimension, 	
	Array2<double> matrix_A,	
	Array1<double> preconditioner_matrix_M 	
	);
     int conjugate_gradient_method(		          // solve linear system with conjugate gradient
        int i_dimension,   			          // method (only SPD systems)
        int j_dimension,   		
        int k_dimension,   		
        Array2<double> matrix_A,   		
        Array1<double> preconditioner_matrix_M,  	
        Array1<double> rhside_vector_b,   		
        Array1<double> solution_vector_x, 		
        double   tolerance,	  		
        int     &iteration_number,  	
        double  &relative_L2_norm_residual, 
        double  &relative_Linfinity_norm_residual,
        int maximum_iterations_allowed	
        );
     void decompress_solution_pressure(                   //  remap one dimensional pressure array
	 Array3<double> pressure,		                  // to multidimensional pressure array
	 Array1<double> pressure_correction,  
	 int number_primary_cells_i,	     
	 int number_primary_cells_j,	     
	 int number_primary_cells_k	     
	 );
     void compress_solution_pressure(                     // map multidimensional pressure array
         Array3<double> pressure,                              // to one dimensional array
         Array1<double> compressed_pressure,         
         int number_primary_cells_i,            
         int number_primary_cells_j,            
         int number_primary_cells_k             
         );
     int export_matrix_matlab(           		   // export the matrix to matlab file
	 int i_dimension,    		
	 int j_dimension, 		
	 int k_dimension, 		
	 int number_matrix_connections,  
	 Array2<double> matrix_A,		
	 Array1<double> rhside_vector,		
	 Array1<double> solution_vector,	
	 string variable_name
	 );
     void  field_extrapolate_boundary(      	   // extrapolate field to virtual cells
        Array3<double> field, 			
        int number_primary_cells_i,	
        int number_primary_cells_j,	
        int number_primary_cells_k	
	);
     void field_neumann_boundary(                     // apply neumann boundary condition to
        Array3<double> field,                               // cell centered field
        int number_primary_cells_i,
        int number_primary_cells_j,
        int number_primary_cells_k
        );
     void set_constant_vector(      	                 // set 1-dimensional array to constant value
	 int vector_length,	
	 Array1<double> vector_to_set,	
	 double constant_value	
     );
     int project_pressure_rhside(			  // project rhside vector on column space
	 int total_number_pressure_points,		
	 Array1<double> pressure_rhside		 	
     );
     int shift_pressure_solution(                   // shift pressure solution
        int total_number_pressure_points,
        Array1<double> compressed_pressure               
     );
    
    

      int  iteration_number;  		  	   // the number of iterations where the iterative
						          // process was terminated
      double   relative_L2_norm_residual; 	          // the L2 norm of the residual, scaled with
						          // the L2 norm of the right hand side
      double   relative_Linfinity_norm_residual;      // the L infinity norm of the residual
						          // scaled with the maximum difference 
						          // between two components of the residual
      Array1<double> preconditioner_matrix_M;		   // preconditioner matrix (only main diagonal)
      Array1<double> compressed_pressure;		          // 1-D array with pressure , with virtual points excluded
      int total_number_pressure_points;		   // total number of points with pressure
      int i;
      
      /* allocate memory for the main diagonal of preconditioner matrix M */
      /* and for the compressed solution vector */
      
      total_number_pressure_points=number_primary_cells_i*number_primary_cells_j*number_primary_cells_k;
     
      preconditioner_matrix_M.create(total_number_pressure_points);
      compressed_pressure.create(total_number_pressure_points);
      
      /* build the preconditioner matrix M */
      
      build_preconditioner( number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
			    pressure_matrix, preconditioner_matrix_M);
      
      /* project rhside on column space of matrix */
      
      project_pressure_rhside( total_number_pressure_points, pressure_rhside);
      
      /* set starting value for pressure correction to zero */
      
      set_constant_vector(total_number_pressure_points, compressed_pressure, 0.0);

      /* set starting value of pressure to previous time step */
      
      compress_solution_pressure( pressure, compressed_pressure, 
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
 
     
      /* solve the linear system for the new pressure */

      if(conjugate_gradient_method( number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				  pressure_matrix, preconditioner_matrix_M, pressure_rhside, compressed_pressure,
				    tolerance_pressure, iteration_number, relative_L2_norm_residual,
				      relative_Linfinity_norm_residual, maximum_iterations_allowed_pressure))
      {
	std::cout << " No convergence in linear solver for pressure equation.\n";
	exit(1);
      }
      else
      {
	std::cout << " The pressure equation converged in " << iteration_number <<" iterations,\n";
	std::cout << " with relative L2 norm of residual " << relative_L2_norm_residual <<" \n";
	
       /* shift the pressure solution */

      shift_pressure_solution(total_number_pressure_points, compressed_pressure);
      
       /* update pressure with pressure correction */

       decompress_solution_pressure(pressure, compressed_pressure,
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
      }

        /* export the matrix to matlab for inspection if requested */

//       export_matrix_matlab(number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
//                       4, pressure_matrix, pressure_rhside, compressed_pressure, "pressure");

      
     
      /* deallocate memory for the main diagonal of preconditioner matrix M */
      
      
      preconditioner_matrix_M.destroy();
      compressed_pressure.destroy();
      
   }    

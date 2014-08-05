#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

#include <string>
#include <sstream>
#include <fstream>
#include "dpcg_rohit_paralution.hpp"
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
EXPORT void solve_pressure_correction_system(
      Array2<double> pressure_matrix, 		          // pressure matrix
      Array1<double> pressure_rhside,		 	  // pressure rhside
      Array3<double> pressure,			          // pressure field
      int number_primary_cells_i,	 	          // number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	  	          // number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,		          // number of primary (pressure) cells in x3 direction
      double   tolerance_pressure,	  	          // the tolerance with which the system is solved	
      int maximum_iterations_allowed_pressure,             // maximum number of iterations allowed for the
						          //conjugate gradient method
      Array3<double> level_set				// levelset for the entire domain with 1 extra point on all sides					          
      )
   {
      int  iteration_number;  		  	          // the number of iterations where the iterative
						          // process was terminated
      double   relative_L2_norm_residual; 	          // the L2 norm of the residual, scaled with
						          // the L2 norm of the right hand side
      double   relative_Linfinity_norm_residual;          // the L infinity norm of the residual
						          // scaled with the maximum difference 
						          // between two components of the residual
      Array1<double> preconditioner_matrix_M;		  // preconditioner matrix (only main diagonal)
      Array1<double> compressed_pressure;		  // 1-D array with pressure , with virtual points excluded
      int total_number_pressure_points;		          // total number of points with pressure
      struct timeval now;
      double tick, tack;
      
      
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

	gettimeofday(&now, NULL);
      tick = now.tv_sec*1000000.0+(now.tv_usec);
//       if(conjugate_gradient_method( number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
// 				  pressure_matrix, preconditioner_matrix_M, pressure_rhside, compressed_pressure,
// 				    tolerance_pressure, iteration_number, relative_L2_norm_residual,
// 				      relative_Linfinity_norm_residual, maximum_iterations_allowed_pressure))
      	
	if(call_to_cg_wrapper(pressure_matrix, compressed_pressure, pressure_rhside, maximum_iterations_allowed_pressure,
				tolerance_pressure, number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				&iteration_number, &relative_L2_norm_residual, level_set))
      {
	std::cout << " No convergence in linear solver for pressure equation.\n";
	exit(1);
      }
      else
      {
	std::cout << " The pressure equation converged in " << iteration_number <<" iterations,\n";
	std::cout << " with relative L2 norm of residual " << relative_L2_norm_residual <<" \n";
	gettimeofday(&now, NULL);
	tack = now.tv_sec*1000000.0+(now.tv_usec);
	std::cout << "Call to DPCG took in total:" << (tack-tick)/1000000 << " sec" << std::endl;
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

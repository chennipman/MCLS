#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<math.h>
   
/********************************************************************************/
/*  Function to compute the surface tension body force, in the CSF approach     */
/*            										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes									*/
/* The computation of the capillary time-step restriction was part of this      */
/* function in Sanders implementation. Now it is moved to a separate function   */
/* that evaluates all different restrictions.                                   */
/********************************************************************************/

EXPORT void compute_surface_tension_body_force(
 	    Array3<double> level_set, 				// level set field
	    Array3<double> surface_tension_body_force_x1,	// x1 component of the body force due to
								// CSF formulation of surface tension model
	    Array3<double> surface_tension_body_force_x2,	// x2 component of the body force due to
								// CSF formulation of surface tension model
	    Array3<double> surface_tension_body_force_x3,	// x3 component of the body force due to
								// CSF formulation of surface tension model
	    Array3<double> curvature,				// interface curvature
	    Array3<double> unsmoothed_curvature,		// interface curvature before smoothing
	    int number_primary_cells_i,			        // number of primary (pressure) cells in x1 direction
	    int number_primary_cells_j,			        // number of primary (pressure) cells in x2 direction
	    int number_primary_cells_k,			        // number of primary (pressure) cells in x3 direction
	    double mesh_width_x1,				// grid spacing in x1 direction (uniform)
	    double mesh_width_x2,				// grid spacing in x2 direction (uniform)
	    double mesh_width_x3,				// grid spacing in x3 direction (uniform)
	    double rho_plus_over_rho_minus, 			// rho_plus / rho_minus indicated 
	    double sigma_over_rho_minus,			// sigma / rho_minus (scaled sigma)
	    double maximum_weighted_curvature,		        // maximum 'active' value of the curvature 
								// used to evaluate the capillary time step
								// restriction 
	    int apply_curvature_smoothing,			// =1, apply curvature smoothing
								// =0, use unsmoothed curvature
	    int number_curvature_smoothing_steps,		// number of iterations applied in the
								// curvature smoothing algorithm
	    int apply_curvature_smoothing_filter,		// =1, apply curvature smoothing filter
	    							// =0, do not apply curvature smoothing filter		
	    double smoothing_distance_factor
	      )								
{
      double maximum_body_force_x1=0;					// maximum body force, x1 component
      double maximum_body_force_x2=0;					// maximum body force, x2 component
      double maximum_body_force_x3=0;					// maximum body force, x3 component
      
      double total_body_force_x1=0;					// total body force, x1 component
      double total_body_force_x2=0;					// total body force, x2 component
      double total_body_force_x3=0;					// total body force, x3 component
      

	std::cerr<<" computing surface tension body forces \n";
  
	/* compute the curvature */
	
	compute_curvature( level_set, unsmoothed_curvature, 
			    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
			      mesh_width_x1, mesh_width_x2, mesh_width_x3);
	
		
	/* smooth the curvature if necessary */
	
	if(apply_curvature_smoothing)
	{
		
	    /* apply smoothing algorithm */
	    
	    smooth_curvature(level_set, curvature, unsmoothed_curvature,
			      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				mesh_width_x1, mesh_width_x2, mesh_width_x3,
				  apply_curvature_smoothing_filter, number_curvature_smoothing_steps);

	}
	else
	{
	    /* don't bother with smoothing */
	    
      	    copy_cell_centered_field(unsmoothed_curvature,  curvature, 
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
			
	}


	maximum_weighted_curvature=0;
	
	/* compute the body force in x1 direction */
	
	compute_body_force_x1(curvature, level_set, surface_tension_body_force_x1, 
			      mesh_width_x1, mesh_width_x2, mesh_width_x3,
			      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
				rho_plus_over_rho_minus, sigma_over_rho_minus,	
				  maximum_body_force_x1, total_body_force_x1, maximum_weighted_curvature,
				  smoothing_distance_factor);						

	/* compute the body force in x2 direction */
	
	compute_body_force_x2(curvature, level_set, surface_tension_body_force_x2, 
			      mesh_width_x1, mesh_width_x2, mesh_width_x3,
			      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
				rho_plus_over_rho_minus, sigma_over_rho_minus,	
				  maximum_body_force_x2, total_body_force_x2, maximum_weighted_curvature,
				  smoothing_distance_factor);						 

   
	/* compute the body force in x3 direction */
	
	compute_body_force_x3(curvature, level_set, surface_tension_body_force_x3, 
			      mesh_width_x1, mesh_width_x2, mesh_width_x3,
			      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
				rho_plus_over_rho_minus, sigma_over_rho_minus,	
				  maximum_body_force_x3, total_body_force_x3, maximum_weighted_curvature,
				  smoothing_distance_factor);						

 
        
}

#include "../headers/array.h"
/********************************************************************************/
/*  Function to initialize the curvature                                        */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* To be able to use the standard function to write the solution to file for    */
/* postprocessing, the curvature and smoothed curvature need to be computed     */
/* to be added to the initial condition.                                        */
/********************************************************************************/
void initialize_curvature(						
      Array3<double> level_set, 				// level-set field
      Array3<double> curvature,				// interface curvature
      Array3<double> unsmoothed_curvature,		// interface curvature without smoothing
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      int apply_curvature_smoothing,			// =1, apply curvature smoothing
							// =0, use unsmoothed curvature
      int number_curvature_smoothing_steps,		// number of iterations applied in the
							// curvature smoothing algorithm
      int apply_curvature_smoothing_filter 		// =1, apply curvature smoothing filter
	    						// =0, do not apply curvature smoothing filter		
	)
{

      void compute_curvature(				// compute the curvature
	 Array3<double> level_set, 				
	 Array3<double> curvature,				
	 int number_primary_cells_i,			
	 int number_primary_cells_j,			
	 int number_primary_cells_k,			
        double mesh_width_x1,				
        double mesh_width_x2,				
        double mesh_width_x3				
	);
      void smooth_curvature(                       // smooth the curvature
  	 Array3<double> level_set, 			
	 Array3<double> curvature_new,			
	 Array3<double> unsmoothed_curvature,			
	 int number_primary_cells_i,			
	 int number_primary_cells_j,			
	 int number_primary_cells_k,			
	 double mesh_width_x1,			
	 double mesh_width_x2,			
	 double mesh_width_x3,			
	 int apply_curvature_smoothing_filter,
	 int number_curvature_smoothing_steps
       );
      void copy_cell_centered_field( 		// copy cell centered field
	    Array3<double> source_field, 		   			
	    Array3<double> target_field,		
	    int number_primary_cells_i,		
	    int number_primary_cells_j,		
	    int number_primary_cells_k		
	    );
      
      
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
}
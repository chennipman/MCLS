#include "../headers/array.h"

#include<cstdlib>
#include<iostream>
/********************************************************************************/
/********************************************************************************/
/*  Function to advect the volume of fluid field                                */
/*  method. 									       */
/*  										       */
/*  Programmer	: Duncan van der Heul       					       */
/*  Date	: 10-03-2013       						       */
/*  Update	:        							       */
/********************************************************************************/
/* Notes									       */
/* The volume of fluid field is advected using an operator splitting technique  */
/* that is explained in both Sander's thesis as in the papers of Sussman        */
/* and is first order accurate. This does not guarantee that the volume of      */
/* fluid value remains within the interval [0,1]. 	 			       */
/* Furthermore, it can occur that a cell is too large to represent a small      */
/* feature. This leads to cells with intermediate values of the volume of 	*/
/* fluid in cells that are NOT interface cells. Of course this can not be 	*/
/* allowed. Both effects:                                                       */
/* -cells that are overfilled or underdrained					*/
/* -cells that contain numerical vapor                                          */
/* Have to handled in a way that does not disturb the global mass conservation. */
/* level set fields.								       */
/********************************************************************************/
      void apply_volume_of_fluid_redistribution(
	    Array3<double> volume_of_fluid, 				// volume of fluid field
	    Array3<double> level_set_star, 				// level set field at new time level
									// after convection and reinitialization
									// not mass conserving
	    Array3<double> level_set_new, 					// level set field at new time level
									// mass conserving
	    int number_primary_cells_i,				// number of primary (pressure) cells in x1 direction
	    int number_primary_cells_j,				// number of primary (pressure) cells in x2 direction
	    int number_primary_cells_k,				// number of primary (pressure) cells in x3 direction
	    double mesh_width_x1,					// grid spacing in x1 direction (uniform)
	    double mesh_width_x2,					// grid spacing in x2 direction (uniform)
	    double mesh_width_x3,					// grid spacing in x3 direction (uniform)
	    double volume_of_fluid_tolerance,			// tolerance on the value of the volume
									// of fluid value
	    double vof_2_level_set_tolerance,			// tolerance in the conversion from volume
									// of fluid value to level-set value
	    double lower_bound_derivatives,				// lower bound for the first partial derivatives
									// to consider it a limiting case of vanishing
									// partial derivatives
	    int number_vof_2_level_set_iterations,			// number of OUTER iterations in the conversion 
									// from volume of fluid to level-set
	    int number_iterations_ridder,				// maximum number of iterations allowed in the
									// nonlinear root finding algorithm
	    double time_step_mass_redistribution,			// time step for the mass redistribution
	    double redistribution_vof_tolerance,			// threshold value of time-derivative 
									// in volume of fluid redistribution equation
      	    int maximum_number_mass_redistribution_iterations 	// number of iterations allowed to make
									// the volume of fluid field valid
									// these are the sweeps on the vof error
	    )
      {
      void match_level_set_to_volume_of_fluid(			// compute a new level-set field
		 Array3<double> level_set_not_mass_conserving, 		// that corresponds to given
		 Array3<double> volume_of_fluid,				// volume of fluid field
		 Array3<double> level_set_mass_conserving,		// using an existing level-set
		 int number_primary_cells_i, 				// field as a starting value
		 int number_primary_cells_j, 
		 int number_primary_cells_k,			
		 double volume_of_fluid_tolerance,		
		 double lower_bound_derivatives,		
		 int number_vof_2_level_set_iterations,		
		 int number_iterations_ridder,			
		 double vof_2_level_set_tolerance		
		);
      void copy_cell_centered_field(					// copy a source field
		 Array3<double> source_field, 				// to a target field 
		 Array3<double> target_field, 				
		 int number_primary_cells_i, 
		 int number_primary_cells_j,
		 int number_primary_cells_k);
       void match_level_set_to_volume_of_fluid(			// compute a new level-set field
		 Array3<double> level_set_not_mass_conserving, 		// that corresponds to given
		 Array3<double> volume_of_fluid,				// volume of fluid field
		 Array3<double> level_set_mass_conserving,		// using an existing level-set
		 int number_primary_cells_i, 				// field as a starting value
		 int number_primary_cells_j, 
		 int number_primary_cells_k,			
		 double volume_of_fluid_tolerance,		
		 double lower_bound_derivatives,		
		 int number_vof_2_level_set_iterations,		
		 int number_iterations_ridder,			
		 double vof_2_level_set_tolerance,
		 double mesh_width_x1,
		 double mesh_width_x2,
		 double mesh_width_x3
		 );
     int modify_volume_of_fluid_values(				// modify volume of fluid
		Array3<double> level_set, 					// values and count number
		Array3<double> volume_of_fluid,			       // of cells with invalid	
		Array3<double> volume_of_fluid_correction,		// volume of fluid valuesluid_correction,
		int number_primary_cells_i,			  
		int number_primary_cells_j,			  
		int number_primary_cells_k,			  
		double volume_of_fluid_tolerance
		  );
      void redistribute_volume_of_fluid_error(			// solve artifical convection				
		Array3<double> level_set, 					// equation to redistribute
		Array3<double> volume_of_fluid, 				// volume of fluid error
		Array3<double> volume_of_fluid_correction,	
		int number_primary_cells_i,		 
		int number_primary_cells_j,		 
		int number_primary_cells_k,		 
		double mesh_width_x2,			
		double mesh_width_x1,			
		double mesh_width_x3,			
		double time_step_mass_redistribution,	
		double volume_of_fluid_tolerance,	
		double redistribution_vof_tolerance,	
		int maximum_number_mass_redistribution_iterations  	
		  );
      int apply_volume_of_fluid_clipping(				// clip the values of		
		 Array3<double> volume_of_fluid, 		              	// the volume of fluid field                     	
		 int number_primary_cells_i,		              	// to bring them in the correct
		 int number_primary_cells_j,		           	// interval [0,1]
		 int number_primary_cells_k,		
		 double volume_of_fluid_tolerance	
		
      );
      int check_volume_of_fluid(
               Array3<double> volume_of_fluid,
               int number_primary_cells_i,
               int number_primary_cells_j,
               int number_primary_cells_k,
               double volume_of_fluid_tolerance

              );
      Array3<double> volume_of_fluid_star;				// volume of fluid field, uncorrected
									// so with possible vapour cells and 
									// under/overfilled cells
      Array3<double> volume_of_fluid_correction;				// correction to the volume of fluid field
								       // to make it valid
      
      int number_cells_vof_out_of_bounds=100;	              // number of control volumes where the volume of fluid
						                     // function is OUTSIDE the interval [0,1]
      int number_cells_numerical_vapor=100;		              // number of control volumes where the volume of fluid
						                     // function is INSIDE the interval [0,1]
						                     // while the cell has 6 neighbours with the same sign:
						                     // the cell is NOT an interface cell
      int number_cells_invalid_volume_of_fluid
				      =200; 	                     // sum of number of vapour cells and number of cells
						                     // with the volume of fluid outside [0,1];
      int index_redistribution_attempt=0;	                     // number of attempts to achieve a
						                     // valid volume of fluid field
						                     // through the redistribution algorithm
						
      /* allocate memory for the volume of fluid correction and the tentative */
      /* volume of fluid field 						      */
      
	volume_of_fluid_correction.create(number_primary_cells_i+2, 
						    number_primary_cells_j+2,
						      number_primary_cells_k+2);
	volume_of_fluid_star.create(number_primary_cells_i+2, 
						    number_primary_cells_j+2,
						      number_primary_cells_k+2);
      
      /*---start the update loop ---*/
      /* it is unclear to me why this should be attempted multiple times */
      /* but it is possible the changes can be locally that strong it is */
      /* necessary.							 */
      
        std::cerr<<"in the redistribution algorithm \n";
      
	while(index_redistribution_attempt<=maximum_number_mass_redistribution_iterations &&
	      number_cells_numerical_vapor>0 &&          
		  number_cells_vof_out_of_bounds>0 )
	{
	  
		      std::cerr<<"iteration number "<< index_redistribution_attempt<<" \n";
		      std::cerr<<" out of "<<maximum_number_mass_redistribution_iterations <<"\n";

      /* copy the original volume of fluid field to the star volume of fluid field */

	 	copy_cell_centered_field(volume_of_fluid, volume_of_fluid_star, 				
		      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
      /* apply clipping to the tentative volume of fluid field */
      
	  	number_cells_vof_out_of_bounds=
	      		apply_volume_of_fluid_clipping(volume_of_fluid_star, 
		    		number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
		                  volume_of_fluid_tolerance);
			
		      std::cerr<<" number_cells_vof_out_of_bounds "<< number_cells_vof_out_of_bounds<<" \n";
                    exit(1);

      /* bring the level-set field in accordance with the clipped volume of fluid field */

	  	match_level_set_to_volume_of_fluid(level_set_star, volume_of_fluid_star, level_set_new,
								number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				  					volume_of_fluid_tolerance, lower_bound_derivatives,		
				    						number_vof_2_level_set_iterations, number_iterations_ridder,			
				      							vof_2_level_set_tolerance,
				      	 						mesh_width_x1, mesh_width_x2, mesh_width_x3);	
      
      /* modify the volume of fluid field: determine the error and the number of */
      /* invalid cells: both vapour cells and cells that have volume of fluid out of bounds */
      
	    	number_cells_invalid_volume_of_fluid=			
			modify_volume_of_fluid_values( level_set_new, volume_of_fluid, volume_of_fluid_correction,	
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			  
						      volume_of_fluid_tolerance);
		
		std::cerr<<" number_cells_invalid_volume_of_fluid "<< number_cells_invalid_volume_of_fluid<<" \n";

      /* if necessary redistribute the volume of fluid correction */
      /* so we end up with a valid field  */
      
      		if(number_cells_invalid_volume_of_fluid>0)
			
		{
      			redistribute_volume_of_fluid_error(level_set_new, volume_of_fluid, volume_of_fluid_correction,
						number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
							mesh_width_x1,	mesh_width_x1,	mesh_width_x3,			
								time_step_mass_redistribution, volume_of_fluid_tolerance,		
	 								redistribution_vof_tolerance,	
						 				maximum_number_mass_redistribution_iterations);
     
      
	  		index_redistribution_attempt++;
		}
	  

      /* copy the new tentative level-set field to the star level-set field */

	 	copy_cell_centered_field(level_set_new, level_set_star, 				
		      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);


	  
	  
	}
	
      /* if there are remaining invalid cells, apply the correction one more time */
      /* and bring the level-set field in accordance with the new volume of fluid */
      /* field 									  */

      	if(number_cells_invalid_volume_of_fluid>0)
	{
	    	number_cells_invalid_volume_of_fluid=			
			modify_volume_of_fluid_values(level_set_new, volume_of_fluid, volume_of_fluid_correction,			  
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			  
								volume_of_fluid_tolerance);

	    	match_level_set_to_volume_of_fluid(level_set_star, volume_of_fluid, level_set_new, 
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				  			   volume_of_fluid_tolerance, lower_bound_derivatives,		
				    				number_vof_2_level_set_iterations, number_iterations_ridder,			
				      					vof_2_level_set_tolerance,
				      	 			mesh_width_x1, mesh_width_x2, mesh_width_x3);
		
		  std::cerr<< "**************************************************** \n";
		  std::cerr<<"WARNING \n";
		  std::cerr<<"Despite the application of the volume of fluid  \n";
		  std::cerr<<"redistribution algorithm, there are still "<<
		  			number_cells_invalid_volume_of_fluid<< "\n";
		  std::cerr<<"cells that have a volume of fluid value      \n";
		  std::cerr<<"outside the physically relevant interval  \n";
		  std::cerr<<"in apply_volume_of_fluid_redistribution line 241. \n";
		  std::cerr<<"**************************************************** \n";
	    
	    
	}
	
	/* deallocate temporary storage */
	  
      volume_of_fluid_star.destroy();
      volume_of_fluid_correction.destroy();
}
  

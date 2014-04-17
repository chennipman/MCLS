#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<math.h>
/***********************************************************************************/
/***********************************************************************************/
/*  Function to redistribute the over and undershoots in the volume of             */
/*  fluid field. Additonally, 'vapour' cells, which contain an intermediate        */
/*  volume of fluid value but NO interface, are removed.                           */
/*                                                                                 */
/*  Programmer : Duncan van der Heul                                               */
/*  Date       : 15-08-2013                                                        */
/*  Update     :                                                                   */
/***********************************************************************************/
/* Notes                                                                           */
/* Due to the operator splitting in the volume of fluid advection equation         */
/* over and undershoots can occur in this field.                                   */
/* This function will advect the over and undershoots towards the interface        */
/* and update the interface position accordingly.                                  */
/***********************************************************************************/

EXPORT void redistribute_volume_of_fluid_error(						
	Array3<double> level_set_old, 				// level set field 
								// before mass redistribution
        Array3<double> level_set,                               // level set field updated after
                                                                // mass redistribution
	Array3<double> volume_of_fluid, 			// volume of fluid field
	Array3<double> volume_of_fluid_correction,		// correction to the volume of fluid field
								// to make it valid
        Array3<double> invalid_vof_cells,                       // indication field to show what is wrong
                                                                // indicator field showing cells that are either
                                                                // within bounds =0
                                                                // underfilled   =-1
                                                                // overfilled    =+1
                                                                // vapour cells  = 5
	int number_primary_cells_i,				// number of primary (pressure) cells in x1 direction
	int number_primary_cells_j,				// number of primary (pressure) cells in x2 direction
	int number_primary_cells_k,				// number of primary (pressure) cells in x3 direction
	double mesh_width_x1,					// grid spacing in x1 direction (uniform)
	double mesh_width_x2,					// grid spacing in x2 direction (uniform)
	double mesh_width_x3,					// grid spacing in x3 direction (uniform)
	double time_step_mass_redistribution,		        // time step for the mass redistribution
								// algorithm
        double volume_of_fluid_tolerance,			// tolerance for volume of fluid value
	double redistribution_vof_tolerance,			// threshold value of time-derivative 
								// in volume of fluid redistribution equation
        int maximum_number_mass_redistribution_iterations	// number of iterations allowed to make
								// the volume of fluid field valid
								// these are the sweeps on the vof error
 
 	)
     {
	Array3<double> redistribution_velocity_x1;			// artificial redistribution velocity x1 direction
	Array3<double> redistribution_velocity_x2;			// artificial redistribution velocity x2 direction
	Array3<double> redistribution_velocity_x3;			// artificial redistribution velocity x3 direction
	Array3<double> time_derivative_volume_of_fluid_correction;	// time derivative in the discretised
									// volume of fluid redistribution equation
	double time_step_vof_error_distribution;			// time step used in the volume of fluid
									// error redistribution algorithm
	
       
	double cfl_number_vof_error_distribution=0.25;       		// cfl number used in the volume of fluid
									// error redistribution algorithm
	double maximum_time_derivative_vof_error=100;			// largest value of time derivative in the volume of fluid
									// error redistribution equation
        int number_cells_vof_out_of_bounds=0;			        // number of control volumes where the volume of fluid
									// function is OUTSIDE the interval [0,1]
        int number_cells_numerical_vapor=0;                             // number of control volumes where the volume of fluid
                                                                        // function is INSIDE the interval [0,1]
                                                                        // while the cell has 6 neighbours with the same sign:
                                                                        // the cell is NOT an interface cell
        int number_cells_invalid_volume_of_fluid=100;                   // sum of number of vapour cells and number of cells
                                                                        // with the volume of fluid outside [0,1];
	int i_index, j_index, k_index;  				// local variables for loop indexing
	int iteration_index=0;						// index of the iteration/sweep in the 
									// vof redistribution algorithm
        bool detailed_output=0;                                         // =1, provide detailed output
                                                                        // =0, provide non output 
							
	/* compute the time step for the error redistribution algorithm */
	
	/* three-dimensional case */
	
	time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution*(mesh_width_x1+mesh_width_x2+mesh_width_x3)/3.0;
	
	/* two-dimensional cases */
	
	if(number_primary_cells_i==1){
	    time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution*(mesh_width_x2+mesh_width_x3)/2.0;
	}
	if(number_primary_cells_j==1){
	    time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution*(mesh_width_x1+mesh_width_x3)/2.0;
	}
	if(number_primary_cells_k==1){
	    time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution*(mesh_width_x1+mesh_width_x2)/2.0;
	}
       
        /* allocate memory for the 'artificial' advection velocities in the update equation */
        /* and the time-derivative of the volume of fluid error                             */
	
      
	 redistribution_velocity_x1.create(number_primary_cells_i+1, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	 redistribution_velocity_x2.create(number_primary_cells_i+2, number_primary_cells_j+1,
				    number_primary_cells_k+2);
	 redistribution_velocity_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+1);
  	 time_derivative_volume_of_fluid_correction.create(number_primary_cells_i+2, number_primary_cells_j+2,
					  number_primary_cells_k+2);

        /* compute_redistribution_velocity_field */
        
        compute_redistribution_velocity_field( level_set,                               
                                                 redistribution_velocity_x1, redistribution_velocity_x2, redistribution_velocity_x3,              
                                                   mesh_width_x1, mesh_width_x2, mesh_width_x3,                                   
                                                     number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);                              
        
        	
	/* Start the iterative process of finding appropriate corrections to the        */
	/* volume of fluid field, such that the over and undershoots are eliminated.    */
	/* Note that the total correction is updated iteratively until it fulfills      */
	/* the requirements, while the original volume of fluid field is left unchanged */
	
	
	while(iteration_index<10 && number_cells_invalid_volume_of_fluid>0 &&
	      maximum_time_derivative_vof_error>redistribution_vof_tolerance)
	{
            /* reset all numbers of invalid cells  */
            /* and the maximum change in the error */
	  
	    number_cells_vof_out_of_bounds=0;
            number_cells_invalid_volume_of_fluid=0;
            number_cells_numerical_vapor=0;
	    maximum_time_derivative_vof_error=0.0;
            
            iteration_index++;
            
            /* compute the time-derivative of the volume of fluid redistribution equation */
            
            maximum_time_derivative_vof_error=
                                compute_redistribution_time_derivative(
                                                redistribution_velocity_x1, redistribution_velocity_x2, redistribution_velocity_x3,                      
                                                   time_derivative_volume_of_fluid_correction, volume_of_fluid_correction,
                                                      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                     
                                                        mesh_width_x1, mesh_width_x2, mesh_width_x3);
                                
            /* update the correction for the volume of fluid field */
                                
            for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
                for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
                    for(k_index=1;k_index<number_primary_cells_k+1;k_index++){
                               
                        volume_of_fluid_correction[i_index][j_index][k_index]+=
                                time_step_vof_error_distribution*
                                    time_derivative_volume_of_fluid_correction[i_index][j_index][k_index];
                   }
                }
            }
 
            /* verify if the correction to the volume of fluid field will bring it into the valid range */
            /* if not, then perform additional sweeps  */
	
            analyse_validity_vof_correction(level_set, volume_of_fluid_correction, volume_of_fluid, invalid_vof_cells,                               
                                              number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                 
                                                volume_of_fluid_tolerance, 
                                                   number_cells_vof_out_of_bounds, number_cells_numerical_vapor,                           
                                                      number_cells_invalid_volume_of_fluid);
	    
 	}
	
	/* use the correction to update the volume of fluid values */

	for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
	    for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
		for(k_index=0;k_index<number_primary_cells_k+1;k_index++){
			  
		    volume_of_fluid[i_index][j_index][k_index]+=
			    volume_of_fluid_correction[i_index][j_index][k_index];
		}
	    }
	}
	
	/* analyse if the volume of fluid field is now valid */
	
         analyse_validity_vof_field( level_set, volume_of_fluid, invalid_vof_cells,                               
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                    
                                          volume_of_fluid_tolerance, number_cells_vof_out_of_bounds, number_cells_numerical_vapor,                             
                                                number_cells_invalid_volume_of_fluid);
         
        /* if detailed output is required, write all fields to file */

         if(detailed_output)
         {
                dump_redistribution_for_debugging(level_set_old, volume_of_fluid,        
                                                    level_set, invalid_vof_cells, volume_of_fluid_correction,          
                                                      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                    
                                                        mesh_width_x1, mesh_width_x2, mesh_width_x3 );
         }
         
        /* deallocate the memory for the 'artificial' advection velocities in the update equation */

        redistribution_velocity_x1.destroy();
	redistribution_velocity_x2.destroy();
	redistribution_velocity_x3.destroy();
        time_derivative_volume_of_fluid_correction.destroy();
     }

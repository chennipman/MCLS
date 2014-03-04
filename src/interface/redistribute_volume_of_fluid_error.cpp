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
	Array3<double> level_set, 				// level set field 
								// mass conserving
	Array3<double> volume_of_fluid, 			// volume of fluid field
	Array3<double> volume_of_fluid_correction,		// correction to the volume of fluid field
								// to make it valid
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
	
       
	double cfl_number_vof_error_distribution=0.5;       		// cfl number used in the volume of fluid
									// error redistribution algorithm
	double flux_minus_x1;						// flux of error advection equation at plus i face
	double flux_minus_x2;						// flux of error advection equation at plus j face
	double flux_minus_x3;						// flux of error advection equation at plus k face
	double flux_pluss_x1;						// flux of error advection equation at minus i face
	double flux_pluss_x2;						// flux of error advection equation at minus j face
	double flux_pluss_x3;						// flux of error advection equation at minus k face
	
	double level_set_value;					// level set field in cell center
	double level_set_minusi;					// level set field at minus i neighbouring cell
	double level_set_minusj;					// level set field at minus j neighbouring cell
	double level_set_minusk;					// level set field at minus k neighbouring cell
	double level_set_plusi;					// level set field at plus i neighbouring cell
	double level_set_plusj;					// level set field at plus j neighbouring cell
	double level_set_plusk;					// level set field at plus k neighbouring cell
	double current_volume_of_fluid;				// uncorrected, cell centered value of volume of fluid
	double maximum_time_derivative_vof_error;			// largest value of time derivative in the volume of fluid
									// error redistribution equation
        int number_cells_vof_out_of_bounds=100;			// number of control volumes where the volume of fluid
									// function is OUTSIDE the interval [0,1]
	int i_index, j_index, k_index;  				// local variables for loop indexing
	int iteration_index;						// index of the iteration/sweep in the 
									// vof redistribution algorithm
	bool pure_cell;						// true: all vertices on one side of interface,
									// false: not all vertices on one side of interface
							
	/* compute the time step for the error redistribution algorithm */
	
	/* three-dimensional case */
	
	time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution/(mesh_width_x1+mesh_width_x2+mesh_width_x3);
	
	/* two-dimensional cases */
	
	if(number_primary_cells_i==1){
	    time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution/(mesh_width_x2+mesh_width_x3);
	}
	if(number_primary_cells_j==1){
	    time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution/(mesh_width_x1+mesh_width_x3);
	}
	if(number_primary_cells_k==1){
	    time_step_vof_error_distribution=
		    cfl_number_vof_error_distribution/(mesh_width_x2+mesh_width_x3);
	}
       
        /* allocate memory for the 'artificial' advection velocities in the update equation */
	
      
	redistribution_velocity_x1.create(number_primary_cells_i+1, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	redistribution_velocity_x2.create(number_primary_cells_i+2, number_primary_cells_j+1,
				    number_primary_cells_k+2);
	redistribution_velocity_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+1);
	
	/* allocate memory for the time-derivative of the volume of fluid error */
	
	time_derivative_volume_of_fluid_correction.create(number_primary_cells_i+2, number_primary_cells_j+2,
					  number_primary_cells_k+2);

	/* compute advection velocity in 1-direction */
	
	for( i_index=0;i_index<number_primary_cells_i+1;i_index++){
	    for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++){
		  
		  redistribution_velocity_x1[i_index][j_index][k_index]=
					      compute_redistribution_velocity(level_set[i_index][j_index][k_index],
						level_set[i_index+1][j_index][k_index],
						  mesh_width_x1);
		}
	    }
	}
	
	/* compute advection velocity in 2-direction */
	
	for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
	    for(j_index=0;j_index<number_primary_cells_j+1;j_index++){
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++){
		  
		  redistribution_velocity_x2[i_index][j_index][k_index]=
					      compute_redistribution_velocity(level_set[i_index][j_index][k_index],
						level_set[i_index][j_index+1][k_index],
						  mesh_width_x2);
		}
	    }
	}
	
	
	/* compute advection velocity in 3-direction */
	
	
	for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
	    for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
		for(k_index=0;k_index<number_primary_cells_k+1;k_index++){
		  
		  redistribution_velocity_x3[i_index][j_index][k_index]=
					      compute_redistribution_velocity(level_set[i_index][j_index][k_index],
						level_set[i_index][j_index][k_index+1],
						  mesh_width_x3);
		}
	    }
	}
	
	
	/* Start the iterative process of finding appropriate corrections to the        */
	/* volume of fluid field, such that the over and undershoots are eliminated.    */
	/* Note that the total correction is updated iteratively until it fulfills      */
	/* the requirements, while the original volume of fluid field is left unchanged */
	
	iteration_index=0;
	number_cells_vof_out_of_bounds=100;
	maximum_time_derivative_vof_error=100;
	
	while(iteration_index<=maximum_number_mass_redistribution_iterations&& number_cells_vof_out_of_bounds>0 &&
	      maximum_time_derivative_vof_error>redistribution_vof_tolerance)
	{
	  
	    number_cells_vof_out_of_bounds=0;
	    maximum_time_derivative_vof_error=0.0;
	    
	    for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
		for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
		    for(k_index=1;k_index<number_primary_cells_k+1;k_index++){
		  
			flux_pluss_x1=
				upwind_flux_mass_redistribution(
				    redistribution_velocity_x1[i_index  ][j_index][k_index],
					volume_of_fluid_correction[i_index  ][j_index][k_index],
					    volume_of_fluid_correction[i_index+1][j_index][k_index]);
			flux_minus_x1=
				upwind_flux_mass_redistribution(
				    redistribution_velocity_x1[i_index-1][j_index][k_index],
					volume_of_fluid_correction[i_index-1][j_index][k_index],
					    volume_of_fluid_correction[i_index  ][j_index][k_index]);
			flux_pluss_x2=
				upwind_flux_mass_redistribution(
				    redistribution_velocity_x2[i_index][j_index  ][k_index],
					volume_of_fluid_correction[i_index][j_index  ][k_index],
					    volume_of_fluid_correction[i_index][j_index+1][k_index]);
			flux_minus_x2=
				upwind_flux_mass_redistribution(
				    redistribution_velocity_x2[i_index][j_index-1][k_index],
					volume_of_fluid_correction[i_index][j_index  ][k_index],
					    volume_of_fluid_correction[i_index][j_index+1][k_index]);
			flux_pluss_x3=
				upwind_flux_mass_redistribution(
				    redistribution_velocity_x3[i_index][j_index][k_index  ],
					volume_of_fluid_correction[i_index][j_index][k_index  ],
					    volume_of_fluid_correction[i_index][j_index][k_index+1]);
			flux_minus_x3=
				upwind_flux_mass_redistribution(
				    redistribution_velocity_x3[i_index][j_index][k_index-1],
					volume_of_fluid_correction[i_index][j_index][k_index  ],
					    volume_of_fluid_correction[i_index][j_index][k_index+1]);
				
			time_derivative_volume_of_fluid_correction[i_index][j_index][k_index]=
					      -1.0*(
					mesh_width_x1*(flux_pluss_x1-flux_minus_x1)+
					mesh_width_x2*(flux_pluss_x2-flux_minus_x2)+
					mesh_width_x3*(flux_pluss_x3-flux_minus_x3)
						   );
			maximum_time_derivative_vof_error=std::max(maximum_time_derivative_vof_error,
			    fabs(time_derivative_volume_of_fluid_correction[i_index][j_index][k_index]));
	    
		    }
		}
	    }
	
	
	/* verify if the correction to the volume of fluid field will bring it into the valid range */
	/* update the volume of fluid field */
	
	
	    
	    for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
		for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
		    for(k_index=0;k_index<number_primary_cells_k+1;k_index++){

                     /* update the volume of fluid correction for this cell */
                     
			volume_of_fluid_correction[i_index][j_index][k_index]+=
				time_step_vof_error_distribution*
				    time_derivative_volume_of_fluid_correction[i_index][j_index][k_index];

                     /* these are the current values for level set and  */
                     /* volume of fluid for this cell, note that the    */
                     /* current volume of fluid cell is not stored in   */
                     /* the array volume of fluid until the algorithm   */
                     /* actually is converged.                          */
                     
			level_set_value =level_set[i_index][j_index][k_index];
                      current_volume_of_fluid=volume_of_fluid[i_index][j_index][k_index]+
                                   volume_of_fluid_correction[i_index][j_index][k_index];

                     /* determine the level set value in all the neighbouring */
                     /* cells to identify vapour cells: cells that have an    */
                     /* intermediate volume of fluid value, but contain NO interface */
                     
			level_set_plusi =level_set[i_index+1][j_index  ][k_index  ];
			level_set_minusi=level_set[i_index-1][j_index  ][k_index  ];
			level_set_plusj =level_set[i_index  ][j_index+1][k_index  ];
			level_set_minusj=level_set[i_index  ][j_index+1][k_index  ];
			level_set_plusk =level_set[i_index  ][j_index  ][k_index+1];
			level_set_minusk=level_set[i_index  ][j_index  ][k_index+1];
			
			
			
			/* two conditions are checked to verify if the volume of fluid value is valid :    */
			/*1) if this is a NON-interface cell: the volume of fluid should also be either 0 or 1 */
			/*				with some tolerance 				   */
			/*2) if this is a mixed cell: the volume of fluid should be in the interval 0 to 1 */
			/*				with some tolerance 				   */
			pure_cell=((level_set_value*level_set_plusi>0)&&
				    (level_set_value*level_set_plusj>0)&&
				      (level_set_value*level_set_plusk>0)&&
					(level_set_value*level_set_minusi>0)&&
					  (level_set_value*level_set_minusj>0)&&
					    (level_set_value*level_set_minusk>0));
			if( (pure_cell && (current_volume_of_fluid<1.0-volume_of_fluid_tolerance ||
					    current_volume_of_fluid>volume_of_fluid_tolerance)) ||
			    (current_volume_of_fluid > 1.0+volume_of_fluid_tolerance||
					    current_volume_of_fluid<volume_of_fluid_tolerance)){
			
			      number_cells_vof_out_of_bounds++;
			}
		    }
		}
	    }
	
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
	
	/* deallocate the memory for the 'artificial' advection velocities in the update equation */

	redistribution_velocity_x1.destroy();
	redistribution_velocity_x2.destroy();
	redistribution_velocity_x3.destroy();
     }

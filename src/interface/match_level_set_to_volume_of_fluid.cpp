#include "../headers/array.h"

#include<cstdlib>
#include<iostream>
#include<stdio.h>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to make the level-set field mass conserving                        */
/*  method. 									*/
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
EXPORT void match_level_set_to_volume_of_fluid(				
		 Array3<double> level_set_star,			// non mass conserving level-set field
		 Array3<double> volume_of_fluid,		// volume of fluid field
		 Array3<double> level_set_mass_conserving,	// corrected, mass conserving level-set field
		 int number_primary_cells_i,			// number of primary (pressure) 
								// cells in x1 direction
		 int number_primary_cells_j,			// number of primary (pressure) 
								// cells in x2 direction
		 int number_primary_cells_k,			// number of primary (pressure) 
								// cells in x3 direction
		 double volume_of_fluid_tolerance,		// tolerance in volume of fluid field
								// for admissable values
		 double lower_bound_derivatives,		// lower bound for the first partial derivatives
								// to consider it a limiting case of vanishing
								// partial derivatives
		 int number_vof_2_level_set_iterations,	        // number of OUTER iterations in the conversion 
								// from volume of fluid to level-set
		 int number_iterations_ridder,		        // maximum number of iterations allowed in the
							        // nonlinear root finding algorithm
								// this is the number of INNER iterations
		 double vof_2_level_set_tolerance,		// tolerance in the conversion from volume
								// of fluid value to level-set value
	  	 double mesh_width_x1,			        // grid spacing in x1 direction (uniform)
	  	 double mesh_width_x2,			        // grid spacing in x2 direction (uniform)
	  	 double mesh_width_x3				// grid spacing in x3 direction (uniform)
		 )
		{
	Array3<double> d_level_set_d_x1; 			// first partial derivative of level-set with 
							        // respect to x1, central approximation
	Array3<double> d_level_set_d_x2; 			// first partial derivative of level-set with 
							        // respect to x2, central approximation 
	Array3<double> d_level_set_d_x3; 			// first partial derivative of level-set with 
							        // respect to x3, central approximation
	Array3<double> volume_of_fluid_star;		        // volume of fluid field, evaluated from the star
							        // level-set field
	Array3<double> level_set_correction;		        // correction to be applied to the level-set field
							        // to make it mass-conserving
	double maximum_level_set_correction=0.0;		// maximum of correction that needs to be applied
							        // to the convected level-set field to make it mass
							        // conserving
	double maximum_volume_of_fluid_deviation=0.0;	        // maximum of volume of fluid deviation
	Array3<double> volume_of_fluid_deviation;		// difference between the converted, advected level
							        // set field and the advected volume of fluid field
	
	int i_index, j_index, k_index;  	                // local variables for loop indexing
	int correction_step_index;		                // index of the iteration number to make the level-set
						                // field mass-conserving
	int number_of_error_cells=0;    	                // Number of cells where the volume of fluid fiels is outside the
						                // allowed interval [0,1] (with proper tolerances)
	int number_of_cells_to_correct=0;	                // number of cells for which the level-set does
						                // not match the mass conservation criterion
	
      
      /* allocate memory for the volume of fluid field corresponding to    */
      /* the star level-set field, the gradient of this level-set field    */
      /* and the correction to the level-set field	and volume of fluid deviation		   */
      
	volume_of_fluid_star.create( number_primary_cells_i+2, number_primary_cells_j+2,
			  number_primary_cells_k+2);
	d_level_set_d_x1.create(number_primary_cells_i+2, number_primary_cells_j+2,
			  number_primary_cells_k+2);
	d_level_set_d_x2.create(number_primary_cells_i+2, number_primary_cells_j+2,
			  number_primary_cells_k+2);
	d_level_set_d_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
			  number_primary_cells_k+2);
	level_set_correction.create(number_primary_cells_i+2, number_primary_cells_j+2,
			  number_primary_cells_k+2);
	volume_of_fluid_deviation.create(number_primary_cells_i+2, number_primary_cells_j+2,
			  number_primary_cells_k+2);
      
	
      /* check the volume of fluid field for containing only 'valid' values */
      
	number_of_error_cells=check_volume_of_fluid(volume_of_fluid, 
				    number_primary_cells_i, number_primary_cells_j,
					number_primary_cells_k,	
					  volume_of_fluid_tolerance);
	
	if(number_of_error_cells)
	{
	  /* there are cells with their volume of fluid outside the interval [0,1] */
	  
	  std::cerr << "***************************************************** \n";
	  std::cerr<<"ERROR \n";
	  std::cerr<<"A number of cells have a VOF value outside the interval [0,1] \n";
	  std::cerr<<"the number of invalid VOF cells is "<< number_of_error_cells << "\n";
	  std::cerr<<"in function match_level_set_to_volume_of_fluid line 152 \n";
	  std::cerr << "***************************************************** \n";
	  exit(1);
	}
	  
       /* start with copying the star level-set field to the mass conservative */
       /* level-set field						       */
	
	  copy_cell_centered_field( level_set_star, level_set_mass_conserving,
				      number_primary_cells_i, number_primary_cells_j,
					number_primary_cells_k);
	  
	/* match the level-set field to the volume of fluid through a set of sweeps */
	/* note that this is an iterative process, because the inversion needs 	    */
	/* information on the gradient of the level-set that will change in this    */
	/* process */
	  
	 for(correction_step_index=1;
		      correction_step_index<=number_vof_2_level_set_iterations;correction_step_index++)
	 {
	  
       /* compute the gradient of the level-set field */
      
	      compute_level_set_gradient(level_set_mass_conserving, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				     number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);

      /* compute the volume of fluid field corresponding to the 'actual' level-set field */ 	
      
	      compute_volume_of_fluid(level_set_mass_conserving, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3, 
				  volume_of_fluid_star, number_primary_cells_i, number_primary_cells_j,
			   number_primary_cells_k, lower_bound_derivatives);

	      maximum_volume_of_fluid_deviation=0.0;
	      number_of_cells_to_correct=0;
	  
	      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	      {
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			  {
			      volume_of_fluid_deviation[i_index][j_index][k_index]=
			      				volume_of_fluid[i_index][j_index][k_index]-
				                     volume_of_fluid_star[i_index][j_index][k_index];
			      maximum_volume_of_fluid_deviation=
			      std::max(fabs(volume_of_fluid_deviation[i_index][j_index][k_index]), maximum_volume_of_fluid_deviation);
			      if(fabs(volume_of_fluid_deviation[i_index][j_index][k_index])>volume_of_fluid_tolerance)
			      {
				number_of_cells_to_correct++;
			      }
			  }
		    }
	      }
	      
// 	      std::cerr<<"number of cells to correct "<< number_of_cells_to_correct<<"\n";
// 	      std::cerr<<"maximum_volume_of_fluid_deviation "<< maximum_volume_of_fluid_deviation<<"\n";
	      
	      if( maximum_volume_of_fluid_deviation < volume_of_fluid_tolerance  )
	      {
	        /* the volume of fluid field corresponding to the level-set field is within */
		/* +/- tolerance   of the volume_of_fluid field at the new time level       */
		/* no further correction is required 					    */
		  break;
	      }
	      else
	      {
		/* the level-set field needs to be corrected to make it mass-conserving */
	  
		  maximum_level_set_correction=0.0;
		  for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		  {
			for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
			{
			    for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			    {
				volume_of_fluid_deviation[i_index][j_index][k_index]=volume_of_fluid[i_index][j_index][k_index]-
				                     volume_of_fluid_star[i_index][j_index][k_index];
				if(fabs(volume_of_fluid_deviation[i_index][j_index][k_index])>volume_of_fluid_tolerance)
				{
					
				   if( volume_of_fluid[i_index][j_index][k_index] <1E-12){
// 					  std::cerr<<"dangerously low value of vof field "<<
// 					  		volume_of_fluid[i_index][j_index][k_index]<<"\n";
// 					  std::cerr<<"index "<< i_index<<" "<<j_index<<" "<<k_index<<"\n";
					   volume_of_fluid[i_index][j_index][k_index]=1E-12;
				   }
				   if( volume_of_fluid[i_index][j_index][k_index] >1-1E-12){
// 			  		std::cerr<<"dangerously high value of "<<
// 					  		volume_of_fluid[i_index][j_index][k_index]<<"\n";
// 					  std::cerr<<"index "<< i_index<<" "<<j_index<<" "<<k_index<<"\n";
// 					   volume_of_fluid[i_index][j_index][k_index]-=1-1E-9;
					   volume_of_fluid[i_index][j_index][k_index]=1-1E-12;
				   }
				  
		/* compute the level-set field corresponding to the volume of fluid value */

				    if(vof_2_level_set(level_set_mass_conserving[i_index][j_index][k_index], 
						    d_level_set_d_x1[i_index][j_index][k_index],
						      d_level_set_d_x2[i_index][j_index][k_index], 
							d_level_set_d_x3[i_index][j_index][k_index], 
							  volume_of_fluid[i_index][j_index][k_index],
							  lower_bound_derivatives,    		
							    vof_2_level_set_tolerance,		
							      number_iterations_ridder))
				    {
					    std::cerr<<"problem occured for cell index ";
				    	    std::cerr<< i_index <<" "<<j_index<<" "<<k_index<<" \n";
					    exit(1);    
				    }
		/* the level-set correction is only required for debugging/visualisation purposes */	
 				    level_set_correction[i_index][j_index][k_index]=
					level_set_mass_conserving[i_index][j_index][k_index]-
					level_set_star[i_index][j_index][k_index];
                                    maximum_level_set_correction=std::max(maximum_level_set_correction,
                                                            fabs(level_set_correction[i_index][j_index][k_index]));
				}
			    }
			}
		  }
          
         /* update the virtual cells with the new level-set field */ 
	 
		  field_neumann_boundary(level_set_mass_conserving, number_primary_cells_i, 
							number_primary_cells_j,
							  number_primary_cells_k);
// 		  field_extrapolate_boundary(level_set_mass_conserving, number_primary_cells_i, 
// 							number_primary_cells_j,
// 							  number_primary_cells_k);
	      }
      }
      
        /* if after the allowed number of iterations the volume of fluid field */
	/* is still not matching the volume of fluid field induced by the level set field */
	/* something is wrong and all quantities involved are written to file  */

      if(maximum_volume_of_fluid_deviation>volume_of_fluid_tolerance)
      {	
	    std::cerr<<" Level-set field could not be corrected properly, aborting...\n";
	    std::cerr<<" in match_level_set_to_volume_of_fluid Line 278. \n";
	
	    if(correction_step_index>=1)
	    {
	      /* the deviation in the volume of fluid field is still too large  */
	      /* and the number of corrections applied is larger than 1         */
	      /* proceed with writing all quantities involved to file           */
	      /* this is purely for debugging purposes 				*/
	      
  		dump_solution_for_debugging( level_set_star, volume_of_fluid,			
						   level_set_mass_conserving,		
						      level_set_correction, volume_of_fluid_deviation,	
							  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			 
							mesh_width_x1,	mesh_width_x2, mesh_width_x3);
	      
	      
	    }
	    exit(1);
      }

	/* deallocate the temporary storage used in this function */
	
	volume_of_fluid_star.destroy();
	d_level_set_d_x1.destroy();
	d_level_set_d_x2.destroy();
	d_level_set_d_x3.destroy();
	level_set_correction.destroy();
       volume_of_fluid_deviation.destroy();
   }

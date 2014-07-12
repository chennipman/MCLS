#include "../headers/array.h"

#include<cstdlib>
#include<iostream> 
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
/* This function applies the correction to the advected level-set field to      */
/* make it mass conserving. To do so the volume of fluid field is advected      */
/* and next this advected volume of fluid field is used to make a correction    */
/* to the level-set field.                                                      */
/* This function should be split up further.                                    */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/
EXPORT void make_level_set_mass_conserving
    (
      	Array3<double> level_set_star, 				// level set field at new time level
								// after convection and reinitialization
								// not mass conserving
      	Array3<double> level_set_new, 				// level set field at new time level
								// mass conserving
      	Array3<double> level_set_old, 				// level set field at old time level
								// mass conserving
      	Array3<double> volume_of_fluid, 			// volume of fluid field
      	Array3<double> u_1_velocity_new, 			// velocity field at new time level x1 direction
      	Array3<double> u_2_velocity_new, 			// velocity field at new time level x2 direction
      	Array3<double> u_3_velocity_new,			// velocity field at new time level x3 direction
      	int number_primary_cells_i,				// number of primary (pressure) cells in x1 direction
      	int number_primary_cells_j,				// number of primary (pressure) cells in x2 direction
      	int number_primary_cells_k,				// number of primary (pressure) cells in x3 direction
      	double actual_time_step_level_set,			// time step used for level-set advection
								// computed from all stability restrictions and 
								// possibly subscycling
      	double mesh_width_x1,					// grid spacing in x1 direction (uniform)
      	double mesh_width_x2,					// grid spacing in x2 direction (uniform)
      	double mesh_width_x3,					// grid spacing in x3 direction (uniform)
      	int apply_mass_distribution_algorithm,    		// =1:apply the mass redistribution algorithm
								// to avoid numerical vapour
      	int apply_mass_conservation_correction,    		// =1:apply the mass conservation correction
								// to the level-set field
      	double volume_of_fluid_tolerance,			// tolerance in volume of fluid field
								// for admissable values
      	double lower_bound_derivatives,			        // lower bound for the first partial derivatives
								// to consider it a limiting case of vanishing
								// partial derivatives
      	int number_vof_2_level_set_iterations,		        // number of OUTER iterations in the conversion 
								// from volume of fluid to level-set
      	int number_iterations_ridder,				// maximum number of iterations allowed in the
								// nonlinear root finding algorithm
								// this is the number of INNER iterations
      	double vof_2_level_set_tolerance,			// tolerance in the conversion from volume
								// of fluid value to level-set value
      	int maximum_number_mass_redistribution_iterations,      // number of iterations allowed to make
								// the volume of fluid field valid
								// these are the sweeps on the vof error
      	double time_step_mass_redistribution,		        // time step for the mass redistribution
								// algorithm
	double redistribution_vof_tolerance, 			// threshold value of time-derivative 
								// in volume of fluid redistribution equation
        double mass_redistribution_diffusion_coefficient        // diffusion coefficient for mass redistribution equation
	    
    )
    {
      Array3<double> vof_after_x1_update;			// intermediate volume of fluid field, 
							        // after update in x1 direction
      Array3<double> vof_after_x2_update;			// intermediate volume of fluid field, 
							        // after update in x2 direction
      Array3<double> vof_after_x3_update;			// intermediate volume of fluid field, 
							        // after update in x3 direction
      Array3<double> d_level_set_d_x1;			        // first partial derivative of
							        // the level-set field wrt x1
							        // second order central approximation
      Array3<double> d_level_set_d_x2;			        // first partial derivative of 
							        // the level-set field wrt x2
							        // second order central approximation
      Array3<double> d_level_set_d_x3;			        // first partial derivative of
 							        // the level-set field wrt x3
							        // second order central approximation
      Array3<double> level_set_after_x1_update;		        // level-set field, converted from 
							        // intermediate	volume of fluid field, 
							        // after update in x1 direction
      Array3<double> level_set_after_x2_update;		        // level-set field, converted from 
							        // intermediate	volume of fluid field, 
							        // after update in x2 direction
      Array3<double> level_set_after_x3_update;		        // level-set field, converted from 
 							        // intermediate	volume of fluid field, 
							        // after update in x3 direction
      Array3<double> invalid_vof_cells;                         // indication field to show what is wrong
                                                                // indicator field showing cells that are either
                                                                // within bounds =0
                                                                // underfilled   =-1
                                                                // overfilled    =+1
                                                                // vapour cells  = 5
     
      static int order_of_updates;			        // indicates the order of the updates
							        // of volume of fluid field due to
							        // fluxes in different directions
							        // =1
							        // order = x1, x2, x3
							        // =2
							        // order = x2, x3, x1
							        // =3
							        // order = x3, x1, x2
							        // =4 (not implemented)
							        // order = x2, x1, x3
							        // =5 (not implemented)
							        // order = x1, x3, x2
							        // =6 (not implemented)
							        // order = x3, x2, x1
      
      int i_index, j_index, k_index;  		                // local variables for loop indexing
      int number_cells_vof_out_of_bounds;                       // number of control volumes where the volume of fluid
                                                                // function is OUTSIDE the interval [0,1]
      int number_cells_numerical_vapor;                         // number of control volumes where the volume of fluid
                                                                // function is INSIDE the interval [0,1]
                                                                // while the cell has 6 neighbours with the same sign:
                                                                // the cell is NOT an interface cell
      int number_cells_invalid_volume_of_fluid;                 // sum of number of vapour cells and number of cells
                                                                // with the volume of fluid outside [0,1];
      int number_clipped_cells;                                 // number of cells that have been clipped
      double one_over_dx1	=    			        // 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      double one_over_dx2	=    			        // 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
      double one_over_dx3	=    			        // 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);
      
/* allocate memory for the temporary storage of the intermediate fields */

	vof_after_x1_update.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	vof_after_x2_update.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	vof_after_x3_update.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	
	level_set_after_x1_update.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	level_set_after_x2_update.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	level_set_after_x3_update.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	
	d_level_set_d_x1.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	d_level_set_d_x2.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
	d_level_set_d_x3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
        
	invalid_vof_cells.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
				    
    /* keep order_of_updates in its allowable range */
    
      if( order_of_updates > 6)
      { 	
	order_of_updates=1;
      }
      
    /* Sander chose to keep the order_of_updates fixed to 1 */    
    
	order_of_updates=1;
	
	/* compute the gradient of the old level-set field, necessary for the level-set/vof conversion */
        
	compute_level_set_gradient(level_set_old, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				      );
	/* compute the current volume of fluid field */
	
	if(!compute_volume_of_fluid(level_set_old, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				volume_of_fluid,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,					  
				  lower_bound_derivatives))
	{
	  
	};

	if(apply_mass_conservation_correction)
	{
		
            std::cerr<<"start mass conserving correction \n";
	  
	    /* the volume of fluid field is advected and used to provide a */
	    /* mass conserving correction to the level set field           */
	    
	    switch (order_of_updates)
	    {
	    case 1:
		  /* the order is x1, x2, x3 */

		  std::cerr<<"initial mass "<< compute_mass_in_domain(
    	   							volume_of_fluid,		
	  							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
	  						mesh_width_x1, mesh_width_x2, mesh_width_x3)<< " \n";

		  
		  /* update the volume of fluid field with the flux in x1 direction */
		  
		  

		  update_volume_of_fluid_x1_flux(level_set_old, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						  volume_of_fluid, vof_after_x1_update, u_1_velocity_new,
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
						      actual_time_step_level_set, mesh_width_x1, lower_bound_derivatives);
		  
		  /* compute the corresponding first intermediate level-set field */
		  /* the field level_set_after_x1_update corresponds exactly to   */
		  /* vof_after_x1_update					  */
		  /* the OLD level-set field is used as starting value for the    */
		  /* iterative process						  */

		  
		  if( match_level_set_to_volume_of_fluid(level_set_old, vof_after_x1_update,
						     level_set_after_x1_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k,			
				  volume_of_fluid_tolerance, lower_bound_derivatives,		
				    number_vof_2_level_set_iterations, number_iterations_ridder,			
				      vof_2_level_set_tolerance,
				      	 mesh_width_x1, mesh_width_x2, mesh_width_x3))
		  {
		  	std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
		  	std::cerr<<" make_level_set_mass_conserving line 202 \n";
		  	exit(1);
		  }
		  

		  /* compute the gradient of the first intermediate level-set field */
		  /* for evaluation of the flux in the x2 direction */
      
		  compute_level_set_gradient(level_set_after_x1_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );

	
 		  /* update the volume of fluid field with the flux in x2 direction */

 
		  update_volume_of_fluid_x2_flux(level_set_after_x1_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						  vof_after_x1_update, vof_after_x2_update, u_2_velocity_new,
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
						      actual_time_step_level_set, mesh_width_x2, lower_bound_derivatives);
	

		  
		  /* compute the corresponding second intermediate level-set field  	*/
		  /* the field level_set_after_x2_update corresponds exactly to     	*/
		  /* vof_after_x2_update					   	   	*/
		  /* the NEW level-set field is used as starting value for the     	*/
		  /* iterative process						   	*/
		  
		  if( match_level_set_to_volume_of_fluid(level_set_after_x1_update, vof_after_x2_update,
						     level_set_after_x2_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k,			
				  volume_of_fluid_tolerance, lower_bound_derivatives,		
				    number_vof_2_level_set_iterations, number_iterations_ridder,			
				      vof_2_level_set_tolerance,
				      	 mesh_width_x1, mesh_width_x2, mesh_width_x3))
		  { 
		  	std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
		  	std::cerr<<" make_level_set_mass_conserving line 242 \n";
		  	exit(1);
		  }
		  

		  /* compute the gradient of the second intermediate level-set field */
		  /* for evaluation of the flux in the x3 direction */
     
		  compute_level_set_gradient(level_set_after_x2_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );


 		  /* update the volume of fluid field with the flux in x2 direction */
 
		  
 		  update_volume_of_fluid_x3_flux(level_set_after_x2_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						  vof_after_x2_update, vof_after_x3_update, u_3_velocity_new,
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
						      actual_time_step_level_set, mesh_width_x3, lower_bound_derivatives);
 	    

		  break;

	    case 2:
		  /* the order is x2, x3, x1 */
		  
		  /* update the volume of fluid field with the flux in x2 direction */
		  update_volume_of_fluid_x2_flux(level_set_old, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						   volume_of_fluid, vof_after_x2_update, u_2_velocity_new,
						     number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
							actual_time_step_level_set, mesh_width_x2, lower_bound_derivatives);
		  
		  /* compute the corresponding first intermediate level-set field */
		  /* the field level_set_after_x2_update corresponds exactly to   */
		  /* vof_after_x2_update					  */
		  /* the OLD level-set field is used as starting value for the    */
		  /* iterative process						  */
		  
		  if( match_level_set_to_volume_of_fluid(level_set_old, vof_after_x2_update,
						     level_set_after_x2_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k,			
				  volume_of_fluid_tolerance, lower_bound_derivatives,		
				    number_vof_2_level_set_iterations, number_iterations_ridder,			
				      vof_2_level_set_tolerance,
				      	 mesh_width_x1, mesh_width_x2, mesh_width_x3))
		  { 
		  	std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
		  	std::cerr<<" make_level_set_mass_conserving line 293 \n";
		  	exit(1);
		  }
		  
		  

		  /* compute the gradient of the first intermediate level-set field */
		  /* for evaluation of the flux in the x3 direction */
      
		  compute_level_set_gradient(level_set_after_x2_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );
	
 		  /* update the volume of fluid field with the flux in x3 direction */
 
		  update_volume_of_fluid_x3_flux(level_set_after_x2_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						    vof_after_x2_update, vof_after_x3_update, u_3_velocity_new,
						      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
							actual_time_step_level_set, mesh_width_x3, lower_bound_derivatives);
	
		  /* compute the corresponding second intermediate level-set field */
		  /* the field level_set_after_x3_update corresponds exactly to    */
		  /* vof_after_x3_update					   */
		  /* the NEW level-set field is used as starting value for the     */
		  /* iterative process						   */
		  
		  if( match_level_set_to_volume_of_fluid(level_set_after_x2_update, vof_after_x3_update,
						     level_set_after_x3_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k,			
				  volume_of_fluid_tolerance, lower_bound_derivatives,		
				    number_vof_2_level_set_iterations, number_iterations_ridder,			
				      vof_2_level_set_tolerance,
				      	 mesh_width_x1, mesh_width_x2, mesh_width_x3))
		  { 
		  	std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
		  	std::cerr<<" make_level_set_mass_conserving line 330 \n";
		  	exit(1);
		  }

		  /* compute the gradient of the second intermediate level-set field */
		  /* for evaluation of the flux in the x1 direction */
     
		  compute_level_set_gradient(level_set_after_x3_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );

 		  /* update the volume of fluid field with the flux in x1 direction */
 
		  update_volume_of_fluid_x1_flux(level_set_after_x3_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						  vof_after_x3_update, vof_after_x1_update, u_1_velocity_new,
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
						      actual_time_step_level_set, mesh_width_x1, lower_bound_derivatives);
 	    
		
		  break;
	    case 3:
		  /* the order is x3, x1, x2 */
		  
		  /* update the volume of fluid field with the flux in x3 direction */
		  update_volume_of_fluid_x3_flux(level_set_old, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						  volume_of_fluid, vof_after_x3_update, u_2_velocity_new,
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
						      actual_time_step_level_set, mesh_width_x3, lower_bound_derivatives);
		  
		  /* compute the corresponding first intermediate level-set field */
		  /* the field level_set_after_x3_update corresponds exactly to   */
		  /* vof_after_x3_update					  */
		  /* the OLD level-set field is used as starting value for the    */
		  /* iterative process						  */
		  
		  if( match_level_set_to_volume_of_fluid(level_set_old, vof_after_x3_update,
						     level_set_after_x3_update,
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				  				volume_of_fluid_tolerance, lower_bound_derivatives,		
				    					number_vof_2_level_set_iterations, number_iterations_ridder,			
				      						vof_2_level_set_tolerance,
				      	 mesh_width_x1, mesh_width_x2, mesh_width_x3))
		  { 
		  	std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
		  	std::cerr<<" make_level_set_mass_conserving line 377 \n";
		  	exit(1);
		  }
		  

		  /* compute the gradient of the first intermediate level-set field */
		  /* for evaluation of the flux in the x1 direction */
      
		  compute_level_set_gradient(level_set_after_x3_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );
	
 		  /* update the volume of fluid field with the flux in x1 direction */
 
		  update_volume_of_fluid_x1_flux(level_set_after_x3_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						    vof_after_x3_update, vof_after_x1_update, u_1_velocity_new,
						      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
							actual_time_step_level_set, mesh_width_x1, lower_bound_derivatives);
	
		  /* compute the corresponding second intermediate level-set field */
		  /* the field level_set_after_x1_update corresponds exactly to    */
		  /* vof_after_x1_update					   */
		  /* the NEW level-set field is used as starting value for the     */
		  
		  if( match_level_set_to_volume_of_fluid(level_set_after_x3_update, vof_after_x1_update,
						     level_set_after_x1_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k,			
				  volume_of_fluid_tolerance, lower_bound_derivatives,		
				    number_vof_2_level_set_iterations, number_iterations_ridder,			
				      vof_2_level_set_tolerance,
				      	 mesh_width_x1, mesh_width_x2, mesh_width_x3))
		  
		  { 
		  	std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
		  	std::cerr<<" make_level_set_mass_conserving line 411 \n";
		  	exit(1);
		  }
		  
		  

		  /* compute the gradient of the second intermediate level-set field */
		  /* for evaluation of the flux in the x2 direction */
     
		  compute_level_set_gradient(level_set_after_x1_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );

 		  /* update the volume of fluid field with the flux in x2 direction */
 
		  update_volume_of_fluid_x2_flux(level_set_after_x1_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
						  vof_after_x1_update, vof_after_x2_update, u_2_velocity_new,
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	
						      actual_time_step_level_set, mesh_width_x2, lower_bound_derivatives);
 	    
		
		  
		  break;
	    case 4:
		  /* the order is x2, x1, x3 */
		  /*left out for the moment */
		  break;
	    case 5:
		  /* the order is x1, x3, x2 */
		  /*left out for the moment */
	    
		  break;
	    case 6:
		  /* the order is x3, x2, x1 */
		  /*left out for the moment */
		
		  break;
	    }  
/* use the different fluxes for the volume of fluid field to  compute */
/* the volume of fluid field at the new time level 		      */

            
	    switch (order_of_updates)
	    {
	    case 1:
		  /* the order is x1, x2, x3 */
		  /* Note index 1 = first nonvirtual cell and  number_primary_cells */
		  /* the last nonvirtual cell					    */
		for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			{
			    volume_of_fluid[i_index][j_index][k_index]=
				vof_after_x3_update[i_index][j_index][k_index]
				-
				actual_time_step_level_set*(
				vof_after_x1_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index-1][j_index][k_index])*one_over_dx1 +
				vof_after_x2_update[i_index][j_index][k_index]*
				  (u_2_velocity_new[i_index][j_index][k_index]-
				    u_2_velocity_new[i_index][j_index-1][k_index])*one_over_dx2 +
				vof_after_x3_update[i_index][j_index][k_index]*
				  (u_3_velocity_new[i_index][j_index][k_index]-
				    u_3_velocity_new[i_index][j_index][k_index-1])*one_over_dx3)
				;
			}
		    }	
		}
		
		/* update virtual cells: apply the boundary condition */
		
// 		std::cerr<<" mass after total update "<< compute_mass_in_domain(
//     	   							volume_of_fluid,		
// 	  							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
// 	  						mesh_width_x1, mesh_width_x2, mesh_width_x3)<< " \n";
		
		field_neumann_boundary(volume_of_fluid, 
				 number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k);
		break;
	    case 2:
		  /* the order is x2, x3, x1 */
		  /* Note index 1 = first nonvirtual cell and  number_primary_cells */
		  /* the last nonvirtual cell					    */
		for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			{
			    volume_of_fluid[i_index][j_index][k_index]=
				vof_after_x1_update[i_index][j_index][k_index]-
				actual_time_step_level_set*(
				vof_after_x1_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index-1][j_index][k_index])*one_over_dx1 +
				vof_after_x2_update[i_index][j_index][k_index]*
				  (u_2_velocity_new[i_index][j_index][k_index]-
				    u_2_velocity_new[i_index][j_index-1][k_index])*one_over_dx2 +
				vof_after_x3_update[i_index][j_index][k_index]*
				  (u_3_velocity_new[i_index][j_index][k_index]-
				    u_3_velocity_new[i_index][j_index][k_index-1])*one_over_dx3);
			}
		    }	
		}
	    
		/* update virtual cells: apply the boundary condition */
		
		field_neumann_boundary(volume_of_fluid, 
				 number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k);
		break;
	    case 3:
		  /* the order is x3, x1, x2 */
		  /* Note index 1 = first nonvirtual cell and  number_primary_cells */
		  /* the last nonvirtual cell					    */
		for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
			{
			    volume_of_fluid[i_index][j_index][k_index]=
				vof_after_x2_update[i_index][j_index][k_index]-
				actual_time_step_level_set*(
				vof_after_x1_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index-1][j_index][k_index])*one_over_dx1 +
				vof_after_x2_update[i_index][j_index][k_index]*
				  (u_2_velocity_new[i_index][j_index][k_index]-
				    u_2_velocity_new[i_index][j_index-1][k_index])*one_over_dx2 +
				vof_after_x3_update[i_index][j_index][k_index]*
				  (u_3_velocity_new[i_index][j_index][k_index]-
				    u_3_velocity_new[i_index][j_index][k_index-1])*one_over_dx3);
			}
		    }	
		}
	    
		field_neumann_boundary(volume_of_fluid, 
				 number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k);
		
		  break;
	    }  
	    

	    /* check if there are any nonvalid volume of fluid values in the cells */
 
                analyse_validity_vof_field( level_set_star, volume_of_fluid, invalid_vof_cells ,                               
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                    
                                          volume_of_fluid_tolerance, number_cells_vof_out_of_bounds, number_cells_numerical_vapor,                             
                                                number_cells_invalid_volume_of_fluid);
              
                std::cerr<<" number_cells_invalid_volume_of_fluid in vof "<< number_cells_invalid_volume_of_fluid<<" \n";
                std::cerr<<" number_cells_numerical_vapor in vof "<< number_cells_numerical_vapor<<" \n";
                std::cerr<<" number_cells_vof_out_of_bounds in vof "<< number_cells_vof_out_of_bounds<<" \n";
             
             if(number_cells_invalid_volume_of_fluid>0)
             {
	    
                            std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
                            std::cerr<<"invalid VOF cels BEFORE redistribution: "<< number_cells_invalid_volume_of_fluid <<" \n";
                            std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n";

	              if(apply_mass_distribution_algorithm)

   /* the mass redistribution algorithm is applied to circumvent the occurence  */ 
   /* of flotsam and jetsam and volume of fluid values that are outside the     */
   /* physically realizable range [0,1] within the function the level-set field */
   /* is adapted to conform to the new volume of fluid field                    */

	               {
		              apply_volume_of_fluid_redistribution(volume_of_fluid, level_set_star, level_set_new, 				       
						      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,				
							mesh_width_x1, mesh_width_x2, mesh_width_x3,  
							 volume_of_fluid_tolerance, vof_2_level_set_tolerance, 
							    lower_bound_derivatives, number_vof_2_level_set_iterations,		
							      number_iterations_ridder, time_step_mass_redistribution,		
	   							redistribution_vof_tolerance, maximum_number_mass_redistribution_iterations, 
                                                                  mass_redistribution_diffusion_coefficient);
	      
                              // cells_out_of_bounds=check_volume_of_fluid(volume_of_fluid, number_primary_cells_i,
                                           // number_primary_cells_j,
                                            // number_primary_cells_k,
                                               // volume_of_fluid_tolerance);
                                               
                              analyse_validity_vof_field( level_set_new, volume_of_fluid, invalid_vof_cells ,                               
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                    
                                          volume_of_fluid_tolerance, number_cells_vof_out_of_bounds, number_cells_numerical_vapor,                             
                                                number_cells_invalid_volume_of_fluid);
                              
                              std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
                              std::cerr<<"invalid VOF cells  AFTER redistribution: "<< number_cells_invalid_volume_of_fluid <<" \n";
                              std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	               }
	               else
	               {
                              
/* the volume of fluid values are clipped, such that they are inside 				  */
/* the physically realizable range [0,1]              						  */      
/* the array vof_after_x1_update is used here as a work array for temporary storage 		  */


                            number_clipped_cells=apply_volume_of_fluid_clipping(volume_of_fluid,
                                                    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                                                        volume_of_fluid_tolerance);
	                     if(number_clipped_cells>0)
                            {
                                   
                                  std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	    		          std::cerr<<"WARNING: clipping algorithm is applied \n";
                                  std::cerr<<"to "<< number_clipped_cells<<" cells .\n";
                                  std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
			 
		              }
		              else
                            {
                                std::cerr<<"No clipping is required \n";
                            }
		
                      }
              }
		            /* compute the corresponding level-set field 			   */

		              if( match_level_set_to_volume_of_fluid(level_set_star, volume_of_fluid, level_set_new,
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				  			  volume_of_fluid_tolerance, lower_bound_derivatives,		
				                          number_vof_2_level_set_iterations, number_iterations_ridder,			
				      				vof_2_level_set_tolerance,
				      	 mesh_width_x1, mesh_width_x2, mesh_width_x3 ))
			      { 
			      	      std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
			      	      std::cerr<<" make_level_set_mass_conserving line 636 \n";
			      	      exit(1);
			      }
		  

		
		              /* compute the gradient of the old level-set field, necessary for the level-set/vof conversion */

		              compute_level_set_gradient(level_set_new, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				      );
                              
		              /* compute the current volume of fluid field */
	
		              if(compute_volume_of_fluid(level_set_new, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				                            volume_of_fluid,
				                                   number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,					  
				                                   lower_bound_derivatives))
		              {
	  		              std::cerr<<" conversion from volume of fluid to level set unsuccesfull\n";
			              std::cerr<<" in make_level_set_mass_conserving line 784 \n";
			              exit(1); 
		              };
		
	                     std::cerr<<" mass after final conversion "<<
	                                   compute_mass_in_domain( volume_of_fluid,
	  							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
	  						mesh_width_x1, mesh_width_x2, mesh_width_x3)<< " \n";

	}

	else
	{	
	    std::cerr<<"No mass conservation applied: copy star to new time level \n";

	    /* NO mass conserving correction to the level set field is applied */
	    /* the new volume of fluid field is simply evaluated from the new  */
	    /* not necessarily mass conserving level-set field 		       */

	    copy_cell_centered_field(level_set_star, level_set_new, number_primary_cells_i,
					    number_primary_cells_j, number_primary_cells_k);

		
	}    
	
	
	
/* if a variable ordering of the updates is applied, then use a different order on the next entry */
/* order_of_updates is a STATIC variable, so its value is preserved */

	order_of_updates+=1;
	
/* deallocate memory for the temporary storage of the intermediate fields */

	vof_after_x1_update.destroy();
	vof_after_x2_update.destroy();
	vof_after_x3_update.destroy();
	level_set_after_x1_update.destroy();
	level_set_after_x2_update.destroy();
	level_set_after_x3_update.destroy();
	d_level_set_d_x1.destroy();
	d_level_set_d_x2.destroy();
	d_level_set_d_x3.destroy();
	invalid_vof_cells.destroy();
    }

 
 
      

   

#include "../headers/array.h"

#include<cstdlib>
#include<iostream>
/********************************************************************************/
/********************************************************************************/
/*  Function to advect the volume of fluid field                                */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The volume of fluid field is advected using an operator splitting technique  */
/* that is explained in both Sander's thesis as in the papers of Sussman        */
/* and is first order accurate. The contributions of fluxes in the different    */
/* directions are added one after the other. To circumvent any directional      */
/* dependence, to order in which the fluxes are added is changed every time     */
/* with 3 spatial directions there are 3!=6 possible orderings.                 */
/* Six orderings are defined, and after ordering 6 is applied the process starts*/
/* over again.									*/
/* Note that after each flux computation, the changes in the intermediate       */
/* volume of fluid fields have to be accounted for in changes in corresponding  */
/* level set fields.								*/
/********************************************************************************/
    void volume_of_fluid_advection_main_control
    (
      Array3<double> level_set_star, 		// level set field at new time level
						// after convection and reinitialization
						// not mass conserving
      Array3<double> level_set_new, 			// level set field at new time level
						// mass conserving
      Array3<double> level_set_old, 			// level set field at old time level
						// mass conserving
      Array3<double> volume_of_fluid, 		// volume of fluid field
      Array3<double> u_1_velocity_new, 		// velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new, 		// velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,		// velocity field at new time level x3 direction
      int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
      double actual_time_step_level_set,	// time step used for level-set advection
						// computed from all stability restrictions and 
						// possibly subscycling
      double mesh_width_x1,			// grid spacing in x1 direction (uniform)
      double mesh_width_x2,			// grid spacing in x2 direction (uniform)
      double mesh_width_x3,			// grid spacing in x3 direction (uniform)
      int apply_mass_distribution_algorithm,    // =1:apply the mass redistribution algorithm
						// to avoid numerical vapour
      int apply_mass_conservation_correction    // =1:apply the mass conservation correction
						// to the level-set field
    
    )
    {
      /* function definitions */
      
      void 	compute_level_set_gradient(				// compute gradient of level-set field
		Array3<double> level_set_star, Array3<double> d_level_set_d_x1, 
		Array3<double> d_level_set_d_x2, Array3<double> d_level_set_d_x3,
		int number_primary_cells_i, int number_primary_cells_j, 
		int number_primary_cells_k);
      int    compute_volume_of_fluid(					// compute volume of fluid field
		 Array3<double> level_set_star, 				// corresponding to a given 
		 Array3<double> d_level_set_d_x1, 				// level-set field
		 Array3<double> d_level_set_d_x2, 
		 Array3<double> d_level_set_d_x3,
		 Array3<double> volume_of_fluid,
		 int number_primary_cells_i, 
		 int number_primary_cells_j, 
		 int number_primary_cells_k);
      int   apply_volume_of_fluid_clipping(				// apply simple clipping
		 Array3<double> volume_of_fluid, 				// to keep the volume of fluid
		 int number_primary_cells_i, 				// field in the interval [0,1] 
		 int number_primary_cells_j,
		 int number_primary_cells_k);                           
      void update_volume_of_fluid_x1_flux(				// compute intermediate 
		 Array3<double> level_set, 					// volume of fluid field 
		 Array3<double> d_level_set_d_x1,				// updated by the flux in x1
		 Array3<double> d_level_set_d_x2, 				// direction
		 Array3<double> d_level_set_d_x3,
	         Array3<double> volume_of_fluid, 
		 Array3<double> vof_after_x1_update,
		 int number_primary_cells_i, int number_primary_cells_j, 
		 int number_primary_cells_k);
      void update_volume_of_fluid_x2_flux(				// compute intermediate
		 Array3<double> level_set, 					// volume of fluid field 
		 Array3<double> d_level_set_d_x1,				// updated by the flux in x2
		 Array3<double> d_level_set_d_x2, 				// direction
		 Array3<double> d_level_set_d_x3,
	         Array3<double> volume_of_fluid, 
		 Array3<double> vof_after_x2_update,
		 int number_primary_cells_i, int number_primary_cells_j, 
		 int number_primary_cells_k);
      void update_volume_of_fluid_x3_flux(				// compute intermediate
		 Array3<double> level_set, 					// volume of fluid field 
		 Array3<double> d_level_set_d_x1,				// updated by the flux in x3
		 Array3<double> d_level_set_d_x2,  				// direction
		 Array3<double> d_level_set_d_x3,
	         Array3<double> volume_of_fluid, 
		 Array3<double> vof_after_x3_update,
		 int number_primary_cells_i, int number_primary_cells_j, 
		 int number_primary_cells_k);
      void match_level_set_to_volume_of_fluid(				// compute a new level-set field
		 Array3<double> level_set_not_mass_conserving, 		// that corresponds to given
		 Array3<double> volume_of_fluid,				// volume of fluid field
		 Array3<double> level_set_mass_conserving,			// using an existing level-set
		 int number_primary_cells_i, 				// field as a starting value
		 int number_primary_cells_j, 
		 int number_primary_cells_k,			
		 double volume_of_fluid_tolerance,		
		 double lower_bound_derivatives,		
		 int number_vof_2_level_set_iterations,		
		 int number_iterations_ridder,			
		 double vof_2_level_set_tolerance		
		);
      void apply_mass_distribution(					// apply mass distribution
		 Array3<double> level_set_star, 				// algorithm to find an 
		 Array3<double> volume_of_fluid, 				// updated level-set and 
		 Array3<double> level_set_new,				// volume of fluid field	
		 int number_primary_cells_i, int number_primary_cells_j, 
		 int number_primary_cells_k);
      void copy_cell_centered_field(					// copy a source field
		 Array3<double> source_field, 				// to a target field 
		 Array3<double> target_field, 				
		 int number_primary_cells_i, int number_primary_cells_j,
		 int number_primary_cells_k);

      
      
      Array3<double> vof_after_x1_update;			// intermediate volume of fluid field, 
							// after update in x1 direction
      Array3<double> vof_after_x2_update;			// intermediate volume of fluid field, 
							// after update in x2 direction
      Array3<double> vof_after_x3_update;			// intermediate volume of fluid field, 
							// after update in x3 direction
      Array3<double> d_level_set_d_x1;			// first partial derivative of
							// the level-set field wrt x1
							// second order central approximation
      Array3<double> d_level_set_d_x2;			// first partial derivative of 
							// the level-set field wrt x2
							// second order central approximation
      Array3<double> d_level_set_d_x3;			// first partial derivative of
 							// the level-set field wrt x3
							// second order central approximation
      Array3<double> level_set_after_x1_update;		// level-set field, converted from 
							// intermediate	volume of fluid field, 
							// after update in x1 direction
      Array3<double> level_set_after_x2_update;		// level-set field, converted from 
							// intermediate	volume of fluid field, 
							// after update in x2 direction
      Array3<double> level_set_after_x3_update;		// level-set field, converted from 
 							// intermediate	volume of fluid field, 
							// after update in x3 direction
     
      static int order_of_updates;			// indicates the order of the updates
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
      
      int i_index, j_index, k_index;  			// local variables for loop indexing
      int number_of_clipped_cells=0;			// the number of cells where clipping was
      double one_over_dx1	=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      double one_over_dx2	=    			// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
      double one_over_dx3	=    			// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);
      
/* allocate memory for the temporary storage of the intermediate fields */

	vof_after_x1_update.create(number_primary_cells_i, number_primary_cells_j,
					   number_primary_cells_k);
	vof_after_x2_update.create(number_primary_cells_i, number_primary_cells_j,
					   number_primary_cells_k);
	vof_after_x3_update.create(number_primary_cells_i, number_primary_cells_j,
					   number_primary_cells_k);
	
	level_set_after_x1_update.create(number_primary_cells_i, number_primary_cells_j,
				    number_primary_cells_k);
	level_set_after_x2_update.create(number_primary_cells_i, number_primary_cells_j,
				    number_primary_cells_k);
	level_set_after_x3_update.create(number_primary_cells_i, number_primary_cells_j,
				    number_primary_cells_k);
	
	d_level_set_d_x1.create(number_primary_cells_i, number_primary_cells_j,
				    number_primary_cells_k);
	d_level_set_d_x2.create(number_primary_cells_i, number_primary_cells_j,
				    number_primary_cells_k);
	d_level_set_d_x3.create(number_primary_cells_i, number_primary_cells_j,
				    number_primary_cells_k);
				    
    /* keep order_of_updates in its allowable range */
    
      if( order_of_updates > 6)
      { 	
	order_of_updates=1;
      }
      
    /* Sander chose to keep the order_of_updates fixed to 1 */    
    
	order_of_updates=1;
	
	
    /* deallocate the used memory */
    
	vof_after_x1_update.destroy();
	vof_after_x2_update.destroy();
	vof_after_x3_update.destroy();
	level_set_after_x1_update.destroy();
	level_set_after_x2_update.destroy();
	level_set_after_x3_update.destroy();

	
	/* compute the gradient of the level-set field, necessary for the level-set/vof conversion */
	compute_level_set_gradient(level_set_star, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );
	/* compute the current volume of fluid field */
	
	if(!compute_volume_of_fluid(level_set_star, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				volume_of_fluid,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  ))
	{
	   cout<< "something went wrong in the level-set to volume of fluid conversion, aborting..."
	   exit(1);
	};
	
	if(apply_mass_conservation_correction)
	{
	  
	    /* the volume of fluid field is advected and used to provide a */
	    /* mass conserving correction to the level set field */
	    
	    switch (order_of_updates)
	    {
	    case 1:
		  /* the order is x1, x2, x3 */
		  
		  /* update the volume of fluid field with the flux in x1 direction */
		  update_volume_of_fluid_x1_flux(level_set_old, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   volume_of_fluid, vof_after_x1_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
		  
		  /* compute the corresponding first intermediate level-set field */
		  /* the field level_set_after_x1_update corresponds exactly to   */
		  /* vof_after_x1_update					  */
		  /* the OLD level-set field is used as starting value for the    */
		  /* iterative process						  */
		  
		  match_level_set_to_volume_of_fluid(level_set_old, vof_after_x1_update,
						     level_set_after_x1_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);

		  /* compute the gradient of the first intermediate level-set field */
		  /* for evaluation of the flux in the x2 direction */
      
		  compute_level_set_gradient(level_set_after_x1_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );
	
 		  /* update the volume of fluid field with the flux in x2 direction */
 
		  update_volume_of_fluid_x2_flux(level_set_after_x1_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   vof_after_x1_update, vof_after_x2_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
	
		  /* compute the corresponding second intermediate level-set field */
		  /* the field level_set_after_x2_update corresponds exactly to    */
		  /* vof_after_x2_update					   */
		  /* the NEW level-set field is used as starting value for the     */
		  /* iterative process						   */
		  
		  match_level_set_to_volume_of_fluid(level_set_star, vof_after_x2_update,
						     level_set_after_x2_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);

		  /* compute the gradient of the second intermediate level-set field */
		  /* for evaluation of the flux in the x3 direction */
     
		  compute_level_set_gradient(level_set_after_x2_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );

 		  /* update the volume of fluid field with the flux in x2 direction */
 
		  update_volume_of_fluid_x3_flux(level_set_after_x2_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   vof_after_x2_update, vof_after_x3_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
 	    
		  break;
	    case 2:
		  /* the order is x2, x3, x1 */
		  
		  /* update the volume of fluid field with the flux in x2 direction */
		  update_volume_of_fluid_x2_flux(level_set_old, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   volume_of_fluid, vof_after_x2_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
		  
		  /* compute the corresponding first intermediate level-set field */
		  /* the field level_set_after_x2_update corresponds exactly to   */
		  /* vof_after_x2_update					  */
		  /* the OLD level-set field is used as starting value for the    */
		  /* iterative process						  */
		  
		  match_level_set_to_volume_of_fluid(level_set_old, vof_after_x2_update,
						     level_set_after_x2_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);

		  /* compute the gradient of the first intermediate level-set field */
		  /* for evaluation of the flux in the x3 direction */
      
		  compute_level_set_gradient(level_set_after_x2_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );
	
 		  /* update the volume of fluid field with the flux in x3 direction */
 
		  update_volume_of_fluid_x3_flux(level_set_after_x2_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   vof_after_x2_update, vof_after_x3_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
	
		  /* compute the corresponding second intermediate level-set field */
		  /* the field level_set_after_x3_update corresponds exactly to    */
		  /* vof_after_x3_update					   */
		  /* the NEW level-set field is used as starting value for the     */
		  /* iterative process						   */
		  
		  match_level_set_to_volume_of_fluid(level_set_star, vof_after_x3_update,
						     level_set_after_x3_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);

		  /* compute the gradient of the second intermediate level-set field */
		  /* for evaluation of the flux in the x1 direction */
     
		  compute_level_set_gradient(level_set_after_x3_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );

 		  /* update the volume of fluid field with the flux in x1 direction */
 
		  update_volume_of_fluid_x1_flux(level_set_after_x3_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   vof_after_x3_update, vof_after_x1_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
 	    
		
		  break;
	    case 3:
		  /* the order is x3, x1, x2 */
		  
		  /* update the volume of fluid field with the flux in x3 direction */
		  update_volume_of_fluid_x3_flux(level_set_old, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   volume_of_fluid, vof_after_x3_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
		  
		  /* compute the corresponding first intermediate level-set field */
		  /* the field level_set_after_x3_update corresponds exactly to   */
		  /* vof_after_x3_update					  */
		  /* the OLD level-set field is used as starting value for the    */
		  /* iterative process						  */
		  
		  match_level_set_to_volume_of_fluid(level_set_old, vof_after_x3_update,
						     level_set_after_x3_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);

		  /* compute the gradient of the first intermediate level-set field */
		  /* for evaluation of the flux in the x1 direction */
      
		  compute_level_set_gradient(level_set_after_x3_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );
	
 		  /* update the volume of fluid field with the flux in x1 direction */
 
		  update_volume_of_fluid_x1_flux(level_set_after_x3_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   vof_after_x3_update, vof_after_x1_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
	
		  /* compute the corresponding second intermediate level-set field */
		  /* the field level_set_after_x1_update corresponds exactly to    */
		  /* vof_after_x1_update					   */
		  /* the NEW level-set field is used as starting value for the     */
		  
		  match_level_set_to_volume_of_fluid(level_set_star, vof_after_x1_update,
						     level_set_after_x1_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);

		  /* compute the gradient of the second intermediate level-set field */
		  /* for evaluation of the flux in the x2 direction */
     
		  compute_level_set_gradient(level_set_after_x1_update, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );

 		  /* update the volume of fluid field with the flux in x2 direction */
 
		  update_volume_of_fluid_x2_flux(level_set_after_x1_update, d_level_set_d_x1,
						 d_level_set_d_x2, d_level_set_d_x3,
				   vof_after_x1_update, vof_after_x2_update,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
 	    
		
		  
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
				vof_after_x3_update[i_index][j_index][k_index]-
				actual_time_step_level_set*(
				vof_after_x1_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index-1][j_index][k_index])*one_over_dx1 +
				vof_after_x2_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index][j_index-1][k_index])*one_over_dx2 +
				vof_after_x3_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index][j_index][k_index-1])*one_over_dx3);
			}
		    }	
		}
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
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index][j_index-1][k_index])*one_over_dx2 +
				vof_after_x3_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index][j_index][k_index-1])*one_over_dx3);
			}
		    }	
		}
	    
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
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index][j_index-1][k_index])*one_over_dx2 +
				vof_after_x3_update[i_index][j_index][k_index]*
				  (u_1_velocity_new[i_index][j_index][k_index]-
				    u_1_velocity_new[i_index][j_index][k_index-1])*one_over_dx3);
			}
		    }	
		}
	    
		
		  break;
	    }  

	    

	    if(apply_mass_distribution_algorithm)

/* the mass redistribution algorithm is applied to circumvent the occurence of flotsam and jetsam */
/* and volume of fluid values that are outside the physically realizable range [0,1]              */
/* within the function the level-set field is adapted to conform to the new volume of fluid field */

	    {
		apply_mass_distribution(level_set_star, volume_of_fluid, level_set_new,
		    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
	      
	    }
	    else
	    {
/* the volume of fluid values are clipped, such that they are inside 				  */
/* the physically realizable range [0,1]              						  */      
/* the array vof_after_x1_update is used here as a work array for temporary storage 		  */
	      
	        number_of_clipped_cells=apply_volume_of_fluid_clipping(volume_of_fluid, 
		    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
		
		  /* compute the corresponding level-set field 			   */
		  /* the last intermediate */
		match_level_set_to_volume_of_fluid(level_set_star, volume_of_fluid,
						     level_set_new,
				number_primary_cells_i, number_primary_cells_j, 
							    number_primary_cells_k	
						);
	    }
	}

	else
	{
	    /* NO mass conserving correction to the level set field is applied */
	    /* the new volume of fluid field is simply evaluated from the new  */
	    /* not necessarily mass conserving level-set field 		       */
	    copy_cell_centered_field(level_set_star, level_set_new, number_primary_cells_i,
					    number_primary_cells_j, number_primary_cells_k);
	}    
	
	
	
/* if a variable ordering of the updates is applied, then use a different order on the next entry */
/* order_of_updates is a STATIC variable, so its value is preserved */

	order_of_updates+=1;
	
/* allocate memory for the temporary storage of the intermediate fields */

	vof_after_x1_update.destroy();
	vof_after_x2_update.destroy();
	vof_after_x3_update.destroy();
	level_set_after_x1_update.destroy();
	level_set_after_x2_update.destroy();
	level_set_after_x3_update.destroy();
	d_level_set_d_x1.destroy();
	d_level_set_d_x2.destroy();
	d_level_set_d_x3.destroy();
    }

 
 
      

   

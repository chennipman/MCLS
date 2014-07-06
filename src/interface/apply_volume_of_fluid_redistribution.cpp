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
/* and is first order accurate. This does not guarantee that the volume of      */
/* fluid value remains within the interval [0,1]. 	 			*/
/* Furthermore, it can occur that a cell is too large to represent a small      */
/* feature. This leads to cells with intermediate values of the volume of 	*/
/* fluid in cells that are NOT interface cells. Of course this can not be 	*/
/* allowed. Both effects:                                                       */
/* -cells that are overfilled or underdrained					*/
/* -cells that contain numerical vapor                                          */
/* Have to handled in a way that does not disturb the global mass conservation. */
/* level set fields.								*/
/********************************************************************************/
EXPORT void apply_volume_of_fluid_redistribution(
	    Array3<double> volume_of_fluid, 				// volume of fluid field
	    Array3<double> level_set_star, 				// level set field at new time level
									// after convection and reinitialization
									// not mass conserving
	    Array3<double> level_set_new, 				// level set field at new time level
									// mass conserving
	    int number_primary_cells_i,				        // number of primary (pressure) cells in x1 direction
	    int number_primary_cells_j,				        // number of primary (pressure) cells in x2 direction
	    int number_primary_cells_k,				        // number of primary (pressure) cells in x3 direction
	    double mesh_width_x1,					// grid spacing in x1 direction (uniform)
	    double mesh_width_x2,					// grid spacing in x2 direction (uniform)
	    double mesh_width_x3,					// grid spacing in x3 direction (uniform)
	    double volume_of_fluid_tolerance,			        // tolerance on the value of the volume
									// of fluid value
	    double vof_2_level_set_tolerance,			        // tolerance in the conversion from volume
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
      	    int maximum_number_mass_redistribution_iterations, 	        // number of iterations allowed to make
									// the volume of fluid field valid
									// these are the sweeps on the vof error
            double mass_redistribution_diffusion_coefficient            // diffusion coefficient for mass redistribution equation
	    )
      {
      Array3<double> volume_of_fluid_star;				// volume of fluid field, uncorrected
									// so with possible vapour cells and 
									// under/overfilled cells
      Array3<double> volume_of_fluid_correction;			// correction to the volume of fluid field
								        // to make it valid
      Array3<double> invalid_vof_cells;                                 // indication field to show what is wrong
                                                                        // indicator field showing cells that are either
                                                                        // within bounds =0
                                                                        // underfilled   =-1
                                                                        // overfilled    =+1
                                                                        // vapour cells  = 5
      int number_cells_vof_out_of_bounds=100;	                        // number of control volumes where the volume of fluid
						                        // function is OUTSIDE the interval [0,1]
      int number_cells_numerical_vapor=100;		                // number of control volumes where the volume of fluid
						                        // function is INSIDE the interval [0,1]
						                        // while the cell has 6 neighbours with the same sign:
						                        // the cell is NOT an interface cell
      int number_cells_invalid_volume_of_fluid=200; 	                // sum of number of vapour cells and number of cells
						                        // with the volume of fluid outside [0,1];
      int index_redistribution_attempt=0;	                        // number of attempts to achieve a
						                        // valid volume of fluid field
						                        // through the redistribution algorithm
						
      double initial_mass;                                              // the mass present in the computational
                                                                        // domain initial value
      double final_mass;                                                // the mass present in the computational
                                                                        // domain final value
      double change_in_mass;                                            // the change in mass present in the 
                                                                        // domain with respect to the initial mass
      bool detailed_output=0;                                           // =1, provide detailed output
                                                                        // =0, provide non output 
                                                                       
      /* allocate memory for the volume of fluid correction, the tentative    */
      /* volume of fluid field, and the indicator field                       */
      
	volume_of_fluid_correction.create(number_primary_cells_i+2, 
						    number_primary_cells_j+2,
						      number_primary_cells_k+2);
	volume_of_fluid_star.create(number_primary_cells_i+2, 
						    number_primary_cells_j+2,
						      number_primary_cells_k+2);
        invalid_vof_cells.create(number_primary_cells_i+2, 
                                                    number_primary_cells_j+2,
                                                      number_primary_cells_k+2);
      
       /* compute the initial mass in the computational domain */
        
        initial_mass= compute_mass_in_domain(volume_of_fluid, number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,         
                                                        mesh_width_x1, mesh_width_x2, mesh_width_x3);
                  
        std::cerr<<"initial mass before redistribution "<< initial_mass<< " \n";
        
        
        /* do a number of redistribution attempts to move the error to the interface      */
        /* this is an iterative process, because the position of the interface is updated */

        while(index_redistribution_attempt<5 &&
	      number_cells_invalid_volume_of_fluid>0)
	{
	  
		
	

                /* copy the original volume of fluid field to the star volume of fluid field */

	 	copy_cell_centered_field(volume_of_fluid, volume_of_fluid_star, 				
		      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
                /* apply clipping to the tentative volume of fluid field */
                /* in this way it can be safely converted to a level-set field */
      
                if(index_redistribution_attempt<1)
                {
	  	        number_cells_vof_out_of_bounds=
	      		     apply_volume_of_fluid_clipping(volume_of_fluid_star, 
		    	 	number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
		                   volume_of_fluid_tolerance);
                //         make_vof_field_valid(level_set_star, volume_of_fluid_star, invalid_vof_cells,                               
                //                 number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                    
                //                    volume_of_fluid_tolerance, number_cells_vof_out_of_bounds,                           
                //                      number_cells_numerical_vapor, number_cells_invalid_volume_of_fluid);
                }
                else
                {
                	/* this function changes the volume of fluid field in the same way */
                	/* as apply_volume_of_fluid_clipping, but also reports on the      */
                	/* different error cells: overfilled, underfilled and vapour cells */
                	
                        make_vof_field_valid(level_set_new, volume_of_fluid_star, invalid_vof_cells,                               
                                number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                    
                                   volume_of_fluid_tolerance, number_cells_vof_out_of_bounds,                           
                                     number_cells_numerical_vapor, number_cells_invalid_volume_of_fluid);
                }



                /* bring the level-set field in accordance with the clipped volume of fluid field */
                /* the level-set field can only be computed from the clipped field, because only  */
                /* when the volume of fluid is in the closed interval [0,1], the function can be  */
                /* inverted.                                                                      */


                if(index_redistribution_attempt<1)
                {

	  	if(match_level_set_to_volume_of_fluid(level_set_star, volume_of_fluid_star, level_set_new,
								number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				  					volume_of_fluid_tolerance, lower_bound_derivatives,		
				    						number_vof_2_level_set_iterations, number_iterations_ridder,			
				      							vof_2_level_set_tolerance,
                		      	 						mesh_width_x1, mesh_width_x2, mesh_width_x3))	
 		  	{ 
 		  		std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
 		  		std::cerr<<" apply_volume_of_fluid_redistribution line 154 \n";
 		  		exit(1);
		  	}
		  

                
                }
                else
                {
                       
                if(match_level_set_to_volume_of_fluid(level_set_new, volume_of_fluid_star, level_set_new,
                                                                number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                 
                                                                        volume_of_fluid_tolerance, lower_bound_derivatives,             
                                                                                number_vof_2_level_set_iterations, number_iterations_ridder,                    
                                                                                        vof_2_level_set_tolerance,
                                                                                        mesh_width_x1, mesh_width_x2, mesh_width_x3))   
  		  	{ 
  		  		std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
  		  		std::cerr<<" apply_volume_of_fluid_redistribution line 171 \n";
  		  		exit(1);
		  	}
               }
               
                /* extend the level-set field to the virtual cells */
    
                field_extrapolate_boundary(level_set_new, number_primary_cells_i, 
                           number_primary_cells_j,number_primary_cells_k);    
                
                /* The interface location is now fixed by level_set_new */
                /* The criterium to decide if a cell is a vapour cell, is based on this level-set field. */
                /* This is the cause for the need to iterate the algorithm.                      */
                /* The over and undershoots for the original volume of fluid field are computed  */
      
      
	    	number_cells_invalid_volume_of_fluid=			
			modify_volume_of_fluid_values( level_set_new, volume_of_fluid, volume_of_fluid_correction, invalid_vof_cells,
                                                   mesh_width_x1, mesh_width_x2, mesh_width_x3,
						    number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			  
						      volume_of_fluid_tolerance);
		
                if(detailed_output)
                {
                        dump_redistribution_for_debugging(level_set_star, volume_of_fluid,        
                                                    level_set_new, invalid_vof_cells, volume_of_fluid_correction,          
                                                      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                    
                                                        mesh_width_x1, mesh_width_x2, mesh_width_x3, index_redistribution_attempt );  
                }
                
                /* if necessary redistribute the volume of fluid correction */
                /* so we end up with a valid field  */
      
      		if(number_cells_invalid_volume_of_fluid>0)
			
		{
      			redistribute_volume_of_fluid_error(level_set_star, level_set_new, volume_of_fluid, 
                                              volume_of_fluid_correction, invalid_vof_cells,
						number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
							mesh_width_x1,	mesh_width_x1,	mesh_width_x3,			
								time_step_mass_redistribution, volume_of_fluid_tolerance,		
	 								redistribution_vof_tolerance,	
						 				maximum_number_mass_redistribution_iterations,
                                                                                  mass_redistribution_diffusion_coefficient,
                                                                                  index_redistribution_attempt);
     
      
		}
		
		/* compute the change in mass that resulted from this redistribution sweep */
		
                change_in_mass= compute_mass_in_domain(volume_of_fluid, number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,         
                                                        mesh_width_x1, mesh_width_x2, mesh_width_x3)-initial_mass;
                std::cerr<<"change in mass after redistribution step  "<< index_redistribution_attempt<<" "<< change_in_mass<< " \n";
                

                index_redistribution_attempt++;
                        
                


        }

        /* if there were cells to be corrected in the last iteration there is a     */
        /* need to make the last level-set field comply with the last volume of     */
        /* fluid field.                                                             */
     
        analyse_validity_vof_field( level_set_new, volume_of_fluid, invalid_vof_cells,                               
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                                    
                                          volume_of_fluid_tolerance, number_cells_vof_out_of_bounds, number_cells_numerical_vapor,                             
                                                number_cells_invalid_volume_of_fluid);
        
        std::cerr<<"The analysis shows there are "<< number_cells_invalid_volume_of_fluid<< "  invalid cells\n";

      	if(number_cells_invalid_volume_of_fluid>0)
	{
	    	  number_cells_invalid_volume_of_fluid=			
			modify_volume_of_fluid_values(level_set_new, volume_of_fluid, volume_of_fluid_correction, invalid_vof_cells,
                                                       mesh_width_x1, mesh_width_x2, mesh_width_x3,
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			  
								volume_of_fluid_tolerance);

	    	  if(match_level_set_to_volume_of_fluid(level_set_star, volume_of_fluid, level_set_new, 
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
				  			   volume_of_fluid_tolerance, lower_bound_derivatives,		
				    				number_vof_2_level_set_iterations, number_iterations_ridder,			
				      					vof_2_level_set_tolerance,
				      	 			mesh_width_x1, mesh_width_x2, mesh_width_x3))
		  		  { 
		  	std::cerr<<" match_level_set_to_volume_of_fluid was called from \n" ;
		  	std::cerr<<" apply_volume_of_fluid_redistribution line 248 \n";  
		  	exit(1);
		  }
		  

                  if(number_cells_invalid_volume_of_fluid>0)
                  {
		  std::cerr<< "**************************************************** \n";
		  std::cerr<<"ERROR \n";
		  std::cerr<<"Despite the application of the volume of fluid  \n";
		  std::cerr<<"redistribution algorithm, there were still "<<
		  			number_cells_invalid_volume_of_fluid<< "\n";
		  std::cerr<<"cells that have a volume of fluid value      \n";
		  std::cerr<<"outside the physically relevant interval  \n";
                  std::cerr<<"Increase the number of iterations, or     \n";
                  std::cerr<<"check your input.                         \n";
		  std::cerr<<"in apply_volume_of_fluid_redistribution line 241. \n";
		  std::cerr<<"**************************************************** \n";
                  exit(1);
	          
                  }
	    
	}
                
        final_mass=compute_mass_in_domain( volume_of_fluid, number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,         
                                                        mesh_width_x1, mesh_width_x2, mesh_width_x3);
        
        std::cerr<<"final mass after redistribution "<< final_mass << " \n";
	
	/* deallocate temporary storage */
	  
        volume_of_fluid_star.destroy();
        volume_of_fluid_correction.destroy();
        invalid_vof_cells.destroy();
}
  

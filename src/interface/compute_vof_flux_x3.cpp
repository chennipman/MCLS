#include "../headers/array.h"

#include<cstdlib>
#include<iostream>
/********************************************************************************/
/********************************************************************************/
/*  Function to  compute the volume of fluid fluxes in the x3 direction       	*/
/*  direction										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
//
     void compute_vof_flux_x3(						
	Array3<double> level_set, 				// level set field 
							// mass conserving
       Array3<double> u_3_velocity_new, 			// velocity field at new time level x1 direction
	Array3<double> d_level_set_d_x1,			// first partial derivative of
							// the level-set field wrt x1
							// second order central approximation
	Array3<double> d_level_set_d_x2,			// first partial derivative of 
							// the level-set field wrt x2
							// second order central approximation
	Array3<double> d_level_set_d_x3,			// first partial derivative of
 							// the level-set field wrt x3
							// second order central approximation
	Array3<double> flux_x3,				// volume of fluid flux in x1 direction
	int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	double actual_time_step_level_set,		// time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
	double mesh_width_x3,				// grid spacing in x1 direction (uniform)
	double lower_bound_derivatives    		// lower bound for the first partial derivatives
							// to consider it a limiting case of vanishing
							// partial derivatives
	)
{
	  /* function definitions */
	double cfl_number;
	double scaled_level_set;
	double scaled_volume_donating_region;
	
	int level_set_2_vof( 
	      double level_set, 			// compute the volume of fluid field value from 
	      double d_level_set_d_x1, 			// a given level-set field value
	      double d_level_set_d_x2, 		
	      double d_level_set_d_x3, 		
	      double &volume_of_fluid,		
	      double lower_bound_derivatives    
      );

	int i_index, j_index, k_index;  			// local variables for loop indexing
	double one_over_dx3	=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x3);
	
	
	/* compute the volume of fluid fluxes in the x3 direction */
	  
	for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
	    for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
		for(k_index=0;k_index<number_primary_cells_k+1;k_index++){
		    
		/* check the local cfl number */
		    cfl_number= u_3_velocity_new[i_index][j_index][k_index]*
		    actual_time_step_level_set*one_over_dx3;
		    if(cfl_number<-1.0 || cfl_number>1.0) 
		    {
			    /* the magnitude of the cfl number exceeds unity   */
			    /* the algorithm is terminated. the time step size */
			    /* has to be reduced			       */
			    std::cerr<< "cfl number exceeds unity";
			    std::cerr<< "in compute_vof_flux_x3 \n";
			    std::cerr<< " line 68 \n";
			    exit(1);
		    }  
		    if( cfl_number >=0.0)
		    {
			/* use is made of the original function that */
			/* computes the volume of fluid from a given */
			/* level-set field to compute the volume of  */
			/* the donating region of the volume of fluid*/
			/* flux. The conversion is sign dependent    */
			/* consider the case cfl>0                   */
			scaled_level_set=
			level_set[i_index][j_index][k_index]+
			0.5*(1.0-cfl_number)*
			  d_level_set_d_x3[i_index][j_index][k_index];
			if(!level_set_2_vof( scaled_level_set,
				       d_level_set_d_x1[i_index][j_index][k_index],
					d_level_set_d_x2[i_index][j_index][k_index],
					  d_level_set_d_x3[i_index][j_index][k_index]*cfl_number,
					    scaled_volume_donating_region,
					      lower_bound_derivatives));
			flux_x3[i_index][j_index][k_index]=
			  cfl_number*scaled_volume_donating_region;
		    }
		    else
		    {
			    
			/* consider the case cfl<0                   */
			scaled_level_set=
			level_set[i_index][j_index][k_index+1]-
			0.5*(1.0+cfl_number)*
			  d_level_set_d_x3[i_index][j_index][k_index+1];
			if(!level_set_2_vof( scaled_level_set,
				       d_level_set_d_x1[i_index][j_index][k_index+1],
					d_level_set_d_x2[i_index][j_index][k_index+1],
					  -d_level_set_d_x3[i_index][j_index][k_index+1]*cfl_number,
					    scaled_volume_donating_region,
					    lower_bound_derivatives));
			flux_x3[i_index][j_index][k_index]=
			  cfl_number*scaled_volume_donating_region;
		    
		    }
		}
	    }
	}
}
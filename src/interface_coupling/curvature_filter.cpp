#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the filter function for the curvature smoothing.        */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* For the smoothing a diffusion equation is solved. The diffusion coefficient  */
/* in the equation is the filter value. This coefficient drops off to zero      */
/* outside a small band around the interface if the filtering is active.        */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/
    double curvature_filter(
	      double level_set_value,			// level-set value 
	      double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	      double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	      double mesh_width_x3			// grid spacing in x3 direction (uniform)
		   )
    
    {
      /* compute the scaling for the filter function 	      */
      /* it is equal to 2 / sqrt(dx1^2 + dx2^2 + dz3^2)	      */
      double one_over_scaling_distance;
      double curvature_filter_value;
      
      one_over_scaling_distance=1.0/(0.5*sqrt( mesh_width_x1*mesh_width_x1+mesh_width_x2*mesh_width_x2+
					  mesh_width_x3*mesh_width_x3));
      curvature_filter_value=1.0-exp(-1.0*
				    (one_over_scaling_distance*level_set_value)*
				     (one_over_scaling_distance*level_set_value)
				      );
				    
      return curvature_filter_value;
    }

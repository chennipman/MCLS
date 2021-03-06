#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the scaled viscosity:mu/rho_minus                       */
/*                                                                 		*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/*  In this approach the viscosity is not discontinuous but smooth.             */
/*  Use is made of a smoothed heaviside function to interpolate between         */
/*  the viscosity of the two phases.                                            */
/********************************************************************************/
EXPORT double compute_scaled_viscosity(			
	  double level_set,
	  double mesh_width_x1,				// grid spacing in x1 direction (uniform)
	  double mesh_width_x2,				// grid spacing in x2 direction (uniform)
	  double mesh_width_x3,				// grid spacing in x3 direction (uniform)
	  double smoothing_distance_factor,		// the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
	  double rho_minus_over_mu_minus,		// this was the 'Reynolds' number
							// in the original implementation of Sander
	  double mu_plus_over_mu_minus			// ratio of the viscosities of both phases
				     )  				
     
     
      {
      double smoothing_distance;			// the heaviside function is smoothed over
							// an interval of width 2*smoothing_distance
      double smallest_meshwidth;			// the smaller of the three meshwidths
      double heaviside_value;				// value of the smoothed heaviside function
      double scaled_viscosity;				// local viscosity scaled with rho_minus
      
      
      /* determine the smallest meshwidth */
      
      smallest_meshwidth=std::min(mesh_width_x1,std::min(mesh_width_x2, mesh_width_x3));
      
      
      /* the smoothing distance is twice the smoothing distance factor times the smallest meshwidth */
      
      smoothing_distance=smoothing_distance_factor*smallest_meshwidth;
      
      if(level_set<=-1.0*smoothing_distance)
      {
	/*on the negative side, outside smoothing interval */
	  heaviside_value=0.0;
      }
      else
      {
	  if(level_set>=-1.0*smoothing_distance)
	  {
	  /*on the positive side, outside smoothing interval */
	      heaviside_value=1.0;
	  }
	  else
	  {
	    /* in the smoothing region */
	    /* use is made of 2*atan(1)=pi/2 to replace pi/2 */
	      heaviside_value=0.5*(1.0+sin(level_set/smoothing_distance*2*atan(1.0)));
	  }
      }
      
      
      /* scale the viscosity */
      /* we need mu/rho_m */

      /* we use						        */
      /*  mu_minus  		mu_plus				*/
      /*----------- (1 + H *  (--------   - 1))   		*/
      /*  rho_minus		mu_minus			*/

      /* which is exactly:  				        */
      /*     1      		mu_plus				*/
      /*----------- (mu_minus + H *  (mu_plus-mu_minus))        */
      /*  rho_minus					        */
      
        scaled_viscosity=1.0/ rho_minus_over_mu_minus*
				      (1.0+heaviside_value*(mu_plus_over_mu_minus-1.0));
	return scaled_viscosity;
}

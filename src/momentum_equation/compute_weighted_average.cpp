#include<cstdlib>
#include<math.h>
#include<algorithm>
/********************************************************************************/
/*  Function to compute the apply linear interpolation to the density           */
/*  based on the level-set values in the neighbouring cells.			*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* In the original implementation this function was called harmonic averaging   */
/* but it seems to be just linear interpolation based on the level-set values   */
/* faces are exactly half way between the centers of the primary cells.         */
/* If the level-set value is positive:						*/
/*				 the scaled density is = rho_plus/rho_minus     */
/* If the level-set value is negative:						*/
/*				 the scaled density is = 1                      */
/********************************************************************************/
EXPORT double compute_weighted_average(
	    double level_set_left, 		// level-set value in left hand adjacent cell
	    double level_set_right,		// level-set value in right hand adjacent cell
	    double density_left,		// density in left hand adjacent cell
	    double density_right		// density in right hand adjacent cell
	    )	
  
  {
    double weight_factor_left;			// weight factor for the left hand cell value
    double weight_factor_right;			// weight factor for the right hand cell value
    double weighted_average;			// weighted average value based on the level-set value
    
    if( fabs(level_set_right-level_set_left) >0)
    {
        /* the density at the interface is obtained using linear interpolation */
	/* based on the level-set field, that is a signed distance function */
	
	/* the weight for the right value is the distance between the interface and the center */
	/* of the left neighbouring cell */
	weight_factor_left = -1.0*level_set_left/(level_set_right-level_set_left); 
	
	/* the following statements were in Sanders code, but seem redundant */
	
	weight_factor_left = std::min(weight_factor_left, 1.0);
	weight_factor_left = std::max(weight_factor_left, 0.0);
    }
    else
    {
	weight_factor_left = 1.0;
      
    }
	/* of course the weights add up to 1 */
	
	weight_factor_right=1.0-weight_factor_left;
	
	
	weighted_average=weight_factor_right*density_right+weight_factor_left*density_left;
	return weighted_average;
  }
    
    

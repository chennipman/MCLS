#include<cstdlib>
#include<math.h>
#include<algorithm>
	
/********************************************************************************/
/********************************************************************************/
/*  Function to compute convection velocity for the mass redistribution         */
/*  algorithm									*/
/*  equation									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/********************************************************************************/
      double compute_redistribution_velocity(
	      double left_hand_value_level_set, 	// left hand value of the level-set field
							// with respect to cell face
	      double right_hand_value_level_set,	// right hand value of the level-set field
							// with respect to cell face
	      double meshwidth				// grid spacing / mesh width
	)
      {
	double sign(double value, double set_sign);    // return the sign of 'value'
	
	double mean_value_level_set;			// average of left and right hand side value
							// of the level-set field
	double redistribution_velocity_threshold=0.1;	// threshold value for the redistribution
							// velocity
	double redistribution_velocity;			// velocity for the mass redistribution
							// algorithm
	
	/* compute the mean value of the level-set value to evaluate the threshold value */
	
	mean_value_level_set=0.5*(fabs(left_hand_value_level_set)+fabs(right_hand_value_level_set));
	
	/* make sure the velocity is directed towards the location of the interface, i.e. the */
	/* level-set=0 contourline							      */
	
	redistribution_velocity= -sign(1.0, right_hand_value_level_set-left_hand_value_level_set)*
			      sign(1, right_hand_value_level_set+left_hand_value_level_set);
	if( fabs(fabs(right_hand_value_level_set)-fabs(left_hand_value_level_set))<
	  redistribution_velocity_threshold*(std::max(mean_value_level_set, meshwidth)))
	{
	    redistribution_velocity=0.0;
	}
	return redistribution_velocity;
	
      }
      

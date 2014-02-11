#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the Dirac delta function.                               */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The smoothed delta function is defined as:					*/
/*										*/
/*										*/
/*	    / 1/(2*alpha) * (1 + cos(phi * pi / alpha)),	|phi| <= alpha	*/
/*  delta = |									*/
/*	    \ 0,						|phi| > alpha	*/
/* where a is the smoothing distance						*/
/********************************************************************************/
      double delta_function(		 
	      double level_set,				// level-set value 
	      double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	      double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	      double mesh_width_x3,			// grid spacing in x3 direction (uniform)
	      double smoothing_distance_factor		// the smoothing distance is smoothing_distance_factor
	  )
 {
      double smoothing_distance;			// the heaviside function is smoothed over
							// an interval of width 2*smoothing_distance
      double minimal_meshwidth;				// smallest mesh width
      double smallest_meshwidth;			// the smaller of the three meshwidths
      double heaviside_left;
      double heaviside_right;
      double delta;
      
      /* determine the smallest meshwidth */
      
      smallest_meshwidth=std::min(mesh_width_x1,std::min(mesh_width_x2, mesh_width_x3));
      
      
      /* the smoothing distance is twice the smoothing distance factor times the smallest meshwidth */
      
      smoothing_distance=smoothing_distance_factor*smallest_meshwidth;
      
      delta = 0.0;
      if (fabs(level_set) < smoothing_distance) 
      {
        delta = (1.0 + cos(atan(1.0)*4.0*level_set/smoothing_distance)) * 0.5/smoothing_distance;
	 
	 delta = 0.5/smoothing_distance*2.0*atan(1.0)*cos(level_set/smoothing_distance*2*atan(1.0));
// 	heaviside_right=0.5*(1.0+sin(level_set_right/smoothing_distance*2*atan(1.0)));
	     }
      return delta;
 }
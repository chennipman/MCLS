#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the derivative of the smoothed heaviside function       */
/*  by numerical approximation. Unclear to me why this is necessary		*/
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The smoothed heaviside function is defined as:				*/
/*										*/
/*	    /  0,  					 phi < -alpha		*/
/* H_a(phi)=|  1/2 * (1 + sin(phi*pi / (2*alpha))),	|phi| <= alpha		*/
/* 	    \  1,					 phi > alpha		*/
/*	             / [H_a(Phi1) - H_a(Phi0)] / [Phi1 - Phi0], if Phi0 /= Phi1	*/
/* Then: H_a'(phi) = |								*/
/*  		     \ Delta_a(Phi0),			         if Phi0 = Phi1 */
/* where a is the smoothing distance						*/
/********************************************************************************/
      double computed_derivative_heaviside_function(		 
	      double level_set_left,			// level-set value in left hand neighbouring cell
	      double level_set_right,	  		// level-set value in right hand neighbouring cell
	      double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	      double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	      double mesh_width_x3,			// grid spacing in x3 direction (uniform)
	      double smoothing_distance_factor		// the smoothing distance is smoothing_distance_factor
	  )
 {
      double delta_function(		 		// Dirac delta function
	      double level_set,				
	      double mesh_width_x1,			
	      double mesh_width_x2,			
	      double mesh_width_x3,			
	      double smoothing_distance_factor		
	  );
      double smoothing_distance;			// the heaviside function is smoothed over
							// an interval of width 2*smoothing_distance
      double largest_meshwidth;			// the larger of the three meshwidths
      double heaviside_left;
      double heaviside_right;
      double computed_derivative_heaviside;
      
      /* determine the smallest meshwidth */
      
      largest_meshwidth=std::max(mesh_width_x1,std::max(mesh_width_x2, mesh_width_x3));
      
      
      /* the smoothing distance is twice the smoothing distance factor times the smallest meshwidth */
      
      smoothing_distance=smoothing_distance_factor*largest_meshwidth;
//       smoothing_distance=1.5*largest_meshwidth;
      
      
      if( fabs( level_set_left -level_set_right)>1E-8)	
//      if( level_set_left != level_set_right)	
      {
	
	if(level_set_left <= -1.0*smoothing_distance)
	{
	/*on the negative side, outside smoothing interval */
	    heaviside_left=0.0;
	}
	else
	{
	    if(level_set_left>= smoothing_distance)
	    {
	  /*on the positive side, outside smoothing interval */
		heaviside_left=1.0;
	    }
	    else   
	    {
	    /* in the smoothing region */
	    /* use is made of 2*atan(1)=pi/2 to replace pi/2 */
		heaviside_left=0.5*(1.0+sin(level_set_left/smoothing_distance*2*atan(1.0)));
	    }
	}
	
	if(level_set_right <= -1.0*smoothing_distance)
	{
	/*on the negative side, outside smoothing interval */
	    heaviside_right=0.0;
	}
	else
	{
	    if(level_set_right>= smoothing_distance)
	    {
	  /*on the positive side, outside smoothing interval */
		heaviside_right=1.0;
	    }
	    else
	    {
	    /* in the smoothing region */
	    /* use is made of 2*atan(1)=pi/2 to replace pi/2 */
		heaviside_right=0.5*(1.0+sin(level_set_right/smoothing_distance*2*atan(1.0)));
	    }
	}
	 
	computed_derivative_heaviside=(heaviside_right-heaviside_left)/(level_set_right-level_set_left);
      }
      else
      {
	computed_derivative_heaviside=delta_function(0.5*(level_set_left+level_set_right),				
						      mesh_width_x1, mesh_width_x2, mesh_width_x3,			
							smoothing_distance_factor	);
      }	      
// 	computed_derivative_heaviside=delta_function(0.5*(level_set_left+level_set_right),				
// 						      mesh_width_x1, mesh_width_x2, mesh_width_x3,			
// 							smoothing_distance_factor	);
      
       return computed_derivative_heaviside;
 }


#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
      int ridders_method(			
	    double level_set_left, 			// left hand initial value for root finding
	    double level_set_right, 			// right hand initial value for root finding
	    double function_g_left, 			// left hand initial function value for root finding
	    double function_g_right,			// right hand initial function value for root findin
	    double volume_of_fluid, 			// volume of fluid value for which we need the 
							// the corresponding level-set value
	    double &level_set,				// the level-set value corresponding to given
							// volume of fluid value
	    double d_level_set_d_x1, 			// first partial derivative in x1 
							// direction of level-set
	    double d_level_set_d_x2, 			// first partial derivative in x2 
							// direction of level-set 
	    double d_level_set_d_x3, 			// first partial derivative in x3 
							// direction of level-set
	    double vof_2_level_set_tolerance,		// tolerance in the conversion from volume
							// of fluid value to level-set value
	    int number_iterations_ridder,		// maximum number of iterations allowed in the
							// nonlinear root finding algorithm
	    int lower_bound_derivatives			// lower bound for the first partial derivatives
							// to consider it a limiting case of vanishing
							// partial derivatives
			)

{
      /* function definitions */
      int level_set_2_vof( 
	      double level_set, 			// compute the volume of fluid field value from 
	      double d_level_set_d_x1, 			// a given level-set field value
						
	      double d_level_set_d_x2, 		
						
	      double d_level_set_d_x3, 		
						
	      double &volume_of_fluid,		
	      double lower_bound_derivatives    
      );
      double sign(					// return value times sign(set_sign)
	    double value, 
	    double set_sign
		 );
      double new_mid_point_beta;			// new estimate of the root location based on
							// midpoint search interval
      double function_g_value_beta;			// function value in the new estimate of the
							// root location based on the midpoint of the
							// search interval
      double function_g_value;				// function value in the new estimate of the
							// root location
      double new_volume_of_fluid;			// volume of fluid at the new estimate of the
							// root location
      double ridders_denominator;			// denominator of Ridders expression for the
							// estimate of the root
      double ridders_enumerator;			// enumerator of Ridders expression for the
							// estimate of the root
      int iteration_index_ridder=1;			// index of the iteration in the root finding
							// algorithm
      
      /* first check if either of the two starting values are sufficiently */
      /* close to the root */
      
      std::cout<<"start with vof "<< volume_of_fluid << " \n";
      
      if( fabs(function_g_left)< vof_2_level_set_tolerance)
      {
	    /* left hand starting value is a root */
	    level_set=level_set_left;
	    return 0;
	 
      }
      else
      {
	  if(fabs(function_g_right)< vof_2_level_set_tolerance)
	  {
   /*	     right hand starting value is a root */
	      level_set=level_set_right;
	      return 0;
	  }
	  else
	  {
	      while( fabs(function_g_value)>vof_2_level_set_tolerance &&
		      iteration_index_ridder<number_iterations_ridder)
	      {
		  new_mid_point_beta=0.5*(level_set_left+level_set_right);
		  if(level_set_2_vof(new_mid_point_beta, 
			    d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3, 
			      new_volume_of_fluid, lower_bound_derivatives))
		  {
		      std::cout<< "error in computation of volume of fluid in ridders_method \n";
		      std::cout<< "new mid point, aborting \n";
		      return 1;
		  }
		  function_g_value_beta=new_volume_of_fluid-volume_of_fluid;
		  ridders_enumerator=0.5*(level_set_right-level_set_left)*
				    sign(1.0,function_g_left)*function_g_value_beta;
		  ridders_denominator=sqrt(function_g_value_beta*function_g_value_beta-
					function_g_right*function_g_left);
	      
		  /* make a new estimate of the root */
	      
		  level_set=new_mid_point_beta+ridders_enumerator/ridders_denominator;
		  if(level_set_2_vof(level_set, 
			    d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3, 
			      new_volume_of_fluid, lower_bound_derivatives))
		  {
		      std::cout<< "error in computation of volume of fluid in ridders_method \n";
		      std::cout<< "new mid point, aborting \n";
		      return 1;
		  }
	      
// 		  std::cerr<<" new_volume_of_fluid "<< new_volume_of_fluid << " \n";
	      /* compute the new residual */
	      
		  function_g_value=new_volume_of_fluid-volume_of_fluid;
		  std::cout<<"g_value "<< function_g_value << " \n";
	      
	      /* depending on the sign of the residual, update left or right point */
	      
		  if(function_g_left*function_g_value<0)
		  {
		      level_set_right=level_set;
		      function_g_right=function_g_value;
		  }
		  else
		  {
		      level_set_left=level_set;
		      function_g_left=function_g_value;
		  }
	      
		  /* if the latest approximation is sufficiently close to the root */
		  /* terminate the iteration  				       */
		  if(fabs(function_g_value)<vof_2_level_set_tolerance)
		  {	
		      return 0;
		  }
		  iteration_index_ridder++;
	      }
	  
	      if(fabs(function_g_value)>vof_2_level_set_tolerance)
	      {
		  /* apparently the algorithm failed to converge in the allowed number of iterations */
		  return 1;
	      }
	    
	  }
      }
}
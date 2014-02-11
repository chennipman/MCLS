
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

      int bisection_method(			
	    double level_set_left, 			// left hand initial value for root finding
	    double level_set_right, 			// right hand initial value for root finding
	    double function_g_left, 			// left hand initial function value for root finding
	    double function_g_right,			// right hand initial function value for root findin
	    double volume_of_fluid, 			// volume of fluid value for which we need the 
							// the corresponding level-set value
	    double &level_set,			// the level-set value corresponding to given
							// volume of fluid value
	    double d_level_set_d_x1, 			// first partial derivative in x1 
							// direction of level-set
	    double d_level_set_d_x2, 			// first partial derivative in x2 
							// direction of level-set 
	    double d_level_set_d_x3, 			// first partial derivative in x3 
							// direction of level-set
	    double vof_2_level_set_tolerance,	// tolerance in the conversion from volume
							// of fluid value to level-set value
	    int number_iterations_bisection,		// maximum number of iterations allowed in the
							// nonlinear root finding algorithm
	    double lower_bound_derivatives		// lower bound for the first partial derivatives
							// to consider it a limiting case of vanishing
							// partial derivatives
			)

{
      /* function definitions */
      int level_set_2_vof( 
	      double level_set, 			// compute the volume of fluid field value from 
	      double d_level_set_d_x1, 		// a given level-set field value
						
	      double d_level_set_d_x2, 		
						
	      double d_level_set_d_x3, 		
						
	      double &volume_of_fluid,		
	      double lower_bound_derivatives    
      );
      double sign(					// return value times sign(set_sign)
	    double value, 
	    double set_sign
		 );
      double new_mid_point;				// new estimate of the root location based on
							// midpoint search interval
      double function_g_value_midpoint=100.0;	// function value in the new estimate of the
							// root location based on the midpoint of the
							// search interval
      double new_volume_of_fluid;			// volume of fluid at the new estimate of the
							// root location
      int iteration_index_bisection=0;		// index of the iteration in the root finding
							// algorithm
      double function_g_left_start;			// starting value for function_g_left
      double function_g_right_start;			// starting value for function_g_right
      double level_set_left_start;			// starting value for level_set_left
      double level_set_right_start;			// starting value for level_set_right
      double iterate[5][200];
  
  	/* save the starting values for the nonlinear root finding algorithm for */
	/* those case the solver does not converge */
  
      function_g_left_start=function_g_left;
      function_g_right_start=function_g_right;
      level_set_left_start=level_set_left;
      level_set_right_start=level_set_right;

      /* first check if either of the two starting values are sufficiently */
      /* close to the root */
      
      
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
	      while( fabs(function_g_value_midpoint)>vof_2_level_set_tolerance &&
		      iteration_index_bisection<number_iterations_bisection)
	      {
		  new_mid_point=0.5*(level_set_left+level_set_right);
		  if(level_set_2_vof(new_mid_point, 
			    d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3, 
			      new_volume_of_fluid, lower_bound_derivatives))
		  {
		      cout<< "error in computation of volume of fluid in ridders_method \n";
		      cout<< "new mid point, aborting \n";
		      return 1;
		  }
		  function_g_value_midpoint=new_volume_of_fluid-volume_of_fluid;
	      
		  /* make a new estimate of the root */
	      
		  if(function_g_left*function_g_value_midpoint<0)
		  {
		      level_set_right=new_mid_point;
		      function_g_right=function_g_value_midpoint;
		  }
		  else
		  {
		      level_set_left=new_mid_point;
		      function_g_left=function_g_value_midpoint;
		  }
	      
		  /* if the latest approximation is sufficiently close to the root */
		  /* terminate the iteration  				       */
		  if(fabs(function_g_value_midpoint)<vof_2_level_set_tolerance)
		  {	
// 		    std::cout<<"new_volume_of_fluid"<< new_volume_of_fluid << " \n";
		      
		      level_set=new_mid_point;
		      return 0;
		  }
                iterate[0][iteration_index_bisection]=function_g_value_midpoint;
                iterate[1][iteration_index_bisection]=level_set_left;
                iterate[2][iteration_index_bisection]=level_set_right;
                iterate[3][iteration_index_bisection]=function_g_left;
                iterate[4][iteration_index_bisection]=function_g_right;
                
		  iteration_index_bisection++;
	      }

	      if(fabs(function_g_value_midpoint)>vof_2_level_set_tolerance)
	      {
		  /* apparently the algorithm failed to converge in the allowed number of iterations */
		  cout<<"nonlinear solver failed to converge with the following \n";
		  cout<<"starting values: \n";
		  cout<<"function_g_left : "<<function_g_left_start<<"\n";
		  cout<<"function_g_right: "<<function_g_right_start<<"\n";
		  cout<<"level_set_left  : "<<level_set_left_start<<"\n";
		  cout<<"level_set_right : "<<level_set_right_start<<"\n";
		  cout<<"d_level_set_d_x1: "<<d_level_set_d_x1<<"\n";
		  cout<<"d_level_set_d_x2: "<<d_level_set_d_x2<<"\n";
		  cout<<"d_level_set_d_x3: "<<d_level_set_d_x3<<"\n";
		  cout<<"and end value   : \n";
		  cout<<"function_g_value_midpoint :"<<function_g_value_midpoint<<"\n";
                cout<<"lower_bound_derivatives :"<<lower_bound_derivatives<<"\n";
               int i;
                double level_set_punt;
                double functie_punt;

                ofstream problem_function_file( "problem_function.dat");
                if(!problem_function_file)
                {
                    /* the contructor returned a 0-pointer :-( */
                    cout << "Cannot open file.\n";
                    cout << "In function bisection_method \n";
                    cout << "Line 168 \n";
                    exit(1);
                }
               for(i=1;i<100;i++)
               {
                    level_set_punt=(level_set_right_start-level_set_left_start)/100*i+level_set_right_start;
                    if(!level_set_2_vof(level_set_punt,
                      d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
                        new_volume_of_fluid, lower_bound_derivatives))
                    {
                    functie_punt=new_volume_of_fluid-volume_of_fluid;
                    }
                    else
                    {
                           cout<<"problem in evaluation of level_set_2_vof \n";
                           cout << "In function bisection_method \n";
                           cout << "Line 174 \n";
                    }
                    problem_function_file<<level_set_punt<<" "<<functie_punt<<"\n";
               }
               cout<<"convergence history \n";
               for(i=0;i<number_iterations_bisection;i++)
               {
                cout<< iterate[0][i] <<"  "<< iterate[1][i] <<"  "<< iterate[2][i] <<"  "<< iterate[3][i] <<"  "<< iterate[4][i] <<" \n";
               }

		  return 1; 
	      }
	    
	  }
      }
}






// 		  int i;
// 		  double level_set_punt;
// 		  double functie_punt;
// 		  for(i=1;i<100;i++)
// 		  {
// 			  level_set_punt=(level_set_right_start-level_set_left_start)/100*i;
// 			  if(!level_set_2_vof(level_set_punt, 
// 			    d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3, 
// 			      new_volume_of_fluid, lower_bound_derivatives))
// 			  {
// 			  functie_punt=new_volume_of_fluid-volume_of_fluid;
// 			  }
// 			  else
// 			  {
// 				  std::cout<<"probleem \n";
// 			  }
// 			  std::cout<<level_set_punt<<" "<<functie_punt<<"\n";
// 		  }
// 
// 




















/*
		  for(iteration_index_bisection=1;iteration_index_bisection<100;iteration_index_bisection++)
		  {
			 std::cout<<"iterate :"<<iteration_index_bisection<<" "<<iterate[iteration_index_bisection]<<"\n";
			 
		  }*/

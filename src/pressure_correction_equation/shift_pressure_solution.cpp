#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/*  Function to shift the the value of the pressure such that the minimum       */
/*   value is zero.               						       */
/*  										       */
/*  Programmer	: Duncan van der Heul       					       */
/*  Date	: 10-03-2013       						       */
/*  Update	:        							       */
/********************************************************************************/
/* Notes									       */
/* Because for a closed container the pressure matrix becomes singular the      */
/* solution is not unique but defined up to a constant. To be able to           */
/* compare solutions at different time levels the pressure is shifted.          */
/********************************************************************************/
EXPORT int shift_pressure_solution(
      int total_number_pressure_points,		// total number of points with pressure
      Array1<double> compressed_pressure		 	// pressure rhside
      )
  {
    int one_dimensional_index;			// index of the point in the 1-d array
    double minimum_pressure;			       // minimum value of pressure vector
    
    /* compute the minimum of all components of the pressure vector */
    
    minimum_pressure=1e6;
    for( one_dimensional_index=0;one_dimensional_index<total_number_pressure_points;one_dimensional_index++)
    {
      
	minimum_pressure=std::min(minimum_pressure, compressed_pressure[one_dimensional_index]);
    }
  
    /* shift the pressure such that the minimum will be zero */

    for(one_dimensional_index=0;one_dimensional_index<total_number_pressure_points;one_dimensional_index++)
    {
	compressed_pressure[one_dimensional_index]-=minimum_pressure;
    }
    return 0; 
  }
  

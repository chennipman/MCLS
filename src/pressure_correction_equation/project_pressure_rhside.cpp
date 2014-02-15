#include "../headers/array.h"
/********************************************************************************/
/*  Function to project the right hand side vector on the column space of the   */
/*    pressure correction matrix  						*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The component of the right hand side that is in the direction of the kernel  */
/* of the pressure correction matrix is removed from the right hand side.       */
/* This function should be in-lined.                                            */
/********************************************************************************/
  int project_pressure_rhside(
      int total_number_pressure_points,		// total number of points with pressure
      Array1<double> pressure_rhside		 	// pressure rhside
      )
  {
    int one_dimensional_index;			// index of the point in the 1-d array
    double sum_all_components;			// sum of all components of the rhside vector
    double average_value;			// average of all components of the rhside vector
    
    /* compute the sum of all components of the rhside vector */
    
    sum_all_components=0;
    for( one_dimensional_index=0;one_dimensional_index<total_number_pressure_points;one_dimensional_index++)
    {
      
	sum_all_components+=pressure_rhside[one_dimensional_index];
    }
  
    average_value=sum_all_components/total_number_pressure_points;
    for(one_dimensional_index=0;one_dimensional_index<total_number_pressure_points;one_dimensional_index++)
    {
	pressure_rhside[one_dimensional_index]-=average_value;
    }
    return 0; 
  }
  

#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the volume of fluid field from a given level-set field  */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* For the case of all first partial derivatives nonzero, the standard          */
/* formulation is used. Otherwise a number of limiting cases are considered     */
/* The algorithm is discussed in paragraph 5.4 of Sander's thesis   		*/
/********************************************************************************/
      int level_set_2_vof( 
	      double level_set, 		// level-set field value in this cell
	      double d_level_set_d_x1, 		// first partial derivative in x1 
						// direction of level-set
	      double d_level_set_d_x2, 		// first partial derivative in x2 
						// direction of level-set 
	      double d_level_set_d_x3, 		// first partial derivative in x3 
						// direction of level-set
	      double &volume_of_fluid,		// volume of fluid value in this cell
	      double lower_bound_derivatives    // lower bound for the first partial derivatives
						// to consider it a limiting case of vanishing
						// partial derivatives
	
      )
{
       double complement_volume_of_fluid;	// (volume of fluid)-1 for -1* level set
	    /* function definitions */
     
      int level_set_2_vof_phi_negative( 	// compute volume of fluid from level-set
	      double level_set, 		// field, when level-set field is less than zero.
	      double d_level_set_d_x1, 		
						
	      double d_level_set_d_x2, 		
						
	      double d_level_set_d_x3, 		
						
	      double &volume_of_fluid,		
	      double lower_bound_derivatives    
	      );
 
      if(level_set<=0.0)
      {
	  if(level_set_2_vof_phi_negative(level_set, 
 		  d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
			  volume_of_fluid, lower_bound_derivatives))
	  {
	    std::cerr << "level-set to vof conversion failed in call from level_set_2_vof. \n";
	    std::cerr << "for level_set< 0 \n";
	    exit(1);
	    return 1;
	  }
      }
      else
      {
	  if(!level_set_2_vof_phi_negative(-1.0*level_set, 
 		  d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
			  complement_volume_of_fluid, lower_bound_derivatives))
	  {  
	      volume_of_fluid=1.0-complement_volume_of_fluid;    
	  }
	  else
	  {
	    std::cerr << "level-set to vof conversion failed in call from level_set_2_vof. \n";
	    std::cerr << "for level_set >0 \n";
	    exit(1);
	    return 1;
	  }  
      }
      if(volume_of_fluid >1.01)
      {
	  std::cerr <<"vof too large \n";
	  exit(1);
       	  return 1;
      }
	  
      return 0;
}
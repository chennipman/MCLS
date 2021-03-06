#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the volume of fluid field from a given level-set field  */
/*  when the level-set function is less than zero 				*/
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
EXPORT int level_set_2_vof_phi_negative( 
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
	/*
C-------------------------------------------------------------------------------
C  Psi = f_std::min(Phi, grad Phi) for Phi <= 0.
C  Section 5.4 in thesis.
C-------------------------------------------------------------------------------
*/
      
      double absolute_d_level_set_d_x1;		// absolute value of first partial derivative in x1 
						// direction of level-set
      double absolute_d_level_set_d_x2;		// absolute value of first partial derivative in x2 
						// direction of level-set
      double absolute_d_level_set_d_x3;		// absolute value of first partial derivative in x3 
						// direction of level-set
      
      double level_set_value_vertex_A;		// value of LINEARIZED level-set field at vertex A
      double level_set_value_vertex_B;		// value of LINEARIZED level-set field at vertex B
      double level_set_value_vertex_C;		// value of LINEARIZED level-set field at vertex C
      double level_set_value_vertex_D;		// value of LINEARIZED level-set field at vertex D
      double level_set_value_vertex_E;		// value of LINEARIZED level-set field at vertex E
      
      double d_level_set_d_xi;			// first partial derivative with respect to logical
						// space coordinate xi of level set field 
      double d_level_set_d_eta;			// first partial derivative with respect to logical
						// space coordinate eta of level set field 
      double d_level_set_d_zeta;		// first partial derivative with respect to logical
						// space coordinate zeta of level set field 
      

      if(level_set>0.0)
      {
	std::cerr<< "level_set value positive in level_set_2_vof_phi_negative.. \n";
	std::cerr<< "Check the source code.\n";
	return 1;
      }
      
      absolute_d_level_set_d_x1=fabs(d_level_set_d_x1);
      absolute_d_level_set_d_x2=fabs(d_level_set_d_x2);
      absolute_d_level_set_d_x3=fabs(d_level_set_d_x3);


      /* d_level_set_d_xi is the largest absolute value of the three first derivatives */
      
      d_level_set_d_xi=std::max(std::max(absolute_d_level_set_d_x1, absolute_d_level_set_d_x2), absolute_d_level_set_d_x3);
      
      /* d_level_set_d_zeta is the smallest absolute value of the three first derivatives */
      
      d_level_set_d_zeta=std::min(std::min(absolute_d_level_set_d_x1, absolute_d_level_set_d_x2), absolute_d_level_set_d_x3);
      
      /* d_level_set_d_eta is the remaining, not yet assigned first derivative */
      
      d_level_set_d_eta= absolute_d_level_set_d_x1+absolute_d_level_set_d_x2+absolute_d_level_set_d_x3-
		      d_level_set_d_xi - d_level_set_d_zeta;

      /* compute the value of the linearized level_set field in the vertices of the control volume */
      
      level_set_value_vertex_A = std::max (level_set - 0.5* ( (-d_level_set_d_xi) - d_level_set_d_eta - d_level_set_d_zeta), 0.0);
      level_set_value_vertex_B = std::max (level_set - 0.5* ( (-d_level_set_d_xi) - d_level_set_d_eta + d_level_set_d_zeta), 0.0);
      level_set_value_vertex_C = std::max (level_set - 0.5* ( (-d_level_set_d_xi) + d_level_set_d_eta - d_level_set_d_zeta), 0.0);
      level_set_value_vertex_D = std::max (level_set - 0.5* (  d_level_set_d_xi - d_level_set_d_eta - d_level_set_d_zeta), 0.0);
      level_set_value_vertex_E = std::max (level_set + 0.5* (  d_level_set_d_xi - d_level_set_d_eta - d_level_set_d_zeta), 0.0);
      
	    
      /* limit cases where one or more of the derivatives in the linearisation vanish have to be taken into account */

      if( d_level_set_d_xi > lower_bound_derivatives )
      {
	    if( d_level_set_d_eta > lower_bound_derivatives )
	    {
		  if( d_level_set_d_zeta > lower_bound_derivatives )
		  { 
		    /* this is the standard case, where all derivatives are nonzero */
		    volume_of_fluid= ( level_set_value_vertex_A*level_set_value_vertex_A*level_set_value_vertex_A-
					 level_set_value_vertex_B*level_set_value_vertex_B*level_set_value_vertex_B-
					   level_set_value_vertex_C*level_set_value_vertex_C*level_set_value_vertex_C-
					     level_set_value_vertex_D*level_set_value_vertex_D*level_set_value_vertex_D+
					       level_set_value_vertex_E*level_set_value_vertex_E*level_set_value_vertex_E)/
					       (6.0*d_level_set_d_xi*d_level_set_d_eta*d_level_set_d_zeta);
		  }
		  else
		  {
		  /* limiting case of vanishing d_level_set_d_zeta */
//                      if( level_set_value_vertex_A>0 && level_set_value_vertex_C>0)
//                      {
//                      volume_of_fluid=level_set/d_level_set_d_xi+0.5;
//                      }
//                      else
//                      {
		      volume_of_fluid= ( level_set_value_vertex_A*level_set_value_vertex_A-
					    level_set_value_vertex_C*level_set_value_vertex_C)/
					    (2.0*d_level_set_d_xi*d_level_set_d_eta);
//                      }
		  }  
		  
	    }
	    else
	    {
	         /* limiting case of vanishing d_level_set_d_zeta AND d_level_set_d_eta */
		  
		  volume_of_fluid= level_set_value_vertex_A/d_level_set_d_xi;
		      
	    }
      }
      else
      {

		/* limiting case of ALL partial derivatives vanishing */
		if(level_set <-1E-07)
		{
		  /* set the volume of fluid to zero */
		    volume_of_fluid=0.0;
		}
		else
		{
		  /* set the volume of fluid to 0.5 */
		  
		    volume_of_fluid=0.5;
		}
      }
      return 0;
      }

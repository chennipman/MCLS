#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the upwind flux for the mass redistribution             */
/*  equation									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/********************************************************************************/
      double upwind_flux_mass_redistribution( 
	      double mass_redistribution_velocity,	// (nonphysical) velocity for the 
							// mass redistribution equation at 
							// the cell face
	      double left_hand_value_correction, 	// left hand value of the correction
							// used to define the upwind flux
	      double right_hand_value_correction)	// right hand value of the correction
							// used to define the upwind flux
	{
	 double sign(double value, double set_sign);    // return the sign of 'value'
	 double upwind_flux_mass_redistribution;  // flux for the mass redistribution equation
						  // based on first order upwind discretisation
	 double upwind_selector;		  // switch to take the upwind value depending
						  // on the sign of the interface velocity
	 
	 upwind_selector=0.5+sign(0.5, mass_redistribution_velocity);
	 upwind_flux_mass_redistribution=mass_redistribution_velocity*
		  (left_hand_value_correction+ 
		      upwind_selector*(left_hand_value_correction-right_hand_value_correction));
	
	  return upwind_flux_mass_redistribution;
      }


      
/********************************************************************************/
/*  Function to compute the scaled density at the center of a cell face.        */
/*  boundary conditions								*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The centers of the cell faces of the primary cells coincide with the       o */
/* location of the velocity unknowns. On the uniform Cartesian grid, the cell   */
/* faces are exactly half way between the centers of the primary cells.         */
/* If the level-set value is positive:						*/
/*				 the scaled density is = rho_plus/rho_minus     */
/* If the level-set value is negative:						*/
/*				 the scaled density is = 1                      */
/********************************************************************************/
       double compute_scaled_density(			
	  double level_set_left, 			// level-set value in left hand neighbouring cell
	  double level_set_right,			// level-set value in right hand neighbouring cell
	  double rho_plus_over_rho_minus		// ratio of density in the two fases
		)
       {
      double sign(double value, double set_sign);    	// return the sign of 'value'=+/- 1.0
      double compute_weighted_average(			// use linear interpolation/weighted averaging
		  double level_set_left, 		// based on the level-set values
		  double level_set_right,
		  double density_left,
		  double density_right
	    );	
      double density_left;
      double density_right;
      double scaled_density;
      
      /* when the level-set value is not equal to zero, it is straightforward to pick the right value */
      /* for the density 									      */
      
      density_left =1.0 + ( sign(0.5, level_set_left ) + 0.5)*(rho_plus_over_rho_minus - 1.0);
      density_right=1.0 + ( sign(0.5, level_set_right) + 0.5)*(rho_plus_over_rho_minus - 1.0);

      /* if one of the two neighbouring points is at the interface, set the density to the average of the */
      /* density of the two phases 									  */
      if(level_set_left ==0.0) density_left =0.5*(rho_plus_over_rho_minus+1.0);
      if(level_set_right==0.0) density_right=0.5*(rho_plus_over_rho_minus+1.0);
      
      scaled_density=compute_weighted_average(level_set_left, level_set_right, density_left, density_right);
      return scaled_density;
       }
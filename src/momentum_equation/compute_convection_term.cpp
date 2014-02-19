      
/********************************************************************************/
/*  Function to compute the convection term for the momentum equation           */
/*  conditions									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function computes the convection term for the momentum equation.        */
/* In the calling function the term is multiplied by the local density.         */
/* In the original function of Sander there was the option to use an upwind     */
/* scheme for this term, but this was rarely, if ever used.                     */
/* For the moment it is left out, but could be included later. Also classic     */
/* higher order kappa-schemes could be used here.                               */
/* This is the reason why it is contained within a separate function		*/
/********************************************************************************/
EXPORT double compute_convection_term( 			
	      double u_1_velocity_cell_center,		// u_1 velocity component at center of
							// control volume of u_alpha
	      double u_2_velocity_cell_center,		// u_2 velocity component at center of
							// control volume of u_alpha		
	      double u_3_velocity_cell_center,		// u_3 velocity component at center of
							// control volume of u_alpha
	      double u_alpha_center,			// u_alpha velocity component at center of
							// control volume of u_alpha 
	      double u_alpha_i_min,			// u_alpha velocity component at i_index+1
							// control volume of u_alpha 
	      double u_alpha_i_plus,			// u_alpha velocity component at i_index-1
							// control volume of u_alpha 
	      double u_alpha_j_min,			// u_alpha velocity component at j_index-1
							// control volume of u_alpha 
	      double u_alpha_j_plus,			// u_alpha velocity component at j_index+1
							// control volume of u_alpha 
	      double u_alpha_k_min,			// u_alpha velocity component at k_index-1
							// control volume of u_alpha 
	      double u_alpha_k_plus,			// u_alpha velocity component at k_index+1
							// control volume of u_alpha 
	      double mesh_width_x1,			// grid spacing in x1 direction (uniform)
	      double mesh_width_x2,			// grid spacing in x2 direction (uniform)
	      double mesh_width_x3			// grid spacing in x3 direction (uniform)
	      )
     {
	      double convection_term;
	      double one_over_dx1=    			// 1/(grid spacing in x1 direction)
		  1.0/(mesh_width_x1);
	      double one_over_dx2=    			// 1/(grid spacing in x2 direction)
		  1.0/(mesh_width_x2);
	      double one_over_dx3=    			// 1/(grid spacing in x3 direction)
		  1.0/(mesh_width_x3);	
	      
	      /* standard centered second-order approximation of the first derivative */
	      
	      convection_term=0.5*(
		    u_1_velocity_cell_center*one_over_dx1*(u_alpha_i_plus-u_alpha_i_min)+
		    u_2_velocity_cell_center*one_over_dx2*(u_alpha_j_plus-u_alpha_j_min)+
		    u_3_velocity_cell_center*one_over_dx3*(u_alpha_k_plus-u_alpha_k_min)
				  );
	      return convection_term;
     }

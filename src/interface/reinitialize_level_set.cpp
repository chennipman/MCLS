#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to reinitialize the level-set field                                */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The level-set field is reinitialized by solving a diffusion type equation    */
/* in the vicinity of the interface.                                            */
/* the reinitialization equation is given in Sander's thesis as equations 6.4   */
/* which should be studied with equation 5.6                                    */
/* The standard Sussman term                                                    */
/* N_h(Phi^k,Phi^0) is weighted with (1-q(Phi^0)) and a term is added           */
/*     Phi^0- Phi^k   								*/
/*    --------------  * q(Phi^0)						*/
/*        Delta t'								*/
/********************************************************************************/
EXPORT void reinitialize_level_set(
      Array3<double> level_set_star, 		// level set field at star time level
      int number_primary_cells_i,		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,		// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,			// grid spacing in x1 direction (uniform)
      double mesh_width_x2,			// grid spacing in x2 direction (uniform)
      double mesh_width_x3,			// grid spacing in x3 direction (uniform)
      double cfl_number_reinitialization,	// courant-friedrichs-lewy number for the reinitialization
						// equation
      int maximum_reinitialization_steps,	// maximum number of time steps in the reinitialization
						// algorithm
      double tolerance_reinitialization	        // stop the reinitialization when the infinite norm of
						// the pseudo time derivative has fallen below this
						// tolerance value
	 )
  

 {
      Array3<double> level_set_0;					// level-set field, after advection
								// but not reinitialized
      Array3<double> level_set_residual_reinit;                      // residual of the reinitialization
                                                                // equation
      Array3<double> d_level_set_d_x1;				// first partial derivative of
								// the level-set field wrt x1
								// second order central approximation
								// defined at the cell faces
      Array3<double> d_level_set_d_x2;				// first partial derivative of 
								// the level-set field wrt x2
								// second order central approximation
								// defined at the cell faces
      Array3<double> d_level_set_d_x3;				// first partial derivative of
 								// the level-set field wrt x3
								// second order central approximation
								// defined at the cell faces
      Array3<double> d_level_set_d_t_prime;			// time derivative in the level-set
 								// reinitialization equation          
								// this defines the complete right hand side
								// of the equation, before used to update the level-set
								// defined at the cell faces
      Array3<double> weighting_function_q0;			// weighting function in the level-set
								// reinitialization equation 
      double scalefactor_bandwidth_reinitialization=1.0;	// scale factor for the width of the reinitialization zone			
      double time_step_reinitialization;			// time step for the reinitialization equation
      double one_over_dx1	=    				// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      double one_over_dx2	=    				// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
      double one_over_dx3	=    				// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);
      int i_index, j_index, k_index;  			        // local variables for loop indexing
      double alpha_q_0;					        // alpha_0
								// term in the weighting function q_0
								// of the reinitialization equation
      double one_over_alpha_q_0;				// reciprocal value of the alpha_0
								// term in the weighting function q_0
								// of the reinitialization equation
      double one_over_time_step_reinitialization;		// reciprocal value of the time
								// step for the reinitialization
      double weight_factor_q_0;				        // weight factor in the weighting
								// function of the reinitialization equation
      double exponent_q_0=2;					// exponent of the scaled level-set field
								// in the weighting function q_0
								// of the reinitialization equation
      double d_x1_left__plus;					// see Eq. 5.7 Sander's thesis
      double d_x1_left__min_;					// see Eq. 5.7 Sander's thesis
      double d_x1_right_plus;					// see Eq. 5.7 Sander's thesis
      double d_x1_right_min_;					// see Eq. 5.7 Sander's thesis
      double d_x2_left__plus;					// see Eq. 5.7 Sander's thesis
      double d_x2_left__min_;					// see Eq. 5.7 Sander's thesis
      double d_x2_right_plus;					// see Eq. 5.7 Sander's thesis
      double d_x2_right_min_;					// see Eq. 5.7 Sander's thesis
      double d_x3_left__plus;					// see Eq. 5.7 Sander's thesis
      double d_x3_left__min_;					// see Eq. 5.7 Sander's thesis
      double d_x3_right_plus;					// see Eq. 5.7 Sander's thesis
      double d_x3_right_min_;					// see Eq. 5.7 Sander's thesis
      double time_derivative_phi_0_part;			// see Eq. 5.7 Sander's thesis
      double inf_norm_residual_reinitialization;		// infinite norm of the residual
								// of the vector of pseudo time
								// derivatives of the reinitialization equation
      int index_reinitialization_step;			        // index of the reinitialization step
      
      std::cerr<< " begin reinitialize_level_set \n";
     
	/* allocate memory to save the un-reinitialized level-set field */
	
	level_set_0.create( number_primary_cells_i+2, 
			number_primary_cells_j+2  , number_primary_cells_k+2);
	
	
	/* allocate memory for the diffusive fluxes */
	
	d_level_set_d_x1.create( number_primary_cells_i+1, 
			number_primary_cells_j+2, number_primary_cells_k+2);
	d_level_set_d_x2.create( number_primary_cells_i+2, 
			number_primary_cells_j+1, number_primary_cells_k+2);
	d_level_set_d_x3.create( number_primary_cells_i+2, 
			number_primary_cells_j+2, number_primary_cells_k+1);

	/* allocate memory for the time derivative part of the equation */

	d_level_set_d_t_prime.create(number_primary_cells_i+2, number_primary_cells_j+2, 
					     number_primary_cells_k+2);

	/* allocate memory for the weighting function in the re-initialization equation*/
        /* and for the residual of the reinitialization equation                       */

	weighting_function_q0.create( number_primary_cells_i+2  , 
			number_primary_cells_j+2  , number_primary_cells_k+2);
	level_set_residual_reinit.create( number_primary_cells_i+2  , 
                        number_primary_cells_j+2  , number_primary_cells_k+2);


	time_step_reinitialization=cfl_number_reinitialization/
				(one_over_dx1+one_over_dx2+one_over_dx3);
      
      /* handle the special case of two-dimensional flow for the determination of the time step size */
	

      if( number_primary_cells_i==1) time_step_reinitialization=cfl_number_reinitialization/
	    (one_over_dx2+one_over_dx3);
      if( number_primary_cells_j==1) time_step_reinitialization=cfl_number_reinitialization/
	    (one_over_dx1+one_over_dx3);
      if( number_primary_cells_k==1) time_step_reinitialization=cfl_number_reinitialization/
	    (one_over_dx1+one_over_dx2);
     

      one_over_time_step_reinitialization= 1.0/time_step_reinitialization;

      /* save the non mass-conserving level-set field */
	    
      copy_cell_centered_field( level_set_star, level_set_0,
			  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);

      

      /* compute the 'band width' for the reinitialization  and the parameters of the function    */
      /* there is not really a band width, because the function approaches 0 in an asymptotic way */
      /* q(Phi_0)= exp( weight_factor_q_0*(Phi_0/alpha_q_0)^exponent_q_0)                         */
	  
	  
 

      /* the two parameters of q_0 */
      
      weight_factor_q_0=-1.0*(exponent_q_0-1.0)/exponent_q_0;
      
      /* to handle the two-dimensional case, make sure that only the relevant meshwidths are 	*/
      /* taken into account in the computation of the bandwidth 			     		*/
      /* this part is documented in Sanders thesis on page 56                                	*/
      /* scalefactor_bandwidth_reinitialization has been added later 			     	*/
      /* take into account that the solution of the level-set equation is not based on       	*/
      /* first principles......								     	*/
      
      alpha_q_0=0.0;
      
      if(number_primary_cells_i>1) alpha_q_0+= mesh_width_x1*mesh_width_x1;
      if(number_primary_cells_j>1) alpha_q_0+= mesh_width_x2*mesh_width_x2;
      if(number_primary_cells_k>1) alpha_q_0+= mesh_width_x3*mesh_width_x3;

      one_over_alpha_q_0=sqrt(3.0/alpha_q_0)/scalefactor_bandwidth_reinitialization;
      
      for(i_index=1;i_index<number_primary_cells_i+1; i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1; j_index++)
	  {
		for(k_index=1;k_index<number_primary_cells_k+1; k_index++)
		{
		  weighting_function_q0[i_index][j_index][k_index]=
			    exp(weight_factor_q_0*
			      pow(one_over_alpha_q_0*level_set_star[i_index][j_index][k_index], exponent_q_0));
		}
	  }
      }
      

       
      for( index_reinitialization_step=1; 
	      index_reinitialization_step<=maximum_reinitialization_steps;
		    index_reinitialization_step++)
      {
	
      /* compute the first derivatives of the level-set field at the cell faces */
      /* that are used to define the fluxes for the reinitialization equation */
	
       /* in each iteration the maximum norm of the time-derivative of the reinitialization equation 	*/
       /* is determined and used as a convergence criterium 					    	*/
       /* initially set to zero									    	*/

	    inf_norm_residual_reinitialization=0.0;
            
            
      /*    extrapolate the level set field to the virtual cells */
	    
            field_extrapolate_boundary( level_set_star,                        
                                        number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
	    
      /*   compute the normal derivative of the level set field at the cell faces */

	    compute_normal_derivative_at_faces(level_set_star, 
						d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,			
						  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,			
						    mesh_width_x1, mesh_width_x2, mesh_width_x3);

      /* details of all terms are explained in equations 5.7-5.10 on page 35 of the thesis */
      /* of Sander									  */
      
	    for(i_index=1;i_index<number_primary_cells_i+1; i_index++)
	    {
		for(j_index=1;j_index<number_primary_cells_j+1; j_index++)
		{
		      for(k_index=1;k_index<number_primary_cells_k+1; k_index++)
		      {
			
			  d_x1_left__plus=std::max(d_level_set_d_x1[i_index-1][j_index  ][k_index  ], 0.0);
			  d_x1_left__min_=std::min(d_level_set_d_x1[i_index-1][j_index  ][k_index  ], 0.0);
			  d_x1_right_plus=std::max(d_level_set_d_x1[i_index  ][j_index  ][k_index  ], 0.0);
			  d_x1_right_min_=std::min(d_level_set_d_x1[i_index  ][j_index  ][k_index  ], 0.0);
		      
			  d_x2_left__plus=std::max(d_level_set_d_x2[i_index  ][j_index-1][k_index  ], 0.0);
			  d_x2_left__min_=std::min(d_level_set_d_x2[i_index  ][j_index-1][k_index  ], 0.0);
			  d_x2_right_plus=std::max(d_level_set_d_x2[i_index  ][j_index  ][k_index  ], 0.0);
			  d_x2_right_min_=std::min(d_level_set_d_x2[i_index  ][j_index  ][k_index  ], 0.0);
		      
			  d_x3_left__plus=std::max(d_level_set_d_x3[i_index  ][j_index  ][k_index-1], 0.0);
			  d_x3_left__min_=std::min(d_level_set_d_x3[i_index  ][j_index  ][k_index-1], 0.0);
			  d_x3_right_plus=std::max(d_level_set_d_x3[i_index  ][j_index  ][k_index  ], 0.0);
			  d_x3_right_min_=std::min(d_level_set_d_x3[i_index  ][j_index  ][k_index  ], 0.0);
		      
	    
			  if(level_set_0[i_index][j_index][k_index] > 0.0) 
			  {
			      d_level_set_d_t_prime[i_index][j_index][k_index]=1.0 -
			      sqrt(
				std::max(d_x1_left__plus*d_x1_left__plus, d_x1_right_min_*d_x1_right_min_)+
				std::max(d_x2_left__plus*d_x2_left__plus, d_x2_right_min_*d_x2_right_min_)+
				std::max(d_x3_left__plus*d_x3_left__plus, d_x3_right_min_*d_x3_right_min_)
				  );
	      
			  }
			  else
			  {
			      if(level_set_0[i_index][j_index][k_index] < 0.0) 
			      {
				  d_level_set_d_t_prime[i_index][j_index][k_index]=-1.0 +
				    sqrt(
				      std::max(d_x1_right_plus*d_x1_right_plus, d_x1_left__min_*d_x1_left__min_)+
				      std::max(d_x2_right_plus*d_x2_right_plus, d_x2_left__min_*d_x2_left__min_)+
				      std::max(d_x3_right_plus*d_x3_right_plus, d_x3_left__min_*d_x3_left__min_)
					);
	      
			      }
			      else
			      {
				  d_level_set_d_t_prime[i_index][j_index][k_index]=0.0;
			      }  
  
			  }  

	                   /* time_derivative_phi_0_part is the following term */
	                   /*    Phi_0-Phi_k     */
	                   /*    ------------    */
	                   /*     delta t'       */
	      
			  time_derivative_phi_0_part=one_over_time_step_reinitialization*
				  (level_set_0[i_index][j_index][k_index]-
					level_set_star[i_index][j_index][k_index]);
				  
                          /* now d_level_set_d_t_prime contains the following terms: 		*/
	                  /*                                                   Phi_0-Phi_k     	*/
	                  /*      N_h(Phi^k,Phi^0)(1-q(Phi_0)) +      q(Phi_0) ------------    	*/
	                  /*                                                    delta t'       	*/
	                         			  
	      
			  d_level_set_d_t_prime[i_index][j_index][k_index]+=
			      weighting_function_q0[i_index][j_index][k_index]*
				(time_derivative_phi_0_part- d_level_set_d_t_prime[i_index][j_index][k_index]);
	      
	      
			  /* check if the level-set field is still changing */

			  inf_norm_residual_reinitialization = 
			     fabs(d_level_set_d_t_prime[i_index][j_index][k_index]);

			  std::max(inf_norm_residual_reinitialization, 
						       fabs(d_level_set_d_t_prime[i_index][j_index][k_index]));

                          /* store the residual */
                          
                          level_set_residual_reinit[i_index][j_index][k_index]=
                             fabs(d_level_set_d_t_prime[i_index][j_index][k_index]);

                          
                          
			  /* update the level-set field */
	      
			  level_set_star[i_index][j_index][k_index]+=
			      time_step_reinitialization*d_level_set_d_t_prime[i_index][j_index][k_index];
 	      
		      }
		 }
	    }

// 	    if(inf_norm_residual_reinitialization < tolerance_reinitialization){
// 		  
// 		  /* the reinitialization is converged */
// 		  
// 		  std::cerr<<"Reinintialization converged... \n";
// 		   
// 		  break;
// 	    }
	    
       }
                  std::cerr<<"The inf norm of the residual is :";
                  std::cerr<<inf_norm_residual_reinitialization<<"\n";
                  std::cerr<<"The tolerance is :";
                  std::cerr<<tolerance_reinitialization<<"\n";
                  std::cerr<<"steps performed  :";
                  std::cerr<< index_reinitialization_step <<"\n";

       
// 	if(inf_norm_residual_reinitialization > tolerance_reinitialization &&
// 			index_reinitialization_step>=maximum_reinitialization_steps)
// 	{
// 		  /* the reinitialization failed to converge */
// 		  
// 		  std::cerr<<"Reinintialization failed to converge. \n";
// 		  std::cerr<<"The inf norm of the residual is :";
// 		  std::cerr<<inf_norm_residual_reinitialization<<"\n";
// 		  std::cerr<<"In reinitialize_level_set line 355. \n";
//                   
//                   /* dump the residuals to file  */
//                   
//                   dump_reinitialization_for_debugging( level_set_star,                         
//                                                         level_set_0, level_set_residual_reinit,    
//                                                            number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,                 
//                                                              mesh_width_x1, mesh_width_x2, mesh_width_x3);
//                   exit(1);
// 	}
       
       /* deallocate all temporary storage */
     
	/* allocate memory to save the un-reinitialized level-set field */
	
	level_set_0.destroy();
	
	
	/* allocate memory for the diffusive fluxes */
	
	d_level_set_d_x1.destroy();
	d_level_set_d_x2.destroy();
	d_level_set_d_x3.destroy();

	/* allocate memory for the time derivative part of the equation */

	d_level_set_d_t_prime.destroy();

	/* allocate memory for the weighting function in the re-initialization equation*/
        /* and the residual of the reinitialization equation                           */

        weighting_function_q0.destroy();
        level_set_residual_reinit.destroy();
	
	std::cerr<< " end reinitialize_level_set \n";

}

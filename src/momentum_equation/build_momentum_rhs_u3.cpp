#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
      
/********************************************************************************/
/*  Function to build momentum right hand side for the u2 velocity, without     */
/*  boundary conditions								*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* In the current implementation the generation of the matrix, right-hand-side  */
/* and the application of the boundary conditions are completely separated.     */
/* In the current function the matrix is build WITHOUT considering the          */
/* boundary conditions. Next the boundary conditions are considered in a        */
/* separate 'matrix folding' function that eliminates all known values  and     */
/* discards any reference to virtual values.                                    */
/********************************************************************************/

    void build_momentum_rhs_u3(
      Array1<double> momentum_rhside_u3,			// momentum matrix velocity x1 direction
      Array3<double> level_set, 				// level-set field
      Array3<double> scaled_density_u3,                 // scaled density for the controlvolumes
                                                   // of the momentum equation in x3 direction
      Array3<double> momentum_source_term_u_3, 		// complete source term for the momentum equation
							// in x1 direction=(-p,1+ g_1 +F1)
      Array3<double> u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction

      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double actual_time_step_navier_stokes,		// time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      double rho_plus_over_rho_minus,			// ratio of the densities of the two phases
      double smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
      double rho_minus_over_mu_minus,			// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus			// ratio of the viscosities of both phases
				      				
       
       )
  {
      double compute_scaled_density(			// compute the density at a cell face
	  double level_set_left, 			// based on the level-set values at the cell
	  double level_set_right,			// centers of the two cells that share that face
	  double rho_plus_over_rho_minus);
	  double compute_convection_term( 		// compute the convection term for the right hand
	      double u_1_velocity_cell_center,		// side of the momentum equation for velocity
	      double u_2_velocity_cell_center,		// component u_alpha
	      double u_3_velocity_cell_center,
	      double u_alpha_center,
	      double u_alpha_i_min,
	      double u_alpha_i_plus,
	      double u_alpha_j_min,
	      double u_alpha_j_plus,
	      double u_alpha_k_min,
	      double u_alpha_k_plus,      
	      double mesh_width_x1,				
	      double mesh_width_x2,				
	      double mesh_width_x3
	      );
      double compute_scaled_viscosity(			// compute the local value of the viscosity
	  double level_set,                     	// from the level-set value
	  double mesh_width_x1,				
	  double mesh_width_x2,				
	  double mesh_width_x3,				
	  double smoothing_distance_factor,		
							
	  double rho_minus_over_mu_minus,		
							
	  double mu_plus_over_mu_minus			
		);  				
      int map_index_u3(                    		// map 3-D array index to 1-D array
	  int i_index,					// index
	  int j_index, 				
	  int k_index,  			
	  int number_primary_cells_i,		
	  int number_primary_cells_j,		
	  int number_primary_cells_k		
      );
		  
      int i_index, j_index, k_index;  			// local variables for loop indexing
      double one_over_dx1	=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      double one_over_dx2=    				// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
      double one_over_dx3=    				// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);	
      double level_set_face_x1_min;                     // level-set at the x1- face of the cell
      double level_set_face_x1_pls;                     // level-set at the x1+ face of the cell
      double level_set_face_x2_min;                     // level-set at the x2- face of the cell
      double level_set_face_x2_pls;                     // level-set at the x2+ face of the cell
      double level_set_face_x3_min;                     // level-set at the x3- face of the cell
      double level_set_face_x3_pls;                     // level-set at the x3+ face of the cell
      double level_set_cell_center;			// level-set at cell center
      double viscosity_x1_min;				// viscosity at the x1- face of the cell
      double viscosity_x1_pls;				// viscosity at the x1+ face of the cell
      double viscosity_x2_min;				// viscosity at the x2- face of the cell
      double viscosity_x2_pls;				// viscosity at the x2+ face of the cell
      double viscosity_x3_min;				// viscosity at the x3- face of the cell
      double viscosity_x3_pls;				// viscosity at the x3+ face of the cell
      double density_cell_center;			// scaled density at the center of the cell
      double convection_term;				// convection term in the momentum equation
      double diffusion_term;				// diffusion term in the momentum equation
      double u_1_velocity_cell_center;			// u_1 velocity at cell center 
      double u_2_velocity_cell_center;			// u_2 velocity at cell center 
      double u_3_velocity_cell_center;			// u_3 velocity at cell center 
      int one_dimensional_index;			// index of point in 1-D array
		  
	    

      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
	      {
		  
		  /* compute the level-set value at the center of the control volume */
		  /* and at the centers of the 6 faces */
		  
		  level_set_cell_center=0.5*(
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index  ][j_index  ][k_index+1]
					         );
		  level_set_face_x1_min=0.25*(
					level_set[i_index-1][j_index  ][k_index  ]+
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index-1][j_index  ][k_index+1]+
					level_set[i_index  ][j_index  ][k_index+1]
					      );
		  level_set_face_x1_pls=0.25*(
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index  ][j_index  ][k_index+1]+
					level_set[i_index+1][j_index  ][k_index  ]+
					level_set[i_index+1][j_index  ][k_index+1]
					      );
		  level_set_face_x2_min=0.25*(
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index  ][j_index-1][k_index  ]+
					level_set[i_index  ][j_index  ][k_index+1]+
					level_set[i_index  ][j_index-1][k_index+1]
					      );
		  level_set_face_x2_pls=0.25*(
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index  ][j_index+1][k_index  ]+
					level_set[i_index  ][j_index  ][k_index+1]+
					level_set[i_index  ][j_index+1][k_index+1]
					      );
		  level_set_face_x3_min=level_set[i_index  ][j_index  ][k_index  ];
		  level_set_face_x3_pls=level_set[i_index  ][j_index  ][k_index+1];
		  
		  /* evaluate the density at the center of the control volume */

                density_cell_center=scaled_density_u3[i_index][j_index][k_index];
		  
// 		  density_cell_center=compute_scaled_density(level_set_cell_center, level_set_cell_center,
// 					 rho_plus_over_rho_minus);
		  
		  /* evaluate the viscosity at the six faces of the cell */
		  
		  viscosity_x1_min=compute_scaled_viscosity(level_set_face_x1_min,
			mesh_width_x1, mesh_width_x2, mesh_width_x3, smoothing_distance_factor,
			  rho_minus_over_mu_minus, mu_plus_over_mu_minus);
		  viscosity_x1_pls=compute_scaled_viscosity(level_set_face_x1_pls,
			mesh_width_x1, mesh_width_x2, mesh_width_x3, smoothing_distance_factor,
			  rho_minus_over_mu_minus, mu_plus_over_mu_minus);
		  viscosity_x2_min=compute_scaled_viscosity(level_set_face_x2_min,
			mesh_width_x1, mesh_width_x2, mesh_width_x3, smoothing_distance_factor,
			  rho_minus_over_mu_minus, mu_plus_over_mu_minus);
		  viscosity_x2_pls=compute_scaled_viscosity(level_set_face_x2_pls,
			mesh_width_x1, mesh_width_x2, mesh_width_x3, smoothing_distance_factor,
			  rho_minus_over_mu_minus, mu_plus_over_mu_minus);
		  viscosity_x3_min=compute_scaled_viscosity(level_set_face_x3_min,
			mesh_width_x1, mesh_width_x2, mesh_width_x3, smoothing_distance_factor,
			  rho_minus_over_mu_minus, mu_plus_over_mu_minus);
		  viscosity_x3_pls=compute_scaled_viscosity(level_set_face_x3_pls,
			mesh_width_x1, mesh_width_x2, mesh_width_x3, smoothing_distance_factor,
			  rho_minus_over_mu_minus, mu_plus_over_mu_minus);

		  /* evaluate the velocity at the cell center for the convective term in */
		  /* the momentum equation, remember: nonconservative approach */
		  
		  u_1_velocity_cell_center=0.25*(
					   u_1_velocity_old[i_index  ][j_index  ][k_index  ]+
					   u_1_velocity_old[i_index  ][j_index  ][k_index+1]+
					   u_1_velocity_old[i_index-1][j_index  ][k_index  ]+
					   u_1_velocity_old[i_index-1][j_index  ][k_index+1]
					     );
		  u_2_velocity_cell_center=0.25*(
					   u_2_velocity_old[i_index  ][j_index  ][k_index  ]+
					   u_2_velocity_old[i_index  ][j_index  ][k_index+1]+
					   u_2_velocity_old[i_index  ][j_index-1][k_index  ]+
					   u_2_velocity_old[i_index  ][j_index-1][k_index+1]
					     );
		  u_3_velocity_cell_center=u_3_velocity_old[i_index][j_index][k_index];
		  
		  if(k_index==0)
		  {
			  /* we are dealing with a half control volume    */
			  /* and replace the normal central approximation */
			  /* with a one-sided approximation               */
			  
			convection_term=density_cell_center*
			      compute_convection_term(u_1_velocity_cell_center,
							u_2_velocity_cell_center,
							  u_3_velocity_cell_center,
							      u_3_velocity_old[i_index  ][j_index  ][k_index  ],
							      u_3_velocity_old[i_index-1][j_index  ][k_index  ],
							      u_3_velocity_old[i_index+1][j_index  ][k_index  ],
							      u_3_velocity_old[i_index  ][j_index-1][k_index  ],
							      u_3_velocity_old[i_index  ][j_index+1][k_index  ],
							      u_3_velocity_old[i_index  ][j_index  ][k_index  ],
							      u_3_velocity_old[i_index  ][j_index  ][k_index+1],
							      mesh_width_x1, mesh_width_x2, mesh_width_x3	
 						    );
		  }
		  else
		  {
		      if(k_index==number_primary_cells_k)
		      {
			  /* we are dealing with a half control volume    */
			  /* and replace the normal central approximation */
			  /* with a one-sided approximation               */
			  
			    convection_term=density_cell_center*
				  compute_convection_term(u_1_velocity_cell_center,
							u_2_velocity_cell_center,
							  u_3_velocity_cell_center,
							      u_3_velocity_old[i_index  ][j_index  ][k_index  ],
							      u_3_velocity_old[i_index-1][j_index  ][k_index  ],
							      u_3_velocity_old[i_index+1][j_index  ][k_index  ],
							      u_3_velocity_old[i_index  ][j_index-1][k_index  ],
							      u_3_velocity_old[i_index  ][j_index+1][k_index  ],
							      u_3_velocity_old[i_index  ][j_index  ][k_index-1],
							      u_3_velocity_old[i_index  ][j_index  ][k_index  ],
							      mesh_width_x1, mesh_width_x2, mesh_width_x3	
 						    );
		      }
		      else
		      {
			  /* regular interior point 		      */
			  /* central approximation of convection term */
		      
			    convection_term=density_cell_center*
				  compute_convection_term(u_1_velocity_cell_center,
							u_2_velocity_cell_center,
							  u_3_velocity_cell_center,
							      u_3_velocity_old[i_index  ][j_index  ][k_index  ],
							      u_3_velocity_old[i_index-1][j_index  ][k_index  ],
							      u_3_velocity_old[i_index+1][j_index  ][k_index  ],
							      u_3_velocity_old[i_index  ][j_index-1][k_index  ],
							      u_3_velocity_old[i_index  ][j_index+1][k_index  ],
							      u_3_velocity_old[i_index  ][j_index  ][k_index-1],
							      u_3_velocity_old[i_index  ][j_index  ][k_index+1],
							      mesh_width_x1, mesh_width_x2, mesh_width_x3	
 						    );
		      }
		  }

		  /* evaluate the explicit part of the diffusion part of the momentum equation */
			      
			      
		  diffusion_term=
				 one_over_dx1*(
				   viscosity_x1_pls*one_over_dx3*(
				      u_1_velocity_old[i_index  ][j_index  ][k_index+1]-
					u_1_velocity_old[i_index  ][j_index  ][k_index  ]) -
					
				   viscosity_x1_min*one_over_dx3*(
				      u_1_velocity_old[i_index-1][j_index  ][k_index+1]-
					u_1_velocity_old[i_index-1][j_index  ][k_index  ]))-
					
				 one_over_dx3*(
				   viscosity_x3_pls*one_over_dx1*(
				      u_1_velocity_old[i_index  ][j_index  ][k_index+1]-
					u_1_velocity_old[i_index-1][j_index  ][k_index+1]) -
					
				   viscosity_x3_min*one_over_dx1*(
				      u_1_velocity_old[i_index  ][j_index  ][k_index  ]-
					u_1_velocity_old[i_index-1][j_index  ][k_index  ]))+		
				 
				 one_over_dx2*(
				   
				   viscosity_x2_pls*one_over_dx3*(
				      u_2_velocity_old[i_index  ][j_index  ][k_index+1]-
					u_2_velocity_old[i_index  ][j_index  ][k_index  ]) -
					
				   viscosity_x2_min*one_over_dx3*(
				      u_2_velocity_old[i_index  ][j_index-1][k_index+1]-
					u_2_velocity_old[i_index  ][j_index-1][k_index  ]))-	
					
				 one_over_dx3*(
				   
				   viscosity_x3_pls*one_over_dx2*(
				      u_2_velocity_old[i_index  ][j_index  ][k_index+1]-
					u_2_velocity_old[i_index  ][j_index-1][k_index+1]) -
					
				   viscosity_x3_min*one_over_dx2*(
				      u_2_velocity_old[i_index  ][j_index  ][k_index  ]-
					u_2_velocity_old[i_index  ][j_index-1][k_index  ]));	  
				 
		  one_dimensional_index=map_index_u3(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
   
		    momentum_rhside_u3[one_dimensional_index]=
			density_cell_center*u_3_velocity_old[i_index][j_index][k_index]+
			  actual_time_step_navier_stokes*(-1.0*convection_term+
				    diffusion_term + 
				      density_cell_center*momentum_source_term_u_3[i_index][j_index][k_index]);
	      }
	  }  
      }      
  }
#include "../headers/array.h"
      
/********************************************************************************/
/*  Function to build momentum matrix for the u3 velocity, without boundary     */
/*  conditions									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* In the current implementation the generation of the matrix, right-hand-side  */
/* and the application of the boundary conditions is completely separated.      */
/* In the current function the matrix is build WITHOUT considering the          */
/* boundary conditions. Next the boundary conditions are considered in a        */
/* separate 'matrix folding' function that eliminates all known values  and     */
/* discards any reference to virtual values.                                    */
/********************************************************************************/

EXPORT void build_momentum_matrix_u3(
      Array2<double> momentum_matrix_u3,		// momentum matrix velocity x1 direction
      Array3<double> level_set, 			// level-set field
      Array3<double> scaled_density_u3,                 // scaled density for the controlvolumes
                                                        // of the momentum equation in x3 direction
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double actual_time_step_navier_stokes,	        // time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      double rho_plus_over_rho_minus,		        // ratio of density in the two fases
      double rho_minus_over_mu_minus,		        // this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      double smoothing_distance_factor		        // the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
       )
  {
      int i_index, j_index, k_index;  			// local variables for loop indexing
      double one_over_dx1_squared=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1*mesh_width_x1);
      double one_over_dx2_squared=    			// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2*mesh_width_x2);
      double one_over_dx3_squared=    			// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3*mesh_width_x3);	
      double level_set_face_x1_min;                     // level-set at the x1- face of the cell
      double level_set_face_x1_pls;                     // level-set at the x1+ face of the cell
      double level_set_face_x2_min;                     // level-set at the x2- face of the cell
      double level_set_face_x2_pls;                     // level-set at the x2+ face of the cell
      double level_set_face_x3_min;                     // level-set at the x3- face of the cell
      double level_set_face_x3_pls;                     // level-set at the x3+ face of the cell
//       double level_set_cell_center;			// level-set at cell center
      double viscosity_x1_min;				// viscosity at the x1- face of the cell
      double viscosity_x1_pls;				// viscosity at the x1+ face of the cell
      double viscosity_x2_min;				// viscosity at the x2- face of the cell
      double viscosity_x2_pls;				// viscosity at the x2+ face of the cell
      double viscosity_x3_min;				// viscosity at the x3- face of the cell
      double viscosity_x3_pls;				// viscosity at the x3+ face of the cell
      double density_cell_center;			// scaled density at the center of the cell
      int one_dimensional_index;			// index of point in 1-D array
		  
	    
	    
      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
	      {
		  /* compute the level-set value at the center of the control volume */
		  /* and at the centers of the 6 faces */
		  
// 		  level_set_cell_center=0.5*(
// 					level_set[i_index  ][j_index  ][k_index  ]+
// 					level_set[i_index  ][j_index  ][k_index+1]
// 					         );
		  level_set_face_x1_min=0.25*(
					level_set[i_index-1][j_index  ][k_index  ]+
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index-1][j_index  ][k_index+1]+
					level_set[i_index  ][j_index  ][k_index+1]
					      );
		  level_set_face_x1_pls=0.25*(
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index  ][j_index  ][k_index+1]+
					level_set[i_index+1][j_index  ][k_index+1]+
					level_set[i_index+1][j_index  ][k_index  ]
					      );
		  level_set_face_x3_min=level_set[i_index  ][j_index  ][k_index  ];
		  level_set_face_x3_pls=level_set[i_index  ][j_index  ][k_index+1];
		  level_set_face_x2_min=0.25*(
					level_set[i_index  ][j_index  ][k_index+1]+
					level_set[i_index  ][j_index-1][k_index+1]+
					level_set[i_index  ][j_index  ][k_index  ]+
 					level_set[i_index  ][j_index-1][k_index  ]
 					      );
		  level_set_face_x2_pls=0.25*(
					level_set[i_index  ][j_index  ][k_index  ]+
					level_set[i_index  ][j_index+1][k_index  ]+
					level_set[i_index  ][j_index  ][k_index+1]+
					level_set[i_index  ][j_index+1][k_index+1]
					      );
                 density_cell_center=scaled_density_u3[i_index][j_index][k_index];
// 		  density_cell_center=compute_scaled_density(level_set_cell_center, level_set_cell_center,
// 					 rho_plus_over_rho_minus);
		  
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
		  
		  /* evaluate the density at the center of the control volume */
		  
		  one_dimensional_index=map_index_u3(i_index,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);
		  
		  momentum_matrix_u3[0][one_dimensional_index]= density_cell_center+
			actual_time_step_navier_stokes*	(
			      (viscosity_x1_min+viscosity_x1_pls)*one_over_dx1_squared+
			      (viscosity_x2_min+viscosity_x2_pls)*one_over_dx2_squared+
			      (viscosity_x3_min+viscosity_x3_pls)*one_over_dx3_squared
							);
			      
		    
		  momentum_matrix_u3[1][one_dimensional_index]=
			    -1.0*actual_time_step_navier_stokes*viscosity_x1_pls*
				one_over_dx1_squared;
		  momentum_matrix_u3[2][one_dimensional_index]=
			    -1.0*actual_time_step_navier_stokes*viscosity_x2_pls*
				one_over_dx2_squared;
		  momentum_matrix_u3[3][one_dimensional_index]=
			    -1.0*actual_time_step_navier_stokes*viscosity_x3_pls*
				one_over_dx3_squared;

	/* later on the momentum equation will be made nonsymmetric */
	/* so now we already use it for convencience */

		  momentum_matrix_u3[4][one_dimensional_index]=
			    -1.0*actual_time_step_navier_stokes*viscosity_x1_min*
				one_over_dx1_squared;
		  momentum_matrix_u3[5][one_dimensional_index]=
			    -1.0*actual_time_step_navier_stokes*viscosity_x2_min*
				one_over_dx2_squared;
		  momentum_matrix_u3[6][one_dimensional_index]=
			    -1.0*actual_time_step_navier_stokes*viscosity_x3_min*
				one_over_dx3_squared;
			    
		
	      }  
	  }  
       }  
  }

#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/*  Function to smooth the computed curvature distribution                      */
/*            										*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* Instead of pre-smoothing the level-set field, before computing the curvature,*/
/* Sander has chosen to smooth the curvature itself, because it is the non-     */
/* smoothness in the curvature that leads to undesirable effects in the         */
/* solution.										*/
/********************************************************************************/
void smooth_curvature(
  	    double ***level_set, 				// level set field
	    double ***curvature_new,				// interface curvature
	    double ***unsmoothed_curvature,			// interface curvature after smoothing
	    int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	    int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	    int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	    double mesh_width_x1,				// grid spacing in x1 direction (uniform)
	    double mesh_width_x2,				// grid spacing in x2 direction (uniform)
	    double mesh_width_x3,				// grid spacing in x3 direction (uniform)
	    int apply_curvature_smoothing_filter,		// =1, apply curvature smoothing filter
	    							// =0, do not apply curvature smoothing filter
	    int number_curvature_smoothing_steps		// number of curvature smoothing steps
  
       )
	{
    
      double ***double_Matrix2(
	    int number_primary_cells_i,		// allocate memory for a three-
	    int number_primary_cells_j, 		// dimensional array of doubles
	    int number_primary_cells_k
	    );
      void free_double_Matrix2( 			// deallocate memory for a three
	    double ***doubleMatrix2, 			// dimensional array of doubles
	    int number_primary_cells_i,	
	    int number_primary_cells_j);
      void copy_cell_centered_field( 		// copy cell centered field
	    double ***source_field, 		   			
	    double ***target_field,		
	    int number_primary_cells_i,		
	    int number_primary_cells_j,		
	    int number_primary_cells_k		
	    );
      double curvature_filter(			// the curvature filter modifies
	    double level_set_value,			// the diffusion coefficient of
	    double mesh_width_x1,			// the filtering equation, so only
	    double mesh_width_x2,			// a small band is affected
	    double mesh_width_x3			
		   );
      double ***curvature_old;			// curvature at the previous time step in the 
							// curvature smoothing algorithm
      double time_step_curvature_smoothing;		// time-step in the curvature smoothing algorithm
      double filter_value_i_min;			// value of the filter, between cell and its i- neighbour
      double filter_value_i_plus;			// value of the filter, between cell and its i+ neighbour		
      double filter_value_j_min;			// value of the filter, between cell and its j- neighbour
      double filter_value_j_plus;			// value of the filter, between cell and its j+ neighbour
      double filter_value_k_min;			// value of the filter, between cell and its k- neighbour
      double filter_value_k_plus;			// value of the filter, between cell and its k+ neighbour
      int i_index, j_index, k_index;  		// local variables for loop indexing
      int curvature_smoothing_step_index;       	// index of the time step of the curvature smoothing
							// algorithm
     /* save the unsmoothed curvature for inspection/debugging */
      
      copy_cell_centered_field(unsmoothed_curvature,  curvature_new, 
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
      /* allocate memory for a second time level for the curvature */
      /* and copy the curvature to this second 'old' time level */
      
      curvature_old=double_Matrix2(number_primary_cells_i+2,
				   number_primary_cells_j+2, number_primary_cells_k+2);
      copy_cell_centered_field(curvature_new, curvature_old, 
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
       
      /* the time step is fixed, because the maximum diffusion coefficient is fixed */
    
      time_step_curvature_smoothing=1.0/6.0;
      
      std::cerr<<" time_step_curvature_smoothing "<<time_step_curvature_smoothing<<" \n";
      
      
      if(!apply_curvature_smoothing_filter)
      {
	      
	/* NO use is made of the additional filter for the curvature smoothing */
	  for(curvature_smoothing_step_index=0; 
		curvature_smoothing_step_index< number_curvature_smoothing_steps;
		      curvature_smoothing_step_index++)
	  {
	      /* shift the field to the old time level */
	      
	      copy_cell_centered_field(curvature_new, curvature_old, 
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
	      
	      /* do one smoothing step */
      
	      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)   
	      {
		  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		  {
		      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		      {
			curvature_new[i_index][j_index][k_index]=
			      curvature_old[i_index][j_index][k_index]+ 
				time_step_curvature_smoothing*(
			(curvature_old[i_index+1][j_index  ][k_index  ]-
			 curvature_old[i_index  ][j_index  ][k_index  ])+
			(curvature_old[i_index-1][j_index  ][k_index  ]-
			 curvature_old[i_index  ][j_index  ][k_index  ])+
			(curvature_old[i_index  ][j_index+1][k_index  ]-
			 curvature_old[i_index  ][j_index  ][k_index  ])+
			(curvature_old[i_index  ][j_index-1][k_index  ]-
			 curvature_old[i_index  ][j_index  ][k_index  ])+
			(curvature_old[i_index  ][j_index  ][k_index+1]-
			 curvature_old[i_index  ][j_index  ][k_index  ])+
			(curvature_old[i_index  ][j_index  ][k_index-1]-
			 curvature_old[i_index  ][j_index  ][k_index  ]));
		      }
		  }
	      }
	  }
      }
      else
      {
	/* use is made of the additional filter for the curvature smoothing */
      }
	  for(curvature_smoothing_step_index=0; 
		curvature_smoothing_step_index< number_curvature_smoothing_steps;
		      curvature_smoothing_step_index++)
	  {
	      /* shift the field to the old time level */
	      
	      copy_cell_centered_field(curvature_new, curvature_old, 
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
	      
	      /* do one smoothing step */
      
	      for(i_index=1;i_index<number_primary_cells_i;i_index++)   
	      {
		  for(j_index=1;j_index<number_primary_cells_j;j_index++)
		  {
		      for(k_index=1;k_index<number_primary_cells_k;k_index++)
		      {
			filter_value_i_plus =curvature_filter(
				    0.5*(level_set[i_index+1][j_index][k_index]-
					level_set[i_index][j_index][k_index]),
						mesh_width_x1, mesh_width_x2, mesh_width_x3);
			filter_value_i_min  =curvature_filter(
				    0.5*(level_set[i_index-1][j_index][k_index]-
					level_set[i_index][j_index][k_index]),
						mesh_width_x1, mesh_width_x2, mesh_width_x3);
			filter_value_j_plus =curvature_filter(
				    0.5*(level_set[i_index][j_index+1][k_index]-
					level_set[i_index][j_index][k_index]),
						mesh_width_x1, mesh_width_x2, mesh_width_x3);
			filter_value_j_min  =curvature_filter(
				    0.5*(level_set[i_index][j_index-1][k_index]-
					level_set[i_index][j_index][k_index]),
						mesh_width_x1, mesh_width_x2, mesh_width_x3);
			filter_value_k_plus =curvature_filter(
				    0.5*(level_set[i_index][j_index][k_index+1]-
					level_set[i_index][j_index][k_index]),
						mesh_width_x1, mesh_width_x2, mesh_width_x3);
			filter_value_k_min  =curvature_filter(
				    0.5*(level_set[i_index][j_index][k_index-1]-
					level_set[i_index][j_index][k_index]),
						mesh_width_x1, mesh_width_x2, mesh_width_x3); 
			curvature_new[i_index][j_index][k_index]=
			      curvature_old[i_index][j_index][k_index]+
				time_step_curvature_smoothing*(
			filter_value_i_plus*
			  (curvature_old[i_index+1][j_index][k_index]-
			    curvature_old[i_index][j_index][k_index])+
			filter_value_i_min*
			  (curvature_old[i_index-1][j_index][k_index]-
			    curvature_old[i_index][j_index+1][k_index])+
			filter_value_j_plus*
			  (curvature_old[i_index][j_index+1][k_index]-
			    curvature_old[i_index][j_index][k_index])+
			filter_value_j_min*
			  (curvature_old[i_index][j_index-1][k_index]-
			    curvature_old[i_index][j_index][k_index])+
			filter_value_k_plus*
			  (curvature_old[i_index][j_index][k_index+1]-
			    curvature_old[i_index][j_index][k_index])+
			filter_value_k_min*
			  (curvature_old[i_index][j_index][k_index-1]-
			    curvature_old[i_index][j_index][k_index]));
		      }
		  }
	      }
	  }
// 	/* check max difference between smoothed and unsmoothed_curvature */
// 	
// 	      double max_curv_difference=0;
// 	      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)   
// 	      {
// 		  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
// 		  {
// 		      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
// 		      {
// 		       max_curv_difference=std::max(fabs(
// 			curvature_new[i_index][j_index][k_index]-
// 			      unsmoothed_curvature[i_index][j_index][k_index]), max_curv_difference); 
// 		      }
// 		  }
// 	      }
// 	
// 		std::cerr<<"max curv difference "<<max_curv_difference<<"\n";
	
	
	free_double_Matrix2(curvature_old, number_primary_cells_i+2,number_primary_cells_j+2);
  }
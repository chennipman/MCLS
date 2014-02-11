#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/*  Function to initialize the volume of fluid field                            */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
      void initialize_volume_of_fluid(
	  double ***level_set, 			// level set field 
	  double ***volume_of_fluid, 			// volume of fluid field
	  int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	  int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	  int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
	  double lower_bound_derivatives 		// lower bound for the first partial derivatives
							// to consider it a limiting case of vanishing
							// partial derivatives
	    )
      {
      /* function definitions */
      
      int    compute_volume_of_fluid(		// compute volume of fluid field
		double ***level_set, 		// corresponding to a given 
		double ***d_level_set_d_x1, 		// level-set field
		double ***d_level_set_d_x2, 
		double ***d_level_set_d_x3,
		double ***volume_of_fluid,
		int number_primary_cells_i, 
		int number_primary_cells_j, 
		int number_primary_cells_k,
		double lower_bound_derivatives
		);
      void 	compute_level_set_gradient(		// compute gradient of level-set field
		double ***level_set_star, 
		double ***d_level_set_d_x1, 
		double ***d_level_set_d_x2, 
		double ***d_level_set_d_x3,
		int number_primary_cells_i, 
		int number_primary_cells_j, 
		int number_primary_cells_k
		);
      double ***double_Matrix2(			// allocate memory for a three-dimensional array of doubles
		int number_primary_cells_i,				
		int number_primary_cells_j, 		
		int number_primary_cells_k
		);
      void   free_double_Matrix2( 			// deallocate memory for a three
		double ***doubleMatrix2, 		// dimensional array of doubles
		int number_primary_cells_i,	
		int number_primary_cells_j
		);
      void field_neumann_boundary(			// apply neumann boundary condition to
	  	double ***field, 			// cell centered field
	  	int number_primary_cells_i,	
	  	int number_primary_cells_j,	
	  	int number_primary_cells_k	
	  	);
      	void  field_extrapolate_boundary(      	// extrapolate field to virtual cells
            double ***field, 			
            int number_primary_cells_i,	
            int number_primary_cells_j,	
            int number_primary_cells_k	
	    );

      double ***d_level_set_d_x1;			// first partial derivative of
							// the level-set field wrt x1
							// second order central approximation
      double ***d_level_set_d_x2;			// first partial derivative of 
							// the level-set field wrt x2
							// second order central approximation
      double ***d_level_set_d_x3;			// first partial derivative of
 							// the level-set field wrt x3
							// second order central approximation



	/* allocate  memory for the gradient of the level-set field */
	
	d_level_set_d_x1=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x2=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x3=double_Matrix2(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	 
	
	/* compute the gradient of the level-set field, necessary for the level-set/vof conversion */
	
	compute_level_set_gradient(level_set, d_level_set_d_x1, d_level_set_d_x2,
									  d_level_set_d_x3,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k					  
				  );


	/* compute the volume of fluid field */
	
	if(compute_volume_of_fluid(level_set, d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
				volume_of_fluid,
				number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
				lower_bound_derivatives
				  ))
	{
	  std::cerr << "**************************************************** \n";
	  std::cerr<<"ERROR \n";
	  std::cerr<<"Something went wrong in the level-set to volume of fluid conversion \n";
	  std::cerr<<"Aborting...\n";
	  std::cerr<<"in function initialize_volume_of_fluid line 97 \n";
	  std::cerr << "**************************************************** \n";
	};
	
    /* extend the volume of fluid field to the virtual cells */
    
//     	field_neumann_boundary(volume_of_fluid, number_primary_cells_i, 
// 			   number_primary_cells_j,number_primary_cells_k);	
     	field_extrapolate_boundary(volume_of_fluid, number_primary_cells_i, 
			   number_primary_cells_j,number_primary_cells_k);	
	
	/* de-allocate  memory for the gradient of the level-set field */
	
	free_double_Matrix2(d_level_set_d_x1, number_primary_cells_i+2, number_primary_cells_j+2);
	free_double_Matrix2(d_level_set_d_x2, number_primary_cells_i+2, number_primary_cells_j+2);
	free_double_Matrix2(d_level_set_d_x3, number_primary_cells_i+2, number_primary_cells_j+2);
      }
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
class vector
{
public:
  double u1,u2,u3;
  vector(double u1, double u2, double u3);
  vector( void);
};
/********************************************************************************/
/********************************************************************************/
/*  Function to set the inhomogeneous boundary condition for the pressure       */
/*                                                                              */
/*                                                                              */
/*  Programmer       : Duncan van der Heul                                      */
/*  Date      : 10-03-2013                                                      */
/*  Update    :                                                                 */
/********************************************************************************/
/* Notes                                                                        */
/* At the moment only an inhomogeneous boundary condition can be specified      */
/* but later this will be extended to more advanced boundary conditions.        */
/********************************************************************************/
  void set_pressure_boundary_condition(
       vector gravity,				       		// gravitational acceleration vector 
       double ***pressure_boundary_condition_x1,          	// inhomogeneous boundary condition for
						              	// the pressure planes with normal in x1 direction
       double ***pressure_boundary_condition_x2,          	// inhomogeneous boundary condition for the pressure
						              	// the pressure planes with normal in x1 direction
       double ***pressure_boundary_condition_x3,          	// inhomogeneous boundary condition for the pressure
						              	// the pressure planes with normal in x1 direction
       int number_primary_cells_i,	                     	// number of primary (pressure) cells in x1 direction
       int number_primary_cells_j,	                     	// number of primary (pressure) cells in x2 direction
       int number_primary_cells_k	                     	// number of primary (pressure) cells in x3 direction
 	)
    {
      void set_constant_matrix(			       		// set 2 dimensional array to constant
	    int first_dimension,		              	// value
	    int second_dimension,	
	    double **matrix_to_set,	
	    double constant_value	
       );

       int i_index, j_index, k_index;  				// local variables for loop indexing
       
       /* set the inhomogeneous neumann boundary condition for the pressure to the x1 */
       /* component of the gravitational acceleration vector */
       /* -1 for the direction of the outward normal vector for i=1 slice*/
       /* +1 for the direction of the outward normal vector for i=number_primary_cells_i slice*/
       
       // set_constant_matrix(number_primary_cells_j+2, number_primary_cells_k+2,  pressure_boundary_condition_x1[0], -1.0*gravity.u1);
       // set_constant_matrix(number_primary_cells_j+2, number_primary_cells_k+2,  pressure_boundary_condition_x1[1],  1.0*gravity.u1);

  	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
  	{
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
			pressure_boundary_condition_x1[0][j_index][k_index]=-1.0*gravity.u1;

      		}  
	}
  	for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
  	{
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
			pressure_boundary_condition_x1[1][j_index][k_index]= 1.0*gravity.u1;

      		}  
	}
       
       /* set the inhomogeneous neumann boundary condition for the pressure to the x2 */
       /* component of the gravitational acceleration vector */
       /* -1 for the direction of the outward normal vector for j=1 slice*/
       /* +1 for the direction of the outward normal vector for j=number_primary_cells_j slice*/

       // set_constant_matrix(number_primary_cells_i+2, number_primary_cells_k+2,  pressure_boundary_condition_x2[0], -1.0*gravity.u2);
       // set_constant_matrix(number_primary_cells_i+2, number_primary_cells_k+2,  pressure_boundary_condition_x2[1],  1.0*gravity.u2);
       
 	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
 	{
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
			pressure_boundary_condition_x2[0][i_index][k_index]=-1.0*gravity.u2;

     		}  
	}
	
 	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
 	{
		for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		{
			pressure_boundary_condition_x2[1][i_index][k_index]= 1.0*gravity.u2;

     		}  
	}
        
       /* set the inhomogeneous neumann boundary condition for the pressure to the x3 */
       /* component of the gravitational acceleration vector */
       /* -1 for the direction of the outward normal vector for k=1 slice*/
       /* +1 for the direction of the outward normal vector for k=number_primary_cells_k slice*/

       // set_constant_matrix(number_primary_cells_i+2, number_primary_cells_j+2,  pressure_boundary_condition_x3[0], -1.0*gravity.u3);
       // set_constant_matrix(number_primary_cells_i+2, number_primary_cells_j+2,  pressure_boundary_condition_x3[1],  1.0*gravity.u3);

	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	{
		for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		{
			pressure_boundary_condition_x3[0][i_index][j_index]=-1.0*gravity.u3;
    		}  
	}

	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	{
		for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		{
			pressure_boundary_condition_x3[1][i_index][j_index]= 1.0*gravity.u3;
    		}  
	}


}
       
       
       
	
#include "../headers/array.h"


#include<cstdlib>
#include<iostream>
#include<algorithm>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the density at the controlvolumes of the u3 velocity    */
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/********************************************************************************/
	void compute_scaled_density_u3(
	Array3<double> scaled_density_u3,			// scaled density for the controlvolumes
							// of the momentum equation in x1 direction
	Array3<double> volume_of_fluid_u3,			// volume of fluid value for the controlvolumes
							// of the momentum equation in x1 direction
       double rho_plus_over_rho_minus,		// ratio of the densities of the two phases
	int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
	int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
	int number_primary_cells_k			// number of primary (pressure) cells in x3 direction
	)
{
      	int i_index, j_index, k_index;  		// local variables for loop indexing


	/* compute the density for the control volumes of the u3 velocity  using local vof values */
	
	
	for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
  	{		
		for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
      		{			
	  		for(k_index=0;k_index<number_primary_cells_k+1;k_index++)
	  		{
				if(k_index==0 || k_index==number_primary_cells_k)
				{
 				scaled_density_u3[i_index][j_index][k_index]=
 						2.0*volume_of_fluid_u3[i_index][j_index][k_index]*rho_plus_over_rho_minus+
 						(1.0-2.0*volume_of_fluid_u3[i_index][j_index][k_index]);
 							
				}
				else
				{
				scaled_density_u3[i_index][j_index][k_index]=
 						volume_of_fluid_u3[i_index][j_index][k_index]*rho_plus_over_rho_minus+
 						(1.0-volume_of_fluid_u3[i_index][j_index][k_index]);
 				}
	  		}  
  
      		}  
     
  	} 
}	

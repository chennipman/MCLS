#include "../headers/array.h"
 
#include<cstdlib>
#include<iostream>
     
/********************************************************************************/
/*  Function to modify the volume of fluid field.					*/
/*  											*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/* If a cell is not an interface cell (it is surrounded by 6 neighbours with    */
/* a level-set value of the same sign, the volume of fluid value can not        */
/* have an intermediate value. It has to be either 1 or 0.                      */
/* Furthermore, if a cell has a volume of fluid value less than 0 or larger     */
/* than 1, this is a nonphysical situation                                      */
/* Both can be corrected, but to have global mass conservation the correction   */
/* has to be added to or removed from cells that allow their volume of fluid    */
/* value to be lowered or raised. Only interface cells qualify.                 */
/* This function determines the errors and the number of cells where errors     */
/* occur. The volume of fluid in the erroneous cells is set to the correct      */
/* value.										*/
/********************************************************************************/
EXPORT int modify_volume_of_fluid_values(
      Array3<double> level_set, 				    	// level-set field
      Array3<double> volume_of_fluid,					// volume of fluid field, uncorrected
									// so with possible vapour cells and 
									// under/overfilled cells
      Array3<double> volume_of_fluid_correction,		    	// correction to the volume of fluid field
									// to make it valid
      int number_primary_cells_i,			    		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			    	    	// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			            	// number of primary (pressure) cells in x3 direction
      double volume_of_fluid_tolerance			    	        // tolerance for volume of fluid value
      )
      {
//       double mean_value_level_set;					// average of left and right hand side value
// 									// of the level-set field
      double current_level_set_value;				        // current level-set field value in the control volume
      double current_volume_of_fluid_value;				// current volume of fluid field value in 
									// the control volume 
      double correct_vof_value;					        // the value a cell with an invalid volume of
									// fluid is set to: 1, 0 dependent on the sign
									// of the level-set function
      int number_cells_vof_out_of_bounds=0;				// number of control volumes where the volume of fluid
									// function is OUTSIDE the interval [0,1]
      int number_cells_numerical_vapor=0;				// number of control volumes where the volume of fluid
									// function is INSIDE the interval [0,1]
									// while the cell has 6 neighbours with the same sign:
									// the cell is NOT an interface cell
      int number_cells_invalid_volume_of_fluid; 			// sum of number of vapour cells and number of cells
      int i_index, j_index, k_index;  				        // local variables for loop indexing
      int adapt_vof_value;						// =1 adapt the volume of fluid value in this
									// cell, because it has an invalid value


         /* initial check on the value of the volume of fluid values */									
  
	  std::cerr<<"in modify_volume_of_fluid_values \n";
         
         number_cells_vof_out_of_bounds=check_volume_of_fluid(volume_of_fluid, number_primary_cells_i,
                                           number_primary_cells_j,
                                            number_primary_cells_k,
                                               volume_of_fluid_tolerance);

         /* all nonvirtual cells are modified */
         /* the required modification is stored in volume_of_fluid_correction */
         
         std::cerr<<"start modifying loop \n";
         
	  for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
	  {
	      for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	      {
		  for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
		  {
		      current_level_set_value=level_set[i_index][j_index][k_index];
		      current_volume_of_fluid_value=volume_of_fluid[i_index][j_index][k_index];

                    
		    adapt_vof_value=0;
		      
		      /* check if the cell is a vapour cell, when this update is applied */
		    if(
			/* is the cell surrounded by neighbours with the same sign ?*/
		      (
		      (level_set[i_index-1][j_index  ][k_index  ]*current_level_set_value>0)&&
		      (level_set[i_index+1][j_index  ][k_index  ]*current_level_set_value>0)&&
		      (level_set[i_index  ][j_index-1][k_index  ]*current_level_set_value>0)&&
		      (level_set[i_index  ][j_index+1][k_index  ]*current_level_set_value>0)&&
		      (level_set[i_index  ][j_index  ][k_index+1]*current_level_set_value>0)&&
		      (level_set[i_index  ][j_index  ][k_index-1]*current_level_set_value>0))
			    &&
			/* is the volume of fluid at the same time in the OPEN interval <0,1>? */
		      (current_volume_of_fluid_value<1.0-volume_of_fluid_tolerance ||
			  current_volume_of_fluid_value>0.0+volume_of_fluid_tolerance)
		    )
		    {
		    /* then the cell is a vapour cell, and the volume_of_fluid_correction is not correct yet */
		    /* additional iterations are required */
		    /* update the number of vapour cells */
		    
			  number_cells_numerical_vapor++;
			  adapt_vof_value=1;
                     
			  std::cerr<<"vapour cell "<< i_index<<" "<<j_index<<" "<<k_index<<"\n";
			  std::cerr<<"value "<< current_volume_of_fluid_value <<"\n";
		    }
		    if( (current_volume_of_fluid_value> 1.0+volume_of_fluid_tolerance)||
			    (current_volume_of_fluid_value < 0.0-volume_of_fluid_tolerance))
		    {
		      /* the volume of fluid value is out of bounds */
		      /* additional iterations are required */
		      /* update the number of out of bounds cells */
			  
			  number_cells_vof_out_of_bounds++;
			  adapt_vof_value=1;
                       
			  std::cerr<<"out of bounds cell "<< i_index<<" "<<j_index<<" "<<k_index<<"\n";
			  std::cerr<<"value "<< current_volume_of_fluid_value <<"\n";
		    }  
		    if(adapt_vof_value)
		    {
		      /* cell is either a vapour cell or an over/under filled cell   */
		      /* the correct value shoul be either 1/0 depending on the sign */
		      /* of the level-set field
                     */
			correct_vof_value=0.5+sign(0.5,current_level_set_value);
			if(current_level_set_value==0.0)
			{
                            
		      /* this is not really very likely to happen but still */

                            correct_vof_value=0.5;
			}
			
			volume_of_fluid_correction[i_index][j_index][k_index]=volume_of_fluid[i_index][j_index][k_index]-
				      correct_vof_value;
			volume_of_fluid[i_index][j_index][k_index]=correct_vof_value;
				      
		    }
		 } 
  
	      }  
     
	  }
	  
	  number_cells_invalid_volume_of_fluid=number_cells_numerical_vapor+number_cells_vof_out_of_bounds;

	  return number_cells_invalid_volume_of_fluid;
    }

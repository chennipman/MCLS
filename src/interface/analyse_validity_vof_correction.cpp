#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
/********************************************************************************/
/********************************************************************************/
/*  Function to establish the validity of the correction to                     */
/*  the volume of fluid field             					*/
/*                                                                              */
/*  Programmer  : Duncan van der Heul                                           */
/*  Date        : 10-03-2013                                                    */
/*  Update      :                                                               */
/********************************************************************************/
/* Notes                                                                        */
/*    Check for cells with a corrected volume of fluid outside the              */
/*    interval [-eps,1+eps] , with eps the volume of fluid tolerance            */
/*    and 'vapor' cells. The correction is 'tried out' to see if it would       */
/*    lead to a valid volume of fluid field. In that case the redistribution    */
/*    algorithm can be terminated.                                              */
/*                                                                              */
/********************************************************************************/
EXPORT void analyse_validity_vof_correction(
        Array3<double> level_set,                                       // level set field 
                                                                        // mass conserving
        Array3<double> volume_of_fluid_correction,                      // correction to the volume of fluid field
                                                                        // to make it valid
        Array3<double> volume_of_fluid,                                 // volume of fluid field
        
        Array3<double> invalid_vof_cells,                               // indication field to show what is wrong
                                                                        // indicator field showing cells that are either
                                                                        // within bounds = 0
                                                                        // underfilled   =-1
                                                                        // overfilled    =+1
                                                                        // vapour cells  = 5
        int number_primary_cells_i,                                     // number of primary (pressure) cells in x1 direction
        int number_primary_cells_j,                                     // number of primary (pressure) cells in x2 direction
        int number_primary_cells_k,                                     // number of primary (pressure) cells in x3 direction
        double volume_of_fluid_tolerance,                               // tolerance for volume of fluid value
        int &number_cells_vof_out_of_bounds,                            // number of control volumes where the volume of fluid
                                                                        // function is OUTSIDE the interval [0,1]
        int &number_cells_numerical_vapor,                              // number of control volumes where the volume of fluid
                                                                        // function is INSIDE the interval [0,1]
                                                                        // while the cell has 6 neighbours with the same sign:
                                                                        // the cell is NOT an interface cell
        int &number_cells_invalid_volume_of_fluid                       // sum of number of vapour cells and number of cells
                                                                        // with the volume of fluid outside [0,1];
)
{        
        double level_set_value;                                         // level set field in cell center
        double current_volume_of_fluid;                                 // corrected, cell centered value of volume of fluid
        int pure_cell;                                                   // true: all vertices on one side of interface,
                                                                        // false: not all vertices on one side of interface
        int i_index, j_index, k_index;                                  // local variables for loop indexing

         /* reset all numbers */
        
        number_cells_invalid_volume_of_fluid=0;
        number_cells_numerical_vapor=0;
        number_cells_vof_out_of_bounds=0;
          
            for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
                for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
                    for(k_index=1;k_index<number_primary_cells_k+1;k_index++){

                     /* start with the assumption all cells have valid vof values */
                     
                        invalid_vof_cells[i_index][j_index][k_index]=0;

                     /* these are the current values for level set and  */
                     /* volume of fluid for this cell, note that the    */
                     /* current volume of fluid cell is not stored in   */
                     /* the array volume of fluid until the algorithm   */
                     /* actually is converged.                          */
                     
                        level_set_value =level_set[i_index][j_index][k_index];
                        current_volume_of_fluid=volume_of_fluid[i_index][j_index][k_index]+
                                   volume_of_fluid_correction[i_index][j_index][k_index];

                        /* two conditions are checked to verify if the volume of fluid value is valid :    */
                        /*1) if this is a NON-interface cell: the volume of fluid should also be either 0 or 1 */
                        /*                              with some tolerance                                */
                        /*2) if this is a mixed cell: the volume of fluid should be in the interval 0 to 1 */
                        /*                              with some tolerance                                */
                        
                        /* first condition is checked */
                        
                        if(current_volume_of_fluid>=volume_of_fluid_tolerance && 
                        	current_volume_of_fluid<=1.0-volume_of_fluid_tolerance)
                        {
                        
                        	/* is the cell a cell removed from the interface ? */
                        	
                        	pure_cell=determine_cell_is_pure(level_set, i_index, j_index, k_index,						
                        					number_primary_cells_i,                             
                        					  number_primary_cells_j,                             
                        					  	number_primary_cells_k); 
                        
                        	if(pure_cell)
                                {

                             /*  this cell is a vapor cell */
                             
                                 	number_cells_numerical_vapor++;
                                        invalid_vof_cells[i_index][j_index][k_index]=5;
                             	}
                        
                        }
                        
                        /* second condition is checked */
                        
                        if(current_volume_of_fluid > 1.0+volume_of_fluid_tolerance||
                                            current_volume_of_fluid<-1.0*volume_of_fluid_tolerance)
                        {
                              /* this cell has an invalid volume of fluid value */
                        
                              number_cells_vof_out_of_bounds++;
                              
                              if(current_volume_of_fluid > 1.0+volume_of_fluid_tolerance)
                              {
                                     
                                        /* mark this cell as an overfilled cell */
                                        
                                     invalid_vof_cells[i_index][j_index][k_index]=1; 
                              }
                              else
                              {
                                        /* mark this cell as an underfilled cell */
                                        
                                    invalid_vof_cells[i_index][j_index][k_index]=-1; 
                              }
                              
                        }
                       
                     }
                }
            }
                number_cells_invalid_volume_of_fluid=number_cells_numerical_vapor+number_cells_vof_out_of_bounds;
                std::cerr<<"After this sweep the correction would lead to:"<<"\n";	
                std::cerr<<"number_cells_numerical_vapor "<< number_cells_numerical_vapor<<"\n" ;
  		std::cerr<<"number_cells_vof_out_of_bounds "<< number_cells_vof_out_of_bounds<<"\n" ;
}
        

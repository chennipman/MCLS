#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
/********************************************************************************/
/********************************************************************************/
/*  Function to make the volume of fluid field valid                            */
/*  method.                                                                     */
/*                                                                              */
/*  Programmer  : Duncan van der Heul                                           */
/*  Date        : 10-03-2013                                                    */
/*  Update      :                                                               */
/********************************************************************************/
/* Notes                                                                        */
/*                                                                              */
/* All the overshoots and undershoots of the volume of fluid field are clipped  */
/* and the vapor cells are removed. This will provide the basis for the         */
/* correct level-set field.                                                     */
/********************************************************************************/
EXPORT void make_vof_field_valid(
        Array3<double> level_set,                                       // level set field 
                                                                        // mass conserving
        Array3<double> volume_of_fluid,                                 // volume of fluid field
        
        Array3<double> invalid_vof_cells,                               // indication field to show what is wrong
                                                                        // indicator field showing cells that are either
                                                                        // within bounds =0
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
        double current_volume_of_fluid;                                 // uncorrected, cell centered value of volume of fluid
        int pure_cell;                                                  // true: all vertices on one side of interface,
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
                     /* volume of fluid for this cell                   */
                      
                        level_set_value =level_set[i_index][j_index][k_index];
                        current_volume_of_fluid=volume_of_fluid[i_index][j_index][k_index];

                        
                        
                        /* two conditions are checked to verify if the volume of fluid value is valid :    */
                        /*1) if this is a NON-interface cell: the volume of fluid should also be either 0 or 1 */
                        /*                              with some tolerance                                */
                        /*2) if this is a mixed cell: the volume of fluid should be in the interval 0 to 1 */
                        /*                              with some tolerance                                */
                        
                        /* first condition is checked */
                        
                        
                        /* a cell can only be a vapor cell if the vof value is in <eps,1-eps> */
                        
                        if(current_volume_of_fluid>=volume_of_fluid_tolerance && 
                        	current_volume_of_fluid<=1.0-volume_of_fluid_tolerance)
                        {
                        
                        	pure_cell=determine_cell_is_pure(level_set, i_index, j_index, k_index,						
                        					number_primary_cells_i,                             
                        					  number_primary_cells_j,                             
                        					  	number_primary_cells_k); 
                        
                        	if(pure_cell)
                                {

                             /*  this cell is a vapor cell */
                             
                                 	number_cells_numerical_vapor++;
                                 	volume_of_fluid[i_index][j_index][k_index]=0.5+sign(0.5,level_set_value);
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
                                     volume_of_fluid[i_index][j_index][k_index]=1;
                             }
                              else
                              {
                                        /* mark this cell as an underfilled cell */
                                        
                                    invalid_vof_cells[i_index][j_index][k_index]=-1; 
                                    volume_of_fluid[i_index][j_index][k_index]=0;
                             }
                              
                        }
                                                
                    }
                }
            }
            number_cells_invalid_volume_of_fluid=number_cells_numerical_vapor+number_cells_vof_out_of_bounds;
            std::cerr<<" number_cells_invalid_volume_of_fluid in validation "<< number_cells_invalid_volume_of_fluid<<" \n";
            std::cerr<<" number_cells_numerical_vapor in validation "<< number_cells_numerical_vapor<<" \n";
            std::cerr<<" number_cells_vof_out_of_bounds in validation"<< number_cells_vof_out_of_bounds<<" \n";
        
}
        

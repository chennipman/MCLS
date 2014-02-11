/********************************************************************************/
/********************************************************************************/
/*  Function to copy a cell centered field including virtual cells              */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* when we switch to more efficient definition for the multi-dimensional arrays */
/* this function can be simplified of even discarded				*/
/********************************************************************************/
      void copy_cell_centered_field( 
	    double ***source_field, 		// original field
	    double ***target_field,		// copy of the original field
	    int number_primary_cells_i,		// number of primary (pressure) 
						// cells in x1 direction
	    int number_primary_cells_j,		// number of primary (pressure) 
						// cells in x2 direction
	    int number_primary_cells_k		// number of primary (pressure) 
						// cells in x3 direction
	   )
     {
     int i_index, j_index, k_index;  		// local variables for loop indexing
     
     
          for(i_index=0;i_index<number_primary_cells_i+2;i_index++)
          {
		for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
           	{
           	      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
           	      {
			    target_field[i_index][j_index][k_index]=
				  source_field[i_index][j_index][k_index];
           	      }
           	}
          }
      }		  
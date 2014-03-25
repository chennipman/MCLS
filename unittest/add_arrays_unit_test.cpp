#include "../src/headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to test the add_arrays function					*/
/*                     								*/
/*  Programmer	: Coen Hennipman	       					*/
/*  Date	: 13-03-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* Unit test for add_arrays		 		                        */
/* Input and output can be compared. 	 		                        */
/********************************************************************************/



EXPORT void add_arrays_unit_test()
{
      Array3<double> output_array;			// output array
      Array3<double> input_array_1;			// input_array_1
      Array3<double> input_array_2;			// input_array_2
	int a = 1;
	int b = 3; 
	int i,j,k;

int number_primary_cells_i,number_primary_cells_j,number_primary_cells_k;
number_primary_cells_i = 2;
number_primary_cells_j = 2;
number_primary_cells_k = 2;

      
      output_array.create(number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+2);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, output_array, -1);
			    
      input_array_1.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, input_array_1, 10.0);	
			    
      input_array_2.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, input_array_2, 15.0);			    

add_arrays(output_array,a,input_array_1,b,input_array_2, number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+2);

  for(i=0;i<number_primary_cells_i+2;i++)
  {
     for(j=0;j<number_primary_cells_j+2;j++)
      {
	  for(k=0;k<number_primary_cells_k+2;k++)
	  {
	  printf("i=%i j=%i k=%i \n", i,j,k);
	  printf("output array = %f \n",output_array[i][j][k]);
	  }  
      }  
  } 
}
	
	


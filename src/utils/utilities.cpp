#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the generic matrix vector product                       */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:       							*/
/********************************************************************************/
/* Notes									*/
/* This function computes the product of a given matrix A and a given vector    */
/* x and returns the product in the vector y                                    */
/* Although the velocity components and the pressure unknowns are stored in     */
/* arrays which have different dimensions, this single function can be applied  */
/* in all cases, because the dimensions are transferred to the function as      */
/* arguments.                                                             	*/
/* Note that the matrix is treated as an heptadiagonal matrix, while it         */
/* really is a block pentadiagonal matrix. In the matrix construction routine   */
/* the necessary coefficients are set to zero so this simplification is possible*/
/********************************************************************************/

EXPORT void matrix_vector_product(
int i_dimension,   // number of unknowns in the system in i-direction
int j_dimension,   // number of unknowns in the system in i-direction
int k_dimension,   // number of unknowns in the system in i-direction
Array2<double> A,    // matrix under consideration
Array1<double> x,     // INPUT vector x
Array1<double> y      // OUTPUT vector y such that y=Ax
)
{
  int number_dof_in_slice= 					// number of degrees of freedom in one slice of the domain
	i_dimension*j_dimension; 				// as considered in the matrix
			   
  int total_number_dof=						// total number of degrees of freedom in the whole domain
	i_dimension*j_dimension*k_dimension;	   	
			   					// as considered in the matrix
  int row_index;	   					// index to indicate the row of the matrix 
  
  
			   
      
/* all rows of the matrix have a nonzero coefficient on the main diagonal */
/* first include that contribution 					  */

  for(row_index=0;row_index<total_number_dof;row_index++)
  {
      y[row_index]=A[0][row_index]*x[row_index];
  }
  
/* all but the first and last row have both a connection to the previous */			   
/* and next unknowns							 */
			   
  for(row_index=0;row_index<total_number_dof-1;row_index++)
  {
      y[row_index]  		+=A[1][row_index]*x[row_index+1];
      y[row_index+1]		+=A[1][row_index]*x[row_index];
  }

  for(row_index=0;row_index<total_number_dof-i_dimension;row_index++)
  {
      y[row_index]  		+=A[2][row_index]*x[row_index+i_dimension];
      y[row_index+i_dimension]	+=A[2][row_index]*x[row_index];
  }
			   
  for(row_index=0;row_index<total_number_dof-number_dof_in_slice;row_index++)
  {
      y[row_index]  			+=A[3][row_index]*x[row_index+number_dof_in_slice];
      y[row_index+number_dof_in_slice]	+=A[3][row_index]*x[row_index];
  }
			   
			   
}


/********************************************************************************/
/********************************************************************************/
/*  Function to compute the generic dot product  	                       	*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function computes the inner product of two vectors of length            */
/* vector_length and returns that value		                                */
/********************************************************************************/
EXPORT double dot_product( 
    int vector_length, 			// length of both input vectors
    Array1<double> x, 				// first input vector
    Array1<double> y				// second input vector
     )
{
    int component_index; 		//index to the components of the two vectors
    double inner_product_value=0; 	// the value of the inner product
    
    for( component_index=0;component_index<vector_length;component_index++)
    {
    inner_product_value+=x[component_index]*y[component_index];  
    }
    return inner_product_value;  
  
}
    
/********************************************************************************/
/********************************************************************************/
/*  Function to copy one vector to another vector	                        */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function copies a vector to another vector elementwise                  */
/*                                     		                                */
/********************************************************************************/
EXPORT void copy_vector( 
  int vector_length, 	// length of both vectors
  Array1<double> original_vector, 	// first input vector
  Array1<double> image_vector	// second input vector
  )
{
  int component_index; //index to the components of the two vectors
  for(component_index=0;component_index<vector_length;component_index++)
  {
  image_vector[component_index]=original_vector[component_index];  
  }
}


/********************************************************************************/
/********************************************************************************/
/*  Function to set a vector to a constant value	                        */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function sets a vector to a constant value elementwise                  */
/*                                     		                                */
/********************************************************************************/
EXPORT void set_constant_vector(
    int vector_length,		// length of the vector
    Array1<double> vector_to_set,	// the name of the vector that has to be set
    double constant_value	// the constant value the vector has to be set to
     )
{
      int component_index; //index to the components of the two vectors
      for( component_index=0;component_index<vector_length;component_index++)
      {
	vector_to_set[component_index]=constant_value;
      }
  
}

/********************************************************************************/
/********************************************************************************/
/*  Function to set a three-dimensional array to a constant value	        */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function sets a vector to a constant value elementwise                  */
/*                                     		                                */
/********************************************************************************/
EXPORT void set_constant_matrix2(
    int first_dimension,	// number of elements in first dimension
    int second_dimension,	// number of elements in second dimension
    int third_dimension,	// number of elements in third dimension
    Array3<double> matrix2_to_set,	// the name of the vector that has to be set
    double constant_value	// the constant value the vector has to be set to
     )
{
  int i_index, j_index, k_index;  // local variables for loop indexing

  for(i_index=0;i_index<first_dimension;i_index++)
  {
     for(j_index=0;j_index<second_dimension;j_index++)
      {
	  for(k_index=0;k_index<third_dimension;k_index++)
	  {
	      matrix2_to_set[i_index][j_index][k_index]=constant_value;
	  }  
  
      }  
     
  } 
  
}
/********************************************************************************/
/********************************************************************************/
/*  Function to add two three-dimensional arrays multiplied by a integer        */
/*  										*/
/*  Programmer	: Coen Hennipman	       					*/
/*  Date	: 13-03-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function adds two arrays and multiplies them with a constant value      */
/*                                     		                                */
/********************************************************************************/
EXPORT void add_arrays(
    Array3<double> output_array,	// the output array that is calculated
    double first_constant_value,	// the constant value the vector has to be set to
    Array3<double> first_array,	// the name of the vector that has to be set
    double second_constant_value,	// the constant value the vector has to be set to
    Array3<double> second_array,	// the name of the vector that has to be set
    int first_dimension,		// number of elements in first dimension
    int second_dimension,		// number of elements in second dimension
    int third_dimension			// number of elements in third dimension
     )
{
  int i_index, j_index, k_index;  // local variables for loop indexing

  for(i_index=0;i_index<first_dimension;i_index++)
  {
     for(j_index=0;j_index<second_dimension;j_index++)
      {
	  for(k_index=0;k_index<third_dimension;k_index++)
	  {
	      output_array[i_index][j_index][k_index]=
	      first_constant_value *first_array[i_index][j_index][k_index]+
	      second_constant_value*second_array[i_index][j_index][k_index];
	  }  
        }  
  } 
}
/********************************************************************************/
/********************************************************************************/
/*  Function to add 3 three-dimensional arrays multiplied by a integer       	*/
/*  										*/
/*  Programmer	: Coen Hennipman	       					*/
/*  Date	: 21-04-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function adds 3 arrays and multiplies them with a constant value      	*/
/*                                     		                                */
/********************************************************************************/
EXPORT void add_3_arrays(
    Array3<double> output_array,	// the output array that is calculated
    double first_constant_value,	// the constant value the vector has to be set to
    Array3<double> first_array,		// the name of the vector that has to be set
    double second_constant_value,	// the constant value the vector has to be set to
    Array3<double> second_array,	// the name of the vector that has to be set
    double third_constant_value,	// the constant value the vector has to be set to
    Array3<double> third_array,		// the name of the vector that has to be set
    int first_dimension,		// number of elements in first dimension
    int second_dimension,		// number of elements in second dimension
    int third_dimension			// number of elements in third dimension
     )
{
  int i_index, j_index, k_index;  // local variables for loop indexing

  for(i_index=0;i_index<first_dimension;i_index++)
  {
     for(j_index=0;j_index<second_dimension;j_index++)
      {
	  for(k_index=0;k_index<third_dimension;k_index++)
	  {
	      output_array[i_index][j_index][k_index]=
	      first_constant_value *first_array[i_index][j_index][k_index]+
	      second_constant_value*second_array[i_index][j_index][k_index]+
	      third_constant_value *third_array[i_index][j_index][k_index];
	  }  
        }  
  } 
}
/********************************************************************************/
/********************************************************************************/
/*  Function to add 4 three-dimensional arrays multiplied by a integer       	*/
/*  										*/
/*  Programmer	: Coen Hennipman	       					*/
/*  Date	: 21-04-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function adds 4 arrays and multiplies them with a constant value      	*/
/*                                     		                                */
/********************************************************************************/
EXPORT void add_4_arrays(
    Array3<double> output_array,	// the output array that is calculated
    double first_constant_value,	// the constant value the vector has to be set to
    Array3<double> first_array,		// the name of the vector that has to be set
    double second_constant_value,	// the constant value the vector has to be set to
    Array3<double> second_array,	// the name of the vector that has to be set
    double third_constant_value,	// the constant value the vector has to be set to
    Array3<double> third_array,		// the name of the vector that has to be set
    double fourth_constant_value,	// the constant value the vector has to be set to
    Array3<double> fourth_array,	// the name of the vector that has to be set
    int first_dimension,		// number of elements in first dimension
    int second_dimension,		// number of elements in second dimension
    int third_dimension			// number of elements in third dimension
     )
{
  int i_index, j_index, k_index;  // local variables for loop indexing

  for(i_index=0;i_index<first_dimension;i_index++)
  {
     for(j_index=0;j_index<second_dimension;j_index++)
      {
	  for(k_index=0;k_index<third_dimension;k_index++)
	  {
	      output_array[i_index][j_index][k_index] = 
	      first_constant_value *first_array[i_index][j_index][k_index]+
	      second_constant_value*second_array[i_index][j_index][k_index]+
	      third_constant_value *third_array[i_index][j_index][k_index]+
	      fourth_constant_value*fourth_array[i_index][j_index][k_index];
	  }  
        }  
  } 
}
/********************************************************************************/
/********************************************************************************/
/*  Function to set a two-dimensional array to a constant value	                */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function sets a vector to a constant value elementwise                  */
/*                                     		                                */
/********************************************************************************/
EXPORT void set_constant_matrix(
    int first_dimension,	// number of elements in first dimension
    int second_dimension,	// number of elements in second dimension
    Array2<double> matrix_to_set,	// the name of the vector that has to be set
    double constant_value	// the constant value the vector has to be set to
     )
{
  int i_index, j_index;  // local variables for loop indexing

  for(i_index=0;i_index<first_dimension;i_index++)
  {
     for(j_index=0;j_index<second_dimension;j_index++)
      {
	   matrix_to_set[i_index][j_index]=constant_value;
  
      }  
     
  } 
  
}


/********************************************************************************/
/********************************************************************************/
/*  Function to compute a linear combination of two input vectors with given    */
/*  weights									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function computes a linear combination of 				*/
/* two vectors with given weights						*/
/*                                     		                                */
/********************************************************************************/
EXPORT void linear_combination(
	int vector_length, 		//length of all vectors 
	Array1<double> input_vector_x, 	//input vector 1, named x
	Array1<double> input_vector_y,		//input vector 2, named y
	Array1<double> output_vector_z,	//output vector, named z
					// such that z=x+alpha*y
	double weight_of_y		//weight alpha in the linear combination above
	)
{
      int component_index; //index to the components of the two vectors
      for( component_index=0;component_index<vector_length;component_index++)
      {
	output_vector_z[component_index]=input_vector_x[component_index]+
			    weight_of_y*input_vector_y[component_index];
      }
     
} 

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the norm of an input vector	  			*/
/*  weights									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function computes a the L2 norm of a vector             		*/
/********************************************************************************/
EXPORT double compute_vector_norm(int vector_length, //length of the vector
		     Array1<double> input_vector  //input vector 
	 )
{  
      int component_index; //index to the components of the two vectors
      double norm_vector=0; //the norm of the vector
      for(component_index=0;component_index<vector_length;component_index++)
      {
	norm_vector+=input_vector[component_index]*input_vector[component_index];
      }
      norm_vector=sqrt(norm_vector);
      return norm_vector;
}	
	
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the minimum element in a vector                         */
/*  weights									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function computes a linear combination of 				*/
/* two vectors with given weights						*/
/*                                     		                                */
/********************************************************************************/
EXPORT double minimum_element(int vector_length, //length of the vector
		     Array1<double> input_vector  //input vector 
	 )
{  
      int component_index; //index to the components of the two vectors
      double minimum_element=10.0E10; //the norm of the vector
      for(component_index=0;component_index<vector_length;component_index++)
      {
	minimum_element=std::min(input_vector[component_index],minimum_element);
      }
      return minimum_element;
}	
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the maximum element in a vector                         */
/*  weights									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function computes a linear combination of 				*/
/* two vectors with given weights						*/
/*                                     		                                */
/********************************************************************************/
EXPORT double maximum_element(int vector_length, //length of the vector
		     Array1<double> input_vector  //input vector 
	 )
{  
      int component_index; //index to the components of the two vectors
      double maximum_element=-10.0E10; //the norm of the vector
      for(component_index=0;component_index<vector_length;component_index++)
      {
	maximum_element=std::max(input_vector[component_index],maximum_element);
      }
      return maximum_element;
}	
	
	
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the sign of a <double>                                  */
/*  weights									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The function returns the absolute value of the first argument if the second  */
/* argument is positive, -1 * the absolute value of the first argument is the   */
/* second argument is negative and zero in all other cases                      */
/********************************************************************************/
EXPORT double sign(double value, double set_sign)
{
	    if(set_sign>0) return fabs(value);
	    if(set_sign<0) return -1.0*fabs(value);
	    return 0;
}

/********************************************************************************/
/********************************************************************************/
/*  Function to compute the divergence of a vector field                        */
/*  										*/
/*  Programmer	: Coen Hennipman					*/
/*  Date	: 20-05-2015       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The first argument is calculated and is the divergence of the vector field.	*/
/********************************************************************************/
EXPORT void divergence_of_vector_field(
    Array3<double> output_array,		// the output array that is calculated
    Array3<double> first_vector_field,		// the first input vector field 
    Array3<double> second_vector_field,		// the second input vector field 
    Array3<double> third_vector_field,		// the third input vector field 
    double mesh_width_x1,			// grid spacing in x1 direction (uniform)
    double mesh_width_x2,			// grid spacing in x2 direction (uniform)
    double mesh_width_x3,			// grid spacing in x3 direction (uniform)
    int number_primary_cells_i,	        	// number of primary (pressure) cells in x1 direction
    int number_primary_cells_j,	        	// number of primary (pressure) cells in x2 direction
    int number_primary_cells_k	        	// number of primary (pressure) cells in x3 direction
)
{
      int i_index, j_index, k_index;  			// local variables for loop indexing
      double one_over_dx1	=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      double one_over_dx2	=    			// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
      double one_over_dx3	=    			// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);
	    
      for( i_index=1; i_index< number_primary_cells_i+1;i_index++)
      {
	  for( j_index=1; j_index< number_primary_cells_j+1 ; j_index++)
	  {
	    for( k_index=1; k_index< number_primary_cells_k+1; k_index++)
	    {
	      output_array[i_index][j_index][k_index]=
	      ( first_vector_field[i_index][j_index][k_index]-first_vector_field[i_index-1][j_index  ][k_index  ] )
	      *one_over_dx1+ 
	      ( second_vector_field[i_index][j_index][k_index]-second_vector_field[i_index  ][j_index-1][k_index  ] )
	      *one_over_dx2+
	      ( third_vector_field[i_index][j_index][k_index]-third_vector_field[i_index  ][j_index  ][k_index-1] )
	      *one_over_dx3;
	    }  
	  }
      }
}

/********************************************************************************/
/********************************************************************************/
/*  Function to apply a precomputed incomplete choleski preconditioner          */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The incomplete choleski preconditioner is applied. It is based on a          */
/* splitting of the original matrix and is therefore not explicitly computed.   */
/* It is sufficient to store only a diagonal matrix. The formulation of the     */
/* preconditioner is described in the Lecture Notes of CFD II Burgers course.   */
/********************************************************************************/
     void build_preconditioner(
	int i_dimension,    	// number of unknowns in the system in i-direction
	int j_dimension, 	// number of unknowns in the system in i-direction
	int k_dimension, 	// number of unknowns in the system in i-direction
	double **matrix_A,			// matrix under consideration
	double *preconditioner_matrix_M 	// preconditioner matrix
	)
      
     
     {
	int number_dof_in_slice= 		// number of degrees of freedom in one slice of the domain
		i_dimension*j_dimension; 	// as considered in the matrix
			   
	int dimension_system=			// dimension of the linear system 				
	  i_dimension*j_dimension*k_dimension;
	int row_index;	   			// index to indicate the row of the matrix 
     
     
        preconditioner_matrix_M[0]=matrix_A[0][0];
	
     
	for(row_index=1;row_index<i_dimension;row_index++){
	    preconditioner_matrix_M[row_index]=matrix_A[0][row_index]-
		matrix_A[1][row_index-1]*
		      matrix_A[1][row_index-1]/
				  preconditioner_matrix_M[row_index-1];
	}
	
	for(row_index=i_dimension;row_index<number_dof_in_slice;row_index++){
	    preconditioner_matrix_M[row_index]=matrix_A[0][row_index]-
		matrix_A[1][row_index-1]*
		      matrix_A[1][row_index-1]/
				  preconditioner_matrix_M[row_index-1]-
		matrix_A[2][row_index-i_dimension]*
		      matrix_A[2][row_index-i_dimension]/
				  preconditioner_matrix_M[row_index-i_dimension];
		      
	}

	for(row_index=number_dof_in_slice;row_index<dimension_system;row_index++){
	    preconditioner_matrix_M[row_index]=matrix_A[0][row_index]-
		matrix_A[1][row_index-1]*
		      matrix_A[1][row_index-1]/
				  preconditioner_matrix_M[row_index-1]-
		matrix_A[2][row_index-i_dimension]*
		      matrix_A[2][row_index-i_dimension]/
				  preconditioner_matrix_M[row_index-i_dimension]-
		matrix_A[3][row_index-number_dof_in_slice]*
		      matrix_A[3][row_index-number_dof_in_slice]/
				  preconditioner_matrix_M[row_index-number_dof_in_slice];
		      
	}
	
     }


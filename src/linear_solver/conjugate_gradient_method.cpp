#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to iteratively solve a linears system using the conjugate gradient */
/*  method. 									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The conjugate gradient method is a linear solver for symmetric positive      */
/* definite systems. In this solver an incomplete choleski preconditioner is    */
/* used.                                       	                         	*/
/********************************************************************************/
      int conjugate_gradient_method(
      int i_dimension,   			// number of unknowns in the system in i-direction
      int j_dimension,   			// number of unknowns in the system in i-direction
      int k_dimension,   			// number of unknowns in the system in i-direction
      Array2<double> matrix_A,   			// matrix of the linear system 
      Array1<double> preconditioner_matrix_M,  	// preconditioner matrix
      Array1<double> rhside_vector_b,   		// right hand side vector
      Array1<double> solution_vector_x, 		// solution vector
      double   tolerance,	  		// the tolerance with which the system is solved	
      int   &iteration_number,  		// the number of iterations where the iterative
						// process was terminated
      double  &relative_L2_norm_residual, 	// the L2 norm of the residual, scaled with
						// the L2 norm of the right hand side
      double  &relative_Linfinity_norm_residual,// the L infinity norm of the residual
						// scaled with the maximum difference 
						// between two components of the residual
      int maximum_iterations_allowed	 	// maximum number of iterations allowed for the
						// conjugate gradient method
     )
      {    
	  
/***************************************************************************************/	
      /* function definitions */	  
 /***************************************************************************************/	
      void copy_vector( int length_vector, 		// length of both vectors
	      Array1<double> original_vector, 			// first input vector
	      Array1<double> image_vector			// second input vector
		      );	  
      double compute_vector_norm(int vector_length, 	//length of the vector
		     Array1<double> input_vector  		//input vector 
	 );   
      double dot_product( 
	      int vector_length, 			// length of both input vectors
	      Array1<double> x, 				// first input vector
	      Array1<double> y					// second input vector
	  );
      double maximum_element(int vector_length, 	//length of the vector
		     Array1<double> input_vector  		//input vector 
	  );
      double minimum_element(int vector_length, 	//length of the vector
		     Array1<double> input_vector  		//input vector 
	  );
      void matrix_vector_product(
	    int i_dimension,   // number of unknowns in the system in i-direction
	    int j_dimension,   // number of unknowns in the system in i-direction
	    int k_dimension,   // number of unknowns in the system in i-direction
	    Array2<double> A,    // matrix under consideration
	    Array1<double> x,     // INPUT vector x
	    Array1<double> y      // OUTPUT vector y such that y=Ax
	  );
      void linear_combination(
	int vector_length, 		//length of all vectors 
	Array1<double> input_vector_x, 	//input vector 1, named x
	Array1<double> input_vector_y,		//input vector 2, named y
	Array1<double> output_vector_z,	//output vector, named z
					// such that z=x+alpha*y
	double weight_of_y		//weight alpha in the linear combination above
	);
      void apply_preconditioner(
	int i_dimension,    	// number of unknowns in the system in i-direction
	int j_dimension, 	// number of unknowns in the system in i-direction
	int k_dimension, 	// number of unknowns in the system in i-direction
	Array2<double> matrix_A,			// matrix under consideration
	Array1<double> preconditioner_matrix_M, 	// preconditioner matrix
	Array1<double> residual_vector, 		// residual vector b-Ax
	Array1<double> vector_z			// vector to store the preconditioned residual	
	);			
/***************************************************************************************/	

      double norm_rhside_vector;			// L2 norm of the right hand side vector
      double dot_product_residual_vector_vector_z;	// r_n*z_n
      double dot_product_p_A_p;				// p'*A*p
      double norm_residual_vector;			// ||r||
      double maximum_element_solution_vector;		// max element of solution vector x
      double minimum_element_solution_vector;		// min element of solution vector x
      double maximum_element_residual_vector;		// max element of residual vector r
      double dot_product_residual_vector_vector_z_old;	// r_(n-1)*z_(n-1)
      Array1<double> residual_vector;				// residual vector b-Ax
      Array1<double> direction_vector_p;			// direction vector for the update p	
      Array1<double> store_matrix_vector_product;  		// help vector to store product of A times 
      Array1<double> vector_z;			   		// vector to store the preconditioned residual					      
      double alpha;					// weight for the update of the solution 
      int dimension_system=				// dimension of the linear system 				
	  i_dimension*j_dimension*k_dimension;
      int convergence=100;				//=0 when convergence is achieved
							//=1 when algorithm failed to converge
	
      int i;
      /* allocate memory for the local vectors needed in the algorithm */
      
      residual_vector.create(dimension_system);			
      direction_vector_p.create(dimension_system);			
      store_matrix_vector_product.create(dimension_system);  	
      vector_z.create(dimension_system);   	
      
      for(i=0;i<dimension_system;i++){
	solution_vector_x[i]=0.0;
      }
      
     
      /* compute A*x */
      matrix_vector_product(i_dimension, j_dimension, k_dimension, matrix_A, 
			      solution_vector_x, store_matrix_vector_product);
      /* compute r=b-A*x */
      linear_combination(dimension_system,rhside_vector_b, store_matrix_vector_product,
			     residual_vector,-1.0);
      
      /* solve Mz=r */
      apply_preconditioner(i_dimension, j_dimension, k_dimension, matrix_A,
			      preconditioner_matrix_M, residual_vector, vector_z);
      
      /* copy z to direction vector p */
      copy_vector(dimension_system, vector_z, direction_vector_p);
      
      /* compute the norm of the right hand side vector */
      norm_rhside_vector=compute_vector_norm(dimension_system,rhside_vector_b);
	     
      /* check if the norm of the right hand side vector is not that small that   */
      /* the homogeneous solution is sufficient to meet the convergence criterium */
      /* because in that case the solution process is already completed here      */
      
      if(norm_rhside_vector<tolerance)
      {
	 /* in this case these norms are norms of the zero vector and we are done */
	 
	 relative_L2_norm_residual=0.0;
	 relative_Linfinity_norm_residual=0.0;
	 iteration_number= 0;
	 convergence=0;
	 std::cerr<<"Convergence for trivial solution of homogeneous system. \n";
      }
      else
      {
	 /* in this case more work is needed and the iteration started */
	 /* compute r_0*z_0 */
	dot_product_residual_vector_vector_z=
			dot_product(dimension_system, residual_vector, vector_z);
			
	for(iteration_number=1;iteration_number<maximum_iterations_allowed;iteration_number++)
	{
	  

	  /* compute A*p  */
	  matrix_vector_product(i_dimension, j_dimension, k_dimension, matrix_A, 
			      direction_vector_p, store_matrix_vector_product);
	  /* compute p * (A*p) */
	  dot_product_p_A_p=dot_product(dimension_system,
					  direction_vector_p, store_matrix_vector_product);
	  alpha=dot_product_residual_vector_vector_z/dot_product_p_A_p;
	  
	  /* update solution vector x */
	    linear_combination(dimension_system,solution_vector_x, 
				    direction_vector_p,
			     solution_vector_x, alpha);	
	  
	  /* update residual vector r */
	    linear_combination(dimension_system, residual_vector,
				    store_matrix_vector_product, 
						  residual_vector, -1.0*alpha);
	  /* compute the norm of the residual vector */
	    norm_residual_vector=compute_vector_norm(dimension_system,residual_vector);
	    relative_L2_norm_residual=norm_residual_vector/norm_rhside_vector;
	    
	  /* identify  the largest element of the solution vector */

	    maximum_element_solution_vector=maximum_element(dimension_system,solution_vector_x);

	  /* identify  the smallest element of the solution vector */

	    minimum_element_solution_vector=minimum_element(dimension_system,solution_vector_x);

	   /* identify  the largest element of the residual vector */

	    maximum_element_residual_vector=maximum_element(dimension_system,residual_vector);

	   /* determine the relative L_infinity norm of the residual vector */
	   
	    relative_Linfinity_norm_residual=maximum_element_residual_vector/
	               (maximum_element_solution_vector-minimum_element_solution_vector);
		       
	  /* use both norms to assess convergence of the solution */
	  
// 	    if( (maximum_element_residual_vector<tolerance*
// 	      (maximum_element_solution_vector-minimum_element_solution_vector)) &&
// 		  relative_L2_norm_residual<tolerance) 
// 	    {
	    if(relative_L2_norm_residual<tolerance) 
	    {
	      /* the solution is converged */
	      convergence=0;
	      break;
	    }
	    else
	      
	      if(iteration_number<maximum_iterations_allowed)  
	      {
		  /* continue the iteration process */ 
		  /* apply the preconditioner */
		  
		    /* solve Mz=r */
		  apply_preconditioner(i_dimension, j_dimension, k_dimension, matrix_A,
			      preconditioner_matrix_M, residual_vector, vector_z);
		   /* save r_{n-1}*z_{n-1} */

		  dot_product_residual_vector_vector_z_old=
				    dot_product_residual_vector_vector_z;
		   /* compute r_{n}*z_{n} */
		  dot_product_residual_vector_vector_z=
			dot_product(dimension_system, residual_vector, vector_z);
		  /* compute the new direction vector p*/
		  linear_combination(dimension_system, vector_z,
				    direction_vector_p, direction_vector_p,
		  dot_product_residual_vector_vector_z/dot_product_residual_vector_vector_z_old);
		  
	      }
	      else
	      {
	      /* no convergence could be achieved */
	       	  convergence=1;  
		  break;
	      }
	}	
	
      }
	      /* deallocate memory for the local vectors needed in the algorithm */
      
      residual_vector.destroy();
      direction_vector_p.destroy();
      store_matrix_vector_product.destroy();
      vector_z.destroy();

      return convergence;
      
}      




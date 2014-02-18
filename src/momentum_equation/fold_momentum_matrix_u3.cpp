#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to include the boundary conditions into the momentum equation      */
/*  conditions									*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* In the current implementation the generation of the matrix, right-hand-side  */
/* and the application of the boundary conditions is completely separated.      */
/* In the current function the matrix is build WITHOUT considering the          */
/* boundary conditions. Next the boundary conditions are considered in a        */
/* separate 'matrix folding' function that eliminates all known values  and     */
/* discards any reference to virtual values.                                    */
/* The tasks are:								*/
/* 1) Determine the domain of unknowns						*/
/* 2) remove any connections outside the domain	using the boundary conditions	*/
/* 
/********************************************************************************/
EXPORT void fold_momentum_matrix_u3(
      boundary_face boundary_faces[6],			// array with all the information
							// for the boundary conditions 
      Array2<double> momentum_matrix_u3,			// momentum matrix velocity x2 direction
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      int number_matrix_connections			// number of connections in momentum matrix				       
	   )
      
{
    int start_index_i;			// lowest index i_index in active part of matrix
    int start_index_j;			// lowest index j_index in active part of matrix
    int start_index_k;			// lowest index k_index in active part of matrix
    int final_index_i;			// highest index i_index in active part of matrix
    int final_index_j;			// highest index j_index in active part of matrix
    int final_index_k;			// highest index k_index in active part of matrix
    int cell_label_boundary;		// index of the matrix element to be manipulated
    int face_index;			// index of the face to be manipulated
    int i_index, j_index, k_index;  	// local variables for loop indexing
    int connection_index;		// local variable for connection indexing
    int one_dimensional_index;		// index of point in 1-D array

    
   /* set default values for start and final index of the unknowns */
   
    start_index_i=1;
    final_index_i=number_primary_cells_i;
    start_index_j=1;
    final_index_j=number_primary_cells_j;
    start_index_k=0;
    final_index_k=number_primary_cells_k;

     
   /* start the loop over all faces for the matrix folding */
   
   
       /******************************************************************/
       /*   +/- I-index faces						 */
       /******************************************************************/
       for(face_index=0;face_index<=1; face_index++)
       {
	  if(face_index==0)
	  {
	      cell_label_boundary=final_index_i;
	  }
	  else
	  {
	      cell_label_boundary=start_index_i;
	  }
	  
	  if(boundary_faces[face_index].boundary_variables[2].boundary_condition_type==dirichlet)
	  {
	      
		  /* DIRICHLET BOUNDARY CONDITION */
			
		if(face_index==0)
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using linear interpolation */
			/* U_virtual=2*U_boundary - U_real */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with opposite sign */
			/* and the nonexistent connection is set to zero */

			    one_dimensional_index=map_index_u3(cell_label_boundary,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]-=
			    momentum_matrix_u3[1][one_dimensional_index];
			    momentum_matrix_u3[1][one_dimensional_index]=0;
			
			}	  
  
		    }  
		}  
		else
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using linear interpolation */
			/* U_virtual=2*U_boundary - U_real */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with opposite sign */

			    one_dimensional_index=map_index_u3(cell_label_boundary,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]-=
			    momentum_matrix_u3[4][one_dimensional_index];
			
			}	  
  
		    }  
		}
	  }
	  else
	  {
	      if(boundary_faces[face_index].boundary_variables[2].boundary_condition_type==neumann)
	      {
		  /* NEUMANN BOUNDARY CONDITION */
			
		if(face_index==0)
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using second order approximation */
			/* to the normal derivative */
			/* U_virtual=U_real + dx* dU/dn */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with equal sign */
			/* and the nonexistent connection is set to zero */

			    one_dimensional_index=map_index_u3(cell_label_boundary,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]+=
			    momentum_matrix_u3[1][one_dimensional_index];
			    momentum_matrix_u3[1][one_dimensional_index]=0;
			
			}	  
  
		    }  
		}  
		else
		{
		    for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using second order approximation */
			/* to the normal derivative */
			/* U_virtual=U_real + dx* dU/dn */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with equal sign */
			/* and the nonexistent connection is set to zero */

			    one_dimensional_index=map_index_u3(cell_label_boundary,j_index,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]+=
			    momentum_matrix_u3[4][one_dimensional_index];
			
			}	  
  
		    }  
		}
			
	      }  
	      else
	      {
	    
		  std::cerr<< "**************************************************** \n";
		  std::cerr<<"ERROR \n";
		  std::cerr<<"you are trying to impose a boundary condition that \n";
		  std::cerr<<"is not implemented yet, check your input";
		  std::cerr<<"in function fold_momentum_matrix_u3 line 232 \n";
		  std::cerr<<"**************************************************** \n";
		  exit(1);
	      }
	  }
  
       }  

       
       /******************************************************************/
       /*   +/- J-index faces						 */
       /******************************************************************/

 
       for(face_index=2;face_index<=3; face_index++)
       {
	  if(face_index==2)
	  {
	      cell_label_boundary=final_index_j;
	  }
	  else
	  {
	      cell_label_boundary=start_index_j;
	  }
	  
	  if(boundary_faces[face_index].boundary_variables[2].boundary_condition_type==dirichlet)
	  {
	      
		  /* DIRICHLET BOUNDARY CONDITION */
			
		if(face_index==2){
		    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using linear interpolation */
			/* U_virtual=2*U_boundary - U_real */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with opposite sign */
			/* and the nonexistent connection is set to zero */

			    one_dimensional_index=map_index_u3(i_index,cell_label_boundary,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]-=
			    momentum_matrix_u3[2][one_dimensional_index];
			    momentum_matrix_u3[2][one_dimensional_index]=0;
			
			}	  
  
		    }  
		}  
		else
		{
		    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using linear interpolation */
			/* U_virtual=2*U_boundary - U_real */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with opposite sign */

			    one_dimensional_index=map_index_u3(i_index,cell_label_boundary,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]-=
			    momentum_matrix_u3[5][one_dimensional_index];
			
			}	  
  
		    }  
		}
	  }
	  else
	  {
	      if(boundary_faces[face_index].boundary_variables[2].boundary_condition_type==neumann)
	      {
		  /* NEUMANN BOUNDARY CONDITION */

		  if(face_index==2){
		    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using second order approximation */
			/* to the normal derivative */
			/* U_virtual=U_real + dx* dU/dn */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with equal sign */
			/* and the nonexistent connection is set to zero */

			    one_dimensional_index=map_index_u3(i_index,cell_label_boundary,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]+=
			    momentum_matrix_u3[2][one_dimensional_index];
			    momentum_matrix_u3[2][one_dimensional_index]=0;
			
			}	  
  
		    }  
		}  
		else
		{
		    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		    {
			for(k_index=1;k_index<number_primary_cells_k;k_index++)
			{
			
			/* the boundary condition is applied using second order approximation */
			/* to the normal derivative */
			/* U_virtual=U_real + dx* dU/dn */
			/* the contribution of the virtual cell is moved to the central */
			/* coefficient with equal sign */
			    one_dimensional_index=map_index_u3(i_index,cell_label_boundary,k_index,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);


			    momentum_matrix_u3[0][one_dimensional_index]+=
			    momentum_matrix_u3[5][one_dimensional_index];
			
			}	  
  
		    }  
		}
			
	      }  
	      else
	      {
	    
		  std::cerr<< "**************************************************** \n";
		  std::cerr<<"ERROR \n";
		  std::cerr<<"you are trying to impose a boundary condition that \n";
		  std::cerr<<"is not implemented yet, check your input";
		  std::cerr<<"in function fold_momentum_matrix_u3 line 358 \n";
		  std::cerr<<"**************************************************** \n";
		  exit(1);
	      }
	  }
  
       }  

      
       /******************************************************************/
       /*   +/- K-index faces						 */
       /******************************************************************/
      
       for(face_index=4;face_index<=5; face_index++)
       {
	  if(face_index==4)
	  {
	      cell_label_boundary=final_index_k;
	  }
	  else
	  {
	      cell_label_boundary=start_index_k;
	  }
	  if(boundary_faces[face_index].boundary_variables[2].boundary_condition_type==dirichlet)
	  {
	      if(face_index==4)
	      {
		    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		    {
			for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
			{
			/* U_real=U_boundary */

			    one_dimensional_index=map_index_u3(i_index,j_index,cell_label_boundary,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]=1.0;
			    for(connection_index=1;connection_index<number_matrix_connections;
							      connection_index++)
			    {
				momentum_matrix_u3[connection_index][one_dimensional_index]=0.0;
			    }
			}	  
  
		    }  
		    
		 /* because only the upper triangular part is stored       */
		 /* the reference to these 'known' values is removed 	   */
		 /* only for the face 4, the other non existing references */
		 /* are in the lower triangular part of the matrix	   */
		 
		    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		    {
			for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
			{
			/* U_real=U_boundary */

			    one_dimensional_index=map_index_u3(i_index,j_index,cell_label_boundary-1,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[3][one_dimensional_index]=0.0;
			}	  
  
		    }  
	      }
	      else
	      {
		 
		    for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
		    {
			for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
			{
			/* U_real=U_boundary */

			    one_dimensional_index=map_index_u3(i_index,j_index,cell_label_boundary,
				      number_primary_cells_i, number_primary_cells_j, 
								  number_primary_cells_k);

			    momentum_matrix_u3[0][one_dimensional_index]=1.0;
			    for(connection_index=1;connection_index<number_matrix_connections;
							      connection_index++)
			    {
				momentum_matrix_u3[connection_index][one_dimensional_index]=0.0;
			    }
			}	  
  
		    }  
	      }
	  }  
	  else
	  {
	      std::cerr<< "**************************************************** \n";
	      std::cerr<<"ERROR \n";
	      std::cerr<<"you are trying to impose a boundary condition that \n";
	      std::cerr<<"is not implemented yet, check your input";
	      std::cerr<<"in function fold_momentum_matrix_u2 line 232 \n";
	      std::cerr<<"**************************************************** \n";
	      exit(1);
	  }
  
       }  

       
  
}

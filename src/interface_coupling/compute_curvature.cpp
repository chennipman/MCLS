#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the interface curvature                                 */
/*                                                                 		*/
/*  Programmer	: Duncan van der Heul       						*/
/*  Date	: 10-03-2013       							*/
/*  Update	:        								*/
/********************************************************************************/
/* Notes										*/
/*  The interface curvature is computed using a second order finite difference  */
/*  approximation of the derivatives of the level-set field.                    */
/*  The curvature is computed in all real cells and extrapolated to the         */
/*  virtual cells.									*/
/********************************************************************************/
EXPORT void compute_curvature(						
      Array3<double> level_set, 				// level-set field
      Array3<double> curvature,				// interface curvature
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3				// grid spacing in x3 direction (uniform)
	)
      {
      Array3<double> d_level_set_d_x1;			// first partial derivative wrt x1 of level-set
      Array3<double> d_level_set_d_x2;			// first partial derivative wrt x2 of level-set
      Array3<double> d_level_set_d_x3;			// first partial derivative wrt x3 of level-set
      Array3<double> d_2_level_set_d_x1_2;		// pure second partial derivative wrt x1 of level-set
      Array3<double> d_2_level_set_d_x2_2;		// pure second partial derivative wrt x2 of level-set
      Array3<double> d_2_level_set_d_x3_2;		// pure second partial derivative wrt x3 of level-set
    
      Array3<double> d_2_level_set_d_x1_d_x2;		// mixed second partial derivative wrt x1 and x2 of level-set
      Array3<double> d_2_level_set_d_x1_d_x3;		// mixed second partial derivative wrt x1 and x3 of level-set
      Array3<double> d_2_level_set_d_x2_d_x3;		// mixed second partial derivative wrt x2 and x3 of level-set
      Array3<double> length_gradient;			// length of the level-set gradient vector
      Array3<double> curvature_error;			// error in the curvature for laplace case

      double one_over_dx1	=    			// 1/(grid spacing in x1 direction)
	    1.0/(mesh_width_x1);
      double one_over_dx2=    			// 1/(grid spacing in x2 direction)
	    1.0/(mesh_width_x2);
      double one_over_dx3=    			// 1/(grid spacing in x3 direction)
	    1.0/(mesh_width_x3);	
      double one_over_dx1_squared=    		// 1/(grid spacing in x1 direction)^2
	    1.0/(mesh_width_x1*mesh_width_x1);
      double one_over_dx2_squared=    		// 1/(grid spacing in x2 direction)^2
	    1.0/(mesh_width_x2*mesh_width_x2);
      double one_over_dx3_squared=    		// 1/(grid spacing in x3 direction)^2
	    1.0/(mesh_width_x3*mesh_width_x3);		
      double curvature_threshold;			// computed curvature should not exceed 
							// the threshold on the curvature
      double domain_center_x1;			// location of bubble center x1 coordinate
      double domain_center_x2;			// location of bubble center x2 coordinate
      double domain_center_x3;			// location of bubble center x3 coordinate
      double exact_curvature;				// exact curvature for the bubble
      double cell_center_x1;				// cell center x1 coordinate
      double cell_center_x2;				// cell center x2 coordinate
      double cell_center_x3;				// cell center x3 coordinate
      double distance_cell_center_domain_center;	// distance from cell center to center of the domain
      int i_index, j_index, k_index;  		// local variables for loop indexing
      

	/* allocate  memory for the derivatives of the level-set field */
	/* and for the weighting function that is applied in the CSF model */
	
	d_level_set_d_x1.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x2.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_level_set_d_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_2_level_set_d_x1_2.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_2_level_set_d_x2_2.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_2_level_set_d_x3_2.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_2_level_set_d_x1_d_x2.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_2_level_set_d_x1_d_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	d_2_level_set_d_x2_d_x3.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	length_gradient.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	curvature_error.create(number_primary_cells_i+2, number_primary_cells_j+2,
				    number_primary_cells_k+2);
	
	
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_level_set_d_x1 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_level_set_d_x2 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_level_set_d_x3 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_2_level_set_d_x1_2 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_2_level_set_d_x2_2 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_2_level_set_d_x3_2 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_2_level_set_d_x1_d_x2 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_2_level_set_d_x1_d_x3 , 0.0);
       set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2,d_2_level_set_d_x2_d_x3 , 0.0);
     
      /* the curvature that is computed from the level-set field should not exceed the threshold */
      /* values -curvature_threshold and +curvature_threshold */
      
      /* curvature_threshold= min(1/dx1, 1/dx2, 1/dx3) + max(1/dx1,1/dx2,1/dx3) */
      
      curvature_threshold= std::min(one_over_dx1, std::min(one_over_dx2, one_over_dx3))+
		    std::max(one_over_dx1, std::max(one_over_dx2, one_over_dx3));



	    /* central differences to approximate the derivatives up to second order */
	    
      domain_center_x1=number_primary_cells_i*mesh_width_x1*0.5;
      domain_center_x2=number_primary_cells_j*mesh_width_x2*0.5;
      domain_center_x3=number_primary_cells_k*mesh_width_x3*0.5;

      for(i_index=1;i_index<number_primary_cells_i+1;i_index++)
      {
	  for(j_index=1;j_index<number_primary_cells_j+1;j_index++)
	  {
	      for(k_index=1;k_index<number_primary_cells_k+1;k_index++)
	      {
	    
	    /* first derivatives */
	    
		d_level_set_d_x1[i_index][j_index][k_index]=one_over_dx1*0.5*(level_set[i_index+1][j_index][k_index]-
						  level_set[i_index-1][j_index][k_index]);
		d_level_set_d_x2[i_index][j_index][k_index]=one_over_dx2*0.5*(level_set[i_index][j_index+1][k_index]-
						  level_set[i_index][j_index-1][k_index]);
		d_level_set_d_x3[i_index][j_index][k_index]=one_over_dx3*0.5*(level_set[i_index][j_index][k_index+1]-
						  level_set[i_index][j_index][k_index-1]);
	    
	    /* pure second derivatives */
	    
		d_2_level_set_d_x1_2[i_index][j_index][k_index]=one_over_dx1_squared*(
				      level_set[i_index+1][j_index][k_index]
					-2.0*level_set[i_index][j_index][k_index]
						  +level_set[i_index-1][j_index][k_index]);
		d_2_level_set_d_x2_2[i_index][j_index][k_index]=one_over_dx2_squared*(
				      level_set[i_index][j_index+1][k_index]
					-2.0*level_set[i_index][j_index][k_index]
						  +level_set[i_index][j_index-1][k_index]);
		d_2_level_set_d_x3_2[i_index][j_index][k_index]=one_over_dx3_squared*(
				      level_set[i_index][j_index][k_index+1]
					-2.0*level_set[i_index][j_index][k_index]
						  +level_set[i_index][j_index][k_index-1]);
	    
	    /* mixed second derivatives */
	    
		d_2_level_set_d_x1_d_x2[i_index][j_index][k_index]=0.25*
				      one_over_dx1*one_over_dx2*
				    ( level_set[i_index+1][j_index+1][k_index]-
					 level_set[i_index+1][j_index-1][k_index] -
					    level_set[i_index-1][j_index+1][k_index]+
						level_set[i_index-1][j_index-1][k_index] );
	    
		d_2_level_set_d_x1_d_x3[i_index][j_index][k_index]=0.25*
				      one_over_dx1*one_over_dx3*
				    ( level_set[i_index+1][j_index][k_index+1]-
					 level_set[i_index+1][j_index][k_index-1] -
					    level_set[i_index-1][j_index][k_index+1]+
						level_set[i_index-1][j_index][k_index-1] );
	    
		d_2_level_set_d_x2_d_x3[i_index][j_index][k_index]=0.25*
				      one_over_dx2*one_over_dx3*
				    ( level_set[i_index][j_index+1][k_index+1]-
					 level_set[i_index][j_index+1][k_index-1] -
					    level_set[i_index][j_index-1][k_index+1]+
						level_set[i_index][j_index-1][k_index-1] );
             
	    /* length of the gradient of the level-set vector */
	    /* this should be very close to one of course */
	    
	    
		length_gradient[i_index][j_index][k_index]=sqrt(  
 			d_level_set_d_x1[i_index][j_index][k_index]*d_level_set_d_x1[i_index][j_index][k_index] + 
			  	d_level_set_d_x2[i_index][j_index][k_index]*d_level_set_d_x2[i_index][j_index][k_index]+
		        		d_level_set_d_x3[i_index][j_index][k_index]*d_level_set_d_x3[i_index][j_index][k_index] );

	     /* check that the length of the gradient does not vanish */
	     
	     	if(length_gradient[i_index][j_index][k_index]<1E-12)
		{
			
			length_gradient[i_index][j_index][k_index]=1E-12;
		}
	    
	    /* next the curvature formula, as presented in Sanders thesis on page ** */
	    
	    	cell_center_x1=(i_index-0.5)*mesh_width_x1;
	    	cell_center_x2=(j_index-0.5)*mesh_width_x2;
	    	cell_center_x3=(k_index-0.5)*mesh_width_x3;
	    	distance_cell_center_domain_center=sqrt(
				(domain_center_x1-cell_center_x1)*(domain_center_x1-cell_center_x1)+
				(domain_center_x2-cell_center_x2)*(domain_center_x2-cell_center_x2)+
				(domain_center_x3-cell_center_x3)*(domain_center_x3-cell_center_x3));
		
	    	exact_curvature=2.0/distance_cell_center_domain_center;
	    
		curvature[i_index][j_index][k_index]=
			  (d_2_level_set_d_x1_2[i_index][j_index][k_index]+d_2_level_set_d_x2_2[i_index][j_index][k_index]+d_2_level_set_d_x3_2[i_index][j_index][k_index])/
				  length_gradient[i_index][j_index][k_index]-
			  (d_level_set_d_x1[i_index][j_index][k_index]*(d_level_set_d_x1[i_index][j_index][k_index]*d_2_level_set_d_x1_2[i_index][j_index][k_index]+
					      d_level_set_d_x2[i_index][j_index][k_index]*d_2_level_set_d_x1_d_x2[i_index][j_index][k_index]+
						d_level_set_d_x3[i_index][j_index][k_index]*d_2_level_set_d_x1_d_x3[i_index][j_index][k_index])+
			   d_level_set_d_x2[i_index][j_index][k_index]*(d_level_set_d_x1[i_index][j_index][k_index]*d_2_level_set_d_x1_d_x2[i_index][j_index][k_index]+
					      d_level_set_d_x2[i_index][j_index][k_index]*d_2_level_set_d_x2_2[i_index][j_index][k_index]+
						d_level_set_d_x3[i_index][j_index][k_index]*d_2_level_set_d_x2_d_x3[i_index][j_index][k_index])+
			   d_level_set_d_x3[i_index][j_index][k_index]*(d_level_set_d_x1[i_index][j_index][k_index]*d_2_level_set_d_x1_d_x3[i_index][j_index][k_index]+
					      d_level_set_d_x2[i_index][j_index][k_index]*d_2_level_set_d_x2_d_x3[i_index][j_index][k_index]+
						d_level_set_d_x3[i_index][j_index][k_index]*d_2_level_set_d_x3_2[i_index][j_index][k_index]))/
					      (length_gradient[i_index][j_index][k_index]*length_gradient[i_index][j_index][k_index]*length_gradient[i_index][j_index][k_index]);
		curvature_error[i_index][j_index][k_index]=exact_curvature-curvature[i_index][j_index][k_index];
			
				
	    /* check that the curvature is within bounds */
	    
		if(curvature[i_index][j_index][k_index]> curvature_threshold) 
		{
		    curvature[i_index][j_index][k_index]=curvature_threshold;
		}
		if(curvature[i_index][j_index][k_index]< -1.0* curvature_threshold) 
		{
		    curvature[i_index][j_index][k_index]= -1.0*curvature_threshold;
		}
	    
	      }
	  }
      }
      
      /* impose a homogeneous neumann boundary condition on the curvature */
      /* this means the curvature is extrapolated to the virtual cells    */
      
	  field_neumann_boundary(curvature, 
				  number_primary_cells_i, number_primary_cells_j,	
				    number_primary_cells_k);
	  
	  /*write all quantities that are involved in the computation of the curvature to a vtk file */
/*	  
	 dump_curvature_for_debugging(d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,		
       						d_2_level_set_d_x1_2,	d_2_level_set_d_x2_2,	d_2_level_set_d_x3_2,	
       							d_2_level_set_d_x1_d_x2, d_2_level_set_d_x1_d_x3, d_2_level_set_d_x2_d_x3,	
 								length_gradient, curvature, curvature_error,			
							number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		
	  						mesh_width_x1,	mesh_width_x2, mesh_width_x3	);*/
	  
	d_level_set_d_x1.destroy();
	d_level_set_d_x2.destroy();
	d_level_set_d_x3.destroy();
	d_2_level_set_d_x1_2.destroy();
	d_2_level_set_d_x2_2.destroy();
	d_2_level_set_d_x3_2.destroy();
	d_2_level_set_d_x1_d_x2.destroy();
	d_2_level_set_d_x1_d_x3.destroy();
       d_2_level_set_d_x2_d_x3.destroy();
	length_gradient.destroy();
	curvature_error.destroy();
}

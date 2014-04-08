#include "../src/headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to test the momentum_predictor					*/
/*                     								*/
/*  										*/
/*  Programmer	: Coen Hennipman	       					*/
/*  Date	: 13-03-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* Unit test for the momentum predictor			                        */
/* Both convection and diffusion can be tested		                        */
/* The ouput polynome should be equal to the output from the 		        */
/* momentum_predictor for linear cases. For higher higher order it is only an   */
/* estimation and the error should go to zero for small enough grid size.	*/
/* The convection and diffusion can be tested seperataly. 			*/
/********************************************************************************/

int main()
{
      Array3<double> pressure;				// pressure field
      Array3<double> level_set;				// level-set field
      Array3<double> u_1_velocity_old;			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old; 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old;			// velocity field at old time level x3 direction
      Array3<double> u_1_momentum; 			// momentum in x1 direction
      Array3<double> u_2_momentum; 			// momentum in x2 direction
      Array3<double> u_3_momentum;			// momentum in x3 direction
      Array3<double> scaled_density_u1;			// scaled density for the controlvolumes
      Array3<double> scaled_density_u2;			// scaled density for the controlvolumes
      Array3<double> scaled_density_u3;			// scaled density for the controlvolumes
      
      Array3<double> momentum_source_term_u_1;		// source term of the momentum equation in x1 direction
					       		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2;		// source term of the momentum equation in x2 direction
					        	// defined on all u2 points (including boundaries)
      Array3<double> momentum_source_term_u_3;		// source term of the momentum equation in x3 direction
					        	// defined on all u3 points (including boundaries)
					        	    
      Array3<double> surface_tension_body_force_x1;	// x1 component of the body force due to
      							// CSF formulation of surface tension model
      Array3<double> surface_tension_body_force_x2;	// x2 component of the body force due to
      							// CSF formulation of surface tension model
      Array3<double> surface_tension_body_force_x3;	// x3 component of the body force due to
      							// CSF formulation of surface tension model

      int number_primary_cells_i;			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j;			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k;			// number of primary (pressure) cells in x3 direction
  
      double mesh_width_x1;				// grid spacing in x1 direction (uniform)
      double mesh_width_x2;				// grid spacing in x2 direction (uniform)
      double mesh_width_x3;				// grid spacing in x3 direction (uniform)
      double smoothing_distance_factor;			// the smoothing distance is smoothing_distance_factor 
      
      double rho_plus_over_rho_minus;			// the density ratio      
      double rho_minus_over_mu_minus;			// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus;			// ratio of the viscosities of both phases
      
      vector initial_velocity; 
      vector gravity;
      int source_terms_in_momentum_predictor = 0;    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation  
      int i,j,k; 	// loop indexing
      double x,y,z;	// coordinates
      
      // size of the test problem
      double x_length = 1; 
      double y_length = 1;
      double z_length = 1;
      
      double a,b,c,d,e,f,g,h,l,m,n,o,p,q,r,s,t,u; 	// input parameters
      
      double u_1, u_2, u_3; 
      double viscosity_over_density;
      double output_polynomial; 
      double scaled_density_constant = 1; 
      
      number_primary_cells_i = pow(2,2); 
      number_primary_cells_j = pow(2,2); 
      number_primary_cells_k = pow(2,2); 

	mesh_width_x1 = x_length/number_primary_cells_i; 
	mesh_width_x2 = y_length/number_primary_cells_j; 
	mesh_width_x3 = z_length/number_primary_cells_k; 
      
      smoothing_distance_factor = 2; 
      rho_plus_over_rho_minus = 1.0;
      rho_minus_over_mu_minus = 1e10;
      mu_plus_over_mu_minus = 1.0; 
      
      initial_velocity.u1 = 00.0;
      initial_velocity.u2 = 0.0;
      initial_velocity.u3 = 0.0;
      
      gravity.u1 = 00.0;
      gravity.u2 = 00.0;
      gravity.u3 = 00.0;

	// input parameters
//	a=1;	b=2;	c=3;	d=4;	e=5;	f=6;	g=7;	h=8;	l=9;	m=1;	n=2;	o=3;	p=4; 	q=5; 	r=6;	s=7;	t=8;	u=9; // full test
	a=0;	b=2;	c=0;	d=4;	e=0;	f=6;	g=0;	h=8;	l=0;	m=1;	n=0;	o=3;	p=0; 	q=5; 	r=0;	s=7;	t=0;	u=9; // linear test
//	a=0;	b=0;	c=1;	d=0;	e=0;	f=0;	g=1;	h=0;	l=0;	m=0;	n=0;	o=0;	p=0; 	q=0; 	r=0;	s=0;	t=0;	u=0; // random test
//	a=0;	b=0;	c=0;	d=0;	e=0;	f=0;	g=0;	h=0;	l=0;	m=0;	n=0;	o=0;	p=0; 	q=0; 	r=0;	s=0;	t=1;	u=0; // random test

	            
      // initialize
      u_1_velocity_old.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_old.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_old.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      momentum_source_term_u_1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      momentum_source_term_u_2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      momentum_source_term_u_3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      u_1_momentum.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_momentum.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_momentum.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      scaled_density_u1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      scaled_density_u2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      scaled_density_u3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      surface_tension_body_force_x1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      surface_tension_body_force_x2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      surface_tension_body_force_x3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      pressure.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      level_set.create(number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+2);
 




  for(i=0;i<number_primary_cells_i+2;i++)
  {
     for(j=0;j<number_primary_cells_j+2;j++)
      {
	  for(k=0;k<number_primary_cells_k+2;k++)
	  {
	  x = x_length * i/number_primary_cells_i -0.5*mesh_width_x1;
	  y = y_length * j/number_primary_cells_j -0.5*mesh_width_x2;
	  z = z_length * k/number_primary_cells_k -0.5*mesh_width_x3;
	  
	  if(i!=number_primary_cells_i+1)
	  	{x += 0.5*mesh_width_x1;
//	  	printf("x = %f \t",x);
	  	u_1_velocity_old[i][j][k]= a*x*x+b*x 	+c*y*y+d*y	+e*z*z+f*z;  // convection
//	  	u_1_velocity_old[i][j][k]= a*x*x+b*y*y+c*z*z 	+d*x*y+e*x*z +f*y*z; // diffusion
//	  	printf("u_1 = %f \n", u_1_velocity_old[i][j][k]);
	  	x -= 0.5*mesh_width_x1;}
	  if(j!=number_primary_cells_j+1)
	  	{y += 0.5*mesh_width_x2;
	  	u_2_velocity_old[i][j][k]= g*x*x+h*x 	+l*y*y+m*y	+n*z*z+o*z; // convection
//	  	u_2_velocity_old[i][j][k]= g*x*x+h*y*y+l*z*z 	+m*x*y+n*x*z +o*y*z; // diffusion	  	
	  	y -= 0.5*mesh_width_x2;}
	  if(k!=number_primary_cells_k+1)
	  	{z += 0.5*mesh_width_x3;
	  	u_3_velocity_old[i][j][k]= p*x*x+q*x 	+r*y*y+s*y	+t*z*z+u*z; // convection
//	  	u_3_velocity_old[i][j][k]= p*x*x+q*y*y+r*z*z 	+s*x*y+t*x*z +u*y*z; // diffusion
	  	z -= 0.5*mesh_width_x3;}
	  }  
  
      }  
     
  } 
  
	
      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, u_1_momentum, 0.0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, u_2_momentum, 0.0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, u_3_momentum, 0.0);

      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, scaled_density_u1, scaled_density_constant);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, scaled_density_u2, scaled_density_constant);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, scaled_density_u3, scaled_density_constant);

      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, surface_tension_body_force_x1, 0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, surface_tension_body_force_x2, 0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, surface_tension_body_force_x3, 0);			    			    			    
			    
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, level_set, 100);
      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2, 
			    number_primary_cells_k+2, pressure, 10000.0);
      
       printf("Start momentum predictor \n");
       momentum_predictor(
       u_1_momentum, u_2_momentum, u_3_momentum,
       u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,
       scaled_density_u1, scaled_density_u2, scaled_density_u3,
       momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
       level_set, 
       number_primary_cells_i, number_primary_cells_j, number_primary_cells_k, 
       mesh_width_x1, mesh_width_x2, mesh_width_x3, smoothing_distance_factor, 
       rho_plus_over_rho_minus, rho_minus_over_mu_minus, mu_plus_over_mu_minus, 
       source_terms_in_momentum_predictor); 

viscosity_over_density = (1/scaled_density_constant)*compute_scaled_viscosity(100,mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus);

// check for u,v,w momentum 
for(i=1;i<number_primary_cells_i+1;i++)
{
  for(j=1;j<number_primary_cells_j+1;j++)
  {
      for(k=1;k<number_primary_cells_k+1;k++)
      {
	  x = x_length * i/number_primary_cells_i -0.5*mesh_width_x1;
	  y = y_length * j/number_primary_cells_j -0.5*mesh_width_x2;
	  z = z_length * k/number_primary_cells_k -0.5*mesh_width_x3;
	
	printf("\n i=%i j=%i k=%i \n", i,j,k);


	if(i!=number_primary_cells_i){
		x += 0.5*mesh_width_x1;
		u_1 = a*x*x+b*x 	+c*y*y+d*y	+e*z*z+f*z;
		u_2 = g*x*x+h*x 	+l*y*y+m*y	+n*z*z+o*z;
		u_3 = p*x*x+q*x 	+r*y*y+s*y	+t*z*z+u*z;
		output_polynomial = -(2*u_1*(2*a*x+b) + u_1*(2*y*l+m) + u_2*(2*y*c+d) + u_1*(2*z*t+u) + u_3*(2*z*e+f))	+viscosity_over_density*(2*2*a+2*c+2*e);
//		output_polynomial = viscosity_over_density*(4*a+2*b+m+2*c+t); // only diffusion
		x -= 0.5*mesh_width_x1;
		printf("polynomial   = %f \n", output_polynomial);
		printf("u_1_momentum = %f \n", u_1_momentum[i][j][k]);

	}	
	
	if(j!=number_primary_cells_j){
		y += 0.5*mesh_width_x2;
		u_1 = a*x*x+b*x 	+c*y*y+d*y	+e*z*z+f*z;
		u_2 = g*x*x+h*x 	+l*y*y+m*y	+n*z*z+o*z;
		u_3 = p*x*x+q*x 	+r*y*y+s*y	+t*z*z+u*z;
		output_polynomial = -(u_1*(2*x*g+h)+u_2*(2*x*a+b) + 2*u_2*(2*y*l+m) + u_3*(2*z*n+o)+u_2*(2*z*t+u))		+viscosity_over_density*(2*2*l+2*g+2*n);
//		output_polynomial = viscosity_over_density*(4*h+2*g+d+2*l+u); // only diffusion
		y -= 0.5*mesh_width_x2;
		printf("polynomial   = %f \n", output_polynomial);
		printf("u_2_momentum = %f \n", u_2_momentum[i][j][k]);
	}	
	if(k!=number_primary_cells_k){
		z += 0.5*mesh_width_x3;
		u_1 = a*x*x+b*x 	+c*y*y+d*y	+e*z*z+f*z;
		u_2 = g*x*x+h*x 	+l*y*y+m*y	+n*z*z+o*z;
		u_3 = p*x*x+q*x 	+r*y*y+s*y	+t*z*z+u*z;
		output_polynomial = -(u_1*(2*x*p+q)+u_3*(2*x*a+b) + u_2*(2*y*r+s)+u_3*(2*y*l+m) + 2*u_3*(2*z*t+u))		+viscosity_over_density*(2*2*t+2*p+2*r);
//		output_polynomial = viscosity_over_density*(4*r+2*p+e+2*q+o); // only diffusion
		z -= 0.5*mesh_width_x3;
		printf("polynomial   = %f \n", output_polynomial);
		printf("u_3_momentum = %f \n", u_3_momentum[i][j][k]);
	}	



}}}
       printf("%s \n", "moment_predictor_unit_test finish");

       return 0;
}

#include "../src/headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>

/********************************************************************************/
/*  Function to test the advance flow field					*/
/*                     								*/
/*  										*/
/*  Programmer	: Coen Hennipman	       					*/
/*  Date	: 13-03-2014       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* Unit test for advance flow field 			                        */
/* Is not working properly because the input velocity is not divergence free    */
/* To test the momentum_predictor, use momentum_predictor_unit_test             */
/********************************************************************************/

EXPORT void advance_flow_field_unit_test()
{
      Array3<double> level_set;				// level-set field
      Array3<double> pressure;				// pressure field
      Array3<double> u_1_velocity_old;			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old;			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old;			// velocity field at old time level x3 direction
      Array3<double> u_1_velocity_new;			// velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new;			// velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new;			// velocity field at new time level x3 direction
      Array3<double> momentum_source_term_u_1;		// source term of the momentum equation in x1 direction
					       		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2;		// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3;		// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x1;	// source term of the momentum equation in x1 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x2;	// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> surface_tension_body_force_x3;	// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> scaled_density_u1;			// scaled density for the controlvolumes
							// of the momentum equation in x1 direction
      Array3<double> scaled_density_u2;			// scaled density for the controlvolumes
							// of the momentum equation in x2 direction
      Array3<double> scaled_density_u3;			// scaled density for the controlvolumes
							// of the momentum equation in x3 direction
      boundary_face boundary_faces[6];		// array with all the information
							// for the boundary conditions 
      double mesh_width_x1;				// grid spacing in x1 direction (uniform)
      double mesh_width_x2;				// grid spacing in x2 direction (uniform)
      double mesh_width_x3;				// grid spacing in x3 direction (uniform)
      int number_primary_cells_i;			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j;			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k;			// number of primary (pressure) cells in x3 direction
      int number_matrix_connections = 7;			// number of connections in momentum matrix				       
      vector gravity;					// gravitational acceleration vector 
      double actual_time_step_navier_stokes;	// time step used for level-set advection
							// computed from all stability restrictions and 
							// possibly subscycling
      double rho_plus_over_rho_minus;		// ratio of the densities of the two phases
      double smoothing_distance_factor;		// the smoothing distance is smoothing_distance_factor
							// times the smallest mesh width
      double rho_minus_over_mu_minus;		// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus;			// ratio of the viscosities of both phases
      double tolerance_pressure = 1e-8;	  		// the tolerance with which the system for the pressure is solved	
      int maximum_iterations_allowed_pressure=1e3;	// maximum number of iterations allowed for the
							// conjugate gradient method
      double tolerance_velocity=1e-7;	  		// the tolerance with which the system for the momentum predictor is solved	
      int maximum_iterations_allowed_velocity=1e2;	// maximum number of iterations allowed for the
							// conjugate gradient method
      int continuous_surface_force_model = 1;       	// =1, the continuous surface force model is applied
					        	// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor =1;   	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation
      vector initial_velocity;				// some initial velocity to fill the arrays


int i, j, k;  		// local variables for loop indexing
double x,y,z;		// coordinates
double x_length = 1;
double y_length = 2;
double z_length = 3; 

      double a,b,c,d,e,f,g,h,l,m,n,o,p,q,r,s,t,u;  // input coefficients
      
      double u_1, u_2, u_3; 
      double viscosity_over_density;
      double output_polynomial; 
      double scaled_density_constant = 1; 

      initial_velocity.u1 = 0.0;
      initial_velocity.u2 = 0.0;
      initial_velocity.u3 = 0.0;
 
      smoothing_distance_factor = 2; 
      rho_plus_over_rho_minus = 10.0;
      rho_minus_over_mu_minus = 1e20;
      mu_plus_over_mu_minus = 1; 
      
      
      gravity.u1 = 00.0;
      gravity.u2 = 00.0;
      gravity.u3 = 00.0;
 	
	actual_time_step_navier_stokes = 0.1; 

set_boundary_conditions(boundary_faces, initial_velocity); 
//	set_boundary_conditions(boundary_face[6], initial_velocity);
	
      number_primary_cells_i = pow(2,1); 
      number_primary_cells_j = pow(2,2); 
      number_primary_cells_k = pow(2,2); 

      	mesh_width_x1 = x_length/number_primary_cells_i; 
	mesh_width_x2 = y_length/number_primary_cells_j; 
	mesh_width_x3 = z_length/number_primary_cells_k; 

viscosity_over_density = (1/scaled_density_constant)*compute_scaled_viscosity(100,mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus);
      
      level_set.create(number_primary_cells_i+2,number_primary_cells_j+2,number_primary_cells_k+2);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, level_set, 0.1);

      pressure.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+2);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+2, pressure, 10000.0);
	// velocity old
      u_1_velocity_old.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_old.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_old.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, u_1_velocity_old, initial_velocity.u1);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, u_2_velocity_old, initial_velocity.u2);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, u_3_velocity_old, initial_velocity.u3);


// multiple input sets of parameters
//	a=1;	b=2;	c=3;	d=4;	e=5;	f=6;	g=7;	h=8;	l=9;	m=1;	n=2;	o=3;	p=4; 	q=5; 	r=6;	s=7;	t=8;	u=9; // full test
	a=0;	b=0;	c=3;	d=4;	e=5;	f=6;	g=7;	h=8;	l=0;	m=0;	n=-3;	o=3;	p=4; 	q=5; 	r=6;	s=7;	t=0;	u=0; // divergence free velocity input b+m+n = 0
//	a=0;	b=2;	c=0;	d=4;	e=0;	f=6;	g=0;	h=8;	l=0;	m=1;	n=0;	o=3;	p=0; 	q=5; 	r=0;	s=7;	t=0;	u=9; // linear test
//	a=0;	b=0;	c=0;	d=0;	e=0;	f=0;	g=0;	h=1;	l=0;	m=0;	n=0;	o=0;	p=0; 	q=0; 	r=0;	s=0;	t=0;	u=0; // random test
//	a=0;	b=1;	c=0;	d=0;	e=0;	f=0;	g=0;	h=0;	l=0;	m=0;	n=0;	o=0;	p=0; 	q=0; 	r=0;	s=0;	t=0;	u=0; // random test

//	a=1;	b=2;	c=3;	d=4;	e=5;	f=6;	g=7;	h=8;	l=9;	m=1;	n=2;	o=3;	p=4; 	q=5; 	r=6;	s=7;	t=8;	u=9;  // test for diffusion, only quadratic terms

  for(i=0;i<number_primary_cells_i+2;i++)
  {
     for(j=0;j<number_primary_cells_j+2;j++)
      {
	  for(k=0;k<number_primary_cells_k+2;k++)
	  {
	  x = x_length * i/number_primary_cells_i -0.5*mesh_width_x1;
//	  printf("x = %f \n",x);
	  y = y_length * j/number_primary_cells_j -0.5*mesh_width_x2;
//	  printf("y = %f \n",y);
	  z = z_length * k/number_primary_cells_k -0.5*mesh_width_x3;
	  
	  if(i!=number_primary_cells_i+1)
	  	{x += 0.5*mesh_width_x1;
//	  	printf("x = %f \t",x);
	  	u_1_velocity_old[i][j][k]= a*x*x+b*x 	+c*y*y+d*y	+e*z*z+f*z; 	// 1st test
//	  	u_1_velocity_old[i][j][k]= a*x*x+b*y*y+c*z*z 	+d*x*y+e*x*z +f*y*z; 	// 2nd test, only quadratic terms 
//	  	printf("u_1 = %f \n", u_1_velocity_old[i][j][k]);
	  	x -= 0.5*mesh_width_x1;}
	  if(j!=number_primary_cells_j+1)
	  	{y += 0.5*mesh_width_x2;
	  	u_2_velocity_old[i][j][k]= g*x*x+h*x 	+l*y*y+m*y	+n*z*z+o*z; 	// 1st test
//	  	u_2_velocity_old[i][j][k]= g*x*x+h*y*y+l*z*z 	+m*x*y+n*x*z +o*y*z; 	// 2nd test, only quadratic terms 	  	
	  	y -= 0.5*mesh_width_x2;}
	  if(k!=number_primary_cells_k+1)
	  	{z += 0.5*mesh_width_x3;
	  	u_3_velocity_old[i][j][k]= p*x*x+q*x 	+r*y*y+s*y	+t*z*z+u*z;	// 1st test
//	  	u_3_velocity_old[i][j][k]= p*x*x+q*y*y+r*z*z 	+s*x*y+t*x*z +u*y*z; 	// 2nd test, only quadratic terms 
	  	z -= 0.5*mesh_width_x3;}
	  }  
  
      }  
     
  } 


	// velocity new
      u_1_velocity_new.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      u_2_velocity_new.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      u_3_velocity_new.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, u_1_velocity_new, initial_velocity.u1);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, u_2_velocity_new, initial_velocity.u2);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, u_3_velocity_new, initial_velocity.u3);

	// momentum source terms
      momentum_source_term_u_1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      momentum_source_term_u_2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      momentum_source_term_u_3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, momentum_source_term_u_1, 0.0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, momentum_source_term_u_2, 0.0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, momentum_source_term_u_3, 0.0);
			    
	// surface_tension_body_force      
      surface_tension_body_force_x1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      surface_tension_body_force_x2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      surface_tension_body_force_x3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, surface_tension_body_force_x1, 0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, surface_tension_body_force_x2, 0);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, surface_tension_body_force_x3, 0);			    			    			    



	// scaled density
      scaled_density_u1.create(number_primary_cells_i+1, number_primary_cells_j+2, number_primary_cells_k+2);
      scaled_density_u2.create(number_primary_cells_i+2, number_primary_cells_j+1, number_primary_cells_k+2);
      scaled_density_u3.create(number_primary_cells_i+2, number_primary_cells_j+2, number_primary_cells_k+1);
      set_constant_matrix2(number_primary_cells_i+1, number_primary_cells_j+2,
			    number_primary_cells_k+2, scaled_density_u1, 1);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+1, 
			    number_primary_cells_k+2, scaled_density_u2, 1);
      set_constant_matrix2(number_primary_cells_i+2, number_primary_cells_j+2, 
			    number_primary_cells_k+1, scaled_density_u3, 1);

 	
	advance_flow_field(level_set, pressure,
                       u_1_velocity_old, u_2_velocity_old, u_3_velocity_old,
                         u_1_velocity_new, u_2_velocity_new, u_3_velocity_new,
                           momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,
                            surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
                             scaled_density_u1, scaled_density_u2, scaled_density_u3,
                              boundary_faces,
                                mesh_width_x1, mesh_width_x2, mesh_width_x3,
                                  number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,
                                   number_matrix_connections, gravity,
                                     actual_time_step_navier_stokes,
                                       rho_plus_over_rho_minus,
                                         smoothing_distance_factor,
                                          rho_minus_over_mu_minus, mu_plus_over_mu_minus,
                                           tolerance_pressure, maximum_iterations_allowed_pressure,
                                             tolerance_velocity, maximum_iterations_allowed_velocity,
                                                continuous_surface_force_model,
                                                 source_terms_in_momentum_predictor
	);
	
	
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
		output_polynomial = 2*u_1*(2*a*x+b) + u_1*(2*y*l+m) + u_2*(2*y*c+d) + u_1*(2*z*t+u) + u_3*(2*z*e+f)	+viscosity_over_density*(2*2*a+2*c+2*e); // 1st test
//		output_polynomial = viscosity_over_density*(4*a+2*b+m+2*c+t); // 2nd test, only quadratic terms 
		x -= 0.5*mesh_width_x1;
		printf("\n u_1_velocity_old = %f \n", u_1);		
		printf("u_1_velocity_new = %f \n", u_1_velocity_new[i][j][k]);		
		printf("new - old = %f \n", u_1_velocity_new[i][j][k]-u_1);			
		printf("polynomial   = %f \n", output_polynomial);
	}	
	
	if(j!=number_primary_cells_j){
		y += 0.5*mesh_width_x2;
		u_1 = a*x*x+b*x 	+c*y*y+d*y	+e*z*z+f*z;
		u_2 = g*x*x+h*x 	+l*y*y+m*y	+n*z*z+o*z;
		u_3 = p*x*x+q*x 	+r*y*y+s*y	+t*z*z+u*z;
		output_polynomial = u_1*(2*x*g+h)+u_2*(2*x*a+b) + 2*u_2*(2*y*l+m) + u_3*(2*z*n+o)+u_2*(2*z*t+u)		+viscosity_over_density*(2*2*l+2*g+2*n); // 1st test
//		output_polynomial = viscosity_over_density*(4*h+2*g+d+2*l+u); // 2nd test, only quadratic terms 
		y -= 0.5*mesh_width_x2;
		printf(" \n u_2_velocity_old = %f \n", u_2);		
		printf("u_2_velocity_new = %f \n", u_2_velocity_new[i][j][k]);		
		printf("new - old = %f \n", u_2_velocity_new[i][j][k]-u_2);			
		printf("polynomial   = %f \n", output_polynomial);	
	}	
	if(k!=number_primary_cells_k){
		z += 0.5*mesh_width_x3;
		u_1 = a*x*x+b*x 	+c*y*y+d*y	+e*z*z+f*z;
		u_2 = g*x*x+h*x 	+l*y*y+m*y	+n*z*z+o*z;
		u_3 = p*x*x+q*x 	+r*y*y+s*y	+t*z*z+u*z;
		output_polynomial = u_1*(2*x*p+q)+u_3*(2*x*a+b) + u_2*(2*y*r+s)+u_3*(2*y*l+m) + 2*u_3*(2*z*t+u)		+viscosity_over_density*(2*2*t+2*p+2*r); // 1st test
//		output_polynomial = viscosity_over_density*(4*r+2*p+e+2*q+o); // 2nd test, only quadratic terms 
		z -= 0.5*mesh_width_x3;
		printf("\n u_3_velocity_old = %f \n", u_3);		
		printf("u_3_velocity_new = %f \n", u_3_velocity_new[i][j][k]);		
		printf("new - old = %f \n", u_3_velocity_new[i][j][k]-u_3);			
		printf("polynomial   = %f \n", output_polynomial);
	}
}}}


//		printf("u_1_velocity_new = %f \n", u_1_velocity_new[i][j][k]);
//		printf("u_2_velocity_old = %f \n", u_2_velocity_old[i][j][k]);
//		printf("u_2_velocity_new = %f \n", u_2_velocity_new[i][j][k]);
//		printf("u_3_velocity_old = %f \n", u_3_velocity_old[i][j][k]);
//		printf("u_3_velocity_new = %f \n", u_3_velocity_new[i][j][k]);
//		printf("pressure = %f \n", pressure[i][j][k]);
	

      printf("%s \n", "advance_flow_field_unit_test finish");
      
      }





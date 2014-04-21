#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
      
/********************************************************************************/
/* Function for momentum prediction in all directions				*/
/*  										*/
/*  Programmer	: Coen Hennipman       						*/
/*  Date	: 06-03-2014    						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* Some notes about the direction and the notation are in my report.		*/
/********************************************************************************/

 EXPORT void convection_diffussion_source_terms(
      Array3<double> u_1_momentum,	 		// momentum in x1 direction
      Array3<double> u_2_momentum, 			// momentum in x2 direction
      Array3<double> u_3_momentum, 			// momentum in x3 direction
      Array3<double> u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction
      Array3<double> scaled_density_u1,			// scaled density for the controlvolumes
      Array3<double> scaled_density_u2,			// scaled density for the controlvolumes
      Array3<double> scaled_density_u3,			// scaled density for the controlvolumes      
      
      Array3<double> momentum_source_term_u_1,		// source term of the momentum equation in x1 direction
					       		// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_2,		// source term of the momentum equation in x2 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> momentum_source_term_u_3,		// source term of the momentum equation in x3 direction
					        	// defined on all u1 points (including boundaries)
      Array3<double> level_set, 			// level-set field
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      double smoothing_distance_factor,			// the smoothing distance is smoothing_distance_factor
      double rho_plus_over_rho_minus,
      double rho_minus_over_mu_minus,			// this was the 'Reynolds' number
							// in the original implementation of Sander
      double mu_plus_over_mu_minus,			// ratio of the viscosities of both phases
      int source_terms_in_momentum_predictor    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation  
       )
{
int i, j, k;
double viscosity_right;
double viscosity_left;
double viscosity_back;
double viscosity_front; 
double viscosity_top;
double viscosity_bottom;  

    // momentum prediction for all u
    // The for loop combined with the if statement is over all the interior velocities
for(i=1;i<number_primary_cells_i+1;i++)
{
  for(j=1;j<number_primary_cells_j+1;j++)
  {
     for(k=1;k<number_primary_cells_k+1;k++)
     {
        if(i!=number_primary_cells_i+1)
	{
		// Moment prediction for u direction        
	      	// calculation of the viscosities, the level_set-value on the faces is the average of the surrounding cells 
	      	// first perpendicular direction is x2 and second perpendicular direction is x3      
		viscosity_right = compute_scaled_viscosity(level_set[i+1][j][k],mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus);	// +1 in tangential
     		viscosity_left 	= compute_scaled_viscosity(level_set[i][j][k],mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus);	// same for every flow 
     		viscosity_back 	= compute_scaled_viscosity((level_set[i][j+1][k]+level_set[i+1][j+1][k]+level_set[i][j][k]+level_set[i+1][j][k])/4, 			// 0,+1 in tangential; 0,+1 in first perpendicular direction
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus);			
     		viscosity_front = compute_scaled_viscosity((level_set[i][j-1][k]+level_set[i+1][j-1][k]+level_set[i][j][k]+level_set[i+1][j][k])/4,			// 0,+1 in tangential; 0,-1 in first perpendicular direction 
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus); 
     		viscosity_top 	= compute_scaled_viscosity((level_set[i][j][k+1]+level_set[i+1][j][k+1]+level_set[i][j][k]+level_set[i+1][j][k])/4,			// 0,+1 in tangential; 0,+1 in second perpendicular direction
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus);
     		viscosity_bottom= compute_scaled_viscosity((level_set[i][j][k-1]+level_set[i+1][j][k-1]+level_set[i][j][k]+level_set[i+1][j][k])/4,			// 0,+1 in tangential; 0,-1 in second perpendicular direction 
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus); 
	      
		// momentum_predictor_cell_face( 
		// u_cell_face,u_left,u_right,u_back, u_front, u_bottom, u_top, 
		// v_left_back,v_right_back,v_left_front,v_right_front, 
		// w_left_top,w_right_top,w_left_bottom,w_right_bottom, 
		// mesh_width_x1,mesh_width_x2,mesh_width_x3, 
		// density,momentum_source_term 
		// viscosity_left,viscosity_right,viscosity_back,viscosity_front,viscosity_top,viscosity_bottom) 
		u_1_momentum[i][j][k] = momentum_predictor_cell_face( 
		u_1_velocity_old[i][j][k],u_1_velocity_old[i-1][j][k],u_1_velocity_old[i+1][j][k],u_1_velocity_old[i][j+1][k],u_1_velocity_old[i][j-1][k],u_1_velocity_old[i][j][k+1],u_1_velocity_old[i][j][k-1], 
		// [i][j][k]             , -1 in tangential          ,+1 in tangential           , +1 in first perpendicular , -1 in first perpendicular , +1 in second perpendicular, -1 in second perpendicular; tangential velocity
		u_2_velocity_old[i][j][k],u_2_velocity_old[i+1][j][k],u_2_velocity_old[i][j-1][k],u_2_velocity_old[i+1][j-1][k], 
		// [i][j][k]             , +1 in tangential          , -1 in first perpendicular , +1 tang -1 first perpendicular; first perpendicular velocity
		u_3_velocity_old[i][j][k],u_3_velocity_old[i+1][j][k],u_3_velocity_old[i][j][k-1],u_3_velocity_old[i+1][j][k-1], 
		// [i][j][k]             , +1 in tangential          , -1 in second perpendicular, +1 tang -1 second perpendicular; second perpendicular velocity
		mesh_width_x1,mesh_width_x2,mesh_width_x3, // tangential, first perpendicular, second perpendicular
		scaled_density_u1[i][j][k], momentum_source_term_u_1[i][j][k],
		viscosity_left,viscosity_right,viscosity_back,viscosity_front,viscosity_top,viscosity_bottom,
		source_terms_in_momentum_predictor);
	}
	
	if(j!=number_primary_cells_j+1)
	{
		// Moment prediction for v direction        
	      	// calculation of the viscosities, the level_set-value on the faces is the average of the surrounding cells      
	      	// first perpendicular direction is x1 and second perpendicular direction is x3      
		viscosity_right = compute_scaled_viscosity(level_set[i][j+1][k],mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus); 	// +1 in tangential
     		viscosity_left 	= compute_scaled_viscosity(level_set[i][j][k],mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus);	// same for every flow
     		viscosity_back 	= compute_scaled_viscosity((level_set[i+1][j][k]+level_set[i+1][j+1][k]+level_set[i][j][k]+level_set[i][j+1][k])/4,			// 0,+1 in tangential; 0,+1 in first perpendicular direction
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus);			
     		viscosity_front = compute_scaled_viscosity((level_set[i-1][j][k]+level_set[i-1][j+1][k]+level_set[i][j][k]+level_set[i][j+1][k])/4,			// 0,+1 in tangential; 0,-1 in first perpendicular direction 
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus); 		
     		viscosity_top 	= compute_scaled_viscosity((level_set[i][j][k+1]+level_set[i][j+1][k+1]+level_set[i][j][k]+level_set[i][j+1][k])/4,			// 0,+1 in tangential; 0,+1 in second perpendicular direction
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus);
     		viscosity_bottom= compute_scaled_viscosity((level_set[i][j][k-1]+level_set[i][j+1][k-1]+level_set[i][j][k]+level_set[i][j+1][k])/4,			// 0,+1 in tangential; 0,-1 in second perpendicular direction 
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus); 
	      
		// momentum_predictor_cell_face( 
		// u_cell_face,u_left,u_right,u_back, u_front, u_bottom, u_top, 
		// v_left_back,v_right_back,v_left_front,v_right_front, 
		// w_left_top,w_right_top,w_left_bottom,w_right_bottom, 
		// mesh_width_x1,mesh_width_x2,mesh_width_x3, 
		// density,momentum_source_term
		// viscosity_left,viscosity_right,viscosity_back,viscosity_front,viscosity_top,viscosity_bottom) 
		u_2_momentum[i][j][k] = momentum_predictor_cell_face( 
		u_2_velocity_old[i][j][k],u_2_velocity_old[i][j-1][k],u_2_velocity_old[i][j+1][k],u_2_velocity_old[i+1][j][k],u_2_velocity_old[i-1][j][k],u_2_velocity_old[i][j][k+1],u_2_velocity_old[i][j][k-1], 
		// [i][j][k]             , -1 in tangential          ,+1 in tangential           , +1 in first perpendicular , -1 in first perpendicular , +1 in second perpendicular, -1 in second perpendicular; tangential velocity 
		u_1_velocity_old[i][j][k],u_1_velocity_old[i][j+1][k],u_1_velocity_old[i-1][j][k],u_1_velocity_old[i-1][j+1][k], 
		// [i][j][k]             , +1 in tangential          , -1 in first perpendicular , +1 tang -1 first perpendicular; first perpendicular velocity
		u_3_velocity_old[i][j][k],u_3_velocity_old[i][j+1][k],u_3_velocity_old[i][j][k-1],u_3_velocity_old[i][j+1][k-1], 
		// [i][j][k]             , +1 in tangential          , -1 in second perpendicular, +1 tang -1 second perpendicular; second perpendicular velocity
		mesh_width_x2,mesh_width_x1,mesh_width_x3, // tangential, first perpendicular, second perpendicular
		scaled_density_u2[i][j][k], momentum_source_term_u_2[i][j][k],
		viscosity_left,viscosity_right,viscosity_back,viscosity_front,viscosity_top,viscosity_bottom,
		source_terms_in_momentum_predictor);
	}

	if(k!=number_primary_cells_k+1)
	{
		// Moment prediction for w direction        
	      	// calculation of the viscosities, the level_set-value on the faces is the average of the surrounding cells      
	      	// first perpendicular direction is x1 and second perpendicular direction is x2      
		viscosity_right = compute_scaled_viscosity(level_set[i][j][k+1],mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus); 	// +1 in tangential
     		viscosity_left 	= compute_scaled_viscosity(level_set[i][j][k],mesh_width_x1,mesh_width_x2,mesh_width_x3, smoothing_distance_factor, rho_minus_over_mu_minus,mu_plus_over_mu_minus);	// same for every flow
     		viscosity_back 	= compute_scaled_viscosity((level_set[i+1][j][k]+level_set[i+1][j][k+1]+level_set[i][j][k]+level_set[i][j][k+1])/4,			// 0,+1 in tangential; 0,+1 in first perpendicular direction
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus);			
     		viscosity_front = compute_scaled_viscosity((level_set[i-1][j][k]+level_set[i-1][j][k+1]+level_set[i][j][k]+level_set[i][j][k+1])/4,			// 0,+1 in tangential; 0,-1 in first perpendicular direction 
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus); 
     		viscosity_top 	= compute_scaled_viscosity((level_set[i][j+1][k]+level_set[i][j+1][k+1]+level_set[i][j][k]+level_set[i][j][k+1])/4,			// 0,+1 in tangential; 0,+1 in second perpendicular direction
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus);
     		viscosity_bottom= compute_scaled_viscosity((level_set[i][j-1][k]+level_set[i][j-1][k+1]+level_set[i][j][k]+level_set[i][j][k+1])/4,			// 0,+1 in tangential; 0,-1 in second perpendicular direction 
     								mesh_width_x1,mesh_width_x2,mesh_width_x3,smoothing_distance_factor,rho_minus_over_mu_minus,mu_plus_over_mu_minus); 
	      
		// momentum_predictor_cell_face( 
		// u_cell_face,u_left,u_right,u_back, u_front, u_bottom, u_top, 
		// v_left_back,v_right_back,v_left_front,v_right_front, 
		// w_left_top,w_right_top,w_left_bottom,w_right_bottom, 
		// mesh_width_x1,mesh_width_x2,mesh_width_x3, 
		// density,momentum_source_term
		// viscosity_left,viscosity_right,viscosity_back,viscosity_front,viscosity_top,viscosity_bottom) 
		u_3_momentum[i][j][k] = momentum_predictor_cell_face( 
		u_3_velocity_old[i][j][k],u_3_velocity_old[i][j][k-1],u_3_velocity_old[i][j][k+1],u_3_velocity_old[i+1][j][k],u_3_velocity_old[i-1][j][k],u_3_velocity_old[i][j+1][k],u_3_velocity_old[i][j-1][k], 
		// [i][j][k]             , -1 in tangential          ,+1 in tangential           , +1 in first perpendicular , -1 in first perpendicular , +1 in second perpendicular, -1 in second perpendicular; tangential velocity 
		u_1_velocity_old[i][j][k],u_1_velocity_old[i][j][k+1],u_1_velocity_old[i-1][j][k],u_1_velocity_old[i-1][j][k+1], 
		// [i][j][k]             , +1 in tangential          , -1 in first perpendicular , +1 tang -1 first perpendicular; first perpendicular velocity
		u_2_velocity_old[i][j][k],u_2_velocity_old[i][j][k+1],u_2_velocity_old[i][j-1][k],u_2_velocity_old[i][j-1][k+1], 
		// [i][j][k]             , +1 in tangential          , -1 in second perpendicular, +1 tang -1 second perpendicular; second perpendicular velocity
		mesh_width_x3,mesh_width_x1,mesh_width_x2, // tangential, first perpendicular, second perpendicular
		scaled_density_u3[i][j][k], momentum_source_term_u_3[i][j][k], 
		viscosity_left,viscosity_right,viscosity_back,viscosity_front,viscosity_top,viscosity_bottom,
		source_terms_in_momentum_predictor);
        }
     }
  }
}

}

      
/********************************************************************************/
/* Function for the predictor of the velocity on the cell face, 		*/
/* without boundary conditions							*/
/*  										*/
/*  Programmer	: Coen Hennipman       						*/
/*  Date	: 06-03-2014    						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* This function calculates the effect of advection and viscosity based on	*/
/* the defined input in the momentum_predictor.					*/
/* In the last step the momentum_source_term is added if: 			*/
/* source_terms_in_momentum_predictor = 1 					*/
/********************************************************************************/

 EXPORT double momentum_predictor_cell_face(
      double u_cell_face,				// the velocity on the cell face of the previous time step 
      double u_left,					// the velocity on the left of the momentum cell
      double u_right, 					// the velocity on the right of the momentum cell 
      double u_back,					// the velocity on the back of the momentum cell
      double u_front, 					// the velocity on the front of the momentum cell 
      double u_top,					// the velocity on the top of the momentum cell
      double u_bottom, 					// the velocity on the bottom of the momentum cell 
      
      double v_left_back,				// the velocity in the first perpendicular direction on the left back
      double v_right_back,				// the velocity in the first perpendicular direction on the right back
      double v_left_front,				// the velocity in the first perpendicular direction on the left front
      double v_right_front,				// the velocity in the first perpendicular direction on the right front
      
      double w_left_top,				// the velocity in the second perpendicular direction on the left top
      double w_right_top,				// the velocity in the second perpendicular direction on the right top
      double w_left_bottom,				// the velocity in the second perpendicular direction on the left bottom
      double w_right_bottom,				// the velocity in the second perpendicular direction on the right bottom      
      
      double mesh_width_x1,				// grid spacing in direction of velocity(uniform)
      double mesh_width_x2,				// grid spacing in first perpendicular direction (uniform)
      double mesh_width_x3,				// grid spacing in second perpendicular direction (uniform)
            
      double density, 					// the density of the momentum cell
      double momentum_source_term, 			// source term of the momentum equation defined on all velocity points (including boundaries) 

      double viscosity_left, 				// the viscosity on the left side of the cell
      double viscosity_right, 				// the viscosity on the right side of the cell
      double viscosity_back, 				// the viscosity on the back side of the cell
      double viscosity_front, 				// the viscosity on the front side of the cell
      double viscosity_top, 				// the viscosity on the top side of the cell
      double viscosity_bottom, 				// the viscosity on the bottom side of the cell
      
      int source_terms_in_momentum_predictor    	// =1, the source terms are applied in the momentum predictor
					        	// equation
					        	// =0, the source terms are applied in the momentum corrector
					        	// equation  				      				
       )
  {
     double advection_term;
     double diffusion_term;
     double advection_and_diffusion; 
     double advection_diffusion_and_momentum;

     double one_over_dx1 = 1/mesh_width_x1;
     double one_over_dx2 = 1/mesh_width_x2;
     double one_over_dx3 = 1/mesh_width_x3;

 	// calculation of the advection term
 	// the discretization is based on Morinishi, Fully Conservative Higher Order Finite Differnce Schemes for Incompressible Flow
	advection_term =  one_over_dx1/4.0*(pow((u_right+u_cell_face),2)-pow((u_cell_face+u_left),2));
	advection_term += one_over_dx2/4.0*((u_cell_face+u_back)*(v_left_back+v_right_back)-(u_front+u_cell_face)*(v_left_front+v_right_front));		
	advection_term += one_over_dx3/4.0*((u_cell_face+u_top)*(w_left_top+w_right_top)-(u_bottom+u_cell_face)*(w_left_bottom+w_right_bottom));
		  
	// calculation of the diffusion term
	diffusion_term = 2.0*pow(one_over_dx1,2)*(viscosity_right*(u_right-u_cell_face)-viscosity_left*(u_cell_face-u_left));
	diffusion_term += one_over_dx2*(viscosity_back*((u_back-u_cell_face)*one_over_dx2+(v_right_back-v_left_back)*one_over_dx1)-viscosity_front*((u_cell_face-u_front)*one_over_dx2+(v_right_front-v_left_front)*one_over_dx1));
	diffusion_term += one_over_dx3*(viscosity_top*((u_top-u_cell_face)*one_over_dx3+(w_right_top-w_left_top)*one_over_dx1)-viscosity_bottom*((u_cell_face-u_bottom)*one_over_dx3+(w_right_bottom-w_left_bottom)*one_over_dx1));
      
  	if(source_terms_in_momentum_predictor)
      	{
      	// u(n+1)[ijk] = u(n)[ijk] + L(u,S[ijk]) -> advection_diffusion_and_momentum = L(u)[ijk], where S is momentum_source_term 
      	advection_diffusion_and_momentum= -advection_term + 1/density*diffusion_term + momentum_source_term; 
	return advection_diffusion_and_momentum;
      	}
      	else
      	{
      	// u(n+1)[ijk] = u(n)[ijk] + L(u) -> advection_diffusion = L(u), only advection_and_diffusion is calculated here, the momentum_source_term should be addded somewhere else  
      	advection_and_diffusion 	= -advection_term + 1/density*(diffusion_term); 
      	return advection_and_diffusion;
      	}
  }      










      	

class vector
{
public:
  double u1,u2,u3;
};
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
enum variable{velocity_u1, velocity_u2, velocity_u3, level_set, pressure};
enum boundary_conditions_type{dirichlet, neumann, periodic};
enum boundary_conditions_rule{constant, function};
enum cell_centerings{cell_centered, vertex_centered};


class boundary_variable
{
public:
  variable variable_name;
  boundary_conditions_type boundary_condition_type;
  boundary_conditions_rule boundary_condition_rule;
  cell_centerings cell_centering;
  double boundary_condition_value;
  boundary_variable(variable varname, boundary_conditions_type bound_type,
				     boundary_conditions_rule bound_rule,
				     cell_centerings  cell_cent,
					double bound_value );
  boundary_variable(variable varname);
};

class boundary_face
{
public:
    boundary_variable boundary_variables[5];
    boundary_face(void);
   
};
/********************************************************************************/
/********************************************************************************/
/*  Function to solve for the correction to the velocity field to project it    */
/*  on the space of divergence free fields					*/
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* Note that in all cases the PRESSURE is computed in the pressure correction   */
/* equation, not the pressure correction.                                       */
/* Currently we assume a Dirichlet boundary condition for all normal velocities.*/
/********************************************************************************/

 void solve_momentum_corrector(
      double ***level_set,					// level-set field
      double ***pressure,					// pressure field
      double ***momentum_source_term_u_1,			// source term of the momentum equation in x1 direction
					        		// defined on all u1 points (including boundaries)
      double ***momentum_source_term_u_2,  			// source term of the momentum equation in x2 direction
					        		// defined on all u1 points (including boundaries)
      double ***momentum_source_term_u_3,			// source term of the momentum equation in x3 direction
					        		// defined on all u1 points (including boundaries)
      double ***surface_tension_body_force_x1,  		// source term of the momentum equation in x1 direction
					        		// defined on all u1 points (including boundaries)
      double ***surface_tension_body_force_x2,  		// source term of the momentum equation in x2 direction
					        		// defined on all u1 points (including boundaries)
      double ***surface_tension_body_force_x3,			// source term of the momentum equation in x3 direction
					        		// defined on all u1 points (including boundaries)
      double ***scaled_density_u1,				// scaled density for the controlvolumes
								// of the momentum equation in x1 direction
      double ***scaled_density_u2,				// scaled density for the controlvolumes
								// of the momentum equation in x2 direction
      double ***scaled_density_u3,				// scaled density for the controlvolumes
								// of the momentum equation in x3 direction
      double ***u_1_velocity_star, 	        		// velocity field at star time level x1 direction
      double ***u_2_velocity_star, 	       		 	// velocity field at star time level x2 direction
      double ***u_3_velocity_star,	        		// velocity field at star time level x3 direction
      double mesh_width_x1,		        		// grid spacing in x1 direction (uniform)
      double mesh_width_x2,		        		// grid spacing in x2 direction (uniform)
      double mesh_width_x3,		        		// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,	        		// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,	        		// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,	        		// number of primary (pressure) cells in x3 direction
      vector gravity,						// gravitational acceleration vector 
      double tolerance,	  					// the tolerance with which the system is solved	
      double actual_time_step_navier_stokes,    		// actual time step for Navier-Stokes solution algorithm 
      double rho_plus_over_rho_minus,	        		// ratio of density where (level set >0) and 
					        		// density where (level set < 0)
      int continuous_surface_force_model,       		// =1, the continuous surface force model is applied
					        		// =0, the exact interface boundary conditions are applied
      int source_terms_in_momentum_predictor,   		// =1, the source terms are applied in the momentum predictor
					        		// equation
					        		// =0, the source terms are applied in the momentum corrector
					        		// equation
      int maximum_iterations_allowed,	 			// maximum number of iterations allowed for the
								// conjugate gradient method
      boundary_face boundary_faces[6]				// array with all the information
								// for the boundary conditions 
	
      )
   
   
   {
  void build_pressure_system(					// build system of equations for the
      double **pressure_matrix,    				// pressure equation
      double *pressure_rhs,	     	
      double ***level_set_new,			
      double ***momentum_source_term_u_1,		        
      double ***momentum_source_term_u_2,  		        
      double ***momentum_source_term_u_3,  		        
      double ***surface_tension_body_force_x1,	 	        
      double ***surface_tension_body_force_x2,	                
      double ***surface_tension_body_force_x3,		        
      double ***scaled_density_u1,
      double ***scaled_density_u2,
      double ***scaled_density_u3,
      double ***u_1_velocity_star, 	        
      double ***u_2_velocity_star, 	        
      double ***u_3_velocity_star,	        
      double ***pressure_boundary_condition_x1,
      double ***pressure_boundary_condition_x2,
      double ***pressure_boundary_condition_x3,
      double mesh_width_x1,		        
      double mesh_width_x2,		        
      double mesh_width_x3,		        
      int number_primary_cells_i,	        
      int number_primary_cells_j,	        
      int number_primary_cells_k,	        
      double actual_time_step_navier_stokes,    
      double rho_plus_over_rho_minus,	        
      int continuous_surface_force_model,       
      int source_terms_in_momentum_predictor,   
      vector gravity				
       );
  void solve_pressure_correction_system(       // solve pressure equation
      double **pressure_matrix, 		  
      double  *pressure_rhside,		  
      double  ***pressure,			  
      int number_primary_cells_i,	 	  
      int number_primary_cells_j,	  	  
      int number_primary_cells_k,		  
      double   tolerance_pressure,	  		  
      int maximum_iterations_allowed_pressure	 	  
      );
   
  void apply_pressure_correction(		// apply pressure correction to velocity field
      double ***level_set,
      double ***pressure,
      double ***u_1_velocity_star,
      double ***u_2_velocity_star,
      double ***u_3_velocity_star,
      double ***surface_tension_body_force_x1,
      double ***surface_tension_body_force_x2,
      double ***surface_tension_body_force_x3,
      double ***momentum_source_term_u_1,
      double ***momentum_source_term_u_2,
      double ***momentum_source_term_u_3,
      double ***scaled_density_u1,
      double ***scaled_density_u2,
      double ***scaled_density_u3,
      int number_primary_cells_i,
      int number_primary_cells_j,
      int number_primary_cells_k,
      double mesh_width_x1,
      double mesh_width_x2,
      double mesh_width_x3,
      double rho_plus_over_rho_minus,
      double actual_time_step_navier_stokes,
      vector gravity
      );
    
  void apply_boundary_conditions_velocity(            // apply boundary conditions to velocity field
      boundary_face boundary_faces[6],
      double ***u_1_velocity,
      double ***u_2_velocity,
      double ***u_3_velocity,
      double mesh_width_x1,
      double mesh_width_x2,
      double mesh_width_x3,
      int number_primary_cells_i,
      int number_primary_cells_j,
      int number_primary_cells_k
      );
  

  double **double_Matrix(					// allocate memory for a two
      int number_primary_cells_i,				// dimensional array
      int number_primary_cells_j
      );
  void free_double_Matrix( 					// deallocate memory for a two
      double **doubleMatrix, 					// dimensional array
      int number_primary_cells_i
      );
  double *double_Vector(					// allocate memory for a one
      int number_primary_cells_i				// dimensional array
      );
  double ***double_Matrix2(                               // allocate memory for a three-
      int number_primary_cells_i,                         // dimensional array of doubles
      int number_primary_cells_j,
      int number_primary_cells_k
      );

  void free_double_Matrix2(                               // deallocate memory for a three
      double ***doubleMatrix2,                            // dimensional array of doubles
      int number_primary_cells_i,
      int number_primary_cells_j
      );
  
  void free_double_Vector(					// deallocate memory for a one
      double *double_Vector					// dimensional array
      );
  void apply_boundary_conditions_pressure(
      double ***pressure,
      double ***pressure_boundary_condition_x1,
      double ***pressure_boundary_condition_x2,
      double ***pressure_boundary_condition_x3,
      double ***scaled_density_u1,
      double ***scaled_density_u2,
      double ***scaled_density_u3,
      double mesh_width_x1,
      double mesh_width_x2,
      double mesh_width_x3,
      int number_primary_cells_i,
      int number_primary_cells_j,
      int number_primary_cells_k
      );
      
      double **pressure_matrix; 					// pressure matrix
      double  *pressure_rhside;					// pressure right hand side
      double ***pressure_boundary_condition_x1;
      double ***pressure_boundary_condition_x2;
      double ***pressure_boundary_condition_x3;
      int total_number_pressure_points;				// total number of points with pressure

    /* allocate memory for the pressure correction matrix and right hand side */
   
      total_number_pressure_points=number_primary_cells_i*number_primary_cells_j*number_primary_cells_k;
    
      pressure_matrix=double_Matrix(4,total_number_pressure_points);
      pressure_rhside=double_Vector(total_number_pressure_points);

    /* allocate memory for the pressure boundary condition */
   
      pressure_boundary_condition_x1=double_Matrix2(2,
                                       number_primary_cells_j+2,
                                         number_primary_cells_k+2);
      pressure_boundary_condition_x2=double_Matrix2(2,
                                       number_primary_cells_i+2,
                                         number_primary_cells_k+2);
      pressure_boundary_condition_x3=double_Matrix2(2,
                                       number_primary_cells_i+2,
                                         number_primary_cells_j+2);
   
   /* build the system of equations for the pressure correction equation */
  
      build_pressure_system(pressure_matrix,   
					pressure_rhside, level_set,			
					  momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3,  		        
					   surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,
					    scaled_density_u1, scaled_density_u2, scaled_density_u3,
					     u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,
                                              pressure_boundary_condition_x1, pressure_boundary_condition_x2,
                                                pressure_boundary_condition_x3,
						  mesh_width_x1, mesh_width_x2, mesh_width_x3,		        
						   number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,	        
						     actual_time_step_navier_stokes, rho_plus_over_rho_minus,	        
						       continuous_surface_force_model, source_terms_in_momentum_predictor,   
						        gravity);

   
      /* solve the system of equations for the pressure correction equation */

      solve_pressure_correction_system(pressure_matrix, pressure_rhside, pressure,			  
					number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,		  
					  tolerance, maximum_iterations_allowed);

     /* apply boundary conditions to pressure field */
   
      apply_boundary_conditions_pressure( pressure,
                                          pressure_boundary_condition_x1,
                                             pressure_boundary_condition_x2,
                                                pressure_boundary_condition_x3,
                                                  scaled_density_u1, scaled_density_u2, scaled_density_u3,
                                                    mesh_width_x1, mesh_width_x2, mesh_width_x3,
                                                      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
   
   /* apply the pressure correction to the velocity */
   
      apply_pressure_correction( level_set, pressure,					
				u_1_velocity_star, u_2_velocity_star, u_3_velocity_star,	     		
				  surface_tension_body_force_x1, surface_tension_body_force_x2, surface_tension_body_force_x3,		
				    momentum_source_term_u_1, momentum_source_term_u_2, momentum_source_term_u_3, 
				     scaled_density_u1, scaled_density_u2, scaled_density_u3,
				      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k,				
					mesh_width_x1, mesh_width_x2, mesh_width_x3, rho_plus_over_rho_minus,			
					  actual_time_step_navier_stokes, gravity);
     
    /* allocate memory for the pressure correction matrix and right hand side */
   
      free_double_Matrix(pressure_matrix,4);
      free_double_Vector(pressure_rhside);

    /* deallocate the memory for the pressure boundary conditions and the */

      free_double_Matrix2(pressure_boundary_condition_x1,2, number_primary_cells_j+2);
      free_double_Matrix2(pressure_boundary_condition_x2,2, number_primary_cells_i+2);
      free_double_Matrix2(pressure_boundary_condition_x3,2, number_primary_cells_i+2);

     /* apply boundary conditions to velocity field */
     
      apply_boundary_conditions_velocity( boundary_faces,		
					  u_1_velocity_star, u_2_velocity_star, u_3_velocity_star, 			
					    mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					      number_primary_cells_i, number_primary_cells_j, number_primary_cells_k);
      
  }

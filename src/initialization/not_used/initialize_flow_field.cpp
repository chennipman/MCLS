#include "../headers/array.h"
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
/*  Function to initialize the flow field                                       */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* For the moment the whole solution is initially set to zero, this should be   */
/* extended to more advanced cases, e.g. a parabolic profile etc.		*/
/********************************************************************************/
void initialize_flow_field(
      Array3<double> u_1_velocity_new, 			// velocity field at new time level x1 direction
      Array3<double> u_2_velocity_new, 			// velocity field at new time level x2 direction
      Array3<double> u_3_velocity_new,			// velocity field at new time level x3 direction
      Array3<double> u_1_velocity_old, 			// velocity field at old time level x1 direction
      Array3<double> u_2_velocity_old, 			// velocity field at old time level x2 direction
      Array3<double> u_3_velocity_old,			// velocity field at old time level x3 direction
      Array3<double> pressure,				// pressure
      boundary_face boundary_faces[6],			// array with all the information
							// for the boundary conditions 
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k			// number of primary (pressure) cells in x3 direction
     )
   /* function definitions */
    {

      void set_constant_vector(				// set vector to constant value 
	int vector_length,	
	Array1<double> vector_to_set,	
	double constant_value	
     );
      void apply_boundary_conditions_velocity(         // apply boundary conditions to velocity field
	  boundary_face boundary_faces[6],		
	  Array3<double> u_1_velocity, 			
	  Array3<double> u_2_velocity, 			
	  Array3<double> u_3_velocity, 			
	  double mesh_width_x1,				
	  double mesh_width_x2,				
	  double mesh_width_x3,				
	  int number_primary_cells_i,			
	  int number_primary_cells_j,			
	  int number_primary_cells_k			
     );
      
      int vector_length_u1;			// number of unknowns in array u_1_velocity 
      int vector_length_u2;			// number of unknowns in array u_2_velocity 
      int vector_length_u3;			// number of unknowns in array u_3_velocity
      int vector_length_pressure;		// number of unknowns in array pressure
      
      /* determine the number of unknowns in the velocity arrays */
      
      vector_length_u1=(number_primary_cells_i+1)*(number_primary_cells_j+2)*(number_primary_cells_k+2);
      vector_length_u2=(number_primary_cells_i+2)*(number_primary_cells_j+1)*(number_primary_cells_k+2);
      vector_length_u3=(number_primary_cells_i+2)*(number_primary_cells_j+2)*(number_primary_cells_k+1);
      vector_length_pressure=(number_primary_cells_i+2)*(number_primary_cells_j+2)*(number_primary_cells_k+2);
       
      std::cout << "tot hier 0 \n";

  
    /* initialize the velocity field */
    
      /* old time level */
    
      set_constant_vector(vector_length_u1, **u_1_velocity_new,0.0);
      set_constant_vector(vector_length_u2, **u_2_velocity_new,0.0);
      set_constant_vector(vector_length_u3, **u_3_velocity_new,0.0);
      
      std::cout << "tot hier 1 \n";

      /* new time level */
       
      set_constant_vector(vector_length_u1, **u_1_velocity_old,0.0);
      set_constant_vector(vector_length_u2, **u_2_velocity_old,0.0);
      set_constant_vector(vector_length_u3, **u_3_velocity_old,0.0);
      
            std::cout << "tot hier 2";

      
      /* apply the boundary conditions to the velocity field */
    
      /* new time level */
       
      apply_boundary_conditions_velocity(boundary_faces,		
					  u_1_velocity_new, u_2_velocity_new, u_3_velocity_new, 			
					    mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					      number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);

      /* old time level */
       
      apply_boundary_conditions_velocity(boundary_faces,		
					  u_1_velocity_old, u_2_velocity_old, u_3_velocity_old, 			
					    mesh_width_x1, mesh_width_x2, mesh_width_x3,				
					      number_primary_cells_i, number_primary_cells_j,number_primary_cells_k);
    
    /* intialize the pressure */
    
      set_constant_vector(vector_length_pressure, **pressure, 0.0);
      
    }
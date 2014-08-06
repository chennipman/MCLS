/*%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*         classes 	    */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#pragma once

#include <string>

class coordinate
{
public:
    double x1,x2,x3;

    coordinate( double x1 = 0, double x2 = 0, double x3 = 0 )
        : x1( x1 ), x2( x2 ), x3( x3 )
    {}
};

class vector
{
public:
    double u1,u2,u3;

    vector( double u1 = 0, double u2 = 0, double u3 = 0)
        : u1( u1 ), u2( u2 ), u3( u3 )
    {}
};

const coordinate default_coordinate;

class bubble
{
public:
    double principle_axis_x1;
    double principle_axis_x2;
    double principle_axis_x3;
    int label;
    coordinate center_location;

    bubble( int label = 1, coordinate center_location = default_coordinate,
            double principle_axis_x1 = 0, double principle_axis_x2 = 0, double principle_axis_x3 = 0 )
        : principle_axis_x1( principle_axis_x1 ),
          principle_axis_x2( principle_axis_x2 ),
          principle_axis_x3( principle_axis_x3 ),
          label( label ),
          center_location( center_location )
    {}
};

enum variable {velocity_u1, velocity_u2, velocity_u3, level_set, pressure_field};
enum boundary_conditions_type {dirichlet, neumann, periodic};
enum boundary_conditions_rule {constant, function};
enum cell_centerings {cell_centered, vertex_centered};
enum geometry {bubbly_flow, wavy_flow, rayleigh_taylor};
enum time_stepping_methods {none, explicit_euler, imex, runge_kutta, two_pres_solve, two_pres_solve_output};


class surface
{
public:
    int active;
    int orientation;
    double height;
};

class boundary_variable
{
public:
    variable variable_name;
    boundary_conditions_type boundary_condition_type;
    boundary_conditions_rule boundary_condition_rule;
    cell_centerings cell_centering;
    double boundary_condition_value;
    double (*custom_boundary_condition_value)( double time, double space_x, double space_y, double space_z );

    boundary_variable(variable varname=velocity_u1,
                         boundary_conditions_type bound_type=neumann,
                         boundary_conditions_rule bound_rule=constant,
                         cell_centerings  cell_cent=cell_centered,
                         double bound_value = 0.0)
//                         _custom_boundary_condition_value = NULL);
    {
        variable_name=varname;
        boundary_condition_type=bound_type;
        boundary_condition_rule=bound_rule;
        cell_centering=cell_cent;
        boundary_condition_value=bound_value;
        custom_boundary_condition_value = NULL;
    };

    boundary_variable(variable varname)
    {
        variable_name=varname;
        boundary_condition_type=neumann;
        boundary_condition_rule=constant;
        cell_centering=cell_centered;
        boundary_condition_value=0.0;
    }
     double get_boundary_condition_value( double t, double x, double y, double z)
     {
     if (this-> boundary_condition_rule == constant)
         return this -> boundary_condition_value;
     else
         return this -> custom_boundary_condition_value(t,x,y,z);
     }
};


class boundary_face
{
public:
    boundary_variable boundary_variables[5];

    boundary_face(void)
    {
        boundary_variables[0].variable_name=velocity_u1;
        boundary_variables[1].variable_name=velocity_u2;
        boundary_variables[2].variable_name=velocity_u3;
        boundary_variables[3].variable_name=level_set;
        boundary_variables[4].variable_name=pressure_field;
    }
};

class restart_parameters
{
public:
      int start_from_restart_file;		
      int write_solution_to_restart_file;
      std::string name_restart_file_to_write;
      std::string name_restart_file_to_read;

    restart_parameters()
    {
          start_from_restart_file=0;
          write_solution_to_restart_file=1;
          name_restart_file_to_read="restart_file_mcls_out";
          name_restart_file_to_write="restart_file_mcls_out";
    }
};

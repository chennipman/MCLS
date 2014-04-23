/*%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*         classes 	    */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%*/

class coordinate
{
public:
  double x1,x2,x3;
  coordinate(double xx1=1, double xx2=2, double xx3=3){x1=xx1;x2=xx2;x3=xx3;}
};

class vector
{
public:
  double u1,u2,u3;
};

const coordinate default_coordinate;

class bubble
{
public:
  double radius;
  int label;
  coordinate center_location;
  bubble(int number, coordinate bubble_center, double bubble_radius);
};
bubble::bubble(int number=1, coordinate bubble_center=default_coordinate, double bubble_radius=0)
{
  label=number;
  center_location=bubble_center;
  radius=bubble_radius;
};

enum variable{velocity_u1, velocity_u2, velocity_u3, level_set, pressure};
enum boundary_conditions_type{dirichlet, neumann, taylor_vortex, periodic};
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
// full constructor
// default is: cell-centered, (constant) homogeneous neumann
boundary_variable::boundary_variable(variable varname=velocity_u1,
				     boundary_conditions_type bound_type=neumann,
				     boundary_conditions_rule bound_rule=constant,
				     cell_centerings  cell_cent=cell_centered,
					double bound_value =0.0)
{
    variable_name=varname;
    boundary_condition_type=bound_type;
    boundary_condition_rule=bound_rule;
    cell_centering=cell_cent;
    boundary_condition_value=bound_value;
};
// constructor with only the name
boundary_variable::boundary_variable(variable varname)
{
variable_name=varname;
boundary_condition_type=neumann;
boundary_condition_rule=constant;
cell_centering=cell_centered;
boundary_condition_value=0.0;
};

class boundary_face
{
public:
    boundary_variable boundary_variables[5];
    boundary_face(void);
   
};
boundary_face::boundary_face(void)
{
    boundary_variables[0].variable_name=velocity_u1;
    boundary_variables[1].variable_name=velocity_u2;
    boundary_variables[2].variable_name=velocity_u3;
    boundary_variables[3].variable_name=level_set;
    boundary_variables[4].variable_name=pressure;

    
}

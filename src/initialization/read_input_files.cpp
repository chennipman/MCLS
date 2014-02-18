/********************************************************************************/
/********************************************************************************/
/*  Function to read all input parameters for the computation from         	*/
/*  input file									*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 03-01       							*/
/*  Update	:        							*/
/********************************************************************************/
//
//
EXPORT void read_input_files()
{
   
/*---------------------------*/  
/*      grid dimensions      */
/*---------------------------*/
int  number_primary_cells_i  /* Number of primary cells in the i-direction*/
int  number_primary_cells_j  /* Number of primary cells in the j-direction*/
int  number_primary_cells_k  /* Number of primary cells in the k-direction*/
 
 

/*------------------------------*/
/* momentum predictor equation  */   
/*------------------------------*/

int use_sources_in_predictor;  /* =1 when the momentum predictor equation contains sources, =0 when not */
			       /* sources are here: pressure gradient, surface tension and gravity   */
			       
/*----------------------------*/
/* curvature computation      */
/* algorithm 		      */
/*----------------------------*/

int number_smoothing_steps; /* number of smoothing steps in curvature smoothing algorithm */
int use_curvature_filter;     /* =1 when a curvature filter is employed, =0 when not */



/*---------------*/
/* time stepping */
/*---------------*/

double maximum_time_step_navier_stokes; /* maximum time step for Navier-Stokes solution algorithm */
double maximum_time_step_level_set;	/* maximum time step for Level-Set conservation equation */
double actual_time_step_navier_stokes;  /* actual time step for Navier-Stokes solution algorithm */
double actual_time_step_level_set;	/* actual time step for Level-Set conservation equation */
double time_step_output_solution;	/* interval of writing the solution to file */

  
 }

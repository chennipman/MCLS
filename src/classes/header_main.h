#include "../headers/array.h"
/*********************************************/
/* include the necessary standard functions  */
/*********************************************/
#include <iostream>
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>


/*********************************************/
/* include the necessary classe definitions  */
/*********************************************/
#include "header_clsses.h"


/************************************************/
/* include the necessary constants definitions  */
/************************************************/
#include "header_constants.h"


/*++++++++++++++++++++++++++++++++++++*/ 
/* Primary variables for Navier-Stokes*/
/*++++++++++++++++++++++++++++++++++++*/ 

Array3<double> u_1_velocity_old; /* 3-array with x1-velocity at old time level at u1 points, old time level*/
Array3<double> u_2_velocity_old; /* 3-array with x2-velocity at old time level at u2 points, new time level*/
Array3<double> u_3_velocity_old; /* 3-array with x3-velocity at old time level at u3 points, old time level*/
Array3<double> u_1_velocity_new; /* 3-array with x1-velocity at new time level at u1 points, new time level*/
Array3<double> u_2_velocity_new; /* 3-array with x2-velocity at new time level at u2 points, old time level */
Array3<double> u_3_velocity_new; /* 3-array with x3-velocity at new time level at u3 points, new time level*/
Array3<double> pressure_old;	    /* 3-array with pressure field at cell centers, old time level */
Array3<double> pressure_new;	    /* 3-array with pressure field at cell centers, new time level */
Array3<double> viscosity;	    /* 3-array with viscosity at cell centers */


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
/* Primary variables for the coupling between interface dynamics and Navier-Stokes */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 

/*++++++++++++++++++++++++++++++++++++++++++++*/  
/* Primary variables for interface capturing  */
/*++++++++++++++++++++++++++++++++++++++++++++*/  
Array3<double> volume_of_fluid_old;  /* 3-array with volume_of_fluid at cell centers, old time level*/
Array3<double> volume_of_fluid_new;  /* 3-array with volume_of_fluid at cell centers, new time level*/
Array3<double> level_set_old;        /* 3-array with level set field at cell centers, old time level*/
Array3<double> level_set_new;        /* 3-array with level set field at cell centers, new time level*/

/* parameters of level-set to  vof conversion function */

double lower_bound_level_set_derivative;  	/* lower limit of the smallest component of the gradient of */
						/* level-set, beyond which the limit formulation is used    */
int    number_iterations_level_set_conversion; 	/* maximum number of iterations allowed for the level-set    */
						/* to vof conversion					 */
double tolerance_level_set_conversion;  	/* tolerance for the level-set to vof conversion */
int    write_level_set_correction;		/* =1, write out the level-set correction, =0 do not */

/* parameters of the vof to level-set conversion function */

double vof_2_level_set_tolerance;
double volume_of_fluid_tolerance;		// the volume of fluid field should be in the interval
						// 0-volume_of_fluid_tolerance and 1+volume_of_fluid_tolerance
						// otherwise conversion to level-set is no longer possible
int    do_vof_2_level_set_relaxation;
int    number_level_set_2_vof_iterations;

/* extra information */
int    write_level_set_timing;
int    write_redistribution_output;
      
/*---------------------------*/  
/* Flow properties           */
/*---------------------------*/  

double rho_minus_over_mu_minus;  	/* The "Reynolds_number", defined as rho_m / mu_m */
double rho_plus_over_rho_minus; 	/* rho_plus / rho_minus indicated by level set */
double mu_plus_over_mu_minus;		/* mu_plus / mu_minus */ 
vector gravity;				/* gravitational acceleration vector */
double sigma_over_rho_minus;		/* sigma / rho_minus */


/* Bubble properties */

bubble gnoef[3];

/*---------------------------*/  
/*      grid dimensions      */
/*---------------------------*/
int  number_primary_cells_i;  /* Number of primary cells in the i-direction*/
int  number_primary_cells_j;  /* Number of primary cells in the j-direction*/
int  number_primary_cells_k;  /* Number of primary cells in the k-direction*/
 
 

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
int time_stepping_method;		/* select time scheme 1:explicit euler 2: imex, 3: runge-kutta */
double actual_time;			// the actual time

 
  
  
  
  
  
  

#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
/********************************************************************************/
/*  Function to initialize the level-set field for the presence of a number     */
/*  of bubbles. The presence of a free surface is treated separately            */
/*  										*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 10-03-2013       						*/
/*  Update	:        							*/
/********************************************************************************/
/* Notes									*/
/* The level-set function is assumed to be <0 inside the bubble, irrespective   */
/* of the fact that the density is < or > than the ambient density  to avoid    */
/* confusion.                                                                   */
/* Currently, the number of bubbles is limited to 10.                           */
/* Note that this function handles both elliptical and spherical bubbles.       */
/********************************************************************************/

EXPORT void initialize_bubbles(
      bubble *bubbles, 					// array with the definition of the bubbles
      int number_of_bubbles, 				// number of bubbles in the domain (<10)
      int number_primary_cells_i,			// number of primary (pressure) cells in x1 direction
      int number_primary_cells_j,			// number of primary (pressure) cells in x2 direction
      int number_primary_cells_k,			// number of primary (pressure) cells in x3 direction
      double mesh_width_x1,				// grid spacing in x1 direction (uniform)
      double mesh_width_x2,				// grid spacing in x2 direction (uniform)
      double mesh_width_x3,				// grid spacing in x3 direction (uniform)
      Array3<double> level_set				// level-set field
	    
    )
    {
    double x1_coordinate_cell_center;			// x1 coordinate cell center
    double x2_coordinate_cell_center;			// x2 coordinate cell center
    double x3_coordinate_cell_center;			// x3 coordinate cell center
    double bubble_radius;				// radius of bubble, if bubble is spherical
    double bubble_center_x1;				// coordinate of center of bubble, x1 component
    double bubble_center_x2;				// coordinate of center of bubble, x2 component
    double bubble_center_x3;				// coordinate of center of bubble, x3 component
    double principle_axis_x1;				// principal axis of bubble, x1 direction
    double principle_axis_x2;				// principal axis of bubble, x2 direction
    double principle_axis_x3;				// principal axis of bubble, x3 direction
    double distance_x1_squared;				// x1 'component' of distance
    double distance_x2_squared;				// x2 'component' of distance
    double distance_x3_squared;				// x3 'component' of distance
    double azimuth_angle;				// azimuth angle of the parametrization of ellipse
    double elevation_angle;				// elevation angle of the parametrization of ellipse
    double distance_to_all_bubbles=1000.0;		// distance of point to this bubble
    double distance_to_bubble;				// distance of point to all bubbles in an array
    double tentative_distance_to_ellipse;		// tentative distance to all interfaces
    int i_index, j_index, k_index;  			// local variables for loop indexing
    int azimuth_angle_index;				// index of azimuth angle in the parametrization of ellipse
    int elevation_angle_index;				// index of elevation angle in the parametrization of ellipse
    int bubble_index;					// index of bubble
    int inside_ellipse;					// (=1 if point is inside the ellipse, otherwise 0)
       

		      
       /* for all points determine the distance to all bubble interfaces */
		      
   
      for(i_index=0;i_index<number_primary_cells_i+2;i_index++)
      {
	  x1_coordinate_cell_center=(i_index-0.5)*mesh_width_x1;
	  
	  for(j_index=0;j_index<number_primary_cells_j+2;j_index++)
	  {
	      x2_coordinate_cell_center=(j_index-0.5)*mesh_width_x2;
	      
	      for(k_index=0;k_index<number_primary_cells_k+2;k_index++)
	      {
		x3_coordinate_cell_center=(k_index-0.5)*mesh_width_x3;
		
		  distance_to_all_bubbles=10000;
		  for(bubble_index=0;bubble_index<number_of_bubbles; bubble_index++)
		  {
		      /* assign the active bubble index and its center location */
		      
		      bubble_center_x1=bubbles[bubble_index].center_location.x1;
		      bubble_center_x2=bubbles[bubble_index].center_location.x2;
		      bubble_center_x3=bubbles[bubble_index].center_location.x3;
		      principle_axis_x1=bubbles[bubble_index].principle_axis_x1;
		      principle_axis_x2=bubbles[bubble_index].principle_axis_x2;
		      principle_axis_x3=bubbles[bubble_index].principle_axis_x3;
		      
		      /* compare the principle axes of the bubble */
		      
		      if((fabs(principle_axis_x1-principle_axis_x2)<0.000001)&&
			(fabs(principle_axis_x2-principle_axis_x3)<0.000001))
		      {
			  /* the three principal axes are equal => */
			  /* the bubble is spherical with radius equal to one of the principal axes */
			  
			  bubble_radius=principle_axis_x1;

			  /*  compute the signed distance to the sphere */
			  /* inside the sign is negative */

			  distance_to_bubble= -1.0*bubble_radius+
			  sqrt( ( x1_coordinate_cell_center-bubble_center_x1)*
				    ( x1_coordinate_cell_center-bubble_center_x1)+
				( x2_coordinate_cell_center-bubble_center_x2)*
				    ( x2_coordinate_cell_center-bubble_center_x2)+
				( x3_coordinate_cell_center-bubble_center_x3)*
				    ( x3_coordinate_cell_center-bubble_center_x3));
			  
			  /* update the distance with the distance to this bubble */
			  /* making sure the correct sign is retained 		  */
			  
			  if( abs(distance_to_bubble) <
				  abs(distance_to_all_bubbles))
			  {
				distance_to_all_bubbles=distance_to_bubble;
			  }
		      }
		      else
		      {
		      	      exit(1);
			 /* the three principal axes are not equal => */
			 /* the bubble is elliptical */
			 /* loop over a fine mesh of points on the ellipsoid to determine */
			 /* distance */
			 /* this is still to be changed to computation of the exact distance */
			 /* because it is very time consuming */
			 
			 
			 /* we first compute the unsigned distance to all points on the ellipse */
			 
			  tentative_distance_to_ellipse=100;
			  for(elevation_angle_index=0; elevation_angle_index<=100; elevation_angle_index++)
			  {
			      elevation_angle=-4.0*atan(1.0) + elevation_angle_index*atan(1.0)/25.0;
			      for(azimuth_angle_index=0; azimuth_angle_index<=200; azimuth_angle_index++)
			      {
				azimuth_angle=-4.0*atan(1.0) + azimuth_angle_index*atan(1.0)/25.0;
				distance_x1_squared=(x1_coordinate_cell_center-bubble_center_x1-
				    principle_axis_x1*cos(elevation_angle)*cos(azimuth_angle))*
					    (x1_coordinate_cell_center-bubble_center_x1-
				    principle_axis_x1*cos(elevation_angle)*cos(azimuth_angle));
				distance_x2_squared=(x2_coordinate_cell_center-bubble_center_x2-
				    principle_axis_x2*cos(elevation_angle)*cos(azimuth_angle))*
					    (x2_coordinate_cell_center-bubble_center_x2-
				    principle_axis_x2*cos(elevation_angle)*cos(azimuth_angle));
				distance_x3_squared=(x3_coordinate_cell_center-bubble_center_x3-
				    principle_axis_x3*cos(elevation_angle)*cos(azimuth_angle))*
					    (x3_coordinate_cell_center-bubble_center_x3-
				    principle_axis_x3*cos(elevation_angle)*cos(azimuth_angle));
				tentative_distance_to_ellipse=std::min(tentative_distance_to_ellipse,
				    sqrt(distance_x1_squared+distance_x2_squared+distance_x3_squared));       
			      }
			  }
			 
			 /* the distance to the bubble is the minimum of the distance to all points */
			 /* on the ellipse */
			 
			  distance_to_bubble=tentative_distance_to_ellipse;
			 
			 /* check if the point is inside the ellipse, to be able to add the correct sign*/
			 
			  inside_ellipse= (int)(x1_coordinate_cell_center-bubble_center_x1)/principle_axis_x1*
					  (x1_coordinate_cell_center-bubble_center_x1)/principle_axis_x1+
					 (x2_coordinate_cell_center-bubble_center_x2)/principle_axis_x2*
					  (x2_coordinate_cell_center-bubble_center_x2)/principle_axis_x2+
					 (x3_coordinate_cell_center-bubble_center_x3)/principle_axis_x3*
					  (x3_coordinate_cell_center-bubble_center_x3)/principle_axis_x3;
			
			/* make the unsigned distance a signed distance */		  
					  
			  if(inside_ellipse)
			  {
			    tentative_distance_to_ellipse=-1.0*tentative_distance_to_ellipse;
			  }		  
					  
			  /* update the distance with the distance to this bubble */
			  /* making sure the correct sign is retained 		  */
			  
			  if( abs(distance_to_bubble) <
				  abs(distance_to_all_bubbles))
			  {
				distance_to_all_bubbles=distance_to_bubble;
			  }

		       }
		  }
		  
		  /* the level-set function is the distance to the set of all interfaces */
		  
		  level_set[i_index][j_index][k_index]=distance_to_all_bubbles;
	      }
	  }
      }
}

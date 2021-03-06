#include "../headers/array.h"
#include<cstdlib>
#include<iostream>
#include<math.h>
/********************************************************************************/
/********************************************************************************/
/*  Function to compute the time derivative of the redistribution               */
/*  equation.                                                                   */
/*                                                                              */
/*  Programmer  : Duncan van der Heul                                           */
/*  Date        : 10-03-2013                                                    */
/*  Update      :                                                               */
/********************************************************************************/
/* Notes                                                                        */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/
EXPORT double compute_redistribution_time_derivative(
        Array3<double> redistribution_velocity_x1,                      // artificial redistribution velocity x1 direction
        Array3<double> redistribution_velocity_x2,                      // artificial redistribution velocity x2 direction
        Array3<double> redistribution_velocity_x3,                      // artificial redistribution velocity x3 direction
        Array3<double> time_derivative_volume_of_fluid_correction,      // time derivative in the discretised
                                                                        // volume of fluid redistribution equation
        Array3<double> volume_of_fluid_correction,                      // correction to the volume of fluid field
                                                                        // to make it valid
        int number_primary_cells_i,                                     // number of primary (pressure) cells in x1 direction
        int number_primary_cells_j,                                     // number of primary (pressure) cells in x2 direction
        int number_primary_cells_k,                                     // number of primary (pressure) cells in x3 direction
        double mesh_width_x1,                                           // grid spacing in x1 direction (uniform)
        double mesh_width_x2,                                           // grid spacing in x2 direction (uniform)
        double mesh_width_x3,                                           // grid spacing in x3 direction (uniform)
        double mass_redistribution_diffusion_coefficient                // diffusion coefficient for mass redistribution equation
        )
{
        double flux_minus_x1;                                           // flux of error advection equation at plus i face
        double flux_minus_x2;                                           // flux of error advection equation at plus j face
        double flux_minus_x3;                                           // flux of error advection equation at plus k face
        double flux_pluss_x1;                                           // flux of error advection equation at minus i face
        double flux_pluss_x2;                                           // flux of error advection equation at minus j face
        double flux_pluss_x3;                                           // flux of error advection equation at minus k face

        
        double one_over_dx1 =                                           // 1/(grid spacing in x1 direction)
                        1.0/(mesh_width_x1);
        double one_over_dx2 =                                           // 1/(grid spacing in x2 direction)
                        1.0/(mesh_width_x1);
        double one_over_dx3 =                                           // 1/(grid spacing in x3 direction)
                        1.0/(mesh_width_x3);
        double maximum_time_derivative_vof_error=0;                     // largest value of time derivative in the volume of fluid
                                                                        // error redistribution equation
        int i_index, j_index, k_index;                                  // local variables for loop indexing
        int max_index_i;						// i_index of cell with largest
        								// time derivative of the vof error
        int max_index_j;						// j_index of cell with largest
        								// time derivative of the vof error
        int max_index_k;						// k_index of cell with largest
        								// time derivative of the vof error
//         double diffusion_coefficient=0.00001;
//        double diffusion_coefficient=0.0000002;
        double diffusion_coefficient=0.0000004;
        mass_redistribution_diffusion_coefficient=diffusion_coefficient;
        
        /* compute time derivative of the redistribution equation */
        
            for( i_index=1;i_index<number_primary_cells_i+1;i_index++){
                for(j_index=1;j_index<number_primary_cells_j+1;j_index++){
                    for(k_index=1;k_index<number_primary_cells_k+1;k_index++){
                  
                        flux_pluss_x1=
                                upwind_flux_mass_redistribution(
                                    redistribution_velocity_x1[i_index  ][j_index][k_index],
                                        volume_of_fluid_correction[i_index  ][j_index][k_index],
                                            volume_of_fluid_correction[i_index+1][j_index][k_index])
                                    -mass_redistribution_diffusion_coefficient*one_over_dx1*
                                     (volume_of_fluid_correction[i_index+1][j_index][k_index]-
                                       volume_of_fluid_correction[i_index  ][j_index][k_index]);
                                
                        flux_minus_x1=
                                upwind_flux_mass_redistribution(
                                    redistribution_velocity_x1[i_index-1][j_index][k_index],
                                        volume_of_fluid_correction[i_index-1][j_index][k_index],
                                            volume_of_fluid_correction[i_index  ][j_index][k_index])
                                   -mass_redistribution_diffusion_coefficient*one_over_dx1*
                                     (volume_of_fluid_correction[i_index][j_index][k_index]-
                                       volume_of_fluid_correction[i_index-1][j_index][k_index]);
                        flux_pluss_x2=
                                upwind_flux_mass_redistribution(
                                    redistribution_velocity_x2[i_index][j_index  ][k_index],
                                        volume_of_fluid_correction[i_index][j_index  ][k_index],
                                            volume_of_fluid_correction[i_index][j_index+1][k_index])
                                    -mass_redistribution_diffusion_coefficient*one_over_dx2*
                                     (volume_of_fluid_correction[i_index][j_index+1][k_index]-
                                       volume_of_fluid_correction[i_index][j_index][k_index]);
                       flux_minus_x2=
                                upwind_flux_mass_redistribution(
                                    redistribution_velocity_x2[i_index][j_index-1][k_index],
                                        volume_of_fluid_correction[i_index][j_index-1][k_index],
                                            volume_of_fluid_correction[i_index][j_index][k_index])
                                    -mass_redistribution_diffusion_coefficient*one_over_dx2*
                                     (volume_of_fluid_correction[i_index][j_index][k_index]-
                                       volume_of_fluid_correction[i_index][j_index-1][k_index]);
                         flux_pluss_x3=
                                upwind_flux_mass_redistribution(
                                    redistribution_velocity_x3[i_index][j_index][k_index  ],
                                        volume_of_fluid_correction[i_index][j_index][k_index  ],
                                            volume_of_fluid_correction[i_index][j_index][k_index+1])
                                    -mass_redistribution_diffusion_coefficient*one_over_dx3*
                                     (volume_of_fluid_correction[i_index][j_index][k_index+1]-
                                       volume_of_fluid_correction[i_index][j_index][k_index]);
                        flux_minus_x3=
                                upwind_flux_mass_redistribution(
                                    redistribution_velocity_x3[i_index][j_index][k_index-1],
                                        volume_of_fluid_correction[i_index][j_index][k_index-1],
                                            volume_of_fluid_correction[i_index][j_index][k_index])
                                    -mass_redistribution_diffusion_coefficient*one_over_dx3*
                                     (volume_of_fluid_correction[i_index][j_index][k_index]-
                                       volume_of_fluid_correction[i_index][j_index][k_index-1]);
                                
                        time_derivative_volume_of_fluid_correction[i_index][j_index][k_index]=
                                              -1.0*(
                                        one_over_dx1*(flux_pluss_x1-flux_minus_x1)+
                                        one_over_dx2*(flux_pluss_x2-flux_minus_x2)+
                                        one_over_dx3*(flux_pluss_x3-flux_minus_x3)
                                                   );
                        // maximum_time_derivative_vof_error=std::max(maximum_time_derivative_vof_error,
                        //     fabs(time_derivative_volume_of_fluid_correction[i_index][j_index][k_index]));
                        
                        if(fabs(time_derivative_volume_of_fluid_correction[i_index][j_index][k_index])>=
                        	maximum_time_derivative_vof_error)
                        {
                        	maximum_time_derivative_vof_error=
                        		fabs(time_derivative_volume_of_fluid_correction[i_index][j_index][k_index]);
                        	max_index_i=i_index;
                        	max_index_j=j_index;
                        	max_index_k=k_index;
                        }
            
                    }
                }
            }
            std::cerr<<" maximum_time_derivative_vof_error : "<<maximum_time_derivative_vof_error<<"\n";
            std::cerr<<" at indices "<< max_index_i<<" "<< max_index_j<<" "<< max_index_k<<"\n";
            i_index=max_index_i;
            j_index=max_index_j;
            k_index=max_index_k;
            flux_pluss_x1=
                    upwind_flux_mass_redistribution(
                        redistribution_velocity_x1[i_index  ][j_index][k_index],
                            volume_of_fluid_correction[i_index  ][j_index][k_index],
                                volume_of_fluid_correction[i_index+1][j_index][k_index])
                        -mass_redistribution_diffusion_coefficient*one_over_dx1*
                         (volume_of_fluid_correction[i_index+1][j_index][k_index]-
                           volume_of_fluid_correction[i_index  ][j_index][k_index]);
                    
            flux_minus_x1=
                    upwind_flux_mass_redistribution(
                        redistribution_velocity_x1[i_index-1][j_index][k_index],
                            volume_of_fluid_correction[i_index-1][j_index][k_index],
                                volume_of_fluid_correction[i_index  ][j_index][k_index])
                       -mass_redistribution_diffusion_coefficient*one_over_dx1*
                         (volume_of_fluid_correction[i_index][j_index][k_index]-
                           volume_of_fluid_correction[i_index-1][j_index][k_index]);
            flux_pluss_x2=
                    upwind_flux_mass_redistribution(
                        redistribution_velocity_x2[i_index][j_index  ][k_index],
                            volume_of_fluid_correction[i_index][j_index  ][k_index],
                                volume_of_fluid_correction[i_index][j_index+1][k_index])
                        -mass_redistribution_diffusion_coefficient*one_over_dx2*
                         (volume_of_fluid_correction[i_index][j_index+1][k_index]-
                           volume_of_fluid_correction[i_index][j_index][k_index]);
            flux_minus_x2=
                    upwind_flux_mass_redistribution(
                        redistribution_velocity_x2[i_index][j_index-1][k_index],
                            volume_of_fluid_correction[i_index][j_index-1][k_index],
                                volume_of_fluid_correction[i_index][j_index][k_index])
                        -mass_redistribution_diffusion_coefficient*one_over_dx2*
                         (volume_of_fluid_correction[i_index][j_index][k_index]-
                           volume_of_fluid_correction[i_index][j_index-1][k_index]);
            flux_pluss_x3=
                    upwind_flux_mass_redistribution(
                        redistribution_velocity_x3[i_index][j_index][k_index  ],
                            volume_of_fluid_correction[i_index][j_index][k_index  ],
                                volume_of_fluid_correction[i_index][j_index][k_index+1])
                        -mass_redistribution_diffusion_coefficient*one_over_dx3*
                         (volume_of_fluid_correction[i_index][j_index][k_index+1]-
                           volume_of_fluid_correction[i_index][j_index][k_index]);
            flux_minus_x3=
                    upwind_flux_mass_redistribution(
                        redistribution_velocity_x3[i_index][j_index][k_index-1],
                            volume_of_fluid_correction[i_index][j_index][k_index-1],
                                volume_of_fluid_correction[i_index][j_index][k_index])
                        -mass_redistribution_diffusion_coefficient*one_over_dx3*
                         (volume_of_fluid_correction[i_index][j_index][k_index]-
                           volume_of_fluid_correction[i_index][j_index][k_index-1]);
            
             std::cerr<<" flux_minus_x1 : "<<    flux_minus_x1 <<"\n";   
             std::cerr<<" flux_minus_x2 : "<<    flux_minus_x2 <<"\n";
             std::cerr<<" flux_minus_x3 : "<<    flux_minus_x3 <<"\n";
             std::cerr<<" flux_pluss_x1 : "<<    flux_pluss_x1 <<"\n";
             std::cerr<<" flux_pluss_x2 : "<<    flux_pluss_x2 <<"\n";
             std::cerr<<" flux_pluss_x3 : "<<    flux_pluss_x3 <<"\n";
             std::cerr<<" check on max  : "<<     -1.0*(
                                        one_over_dx1*(flux_pluss_x1-flux_minus_x1)+
                                        one_over_dx2*(flux_pluss_x2-flux_minus_x2)+
                                        one_over_dx3*(flux_pluss_x3-flux_minus_x3)
                                                   )<< "\n";
                           
          return maximum_time_derivative_vof_error;
}
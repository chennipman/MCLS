#include<cstdlib>
#include<iostream>
#include<algorithm>
#include<math.h>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

/************************************************************************************/
/************************************************************************************/
/*  Function to compute the volume of fluid field from a given level-set field      */
/*  when the level-set function is less than zero                                   */
/*                                                                                  */
/*  Programmer       : Duncan van der Heul                                          */
/*  Date      : 10-03-2013                                                          */
/*  Update    :                                                                     */
/************************************************************************************/
/* Notes                                                                            */
/* For the case of all first partial derivatives nonzero, the standard              */
/* formulation is used. Otherwise a number of limiting cases are considered         */
/* The algorithm is discussed in paragraph 5.4 of Sander's thesis                   */
/************************************************************************************/
      int level_set_2_vof_phi_negative(
             double level_set,            // level-set field value in this cell
             double d_level_set_d_x1,            // first partial derivative in x1
                                          // direction of level-set
             double d_level_set_d_x2,            // first partial derivative in x2
                                          // direction of level-set
             double d_level_set_d_x3,            // first partial derivative in x3
                                          // direction of level-set
             double &volume_of_fluid,            // volume of fluid value in this cell
             double lower_bound_derivatives    // lower bound for the first partial derivatives
                                          // to consider it a limiting case of vanishing
                                          // partial derivatives

      )
      {
       /*
C-------------------------------------------------------------------------------
C  Psi = f_std::min(Phi, grad Phi) for Phi <= 0.
C  Section 5.4 in thesis.
C-------------------------------------------------------------------------------
*/

      double absolute_d_level_set_d_x1;          // absolute value of first partial derivative in x1
                                          // direction of level-set
      double absolute_d_level_set_d_x2;          // absolute value of first partial derivative in x2
                                          // direction of level-set
      double absolute_d_level_set_d_x3;          // absolute value of first partial derivative in x3
                                          // direction of level-set

      double level_set_value_vertex_A;           // value of LINEARIZED level-set field at vertex A
      double level_set_value_vertex_B;           // value of LINEARIZED level-set field at vertex B
      double level_set_value_vertex_C;           // value of LINEARIZED level-set field at vertex C
      double level_set_value_vertex_D;           // value of LINEARIZED level-set field at vertex D
      double level_set_value_vertex_E;           // value of LINEARIZED level-set field at vertex E

      double level_set_value_vertex_A_signed;
      double level_set_value_vertex_B_signed;
      double level_set_value_vertex_C_signed;
      double level_set_value_vertex_D_signed;
      double level_set_value_vertex_E_signed;
      double level_set_value_vertex_F_signed;

      double d_level_set_d_xi;                   // first partial derivative with respect to logical
                                          // space coordinate xi of level set field
      double d_level_set_d_eta;                  // first partial derivative with respect to logical
                                          // space coordinate eta of level set field
      double d_level_set_d_zeta;          // first partial derivative with respect to logical
                                          // space coordinate zeta of level set field


      if(level_set>0.0)
      {
       std::cerr<< "level_set value positive in level_set_2_vof_phi_negative.. \n";
       std::cerr<< "Check the source code.\n";
       return 1;
      }

      absolute_d_level_set_d_x1=fabs(d_level_set_d_x1);
      absolute_d_level_set_d_x2=fabs(d_level_set_d_x2);
      absolute_d_level_set_d_x3=fabs(d_level_set_d_x3);


      /* d_level_set_d_xi is the largest absolute value of the three first derivatives */

      d_level_set_d_xi=std::max(std::max(absolute_d_level_set_d_x1, absolute_d_level_set_d_x2), absolute_d_level_set_d_x3);

      /* d_level_set_d_zeta is the smallest absolute value of the three first derivatives */

      d_level_set_d_zeta=std::min(std::min(absolute_d_level_set_d_x1, absolute_d_level_set_d_x2), absolute_d_level_set_d_x3);

      /* d_level_set_d_eta is the remaining, not yet assigned first derivative */

      d_level_set_d_eta= absolute_d_level_set_d_x1+absolute_d_level_set_d_x2+absolute_d_level_set_d_x3-
                    d_level_set_d_xi - d_level_set_d_zeta;

      /* compute the value of the linearized level_set field in the vertices of the control volume */

      std::cerr<<"d_level_set_d_xi : "<< d_level_set_d_xi<<" \n";
      std::cerr<<"d_level_set_d_eta : "<< d_level_set_d_eta<<" \n";
      std::cerr<<"d_level_set_d_zeta : "<< d_level_set_d_zeta<<" \n";
      std::cerr<<"lower_bound_derivatives : "<< lower_bound_derivatives<<" \n";

      level_set_value_vertex_A = std::max (level_set - 0.5* ( (-d_level_set_d_xi) - d_level_set_d_eta - d_level_set_d_zeta), 0.0);
      level_set_value_vertex_B = std::max (level_set - 0.5* ( (-d_level_set_d_xi) - d_level_set_d_eta + d_level_set_d_zeta), 0.0);
      level_set_value_vertex_C = std::max (level_set - 0.5* ( (-d_level_set_d_xi) + d_level_set_d_eta - d_level_set_d_zeta), 0.0);
      level_set_value_vertex_D = std::max (level_set - 0.5* (  d_level_set_d_xi - d_level_set_d_eta - d_level_set_d_zeta), 0.0);
      level_set_value_vertex_E = std::max (level_set + 0.5* (  d_level_set_d_xi - d_level_set_d_eta - d_level_set_d_zeta), 0.0);

      level_set_value_vertex_A_signed=level_set - 0.5* ( (-d_level_set_d_xi) - d_level_set_d_eta - d_level_set_d_zeta);
      level_set_value_vertex_B_signed=level_set - 0.5* ( (-d_level_set_d_xi) - d_level_set_d_eta + d_level_set_d_zeta);
      level_set_value_vertex_C_signed=level_set - 0.5* ( (-d_level_set_d_xi) + d_level_set_d_eta - d_level_set_d_zeta);
      level_set_value_vertex_D_signed=level_set - 0.5* (  d_level_set_d_xi - d_level_set_d_eta - d_level_set_d_zeta);
      level_set_value_vertex_E_signed=level_set + 0.5* (  d_level_set_d_xi - d_level_set_d_eta - d_level_set_d_zeta);
      level_set_value_vertex_F_signed=level_set + 0.5* ( (-d_level_set_d_xi) + (-d_level_set_d_eta) + d_level_set_d_zeta);

      /* limit cases where one or more of the derivatives in the linearisation vanish have to be taken into account */

      if( d_level_set_d_xi > lower_bound_derivatives )
      {
           if( d_level_set_d_eta > lower_bound_derivatives )
           {
                if( d_level_set_d_zeta > lower_bound_derivatives )
                {
                  /* this is the standard case, where all derivatives are nonzero */
                  volume_of_fluid= ( level_set_value_vertex_A*level_set_value_vertex_A*level_set_value_vertex_A-
                                    level_set_value_vertex_B*level_set_value_vertex_B*level_set_value_vertex_B-
                                      level_set_value_vertex_C*level_set_value_vertex_C*level_set_value_vertex_C-
                                        level_set_value_vertex_D*level_set_value_vertex_D*level_set_value_vertex_D+
                                          level_set_value_vertex_E*level_set_value_vertex_E*level_set_value_vertex_E)/
                                          (6.0*d_level_set_d_xi*d_level_set_d_eta*d_level_set_d_zeta);
                    std::cerr<<"case 0\n";
                }
                else
                {
                /* limiting case of vanishing d_level_set_d_zeta */
                     if( level_set_value_vertex_A_signed>0 && level_set_value_vertex_C_signed>0)
                     {
                     volume_of_fluid=0.5+level_set*sqrt(1.0+(d_level_set_d_eta/d_level_set_d_xi)*(d_level_set_d_eta/d_level_set_d_xi));

                     std::cerr<<"case 1\n";
                     }
                     else
                     {
                       volume_of_fluid=0.5*(fabs(level_set_value_vertex_C_signed)/fabs(level_set_value_vertex_C_signed)+
                                                                                         fabs(level_set_value_vertex_F_signed))*
                                            (fabs(level_set_value_vertex_A_signed)/(fabs(level_set_value_vertex_A_signed)+fabs(level_set_value_vertex_D_signed)));
                      std::cerr<<"case 2\n";
                    }
 /*                    if( level_set_value_vertex_A>0 && level_set_value_vertex_C>0)
                     {
                     volume_of_fluid=level_set/d_level_set_d_xi+0.5;
                     }
                     else
                     {
 */ /*                  volume_of_fluid= ( level_set_value_vertex_A*level_set_value_vertex_A-
                                       level_set_value_vertex_C*level_set_value_vertex_C)/
                                       (2.0*d_level_set_d_xi*d_level_set_d_eta);
//                      }*/
           }

           }
           else
           {
                /* limiting case of vanishing d_level_set_d_zeta AND d_level_set_d_eta */

                volume_of_fluid= level_set_value_vertex_A/d_level_set_d_xi;
                    std::cerr<<"case 4\n";

           }
      }
      else
      {

              /* limiting case of ALL partial derivatives vanishing */
              if(level_set <-1E-07)
              {
                /* set the volume of fluid to zero */
                  volume_of_fluid=0.0;
              }
              else
              {
                /* set the volume of fluid to 0.5 */

                  volume_of_fluid=0.5;
              }
      }
      return 0;
      }


/********************************************************************************/
/********************************************************************************/
/*  Function to compute the volume of fluid field from a given level-set field  */
/*                                                                    */
/*  Programmer       : Duncan van der Heul                                          */
/*  Date      : 10-03-2013                                            */
/*  Update    :                                                       */
/********************************************************************************/
/* Notes                                                              */
/* For the case of all first partial derivatives nonzero, the standard          */
/* formulation is used. Otherwise a number of limiting cases are considered     */
/* The algorithm is discussed in paragraph 5.4 of Sander's thesis            */
/********************************************************************************/
      int level_set_2_vof(
             double level_set,            // level-set field value in this cell
             double d_level_set_d_x1,            // first partial derivative in x1
                                          // direction of level-set
             double d_level_set_d_x2,            // first partial derivative in x2
                                          // direction of level-set
             double d_level_set_d_x3,            // first partial derivative in x3
                                          // direction of level-set
             double &volume_of_fluid,            // volume of fluid value in this cell
             double lower_bound_derivatives    // lower bound for the first partial derivatives
                                          // to consider it a limiting case of vanishing
                                          // partial derivatives

      )
{
       double complement_volume_of_fluid; // (volume of fluid)-1 for -1* level set
           /* function definitions */

      int level_set_2_vof_phi_negative(   // compute volume of fluid from level-set
             double level_set,            // field, when level-set field is less than zero.
             double d_level_set_d_x1,

             double d_level_set_d_x2,

             double d_level_set_d_x3,

             double &volume_of_fluid,
             double lower_bound_derivatives
             );

      if(level_set<=0.0)
      {
         if(level_set_2_vof_phi_negative(level_set,
                d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
                       volume_of_fluid, lower_bound_derivatives))
         {
           std::cerr << "level-set to vof conversion failed in call from level_set_2_vof. \n";
           std::cerr << "for level_set< 0 \n";
           exit(1);
           return 1;
         }
      }
      else
      {
         if(!level_set_2_vof_phi_negative(-1.0*level_set,
                d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
                       complement_volume_of_fluid, lower_bound_derivatives))
         {
             volume_of_fluid=1.0-complement_volume_of_fluid;
         }
         else
         {
           std::cerr << "level-set to vof conversion failed in call from level_set_2_vof. \n";
           std::cerr << "for level_set >0 \n";
           exit(1);
           return 1;
         }
      }
      if(volume_of_fluid >1.01)
      {
         std::cerr <<"vof too large \n";
         exit(1);
                return 1;
      }

      return 0;
}
      
/* program to test formulations of the function level-set to volume of fluid */


int main()
{
           double level_set_left;                // left hand initial value for root finding
           double level_set_right;               // right hand initial value for root finding
           double volume_of_fluid;               // volume of fluid value for which we need the
           double d_level_set_d_x1;              // first partial derivative in x1
                                                 // direction of level-set
           double d_level_set_d_x2;                     // first partial derivative in x2
                                                 // direction of level-set
           double d_level_set_d_x3;                     // first partial derivative in x3
                                                 // direction of level-set
           double lower_bound_derivatives;           // lower bound for the first partial derivatives
                                                 // to consider it a limiting case of vanishing
                                                 // partial derivatives
           double level_set_punt;
           double functie_punt;
           double new_volume_of_fluid;

           int i;

           ofstream problem_function_file( "problem_function.m");

           problem_function_file << "clear all ;\n";
           problem_function_file << "function_g=[ \n";

           volume_of_fluid=0.5;
// function_g_left : -0.167875
// function_g_right: 0.832125
// level_set_left  : -0.000272041
// level_set_right : 0.000272041
// d_level_set_d_x1: 2.16768e-08
// d_level_set_d_x2: 5.58434e-09
// d_level_set_d_x3: -0.000544055
           d_level_set_d_x1= 2.16768e-08;
           d_level_set_d_x2= 5.58434e-09;
           d_level_set_d_x3= -0.000544055;
           d_level_set_d_x1= 1e-12;
           d_level_set_d_x2= 0.1;
           d_level_set_d_x3= 0.2;
           lower_bound_derivatives=0.000000001;
      level_set_left= -0.5*(fabs(d_level_set_d_x1)+
                            fabs(d_level_set_d_x2)+
                              fabs(d_level_set_d_x3));
      level_set_right=-1.0*level_set_left;

                if(!problem_function_file)
                {
                    /* the contructor returned a 0-pointer :-( */
                    cout << "Cannot open file.\n";
                    cout << "In function bisection_method \n";
                    cout << "Line 168 \n";
                    exit(1);
                }
               for(i=1;i<100;i++)
               {
                    level_set_punt=(level_set_right-level_set_left)/100*i+level_set_left;
                    if(!level_set_2_vof(level_set_punt,
                      d_level_set_d_x1, d_level_set_d_x2, d_level_set_d_x3,
                        new_volume_of_fluid, lower_bound_derivatives))
                    {
                    functie_punt=new_volume_of_fluid-volume_of_fluid;
                    }
                    else
                    {
                           cout<<"problem in evaluation of level_set_2_vof \n";
                           cout << "In function bisection_method \n";
                           cout << "Line 174 \n";
                    }
                    problem_function_file<<level_set_punt<<" "<<functie_punt<<";"<<"\n";
               }

             problem_function_file <<"]; \n";
             problem_function_file <<"clf; \n";
             problem_function_file <<" plot(function_g(:,1), function_g(:,2)); \n";

}

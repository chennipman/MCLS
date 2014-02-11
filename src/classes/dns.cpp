/********************************************************************************/
/********************************************************************************/
/*  Main program for the Cartesian Mass Conserving Level Set Method         	*/
/*  Programmer	: Duncan van der Heul       					*/
/*  Date	: 03-01       							*/
/*  Update	:        							*/
/********************************************************************************/
//
//
/********************************************************************************/
#include "header_main.h" 
/********************************************************************************/
int main ()
{
//
    
    std::cout << "bubble " << gnoef[1].center_location.x2 << "\n";
    
    
    boundary_face boundary_faces[6];
    
    std::cout << boundary_faces[0].boundary_variables[0].boundary_condition_value;

/* intialization */
   
//     initialize_computation();
    
    std::cout << "Initialization completed \n";

/* time stepping sequence */

    std::cout << "Time stepping completed \n";

/* post processing */

    std::cout << "Post processing completed \n";

    return 0;
}













/********************************************************************************/
/********************************************************************************/
/*  Function to compute the WENO flux for a cell centered field                 */
/*                                                                              */
/*                                                                              */
/*  Programmer  : Duncan van der Heul                                           */
/*  Date        : 10-03-2013                                                    */
/*  Update      :                                                               */
/********************************************************************************/
/* Notes:                                                                       */
/* The notation of Shu is adopted for all the terms that comprise the flux.     */
/********************************************************************************/
  double weno_flux_computation(
        double field_value_minus_3,             // field at i-3
        double field_value_minus_2,             // field at i-2
        double field_value_minus_1,             // field at i-1
        double field_value_0,                   // field at i
        double field_value_plus_1,              // field at i+1
        double field_value_plus_2,              // field at i+2
        double mesh_width                       // mesh width
               )
{
       
    double Flux_one_minus_2;            // first candidate first order ENO
    double Flux_one_minus_1;            // second candidate first order ENO
    double Flux_one_0;                  // third candidate first order ENO
    double Flux_one_plus_1;             // fourth candidate first order ENO
    double Flux_one_plus_2;             // fifth candidate first order ENO
    double Flux_three_0;                // first candidate third order ENO
    double Flux_three_1;                // second candidate third order ENO
    double Flux_three_2;                // third candidate third order ENO
    double  g_0=0.1, g_1=0.6, g_2=0.3;  // linear weights
    double  b_0, b_1, b_2;              // smoothness indicators
    double  w_bar_0, w_bar_1, w_bar_2;  // WENO weights (notation of Shu)
    double  w_0, w_1, w_2;              // WENO weights (notation of Shu)
    double eps=10E-6;                   // small number to avoid division by zero
                                        // value taken from Shu paper.
    double  weno_flux;                  // flux based on WENO interpolation
    
    /* define the first order ENO candidates */
    
//     v_im2 = (f2 - f1) / dx1
//     v_im1 = (f3 - f2) / dx2
//     v_i   = (f4 - f3) / dx3
//     v_ip1 = (f5 - f4) / dx4
//     v_ip2 = (f6 - f5) / dx5

    Flux_one_minus_2= (field_value_minus_2-field_value_minus_3)/mesh_width;
    Flux_one_minus_1= (field_value_minus_1-field_value_minus_2)/mesh_width;
    Flux_one_0      = (field_value_0      -field_value_minus_1)/mesh_width;
    Flux_one_plus_1 = (field_value_plus_1 -field_value_0      )/mesh_width;
    Flux_one_plus_2 = (field_value_plus_2 -field_value_plus_1 )/mesh_width;
    
    /* define third order ENO candidates  */
    
    Flux_three_0= 1./3.* Flux_one_minus_2 - 7./6.* Flux_one_minus_1 + 11./6.* Flux_one_0;
    Flux_three_1=-1./6.* Flux_one_minus_1 + 5./6.* Flux_one_0       +  1./3.* Flux_one_plus_1;
    Flux_three_2= 1./3.* Flux_one_0       + 5./6.* Flux_one_plus_1  -  1./6.* Flux_one_plus_2;
    
//     p0 = v_im2/3. - 7.*v_im1/6. + 11.*v_i/6.
//     p1 = -v_im1/6. + 5.*v_i/6. + v_ip1/3.
//     p2 = v_i/3. + 5.*v_ip1/6. - v_ip2/6.
    
    
    /*  Smoothness indicators */
    
    
    
//     b0 = 13.*(v_im2 - 2.*v_im1 + v_i)**2. / 12. + (3.*v_im2 - 4.*v_im1 + v_i)**2. / 4.
//     b1 = 13.*(v_im1 - 2.*v_i + v_ip1)**2. / 12. + (3.*v_im1 - v_ip1)**2. / 4.
//     b2 = 13.*(v_i - 2.*v_ip1 + v_ip2)**2. / 12. + (v_i - 4.*v_ip1 + v_ip2)**2. / 4.
    
//     b_0 = 13./12*(Flux_one_minus_2 - 2.*Flux_one_minus_1 + Flux_one_0)*(Flux_one_minus_2 - 2.*Flux_one_minus_1 + Flux_one_0)+ 
//                         1./4.*(3.*Flux_one_minus_2 - 4.*Flux_one_minus_1 + Flux_one_0)*(3.*Flux_one_minus_2 - 4.*Flux_one_minus_1 + Flux_one_0);
//     b_1 = 13./12*(Flux_one_minus_1 - 2.*Flux_one_0 + Flux_one_plus_1)*(Flux_one_minus_1 - 2.*Flux_one_0 + Flux_one_plus_1)  + 
//                         1./4.*(3.*Flux_one_minus_1 - Flux_one_plus_1)*(3.*Flux_one_minus_1 - Flux_one_plus_1);
//     b_2 = 13./12*(Flux_one_0 - 2.*Flux_one_plus_1 + Flux_one_plus_2)*(Flux_one_0 - 2.*Flux_one_plus_1 + Flux_one_plus_2) +
//                         1./4.*(Flux_one_0 - 4.*Flux_one_plus_1 + Flux_one_plus_2)*(Flux_one_0 - 4.*Flux_one_plus_1 + Flux_one_plus_2);
    b_0 = 13./12*(Flux_one_minus_2 - 2.*Flux_one_minus_1 + Flux_one_0)*(Flux_one_minus_2 - 2.*Flux_one_minus_1 + Flux_one_0)+ 
                        1./4.*(Flux_one_minus_2 - 4.*Flux_one_minus_1 + 3.*Flux_one_0)*(Flux_one_minus_2 - 4.*Flux_one_minus_1 + 3.*Flux_one_0);
    b_1 = 13./12*(Flux_one_minus_1 - 2.*Flux_one_0 + Flux_one_plus_1)*(Flux_one_minus_1 - 2.*Flux_one_0 + Flux_one_plus_1)  + 
                        1./4.*(Flux_one_minus_1 - Flux_one_plus_1)*(Flux_one_minus_1 - Flux_one_plus_1);
    b_2 = 13./12*(Flux_one_0 - 2.*Flux_one_plus_1 + Flux_one_plus_2)*(Flux_one_0 - 2.*Flux_one_plus_1 + Flux_one_plus_2) +
                        1./4.*(3*Flux_one_0 - 4.*Flux_one_plus_1 + Flux_one_plus_2)*(3*Flux_one_0 - 4.*Flux_one_plus_1 + Flux_one_plus_2);
    
    /* Nonlinear weights */
    /* eps is added to the denominator to avoid division by zero */
    /* the value is taken from the paper by Shu  */
    
    eps = 1.e-6;         
    
    w_bar_0 = g_0 / ( (eps + b_0)*(eps + b_0));
    w_bar_1 = g_1 / ( (eps + b_1)*(eps + b_1));
    w_bar_2 = g_2 / ( (eps + b_2)*(eps + b_2));
    
    w_0 = w_bar_0 / (w_bar_0 + w_bar_1 + w_bar_2);
    w_1 = w_bar_1 / (w_bar_0 + w_bar_1 + w_bar_2);
    w_2 = w_bar_2 / (w_bar_0 + w_bar_1 + w_bar_2);
    
    /* the flux */
    
    weno_flux = w_0*Flux_three_0 + w_1*Flux_three_1 + w_2*Flux_three_2;
    
    
    return weno_flux ;

}


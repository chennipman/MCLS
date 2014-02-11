
    ! First order derivatives.

    /* define the first order ENO stencils */
    
//     v_im2 = (f2 - f1) / dx1
//     v_im1 = (f3 - f2) / dx2
//     v_i   = (f4 - f3) / dx3
//     v_ip1 = (f5 - f4) / dx4
//     v_ip2 = (f6 - f5) / dx5

    Flux_one_minus_2= (field_value_minus_2-field_value_minus_3)/mesh_width;
    Flux_one_minus_1= (field_value_minus_1-field_value_minus_2)/mesh_width;
    Flux_one_0      = (field_value_0      -field_value_plus_1 )/mesh_width;
    Flux_one_plus_1 = (field_value_plus_1 -field_value_plus_0 )/mesh_width;
    Flux_one_plus_2 = (field_value_plus_2 -field_value_plus_1 )/mesh_width;
    
    /* define third order ENO stencils  */
    
    Flux_three_0= 1/3* Flux_one_minus_2 - 7/6* Flux_one_minus_1 + 11/6* Flux_one_0;
    Flux_three_1=-1/6* Flux_one_minus_1 + 5/6* Flux_one_0       +  1/3* Flux_one_plus_1;
    Flux_three_2= 1/3* Flux_one_0       + 5/6* Flux_one_plus_1  -  1/6* Flux_one_plus_2;
    
//     p0 = v_im2/3. - 7.*v_im1/6. + 11.*v_i/6.
//     p1 = -v_im1/6. + 5.*v_i/6. + v_ip1/3.
//     p2 = v_i/3. + 5.*v_ip1/6. - v_ip2/6.
    
    /* Linear weights, again following notation of Shu */
    
    g0 = 0.1;
    g1 = 0.6;
    g2 = 0.3;
    
    /*  Smoothness indicators */
    
    
    
//     b0 = 13.*(v_im2 - 2.*v_im1 + v_i)**2. / 12. + (3.*v_im2 - 4.*v_im1 + v_i)**2. / 4.
//     b1 = 13.*(v_im1 - 2.*v_i + v_ip1)**2. / 12. + (3.*v_im1 - v_ip1)**2. / 4.
//     b2 = 13.*(v_i - 2.*v_ip1 + v_ip2)**2. / 12. + (v_i - 4.*v_ip1 + v_ip2)**2. / 4.
    
    b_0 = 13./12*(Flux_one_minus_2 - 2.*Flux_one_minus_1 + Flux_one_0)*(Flux_one_minus_2 - 2.*Flux_one_minus_1 + Flux_one_0)+ 
                        1./4.*(3.*Flux_one_minus_2 - 4.*Flux_one_minus_1 + Flux_one_0)*(3.*Flux_one_minus_2 - 4.*Flux_one_minus_1 + Flux_one_0);
    b_1 = 13./12*(Flux_one_minus_1 - 2.*Flux_one_0 + Flux_one_plus_1)*(Flux_one_minus_1 - 2.*Flux_one_0 + Flux_one_plus_1)  + 
                        1./4.*(3.*Flux_one_minus_1 - Flux_one_plus_1)*(3.*Flux_one_minus_1 - Flux_one_plus_1);
    b_2 = 13./12*(Flux_one_0 - 2.*Flux_one_plus_1 + Flux_one_plus_2)*(Flux_one_0 - 2.*Flux_one_plus_1 + Flux_one_plus_2) +
                        1./4.*(Flux_one_0 - 4.*Flux_one_plus_1 + Flux_one_plus_2)*(Flux_one_0 - 4.*Flux_one_plus_1 + Flux_one_plus_2);
    
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
    
    !-------------------------------------------------------------------------------

    end function weno_flux


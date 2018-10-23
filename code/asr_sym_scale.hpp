#ifndef ASRSYMSC_H
#define ASRSYMSC_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include<cmath>
#include <complex>
#include "getHtotOMEN.hpp"



/* ==========================================================================

                    FUNCTION  :  asr  (Acoustic Sum Rule) 
                                 symmetrizing
                                 scale

=============================================================================
in this function, we
1)impose the acoustic sum rule 
and 
2)symmetrizing the different parts of the force-constant matrix 
3) Scale the elements of force constant matrices with the corresponding masses 
    --> dynamical matrix.

==============================================================================
1) ASR

1a) getHtotOMEN  at k = (0.0,0.0,0.0)
    this function generates the total Hamiltonian HtotOMEN (3*Na_uc_new, 3*Na_uc_new)  
    at a given k-point out of the tensor Hs_mod 
    (NFold_bs_x,3*Na_uc_new,NFold_bs_y*NFold_bs_z*3*Na_uc_new)
    as it is done in OMEN.

    Note that HtotOMEN will be used to calculate the new phonon dispersion

    input: Hs_mod, NFold_bs, k, Na_uc_new
    output: HtotOMEN 

1b) the Acoustic Sum Rule is imposed

    C(I,J,a,b) = IFC of (atom I direction a) & (atom J direction b)
    I,J = .., Na_uc_new
    a,b = x,y,z
    
    ASR : C(I,I,a,b) = - SUM_{J != I} C(I,J,a,b)
    
*/

Eigen::Tensor<double, 3> asr_sym_scale(Eigen::Tensor<double, 3> Hs_mod_input,
                                       Eigen::MatrixXd NFold_bs,int Na_uc_new, 
                                       Eigen::MatrixXd LMs_uc_new, Eigen::MatrixXd mr);


#endif

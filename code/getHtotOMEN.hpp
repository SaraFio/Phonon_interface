#ifndef GETHOMEN_H
#define GETHOMEN_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include<cmath>
#include <complex>

#include"find_vec.hpp"

/* ==========================================================================

                    FUNCTION  :  getHtotOMEN

=============================================================================

this function generates the total Hamiltonian HtotOMEN (3*Na_uc_new, 3*Na_uc_new)  
at a given k-point out of the tensor Hs_mod 
(NFold_bs_x,3*Na_uc_new,NFold_bs_y*NFold_bs_z*3*Na_uc_new)
 as it is done in OMEN.

input: Hs_mod, NFold_bs, k, Na_uc_new
output: HtotOMEN 
*/

Eigen::MatrixXcd getHtotOMEN(Eigen::Tensor<double, 3> Hs_mod,Eigen::MatrixXd NFold_bs,
                             Eigen::MatrixXd k,int Na_uc_new);


#endif


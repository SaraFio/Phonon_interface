#ifndef EXTSUB_H
#define EXTSUB_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include"find_vec.hpp"

/* ==========================================================================

                    FUNCTION  :  extract_submatrices

=============================================================================

this function calculates  matrices in OMEN style, 
using the updated structure information from check_slabs function 
and extend uc section 
(new unit cell = Na_uc_new, Lxyz_uc_new, LMs_uc_new)

Note that we are studying transport properties -> 
at least trasport direction (x) must be repeated = NFold_bs(0)=3

input: FCs_cut_sym, LMs_ord, LMs_uc_new, Lxyz,Lxyz_uc_new, NFold_bs

output : Hs_mod Tensor (NFold_bs(1), NFold_bs(2)*3*Na_uc_new, NFold_bs(1)*3*Na_uc_new) 
         it is composed by
         NFold_bs(0) matrices = 3 matrices: H00,H10,H01 
         matrices dimensions: ( NFold_bs(2)*3*Na_uc_new, NFold_bs(1)*3*Na_uc_new)
         elements : FCs_cut_sym

Note that:
1) 
NFold_bs=[3,3,3] = 3D -> 3 matrices of dim (3*Na_uc_new,3*3*3*Na_uc_new)
NFold_bs=[3,1,3] = 2D -> 3 matrices of dim (3*Na_uc_new,1*3*3*Na_uc_new)
NFold_bs=[3,1,1] = 1D -> 3 matrices of dim (3*Na_uc_new,1*1*3*Na_uc_new)

2)
to fill the matrices, artificially repeated new unit cells are created:
NFold_bs(0) new unit cells along x
NFold_bs(1) new unit cells along y
NFold_bs(2) new unit cells along z
Clearly the numbers of atoms artificially created is > Na
(Remember we only have FCs_cut_sym elements for the original
Na atoms, thus not all the artificially generated atoms 
will be used.)
*/

Eigen::Tensor<double, 3> extract_submatrices(Eigen::MatrixXd FCs_cut_sym,int Na, int Na_uc_new,
                                      Eigen::MatrixXd LMs_ord, Eigen::MatrixXd LMs_uc_new, 
                                             Eigen::MatrixXd Lxyz,Eigen::MatrixXd Lxyz_uc_new,Eigen::MatrixXd NFold_bs);


#endif

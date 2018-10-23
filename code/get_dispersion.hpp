#ifndef GETDISP_H
#define GETDISP_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include<string>
#include<cmath>
#include <complex>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>
#include<algorithm>

#include"getHtotOMEN.hpp"


using Eigen::MatrixXd;
/* 
 #############################################################################

                      GET PHONON BAND STRUCTURE

 #############################################################################

The phonon band structure is calculated using the Hs_mod tensor.



-    for each k-point ka(i) a matrix Hamiltonian is created H(ka(i))
     by the function getHtotOMEN .

-    the corresponding eigenvalues E(ka(i)) are calculated using the 
     Eigen attribute .eigenvalues()
     all the  E(ka(i)) are stored in matrix E -> dim(E)=(3*Na_uc_new,ntotk)
     where ntotk = (Nkx+1)+(Nky+1)+(Nkz+1)
     Note that E(ka(i)) is a vector of length 3*Na_uc_new. 
     The first Nk(0)+1 columns of E contain E(kx), 
     the second Nk(1)+1 columns contain E(ky),
     the last Nk(2)+1 columns contain E(kz)
     
-    in order to obtain the phonon energies 
     E(i,j)=sign*hq*sqrt(sign*E(i,j));
     where hq = hbar/charge;
           sign = -1 if E(i,j)<0
                   1 otherwise

     The E matrix columns are returned in ascending order

*/ 


Eigen::MatrixXd get_dispersion(Eigen::Tensor<double, 3> Hs_mod, Eigen::MatrixXd NFold_bs,
                               Eigen::MatrixXd kx, Eigen::MatrixXd ky, Eigen::MatrixXd kz,
                               int indE, int Na_uc, double hbar, double charge);


#endif

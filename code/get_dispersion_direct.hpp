#ifndef GETDISPDIR_H
#define GETDISPDIR_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include<string>
#include<cmath>
#include <complex>
#include <Eigen/Eigenvalues>
#include<algorithm>

#include"getHtot.hpp"
#include"find_vec.hpp"



using Eigen::MatrixXd;
/* 
 #############################################################################

                      GET PHONON BAND STRUCTURE

 #############################################################################

The phonon band structure is calculated using the Force Constant matrix.

1st) the k-space grid of a rectangular brillouin zone is generated
     by the function get_k(Na). 
     Input:   Nka: number k-points in ka direction (int)  (a = x,y,z)
     Output: ka: vector containing the Nka ka-points between -pi and pi     
     nb. actually ka vector is a  (Nka,1) matrix

2nd) for each k-point ka(i) a matrix Hamiltonian is created H(ka(i))
     by the function getHtot.

3rd) the corresponding eigenvalues E(ka(i)) are calculated using the 
     Eigen attribute .eigenvalues()
     all the  E(ka(i)) are stored in matrix E -> dim(E)=(3*Na_uc,ntotk)
     where ntotk = (Nkx+1)+(Nky+1)+(Nkz+1)
     Note that E(ka(i)) is a vector of length 3*Na_uc. 
     The first Nk(0)+1 columns of E contain E(kx), 
     the second Nk(1)+1 columns contain E(ky),
     the last Nk(2)+1 columns contain E(kz)
     
4th) in order to obtain the phonon energies 
     E(i,j)=sign*hq*sqrt(sign*E(i,j));
     where hq = hbar/charge;
           sign = -1 if E(i,j)<0
                   1 otherwise

     The E matrix columns are returned in ascending order

*/ 


Eigen::MatrixXd get_dispersion_direct(Eigen::MatrixXd FCs_m, Eigen::MatrixXd NFold,
                                      Eigen::MatrixXd  kx, Eigen::MatrixXd  ky, Eigen::MatrixXd kz,
                                      int indE, int Na, int Na_uc,
                                      Eigen::MatrixXd LMs_ord,Eigen::MatrixXd LMs_uc_ord, Eigen::MatrixXd Lxyz_uc,
                                      Eigen::MatrixXd rvec_short, Eigen::MatrixXd Nmulti,double comp_diff, 
                                      double hbar, double charge);


#endif

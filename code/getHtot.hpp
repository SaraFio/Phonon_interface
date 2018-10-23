#ifndef GETHTOT_H
#define GETHTOT_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include<string>
#include<cmath>
#include <complex>

#include"find_vec.hpp"


/* ==========================================================================

              FUNCTION :  getHtot

=============================================================================

this function generates the  total Hamiltonian Htot out of the force constant matrix 
as it is done in PHONOPY at a given k-point

*/


/*entries:
FCs_m, NFold, kx, ky, kz, LMs_ord, LMs_uc_ord, Lxyz_uc, rvec_short, Nmulti, num_SCELL,comp_diff*/

Eigen::MatrixXcd getHtot(Eigen::MatrixXd FCs_m, Eigen::MatrixXd NFold, double kx, double  ky,double kz,
                          Eigen::MatrixXd LM_coord,Eigen::MatrixXd LM_uc_coord, Eigen::MatrixXd Lxyz_uc,
                         Eigen::MatrixXd rvec_short, Eigen::MatrixXd Nmulti, int num_SCELL,double comp_diff);

#endif

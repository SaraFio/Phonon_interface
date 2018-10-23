#ifndef EXTUC_H
#define EXTUC_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

#include"find_vec.hpp"

/* ==========================================================================

                    FUNCTION  :  check_slab

=============================================================================

this function checks how many neighbor slabs in x,y,z directions need to be included
to fully cover the interaction in the FCs_cut_sym

input: FCs_cut_sym, LMs_ord, LMs_uc_ord, Lxyz_uc, NFold, comp_diff

output : (3,1)matrix nneigh:(num. slabs in x-dir, num. slabs in y-dir, num. slabs in z-dir, )
*/

Eigen::MatrixXd check_slabs(Eigen::MatrixXd FCs_cut_sym, Eigen::MatrixXd LMs_ord,
                            Eigen::MatrixXd LMs_uc_ord, Eigen::MatrixXd Lxyz_uc,
                            Eigen::MatrixXd NFold,double comp_diff );

#endif

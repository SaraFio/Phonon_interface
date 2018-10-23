#ifndef FINDNN_H
#define FINDNN_H

#include <fstream>    
#include <iomanip>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include<string>
#include<cmath>


/* ==========================================================================

              C L A S S :  FindNN

=============================================================================

this class looks for the closest neighbors in the reference supercell as well 
as in the artificially generated neighboring supercells. 
The information is needed to preserve the symmetry properties of the 
bandstructure of the  unit-cell
*/

class FindNN{
 private: 
  int Na_; // number of the atoms in the supercell (from sc_file)
  int Nauc_; // number of the atoms in the unitcell (from uc_file)
  Eigen::MatrixXd LMs_ord_; //sorted LM (from sc_file)
  Eigen::MatrixXd LMs_uc_ord_; //sorted LM (from uc_file)
  Eigen::MatrixXd Lxyz_; //x, y, and z extension of the supercell. (from ms_file)
  Eigen::MatrixXd Lxyzuc_; //x, y, and z extension of the unitcell.(from ms_file)
  Eigen::MatrixXd ixyz1_ref_; /* lower left corner of the reference unit-cell.
                                dim (1,3) = (x_ref, y_ref, z_ref)*/
  Eigen::MatrixXd NSCELL_; /* number of neighboring supercells in x, y, and z direction.
                            dim(1,3) = (NSCELL_x,NSCELL_y,NSCELL_z)*/
  double  compdiff_;// numerical criteria.
  Eigen::MatrixXd rvec_short_;/*Containing the bond vectors  and the length of the bond vector 
                                between  an atom in the  unit-cell  and an atom in a supercell.
                                dim (Na_uc * Na * (NSCELL_x*2+1) * (NSCELL_y*2+1) * (NSCELL_z*2+1),11) =
                                (atom in uc,  atom in sc, NSCELL_x , NSCELL_y , NSCELL_z,n_SCELL ,x12, y12, z12, r12, nn_logical)
                                where NSCELL_i = coordinates of the supercell  (ex. 1,-1,0)
                                      n_SCELL = supercell index (ex 23)
                                      x12 = x2 - x1, y12 = y2 - y1, z12 = z2 - z1
                                      r12 = sqrt(x12*x12 + y12*y12 + z12*z12)
          nn_logical is 1 if ( r12-(min dist 1-2) ) < comp_diff, otherwise 0*/
 Eigen::MatrixXd Nmulti_; /* Counts the closest neighbors considering all neighboring supercells
                          Dimension (Na_uc,Na).
                          Is obtained summing up all column 10 of rvec_short, for fixed atmo1 and atom2*/

   public:
 FindNN(int, int,Eigen::MatrixXd ,Eigen::MatrixXd,Eigen::MatrixXd ,Eigen::MatrixXd,
Eigen::MatrixXd ,Eigen::MatrixXd, double); //constructor
 
 void find_nn(); //nearest neighbors are found and stored.
 
 Eigen::MatrixXd rvec_short(){return rvec_short_;}

 Eigen::MatrixXd Nmulti(){return Nmulti_;}

};
#endif

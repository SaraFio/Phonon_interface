#ifndef MATERIALSTRUCTUREFILE_H
#define MATERIALSTRUCTUREFILE_H

#include "InputFile.hpp"

#include <fstream>    
#include <iomanip>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include<string>


/* ==========================================================================

                    CLASS  :  MaterialStructureFile

=============================================================================

this class reads the informations in the ms_file and store them

*/


class MaterialStructureFile{
private:
  int nas_; // number of atomic species
  std::string  filename_matstr; // name of the file with material parameters
  std::string dummy_matstr_; /* this string is used to read indication before the data
                         it's a dummy variable. not used.  */  
  double  hbar_; // h_bar in Joule second
  double charge_; //charge in Coulomb
  double mscalefac_; // mass scale factor
  double tollim_;  //numerical criteria  for tolerance
  Eigen::MatrixXd mr_;  //reduced masses of atomic species matrix(1,n_as)
  Eigen::MatrixXd mt_;  /*total masses of atomic species matrix(1,n_as) 
                          = mr*mscalefac      */
  Eigen::MatrixXd NFold_;     /*Defines the supercell size based on the basic unit-cell.
                                [x-dir,y-dir,z-dir]. (Rectangular unit-cell!) */
  Eigen::MatrixXd NFoldbs_;   /* Periodically extended directions get a 3, confined
                                 directions get a 1. [x-dir,y-dir,z-dir],eg. wire: [3 1 1]*/


  Eigen::MatrixXd ixyz1_ref_; /*Lower left corner of the reference unitcell
                              dim(ixyz1_ref) = (1,3)*/ 

  Eigen::MatrixXd NSCELL_; /* number of artificially genereated repeated supercells
                                along the x,y,z direction
                              dim(NSCELL) = (1,3)*/ 

  Eigen::MatrixXd Nk_; /* Number k-points in ka direction (a=x,y,z)
                          Nk = [Nkx, Nky, Nkz]  dim(NSCELL) = (1,3)*/ 

  double rc_; /*radius cutoff (used to cut the FC) Angstrom*/


public:
  MaterialStructureFile(int); // constructor
  void read(std::string); // read the material structure file
  double hbar() {return hbar_;}
  double charge() {return charge_;}
  double m_scalefac() {return mscalefac_;}
  double tol_lim() {return tollim_;}
  Eigen::MatrixXd mr() {return mr_;}
  Eigen::MatrixXd mt() {return mt_;}
  Eigen::MatrixXd NFold() {return NFold_;}
  Eigen::MatrixXd NFold_bs() {return NFoldbs_;}

  Eigen::MatrixXd ixyz1_ref() {return ixyz1_ref_;}
  Eigen::MatrixXd NSCELL() {return NSCELL_;}
  Eigen::MatrixXd Nk() {return Nk_;}
  double rc() {return rc_;}
};


#endif

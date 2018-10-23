#ifndef POSCARFILE_H
#define POSCARFILE_H

#include "InputFile_tosort.hpp"

#include <fstream>    
#include <iomanip>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include<string>








/* ==========================================================================

              C L A S S :  POSCAR

=============================================================================

this class reads the data contained in the POSCAR  file

*/

class PoscarFile{
private:
  //int Nabuc_; // number of atom in the basic unit cell (this value is get from BucNamesFile.read())

  int nas_; // number of atomic species (this value is get from BucNamesFile.read())

  int Na_; // number of atom in the cell
  
  std::string compoundname_; // name of the compound under investigation

  double scalefac_; // scale factor

  Eigen::MatrixXd lattvec_; // matrix(3,3) of the coomponents of lattice vectors

  std::string atomicspecies_; // atomic species names (dummy variable, not used)

  Eigen::MatrixXd Naps_; // Matrix Naps(1,n_as) : Number of Atoms Per atomic Species

  //std::string coordunit_; 
  char coordunit_[40];// Direct or Cartesian coordinates

  Eigen::MatrixXd LM_onlycoord_; /*matrix(Na,3) : coordinates of the atoms
                         col of LM = (x-coord,y-coord,z-coord)*/

  Eigen::MatrixXd LM_; /*matrix(Na,4) : coordinates of the atoms + atomic species
                         col of LM = (x-coord,y-coord,z-coord, at. spec.)*/

  Eigen::MatrixXd LMs_; /*matrix(Na,4) :sorted according to OMEN convention
                          coordinates of the atoms + atomic species
                          col of LMs = (x-coord,y-coord,z-coord, at. spec.)*/

  Eigen::MatrixXd LMs_ord_; /*matrix(Na,5) :sorted according to OMEN convention
                          coordinates of the atoms + atomic species + 
                          correspondance original LM and LMs
                          col of LMs = (x-coord,y-coord,z-coord, at. spec., correspondance)
                          ex. if a coordinate in the original matrix LM is in position 3 and 
                          in the sorted matrix LMs is in position 9, we have LMs_ord(3,5) = 9*/

public:
  PoscarFile(int); // constructor

  void read(std::string); // read the POSCAR file

  std::string compound_name() {return compoundname_;} 
  double scale_fac() {return scalefac_;}
  Eigen::MatrixXd latt_vec(){return lattvec_;} 
  int Na() {return Na_;}
  Eigen::MatrixXd Naps() {return Naps_;}
  Eigen::MatrixXd LMcoord() {return LM_onlycoord_;}
  Eigen::MatrixXd LM() {return LM_;}

  void sort(); //sort LM matrix 

  Eigen::MatrixXd LMs() {return LMs_;}
  Eigen::MatrixXd LMs_ord() {return LMs_ord_;}
  
};



#endif

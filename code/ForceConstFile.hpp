#ifndef FCFILE_H
#define FCFILE_H

#include "InputFile_tosort.hpp"

#include <fstream>    
#include <iomanip>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include<string>
#include<cmath>

/* ==========================================================================

              C L A S S :  ForceConstFile

=============================================================================

this class 
- reads the data contained in the FC_file (output Phonopy)
  and write them into a matrix FC_eVA
- convert FC_eVA  from eV/A^2 -> J/m^2: FC
- sort FC  according to OMEN -> FC_s 
- scale FC_s scaled by the square root of the corresponding masses
  FC_s(i,j)/sqrt(m_i*m_j)

*/

class ForceConstFile{
private:
  int Na_; // number of the atoms in the supercell (written in the FC_file)
  int Na_sc_; // number of the atoms in the supercell (from sc_file)
  Eigen::MatrixXd LMs_ord_; //sorted LM (from sc_file)
  Eigen::MatrixXd mr_; //sorted LM (from ms_file)
  Eigen::MatrixXd FC_eVA_; // force constants matrix : (eV/A^2)
  Eigen::MatrixXd FC_; // force constants matrix : (J/m^2)
  Eigen::MatrixXd FCs_; // sorted force constants matrix (J/m^2)
  Eigen::MatrixXd FCs_m_; // sorted and mass scaled force constants matrix (J/m^2)
 public:
  ForceConstFile(int,Eigen::MatrixXd ,Eigen::MatrixXd); //constructor

  void read(std::string); //read the FC_file -> FC_eVA

  Eigen::MatrixXd FC_eVA() {return FC_eVA_;}

  void convert(); //convert from eV/A^2 -> J/m^2: FC

  Eigen::MatrixXd FC() {return FC_;}

  void sort(); //sort FC -> FC_s

  Eigen::MatrixXd FCs() {return FCs_;}

  void scale(); /*FCs scaled by the square root 
                  of the corresponding masses*/

  Eigen::MatrixXd FCs_m() {return FCs_m_;}
  
};

#endif

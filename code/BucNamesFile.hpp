#ifndef BUCNAMESFILE_H
#define BUCNAMESFILE_H

#include "InputFile.hpp"

#include<string>
#include <iostream>     // std::cout
#include <fstream>

using namespace std;

/* ==========================================================================

              C L A S S :  BUCNAMES

=============================================================================

*/

class BucNamesFile : public InputFile
{
private:
  std::string compoundname_; //name of the compound under investigation
  int  Nabuc_; // number of atom in the basic unit cell
  int  nas_; // number of atomic species
  std::string  scfile_; // name of the POSCAR file for the supercell
  std::string ucfile_; // name of the POSCAR file for the unitcell
  std::string FCfile_; // name of the FORCE_CONSTANTS file
  std::string msfile_; // name of the material strucuture file
  std::string dummy_; /* dummy variable, read and skip the lines that
                    give you indication about the field to fill */
public:
  BucNamesFile() {;} // constructor
  void read(std::string);
  std::string compound_name() {return compoundname_;} 
  int Na_buc() {return Nabuc_;}
  int n_as() {return nas_;}
  std::string sc_file() {return scfile_;}
  std::string uc_file() {return ucfile_;}
  std::string FC_file() {return FCfile_;}
  std::string ms_file() {return msfile_;}
};


#endif

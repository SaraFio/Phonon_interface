
#include "BucNamesFile.hpp"

/* ________________________________________________________________________

              C L A S S : BUCNAMES : read

   ________________________________________________________________________
*/
void BucNamesFile::read(std::string buc_file) {
    /* open the file with the data
  ================================*/
  std::ifstream is_bucnames(buc_file);
    if (!is_bucnames)
    {
      std::cout << "Unable to open file: " <<buc_file << std::endl;
      exit(1);
    }
    // ----------------------------------------------------------------


    /* read the name of the compound 
    ================================*/
    is_bucnames >> dummy_;
    is_bucnames >> compoundname_;

    // ---------------------------------------------------------------


    /* read the  number of atom in the basic unit cell
    ==================================================*/
    is_bucnames >> dummy_;
    is_bucnames >> Nabuc_;

    // ---------------------------------------------------------------

    /* read the  number of atomic species 
    ======================================*/
    is_bucnames >> dummy_;
    is_bucnames >> nas_;

    // ---------------------------------------------------------------

    /* read the name of the  POSCAR file for the supercell
    ======================================================*/
    is_bucnames >> dummy_;
    is_bucnames>> scfile_;

    // ---------------------------------------------------------------

    /* read the name of the  POSCAR file for the unitcell
    ======================================================*/
     is_bucnames >> dummy_;
    is_bucnames >> ucfile_;

    // ---------------------------------------------------------------

    /* read the name of the FORCE_CONSTANTS file
    ======================================================*/
    is_bucnames >> dummy_;
    is_bucnames >> FCfile_;

    // ---------------------------------------------------------------


    /* read the name of the material structure file
    ======================================================*/
    is_bucnames >> dummy_;
    is_bucnames >> msfile_;

    // ---------------------------------------------------------------

}


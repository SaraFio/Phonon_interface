#ifndef SYMFC_H
#define SYMFC_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

/* ==========================================================================

              FUNCTION :  sym_FC

=============================================================================

this function makes sure that the force constant matrix is invariant under the 
inter-change of the two derivatives.


Input:  FCs_cut: cut Force constant matrix. 
        Na : number of atoms in the supercell

*/


Eigen::MatrixXd sym_FC( Eigen::MatrixXd FCs_cut, int Na);


#endif

#ifndef CUTFC_H
#define CUTFC_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

/* ==========================================================================

              FUNCTION :  cut_FC

=============================================================================

this function cuts the interactions in the force-constant matrix FC after a defined 
cut-off radius rc.


Input:  FCs: Force constant matrix. 
        LMs_ord: layer matrix 
        rc:cut-off radius. Units [Angstrom].

Note that FCs is the sorted FC matrix. (not scaled by the mass)
*/


Eigen::MatrixXd cut_FC(Eigen::MatrixXd LMs_ord, Eigen::MatrixXd FCs, double rc);


#endif

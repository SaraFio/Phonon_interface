/* ==========================================================================

              FUNCTION :  find_vec

=============================================================================
looks if a row defined by vec exists in a matrix and return the index of the row of the matrix
Note that vec is a (1,n)matrix, and matrix is a (m,n) matrix
*/

#ifndef FINDVEC_H
#define FINDVEC_H

#include <fstream>    
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

int find_vec(Eigen::MatrixXd vec, Eigen::MatrixXd matrix, double comp_diff);

#endif 

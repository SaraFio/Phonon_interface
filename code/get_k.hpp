/* ==========================================================================

              FUNCTION :  get_k

============================================================================= 
generates the k-space grid of a rectangular brillouin zone.
 Input:   Nka: number k-points in ka direction (int)
          (a = x,y,z)

 Output: ka: vector containing the Nka ka-points between -pi and pi        
 
nb. actually ka vector is a  (Nka,1) matrix */

#ifndef GETK
#define GETK

#include <Eigen/Dense>
#include<cmath>
#include<math.h>

Eigen::MatrixXd getk(int);

#endif


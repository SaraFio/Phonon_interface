/* ==========================================================================

              FUNCTION :  get_k

============================================================================= 
generates the k-space grid of a rectangular brillouin zone.
 Input:   Nka: number k-points in ka direction (int)
          (a = x,y,z)

 Output: ka: vector containing the Nka ka-points between -pi and pi        
 
nb. actually ka vector is a  (Nka,1) matrix */



/*=====================================================*/

#include "get_k.hpp"

//Here the function is explicited
Eigen::MatrixXd getk(int Nka)
{
const double pi = 3.14159265358979323846; 


     double fac = 2*pi/Nka;
     int colsize = Nka+1; // length of the matrix
     Eigen::MatrixXd ka(colsize,1); // ka matrix definition
     
     if(Nka==0)
       {
         ka(0,0) = 0.0;
       }
     else
       {
         for(int i =0; i<=Nka ; i++)
           {
             // filling the ka matrix
             ka(i,0) = -pi +fac*i; 
           }
       }

  return ka;
}



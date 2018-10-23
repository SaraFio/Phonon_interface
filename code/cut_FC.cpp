
#include "cut_FC.hpp"    

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


Eigen::MatrixXd cut_FC(Eigen::MatrixXd LMs_ord, Eigen::MatrixXd FCs, double rc)
{

// -------------- local variables --------------------

int Na; // number of atom in LMs_ord (supercell)
Na = LMs_ord.rows();

 double x =0;
 double y =0;
 double z =0;
 double r =0;

 // --------------------------------------------------------------

Eigen::MatrixXd FC_cut(3*Na,3*Na); //FC matrix to return after the cutting

FC_cut = Eigen::MatrixXd::Zero(3*Na,3*Na); // inizialization FC_cut

//-----------------------------------------------------------------

for(int IA1=0; IA1 < Na; IA1++)
  {
    for(int IA2=0; IA2 < Na; IA2++)
      {
        x = LMs_ord(IA2,0)-LMs_ord(IA1,0);
        y = LMs_ord(IA2,1)-LMs_ord(IA1,1);
        z = LMs_ord(IA2,2)-LMs_ord(IA1,2);
        
        r = sqrt(x*x +y*y +z*z);
        
        if(r < rc)
          {

            for (int r=0; r<3 ; r++)
              {
                for(int c=0; c<3; c++)
                  {
                    
                    FC_cut(3*IA1+r,3*IA2+c) = FCs(3*IA1+r,3*IA2+c);
                    
                  }
              }
            
          }
      }
 }

 return FC_cut;
}


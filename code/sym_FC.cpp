
#include "sym_FC.hpp"

/* ==========================================================================

              FUNCTION :  sym_FC

=============================================================================

this function makes sure that the force constant matrix is invariant under the 
inter-change of the two derivatives.


Input:  FCs_cut: cut Force constant matrix. 
        Na : number of atoms in the supercell

*/


Eigen::MatrixXd sym_FC( Eigen::MatrixXd FCs_cut, int Na)
{

// -------------- local variables --------------------

  double FC_tmp1, FC_tmp2;

 // --------------------------------------------------------------

Eigen::MatrixXd FCs_cut_sym(3*Na,3*Na); //FCs_cut matrix to return after the sym

FCs_cut_sym = Eigen::MatrixXd::Zero(3*Na,3*Na); // inizialization FC_cut_sym

//-----------------------------------------------------------------

for(int IA1=0; IA1 < Na; IA1++)
  {
    for(int IA2=IA1; IA2 < Na; IA2++)
      {

            for (int r=0; r<3 ; r++)
              {
                for(int c=0; c<3; c++)
                  {
                    
                   FC_tmp1 = FCs_cut(3*IA1+r,3*IA2+c);
                   FC_tmp2 = FCs_cut(3*IA2+c,3*IA1+r);

                   FCs_cut_sym(3*IA1+r,3*IA2+c) = (FC_tmp1+FC_tmp2)/2;
                   FCs_cut_sym(3*IA2+c,3*IA1+r) = (FC_tmp1+FC_tmp2)/2;
                  }
              }
      }
 }

 return FCs_cut_sym;
}


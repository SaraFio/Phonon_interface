
#include"getHtotOMEN.hpp"

/* ==========================================================================

                    FUNCTION  :  getHtotOMEN

=============================================================================

this function generates the total Hamiltonian HtotOMEN (3*Na_uc_new, 3*Na_uc_new)  
at a given k-point out of the tensor Hs_mod 
(NFold_bs_x,3*Na_uc_new,NFold_bs_y*NFold_bs_z*3*Na_uc_new)
 as it is done in OMEN.

input: Hs_mod, NFold_bs, k, Na_uc_new
output: HtotOMEN 
*/

Eigen::MatrixXcd getHtotOMEN(Eigen::Tensor<double, 3> Hs_mod,Eigen::MatrixXd NFold_bs,
                            Eigen::MatrixXd k,int Na_uc_new)
{

Eigen::MatrixXcd HtotOMEN(3*Na_uc_new,3*Na_uc_new);
HtotOMEN.setZero();


//----- local variables -----------

std::complex<double> PhaseH ;
const std::complex<double> Im(0,1);
 int IA1, IA2,ix,iy,iz;

//-------------------------------------------------------------------------

for(int ixt =0; ixt < NFold_bs(0); ixt ++)
  { 
    ix = (ixt+1) -(NFold_bs(0)+1)/2;
     for(int iyt =0; iyt < NFold_bs(1); iyt ++)
       {
         iy = (iyt+1) -(NFold_bs(1)+1)/2;
         for(int izt =0; izt < NFold_bs(2); izt ++)
           {
             iz = (izt+1) -(NFold_bs(2)+1)/2;

             PhaseH = std::exp(Im*k(0,0)*(double)ix)*std::exp(Im*k(0,1)*(double)iy)*std::exp(Im*k(0,2)*(double)iz);
           
             
             // fill HtotOMEN
             for(int IAnuc1=0; IAnuc1 < Na_uc_new; IAnuc1++)
               {
                 for(int IAnuc2=0; IAnuc2 < Na_uc_new; IAnuc2++)
                   {
                     for(int r=0; r < 3; r++)
                       {
                         for(int c=0; c < 3; c++)
                           {
                             HtotOMEN(3*IAnuc1+r, 3*IAnuc2+c)=HtotOMEN(3*IAnuc1+r, 3*IAnuc2+c)+
                               PhaseH*Hs_mod(ixt,3*IAnuc1+r,3*Na_uc_new*(izt+3*iyt) + 3*IAnuc2+c);
                           }
                       }
                   }
               }
           }
       }
  }
 
 
  return HtotOMEN;
}




#include"getHtot.hpp"


/* ==========================================================================

              FUNCTION :  getHtot

=============================================================================

this function generates the  total Hamiltonian Htot out of the force constant matrix 
as it is done in PHONOPY at a given k-point

*/


/*entries:
FCs_m, NFold, kx, ky, kz, LMs_ord, LMs_uc_ord, Lxyz_uc, rvec_short, Nmulti, num_SCELL,comp_diff*/

Eigen::MatrixXcd getHtot(Eigen::MatrixXd FCs_m, Eigen::MatrixXd NFold, double kx, double  ky,double kz,
                          Eigen::MatrixXd LM_coord,Eigen::MatrixXd LM_uc_coord, Eigen::MatrixXd Lxyz_uc,
                          Eigen::MatrixXd rvec_short, Eigen::MatrixXd Nmulti, int num_SCELL,double comp_diff){

  //--- local variables ---------

  int IA1;
  int IA2;
  int L;
  L = LM_coord.cols();
  int Na;
  int Na_uc;
  Na  = LM_coord.rows();
  Na_uc = LM_uc_coord.rows();
  int row_ind;
  Eigen::MatrixXd xyz1(1,L);
  Eigen::MatrixXd xyz2(1,L);
  std::complex<double> phaseH ;
  const std::complex<double> Im(0,1);
  // ----------------------------

  // -- matrix to return ------

    Eigen::MatrixXcd Htot(3*Na_uc,3*Na_uc);
    Htot.setZero();
  
  //-------------------------

  // ---check if LMs_ord and LMs_uc_ord have the same number ---
  int L1;
  L1=  LM_uc_coord.cols();
  if (L1!=L) 
    {
      std::cout << "\n ERROR: LMs_ord and LMs_uc_ord have != numb of columns" << std::endl;
      exit(1);
    }
  //------------------------------------



  for(int IAuc1 = 0; IAuc1 < Na_uc; IAuc1 ++)
    {
      
      xyz1 = LM_uc_coord.row(IAuc1);

      /* find the position of the IAuc1 atom in LMs_ord 
       you can read this function in "find_vec.h"*/
      IA1 = find_vec(xyz1,LM_coord,comp_diff);


      for(int IAuc2=0; IAuc2 < Na_uc; IAuc2 ++)
        {
           for(int ix = 0; ix < NFold(0); ix++)
             {
               for(int iy = 0; iy < NFold(1); iy++)
                 {
                   for(int iz =0; iz < NFold(2); iz++)
                     {                      
                       
                       xyz2(0,0) = LM_uc_coord(IAuc2,0) + ix*Lxyz_uc(0);
                       xyz2(0,1) = LM_uc_coord(IAuc2,1) + iy*Lxyz_uc(1);
                       xyz2(0,2) = LM_uc_coord(IAuc2,2) + iz*Lxyz_uc(2);
                       
                       IA2 = find_vec(xyz2,LM_coord, comp_diff);
                     
                      
                       
                       for(int j =0; j < num_SCELL; j++)
                         {
                           row_ind = IAuc1*Na*num_SCELL + IA2*num_SCELL + (j);
                           
                           if(rvec_short(row_ind,10)!=0)
                             {

                               phaseH= (std::exp(Im*kx/Lxyz_uc(0)*rvec_short(row_ind,6))*std::exp(Im*ky/Lxyz_uc(1)*rvec_short(row_ind,7))
                                        *std::exp(Im*kz/Lxyz_uc(2)*rvec_short(row_ind,8)))/Nmulti(IAuc1,IA2);
                               // std::cout<< phaseH << std::endl;
                               
                               for (int r=0; r<3 ; r++)
                                 {
                                   for(int c=0; c<3; c++)
                                     {
  
                                        Htot(3*IAuc1+r,3*IAuc2+c) =Htot(3*IAuc1+r,3*IAuc2+c) +phaseH*FCs_m(3*IA1+r,3*IA2+c);

                                     }
                                 }
                            }
                         }     
                     }
                 }      
             }
        } 
    }



  return Htot;

}


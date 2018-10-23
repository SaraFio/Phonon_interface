#include "FindNN.hpp"

#include <fstream>    
#include <iomanip>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include<string>
#include<cmath>


/* ________________________________________________________________________

              C L A S S : FindNN

   ________________________________________________________________________
*/
// Constructor
FindNN:: FindNN(int Na_uc, int Na,Eigen::MatrixXd LMs_ord,Eigen::MatrixXd LMs_uc_ord,Eigen::MatrixXd Lxyz,
Eigen::MatrixXd Lxyz_uc,Eigen::MatrixXd ixyz1_ref,Eigen::MatrixXd NSCELL, double comp_diff):
  Nauc_(Na_uc),  Na_(Na), LMs_ord_(LMs_ord), LMs_uc_ord_(LMs_uc_ord),Lxyz_(Lxyz),
  Lxyzuc_(Lxyz_uc), ixyz1_ref_(ixyz1_ref),NSCELL_(NSCELL), compdiff_(comp_diff) {}


//##########################################################################

void FindNN::find_nn()
{

  /* to find the nearest neighbors, We need to evaluate the distance between atoms: r12
    the 1st atom is in the unit cell
    the 2nd atom is in one of the supercell
    
    The information will be stored in the matrix  rvec_short_
    dim rvec_short(Na_uc*Na*(NSCELL_x*2+1)*(NSCELL_y*2+1)*(NSCELL_z*2+1),11)
    Columns of rvec_short =
    (atom in uc,  atom in sc, NSCELL_x , NSCELL_y , NSCELL_z,n_SCELL ,x12, y12, z12, r12, nn_logical)
    where NSCELL_i = coordinates of the supercell  (ex. 1,-1,0)
          n_SCELL = supercell index (ex 23)
          x12 = x2 - x1, y12 = y2 - y1, z12 = z2 - z1
          r12 = sqrt(x12*x12 + y12*y12 + z12*z12)
          nn_logical is 1 if ( r12-(min dist 1-2) ) < comp_diff, otherwise 0

   The filling of rvec_short  will procede in 2 steps.
  
  I  ) fill cols (1,..,10) =  (atom in uc,  atom in sc, NSCELL_x , NSCELL_y , NSCELL_z, n_SCELL ,x12, y12, z12, r12)
  
  II ) fill col(11) = (nn_logical)

  Then the matrix Nmulti (dim (Na_uc , Na)) is filled

  III)  for each fixed couple atom1 -atom2 the sum of the nn_logical is performed
*/

 

  /*=======================================================
    part I : fill rvec_short cols(1,..,9)
  */

  //-----  local variables --------------------------------

  Eigen::MatrixXd xyz1(1,3); // coordinates of the 1st atom 
  Eigen::MatrixXd xyz2(1,3); // coordinates of the 2nd atom 

  int row_size; /*is used to set the row-dim of rvec_short_ 
                remember 
                dim rvec_short(Na_uc*Na*(NSCELL_x*2+1)*(NSCELL_y*2+1)*(NSCELL_z*2+1),11)*/

  row_size = Nauc_*Na_*(NSCELL_(0,0)*2+1)*(NSCELL_(0,1)*2+1)*(NSCELL_(0,2)*2+1);
  
  int num_SCELL; //number of supercells: (NSCELL_x*2+1)*(NSCELL_y*2+1)*(NSCELL_z*2+1) 

  num_SCELL = (NSCELL_(0,0)*2+1)*(NSCELL_(0,1)*2+1)*(NSCELL_(0,2)*2+1);

  //----- set dim of matrix (private attribute) ----------  

  rvec_short_.resize(row_size,11);
 

  //------------------------------------------------------
  

 int ia = 0; 
    
 int I=0 ;

      // --- 1st atom : unit cell ------------
      for(int ia1 = 0; ia1 < Nauc_; ia1++)
        {
          for(int j = 0; j < 3; j++)
            {
              // coordinates of the 1st atom 
              xyz1(0,j) = LMs_uc_ord_(ia1,j)+ixyz1_ref_(0,j)*Lxyzuc_(0,j);
            }
          // --- 2nd atom : super cell ------------          
          for(int ia2 = 0; ia2 < Na_; ia2++) 
            {
             //----in which supercell the 2nd atom is --------
              for(int ix = -NSCELL_(0,0); ix < NSCELL_(0,0)+1; ix ++)
                {
  
                   for(int iy = -NSCELL_(0,1); iy < NSCELL_(0,1)+1; iy ++)
                    {
                      for(int iz = -NSCELL_(0,2); iz < NSCELL_(0,2)+1; iz ++)
                        {

                          //--- supercell index I -----
                          if(I==num_SCELL)
                            {
                              I=1;
                            }
                          else
                            {
                              I++; 
                            }
                          //----------------------------
                          rvec_short_(ia,0) = ia1; // col 1: atom in uc 
                          rvec_short_(ia,1) = ia2; // col 2: atom in sc 
                          rvec_short_(ia,2) = ix; // col 3: which sc (x-coord)
                          rvec_short_(ia,3) = iy; // col 4: which sc (y-coord)
                          rvec_short_(ia,4) = iz; // col 5: which sc (z-coord)
                          rvec_short_(ia,5) = I; // col 6: n_SCELL
                         
                          

                          for(int j=0; j < 3; j++)
                            {  
                              /* coordinates of the 2nd atom 
                               note that is give by the coordinates of the atom in the supercell
                               + the coordinates of the supercell*/
                          xyz2(0,j) = LMs_ord_(ia2,j)+rvec_short_(ia,2+j)*Lxyz_(0,j);
                              
                              // 1st - 2nd atom coordinates difference : col(7,8,9)=(x12,y12,z12)
                              rvec_short_(ia, 6+j)=xyz2(0,j)-xyz1(0,j);
                            }

                          
                          // 1st - 2nd atom distance:  col 10= r12
                          double tmp = 0.0;

                          tmp = rvec_short_(ia,6)*rvec_short_(ia,6) + rvec_short_(ia,7)*rvec_short_(ia,7) + rvec_short_(ia,8)*rvec_short_(ia,8);
                          
                          rvec_short_(ia,9)=sqrt(tmp);
                          
                          ia ++;
                        }
                    }
                  
                }
              
            }
        }
  /*=======================================================
    part II : fill rvec_short col(11)

    Here we find, for each couple atom 1-2,  which is the minimum distance, along all the cells.
    for fixed atom 1 and fixed atom 2, r12 is studied  along all the cells. 
    At the end we'll have a matrix(Natoms,1) where all the minimum r12 are stored*/

  //-----  local variables --------------------------------

      int Natoms = Nauc_*Na_;
      int num_NSCELL = (NSCELL_(0,0)*2+1)*(NSCELL_(0,1)*2+1)*(NSCELL_(0,2)*2+1);
      int start = 0;
      int finish = num_NSCELL;

      int min_ind;
      double min;
      Eigen::MatrixXd min_dist(Natoms,1); // min dist atom 1-2
      
      //--- find the min distances ---------------------
      for(int j = 0; j < Natoms; j++ )
        {
          min_ind = start;
          min = rvec_short_(start,9);
         
          for(int i = start; i < finish; i++)
            {
              if(rvec_short_(i,9)<min) 
                {
                  min=rvec_short_(i,9);
                }
            }
          min_dist(j,0)=min;
          start = finish;
          finish += num_NSCELL;
        }

      //---- fill col 11 (true =1, false =0)-----

     int  iat =0;
     int im =0;
      //  atom 1
      for(int ia1 = 0; ia1 < Nauc_; ia1++)
        {
          // atom 2
          for(int ia2 = 0; ia2 < Na_; ia2++) 
            {
              //----in which supercell the 2nd atom is --------
              for(int ix = -NSCELL_(0,0); ix < NSCELL_(0,0)+1; ix ++)
                {
                  for(int iy = -NSCELL_(0,1); iy < NSCELL_(0,1)+1; iy ++)
                    {
                      for(int iz = -NSCELL_(0,2); iz < NSCELL_(0,2)+1; iz ++)
                        {
                          
                          if(abs(rvec_short_(iat,9)-min_dist(im,0))< compdiff_)
                            {
                              rvec_short_(iat,10) = 1;
                            }
                          else 
                            {
                              rvec_short_(iat,10) = 0; 
                            }
                          iat ++;
                        }
                    }
                }
              im++;
            }

        }



  /*=======================================================
    part III : fill Nmulti_
  */

  
  //----- set dim of matrix (private attribute) ----------  

  Nmulti_.resize(Nauc_,Na_);
 
  //------------------------------------------------------

  start = 0;
  finish=num_NSCELL;

  int tmp=0;


  for(int i = 0; i < Nauc_; i++)
    {

      for(int j = 0; j<Na_; j++)
        {

          for(int k = start; k < finish; k++)
            {
              tmp += rvec_short_(k,10);
            }
          
          Nmulti_(i,j) =tmp;
          tmp = 0;         
          start = finish;
          finish += num_NSCELL;
        }
    }

}

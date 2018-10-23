
#include "asr_sym_scale.hpp"



/* ==========================================================================

                    FUNCTION  :  asr  (Acoustic Sum Rule) 
                                 symmetrizing
                                 scale

=============================================================================
in this function, we
1)impose the acoustic sum rule 
and 
2)symmetrizing the different parts of the force-constant matrix 
3) Scale the elements of force constant matrices with the corresponding masses 
    --> dynamical matrix.

==============================================================================
1) ASR

1a) getHtotOMEN  at k = (0.0,0.0,0.0)
    this function generates the total Hamiltonian HtotOMEN (3*Na_uc_new, 3*Na_uc_new)  
    at a given k-point out of the tensor Hs_mod 
    (NFold_bs_x,3*Na_uc_new,NFold_bs_y*NFold_bs_z*3*Na_uc_new)
    as it is done in OMEN.

    Note that HtotOMEN will be used to calculate the new phonon dispersion

    input: Hs_mod, NFold_bs, k, Na_uc_new
    output: HtotOMEN 

1b) the Acoustic Sum Rule is imposed

    C(I,J,a,b) = IFC of (atom I direction a) & (atom J direction b)
    I,J = .., Na_uc_new
    a,b = x,y,z
    
    ASR : C(I,I,a,b) = - SUM_{J != I} C(I,J,a,b)
    
*/

Eigen::Tensor<double, 3> asr_sym_scale(Eigen::Tensor<double, 3> Hs_mod_input,
                                       Eigen::MatrixXd NFold_bs,int Na_uc_new, 
                                       Eigen::MatrixXd LMs_uc_new, Eigen::MatrixXd mr)
{

  Eigen::Tensor<double, 3> Hs_mod;
  Hs_mod = Hs_mod_input;



//1a)--------------------------------------------------------------------------------------------
    Eigen::MatrixXcd HtotOMEN_zero(3*Na_uc_new,3*Na_uc_new); // matrix (3*Na_uc_new, 3*Na_uc_new) obtained from Hs_mod
    HtotOMEN_zero.setZero();

    Eigen::MatrixXd k(1,3); // k point (kx,ky,kz)
    k << 0.0,0.0,0.0;
    
    HtotOMEN_zero= getHtotOMEN( Hs_mod, NFold_bs, k, Na_uc_new);

   

//1b)------   A S R  ---------------------------------------------------------------


    Eigen::MatrixXd temp_Hs(3*Na_uc_new*3,1); 
    temp_Hs.setZero(3*Na_uc_new*3,1);

    int y,z;

   

    if(NFold_bs(1)==1) y = 0;
    else y = 1;
    
    if(NFold_bs(2)==1) z = 0;
    else z = 1;

    double sum;
    sum = 0.0;
int N,ix,iy,iz;
    
   
    
    for (int IA1 = 0; IA1 < Na_uc_new; IA1 ++)
      {
        for(int r= 0; r < 3; r++)
          {
            for(int c= 0; c < 3; c++)
              {
                if(IA1 == Na_uc_new -1) N = Na_uc_new-1;
                else  N = Na_uc_new;
                
                for(int IA2 = 0; IA2 < N; IA2++)
                  {
                
                    if(IA2 == IA1) IA2 ++; // exclude it self
                 
                    // compute the sum off-site (correction)
                    sum = sum + HtotOMEN_zero(3*IA1+r,3*IA2+c).real();
                    
                  }

             
                Hs_mod(1,3*IA1+r + 3*Na_uc_new*y, 3*Na_uc_new*z+3*IA1+c) = -sum;

   
                
                sum = 0.0;


                /* note that exluding it self we are missing the 
                   contribution of the "same" atom in the repeated
                   cell (due to the folding in HtotOMEN_zero)
                   so, here we add this contribution*/
                
                for(int ixt =0; ixt < NFold_bs(0); ixt ++)
                  { 
                    ix = (ixt+1) -(NFold_bs(0)+1)/2;
                    for(int iyt =0; iyt < NFold_bs(1); iyt ++)
                      { 
                       
                        iy = (iyt+1) -(NFold_bs(1)+1)/2;
                        
                        for(int izt =0; izt < NFold_bs(2); izt ++)
                          {
                            iz = (izt+1) -(NFold_bs(2)+1)/2;
                            

                            // excluding the central matrix (0,0,0)
                            if (((ix==0)&&(iy==0)&&(iz==0))==0)
                              {
                                // std::cout << ix << "\t" << iy << "\t" << iz << std::endl;
                              Hs_mod(1,3*IA1+r + 3*Na_uc_new*y, 3*Na_uc_new*z+3*IA1+c) = 
                                Hs_mod(1,3*IA1+r + 3*Na_uc_new*y, 3*Na_uc_new*z+3*IA1+c) 
                              - Hs_mod(ixt,3*IA1+r +3*Na_uc_new*iyt,3*Na_uc_new*izt + 3*IA1+c);
                              
                              }
                            
                          }
                      }
                  }

                
                
              }
          }      
      }



    //======================================================================================

   // 2) ------------- symmetrizing ---------------------

    Eigen::MatrixXd temp_1(3,3);
    temp_1.setZero(3,3);

    int Im1 = 0; /* this index reads in LMs_uc_new matrix the atomic 
                    species index, such as 1, 2, .. , n_as
                    it is used to assign the value of the mass
                    Note. the  atomic species indeces starts from 1,
                    while indeces in matrices from zero
                    */
    int Im2 = 0; /* this index reads in LMs_uc_new matrix the atomic 
                    species index
                    it is used to assign the value of the mass.*/
    
    double m1=0.0; /* mass of the atomic species i*/
    double m2=0.0; /* mass of the atomic species j*/

    double denom=0.0 ; /*scale factor sqrt(m1m2)*/
    
    int ind1,ind2;

    for (int IA1 = 0; IA1 < Na_uc_new; IA1 ++)
      {

        Im1 = LMs_uc_new(IA1,3); // index of atom1
        m1 = mr(Im1-1); // mass of atom1

        for(int r= 0; r < 3; r++)
          {
            for(int c= 0; c < 3; c++)
              {

                temp_1(r,c) = Hs_mod(1,3*IA1+r + 3*Na_uc_new*y, 3*Na_uc_new*z+3*IA1+c);

                Hs_mod(1,3*IA1+r + 3*Na_uc_new*y, 3*Na_uc_new*z+3*IA1+c) = (temp_1(r,c)+temp_1(c,r))/2.0;
                
                for (int IA2 = 0; IA2 < Na_uc_new; IA2 ++)
                  {
                    
                    Im2 = LMs_uc_new(IA2,3); // index of atom2
                    m2 = mr(Im2-1); // mass of atom2

                    denom =1./ sqrt(m1*m2);  

                    // --- mass scaling -----
                    
                    for(int ixt =0; ixt < NFold_bs(0); ixt ++)
                      { 
                        for(int iyt =0; iyt < NFold_bs(1); iyt ++)
                          { 
                            for(int izt =0; izt < NFold_bs(2); izt ++)
                              {
                                
                                ind1 =3*IA1+r + 3*Na_uc_new*iyt;
                                ind2 =3*IA2+c + 3*Na_uc_new*izt;
                                
                                
                                Hs_mod(ixt,ind1,ind2) =  Hs_mod(ixt,ind1,ind2)*denom;
                                
                              }
                          }
                      }
                                
                    
                  }
                
                
              }
          }      
      }
    
 




    return Hs_mod;

}


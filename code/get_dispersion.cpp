

#include"get_dispersion.hpp"


using Eigen::MatrixXd;
/* 
 #############################################################################

                      GET PHONON BAND STRUCTURE

 #############################################################################

The phonon band structure is calculated using the Hs_mod tensor.



-    for each k-point ka(i) a matrix Hamiltonian is created H(ka(i))
     by the function getHtotOMEN .

-    the corresponding eigenvalues E(ka(i)) are calculated using the 
     Eigen attribute .eigenvalues()
     all the  E(ka(i)) are stored in matrix E -> dim(E)=(3*Na_uc_new,ntotk)
     where ntotk = (Nkx+1)+(Nky+1)+(Nkz+1)
     Note that E(ka(i)) is a vector of length 3*Na_uc_new. 
     The first Nk(0)+1 columns of E contain E(kx), 
     the second Nk(1)+1 columns contain E(ky),
     the last Nk(2)+1 columns contain E(kz)
     
-    in order to obtain the phonon energies 
     E(i,j)=sign*hq*sqrt(sign*E(i,j));
     where hq = hbar/charge;
           sign = -1 if E(i,j)<0
                   1 otherwise

     The E matrix columns are returned in ascending order

*/ 


Eigen::MatrixXd get_dispersion(Eigen::Tensor<double, 3> Hs_mod, Eigen::MatrixXd NFold_bs,
                               Eigen::MatrixXd kx, Eigen::MatrixXd ky, Eigen::MatrixXd kz,
                               int indE, int Na_uc, double hbar, double charge)

{


       
// =============  get Htot(ka(i)) =====================================



Eigen::MatrixXcd Htot(3*Na_uc,3*Na_uc); /* here the Hamiltonian created by getHtot
                                        will be stored*/

 //inizialize the Htot: zeros matrix
 Htot.setZero();

Eigen::MatrixXd E(3*Na_uc, indE); /* Eigenvlues matrix 
                                     = phonon energies for each ka(i)*/
 E.setZero();

 Eigen::MatrixXd k(1,3); // k-point vector in which calculate the Htot




 //--------------------
          k(0,0)=kx(0);
          k(0,1)=ky(0);
          k(0,2)=kz(0);
          
 
 Htot =  getHtotOMEN( Hs_mod, NFold_bs, k,Na_uc);    
 // std::cout << "\n" << std::setprecision(4)<<Htot.real()/1e27 << std::endl;

/* For each k point the Htot is calculated and the corresponding eigenvalues E(ka(i))
are stored. Note that E(ka(i)) is a vector of length 3*Na_uc. 
The first Nk(0)+1 columns contain E(kx), the second Nk(1)+1 columns contain E(ky),
the last Nk(2)+1 columns contain E(kz)*/

for(int ikx = 0; ikx < kx.rows(); ikx++)
  {
    for(int iky = 0; iky < ky.rows(); iky++)
      {
        for(int ikz = 0; ikz < kz.rows(); ikz++)
        {
          // -- fill the k-point vector 
          k(0,0)=kx(ikx);
          k(0,1)=ky(iky);
          k(0,2)=kz(ikz);
          
          //---------- Htot(ka(i)) is created--------------------------------
            Htot =  getHtotOMEN( Hs_mod, NFold_bs, k,Na_uc);           



//3rd)======================== E(ka(i)) is computed -> E   ==================================
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
            ces.compute(Htot);

            E.col(indE-1) = ces.eigenvalues().real();
            // note that E is filled from the last column to the first           

            indE--; 

            //----------re-set to zeros Htot for the next ka(i)------------ 

            Htot.setZero();
          }
      }
  }



//4th) ========= scale and sort the eigenvalues = phonon energies  ==================
            
 double hq = hbar/charge;
 int sign =1;


 for(int j=0; j<E.cols();j++)
   {
     for(int i =0; i < E.rows(); i++)
       {
         if(E(i,j)<0) sign =-1;
         else sign =1;
         E(i,j)=sign*hq*sqrt(sign*E(i,j));
       }     
     // --- sort in ascending order -------
     std::sort(E.col(j).data(), E.col(j).data() + E.rows());

   }
return E;
}

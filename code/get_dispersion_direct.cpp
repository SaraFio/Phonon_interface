
#include"get_dispersion_direct.hpp"




using Eigen::MatrixXd;
/* 
 #############################################################################

                      GET PHONON BAND STRUCTURE

 #############################################################################

The phonon band structure is calculated using the Force Constant matrix.

1st) the k-space grid of a rectangular brillouin zone is generated
     by the function get_k(Na). 
     Input:   Nka: number k-points in ka direction (int)  (a = x,y,z)
     Output: ka: vector containing the Nka ka-points between -pi and pi     
     nb. actually ka vector is a  (Nka,1) matrix

2nd) for each k-point ka(i) a matrix Hamiltonian is created H(ka(i))
     by the function getHtot.

3rd) the corresponding eigenvalues E(ka(i)) are calculated using the 
     Eigen attribute .eigenvalues()
     all the  E(ka(i)) are stored in matrix E -> dim(E)=(3*Na_uc,ntotk)
     where ntotk = (Nkx+1)+(Nky+1)+(Nkz+1)
     Note that E(ka(i)) is a vector of length 3*Na_uc. 
     The first Nk(0)+1 columns of E contain E(kx), 
     the second Nk(1)+1 columns contain E(ky),
     the last Nk(2)+1 columns contain E(kz)
     
4th) in order to obtain the phonon energies 
     E(i,j)=sign*hq*sqrt(sign*E(i,j));
     where hq = hbar/charge;
           sign = -1 if E(i,j)<0
                   1 otherwise

     The E matrix columns are returned in ascending order

*/ 


Eigen::MatrixXd get_dispersion_direct(Eigen::MatrixXd FCs_m, Eigen::MatrixXd NFold,
                                      Eigen::MatrixXd  kx, Eigen::MatrixXd  ky, Eigen::MatrixXd kz,
                                      int indE, int Na, int Na_uc,
                                      Eigen::MatrixXd LMs_ord,Eigen::MatrixXd LMs_uc_ord, Eigen::MatrixXd Lxyz_uc,
                                      Eigen::MatrixXd rvec_short, Eigen::MatrixXd Nmulti,double comp_diff, 
                                      double hbar, double charge)

{

       
// =============  get Htot(ka(i)) =====================================


//--- temporary LM matrixes are created --------------------------------------- 

/* LM_coord is  LMs_ord matrix, reduced to only coordinates 
   (rid off col 4 = atomic species, col 5 = original order) */    
 Eigen::MatrixXd LM_coord(Na,3); 

for(int i = 0; i < 3 ; i++)
  {
LM_coord.col(i) = LMs_ord.col(i);
  } 
            
/* LM_uc_coord is  LMs_uc_ord matrix, reduced to only coordinates 
   (rid off col 4 = atomic species, col 5 = original order)*/
Eigen::MatrixXd LM_uc_coord(Na_uc,3); 

for(int i = 0; i < 3 ; i++)
  {
LM_uc_coord.col(i) = LMs_uc_ord.col(i);
  }

//-------------------------------------------------------------------------



 int num_SCELL = rvec_short.rows()/(Na_uc*Na); //number of artificially repeated supercell


Eigen::MatrixXcd Htot(3*Na_uc,3*Na_uc); /* here the Hamiltonian created by getHtot
                                        will be stored*/

Htot = Eigen::MatrixXd::Zero(3*Na_uc,3*Na_uc); //inizialize the Htot: zeros matrix
 Htot.setZero();

Eigen::MatrixXd E(3*Na_uc, indE); /* Eigenvlues matrix 
                                     = phonon energies for each ka(i)*/



/* For each k point the Htot is calculated and the corresponding eigenvalues E(ka(i))
are stored. Note that E(ka(i)) is a vector of length 3*Na_uc. 
The first Nk(0)+1 columns contain E(kx), the second Nk(1)+1 columns contain E(ky),
the last Nk(2)+1 columns contain E(kz)*/

 Htot = getHtot( FCs_m,NFold, -3.1416 ,0.0 ,0.0, 
                            LM_coord, LM_uc_coord, Lxyz_uc, rvec_short, 
                            Nmulti, num_SCELL ,comp_diff);


for(int ikx = 0; ikx < kx.rows(); ikx++)
  {
    for(int iky = 0; iky < ky.rows(); iky++)
      {
        for(int ikz = 0; ikz < kz.rows(); ikz++)
        {

          //---------- Htot(ka(i)) is created--------------------------------
            Htot = getHtot( FCs_m,NFold, kx(ikx), ky(iky), kz(ikz), 
                            LM_coord, LM_uc_coord, Lxyz_uc, rvec_short, 
                            Nmulti, num_SCELL ,comp_diff);  



//)======================== E(ka(i)) is computed -> E   ==================================
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
            ces.compute(Htot);

            E.col(indE-1) = ces.eigenvalues().real();
            // note that E is filled from the last column to the first           

            indE--; 

            //----------re-set to zeros Htot for the next ka(i)------------ 
            //Htot = Eigen::MatrixXd::Zero(3*Na_uc,3*Na_uc);
            Htot.setZero();
          }
      }
  }

//) ========= scale and sort the eigenvalues = phonon energies  ==================
            
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

//--- delete LM_coord and LM_uc_coord----------------------------------

LM_coord.resize(0,0);
LM_uc_coord.resize(0,0);

//----------------------------------------------------------------------

return E;
}

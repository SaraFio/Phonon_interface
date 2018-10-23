#include <fstream>    
#include <iomanip>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <unsupported/Eigen/CXX11/Tensor>


#include "InputFile.hpp"
#include "BucNamesFile.hpp"
#include "MaterialStructureFile.hpp"
#include "InputFile_tosort.hpp"
#include "PoscarFile.hpp"
#include "ForceConstFile.hpp"
#include "FindNN.hpp"
#include "get_k.hpp"
#include "getHtot.hpp"
#include "get_dispersion_direct.hpp"
#include "cut_FC.hpp"
#include "sym_FC.hpp"
#include "check_slabs.hpp"
#include "extract_submatrices.hpp"
#include "getHtotOMEN.hpp"
#include "asr_sym_scale.hpp"
#include"get_dispersion.hpp"


using Eigen::MatrixXd;
using namespace std;






int main(){

 

/* 
 #############################################################################

                    READ INPUT FILE + SORT ACCORDING TO OMEN

 #############################################################################

classes used: BucFileNames
              MaterialStructureFile
              PoscarFile
              ForceConstFile

Here all the input files are read and their quantities stored. 
Input files: buc_file, sc_file, uc_file, FC_file, ms_file
Note that buc_file and ms_file are just read,
while some operation are performed on 
=sc_file, uc_file : 
  + if sc_file and uc_file are written in Direct coord -> converted to cartesian
  + the LM and LM_uc matrixes (matrix of coordiantes) are sorted according to OMEN
= FC_file:
  + is converted from eV/A^2 -> J/m^2
  + sorted according to OMEN
  + scaled by mass

*/ 

  /*=======================================================================*/

  // buc_file  BASIC UNIT CELL + NAMES

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Reading BUC NAMES file " << endl;
  cout << "-----------------------------------------------" << endl;


  BucNamesFile bucnames; 

  bucnames.read("Input_files/buc_names.txt");

  string compound_name; //name of the compound under investigation
  int Na_buc; // number of atom in the basic unit cell
  int n_as; // number of atomic species
  string sc_file; // name of the POSCAR file for the supercell
  string uc_file; // name of the POSCAR file for the unitcell
  string FC_file; // name of the FORCE_CONSTANTS file
  string ms_file; // name of the material strucuture file

  compound_name = bucnames.compound_name();
  Na_buc = bucnames.Na_buc();
  n_as = bucnames.n_as();
  sc_file = bucnames.sc_file();
  uc_file = bucnames.uc_file();
  FC_file = bucnames.FC_file();
  ms_file = bucnames.ms_file();

  cout << "\n basic-unit cell data & file names read \n" << endl; 


  /*=======================================================================*/

  // ms_file: MATERIAL STRUCTURE INFOS
 
  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Reading MATERIAL STRUCTURE file " << endl;
  cout << "-----------------------------------------------" << endl;

 
  MaterialStructureFile mat_str(n_as);//constructor

  double hbar; // h_bar in Joule second
  double charge;//charge in Coulomb
  double m_scalefac;// mass scale factor
  double tol_lim; //numerical criteria  for tolerance
  Eigen::MatrixXd mr; //reduced masses of atomic species matrix(1,n_as)
  Eigen::MatrixXd mt; /*total masses of atomic species matrix(1,n_as) 
                        mt  = mr*mscalefac      */
  Eigen::MatrixXd NFold; /*Defines the supercell size based on the  unit-cell.
                           [x-dir,y-dir,z-dir]. (Rectangular unit-cell!) */
  Eigen::MatrixXd NFold_bs; /* Periodically extended directions get a 3, confined
                               directions get a 1. [x-dir,y-dir,z-dir],eg. wire: [3 1 1]*/


  Eigen::MatrixXd ixyz1_ref; /*Lower left corner of the reference unitcell
                              dim(ixyz1_ref) = (1,3)*/ 
  Eigen::MatrixXd NSCELL; /* number of artificially genereated repeated supercells
                                along the x,y,z direction
                              dim(NSCELL) = (1,3)*/ 
  Eigen::MatrixXd Nk; /* Number k-points in ka direction (a=x,y,z)
                          Nk = [Nkx, Nky, Nkz]  dim(Nk) = (1,3)*/ 

  double rc;  /*radius cutoff (used to cut the FC) Angstrom*/


  mat_str.read(ms_file); // read the material structure file

  hbar= mat_str.hbar();
  charge=mat_str.charge();
  m_scalefac= mat_str.m_scalefac();
  tol_lim=mat_str.tol_lim();
  mr=mat_str.mr();
  mt=mat_str.mt();
  NFold=mat_str.NFold();
  NFold_bs=mat_str.NFold_bs();
  ixyz1_ref=mat_str.ixyz1_ref();
  NSCELL = mat_str.NSCELL();
  Nk = mat_str.Nk();
  rc = mat_str.rc();

  cout << "\n material strucure data read \n" << endl;

  /*=======================================================================*/

  // sc_file SUPERCELL DATA

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Reading SUPERCELL file " << endl;
  cout << "-----------------------------------------------" << endl;


  PoscarFile pos(n_as);//constructor


    int Na; // number of atom in the supercell
  
    double scale_fac; // scale factor supercell 

    Eigen::MatrixXd latt_vec; // matrix(Na_buc,3) of the components of lattice vectors (supercell)
    
    Eigen::MatrixXd Naps; // Matrix Naps(1,n_as) : Number of Atoms Per atomic Species (supercell)
    
    std::string coordunit; // Direct or Cartesian coordinates (supercell)
    
    Eigen::MatrixXd LM_onlycoord; /*matrix(Na,3) : coordinates of the atoms (supercell)
                                    col of LM = (x-coord,y-coord,z-coord)*/
    
    Eigen::MatrixXd LM; /*matrix(Na,4) :(supercell) 
                          coordinates of the atoms + atomic species
                          col of LM = (x-coord,y-coord,z-coord, at. spec.)*/
    
    Eigen::MatrixXd LMs; /*matrix(Na,4) :(supercell)
                         sorted according to OMEN convention
                         coordinates of the atoms + atomic species
                         col of LMs = (x-coord,y-coord,z-coord, at. spec.)*/
    
    Eigen::MatrixXd LMs_ord; /*matrix(Na,5) :(supercell)
                               sorted according to OMEN convention
                               coordinates of the atoms + atomic species + 
                               correspondance original LM and LMs
                               col of LMs = (x-coord,y-coord,z-coord, at. spec., correspondance)
                               ex. if a coordinate in the original matrix LM is in position 3 and 
                               in the sorted matrix LMs is in position 9, we have LMs_ord(3,5) = 9*/
    Eigen::MatrixXd Lxyz(1,3);/* Extension of the supercell in x, y, and z directions. 
                             (Rectangular unit-cell) */


    pos.read(sc_file); // read the supercell file

    scale_fac = pos.scale_fac();
    latt_vec = pos.latt_vec();
    Na = pos.Na();
    Naps = pos.Naps();
    LM_onlycoord = pos.LMcoord();
    LM = pos.LM();

    cout << "\n supercell poscar file read \n" << endl;

    pos.sort(); /* LM matrix is sorted according to OMEN conventions 
                   Note that if the coordinates are written in Direct
                  coordinates, they are converted in Cartesian*/
    
    LMs = pos.LMs();
    LMs_ord = pos.LMs_ord();

    Lxyz(0,0) = latt_vec(0,0)*scale_fac;
    Lxyz(0,1) = latt_vec(1,1)*scale_fac;
    Lxyz(0,2) = latt_vec(2,2)*scale_fac;

    //The LMs_ord  matrix  is printed out in a file
     ofstream lms("Output_data/LMs_ord.dat");

    if (lms.is_open())
      lms<< LMs_ord << endl;

    cout << "\n LM sorted and written in LMs_ord.dat \n" << endl;

  /*=======================================================================*/

  // uc_file UNIT CELL DATA

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Reading UNIT CELL file " << endl;
  cout << "-----------------------------------------------" << endl;


    PoscarFile posuc(n_as);//constructor


    int Na_uc; // number of atom in the unitcell
  
    double scale_fac_uc; // scale factor unitcell 

    Eigen::MatrixXd latt_vec_uc; // matrix(Na_buc,3) of the coomponents of lattice vectors (unitcell)
    
    Eigen::MatrixXd Naps_uc; // Matrix Naps(1,n_as) : Number of Atoms Per atomic Species (unitcell)
    
    std::string coordunit_uc; // Direct or Cartesian coordinates (unitcell)
    
    Eigen::MatrixXd LM_uc_onlycoord; /*matrix(Na_uc,3) : coordinates of the atoms (unitcell)
                                    col = (x-coord,y-coord,z-coord)*/
    
    Eigen::MatrixXd LM_uc; /*matrix(Na_uc,4) :(unitcell) 
                          coordinates of the atoms + atomic species
                          col of LM_uc = (x-coord,y-coord,z-coord, at. spec.)*/
    
    Eigen::MatrixXd LMs_uc; /*matrix(Na_uc,4) :(unitcell)
                         sorted according to OMEN convention
                         coordinates of the atoms + atomic species
                         col of LMs_uc = (x-coord,y-coord,z-coord, at. spec.)*/
    
    Eigen::MatrixXd LMs_uc_ord; /*matrix(Na_uc,5) :(unitcell)
                               sorted according to OMEN convention
                               coordinates of the atoms + atomic species + 
                               correspondance original LM_uc and LMs_uc
                               col of LMs_uc = (x-coord,y-coord,z-coord, at. spec., correspondance)*/

    Eigen::MatrixXd Lxyz_uc(1,3);/* Extension of the unit cell in x, y, and z directions. 
                              (Rectangular unit-cell) */


    posuc.read(uc_file); // read the unitcell file

    scale_fac_uc = posuc.scale_fac();
    latt_vec_uc = posuc.latt_vec();
    Na_uc = posuc.Na();
    Naps_uc = posuc.Naps();
    LM_uc_onlycoord = posuc.LMcoord();
    LM_uc = posuc.LM();


    cout << "\n unitcell poscar file read \n" << endl;

    posuc.sort(); /* LM matrix is sorted according to OMEN conventions 
                   Note that if the coordinates are written in Direct
                  coordinates, they are converted in Cartesian*/
    
    LMs_uc = posuc.LMs();
    LMs_uc_ord = posuc.LMs_ord();

    Lxyz_uc(0,0) = latt_vec_uc(0,0)*scale_fac_uc;
    Lxyz_uc(0,1) = latt_vec_uc(1,1)*scale_fac_uc;
    Lxyz_uc(0,2) = latt_vec_uc(2,2)*scale_fac_uc;
    
    //The LMs_ord  matrix  is printed out in a file
     ofstream lms_uc("Output_data/LMs_uc_ord.dat");

    if (lms_uc.is_open())
      lms_uc<< LMs_uc_ord << endl;

    cout << "\n LM_uc sorted and written in LMs_uc_ord.dat \n" << endl;

  /*=======================================================================*/

  // FC_file FORCE_CONSTANTS 

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Reading FORCE CONSTANT file " << endl;
  cout << "-----------------------------------------------" << endl;

    ForceConstFile fc(Na,LMs_ord, mr); //constructor

    Eigen::MatrixXd FC_eVA; // force constants matrix : (eV/A^2)
    Eigen::MatrixXd FC; // force constants matrix : (J/m^2)
    Eigen::MatrixXd FCs; // sorted force constants matrix (J/m^2)
    Eigen::MatrixXd FCs_m; // sorted and mass scaled force constants matrix (J/m^2)

    fc.read(FC_file); // read FC_file
    FC_eVA =fc.FC_eVA();

    //The FCs matrix (eV/A^2) is printed out in a file
    ofstream f("Output_data/FC_me");
    
    if (f.is_open())
      f << fixed << setprecision(17) << endl;
      f<< FC_eVA << endl;

  /*---------------------------------------------------------------------------*/

    fc.convert(); //convert FC (eV/A^2) -> (J/m^2)
    FC=fc.FC();

  /*---------------------------------------------------------------------------*/

    fc.sort(); // sort FC 
    FCs = fc.FCs();

    cout << "sorted Force constant matrix -> written in file \n" << endl;
    //The FCs matrix (J/m^2) is printed out in a file
    ofstream fs("Output_data/FCs_me");
    
    if (fs.is_open())
      fs << fixed << setprecision(17) << endl;
      fs<< FCs << endl;
    

  /*---------------------------------------------------------------------------*/

      fc.scale(); //mass scaled FCs -> used to calculate phonon dispersion
      FCs_m = fc.FCs_m();
      

       cout << "scaled sorted Force constant matrix -> written in file" << endl;
      //The FCs matrix (J/m^2) is printed out in a file
             ofstream fsm("Output_data/FCs_m_me");     
      if (fsm.is_open())
        fsm << fixed << setprecision(17) << endl;
      fsm<< FCs_m << endl;
       

/* 
 #############################################################################

                      CUTTING FCs USING A CUT-OFF RADIUS

 #############################################################################

Function used : cut_FC

Be careful: the FC matrix to cut is FCs (sorted FC), 
not FCs_m (sorted and scaled by the mass FC matrix)
*/    

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Introducing CUT OFF RADIUS " << endl;
  cout << "-----------------------------------------------" << endl;


      Eigen::MatrixXd FCs_cut(3*Na,3*Na); //cut FCs
      
      FCs_cut = Eigen::MatrixXd::Zero(3*Na,3*Na); // inizialization FC_cut

      FCs_cut = cut_FC( LMs_ord, FCs, rc);
      
      //The FCs_cut matrix  is printed out in a file
      ofstream ffccut("Output_data/FCs_cut");     
      if (ffccut.is_open())
       {
          ffccut<< fixed << setprecision(17) << endl;
          ffccut<< FCs_cut  << endl;
          
         cout << "\n FCs_cut -> written in file" << endl;
       }



/* 
 #############################################################################

                     SYMMETRIZING FCs_cut

 #############################################################################

Function used : sym_FC
*/    

      Eigen::MatrixXd FCs_cut_sym(3*Na,3*Na); //sym FCs_cut
      
      FCs_cut_sym = Eigen::MatrixXd::Zero(3*Na,3*Na); // inizialization FCs_cut_sym

      FCs_cut_sym = sym_FC(  FCs_cut, Na);
      
      //The FCs_cut_sym  matrix  is printed out in a file
      ofstream ffcsym("Output_data/FCs_cut_sym");     
      if (ffcsym.is_open())
       {
          ffcsym<< fixed << setprecision(17) << endl;
          ffcsym<< FCs_cut_sym  << endl;
          
         cout << "\n FCs_cut_sym -> written in file" << endl;
       }




/* 
 #############################################################################

                        FIND THE UNIT CELL
                        (recursive process)

 #############################################################################

check how many neighbor slabs in x, y, and z directions need to
be included to fully cover the interactions in the force constant matrix. 
get as a return a (3,1)matrix nneigh:
(num. slabs in x-dir, num. slabs in y-dir, num. slabs in z-dir, )

Function used : check_slabs
*/  

  Eigen::MatrixXd nneigh(3,1); /*Defines how many slabs need to be included 
                                  in the corresponding direction to cover 
                                  all interactions. dim(3,1) 
                                  (num. slabs in x-dir, num. slabs in y-dir, 
                                  num. slabs in z-dir, )*/

  // if no redefinition is needed, we get nneigh = 1,1,1
  // the vector No_redefinition is used as check for the do..while loop
  Eigen::MatrixXd No_redefinition(3,1);  
  No_redefinition << 1,1,1;




  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "FIND UNIT CELL " << endl;
  cout << "-----------------------------------------------" << endl;




  nneigh=check_slabs(FCs_cut_sym, LMs_ord, LMs_uc_ord, Lxyz_uc, NFold, tol_lim );


  cout<<"\n number of  neighbor slabs in x, y, and z directions needed to be included" << endl;
  cout << " to fully cover the interactions in the force constant matrix\n" << nneigh << endl;

  if(nneigh!=No_redefinition)
    {
      cout << "\nThe unit cell is not big enough.\n Start redefinition...\n"<< endl;  
    }

do
{
  



/* 
 #############################################################################

                         EXTEND UC 

 #############################################################################

if in the previous section nneigh(i)>1 (i=x,y,z) -> the uc must be resized

1)
  Na_uc      -> Na_uc_new = nneigh(x)*nneigh(y)*nneigh(z)*Na_uc 

2)
  LMs_uc_ord -> LMs_uc_ord_new

3)Lxyz_uc = (Lx_uc,Ly_uc,Lz_uc) -> Lxyz_uc_new = (Lx_uc_new,Ly_uc_new,Lz_uc_new)

   where Li_uc_new = nneigh(i)*Li_uc

4) 
  NFold -> NFold/nneigh = NFold(i)/nneigh(i)

5) Check if the redefined unit cell is big enough (nneigh = 1,1,1)
   (otherwise, the operation is repeated)
*/


// ==========  1) ================================

int Naold;
  Naold = Na_uc; // keep trak of the original Na_uc
  Na_uc= Na_uc*nneigh.prod(); // redefine Na_uc

  cout<< "\n Na_uc_new :  " << Na_uc << endl;


// ==========  2) ===============================

  //enlarge the LMs_uc_ord matrix, we add Na_uc - Naold zero rows
      LMs_uc_ord.conservativeResizeLike(Eigen::MatrixXd::Zero(Na_uc,5));

  //Fill the zero rows with the atom in the enlarged cell

      for(int IAuc = 0; IAuc < Naold; IAuc++)
        { 
          for(int ix = 0; ix < nneigh(0); ix++)
            {
              for(int iy = 0; iy < nneigh(1); iy++)
                {
                  for(int iz =0; iz < nneigh(2); iz++)
                    {                      

                      LMs_uc_ord(IAuc+Naold,0) = LMs_uc_ord(IAuc,0) + ix*Lxyz_uc(0);
                      LMs_uc_ord(IAuc+Naold,1) = LMs_uc_ord(IAuc,1) + iy*Lxyz_uc(1);
                      LMs_uc_ord(IAuc+Naold,2) = LMs_uc_ord(IAuc,2) + iz*Lxyz_uc(2);
                      LMs_uc_ord(IAuc+Naold,3) = LMs_uc_ord(IAuc,3);
                    }
                }
            }
        }
     
// ==========  3) ================================

  for(int i = 0 ; i < 3; i++)
    {
      Lxyz_uc(0,i) = nneigh(i,0)*Lxyz_uc(0,i);
    }

  cout << "\n Lxyz_uc_new : " << Lxyz_uc << endl;

// ==========  4) ===============================

       
  for(int i = 0 ; i < 3; i++)
    {
      NFold(i) = NFold(i)/nneigh(i,0);
    }

  cout<< "\n NFold_new :  " << NFold << endl;



// ==========  5) ===============================

      //Check if the redefined unit cell is big enough (nneigh = 1,1,1)
      // (otherwise, the operation is repeated)

      nneigh=check_slabs(FCs_cut_sym, LMs_ord, LMs_uc_ord, Lxyz_uc, NFold, tol_lim );

      cout<<"\n number of  neighbor slabs in x, y, and z directions needed to be included" << endl;
      cout << " to fully cover the interactions in the force constant matrix\n" << nneigh << endl;
      

 }while(nneigh!=No_redefinition);

 cout << "\n \n...end redefinition\n"<< endl; 
 cout << "-----------------------------------------------" << endl;
 cout << "\nThe new unit cell is\n"<< endl;
 cout<< "\n Na_uc :  " << Na_uc << endl;
 cout <<  setprecision(17) <<"\n Lxyz_uc : " << Lxyz_uc << endl;
 cout<< "\n NFold :  " << NFold << endl;

 //The LMs_uc_ord(new)  matrix  is printed out in a file
 ofstream lms_nuc("Output_data/LMs_uc_ord_new.dat");
 
 if (lms_nuc.is_open())
   lms_nuc<< fixed << setprecision(17) << endl;
   lms_nuc<< LMs_uc_ord << endl;
 
 cout << "\n coordiantes of the new unit cell written in LMs_uc_ord_new.dat \n" << endl;

/* 
 #############################################################################

                       FIND NEAREST NEIGHBORS

 #############################################################################

class used: FindNN

Here the matrixes rvec_short and Nmulti (see below for details)
are created by the class  FindNN.

 */ 

  cout << "================================================" << endl;
  cout << "================================================" << endl;


      Eigen::MatrixXd rvec_short;/*Containing the bond vectors  and the length of the bond vector 
                                   between  an atom in the  unit-cell  and an atom in a supercell.
                                   dim (Na_uc * Na * (NSCELL_x*2+1) * (NSCELL_y*2+1) * (NSCELL_z*2+1),10) =
                                   = (atom in uc,  atom in sc, NSCELL_x , NSCELL_y , NSCELL_z ,
                                       x12, y12, z12, r12, nn_logical)
                                   where x12 = x2 - x1, y12 = y2 - y1, z12 = z2 - z1
                                         r12 = sqrt(x12*x12 + y12*y12 + z12*z12)  
                                         nn_logical is 1 if ( r12-(min dist 1-2) ) < compdiff, otherwise 0*/
     

      Eigen::MatrixXd Nmulti;/* Counts the closest neighbors considering all neighboring supercells
                                 Dimension (Na_uc,Na).
                              Is obtained summing up all column 10 of rvec_short, for fixed atmo1 and atom2*/


      FindNN nn(Na_uc, Na, LMs_ord, LMs_uc_ord, Lxyz, Lxyz_uc, ixyz1_ref, NSCELL, tol_lim);
      
      nn.find_nn();

      rvec_short=nn.rvec_short();
     

      //The rvec_short matrix  is printed out in a file
             ofstream fnn("Output_data/rvec_short");     
      if (fnn.is_open())
        //fnn << fixed << setprecision(17) << endl;
      fnn<<rvec_short  << endl;

       cout << "\n rvec_shor -> written in file" << endl;

       Nmulti = nn.Nmulti();
       
      //The Nmulti matrix  is printed out in a file
             ofstream fnm("Output_data/Nmult");     
      if (fnm.is_open())
      fnm<< Nmulti  << endl;

       cout << "\n Nmulti -> written in file" << endl;
       
       
/* 
 #############################################################################

                      REFERENCE PHONON BAND STRUCTURE

 #############################################################################

Function used: get_dispersion_direct

(nested function: get_k, getHtot, find_vec)

The phonon band structure is calculated using the FCs_m (sorted and scaled by the mass FC matrix)
via the function get_dispersion_direct.
 This will be the reference phonon band structure for the subsequent comparison (after the cut)

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
     where ntotk =(Nkx+1)+(Nky+1)+(Nkz+1)
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

 //) ===============  k-space grid generation ============================ 

Eigen::MatrixXd kx; //matrix dim(Nk(0)+1,1) = dim(Nkx+1,1)
Eigen::MatrixXd ky; //matrix dim(Nk(1)+1,1) = dim(Nky+1,1)
Eigen::MatrixXd kz; //matrix dim(Nk(2)+1,1) = dim(Nkz+1,1)

kx =getk(Nk(0));
ky =getk(Nk(1));
kz =getk(Nk(2));
      

 //) ===============  compute H and eigenvalues ============================ 

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Computing PHONON REFERNCE DISPERSION " << endl;
  cout << "-----------------------------------------------" << endl;


       int indE = (Nk(0)+1)+(Nk(1)+1)+(Nk(2)+1);
       /* if in one direction the k vectory is empty (= just 1 point 
          equal to 0.0, the E columns dimension is reduced by 1, since the
          eigenvalues for ka=(0.0) will not be computed)*/
       if(Nk(0)+1==1) indE = indE -1;
       if(Nk(1)+1==1) indE = indE -1;
       if(Nk(2)+1==1) indE = indE -1;

       Eigen::MatrixXd E_test(3*Na_uc,indE); /* Eigenvlues matrix 
                                                 = phonon energies for each ka(i)*/
       

       E_test = get_dispersion_direct( FCs_m, NFold,kx,ky,kz, indE,
                                       Na, Na_uc, LMs_ord, LMs_uc_ord,
                                       Lxyz_uc,rvec_short, Nmulti, 
                                       tol_lim, hbar, charge);


      //The E_test  matrix  is printed out in a file
             ofstream fet("Output_data/Etest.dat");     
      if (fet.is_open())
        {
          fet << fixed << setprecision(17) << endl;
          fet<< E_test  << endl;
      
          cout << "\n E_test -> written in file" << endl;
        }

/* 
 #############################################################################

                        EXTRACT SUBMATRICES (TENSOR)

 #############################################################################
this function calculates  matrices in OMEN style, 
using the updated structure information from check_slabs function 
and extend uc section 
(new unit cell = Na_uc_new, Lxyz_uc_new, LMs_uc_new)

Note that we are studying transport properties -> 
at least trasport direction (x) must be repeated = NFold_bs(0)=3

input: FCs_cut_sym, LMs_ord, LMs_uc_new, Lxyz,Lxyz_uc_new, NFold_bs

output : Hs_mod Tensor (NFold_bs(1), NFold_bs(2)*3*Na_uc_new, NFold_bs(1)*3*Na_uc_new) 
         it is composed by
         NFold_bs(0) matrices = 3 matrices: H00,H10,H01 
         matrices dimensions: ( NFold_bs(2)*3*Na_uc_new, NFold_bs(1)*3*Na_uc_new)
         elements : FCs_cut_sym


*/

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "\nComputing tensor for the phonon dispersion " << endl;
  cout << "-----------------------------------------------" << endl;


  Eigen::Tensor<double, 3> Hs_mod(NFold_bs(0),3*Na_uc*NFold_bs(1),NFold_bs(2)*3*Na_uc);
    
    Hs_mod.setZero();


     Hs_mod = extract_submatrices(FCs_cut_sym,Na, Na_uc, LMs_ord, LMs_uc_ord, Lxyz,
                                 Lxyz_uc, NFold_bs);

    //The Hs_mod  tensor  is printed out in a file
     ofstream hs("Output_data/Hs_mod.dat");

    if (hs.is_open())
      hs<< fixed << setprecision(17) << endl;
      hs<< Hs_mod << endl;

    cout << "\n Hs_mod is  written in Hs_mod.dat \n" << endl;
    


/* 
 #############################################################################

                         ASR + SYMMETRIZING + SCALE  

 #############################################################################
In this section we  
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

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Applying ASR " << endl;
  cout << "-----------------------------------------------" << endl;

    Hs_mod = asr_sym_scale( Hs_mod, NFold_bs, Na_uc, LMs_uc_ord,mr);

    int Ny,Nz;
    Ny = (int)NFold_bs(1);
    Nz = (int)NFold_bs(2);


    Eigen::MatrixXd asrHs0(3*Na_uc*Ny,3*Na_uc*Nz); 
    asrHs0.setZero(3*Na_uc*Ny,3*Na_uc*Nz); 
     Eigen::MatrixXd asrHs1(3*Na_uc*Ny,3*Na_uc*Nz); 
    asrHs1.setZero(3*Na_uc*Ny,3*Na_uc*Nz); 
     Eigen::MatrixXd asrHs2(3*Na_uc*Ny,3*Na_uc*Nz); 
    asrHs2.setZero(3*Na_uc*Ny,3*Na_uc*Nz); 


    for (int IA1 = 0; IA1 < Na_uc; IA1 ++)
      {
        for (int IA2 = 0; IA2 < Na_uc; IA2 ++)
          {
            for(int r= 0; r < 3; r++)
              {
                for(int c= 0; c < 3; c++)
                  {  

                   
                    for(int iyt =0; iyt < NFold_bs(1); iyt ++)
                      { 
                        for(int izt =0; izt < NFold_bs(2); izt ++)
                          {
                            
                                
                            asrHs0(3*IA1+r + 3*Na_uc*iyt,3*IA2+c + 3*Na_uc*izt) = 
                              Hs_mod(0,3*IA1+r + 3*Na_uc*iyt,3*IA2+c + 3*Na_uc*izt);
                            
                            asrHs1(3*IA1+r + 3*Na_uc*iyt,3*IA2+c + 3*Na_uc*izt) = 
                              Hs_mod(1,3*IA1+r + 3*Na_uc*iyt,3*IA2+c + 3*Na_uc*izt);
                            
                            asrHs2(3*IA1+r + 3*Na_uc*iyt,3*IA2+c + 3*Na_uc*izt) = 
                              Hs_mod(2,3*IA1+r + 3*Na_uc*iyt,3*IA2+c + 3*Na_uc*izt);
                            
                          }
                      }
                    
                  }
              }
          }
      }


      ofstream hsasr0("Output_data/Hs0.dat");

    if (hsasr0.is_open())
      hsasr0<< fixed << setprecision(17) << endl;
    hsasr0<< asrHs0<<endl;


      ofstream hsasr1("Output_data/Hs1.dat");

    if (hsasr1.is_open())
      hsasr1<< fixed << setprecision(17) << endl;
    hsasr1<< asrHs1<<endl;


      ofstream hsasr2("Output_data/Hs2.dat");

    if (hsasr2.is_open())
      hsasr2<< fixed << setprecision(17) << endl;
    hsasr2<< asrHs2<<endl;


 
    
    cout << "\n ASR + symmetrization + scale  \n -> Hs_mod is  written in Hs_mod_asr.dat \n" << endl;


/* 
 #############################################################################

                     PHONON BAND STRUCTURE

 #############################################################################

Function used: get_dispersion

(nested function:  getHtotOMEN)

The phonon band structure is calculated using the Hs_mod (asr, sym and scaling imposed)
via the function get_dispersion_direct.
 This will be campared with E_test

1st) the k-space grid of a rectangular brillouin zone is generated
     by the function get_k(Na). 
     Input:   Nka: number k-points in ka direction (int)  (a = x,y,z)
     Output: ka: vector containing the Nka ka-points between -pi and pi     
     nb. actually ka vector is a  (Nka,1) matrix
     
     NB. we already computed the grid! see section 
         REFERENCE PHONON BAND STRUCTURE

2nd) for each k-point ka(i) a matrix Hamiltonian is created H(ka(i))
     by the function getHtotOMEN.

3rd) the corresponding eigenvalues E(ka(i)) are calculated using the 
     Eigen attribute .eigenvalues()
     all the  E(ka(i)) are stored in matrix E -> dim(E)=(3*Na_uc,ntotk)
     where ntotk =(Nkx+1)+(Nky+1)+(Nkz+1)
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

  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "Computing PHONON DISPERSION " << endl;
  cout << "-----------------------------------------------" << endl;

       Eigen::MatrixXd E(3*Na_uc,indE); /* Eigenvlues matrix 
                                                = phonon energies for each ka(i)*/
       

       E = get_dispersion( Hs_mod, NFold_bs, kx,ky,kz, 
                                indE,Na_uc, hbar, charge);


      //The E  matrix  is printed out in a file
       ofstream fe("Output_data/E.dat");     
      if (fe.is_open())
        {
          fe << fixed << setprecision(17) << endl;
          fe<< E << endl;
      
          cout << "\n E -> written in file" << endl;
        }
           
/* 
 #############################################################################

                   COMPARISON  PHONON BAND STRUCTURES

 #############################################################################
*/

// ================  plot ==================================================

      /* This part only works if one chose Nk = (M,0,0) 
       it generates the data file to be used in the 
       already preprared gnuplot file "Eplot.gnu"*/
      
      if(Nk(0)!= 0 && Nk(1)==0 && Nk(2) == 0)
        {
          Eigen::MatrixXd Eplot(3*Na_uc*((int)Nk(0)+1),3);
          Eplot.setZero(3*Na_uc*(Nk(0)+1),3);
          
          
          for(int i=0; i < 3*Na_uc; i++)
            {
              for(int j=0; j < ((int)Nk(0)+1); j++)
                {
                  Eplot(j+i*((int)Nk(0)+1),0)= kx(j,0);
                  Eplot(j+i*((int)Nk(0)+1),1)= E_test(i,j);
                  Eplot(j+i*((int)Nk(0)+1),2)= E(i,j);
                  
                }
            }
          
          //The Eplot  matrix  is printed out in a file
          ofstream fep("Output_data/plotE/Eplot.dat");     
          if (fep.is_open())
            {
              fep << fixed << setprecision(17) << endl;
              fep<< Eplot << endl;
              
          cout << "\n Eplot -> written in file" << endl;
            }
          
        }

// ================   difference ============================================

      Eigen::MatrixXd Ecomp(3*Na_uc,indE); /* (Etest - E)/E */


      for(int i=0; i < 3*Na_uc; i++)
        {
          for(int j=0; j < indE; j++)
            {
              Ecomp(i,j) = (E_test(i,j)- E(i,j))/E(i,j);
            }
        }

      std::ptrdiff_t im, jm;// index of the maximum diff

      cout << setprecision(2)<< "\n\nThe maximum difference between E_test and E is  " <<abs(Ecomp.maxCoeff(&im,&jm))*100<< "%"<< endl;
      cout << "at position " << im << "," <<jm << endl;

      cout << "\nExcluding the acoustic phonons "<< endl;
      

      Ecomp.row(0).setZero();
      Ecomp.row(1).setZero();
      Ecomp.row(2).setZero();



      std::ptrdiff_t im2, jm2;// index of the maximum diff no acoustic ph

      cout << setprecision(2)<< "The maximum difference between E_test and E is " <<abs(Ecomp.maxCoeff(&im2,&jm2))*100<< "%"<< endl;
      cout << "at position " << im2 << "," <<jm2 << endl;
 
  cout << "================================================" << endl;
  cout << "================================================" << endl;
  cout << "================== J O B =======================" << endl;
  cout << "============ C O M P L E T E D==================" << endl;
  cout << "================================================" << endl;
  cout << "============  G O O D B Y E  ===================" << endl;
  cout << "================================================" << endl;







  return 0;
}

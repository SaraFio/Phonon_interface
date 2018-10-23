
#include "PoscarFile.hpp"




/* ________________________________________________________________________

              C L A S S : PoscarFile

   ________________________________________________________________________
*/
// Constructor
PoscarFile::PoscarFile(int nas) : nas_(nas) {}

//#############################################################################

//read
void PoscarFile::read(std::string filename_POSCAR)
{
  /* open the file with the data
  ================================*/
  std::ifstream is_poscar(filename_POSCAR);
    if (!is_poscar)
    {
      std::cout << "Unable to open file: " << filename_POSCAR << std::endl;
      exit(1);
    }
    // ----------------------------------------------------------------

    /* read the name of the compound 
    ================================*/
    is_poscar >> compoundname_;

    // ----------------------------------------------------------------

    /* read the scale factor 
    ================================*/
    is_poscar >> scalefac_;

    // ----------------------------------------------------------------

    /* read the lattice vector matrix 
    ================================*/
    lattvec_.resize(3,3); // set the dim of lattvec

    for (int i = 0; is_poscar && i < 3; ++i)
      for (int j = 0; j < 3; ++j)
          is_poscar >> lattvec_(i,j); 

    // ----------------------------------------------------------------

    /* read the atomic species name (dummy - not used)
    ================================*/
    is_poscar >> atomicspecies_;

    // ----------------------------------------------------------------

    /* read the Number of Atomos Per atomic Species matrix 
    ======================================================*/
    Na_ = 0;
    Naps_.resize(1,nas_); // set the dim of Naps

        for(int j = 0;  j < nas_; j++)
          {
              is_poscar >> Naps_(0,j);
              Na_ += Naps_(0,j); // sum over all the atom per atomic species
          }


    // ----------------------------------------------------------------


    /* read if the coordinates are given in direct or cartesian units
    ===================================================================*/
        is_poscar >> coordunit_;

        char direct[] = "Direct";
        char cartesian[] = "Cartesian";

        if(strcmp(coordunit_,direct)!=0 || strcmp(coordunit_,direct)!=0)
          {
            std::cout << "Not valid option. Direct or Cartesian coordinates" << std::endl;
            exit(1);
          }
       
       
    // ----------------------------------------------------------------

    /* read the matrix with the coordinates of the atoms
    =====================================================
   1) First a matrix LM_onlycoord(Na,3) is created 
   and filled with the coordinates of the atoms
   directly read from the POSCAR file
   2)  Then a new matrix LM(Na,4) is created
   col 1-2-3 of LM = col 1-2-3 LM_onlycoord
   col 4 pf LM contains the atomic species 
   -> LM = (x-coord, y-coord,z-coord, atomic species)
    Remember that atomic species are labelled by indeces.
    eg. 2 atomic species -> can be (x-coord, y-coord,z-coord, 1)
                            or     (x-coord, y-coord,z-coord, 2)
    Note that the POSCAR file has been written using 
    first all the coordinates of the atomic species 1,
    then all the coordinates of the atomic species 2,
    ...
    up to all the coordinates of the atomic species n_as*/ 

    // (1)
    LM_onlycoord_.resize(Na_,3); // set the dim of LM_onlycoord

    for (int i = 0; is_poscar && i < Na_; ++i)
      for (int j = 0; j < 3; ++j)
        is_poscar >> LM_onlycoord_(i,j);

    //(2)
    LM_.resize(Na_,4); // set the dim of LM

    int start = 0;
    int end = Naps_(0,0);

    for (int k = 0; k<nas_; k++)
      {
        for (int i = start; i < end; i++)
          {
            LM_(i,3) = k+1;
          }
        start = Naps_(0,k) ;
        if(k!=(nas_-1))
          end = Naps_(0,k)+Naps_(0,k+1);
          } 

   

      /*---------- Direct -> Cartesian--------------*/

        if(strcmp(coordunit_,direct)==0)
          {
            std::cout << "changing from direct to cartesian" << std::endl;

               LM_onlycoord_ =  scalefac_*LM_onlycoord_*lattvec_;
          }

        
        /*-------------------------------------------*/
        
         for (int i = 0;  i < Na_; ++i)
        for (int j = 0; j < 3; ++j)
          LM_(i,j) = LM_onlycoord_(i,j);    


}

//#############################################################################

void PoscarFile::sort(){

  LMs_.resize(Na_,4); // set the dim of sorted LM

  LMs_ = LM_; // The original LM is copied, in order to be sorted
 
  Eigen::MatrixXd tmp(1,4); // dummy matrix used to swap


   bool swapped = false;
   int xx = Na_;      // <=== added a variable for the end of the scan
   // Here the row are sorted (ascending order) with respect to x (1st column)
   do
   {
      swapped = false; 
 
      for (int j = 1; j < xx; j++) //  <=== note xx not x 
      {
        if (LMs_(j-1,0) > LMs_(j,0))  
         { 
           tmp = LMs_.row(j);
           LMs_.row(j) =  LMs_.row(j-1);
            LMs_.row(j-1) = tmp;
            swapped = true;
         }
      }
   xx--; // <= end element finalised in position, so don't need to go as far next time
   } while ( swapped );
   

   // if 2 x values are identical -> they are sorted with respect to y
          for(int i =0; i < Na_ ; i++)
            {
              for (int j=i+1; j < Na_; j++)
                {
                  if(LMs_(i,0)==LMs_(j,0) && i!=j && LMs_(i,1)>LMs_(j,1))
                    {
                      tmp = LMs_.row(i);
                      LMs_.row(i) = LMs_.row(j);
                      LMs_.row(j) = tmp;
                    }
                }
            }


          //  if 2 (x,y) values are identical -> they are sorted with respect to z
         for(int i =0; i < Na_ ; i++)
            {
              for (int j=i+1; j < Na_; j++)
                {
                  if(LMs_(i,0)==LMs_(j,0)&& LMs_(i,1)==LMs_(j,1)  && i!=j && LMs_(i,2)>LMs_(j,2))
                    {
                      tmp = LMs_.row(i);
                      LMs_.row(i) = LMs_.row(j);
                      LMs_.row(j) = tmp;
                    }
                }
            }     


 

  LMs_ord_.resize(Na_,5);//set the dim of LMs_ord

  //the first 4 col of LMs_ord are the same of LMs
   for(int i =0; i < Na_ ; i++)
    {
      for (int j=0; j < 4; j++)
        {
          LMs_ord_(i,j) = LMs_(i,j);
        }
    }

   /* the 5th col of LMs_ord is the one of the 
      correspondance LM-LMs*/
  for(int i =0; i < Na_ ; i++)
    {
      for (int j=0; j < Na_; j++)
        {
          if(LMs_.row(i)==LM_.row(j))
            {
              LMs_ord_(i,4) = j;

            }
        }
      
    }

}

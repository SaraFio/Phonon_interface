
#include "ForceConstFile.hpp"


/* ________________________________________________________________________

              C L A S S : ForceConstFile

   ________________________________________________________________________
*/
// Constructor
ForceConstFile::ForceConstFile(int Na, Eigen::MatrixXd LMs_ord, Eigen::MatrixXd mr): Na_sc_(Na), LMs_ord_(LMs_ord), mr_(mr) {}

//##########################################################################
 
//read
void ForceConstFile::read(std::string FC_file)
{
      

  /* open the file with the data
  ================================*/
  std::ifstream is_fc(FC_file);
    if (!is_fc)
    {
      std::cout << "Unable to open file: " << FC_file << std::endl;
      exit(1);
    }
    // ----------------------------------------------------------------

    if (!(is_fc >> Na_ )  )
      std::cout << "Error reading Na" << std::endl;

    if(Na_!=Na_sc_)
      {
        std::cout << "different Na in sc_file and FC_file" << std::endl;
        exit(1);
      }

    /* local variables to fill FC 
     ============================*/

  int Na_tot = Na_*Na_; // is used to read all the combination between the Na atoms
  Eigen::MatrixXd dummy(1,2); // here will be temporary stored the atom indexes of IFC (eg 1 2)
  Eigen::MatrixXd tmp(3,3); // here are temporary stored the IFC for atom i and j 
  int c = 0; // this constant is used to fill the FC matrix
  int d = 0; // this constant is used to fill the FC matrix   
  // ----------------------------------------------------------------


  FC_eVA_.resize(3*Na_,3*Na_); //set dim of FC : (eV/A^2)

  // ----------------------------------------------------------------

  /*   fill the FC matrix 
    ==============================================================================*/

    for(int k=0; k < Na_tot; k++)
      {
       
        /* the dummy variables ot the atoms
           are stored in the matrix dummy
           and deleted
         */
        for(int j = 0; j < 2; j++)
          {
            is_fc >> dummy(0,j);
          }

        /* The FC matrix is here filled*/
        for(int i = 0; i<3; i++)
          {
            for(int j= 0; j < 3; j++)
              {
                is_fc >> tmp(i,j);

                FC_eVA_(i+c,j+d)=tmp(i,j); // FC (eV/A^2)
              }
          }

        d+=3; // in this way we move to the next block (in line)
 
        // in this way we move a block below
        if((k+1)%48==0)
          {
            c+=3;
            d = 0;
          }

      }

}



//##########################################################################
 
//convert the  FC matrix (eV/A^2) --> (J/m^2) 
void ForceConstFile::convert(){

  double q;
  double m_scale;
  q = 1.6022e-19;
  m_scale = 1.0/(6.02214e26);

  FC_.resize(3*Na_,3*Na_); // set dim of FC (J/m^2)

  FC_ = FC_eVA_*q/((1e-10)*(1e-10))*(1/m_scale); // FC (J/m^2)

}

//##########################################################################
 
//sort  FC matrix according to OMEN conventions 
void ForceConstFile::sort(){
  
  FCs_.resize(3*Na_,3*Na_); // set dim of FCs (J/m^2)

    for(int i = 0; i< Na_; i++)
      {
        for(int j= 0; j < Na_; j++)
          {
            for(int l = 0; l < 3; l++)
              {
                for(int m = 0; m < 3 ; m++)
                  {
                    int a, b,c,d;
                    a = 3*LMs_ord_(i,4) + l;
                    b = 3*LMs_ord_(j,4) + m;
                    c = 3*i+l;
                    d = 3*j+m;
                    
                    FCs_(c,d) = FC_(a,b);
                    

                  }
              }
          }
      }

}


//##########################################################################
 
//scale  FCs matrix 

  /* We want to scale FCs by the square root of the corresponding masses
     -> FCs(i,j)/sqrt(m_i,m_j)
     the element FCs(i,j) is related to certain atom i and j 
     so m_i correspond to the mass of the atomic species 
     of atom i
     and the same for m_j
*/

void ForceConstFile::scale(){

   FCs_m_.resize(3*Na_,3*Na_); // set dim of FCs_m (J/m^2)

    /* local variables to fill FC 
     ============================*/
  int Im1 = 0; /* this index reads in LMs_ord matrix the atomic 
                  species index, such as 1, 2, .. , n_as
                  it is used to assign the value of the mass
                  Note that atomic species indeces starts from 1,
                  while indeces in matrices from zero.*/
  int Im2 = 0; /* this index reads in LMs_ord matrix the atomic 
                  species index
                  it is used to assign the value of the mass.*/

  double m1=0.0; /* mass of the atomic species i*/
  double m2=0.0; /* mass of the atomic species j*/

  double denom=0.0 ; /*scale factor sqrt(m1m2)*/

  // ----------------------------------------------------------------

 
    for(int i = 0; i< Na_; i++)
      {
        for(int j= 0; j < Na_; j++)
          {
            for(int l = 0; l < 3; l++)
              {
                for(int m = 0; m < 3 ; m++)
                  {
                    int a,b;
                    a = 3*i+l;
                    b = 3*j+m;
                     Im1 = LMs_ord_(i,3);
                    Im2 = LMs_ord_(j,3);
                    m1 = mr_(Im1-1);
                    m2 = mr_(Im2-1);
                    denom = sqrt(m1*m2);
                    
                    //std::cout << a  << "\t" << Im1<< "\t" << m1 << "\t"<< b << "\t" << Im2 << "\t" << m2 << "\t" << denom<< std::endl;
                 
                    FCs_m_(a,b) = FCs_(a,b)/denom;
                  }
              }
          }
      }

}

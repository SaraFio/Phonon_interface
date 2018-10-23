
#include "MaterialStructureFile.hpp"


/* ________________________________________________________________________

              C L A S S : MaterialStructureFile

   ________________________________________________________________________
*/
// Constructor
MaterialStructureFile::MaterialStructureFile(int nas) : nas_(nas) {}

//read
void MaterialStructureFile::read(std::string ms_file) {
      // open the file with the data
  std::ifstream is_matstr(ms_file);
    if (!is_matstr)
    {
      std::cout << "Unable to open file: " <<  ms_file<< std::endl;

      exit(1);
    }

    // ----------------------------------------------------------------

    //header: constants
    is_matstr >> dummy_matstr_; 

    // ----------------------------------------------------------------

    // h bar
    is_matstr >> dummy_matstr_;

    if (!(is_matstr >> hbar_ )  )
      {
      std::cout << "ERROR reading h bar" << std::endl;
      exit(1);
      }
      
      // ----------------------------------------------------------------

    //  charge
    is_matstr >> dummy_matstr_;

    if (!(is_matstr >> charge_ )  )
      {
      std:: cout << "ERROR reading charge" << std::endl;
      exit(1);
      }

    // ----------------------------------------------------------------

    //  mass scale factor
    is_matstr >> dummy_matstr_;

    if (!(is_matstr >> mscalefac_ )  )
      {
      std::cout << "ERROR reading mass scale factor" << std::endl;
      exit(1);
      }

    // ----------------------------------------------------------------

    //  numerical criteria  for tolerance
    is_matstr >> dummy_matstr_;

    if (!(is_matstr >> tollim_ )  )
      {
      std::cout << "ERROR reading numerical criteria for tolerance" << std::endl;
      exit(1);
      }

    // ----------------------------------------------------------------

    //header: data
    is_matstr >> dummy_matstr_; 

    // ----------------------------------------------------------------
    
    //  reduced masses of atomic species
    is_matstr >> dummy_matstr_;

    mr_.resize(1,nas_); // assign the matrix dimensions

      for(int j = 0;  j < nas_; j++)
          {
            is_matstr >> mr_(0,j);
          }
      

    // ----------------------------------------------------------------

    //  full  masses of atomic species

    mt_.resize(1,nas_); // assign the matrix dimensions

      for(int j = 0;  j < nas_; j++)
          {
            mt_(0,j) =  mscalefac_ *mr_(0,j);
          }

      

    // ----------------------------------------------------------------

    // NFold

    is_matstr >> dummy_matstr_;

    NFold_.resize(1,3); // assign the matrix dimensions

      for(int j = 0;  j < 3; j++)
          {
            is_matstr >> NFold_(0,j);
          }

    // ----------------------------------------------------------------


    // NFold_bs

    is_matstr >> dummy_matstr_;
    Eigen::MatrixXd bulk(1,3);
    Eigen::MatrixXd well(1,3);
    Eigen::MatrixXd wire(1,3);
    NFoldbs_.resize(1,3);  // assign the matrix dimensions

      bulk << 3,3,3;
      well << 3,1,3;
      wire << 3,1,1;

      for(int j = 0;  j < 3; j++)
          {
            is_matstr >> NFoldbs_(0,j);
          }
	
	
      
      if( (NFoldbs_ != bulk) && (NFoldbs_ != well) && (NFoldbs_ !=wire))
        {
          std::cout << "ERROR: the transort direction is not repeated" << std::endl;
          std::cout << " you are not simulating transport." << std::endl;
          std::cout << "Please assign NFold_bs = [3 x x] " << std::endl;
          std::cout << "Remember \n [3 3 3] = bulk \n [3 1 3]= 2D\n [3 1 1] = nanowire" << std::endl;
          exit(1);
        }


    // ----------------------------------------------------------------

 
    // ixyz1_ref

    is_matstr >> dummy_matstr_;
     ixyz1_ref_.resize(1,3); // assign the matrix dimensions

      for(int j = 0;  j < 3; j++)
          {
            is_matstr >> ixyz1_ref_(0,j);
          }
    // ----------------------------------------------------------------
    // NSCELL

    is_matstr >> dummy_matstr_;
     NSCELL_.resize(1,3); // assign the matrix dimensions

      for(int j = 0;  j < 3; j++)
          {
            is_matstr >> NSCELL_(0,j);
          }
    // ----------------------------------------------------------------

    // Nk

    is_matstr >> dummy_matstr_;
     Nk_.resize(1,3); // assign the matrix dimensions

      for(int j = 0;  j < 3; j++)
          {
            is_matstr >> Nk_(0,j);
          }
    // ----------------------------------------------------------------

    //  radius cutoff
    is_matstr >> dummy_matstr_;

    if (!(is_matstr >> rc_ )  )
      std::cout << "Error reading radius cutoff" << std::endl;


}


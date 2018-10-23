/* ==========================================================================

              FUNCTION :  find_vec

=============================================================================
looks if a row defined by vec exists in a matrix and return the index of the row of the matrix
Note that vec is a (1,n)matrix, and matrix is a (m,n) matrix
*/


#include "find_vec.hpp"

int find_vec(Eigen::MatrixXd vec, Eigen::MatrixXd matrix, double comp_diff){


  int N; // number of column of the matrix
  N=matrix.cols();

  int M; // number of rows of the matrix
  M=matrix.rows();

  /* check if the number columns of the matrix and of vec os equal*/
  if(vec.rows()!=1) 
    {
      std::cout << "\n ERROR: vec must be a (1,n) matrix. you gave : vec (m,n), with m!=1" << std::endl;
    
      exit(1);
    }

  /* check if the number columns of the matrix and of vec os equal*/
  if(N!=vec.cols()) 
    {
      std::cout << "\n ERROR: vec must be a (1,n) matrix and matrix must be a (m,n) matrix. you gave : vec (1,n) and matrix (m,N)" << std::endl;
    
      exit(1);
    }

  int index = -1; // this is the index that will be returned
  bool found = 0; /* if =1 : value found -> return
                     if =0 : No value found -> error
                     if >1 : more than 1 value found -> error*/ 
  
  for(int i = 0; i < M ; i++ )
    {
      int j=0;
      while((abs(matrix(i,j)-vec(0,j))<comp_diff)&&(j<N-1))
        {
          j++;
          /*In this way if matrix(i,0)==vec(0,0) it proceeds to check matrix(i,1) and vec(0,1). 
            if matrix(i,1)==vec(0,1) it proceeds to check matrix(i,2) and vec(0,2)
            and so on, up to j = N-1. note that j=N-1 will only be reached if all the elements are equal*/

        }
      if ((j==N-1)&&(abs(matrix(i,N-1)-vec(0,N-1))<comp_diff)) // if matrix(i,j)==vec(0,j), for all j
        {
          /* if this condition has never been fulfilled before,
           index is set equal to 1 and found is now true (=1)*/

          if (found==0)
            {
              index =i;
              found = 1;
              
            }
          /*if  more than 1 row fulfil the condition (thus found is already true, found=1): error!*/
          else
            {
              std::cout << "\n ERROR: More than one value found in the find_vec function" << std::endl;
              exit(1);
            }
          
        }
      
    }
  /*if  more than any row fulfils the condition (thus found is still  false, found=0): error!*/  
  if(found==0) 
    {
       std::cout << "\n ERROR: Anything found in the find_vec function" << std::endl;
       exit(1);
    }

  /*if only 1 row fulfils the condition (thus found is true, found=1): return index*/  
  else return index;
  
}


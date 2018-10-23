#ifndef INPUTFILE_tosort_H
#define INPUTFILE_tosort_H


#include "InputFile.hpp"

class InputFile_tosort{
public:
  virtual void sort() =0; 
  /* the matrix is read and sorted according to OMEN conventions
OMEN requires coordintes in ascending order (smallest one = 1st, biggest one = last)
NB. coordinates are given as (x,y,z)
if there are identical x-values -> coord are sorted with respect to y-values
if there are identical x- and y-values -> coord are sorted with respect to z-values*/

};

#endif

#ifndef INPUTFILE_H
#define INPUTFILE_H

#include<string>

class InputFile{
public:
  virtual void read(std::string) =0; // the file is read and all the infos stored

};

#endif

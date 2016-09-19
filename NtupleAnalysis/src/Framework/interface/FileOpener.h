#ifndef FileOpener_h
#define FileOpener_h

// User
#include "../../Auxiliary/src/Datasets.C"

class FileOpener 
{

 public:
  TChain* OpenFile(const std::string SamplePath, const std::string SampleName, TChain* tree = 0);
    
 private:
  Datasets datasets_;

};

#endif // FileOpener_h

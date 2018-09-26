#ifndef FileOpener_h
#define FileOpener_h

// User
#include "../../Auxiliary/src/Datasets.C"

class FileOpener 
{

 public:
  TChain* OpenFiles(const std::string multicrabPath, const std::string dataset, TChain* tree = 0);
  std::string GetFirstFileName(const std::string multicrabPath, const std::string dataset);
  vector<string> GetListOfFiles(const string inDir);
    
 private:
  Datasets datasets_;

};

#endif // FileOpener_h

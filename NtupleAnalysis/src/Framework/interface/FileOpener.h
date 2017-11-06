#ifndef FileOpener_h
#define FileOpener_h

// User
#include "../../Auxiliary/src/Datasets.C"

class FileOpener 
{

 public:
  TChain* OpenFiles(const std::string multicrabPath, const std::string dataset, TChain* tree = 0);
  vector<string> GetListOfFiles(const string inDir, const string filePrefix = "raw2TTree_");
  bool EndsWith(const std::string &str, const std::string &suffix);
    
 private:
  Datasets datasets_;

};

#endif // FileOpener_h

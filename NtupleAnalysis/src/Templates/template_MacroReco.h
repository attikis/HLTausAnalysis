#ifndef MACRONAMERECO_h
#define MACRONAMERECO_h

#include <iostream>

#include "../utilities/TreeAnalyserReco.C"
#include "../utilities/AuxTools.C"

class MACRONAMERECO : public TreeAnalyserReco
{
  public:
     MACRONAMERECO(const std::string SampleName, 
	   const std::string text_, 
	   const int maxEvents_ = -1, 
	   TTree* tree=0) : 
     TreeAnalyserReco("MACRONAMERECO", SampleName, text_, maxEvents_, tree) {};
     virtual void Loop();

  private:
     AuxTools tools;
    // Add your private variables/methods here
};

#endif // MACRONAMERECO_h

#ifndef MACRONAMEMC_h
#define MACRONAMEMC_h

#include <iostream>

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/AuxTools.C"
#include "../utilities/Table.C"
#include "../utilities/Datasets.C"

class MACRONAMEMC : public TreeAnalyserMC
{
  public:
    MACRONAMEMC(const std::string SamplePath,
	const std::string SampleName,
	const std::string text_, 
	const int maxEvents_ = -1, 
	TTree* tree=0) : 
    TreeAnalyserMC("MACRONAMEMC", SamplePath, SampleName, text_, maxEvents_, tree) {
      auxTools_.StopwatchStart();
      mcSample = SampleName;
    };

    ~MACRONAMEMC() {};

    virtual void Loop();

  private:
    C11Func c11;
    std::string mcSample;
    AuxTools auxTools_;
    Datasets datasets_;

    // Add your private variables/methods here
};

#endif // MACRONAMEMC_h

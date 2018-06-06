#ifndef TreeAnalyserMC_h
#define TreeAnalyserMC_h

// User
#include "../../Auxiliary/src/MCTools.C"
#include "../src/TreeAnalyserBase.C"
#include "../src/TreeReaderMC.C"
#include "../src/ObjectMCAssociator.C"

class L1Tracks;

// class TrackingParticles;

class TreeAnalyserMC : public TreeAnalyserBase, public TreeReaderMC, 
                       public virtual MCTools,  public virtual ObjectMCAssociator
{
 public: 
   TreeAnalyserMC() {
     std::cout << " TreeAnalyserMC : The constructor needs some arguments \n";
   };
   TreeAnalyserMC(const std::string MyName_, 
		  const std::string SamplePath, 
		  const std::string SampleName, 
		  const std::string text_, 
		  const int maxEvents_ = -1, 
		  //TTree* tree=0) : 
		  TChain* chain=0) : 
                  TreeAnalyserBase(MyName_, SamplePath, SampleName, text_, maxEvents_), 
                  //TreeReaderMC(SamplePath, SampleName, tree) { InitSelector(); };
		  TreeReaderMC(SamplePath, SampleName, chain) { InitSelector(); };

   friend class L1Tracks;
   // friend class TrackingParticles;

 private:
   void InitSelector();
   L1Tracks* s;
   //TrackingParticles* tp;

};

#endif // TreeAnalyserMC_h

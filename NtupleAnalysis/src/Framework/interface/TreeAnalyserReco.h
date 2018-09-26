#ifndef TreeAnalyserReco_h
#define TreeAnalyserReco_h

// User
#include "../src/TreeAnalyserBase.C"
#include "../src/TreeReaderReco.C"

class L1Tracks;
class L1PixelTrackFit;

class TreeAnalyserReco : public TreeAnalyserBase, public TreeReaderReco
{
 public:
    TreeAnalyserReco() { 
      std::cout << " TreeAnalyserReco : Constructor needs some arguments \n"; 
    }; 
    TreeAnalyserReco(const std::string MyName_, 
		     const std::string SampleName, 
		     const std::string text_, 
		     const int maxEvents_ = -1,
		     //TTree* tree=0) : 
                     TChain* chain=0) :
                     TreeAnalyserBase(MyName_, SampleName, text_, maxEvents_), 
		     //TreeReaderReco(SampleName, tree) { InitSelector(); };
		     TreeReaderReco(SampleName, chain) { InitSelector(); };
  
    friend class L1Tracks;
    friend class L1PixelTrackFit;

 private:
    void InitSelector();
    L1Tracks* s;
    L1PixelTrackFit* f;
    
};

 
#endif // TreeAnalyserFull_h

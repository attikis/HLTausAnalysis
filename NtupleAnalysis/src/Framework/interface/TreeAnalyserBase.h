#ifndef TreeAnalyserBase_h
#define TreeAnalyserBase_h

// ROOT
#include <TFile.h>

class TreeAnalyserBase
{
 public: 
  TreeAnalyserBase() {
    std::cout << " TreeAnalyserBase : The constructor needs some arguments \n";
  };
  TreeAnalyserBase(const std::string MyName_, 
		   const std::string SamplePath_, 
		   const std::string SampleName_, 
		   const std::string Text_, 
		   const int MaxEvents_ = -1) {
    Initialize(MyName_, SamplePath_, SampleName_, Text_, MaxEvents_);
  };
  ~TreeAnalyserBase();
  
  // Variables
  TFile* outFile;
  int MaxEvents;

 private:
  void Initialize(const std::string MyName_, 
		  const std::string SamplePath_, 
		  const std::string SampleName_, 
		  const std::string Text_,
		  const int MaxEvents_);

  // Variables
  // TFile* outFile;
  // int MaxEvents;
  std::string SamplePath;
  std::string SampleName;
  std::string MyName;
  std::string Text;
};

#endif // TreeAnalyserBase_h

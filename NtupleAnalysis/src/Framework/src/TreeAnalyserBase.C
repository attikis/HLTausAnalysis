#ifndef TreeAnalyserBase_cxx
#define TreeAnalyserBase_cxx

// System
#include <typeinfo>

// User
#include "../interface/TreeAnalyserBase.h"

void TreeAnalyserBase::Initialize(const std::string MyName_, 
				  const std::string SamplePath_, 
				  const std::string SampleName_, 
				  const std::string Text_, 
				  const int MaxEvents_)
{

  MyName     = MyName_; // passed from TreeAnalyserMC(, in each analyzer header file
  SamplePath = SamplePath_;
  SampleName = SampleName_;
  MaxEvents  = MaxEvents_;
  Text       = Text_;

  // Create filename
  std::string prefix;
  std::string postfix;
  std::string outFileName;
  if (MyName == "") prefix = "";
  else prefix = MyName + "_";

  if (Text_ == "") postfix = "";
  else postfix = "_" + Text_;

  // Create the output file name
  outFileName = prefix + "histograms-" + SampleName_ + postfix + ".root";
  outFile     = new TFile(outFileName.c_str(),"Recreate");

  return;
}

TreeAnalyserBase::~TreeAnalyserBase()
{
  outFile->Close();
}

#endif // TreeAnalyserBase_cxx

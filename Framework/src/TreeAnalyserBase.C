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

  MyName     = MyName_;
  SamplePath = SamplePath_;
  SampleName = SampleName_;
  MaxEvents  = MaxEvents_;
  Text       = Text_;

  std::string outFileName;
  if (Text == "") outFileName = MyName + "_Histograms_" + SampleName + ".root";
  else outFileName = MyName + "_Histograms_" + SampleName + "_" + Text + ".root";

  outFile = new TFile(outFileName.c_str(),"Recreate");
}

TreeAnalyserBase::~TreeAnalyserBase()
{
  outFile->Close();
}

#endif // TreeAnalyserBase_cxx

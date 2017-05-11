//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage:
// .x run.cc(<multicrab_dir>, <sample_alias>, <text_to_append_to_output>, <numberOfEvents>);
//
// Example:
// root -l
// root[0] .x run.cc("multicrab_CaloTk_v910p2_test", "SingleTau-PU140", "", -1);
//
// WARNING: Does NOT work with ROOT v6! Works fine with ROOT5
// source /Users/attikis/ROOT/v5-34-00-patches/bin/thisroot.csh
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#include "../Tracking.h"
// R__LOAD_LIBRARY(../Tracking.C)

//R__LOAD_LIBRARY(../Tracking.C)
void run(const std::string MulticrabDir = "", 
	 const std::string SampleName = "", 
	 const std::string text = "", 
	 const int maxEvents = -1)
{
    
  //gSystem->CompileMacro("../Tracking.C");

  // const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/P2L1T_HLTaus_91X";
  const std::string absolutePath = "/afs/cern.ch/user/a/attikis/workspace/multicrab/";
  Tracking macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

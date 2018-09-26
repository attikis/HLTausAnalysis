//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage:
// .x run.cc(<multicrab_dir>, <sample_alias>, <text_to_append_to_output>, <numberOfEvents>);
//
// Example:
// root -l
// root[0].x run.cc("multicrab_CaloTk_v910p2_20170613T2048", "SingleTau_PU140", "", -1)
// root[0].x run.cc("multicrab_CaloTk_v910p2_test", "SingleTau-PU140", "", -1);
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../Sanity.C+" // how to compile macro in ROOT6 (compatible with ROOT5 as well)

void run(const std::string MulticrabDir = "", 
	 const std::string SampleName = "", 
	 const std::string text = "", 
	 const int maxEvents = -1)
{

  const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/P2L1T_HLTaus_91X";
  // const std::string absolutePath = "/afs/cern.ch/user/a/attikis/workspace/multicrab/";
  // const std::string absolutePath = "/Users/attikis/disk/hltaus/rootFiles/TTrees/P2L1T_HLTaus_91X/";
  
  Sanity macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

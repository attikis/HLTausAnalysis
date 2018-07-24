//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage:
// .x run.cc(<multicrab_dir>, <sample_alias>, <text_to_append_to_output>, <numberOfEvents>);
//
// Example:
// root -l
// root[0] .x run.cc("multicrab_CaloTk_v910p2_20170520T1958", "SingleTau-PU140", "", -1);
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../CaloTk.C+" // how to compile macro in ROOT6 (compatible with ROOT5 as well)

void run(const std::string MulticrabDir = "", 
	 const std::string SampleName = "", 
	 const std::string text = "", 
	 const int maxEvents = -1)
{

  // Alexandros' files
  // const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/P2L1T_HLTaus_91X";
  // const std::string absolutePath = "/afs/cern.ch/user/a/attikis/workspace/multicrab/";
  // const std::string absolutePath = "/Users/attikis/disk/hltaus/rootFiles/TTrees/P2L1T_HLTaus_91X/";

  // Marina's files 
  //const std::string absolutePath = "/afs/cern.ch/user/m/mtoumazo/workspace/multicrab";
  const std::string absolutePath = "/eos/user/m/mtoumazo";
  
  // Mikko's files
  //const std::string absolutePath = "/eos/user/m/mlotti";
  
  // Santeri's files
  // const std::string absolutePath = "/afs/cern.ch/work/s/slaurila/public/HLTaus/CMSSW_9_1_0_pre2/src/HLTausAnalysis/NtupleAnalysis/src/Tracking/work";
  
  CaloTk macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

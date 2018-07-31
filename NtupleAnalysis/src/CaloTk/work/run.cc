//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage:
// .x run.cc(<multicrab_dir>, <sample_alias>, <text_to_append_to_output>, <numberOfEvents>);
//
// Example:
// root -l
// root[0] .x run.cc("multicrab_HLTaus_v1015_20180710T1650", "SingleNeutrino_L1TPU200", "", -1)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../CaloTk.C+" // how to compile macro in ROOT6 (compatible with ROOT5 as well)

void run(const std::string MulticrabDir = "", 
	 const std::string SampleName = "", 
	 const std::string text = "", 
	 const int maxEvents = -1)
{

  // Files on LXPLUS
  const std::string absolutePath = "/eos/user/m/mtoumazo";
  //const std::string absolutePath = "/eos/user/m/mlotti";
  
  CaloTk macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

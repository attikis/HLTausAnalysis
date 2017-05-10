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
void run(const std::string MulticrabDir = "", 
	 const std::string SampleName = "", 
	 const std::string text = "", 
	 const int maxEvents = -1)
{
  
  gSystem->CompileMacro("../TTTracks.C");
  const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/P2L1T_HLTaus_91X";
  TTTracks macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

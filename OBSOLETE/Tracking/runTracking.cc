void runTracking(const std::string MulticrabDir = "", 
		 const std::string SampleName = "", 
		 const std::string text = "", 
		 const int maxEvents = -1)
{

  gSystem->CompileMacro("Tracking.C");

  const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v6";
  
  Tracking macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

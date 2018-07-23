// root[0] .x runKinematics.cc("TTracks_TTPixelTracks_CaloTausCorr_FitParam5_PtMin2_TTI2023Upg14D_TestVersionWithPixelBugs_28Jan2015_144735", "VBF", "all", -1);
// root[0] .x runKinematics.cc("TTracks_TTPixelTracks_CaloTausCorr_FitParam5_PtMin2_TTI2023Upg14D_TestVersionWithPixelBugs_28Jan2015_144735", "VBF", "1p", -1);
// root[0] .x runKinematics.cc("TTracks_TTPixelTracks_CaloTausCorr_FitParam5_PtMin2_TTI2023Upg14D_TestVersionWithPixelBugs_28Jan2015_144735", "VBF", "3p", -1);
void runKinematics(const std::string MulticrabDir = "", 
		   const std::string SampleName = "", 
		   const std::string text = "", 
		   const int maxEvents = -1)
{

  gSystem->CompileMacro("TauTriggerKinematics.C");

  const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v3";

  TauTriggerKinematics macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage:
// .x runPixelTracking.cc(<multicrab_dir>, <sample_alias>, <text_to_append_to_output>, <numberOfEvents>);
//
// Example-1:
// root -l
// root[0] .x runPixelRefitting.cc ("L1PixTks_CandPixHits_TPs_GenPs_NPixHitsMin3_rphi3mm_rz5mm_PrivateProduction_29July2015", "SingleMuon_NoPU", "", 1);
// (Run on all events to produce the file "PixelRefitting_Histograms_SingleMuon_NoPU.root")
//
// Example-2:
// root -l
// root[0] .x runTauTrigger.cc("L1PixTks_CandPixHits_TPs_GenPs_NPixHitsMin3_rphi3mm_rz5mm_TTI2023Upg14D_24Jul2015_154331", "VBF", "alex", 10);
// (Run on 10 events to produce the file "PixelRefitting_Histograms_SingleMuon_NoPU_alex.root")
//
// Multicrab directories previously used:
// "AllTracks_CaloTausCorr_FitParam4_TTI2023Upg14D_17Feb2015_165244"
// "AllTracks_CaloTausCorr_FitParam5_TTI2023Upg14D_16Feb2015_113104"
//
// WARNING: Does NOT work with ROOT v6! Works fine with ROOT v5-34-00-patches/
// source /Users/attikis/ROOT/v5-34-00-patches/bin/thisroot.csh
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void runPixelTracking(const std::string MulticrabDir = "", 
		      const std::string SampleName = "", 
		      const std::string text = "", 
		      const int maxEvents = -1)
{
  
  gSystem->CompileMacro("PixelTracking.C");

  const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v6";
  
  PixelTracking macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

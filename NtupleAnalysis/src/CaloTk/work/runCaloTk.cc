//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage:
// .x runCaloTk.cc(<multicrab_dir>, <sample_alias>, <text_to_append_to_output>, <numberOfEvents>);
//
// Example-1:
// root -l
// root[0] .x runCaloTk.cc("L1CaloTaus_CaloCorr_TTTracks_Stubs_TTPixelTracks_CandPixHits_TPs_GenPs_v620SLHC12p1_07Nov2016", "VBF", "", -1);
// (Run on all events to produce the file "CaloTk_Histograms_VBF.root")
//
// Example-2:
// root -l
// root[0] .x runCaloTk.cc("L1PixTks_CandPixHits_TPs_GenPs_NPixHitsMin3_rphi3mm_rz5mm_PrivateProduction_29July2015", "VBF", "alex", 10);
// (Run on 10 events to produce the file "CaloTk_Histograms_VBF_alex.root")
//
// root -l
// root[0] .x runCaloTk.cc("L1CaloTaus_CaloCorr_TTTracks_Stubs_TTPixelTracks_CandPixHits_TPs_GenPs_v620SLHC12p1_07Nov2016/", "VBF", ", 10);
//
// WARNING: Does NOT work with ROOT v6! Works fine with ROOT v5-34-00-patches/
// source /Users/attikis/ROOT/v5-34-00-patches/bin/thisroot.csh
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void runCaloTk(const std::string MulticrabDir = "", 
		   const std::string SampleName = "", 
		   const std::string text = "", 
		   const int maxEvents = -1)
{

  gSystem->CompileMacro("../CaloTk.C");

  // const std::string absolutePath = "/Users/attikis/my_work/cms/lxplus/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v3";
  // const std::string absolutePath = "/Users/attikis/my_work/cms/lxplus/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v6";
  // const std::string absolutePath = "/Users/attikis/my_work/cms/lxplus/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v7";
  // const std::string absolutePath = "/Users/attikis/disk/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v7";
  const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v7";
  // const std::string absolutePath = "/Users/attikis/disk/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v7";

   
  CaloTk macro(absolutePath + "/" + MulticrabDir, SampleName, text, maxEvents);
  macro.Loop();
}

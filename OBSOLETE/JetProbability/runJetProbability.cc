void runJetProbability(const std::string MulticrabDir = "", 
                    const std::string SampleName = "", 
                    const std::string text = "",
                    const int MaxEvents = -1)
{
  gSystem->CompileMacro("JetProbability.C");

  //const std::string absolutePath = "/home/athermw/Research/rootFiles";
  const std::string absolutePath = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v5";
  JetProbability macro(absolutePath + "/" + MulticrabDir,SampleName, text, MaxEvents);
  macro.Loop();
}

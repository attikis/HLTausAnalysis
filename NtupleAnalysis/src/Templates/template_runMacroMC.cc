void runMACRONAMEMC(const std::string MulticrabDir = "", 
                    const std::string SampleName = "", 
                    const std::string text = "",
                    const int MaxEvents = -1)
{
  gSystem->CompileMacro("MACRONAMEMC.C");
  const std::string absolutePath = "/raid2/data/cms/TauSLHC/NtupleSamples";
  MACRONAMEMC macro(absolutePath + "/" + MulticrabDir,SampleName, text, MaxEvents);
  macro.Loop();
}

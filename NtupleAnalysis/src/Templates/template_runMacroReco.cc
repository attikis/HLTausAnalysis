void runMACRONAMERECO(const std::string SampleName = "", const std::string text = "", const int maxEvents = -1)
{
  gSystem->CompileMacro("MACRONAMERECO.C");
  MACRONAMERECO macro(SampleName, text, maxEvents);
  macro.Loop();
}

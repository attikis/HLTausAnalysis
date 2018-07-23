void runXSections(void)
{

  gSystem->CompileMacro("XSections.C");
  
  int Precision     = 6;     // for converting numbers to strings
  double Higgs_Mass = 125.0; // in GeV
  double Int_Lumi   = 100;   // in fb
  bool bLatexStyle  = false;
  XSections xs(Precision, Higgs_Mass, Int_Lumi);
 
  xs.PrintHiggsFinalStates("ttH", "HToWW"    , bLatexStyle);
  xs.PrintHiggsFinalStates("ttH", "HToZZ"    , bLatexStyle);
  xs.PrintHiggsFinalStates("ttH", "HToTauTau", bLatexStyle);

  return;
}

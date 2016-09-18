#ifndef XSections_cxx
#define XSections_cxx

#include "XSections.h"
#include "../utilities/constants.h"

//****************************************************************************
XSections::XSections(int Precision,
		     const double Higgs_Mass,
		     const double Int_Lumi)
//****************************************************************************
{
  _Precision  = Precision;  // for converting numbers to strings
  _Int_Lumi   = Int_Lumi;   // in fb-1
  _Higgs_Mass = Higgs_Mass; // in GeV
  _CheckInputParams();
  _InitVars();
  
}


//****************************************************************************
XSections::~XSections()
//****************************************************************************
{

}


//****************************************************************************
void XSections::_CheckInputParams(void)
//****************************************************************************
{

  std::vector<double> v_Higgs_Masses_Allowed;
  // v_Higgs_Masses_Allowed.push_back(123.0); // xSections missing
  // v_Higgs_Masses_Allowed.push_back(123.5); // xSections missing
  // v_Higgs_Masses_Allowed.push_back(124.0); // xSections missing
  // v_Higgs_Masses_Allowed.push_back(124.5); // xSections missing
  v_Higgs_Masses_Allowed.push_back(125.0);
  v_Higgs_Masses_Allowed.push_back(125.5);
  v_Higgs_Masses_Allowed.push_back(126.0);
  // v_Higgs_Masses_Allowed.push_back(126.5); // xSections missing
  // v_Higgs_Masses_Allowed.push_back(127.0); // xSections missing
  
  bool bFoundMatch = false;
  // For-loop: All higgs masses supported
  for (std::vector<double>::iterator it = v_Higgs_Masses_Allowed.begin(); it!=v_Higgs_Masses_Allowed.end(); it++){
    if ( _Higgs_Mass == (*it) ) bFoundMatch = true;
  }
  
  if(!bFoundMatch){
    cout << "E R R O R ! XSections::XSections(...) - Unsupported value for _Higgs_Mass = \"" << _Higgs_Mass << "\" GeV selected. EXIT" << endl;
    auxTools.PrintVector(v_Higgs_Masses_Allowed, "\nCurrently supported values for mH (GeV):");    
    exit(1);
  }

  return;
}



//****************************************************************************
void XSections::_InitVars(void)
//****************************************************************************
{

  _IntLumi_Units  = "fb-1";
  _XSection_Units = "fb";
  _Mass_Units     = "GeV";
  _Width_Units    = "GeV";

  _SetProductionChannels();
  _SetDecayChannels();
  _SetStableFinalState();
  _SetHiggsWidthMaps();
  _SetHiggsDaughterDecayChannelsMap();
  _SetBranchingRatioMaps();     

  return;
}



//****************************************************************************
void XSections::_SetDecayChannels(void)
//****************************************************************************
{

  // W Boson
  _channels_WToAll.push_back("WToLNu");
  _channels_WToAll.push_back("WToTauNu");
  _channels_WToAll.push_back("WToHadrons");
  
  // Higgs Boson (Order matters!)
  _channels_HToAll.push_back("HToWW");
  _channels_HToAll.push_back("HToZZ");
  _channels_HToAll.push_back("HToTauTau");
  
  // Z Boson
  _channels_ZToAll.push_back("ZToHadrons");
  _channels_ZToAll.push_back("ZToLL");
  _channels_ZToAll.push_back("ZToTauTau"); 
  _channels_ZToAll.push_back("ZToNuNu");
  
  // Tau Lepton
  _channels_TauToAll.push_back("TauToHadrons");
  _channels_TauToAll.push_back("TauToLNuNu"); 

  // Tau Lepton
  _channels_Empty.push_back("-");

  return;
}



//****************************************************************************
void XSections::_SetProductionChannels(void)
//****************************************************************************
{

  // Higgs Production Boson
  _channels_HiggsProduction.push_back("ttH");

  return;
}



//****************************************************************************
void XSections::_SetHiggsWidthMaps(void)
//****************************************************************************
{

  m_Channel_To_BR m_Channel_To_BR_m123p0;
  m_Channel_To_BR m_Channel_To_BR_m123p5;
  m_Channel_To_BR m_Channel_To_BR_m124p0;
  m_Channel_To_BR m_Channel_To_BR_m124p5;
  m_Channel_To_BR m_Channel_To_BR_m125p0;
  m_Channel_To_BR m_Channel_To_BR_m125p5;
  m_Channel_To_BR m_Channel_To_BR_m126p0;
  m_Channel_To_BR m_Channel_To_BR_m126p5;
  m_Channel_To_BR m_Channel_To_BR_m127p0;
  
  // BR values copied from ............: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2
  // ttH XSection values copied from...: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeVa
  // Width is in GeV. XSection is in fb
  m_Channel_To_BR_m123p0["HToBB"]         = 6.07E-01;// +3.0	-3.1
  m_Channel_To_BR_m123p0["HToTauTau"]     = 6.63E-02;// +5.9	-5.8
  m_Channel_To_BR_m123p0["HToMuMu"]       = 2.31E-04;// +6.2	-6.0
  m_Channel_To_BR_m123p0["HToCC"]         = 3.07E-02;// +12.2	-12.2
  m_Channel_To_BR_m123p0["HToTT"]         = 0.00E+00;// +0.0	-0.0
  m_Channel_To_BR_m123p0["HToGluGlu"]     = 8.71E-02;// +10.4	-10.1
  m_Channel_To_BR_m123p0["HToGammaGamma"] = 2.27E-03;// +5.2	-5.0
  m_Channel_To_BR_m123p0["HToZGamma"]     = 1.36E-03;// +9.2	-9.0
  m_Channel_To_BR_m123p0["HToWW"]         = 1.83E-01;// +4.5	-4.4
  m_Channel_To_BR_m123p0["HToZZ"]         = 2.18E-02;// +4.5	-4.4
  m_Channel_To_BR_m123p0["Width"]         = 3.82E-03;// +4.1	-4.1
  m_Channel_To_BR_m123p0["XSection"]      = 0.0;
  m_Channel_To_BR_m123p0["-"]             = 1.0;
  
  m_Channel_To_BR_m123p5["HToBB"]         = 6.00E-01;//	+3.1	-3.1
  m_Channel_To_BR_m123p5["HToTauTau"]     = 6.56E-02;//	+5.8	-5.8
  m_Channel_To_BR_m123p5["HToMuMu"]       = 2.28E-04;//	+6.1	-6.0
  m_Channel_To_BR_m123p5["HToCC"]         = 3.03E-02;//	+12.2	-12.2
  m_Channel_To_BR_m123p5["HToTT"]         = 0.00E+00;//	+0.0	-0.0
  m_Channel_To_BR_m123p5["HToGluGlu"]     = 8.68E-02;//	+10.3	-10.1
  m_Channel_To_BR_m123p5["HToGammaGamma"] = 2.28E-03;//	+5.1	-5.0
  m_Channel_To_BR_m123p5["HToZGamma"]     = 1.41E-03;//	+9.1	-9.0
  m_Channel_To_BR_m123p5["HToWW"]         = 1.91E-01;//	+4.4	-4.3
  m_Channel_To_BR_m123p5["HToZZ"]         = 2.29E-02;//	+4.5	-4.3
  m_Channel_To_BR_m123p5["Width"]         = 3.88E-03;//	+4.1	-4.0  
  m_Channel_To_BR_m123p5["XSection"]      = 0.0;
  m_Channel_To_BR_m123p5["-"]             = 1.0;
  
  m_Channel_To_BR_m124p0["HToBB"]         = 5.92E-01;//	+3.1	-3.2
  m_Channel_To_BR_m124p0["HToTauTau"]     = 6.48E-02;//	+5.8	-5.7
  m_Channel_To_BR_m124p0["HToMuMu"]       = 2.25E-04;//	+6.1	-6.0
  m_Channel_To_BR_m124p0["HToCC"]         = 2.99E-02;//	+12.2	-12.2
  m_Channel_To_BR_m124p0["HToTT"]         = 0.00E+00;//	+0.0	-0.0
  m_Channel_To_BR_m124p0["HToGluGlu"]     = 8.65E-02;//	+10.3	-10.1
  m_Channel_To_BR_m124p0["HToGammaGamma"] = 2.28E-03;//	+5.1	-5.0
  m_Channel_To_BR_m124p0["HToZGamma"]     = 1.45E-03;//	+9.1	-8.9
  m_Channel_To_BR_m124p0["HToWW"]         = 1.98E-01;//	+4.4	-4.3
  m_Channel_To_BR_m124p0["HToZZ"]         = 2.40E-02;//	+4.4	-4.3
  m_Channel_To_BR_m124p0["Width"]         = 3.94E-03;//	+4.0	-4.0    
  m_Channel_To_BR_m124p0["XSection"]      = 0.0;
  m_Channel_To_BR_m124p0["-"]             = 1.0;
  
  m_Channel_To_BR_m124p5["HToBB"]         = 5.84E-01;//	+3.2	-3.2
  m_Channel_To_BR_m124p5["HToTauTau"]     = 6.39E-02;//	+5.8	-5.7
  m_Channel_To_BR_m124p5["HToMuMu"]       = 2.22E-04;//	+6.0	-5.9
  m_Channel_To_BR_m124p5["HToCC"]         = 2.95E-02;//	+12.2	-12.2
  m_Channel_To_BR_m124p5["HToTT"]         = 0.00E+00;//	+0.0	-0.0
  m_Channel_To_BR_m124p5["HToGluGlu"]     = 8.61E-02;//	+10.3	-10.0
  m_Channel_To_BR_m124p5["HToGammaGamma"] = 2.28E-03;//	+5.0	-4.9
  m_Channel_To_BR_m124p5["HToZGamma"]     = 1.49E-03;//	+9.1	-8.9
  m_Channel_To_BR_m124p5["HToWW"]         = 2.07E-01;//	+4.3	-4.2
  m_Channel_To_BR_m124p5["HToZZ"]         = 2.52E-02;//	+4.4	-4.2
  m_Channel_To_BR_m124p5["Width"]         = 4.00E-03;//	+4.0	-4.0
  m_Channel_To_BR_m124p5["XSection"]      = 0.0;
  m_Channel_To_BR_m124p5["-"]             = 1.0;
  
  m_Channel_To_BR_m125p0["HToBB"]         = 5.77E-01;//	+3.2	-3.3
  m_Channel_To_BR_m125p0["HToTauTau"]     = 6.32E-02;//	+5.7	-5.7
  m_Channel_To_BR_m125p0["HToMuMu"]       = 2.20E-04;//	+6.0	-5.9
  m_Channel_To_BR_m125p0["HToCC"]         = 2.91E-02;//	+12.2	-12.2
  m_Channel_To_BR_m125p0["HToTT"]         = 0.00E+00;//	+0.0	-0.0  
  m_Channel_To_BR_m125p0["HToGluGlu"]     = 8.57E-02;//	+10.2	-10.0
  m_Channel_To_BR_m125p0["HToGammaGamma"] = 2.28E-03;//	+5.0	-4.9
  m_Channel_To_BR_m125p0["HToZGamma"]     = 1.54E-03;//	+9.0	-8.8
  m_Channel_To_BR_m125p0["HToWW"]         = 2.15E-01;//	+4.3	-4.2
  m_Channel_To_BR_m125p0["HToZZ"]         = 2.64E-02;//	+4.3	-4.2
  m_Channel_To_BR_m125p0["Width"]         = 4.07E-03;//	+4.0	-3.9
  m_Channel_To_BR_m125p0["XSection"]      = 508.5;   //	+5.7	-9.3	+8.8	-8.8
  m_Channel_To_BR_m125p0["-"]             = 1.0;
  
  m_Channel_To_BR_m125p5["HToBB"]         = 5.69E-01;//	+3.3	-3.3
  m_Channel_To_BR_m125p5["HToTauTau"]     = 6.24E-02;//	+5.7	-5.6
  m_Channel_To_BR_m125p5["HToMuMu"]       = 2.17E-04;//	+6.0	-5.8
  m_Channel_To_BR_m125p5["HToCC"]         = 2.87E-02;//	+12.2	-12.2
  m_Channel_To_BR_m125p5["HToTT"]         = 0.00E+00;//	+0.0	-0.0  
  m_Channel_To_BR_m125p5["HToGluGlu"]     = 8.52E-02;//	+10.2	-9.9
  m_Channel_To_BR_m125p5["HToGammaGamma"] = 2.28E-03;//	+4.9	-4.8
  m_Channel_To_BR_m125p5["HToZGamma"]     = 1.58E-03;//	+8.9	-8.8
  m_Channel_To_BR_m125p5["HToWW"]         = 2.23E-01;//	+4.2	-4.1
  m_Channel_To_BR_m125p5["HToZZ"]         = 2.76E-02;//	+4.3	-4.1
  m_Channel_To_BR_m125p5["Width"]         = 4.14E-03;//	+3.9	-3.9
  m_Channel_To_BR_m125p5["XSection"]      = 502.7;   //	+5.7	-9.3	+8.8	-8.8
  m_Channel_To_BR_m125p5["-"]             = 1.0;
  
  m_Channel_To_BR_m126p0["HToBB"]         = 5.61E-01;//	+3.3	-3.4
  m_Channel_To_BR_m126p0["HToTauTau"]     = 6.15E-02;//	+5.6	-5.6
  m_Channel_To_BR_m126p0["HToMuMu"]       = 2.14E-04;//	+5.9	-5.8
  m_Channel_To_BR_m126p0["HToCC"]         = 2.83E-02;//	+12.2	-12.2
  m_Channel_To_BR_m126p0["HToTT"]         = 0.00E+00;//	+0.0	-0.0  
  m_Channel_To_BR_m126p0["HToGluGlu"]     = 8.48E-02;//	+10.1	-9.9
  m_Channel_To_BR_m126p0["HToGammaGamma"] = 2.28E-03;//	+4.9	-4.8
  m_Channel_To_BR_m126p0["HToZGamma"]     = 1.62E-03;//	+8.9	-8.8
  m_Channel_To_BR_m126p0["HToWW"]         = 2.31E-01;//	+4.1	-4.1
  m_Channel_To_BR_m126p0["HToZZ"]         = 2.89E-02;//	+4.2	-4.0
  m_Channel_To_BR_m126p0["Width"]         = 4.21E-03;//	+3.9	-3.8
  m_Channel_To_BR_m126p0["XSection"]      = 496.6;   //	+5.7	-9.3	+8.8	-8.8
  m_Channel_To_BR_m126p0["-"]             = 1.0;
  
  m_Channel_To_BR_m126p5["HToBB"]         = 5.53E-01;//	+3.4	-3.4
  m_Channel_To_BR_m126p5["HToTauTau"]     = 6.08E-02;//	+5.6	-5.5
  m_Channel_To_BR_m126p5["HToMuMu"]       = 2.11E-04;//	+5.9	-5.7
  m_Channel_To_BR_m126p5["HToCC"]         = 2.79E-02;//	+12.2	-12.2
  m_Channel_To_BR_m126p5["HToTT"]         = 0.00E+00;//	+0.0	-0.0    
  m_Channel_To_BR_m126p5["HToGluGlu"]     = 8.42E-02;//	+10.1	-9.8
  m_Channel_To_BR_m126p5["HToGammaGamma"] = 2.28E-03;//	+4.8	-4.7
  m_Channel_To_BR_m126p5["HToZGamma"]     = 1.66E-03;//	+8.8	-8.7
  m_Channel_To_BR_m126p5["HToWW"]         = 2.39E-01;//	+4.1	-4.0
  m_Channel_To_BR_m126p5["HToZZ"]         = 3.02E-02;//	+4.1	-4.0
  m_Channel_To_BR_m126p5["Width"]         = 4.29E-03;//	+3.8	-3.8  
  m_Channel_To_BR_m126p5["XSection"]      = 0.0;
  m_Channel_To_BR_m126p5["-"]             = 1.0;
  
  m_Channel_To_BR_m127p0["HToBB"]         = 5.45E-01;//	+3.4	-3.5
  m_Channel_To_BR_m127p0["HToTauTau"]     = 5.99E-02;//	+5.5	-5.5
  m_Channel_To_BR_m127p0["HToMuMu"]       = 2.08E-04;//	+5.8	-5.7
  m_Channel_To_BR_m127p0["HToCC"]         = 2.75E-02;//	+12.2	-12.2
  m_Channel_To_BR_m127p0["HToTT"]         = 0.00E+00;//	+0.0	-0.0  
  m_Channel_To_BR_m127p0["HToGluGlu"]     = 8.37E-02;//	+10.1	-9.8
  m_Channel_To_BR_m127p0["HToGammaGamma"] = 2.28E-03;//	+4.8	-4.7
  m_Channel_To_BR_m127p0["HToZGamma"]     = 1.70E-03;//	+8.8	-8.7
  m_Channel_To_BR_m127p0["HToWW"]         = 2.48E-01;//	+4.0	-4.0
  m_Channel_To_BR_m127p0["HToZZ"]         = 3.15E-02;//	+4.1	-3.9
  m_Channel_To_BR_m127p0["Width"]         = 4.37E-03;//	+3.8	-3.7
  m_Channel_To_BR_m127p0["-"]             = 1.0;
  
  // typedef map<double, m_ChannelVsBR>
  m_HiggsMass_To_BR[123.0] = m_Channel_To_BR_m123p0;      
  m_HiggsMass_To_BR[123.5] = m_Channel_To_BR_m123p5;
  m_HiggsMass_To_BR[124.0] = m_Channel_To_BR_m124p0;
  m_HiggsMass_To_BR[124.5] = m_Channel_To_BR_m124p5;
  m_HiggsMass_To_BR[125.0] = m_Channel_To_BR_m125p0;
  m_HiggsMass_To_BR[125.5] = m_Channel_To_BR_m125p5;
  m_HiggsMass_To_BR[126.0] = m_Channel_To_BR_m126p0;
  m_HiggsMass_To_BR[126.5] = m_Channel_To_BR_m126p5;
  m_HiggsMass_To_BR[127.0] = m_Channel_To_BR_m127p0;

  return;
}


//****************************************************************************
void XSections::_SetHiggsDaughterDecayChannelsMap(void)
//****************************************************************************
{

  // m_Higgs_To_DecayModes["HToBB"]         = ;
  // m_Higgs_To_DecayModes["HToMuMu"]       = ;
  // m_Higgs_To_DecayModes["HToCC"]         = ;
  // m_Higgs_To_DecayModes["HToTT"]         = ;
  // m_Higgs_To_DecayModes["HToGluGlu"]     = ;
  // m_Higgs_To_DecayModes["HToGammaGamma"] = ;
  // m_Higgs_To_DecayModes["HToZGamma"]     = ;
  m_Higgs_To_DecayModes["HToTauTau"]     = GetTauDecayChannels();
  m_Higgs_To_DecayModes["HToWW"]         = GetWDecayChannels();
  m_Higgs_To_DecayModes["HToZZ"]         = GetZDecayChannels();
  m_Higgs_To_DecayModes["-"]             = GetEmptyDecayChannels();
  
  return;
}


//****************************************************************************
bool XSections::IsStableFinalState(string finalState)
//****************************************************************************
{

  // Ensure final state is allowed first.
  bool bNotFound = m_FS_To_IsStable.find(finalState) == m_FS_To_IsStable.end();

  if (bNotFound)
    {
      std::cout << "E R R O R ! XSections::IsStableFinalState(...) - Unknown final state \"" << finalState << "\". EXIT" << std::endl;
      std::cout << "m_FS_To_IsStable["<<finalState<<"] = " << m_FS_To_IsStable[finalState] << std::endl;
      exit(1);
    }
  
  return m_FS_To_IsStable[finalState];
}



//****************************************************************************
void XSections::_SetStableFinalState(void)
//****************************************************************************
{

  m_FS_To_IsStable["HToBB"]         = true;
  m_FS_To_IsStable["HToTauTau"]     = true;
  m_FS_To_IsStable["HToMuMu"]       = true;
  m_FS_To_IsStable["HToCC"]         = true;
  m_FS_To_IsStable["HToTT"]         = false;
  m_FS_To_IsStable["HToGluGlu"]     = true;
  m_FS_To_IsStable["HToGammaGamma"] = true;
  m_FS_To_IsStable["HToZGamma"]     = false;
  m_FS_To_IsStable["HToWW"]         = false;
  m_FS_To_IsStable["HToZZ"]         = false;

  m_FS_To_IsStable["WToENu"]        = true;
  m_FS_To_IsStable["WToMuNu"]       = true;
  m_FS_To_IsStable["WToLNu"]        = true;
  m_FS_To_IsStable["WToTauNu"]      = false;
  m_FS_To_IsStable["WToHadrons"]    = true;

  m_FS_To_IsStable["ZToEE"]         = true;
  m_FS_To_IsStable["ZToMuMu"]       = true;
  m_FS_To_IsStable["ZToLL"]         = true;
  m_FS_To_IsStable["ZToTauTau"]     = false;
  m_FS_To_IsStable["ZToNuNu"]       = true;
  m_FS_To_IsStable["ZToHadrons"]    = true;
  m_FS_To_IsStable["ZToBB"]         = true;
  
  m_FS_To_IsStable["TauToENuNu"]       = true;
  m_FS_To_IsStable["TauToMuNuNu"]      = true;
  m_FS_To_IsStable["TauToLNuNu"]       = true;
  m_FS_To_IsStable["TauToPiNu"]        = true;
  m_FS_To_IsStable["TauToPiPizeroNu"]  = true;
  m_FS_To_IsStable["TauToPi2PizeroNu"] = true;
  m_FS_To_IsStable["TauToPi3PizeroNu"] = true;
  m_FS_To_IsStable["TauToHadrons"]     = true;

  return;
}



//****************************************************************************
void XSections::_SetBranchingRatioMaps(void)
//****************************************************************************
{


  BRContainer = m_HiggsMass_To_BR[_Higgs_Mass];

  // http://pdg.lbl.gov/2012/listings/rpp2012-list-w-boson.pdf
  BRContainer["WToENu"]        = 0.1075;
  BRContainer["WToMuNu"]       = 0.1057;
  BRContainer["WToLNu"]        = 0.1075 + 0.1057;
  BRContainer["WToTauNu"]      = 0.125;
  BRContainer["WToHadrons"]    = 0.6760; 
  BRContainer["WToAll"]        = 2.085; // GeV

  // http://pdg.lbl.gov/2012/listings/rpp2012-list-z-boson.pdf
  BRContainer["ZToEE"]         = 0.03363;
  BRContainer["ZToMuMu"]       = 0.03366;
  BRContainer["ZToLL"]         = 0.03363 + 0.03366;
  BRContainer["ZToTauTau"]     = 0.03370;
  BRContainer["ZToNuNu"]       = 0.20000;
  BRContainer["ZToHadrons"]    = 0.6991;
  BRContainer["ZToBB"]         = 0.1512; 
  BRContainer["ZToAll"]        = 2.4952; // GeV

  // http://pdglive.lbl.gov/Particle.action?node=S035
  BRContainer["TauToENuNu"]       = 0.1783;
  BRContainer["TauToMuNuNu"]      = 0.1741;
  BRContainer["TauToLNuNu"]       = 0.1783 + 0.1741;
  BRContainer["TauToPiNu"]        = 0.1083;
  BRContainer["TauToPiPizeroNu"]  = 0.2552;
  BRContainer["TauToPi2PizeroNu"] = 0.0930;
  BRContainer["TauToPi3PizeroNu"] = 0.0105; 
  BRContainer["TauToHadrons"]     = 0.6476;
  BRContainer["TauToAll"]         = -999.9; // GeV

  // Other
  BRContainer["-"] = 1.0;
  _Higgs_XSection  = BRContainer["XSection"];
    
  return;
}


//****************************************************************************
void XSections::_ConvertChannelToFinalState(void)
//****************************************************************************
{

  m_Channel_To_FS["HToBB"]         = "b,b";
  m_Channel_To_FS["HToTauTau"]     = "Tau,Tau";
  m_Channel_To_FS["HToMuMu"]       = "Mu,Mu";
  m_Channel_To_FS["HToCC"]         = "c,c";
  m_Channel_To_FS["HToTT"]         = "t,t";
  m_Channel_To_FS["HToGluGlu"]     = "g,g";
  m_Channel_To_FS["HToGammaGamma"] = "gamma,gamma";
  m_Channel_To_FS["HToZGamma"]     = "Z,gamma";
  m_Channel_To_FS["HToWW"]         = "W,W";
  m_Channel_To_FS["HToZZ"]         = "Z,Z";

  m_Channel_To_FS["WToENu"]        = "e,nu";
  m_Channel_To_FS["WToMuNu"]       = "mu,nu";
  m_Channel_To_FS["WToLNu"]        = "l,nu";
  m_Channel_To_FS["WToTauNu"]      = "tau,nu";
  m_Channel_To_FS["WToHadrons"]    = "q,q";

  m_Channel_To_FS["ZToEE"]         = "e,e";
  m_Channel_To_FS["ZToMuMu"]       = "mu,mu";
  m_Channel_To_FS["ZToLL"]         = "l,l";
  m_Channel_To_FS["ZToTauTau"]     = "tau,tau";
  m_Channel_To_FS["ZToNuNu"]       = "nu,nu";
  m_Channel_To_FS["ZToHadrons"]    = "q,q";
  m_Channel_To_FS["ZToBB"]         = "b,b";

  m_Channel_To_FS["TauToENuNu"]    = "e,nu";
  m_Channel_To_FS["TauToMuNuNu"]   = "mu,nu";
  m_Channel_To_FS["TauToLNuNu"]    = "l,nu";
  m_Channel_To_FS["TauToHadrons"]  = "tau_h";

  return;
}



//****************************************************************************
vector<string> XSections::GetDecayChannels(string motherParticle)
//****************************************************************************
{

  vector<string> decayChannels;
  if (motherParticle.compare("W")        == 0) decayChannels = GetWDecayChannels();
  else if (motherParticle.compare("Z")   == 0) decayChannels = GetZDecayChannels();
  else if (motherParticle.compare("H")   == 0) decayChannels = GetHDecayChannels();
  else if (motherParticle.compare("Tau") == 0) decayChannels = GetTauDecayChannels();
  else
    {
      std::cout << "E R R O R ! XSections::GetDecayChannels(...) - The mother particle selected \"" << motherParticle << "\""
		<< " is not a supported. EXIT" << std::endl;
      exit(1);
    }

  return decayChannels;
}



//****************************************************************************
void XSections::PrintHiggsFinalStates(string HiggsProductionMode,
				      string HiggsDecayMode,
				      bool bLatex)
//****************************************************************************
{

  // Check that Higgs production mode is supported
  IsAllowedHiggsProduction(HiggsProductionMode);

  // Check that Higgs decay mode is supported
  IsAllowedHiggsDecay(HiggsDecayMode);

  // Get the final states table
  Table events_FS = GetAssociatedHiggsToXTable(HiggsDecayMode, bLatex);
  events_FS.Print();

  // Get the final states table (with identical final states merged)
  Table events_FS_merged = _ConvertToFinalStateTable(events_FS, bLatex);
  events_FS_merged.Print();
    
  return;
}


//****************************************************************************
void XSections::IsAllowedHiggsDecay(string HiggsDecayMode)
//****************************************************************************
{
  
  bool bIsAllowed = find(_channels_HToAll.begin(), _channels_HToAll.end(), HiggsDecayMode) != _channels_HToAll.end();

  if(!bIsAllowed)
    {
      std::cout << "E R R O R ! XSections::IsAllowedHiggsDecay(...) - The decay mode selected \"" << HiggsDecayMode << "\""
		<< " is not a supported Higgs decay mode. EXIT" << std::endl;
      exit(1);
    }
    
  return; 
}



//****************************************************************************
void XSections::IsAllowedHiggsProduction(string HiggsProductionMode)
//****************************************************************************
{
  
  bool bIsAllowed = find(_channels_HiggsProduction.begin(), _channels_HiggsProduction.end(), HiggsProductionMode) != _channels_HiggsProduction.end();

  if(!bIsAllowed)
    {
      std::cout << "E R R O R ! XSections::IsAllowedHiggsProduction(...) - The decay mode selected \"" << HiggsProductionMode << "\""
		<< " is not a supported Higgs production mode. EXIT" << std::endl;
      exit(1);
    }

    
  return; 
}


  
//****************************************************************************
Table XSections::GetAssociatedHiggsToXTable(string HiggsDecayMode,
					    bool bLatex)
//****************************************************************************
{
      
  // Initialise local variables
  int counter = 0;
  string HToX   = "-";
  string XToY_a = "-";
  string XToY_b = "-";
  string YToZ_a = "-";
  string YToZ_b = "-";
  string WToA_a = "-";
  string AToB_a = "-";
  string WToA_b = "-";
  string AToB_b = "-";
  double BR_HToX   = 1.0;
  double BR_XToY_a = 1.0;
  double BR_XToY_b = 1.0;
  double BR_YToZ_a = 1.0;
  double BR_YToZ_b = 1.0;
  double BR_WToA_a = 1.0;
  double BR_AToB_a = 1.0;
  double BR_WToA_b = 1.0;
  double BR_AToB_b = 1.0;
  double BR_Total  = 1.0;

  // Determine decay modes according to Higgs Decay mode
  std::vector<string> br_HToX (1, HiggsDecayMode);
  std::vector<string> br_XToY   = GetHiggsDaughterDecayChannels(HiggsDecayMode);
  std::vector<string> br_YToZ   = GetTauDecayChannels(); // in case we have WToTauNu.
  std::vector<string> br_XToY_a = br_XToY;
  std::vector<string> br_XToY_b = br_XToY;
  std::vector<string> br_YToZ_a = br_YToZ;
  std::vector<string> br_YToZ_b = br_YToZ;
  std::vector<string> br_WToA_a = GetWDecayChannels();
  std::vector<string> br_AToB_a;
  std::vector<string> br_WToA_b = GetWDecayChannels();
  std::vector<string> br_AToB_b;
 
  // Setup Table
  string tableType;
  if (bLatex) tableType = "LaTeX";
  else tableType = "Text";
  
  // Create Table Title
  string title = "  | ";
  title += _GetBRInLatexFormat("H" , "X", bLatex) + " | ";
  title += _GetBRInLatexFormat("X" , "Y", bLatex) + " | ";
  title += _GetBRInLatexFormat("Y" , "Z", bLatex) + " | ";
  title += _GetBRInLatexFormat("X" , "Y", bLatex) + " | ";
  title += _GetBRInLatexFormat("Y" , "Z", bLatex) + " | ";
  title += _GetBRInLatexFormat("W+", "A", bLatex) + " | ";
  title += _GetBRInLatexFormat("A" , "B", bLatex) + " | ";
  title += _GetBRInLatexFormat("W-", "A", bLatex) + " | ";
  title += _GetBRInLatexFormat("A" , "B", bLatex) + " | ";
  title += _GetInLatexFormat("BR", bLatex) + " | ";
  title += _GetInLatexFormat("XS (" + _GetInLatexFormat(_XSection_Units, bLatex) + ")", bLatex) + " | ";
  title += _GetInLatexFormat("Events @ " + auxTools.ToString(_Int_Lumi) + " (" + _IntLumi_Units + ")", bLatex);
  Table xs_BR(title, tableType);
  
  // For-loop: All combinations
  for (int i = 0; i < int(br_HToX.size()); i++){
    HToX    = br_HToX.at(i);
    BR_HToX = BRContainer[HToX];    
    

    for (int j = 0; j < int(br_XToY_a.size()); j++){
      XToY_a    = br_XToY_a.at(j);
      BR_XToY_a = BRContainer[ XToY_a ];
      if (!IsStableFinalState(XToY_a) ) br_YToZ_a = GetTauDecayChannels();
      else br_YToZ_a = GetEmptyDecayChannels();
	

      for (int k = 0; k < int(br_YToZ_a.size()); k++){
	YToZ_a    = br_YToZ_a.at(k);
	BR_YToZ_a = BRContainer[ YToZ_a ];

	
	for (int l = 0; l < int(br_XToY_b.size()); l++){
	  XToY_b    = br_XToY_b.at(l);
	  BR_XToY_b = BRContainer[ XToY_b ];
	  if (!IsStableFinalState(XToY_b) ) br_YToZ_b = GetTauDecayChannels();
	  else br_YToZ_b = GetEmptyDecayChannels();

	  
	  for (int m = 0; m < int(br_YToZ_b.size()); m++){
	    YToZ_b    = br_YToZ_b.at(m);
	    BR_YToZ_b = BRContainer[ YToZ_b ];

	    
	    for (int n = 0; n < int(br_WToA_a.size()); n++){
	      WToA_a    = br_WToA_a.at(n);
	      BR_WToA_a = BRContainer[ WToA_a ];
	      if (!IsStableFinalState(WToA_a) ) br_AToB_a = GetTauDecayChannels();
	      else br_AToB_a = GetEmptyDecayChannels();

	      
	      for (int o = 0; o < int(br_AToB_a.size()); o++){
		AToB_a    = br_AToB_a.at(o);
		BR_AToB_a = BRContainer[ AToB_a ];

		
		for (int p = 0; p < int(br_WToA_b.size()); p++){
		  WToA_b    = br_WToA_b.at(p);
		  BR_WToA_b = BRContainer[ WToA_b ];
		  if (!IsStableFinalState(WToA_b) ) br_AToB_b = GetTauDecayChannels();
		  else br_AToB_b = GetEmptyDecayChannels();

		  
		  for (int q = 0; q < int(br_AToB_b.size()); q++, counter++){
		    AToB_b    = br_AToB_b.at(q);
		    BR_AToB_b = BRContainer[ AToB_b ];

		    double BR_Total  = BR_HToX * BR_XToY_a * BR_XToY_b * BR_YToZ_a * BR_YToZ_b * BR_WToA_a * BR_WToA_b * BR_AToB_a * BR_AToB_b;
   		    double XSection  = GetHiggsXSection();
		    double IntLumi   = GetIntLumi();
   		    double NEvents   = BR_Total * XSection * IntLumi;

 		    // Table with decays in text
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( counter  , _Precision) );
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( HToX     , _Precision) );
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( XToY_a   , _Precision) );
		    xs_BR.AddRowColumn(counter, auxTools.ToString( YToZ_a   , _Precision) );
		    xs_BR.AddRowColumn(counter, auxTools.ToString( XToY_b   , _Precision) );
		    xs_BR.AddRowColumn(counter, auxTools.ToString( YToZ_b   , _Precision) );
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( WToA_a   , _Precision) );
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( AToB_a   , _Precision) );
		    xs_BR.AddRowColumn(counter, auxTools.ToString( WToA_b   , _Precision) );
		    xs_BR.AddRowColumn(counter, auxTools.ToString( AToB_b   , _Precision) );
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( BR_Total , _Precision) );
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( XSection , _Precision) );
 		    xs_BR.AddRowColumn(counter, auxTools.ToString( NEvents  , _Precision) );

 		    // Table with decays in numbers
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( counter   , _Precision) );
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_HToX   , _Precision) );
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_XToY_a , _Precision) );
		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_YToZ_a , _Precision) );
		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_XToY_b , _Precision) );
		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_YToZ_b , _Precision) );
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_WToA_a , _Precision) );
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_AToB_a , _Precision) );
		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_WToA_b , _Precision) );
		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_AToB_b , _Precision) );
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( BR_Total  , _Precision) );
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( XSection  , _Precision) );
 		    // xs_BR.AddRowColumn(counter, auxTools.ToString( NEvents   , _Precision) );

		  }
 		}
 	      }
 	    }
 	  }
 	}
      }
    }
  } 

  return xs_BR;
}



//****************************************************************************
Table XSections::_ConvertToFinalStateTable(Table myTable,
					   bool bLatex)
//****************************************************************************
{

  // Initialise local variables
  Table auxTable = myTable;
  string tableType;
  string title = "  | ";

  // Setup table style
  if (bLatex) tableType = "LaTeX";
  else tableType = "Text";

  // Setup title
  title += "Final State | ";
  title += _GetInLatexFormat("nPr", bLatex) + " | ";
  title += _GetInLatexFormat("Events @ " + auxTools.ToString(_Int_Lumi) + " (" + _IntLumi_Units + ")", bLatex);

  cout << "bLatex = " << bLatex << endl;
  
  // Re,place strings
  auxTable.ReplaceStringInTable("HToTauTau"   , _GetInLatexFormat(""     , bLatex) );
  auxTable.ReplaceStringInTable("TauToHadrons", _GetInLatexFormat("tau_h", bLatex) );
  auxTable.ReplaceStringInTable("TauToLNuNu"  , _GetInLatexFormat("l"    , bLatex) );
  auxTable.ReplaceStringInTable("WToLNu"      , _GetInLatexFormat("l"    , bLatex) );
  auxTable.ReplaceStringInTable("WToTauNu"    , _GetInLatexFormat(""     , bLatex) );
  auxTable.ReplaceStringInTable("WToHadrons"  , _GetInLatexFormat("qq"   , bLatex) );
  auxTable.ReplaceStringInTable("ZToLL"       , _GetInLatexFormat("ll"   , bLatex) );
  auxTable.ReplaceStringInTable("ZToNuNu"     , _GetInLatexFormat(""     , bLatex) );
  auxTable.ReplaceStringInTable("ZToHadrons"  , _GetInLatexFormat("qq"   , bLatex) );
  auxTable.ReplaceStringInTable("ZToTauTau"   , _GetInLatexFormat(""     , bLatex) );
  
  // auxTable.ReplaceStringInTable("HToTauTau"   , ""); //\\lTau{} \\lTau{} 
  // auxTable.ReplaceStringInTable("TauToHadrons", "\\lTauJet");
  // auxTable.ReplaceStringInTable("TauToLNuNu"  , "\\lL{}");
  // auxTable.ReplaceStringInTable("WToLNu"      , "\\lL{}");
  // auxTable.ReplaceStringInTable("WToTauNu"    , ""); //\\lTau{}
  // auxTable.ReplaceStringInTable("WToHadrons"  , "\\qQ \\qQ");
  // auxTable.ReplaceStringInTable("ZToLL"       , "\\lL{}\\lL{}");
  // auxTable.ReplaceStringInTable("ZToNuNu"     , "");
  // auxTable.ReplaceStringInTable("ZToHadrons"  , "\\qQ \\qQ");
  // auxTable.ReplaceStringInTable("ZToTauTau"   , ""); // \\lTau{} \\lTau{}
  // // auxTable.Print();
  // auxTable.Print();
  
  
  // Delete unwanted columns
  auxTable.DeleteColumn(8);
  auxTable.DeleteColumn(9);

  // Create a NEW table from the input table
  Table fsTable(title, "LaTeX");
  for (int iRow = 0; iRow < auxTable.GetNumberOfRows()-2; iRow++)
    {
      fsTable.AddRowColumn(iRow, auxTools.ToString(iRow) );
      fsTable.AddRowColumn(iRow, auxTable.GetMergedColumnsInRow(iRow, 1, 9) ); //iRow, 1, 5
      fsTable.AddRowColumn(iRow, "1" );
      fsTable.AddRowColumn(iRow, auxTable.GetColumnForRow(iRow, 12) );
    }
    

  // Replace keywords with final state particles
  vector<string> keyWords;
  keyWords.push_back(_GetInLatexFormat("l"    , bLatex) );
  keyWords.push_back(_GetInLatexFormat("tau_h", bLatex) );
  keyWords.push_back(_GetInLatexFormat("q"    , bLatex) );
  fsTable.ReplaceStringWithOccurences(0, auxTable.GetNumberOfRows()-2, 1, 1, keyWords);
  // fsTable.Print();

  // Create another NEW table from the input table
  Table fsTable_2(title, "LaTeX"); //, "c c l");
  int rowCounter = 0;
  vector<int> usedRows;
  
  // For-loop: All rows
  for (int iRow = 0; iRow < fsTable.GetNumberOfRows()-2; iRow++)
    {
      string columnRow_i = fsTable.GetColumnForRow(iRow, 1);
      double events_sum  = stod(fsTable.GetColumnForRow(iRow, 3));
      int occurrences    = 1;
      
      // Skip row if used already
      if ( find(usedRows.begin(), usedRows.end(), iRow) != usedRows.end() ) continue;
      
      // For-loop: All rows
      for (int jRow = 0; jRow < fsTable.GetNumberOfRows()-2; jRow++)
	{
	  
	  // Skip deleted or identical rows
	  if (iRow == jRow) continue;

	  // Skip row if used already
	  if ( find(usedRows.begin(), usedRows.end(), jRow) != usedRows.end() ) continue;

	  // Find duplicate column in for row with index "iRow" in row with index "jRow"
	  string columnRow_j = fsTable.GetColumnForRow(jRow, 1);
	  if (columnRow_i.compare(columnRow_j) != 0) continue;

	  // If duplicate is found add the contents of columns index 3 to the that of row "iRow". Then mark "jRow" as deleted
	  double events_j = stod(fsTable.GetColumnForRow(jRow, 3));
	  events_sum += events_j;
	  
	  // Increment occurrunces
	  occurrences++;	  

	  // Keep track of used rows
	  usedRows.push_back(jRow);
	  
	}
      // Keep track of used rows
      usedRows.push_back(iRow);      
      
      // Create new table
      fsTable_2.AddRowColumn(rowCounter, auxTools.ToString(rowCounter) );
      fsTable_2.AddRowColumn(rowCounter, columnRow_i );
      fsTable_2.AddRowColumn(rowCounter, auxTools.ToString(occurrences) );
      if(occurrences > 1) fsTable_2.AddRowColumn(rowCounter, auxTools.ToString(events_sum) );
      else fsTable_2.AddRowColumn(rowCounter, fsTable.GetColumnForRow(iRow, 3));
      rowCounter++;
    }
  
  return fsTable_2;
}



//****************************************************************************
string XSections::_GetBRInLatexFormat(string mother,
				      string decayMode,
				      bool bLatex)
//****************************************************************************
{

  string BR_Mother_To_DecayMode;
  if(bLatex) mother = _GetInLatexFormat(mother, bLatex);  
  if (decayMode.compare("-") !=0) BR_Mother_To_DecayMode = "\\BRalt{" + mother + " \\to "  + decayMode + "}";
  if (!bLatex) BR_Mother_To_DecayMode = "BR(" + mother + " -> "  + decayMode + ")";
  
  return BR_Mother_To_DecayMode;
}



//****************************************************************************
string XSections::_GetInLatexFormat(string expression,
				    bool bLatex)
//****************************************************************************
{

  if(!bLatex) return expression;
		
  if (expression.compare("fb-1")         == 0) expression = "\\sInvFb";
  else if (expression.compare("")        == 0) expression = "";
  else if (expression.compare("fb")      == 0) expression = "\\sFb";
  else if (expression.compare("H0")      == 0) expression = "\\bH{0}";
  else if (expression.compare("W+")      == 0) expression = "\\bW{+}";
  else if (expression.compare("W-")      == 0) expression = "\\bW{-}";
  else if (expression.compare("nPr")     == 0) expression = "$_{n} P ^{r}$";
  else if (expression.compare("BR")      == 0) expression = "$\\mathcal{B}$";
  else if (expression.compare"tau_h"     == 0) expression = "\\lTauJet";
  else if (expression.compare"l"         == 0) expression = "\\lL{}";
  else if (expression.compare"q"         == 0) expression = "\\qQ";
  else if (expression.compare"qq"        == 0) expression = "\\qQ \\qQ";
  else if (expression.compare"ll"        == 0) expression = "\\lL{}\\lL{}";
  else if (expression.compare"nu"        == 0) expression = "\\lNu{}";
  else if (expression.compare"nunu"      == 0) expression = "\\lNu{} \\lNu{}";  
  else if (expression.compare("XS (" + _GetInLatexFormat(_XSection_Units, bLatex) + ")") == 0) expression = "$\\sigma$ (" + _GetInLatexFormat(_XSection_Units, bLatex) + ")";
  else if (expression.compare("Events @ " + auxTools.ToString(_Int_Lumi) + " (" + _IntLumi_Units + ")") == 0) expression = "N (" + auxTools.ToString(_Int_Lumi) + " " + _GetInLatexFormat(_IntLumi_Units) + ")";
  else{}

  return expression;
}

#endif // XSections_cxx

#ifndef XSections_h
#define XSections_h

#include <iostream>

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/MCTools.C"
#include "../utilities/AuxTools.C"
#include "../utilities/HistoTools.C"

typedef map<string, vector <string> > m_Channel_To_DecayModes;
typedef map<string, vector <string> >::iterator it_Channel_To_DecayModes;

typedef map<string, string> m_Channel_To_FinalState;
typedef map<string, string>::iterator it_Channel_To_FinalState;

typedef map<string, bool> m_FinalState_To_IsStable;
typedef map<string, bool>::iterator it_FinalState_To_IsStable;

typedef map<string, double> m_Channel_To_BR;
typedef map<string, double>::iterator it_ChannelVsBR;

typedef map<double, m_Channel_To_BR> m_Mass_To_BR;
typedef map<double, m_Channel_To_BR>::iterator it_Mass_To_BR;

class XSections{

 public:
  XSections(int Precision=4, const double Higss_Mass=125.0, const double Int_Lumi=100.0);  
  ~XSections();
  vector<string> GetDecayChannels(string motherParticle);
  vector<string> GetEmptyDecayChannels(void){return _channels_Empty;}  
  vector<string> GetHDecayChannels(void){return _channels_HToAll;}
  vector<string> GetHiggsDaughterDecayChannels(string HiggsDecayMode){return m_Higgs_To_DecayModes[HiggsDecayMode];}
  vector<string> GetTauDecayChannels(void){return _channels_TauToAll;}
  vector<string> GetWDecayChannels(void){return _channels_WToAll;}
  vector<string> GetZDecayChannels(void){return _channels_ZToAll;}
  bool IsStableFinalState(string finalState);
  Table GetAssociatedHiggsToXTable(string HiggsDecayMode, bool bLatex=false);
  void IsAllowedHiggsDecay(string HiggsDecayMode);
  void IsAllowedHiggsProduction(string HiggsDecayMode);
  void PrintHiggsFinalStates(string HiggsProductionMode, string HiggsDecayMode, bool bLatex=false);
  double GetHiggsMass(void){ return _Higgs_Mass;}
  double GetHiggsWidth(void){ return _Higgs_Width;}
  double GetHiggsXSection(void){ return _Higgs_XSection;}
  double GetIntLumi(void){ return _Int_Lumi;}
  
 private:
  string _GetBRInLatexFormat(string mother, string decayMode, bool bLatex);
  string _GetInLatexFormat(string expression, bool bLatex=true);
  Table _ConvertToFinalStateTable(Table myTable, bool bLatex=false);
  void _CheckInputParams(void);
  void _ConvertChannelToFinalState(void);
  void _InitVars(void);
  void _SetBranchingRatioMaps(void);
  void _SetDecayChannels(void);  
  void _SetHiggsDaughterDecayChannelsMap(void);
  void _SetHiggsWidthMaps(void);
  void _SetProductionChannels(void);
  void _SetStableFinalState(void);
  
  // Objects & Variables
  AuxTools auxTools;
  m_Channel_To_BR BRContainer;
  string _IntLumi_Units;
  string _XSection_Units;
  string _Mass_Units;
  string _Width_Units;
  double _Int_Lumi;
  double _Higgs_XSection;
  double _Higgs_Mass;
  double _Higgs_Width;
  int _Precision;
  vector<string> _channels_WToAll;
  vector<string> _channels_HToAll;
  vector<string> _channels_ZToAll;
  vector<string> _channels_TauToAll;
  vector<string> _channels_HiggsProduction;
  vector<string> _channels_Empty;
  
  m_Mass_To_BR m_HiggsMass_To_BR;
  m_Channel_To_FinalState m_Channel_To_FS;
  m_FinalState_To_IsStable m_FS_To_IsStable;
  m_Channel_To_DecayModes m_Higgs_To_DecayModes;

  it_Mass_To_BR it_HiggsMass_To_BR;
  it_Channel_To_FinalState it_Channel_To_FS;
  it_FinalState_To_IsStable it_FS_To_IsStable;
  it_Channel_To_DecayModes it_Higgs_To_DecayModes;

};

#endif // XSections_h


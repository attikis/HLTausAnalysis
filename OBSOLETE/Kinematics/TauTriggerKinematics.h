#ifndef TauTriggerKinematics_h
#define TauTriggerKinematics_h

#include <iostream>
#include <stdlib.h>     // exit, EXIT_FAILURE

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/MCTools.C"
#include "../utilities/AuxTools.C"
#include "../utilities/HistoTools.C"
#include "../utilities/L1Tracks.h"

class TauTriggerKinematics : public TreeAnalyserMC{
 public:
 TauTriggerKinematics(const std::string SamplePath,
		      const std::string SampleName,
		      const std::string text_, 
		      const int maxEvents_ = -1, 
		      TTree* tree=0) : 
  TreeAnalyserMC("TauTriggerKinematics", SamplePath, SampleName, text_, maxEvents_, tree) { 
    mcSample     = SampleName; 
    tauDecayMode = text_;
  };
  virtual void Loop();
  
  // Variable declaration
  int nMaxNumOfHTausPossible;
  int realTauMom;
  std::string mcSample;
  std::string tauDecayMode;
  double tauVisEtCut;
  
 private:
  void InitVars(void);

  void BookHistos(void);

  std::vector<int> GetTriggerHadronicTaus(const int tauMom,
					  const double hTauVisEtCut);
  
  TLorentzVector GetHadronicTauVisP4(const int genP_index);
  
  // Objects/Variables
  AuxTools tools; 
  HistoTools histos; 

  // Tau
  TH1D* hHadronicTau_N;
  TH1D* hHadronicTau_Pt;
  TH1D* hHadronicTau_Eta;
  TH1D* hHadronicTau_Phi;
  TH1D* hHadronicTau_Mass;
  TH1D* hHadronicTau_Charge;
  TH1D* hHadronicTau_PdgId;
  TH1D* hHadronicTau_Status;
  TH1D* hHadronicTau_VertexX;
  TH1D* hHadronicTau_VertexY;
  TH1D* hHadronicTau_VertexZ;

  // Visible-Tau
  TH1D* hHadronicTau_VisEt;
  TH1D* hHadronicTau_VisEta;
  TH1D* hHadronicTau_VisPhi;
  TH1D* hHadronicTau_VisMass;
  TH1D* hHadronicTau_DecayMode;
  TH1D* hHadronicTau_EcalFraction;
  TH1D* hHadronicTau_HcalFraction;

  // Tau Charged Pions 
  TH1D* hHadronicTau_ChargedPion_N;
  TH1D* hHadronicTau_ChargedPion_Pt;
  TH1D* hHadronicTau_ChargedPion_Eta;
  TH1D* hHadronicTau_ChargedPion_Phi;
  TH1D* hHadronicTau_ChargedPion_Mass;
  TH1D* hHadronicTau_ChargedPion_Charge;
  TH1D* hHadronicTau_ChargedPion_PdgId;
  TH1D* hHadronicTau_ChargedPion_Status;
  TH1D* hHadronicTau_ChargedPion_VertexX;
  TH1D* hHadronicTau_ChargedPion_VertexY;
  TH1D* hHadronicTau_ChargedPion_VertexZ;

  // Tau Neutral Pions
  TH1D* hHadronicTau_NeutralPion_N;
  TH1D* hHadronicTau_NeutralPion_Pt; 
  TH1D* hHadronicTau_NeutralPion_Eta; 
  TH1D* hHadronicTau_NeutralPion_Phi;
  TH1D* hHadronicTau_NeutralPion_Mass;
  TH1D* hHadronicTau_NeutralPion_Charge;
  TH1D* hHadronicTau_NeutralPion_PdgId;
  TH1D* hHadronicTau_NeutralPion_Status;
  TH1D* hHadronicTau_NeutralPion_VertexX;
  TH1D* hHadronicTau_NeutralPion_VertexY;
  TH1D* hHadronicTau_NeutralPion_VertexZ;

  // Leading Charged Pion
  TH1D*  hHadronicTau_LdgChPion_Pt;
  TH1D*  hHadronicTau_LdgChPion_Eta;
  TH1D*  hHadronicTau_LdgChPion_Phi;
  TH1D*  hHadronicTau_LdgChPion_Mass;
  TH1D*  hHadronicTau_LdgChPion_Charge;
  TH1D*  hHadronicTau_LdgChPion_PdgId;
  TH1D*  hHadronicTau_LdgChPion_Status;
  TH1D*  hHadronicTau_LdgChPion_VertexX;
  TH1D*  hHadronicTau_LdgChPion_VertexY;
  TH1D*  hHadronicTau_LdgChPion_VertexZ;
  TH1D*  hHadronicTau_LdgChPion_Rtau;
  TH2D*  hHadronicTau_LdgChPion_Pt_VisEt;

  // 3-prong
  TH1D*  hHadronicTau_LdgChPion_DeltaRMax;
  TH2D*  hHadronicTau_LdgChPion_DeltaRMax_Pt;
  TH2D*  hHadronicTau_LdgChPion_DeltaRMax_TauEt;
  TH2D*  hHadronicTau_LdgChPion_DeltaRMax_VisTauEt;

  TH1D*  hHadronicTau_LdgChPion_DeltaRMax_MinPt;
  TH2D*  hHadronicTau_LdgChPion_DeltaRMax_Pt_MinPt;
  TH2D*  hHadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt;
  TH2D*  hHadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt;
  
  // DiTau System
  TH1D* hDiTau_Pt;
  TH1D* hDiTau_InvMass;
  TH1D* hDiTau_VisMass;

};

#endif // TauTriggerKinematics_h

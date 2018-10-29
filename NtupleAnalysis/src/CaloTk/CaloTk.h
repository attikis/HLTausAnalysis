#ifndef CaloTk_h
#define CaloTk_h

// System 
#include <iostream>
#include <stdlib.h> 

// User
#include "../Framework/src/TreeAnalyserMC.C"
#include "../Framework/src/TreeReaderMC.C"


#include "../Auxiliary/src/AuxTools.C"
#include "../Auxiliary/src/Table.C"
#include "../Auxiliary/src/MCTools.C"
#include "../Auxiliary/src/HistoTools.C"
#include "../Auxiliary/src/Datasets.C" 

#include "../DataFormat/src/L1TkTauParticle.C"
#include "../DataFormat/src/GenParticle.C"
#include "../DataFormat/src/TrackingParticle.C"
#include "../DataFormat/interface/TTTrack.h"
//#include "../DataFormat/interface/TTPixelTrack.h"
//#include "../DataFormat/src/L1EG.C"
#include "../DataFormat/src/L1Jet.C"
#include "../DataFormat/src/L1Tau.C"
#include "../DataFormat/src/L1Sum.C"
#include "../DataFormat/src/L1CaloTP.C"

// #include "../Plugins/src/L1TkPrimaryVertex.C"
//#include "../Plugins/src/L1PixelTrackFit.C"

// ROOT
#include "TEfficiency.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

//marina
#include <vector>
#include <iostream>

using namespace std;

class CaloTk : public TreeAnalyserMC{
 public:
  // Constructors/Destructors
  ~CaloTk(){};
 CaloTk(const string SamplePath,
	const string SampleName,
	const string text_, 
	const int maxEvents_ = -1, 
	TChain* chain=0) :

  TreeAnalyserMC("", SamplePath, SampleName, text_, maxEvents_, chain)
    { 
      auxTools_.StopwatchStart();
      mcSample = SampleName;
      InitObjects();
    };
  
  // Public Variables
  virtual void Loop();

  void PrintSettings(void);

  void ApplyDiTauZMatching(vector<L1TkTauParticle> &L1Taus);

  void GetShrinkingConeSizes(double tk_pt,
			     double sigCone_Constant,
			     double isoCone_Constant,
			     const double sigCone_dRCutoff,
			     double &sigCone_dRMin,
			     double &sigCone_dRMax,
			     double &isoCone_dRMin,
			     double &isoCone_dRMax);

  double GetDonutRatio(L1TkTauParticle &L1TkTau, 
		       vector<TTTrack> isoTTTracks,
		       bool bUseCone);


  double GetJetWidth(vector<TTTrack> sigTks, vector<TTTrack> isoTks,
		     TLorentzVector sigTks_p4, TLorentzVector isoTks_p4);

  void GetMatchingTrack(L1TkTauParticle &L1TkTau,
			L1Tau L1CaloTau,
			vector<TTTrack> TTTracks);
  
  void GetSigConeTracks(L1TkTauParticle &L1TkTau,
			vector<TTTrack> TTTracks,
			double sigConeTks_dPOCAz,
			double sigConeTks_invMass);
  
  void GetIsolationTracks(L1TkTauParticle &L1TkTau,
			vector<TTTrack> isoTTTracks,
			double isoConeTks_dPOCAz);

  void GetIsolationValues(L1TkTauParticle &L1TkTau, bool bUseCone);
  
  void GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
			      vector<GenParticle> hadGenTaus);			    

  void GetLdgAndSubldgIndices(vector<L1TkTauParticle> myTaus,
			      int &iLdg,
			      int &iSubldg);
  
  // Public Variables
  bool DEBUG;
  bool mcMatching_unique;
  double diTau_deltaPOCAz;

  // L1TkTau - Matching track
  string seedTk_Collection;
  int seedTk_nFitParams;
  double seedTk_minPt;
  double seedTk_minEta;
  double seedTk_maxEta;
  double seedTk_maxChiSq;
  double seedTk_minStubs;

  // Signal Cone Tracks
  string sigConeTks_Collection;
  int sigConeTks_nFitParams;
  double sigConeTks_minPt;
  double sigConeTks_minEta;
  double sigConeTks_maxEta;
  double sigConeTks_maxChiSq;
  unsigned int sigConeTks_minStubs;
  double sigConeTks_dPOCAz;
  double sigConeTks_maxInvMass;

  // Isolation Cone Tracks
  string isoConeTks_Collection;
  int isoConeTks_nFitParams;
  double isoConeTks_minPt;
  double isoConeTks_minEta;
  double isoConeTks_maxEta;
  double isoConeTks_maxChiSq;
  unsigned int isoConeTks_minStubs;
  double isoConeTks_dPOCAz;

  double mcMatching_dRMax;
  double pv_deltaZMax;
  double pv_z;
  
  double _eta_C;
  double _eta_F;

  //
  double sigCone_Constant;
  double sigCone_cutoffDeltaR;
  double sigCone_dRMax;
  double sigCone_dRMin;
  double isoCone_Constant;
  double isoCone_dRMax;
  double isoCone_dRMin;
  bool   isoCone_useCone; //instead of annulus

  double tau_jetWidth;
  double tau_vtxIsoWP;
  double tau_relIsoWP;
  double tau_relIsodZ0;
  //
  int nMaxNumOfHTausPossible;
  int realTauMom;
  string caloCorrectionType;
  string mcSample;
  string pv_producer;
  string tauDecayMode;
  struct SortAscendingAbs{ bool operator() (double a, double b) const { return abs(a) > abs(b); } }; 
  struct SortDescendingAbs{ bool operator() (double a, double b) const { return abs(a) < abs(b); } }; 
  
 private:
  // Function declaration
  void BookHistos_(void);
  void WriteHistos_(void);
  void InitObjects(void);
  void InitVars_(void);

  void FillEfficiency_(TH1D *hSignalEfficiency,
		       const double ldgEt);

  void FillEfficiency_(TH2D *hSignalEfficiency,
		       const double ldgEt1,
		       const double ldgEt2);  

  void FillRate_(TH1D *hRate,
		 const double ldgEt);

  void FillRate_(TH2D *hRate,
		 const double ldgEt1,
		 const double ldgEt2);  

  void FinaliseEffHisto_(TH1D *histo, 
			 const int nEvtsTotal);

  void FinaliseEffHisto_(TH2D *histo, 
			 const int nEvtsTotal);  

  void FillTurnOn_Numerator_(vector<L1TkTauParticle> L1Taus,
			     const double minEt,
			     TH1D *hTurnOn, 
			     TH1D *hTurnOn_1pr, 
			     TH1D *hTurnOn_3pr, 
			     TH1D *hTurnOn_withNeutrals, 
			     TH1D *hTurnOn_noNeutrals);

  void FillSingleTau_(vector<L1TkTauParticle> L1Taus,
		      TH1D *hRate,
		      TH1D *hEfficiency,
		      double minEta=0.0,
		      double maxEta=999.9);

  void FillDiTau_(vector<L1TkTauParticle> L1Taus, 
		  TH1D *hRate,
		  TH1D *hEfficiency,
		  double minEta=0.0,
		  double maxEta=999.9);

  void FillDiTau_(vector<L1TkTauParticle> L1Taus1,
		  vector<L1TkTauParticle> L1Taus2,
		  TH2D *hRate,
		  TH2D *hEfficiency);

  vector<L1TkTauParticle> GetMcMatchedL1Taus(vector<L1TkTauParticle> L1Taus);

  bool IsWithinEtaRegion(string etaRegion,
			 double eta);

  
  // Variable declaration
  // L1TkPrimaryVertex *pvProducer;
  AuxTools auxTools_;
  MCTools mcTools_;
  Datasets datasets_;
  HistoTools histoTools_;
  bool bFoundAllTaus_;

  // GenParticles Histograms
  TH2D* hGenP_VisEt_Vs_dRMaxLdgPion;
  TH2D* hGenP_PtLdg_Vs_dRMaxLdgPion;

  // Counters
  TH1D* hCounters;

  // L1CaloTaus
  // http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_3_7/doc/html/d9/d40/L1Trigger_2interface_2Tau_8h_source.html#l00019
  // http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_3_7/doc/html/d1/dc8/classl1t_1_1CaloTools.html
  TH1D* hL1CaloTau_Et; 
  TH1D* hL1CaloTau_Eta;
  TH1D* hL1CaloTau_Phi;
  TH1D* hL1CaloTau_IEt;
  TH1D* hL1CaloTau_IEta; // ieta of seed tower
  TH1D* hL1CaloTau_IPhi; // iphi of seed tower
  TH1D* hL1CaloTau_Iso;
  // TH1D* hL1CaloTau_Bx;
  TH1D* hL1CaloTau_TowerIEta;
  TH1D* hL1CaloTau_TowerIPhi;
  TH1D* hL1CaloTau_RawEt; // raw (uncalibrated) cluster sum
  TH1D* hL1CaloTau_IsoEt; // raw isolation sum - cluster sum
  TH1D* hL1CaloTau_NTT;   // n towers above threshold
  TH1D* hL1CaloTau_HasEM;
  TH1D* hL1CaloTau_IsMerged;
  // TH1D* hL1CaloTau_HwQual; //integer hardware (hw) value

  // L1Taus
  TH1D* hL1Tau_SeedTk_DeltaR;
  TH1D* hL1Tau_SeedTk_PtRel;
  TH1D* hL1Tau_SeedTk_Pt;
  TH1D* hL1Tau_SeedTk_Eta;
  TH1D* hL1Tau_SeedTk_POCAz;
  TH1D* hL1Tau_SeedTk_NStubs;
  TH1D* hL1Tau_SeedTk_NPsStubs;
  TH1D* hL1Tau_SeedTk_NBarrelStubs;
  TH1D* hL1Tau_SeedTk_NEndcapStubs;
  TH1D* hL1Tau_SeedTk_ChiSquared;
  TH1D* hL1Tau_SeedTk_RedChiSquared;
  TH1D* hL1Tau_SeedTk_IsGenuine;
  TH1D* hL1Tau_SeedTk_IsUnknown;
  TH1D* hL1Tau_SeedTk_IsCombinatoric;

  TH1D* hL1Tau_SigTks_Pt;
  TH1D* hL1Tau_SigTks_PtRel;
  TH1D* hL1Tau_SigTks_Eta;
  TH1D* hL1Tau_SigTks_POCAz;
  TH1D* hL1Tau_SigTks_DeltaPOCAz;
  TH1D* hL1Tau_SigTks_DeltaR;
  TH1D* hL1Tau_SigTks_NStubs;
  TH1D* hL1Tau_SigTks_NPsStubs;
  TH1D* hL1Tau_SigTks_NBarrelStubs;
  TH1D* hL1Tau_SigTks_NEndcapStubs;
  TH1D* hL1Tau_SigTks_ChiSquared;
  TH1D* hL1Tau_SigTks_RedChiSquared;

  TH1D* hL1Tau_IsoTks_Pt;
  TH1D* hL1Tau_IsoTks_PtRel;
  TH1D* hL1Tau_IsoTks_Eta;
  TH1D* hL1Tau_IsoTks_POCAz;
  TH1D* hL1Tau_IsoTks_DeltaPOCAz;
  TH1D* hL1Tau_IsoTks_DeltaR;
  TH1D* hL1Tau_IsoTks_NStubs;
  TH1D* hL1Tau_IsoTks_NPsStubs;
  TH1D* hL1Tau_IsoTks_NBarrelStubs;
  TH1D* hL1Tau_IsoTks_NEndcapStubs;
  TH1D* hL1Tau_IsoTks_ChiSquared;
  TH1D* hL1Tau_IsoTks_RedChiSquared;

  TH1D* hL1Tau_Multiplicity;
  TH1D* hL1Tau_Multiplicity_MC;
  TH1D* hL1Tau_JetWidth;
  TH1D* hL1Tau_DonutRatio;
  TH1D* hL1Tau_NSigTks;
  TH1D* hL1Tau_SigTksEt;
  TH1D* hL1Tau_SigTksEta;
  TH1D* hL1Tau_NIsoTks;
  TH1D* hL1Tau_IsoTksEt;
  TH1D* hL1Tau_IsoTksEta;
  TH1D* hL1Tau_InvMass;
  TH1D* hL1Tau_IsoConeMass;
  TH1D* hL1Tau_IsoAnnulusMass;
  TH1D* hL1Tau_SigConeRMin;
  TH1D* hL1Tau_SigConeRMax;
  TH1D* hL1Tau_IsoConeRMin;
  TH1D* hL1Tau_IsoConeRMax;
  TH1D* hL1Tau_Charge;
  TH1D* hL1Tau_RelIso;
  TH1D* hL1Tau_VtxIso;
  TH2D* hL1Tau_VtxIso_Vs_RelIso;
  TH1D* hL1Tau_DeltaRGenP;

  // L1IsoTaus
  TH1D* hL1IsoTau_SeedTk_DeltaR;
  TH1D* hL1IsoTau_SeedTk_PtRel;
  TH1D* hL1IsoTau_SeedTk_Pt;
  TH1D* hL1IsoTau_SeedTk_Eta;
  TH1D* hL1IsoTau_SeedTk_POCAz;
  TH1D* hL1IsoTau_SeedTk_NStubs;
  TH1D* hL1IsoTau_SeedTk_NPsStubs;
  TH1D* hL1IsoTau_SeedTk_NBarrelStubs;
  TH1D* hL1IsoTau_SeedTk_NEndcapStubs;
  TH1D* hL1IsoTau_SeedTk_ChiSquared;
  TH1D* hL1IsoTau_SeedTk_RedChiSquared;
  TH1D* hL1IsoTau_SeedTk_IsGenuine;
  TH1D* hL1IsoTau_SeedTk_IsUnknown;
  TH1D* hL1IsoTau_SeedTk_IsCombinatoric;
  TH1D* hL1IsoTau_SigTks_Pt;
  TH1D* hL1IsoTau_SigTks_PtRel;
  TH1D* hL1IsoTau_SigTks_Eta;
  TH1D* hL1IsoTau_SigTks_POCAz;
  TH1D* hL1IsoTau_SigTks_DeltaPOCAz;
  TH1D* hL1IsoTau_SigTks_DeltaR;
  TH1D* hL1IsoTau_SigTks_NStubs;
  TH1D* hL1IsoTau_SigTks_NPsStubs;
  TH1D* hL1IsoTau_SigTks_NBarrelStubs;
  TH1D* hL1IsoTau_SigTks_NEndcapStubs;
  TH1D* hL1IsoTau_SigTks_ChiSquared;
  TH1D* hL1IsoTau_SigTks_RedChiSquared;
  TH1D* hL1IsoTau_IsoTks_Pt;
  TH1D* hL1IsoTau_IsoTks_PtRel;
  TH1D* hL1IsoTau_IsoTks_Eta;
  TH1D* hL1IsoTau_IsoTks_POCAz;
  TH1D* hL1IsoTau_IsoTks_DeltaPOCAz;
  TH1D* hL1IsoTau_IsoTks_DeltaR;
  TH1D* hL1IsoTau_IsoTks_NStubs;
  TH1D* hL1IsoTau_IsoTks_NPsStubs;
  TH1D* hL1IsoTau_IsoTks_NBarrelStubs;
  TH1D* hL1IsoTau_IsoTks_NEndcapStubs;
  TH1D* hL1IsoTau_IsoTks_ChiSquared;
  TH1D* hL1IsoTau_IsoTks_RedChiSquared;

  TH1D* hL1IsoTau_Multiplicity;
  TH1D* hL1IsoTau_Multiplicity_MC;
  TH1D* hL1IsoTau_JetWidth;
  TH1D* hL1IsoTau_DonutRatio;
  TH1D* hL1IsoTau_NSigTks;
  TH1D* hL1IsoTau_SigTksEt;
  TH1D* hL1IsoTau_SigTksEta;
  TH1D* hL1IsoTau_NIsoTks;
  TH1D* hL1IsoTau_IsoTksEt;
  TH1D* hL1IsoTau_IsoTksEta;
  TH1D* hL1IsoTau_InvMass;
  TH1D* hL1IsoTau_IsoConeMass;
  TH1D* hL1IsoTau_IsoAnnulusMass;
  TH1D* hL1IsoTau_SigConeRMin;
  TH1D* hL1IsoTau_SigConeRMax;
  TH1D* hL1IsoTau_IsoConeRMin;
  TH1D* hL1IsoTau_IsoConeRMax;
  TH1D* hL1IsoTau_Charge;
  TH1D* hL1IsoTau_RelIso;
  TH1D* hL1IsoTau_VtxIso;
  TH2D* hL1IsoTau_VtxIso_Vs_RelIso;
  TH1D* hL1IsoTau_DeltaRGenP;

  // Resolutions
  TH1D* hL1IsoTau_ResolutionEt;
  TH1D* hL1IsoTau_ResolutionEt_1pr;
  TH1D* hL1IsoTau_ResolutionEt_3pr;
  TH1D* hL1IsoTau_ResolutionEt_withNeutrals;
  TH1D* hL1IsoTau_ResolutionEt_noNeutrals;

  TH1D* hL1IsoTau_ResolutionEta;
  TH1D* hL1IsoTau_ResolutionEta_1pr;
  TH1D* hL1IsoTau_ResolutionEta_3pr;
  TH1D* hL1IsoTau_ResolutionEta_withNeutrals;
  TH1D* hL1IsoTau_ResolutionEta_noNeutrals;

  TH1D* hL1IsoTau_ResolutionPhi;
  TH1D* hL1IsoTau_ResolutionPhi_1pr;
  TH1D* hL1IsoTau_ResolutionPhi_3pr;
  TH1D* hL1IsoTau_ResolutionPhi_withNeutrals;
  TH1D* hL1IsoTau_ResolutionPhi_noNeutrals;

  TH1D* hL1IsoTau_ResolutionEt_C;
  TH1D* hL1IsoTau_ResolutionEta_C;
  TH1D* hL1IsoTau_ResolutionPhi_C;
  TH1D* hL1IsoTau_ResolutionEt_I;
  TH1D* hL1IsoTau_ResolutionEta_I;
  TH1D* hL1IsoTau_ResolutionPhi_I;
  TH1D* hL1IsoTau_ResolutionEt_F;
  TH1D* hL1IsoTau_ResolutionEta_F;
  TH1D* hL1IsoTau_ResolutionPhi_F;
  
  // SingleTau: Rates
  TH1D* hCalo_Rate; // Inclusive = C+I+F
  TH1D* hCalo_Rate_C;
  TH1D* hCalo_Rate_I;
  TH1D* hCalo_Rate_F;
  TH1D* hTk_Rate;
  TH1D* hTk_Rate_C;
  TH1D* hTk_Rate_I;
  TH1D* hTk_Rate_F;
  TH1D* hVtxIso_Rate;
  TH1D* hVtxIso_Rate_C;
  TH1D* hVtxIso_Rate_I;
  TH1D* hVtxIso_Rate_F;
  TH1D* hRelIso_Rate;
  TH1D* hRelIso_Rate_C;
  TH1D* hRelIso_Rate_I;
  TH1D* hRelIso_Rate_F;
  TH1D* hVtxIsoLoose_Rate;
  TH1D* hVtxIsoLoose_Rate_C;
  TH1D* hVtxIsoLoose_Rate_I;
  TH1D* hVtxIsoLoose_Rate_F;
  TH1D* hVtxIsoTight_Rate;
  TH1D* hVtxIsoTight_Rate_C;
  TH1D* hVtxIsoTight_Rate_I;
  TH1D* hVtxIsoTight_Rate_F;
  TH1D* hRelIsoLoose_Rate;
  TH1D* hRelIsoLoose_Rate_C;
  TH1D* hRelIsoLoose_Rate_I;
  TH1D* hRelIsoLoose_Rate_F;
  TH1D* hRelIsoTight_Rate;
  TH1D* hRelIsoTight_Rate_C;
  TH1D* hRelIsoTight_Rate_I;
  TH1D* hRelIsoTight_Rate_F;

  // SingleTau: Efficiencies
  TH1D* hTk_Eff;
  TH1D* hTk_Eff_C;
  TH1D* hTk_Eff_I;
  TH1D* hTk_Eff_F;
  TH1D* hVtxIso_Eff;
  TH1D* hVtxIso_Eff_C;
  TH1D* hVtxIso_Eff_I;
  TH1D* hVtxIso_Eff_F;      
  TH1D* hRelIso_Eff;
  TH1D* hRelIso_Eff_C;
  TH1D* hRelIso_Eff_I;
  TH1D* hRelIso_Eff_F;      
  TH1D* hVtxIsoLoose_Eff;
  TH1D* hVtxIsoLoose_Eff_C;
  TH1D* hVtxIsoLoose_Eff_I;
  TH1D* hVtxIsoLoose_Eff_F;      
  TH1D* hVtxIsoTight_Eff;
  TH1D* hVtxIsoTight_Eff_C;
  TH1D* hVtxIsoTight_Eff_I;
  TH1D* hVtxIsoTight_Eff_F;      
  TH1D* hRelIsoLoose_Eff;
  TH1D* hRelIsoLoose_Eff_C;
  TH1D* hRelIsoLoose_Eff_I;
  TH1D* hRelIsoLoose_Eff_F;      
  TH1D* hRelIsoTight_Eff;
  TH1D* hRelIsoTight_Eff_C;
  TH1D* hRelIsoTight_Eff_I;
  TH1D* hRelIsoTight_Eff_F;      

  // DiTau: Rates
  TH1D* hDiTau_Rate_Tk;
  TH1D* hDiTau_Rate_Tk_C;
  TH1D* hDiTau_Rate_Tk_I;
  TH1D* hDiTau_Rate_Tk_F;
  TH1D* hDiTau_Rate_VtxIso;
  TH1D* hDiTau_Rate_VtxIso_C;
  TH1D* hDiTau_Rate_VtxIso_I;
  TH1D* hDiTau_Rate_VtxIso_F;
  TH1D* hDiTau_Rate_RelIso;
  TH1D* hDiTau_Rate_RelIso_C;
  TH1D* hDiTau_Rate_RelIso_I;
  TH1D* hDiTau_Rate_RelIso_F;
  TH1D* hDiTau_Rate_VtxIsoLoose;
  TH1D* hDiTau_Rate_VtxIsoLoose_C;
  TH1D* hDiTau_Rate_VtxIsoLoose_I;
  TH1D* hDiTau_Rate_VtxIsoLoose_F;
  TH1D* hDiTau_Rate_VtxIsoTight;
  TH1D* hDiTau_Rate_VtxIsoTight_C;
  TH1D* hDiTau_Rate_VtxIsoTight_I;
  TH1D* hDiTau_Rate_VtxIsoTight_F;
  TH1D* hDiTau_Rate_RelIsoLoose;
  TH1D* hDiTau_Rate_RelIsoLoose_C;
  TH1D* hDiTau_Rate_RelIsoLoose_I;
  TH1D* hDiTau_Rate_RelIsoLoose_F;
  TH1D* hDiTau_Rate_RelIsoTight;
  TH1D* hDiTau_Rate_RelIsoTight_C;
  TH1D* hDiTau_Rate_RelIsoTight_I;
  TH1D* hDiTau_Rate_RelIsoTight_F;

  // DiTau: Efficiencies
  TH1D* hDiTau_Eff_Tk;
  TH1D* hDiTau_Eff_Tk_C;
  TH1D* hDiTau_Eff_Tk_I;
  TH1D* hDiTau_Eff_Tk_F;
  TH1D* hDiTau_Eff_VtxIso;
  TH1D* hDiTau_Eff_VtxIso_C;
  TH1D* hDiTau_Eff_VtxIso_I;
  TH1D* hDiTau_Eff_VtxIso_F;
  TH1D* hDiTau_Eff_RelIso;
  TH1D* hDiTau_Eff_RelIso_C;
  TH1D* hDiTau_Eff_RelIso_I;
  TH1D* hDiTau_Eff_RelIso_F;
  TH1D* hDiTau_Eff_VtxIsoLoose;
  TH1D* hDiTau_Eff_VtxIsoLoose_C;
  TH1D* hDiTau_Eff_VtxIsoLoose_I;
  TH1D* hDiTau_Eff_VtxIsoLoose_F;
  TH1D* hDiTau_Eff_VtxIsoTight;
  TH1D* hDiTau_Eff_VtxIsoTight_C;
  TH1D* hDiTau_Eff_VtxIsoTight_I;
  TH1D* hDiTau_Eff_VtxIsoTight_F;
  TH1D* hDiTau_Eff_RelIsoLoose;
  TH1D* hDiTau_Eff_RelIsoLoose_C;
  TH1D* hDiTau_Eff_RelIsoLoose_I;
  TH1D* hDiTau_Eff_RelIsoLoose_F;
  TH1D* hDiTau_Eff_RelIsoTight;
  TH1D* hDiTau_Eff_RelIsoTight_C;
  TH1D* hDiTau_Eff_RelIsoTight_I;
  TH1D* hDiTau_Eff_RelIsoTight_F;

  // DiTau (Tk-Other)
  TH2D* hDiTau_Rate_Tk_VtxIso;
  TH2D* hDiTau_Rate_Tk_RelIso;
  TH2D* hDiTau_Rate_Tk_VtxIsoLoose;
  TH2D* hDiTau_Rate_Tk_VtxIsoTight;
  TH2D* hDiTau_Rate_Tk_RelIsoLoose;
  TH2D* hDiTau_Rate_Tk_RelIsoTight;

  TH2D* hDiTau_Eff_Tk_VtxIso;
  TH2D* hDiTau_Eff_Tk_RelIso;
  TH2D* hDiTau_Eff_Tk_VtxIsoLoose;
  TH2D* hDiTau_Eff_Tk_VtxIsoTight;
  TH2D* hDiTau_Eff_Tk_RelIsoLoose;
  TH2D* hDiTau_Eff_Tk_RelIsoTight;

  // Turn-Ons
  // TEfficiency* pEff; //fixme: convert all turn-ons
  TH1D* hMcHadronicTau_VisEt;
  TH1D* hMcHadronicTau_VisEt_1pr;
  TH1D* hMcHadronicTau_VisEt_3pr;
  TH1D* hMcHadronicTau_VisEt_withNeutrals;
  TH1D* hMcHadronicTau_VisEt_noNeutrals;

  TH1D* hTk_TurnOn25;
  TH1D* hTk_TurnOn25_1pr;
  TH1D* hTk_TurnOn25_3pr;
  TH1D* hTk_TurnOn25_withNeutrals;
  TH1D* hTk_TurnOn25_noNeutrals;

  TH1D* hVtxIso_TurnOn25;
  TH1D* hVtxIso_TurnOn25_1pr;
  TH1D* hVtxIso_TurnOn25_3pr;
  TH1D* hVtxIso_TurnOn25_withNeutrals;
  TH1D* hVtxIso_TurnOn25_noNeutrals;

  TH1D* hRelIso_TurnOn25;
  TH1D* hRelIso_TurnOn25_1pr;
  TH1D* hRelIso_TurnOn25_3pr;
  TH1D* hRelIso_TurnOn25_withNeutrals;
  TH1D* hRelIso_TurnOn25_noNeutrals;

  TH1D* hVtxIsoLoose_TurnOn25;
  TH1D* hVtxIsoLoose_TurnOn25_1pr;
  TH1D* hVtxIsoLoose_TurnOn25_3pr;
  TH1D* hVtxIsoLoose_TurnOn25_withNeutrals;
  TH1D* hVtxIsoLoose_TurnOn25_noNeutrals;

  TH1D* hVtxIsoTight_TurnOn25;
  TH1D* hVtxIsoTight_TurnOn25_1pr;
  TH1D* hVtxIsoTight_TurnOn25_3pr;
  TH1D* hVtxIsoTight_TurnOn25_withNeutrals;
  TH1D* hVtxIsoTight_TurnOn25_noNeutrals;

  TH1D* hRelIsoLoose_TurnOn25;
  TH1D* hRelIsoLoose_TurnOn25_1pr;
  TH1D* hRelIsoLoose_TurnOn25_3pr;
  TH1D* hRelIsoLoose_TurnOn25_withNeutrals;
  TH1D* hRelIsoLoose_TurnOn25_noNeutrals;

  TH1D* hRelIsoTight_TurnOn25;
  TH1D* hRelIsoTight_TurnOn25_1pr;
  TH1D* hRelIsoTight_TurnOn25_3pr;
  TH1D* hRelIsoTight_TurnOn25_withNeutrals;
  TH1D* hRelIsoTight_TurnOn25_noNeutrals;

  TH1D* hTk_TurnOn50;
  TH1D* hTk_TurnOn50_1pr;
  TH1D* hTk_TurnOn50_3pr;
  TH1D* hTk_TurnOn50_withNeutrals;
  TH1D* hTk_TurnOn50_noNeutrals;

  TH1D* hVtxIso_TurnOn50;
  TH1D* hVtxIso_TurnOn50_1pr;
  TH1D* hVtxIso_TurnOn50_3pr;
  TH1D* hVtxIso_TurnOn50_withNeutrals;
  TH1D* hVtxIso_TurnOn50_noNeutrals;

  TH1D* hRelIso_TurnOn50;
  TH1D* hRelIso_TurnOn50_1pr;
  TH1D* hRelIso_TurnOn50_3pr;
  TH1D* hRelIso_TurnOn50_withNeutrals;
  TH1D* hRelIso_TurnOn50_noNeutrals;

  TH1D* hVtxIsoLoose_TurnOn50;
  TH1D* hVtxIsoLoose_TurnOn50_1pr;
  TH1D* hVtxIsoLoose_TurnOn50_3pr;
  TH1D* hVtxIsoLoose_TurnOn50_withNeutrals;
  TH1D* hVtxIsoLoose_TurnOn50_noNeutrals;

  TH1D* hVtxIsoTight_TurnOn50;
  TH1D* hVtxIsoTight_TurnOn50_1pr;
  TH1D* hVtxIsoTight_TurnOn50_3pr;
  TH1D* hVtxIsoTight_TurnOn50_withNeutrals;
  TH1D* hVtxIsoTight_TurnOn50_noNeutrals;

  TH1D* hRelIsoLoose_TurnOn50;
  TH1D* hRelIsoLoose_TurnOn50_1pr;
  TH1D* hRelIsoLoose_TurnOn50_3pr;
  TH1D* hRelIsoLoose_TurnOn50_withNeutrals;
  TH1D* hRelIsoLoose_TurnOn50_noNeutrals;

  TH1D* hRelIsoTight_TurnOn50;
  TH1D* hRelIsoTight_TurnOn50_1pr;
  TH1D* hRelIsoTight_TurnOn50_3pr;
  TH1D* hRelIsoTight_TurnOn50_withNeutrals;
  TH1D* hRelIsoTight_TurnOn50_noNeutrals;

};

#endif

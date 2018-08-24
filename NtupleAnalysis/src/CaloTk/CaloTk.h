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
// #include "../Auxiliary/src/L1Tracks.C" // needed?
#include "../Auxiliary/src/Datasets.C" 

#include "../DataFormat/src/L1TkTauParticle.C"
#include "../DataFormat/src/GenParticle.C"
#include "../DataFormat/src/TrackingParticle.C"
#include "../DataFormat/interface/TTTrack.h"
#include "../DataFormat/interface/TTPixelTrack.h"
#include "../DataFormat/src/L1EG.C"
#include "../DataFormat/src/L1Jet.C"
#include "../DataFormat/src/L1Tau.C"
#include "../DataFormat/src/L1Sum.C"
#include "../DataFormat/src/L1CaloTP.C"

// #include "../Plugins/src/L1TkPrimaryVertex.C"
#include "../Plugins/src/L1PixelTrackFit.C"

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

  void ApplyDiTauZMatching(string tkCollectionType,
			   vector<L1TkTauParticle> &L1TkTaus);

  void GetShrinkingConeSizes(double calo_et,
			     double sigCone_Constant,
			     double isoCone_Constant,
			     const double sigCone_dRCutoff,
			     double &sigCone_dRMin,
			     double &sigCone_dRMax,
			     double &isoCone_dRMin,
			     double &isoCone_dRMax);

  void GetMatchingTrack(L1TkTauParticle &L1TkTau,
			L1Tau L1CaloTau,
			vector<TTTrack> TTTracks);
  
  void GetSigConeTracks(L1TkTauParticle &L1TkTau,
			vector<TTTrack> TTTracks);
  
  void GetIsoConeTracks(L1TkTauParticle &L1TkTau,
			vector<TTTrack> TTTracks);
    
  void GetIsolationValues(L1TkTauParticle &L1TkTau);
  
  void GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
			      vector<GenParticle> hadGenTaus);			    

  // Public Variables
  bool DEBUG;
  bool mcMatching_unique;
  double diTau_deltaPOCAz;

  // L1TkTau - Matching track
  string matchTk_Collection;
  int matchTk_nFitParams;
  double matchTk_minPt;
  double matchTk_minEta;
  double matchTk_maxEta;
  double matchTk_maxChiSqRed;
  double matchTk_minStubs;
  double matchTk_caloDeltaR;
  // Signal Cone Tracks
  string sigConeTks_Collection;
  int sigConeTks_nFitParams;
  double sigConeTks_minPt;
  double sigConeTks_minEta;
  double sigConeTks_maxEta;
  double sigConeTks_maxChiSqRed;
  unsigned int sigConeTks_minStubs;
  // Isolation Cone Tracks
  string isoConeTks_Collection;
  int isoConeTks_nFitParams;
  double isoConeTks_minPt;
  double isoConeTks_minEta;
  double isoConeTks_maxEta;
  double isoConeTks_maxChiSqRed;
  unsigned int isoConeTks_minStubs;

  double mcMatching_dRMax;
  double pv_deltaZMax;
  double pv_z;
  
  double _eta_C;
  double _eta_F;

  //
  double sigCone_Constant;
  double sigCone_cutoffDeltaR;
  double sigCone_dRMax;
  double sigCone_maxTkDeltaPOCAz;
  double sigCone_maxTkInvMass;
  double sigCone_dRMin;
  //
  double isoCone_Constant;
  double isoCone_VtxIsoWP;
  double isoCone_RelIsoWP;
  double isoCone_dRMax;
  double isoCone_dRMin;
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

  void FillTurnOn_Numerator_(vector<L1TkTauParticle> L1TkTaus,
			     const double minEt,
			     TH1D *hTurnOn);

  void FillSingleTau_(vector<L1TkTauParticle> L1TkTaus,
		      TH1D *hRate,
		      TH1D *hEfficiency,
		      double minEta=0.0,
		      double maxEta=999.9);

  void FillDiTau_(vector<L1TkTauParticle> L1TkTaus, 
		  TH1D *hRate,
		  TH1D *hEfficiency,
		  double minEta=0.0,
		  double maxEta=999.9);

  void FillDiTau_(vector<L1TkTauParticle> L1TkTaus1,
		  vector<L1TkTauParticle> L1TkTaus2,
		  TH2D *hRate,
		  TH2D *hEfficiency);

  vector<L1TkTauParticle> GetMcMatchedL1TkTaus(vector<L1TkTauParticle> L1TkTaus);

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
  TH2D* h_GenP_VisET_dRMaxLdgPion;
  
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

  // L1TkTaus: Matching track
  TH1D* hL1TkTau_MatchTk_DeltaR;
  TH1D* hL1TkTau_MatchTk_PtRel;
  TH1D* hL1TkTau_MatchTk_Pt;
  TH1D* hL1TkTau_MatchTk_Eta;
  TH1D* hL1TkTau_MatchTk_POCAz;
  // TH1D* hL1TkTau_MatchTk_d0;
  // TH1D* hL1TkTau_MatchTk_d0Abs;
  TH1D* hL1TkTau_MatchTk_NStubs;
  TH1D* hL1TkTau_MatchTk_NPsStubs;
  TH1D* hL1TkTau_MatchTk_NBarrelStubs;
  TH1D* hL1TkTau_MatchTk_NEndcapStubs;
  TH1D* hL1TkTau_MatchTk_ChiSquared;
  TH1D* hL1TkTau_MatchTk_RedChiSquared;
  TH1D* hL1TkTau_MatchTk_IsGenuine;
  TH1D* hL1TkTau_MatchTk_IsUnknown;
  TH1D* hL1TkTau_MatchTk_IsCombinatoric;
  TH1D* hL1TkTau_MatchTk_PtMinusCaloEt;

  // L1TkTaus: Signal cone tracks
  TH1D* hL1TkTau_SigTks_Pt;
  TH1D* hL1TkTau_SigTks_PtRel;
  TH1D* hL1TkTau_SigTks_Eta;
  TH1D* hL1TkTau_SigTks_POCAz;
  TH1D* hL1TkTau_SigTks_DeltaPOCAz;
  // TH1D* hL1TkTau_SigTks_d0;
  // TH1D* hL1TkTau_SigTks_d0Abs;
  // TH1D* hL1TkTau_SigTks_d0Sig;
  // TH1D* hL1TkTau_SigTks_d0SigAbs;
  TH1D* hL1TkTau_SigTks_DeltaR;
  TH1D* hL1TkTau_SigTks_NStubs;
  TH1D* hL1TkTau_SigTks_NPsStubs;
  TH1D* hL1TkTau_SigTks_NBarrelStubs;
  TH1D* hL1TkTau_SigTks_NEndcapStubs;
  TH1D* hL1TkTau_SigTks_ChiSquared;
  TH1D* hL1TkTau_SigTks_RedChiSquared;
  TH1D* hL1TkTau_SigTks_PtMinusCaloEt;

  // L1TkTaus: Isolation cone tracks
  TH1D* hL1TkTau_IsoTks_Pt;
  TH1D* hL1TkTau_IsoTks_PtRel;
  TH1D* hL1TkTau_IsoTks_Eta;
  TH1D* hL1TkTau_IsoTks_POCAz;
  TH1D* hL1TkTau_IsoTks_DeltaPOCAz;
  // TH1D* hL1TkTau_IsoTks_d0;
  // TH1D* hL1TkTau_IsoTks_d0Abs;
  // TH1D* hL1TkTau_IsoTks_d0Sig;
  // TH1D* hL1TkTau_IsoTks_d0SigAbs;
  TH1D* hL1TkTau_IsoTks_DeltaR;
  TH1D* hL1TkTau_IsoTks_NStubs;
  TH1D* hL1TkTau_IsoTks_NPsStubs;
  TH1D* hL1TkTau_IsoTks_NBarrelStubs;
  TH1D* hL1TkTau_IsoTks_NEndcapStubs;
  TH1D* hL1TkTau_IsoTks_ChiSquared;
  TH1D* hL1TkTau_IsoTks_RedChiSquared;
  TH1D* hL1TkTau_IsoTks_PtMinusCaloEt;

  // L1TkTaus: VtxIsolated
  TH1D* hL1TkTau_Multiplicity;
  TH1D* hL1TkTau_CaloEt; 
  TH1D* hL1TkTau_CaloEta;
  TH1D* hL1TkTau_CaloPhi;
  TH1D* hL1TkTau_CaloIEt;
  TH1D* hL1TkTau_CaloIEta; // ieta of seed tower
  TH1D* hL1TkTau_CaloIPhi; // iphi of seed tower
  TH1D* hL1TkTau_CaloIso;
  TH1D* hL1TkTau_CaloTowerIEta;
  TH1D* hL1TkTau_CaloTowerIPhi;
  TH1D* hL1TkTau_CaloRawEt; // raw (uncalibrated) cluster sum
  TH1D* hL1TkTau_CaloIsoEt; // raw isolation sum - cluster sum
  TH1D* hL1TkTau_CaloNTT;   // n towers above threshold
  TH1D* hL1TkTau_CaloHasEM;
  TH1D* hL1TkTau_CaloIsMerged;

  TH1D* hL1TkTau_Rtau;
  TH1D* hL1TkTau_CHF;
  TH1D* hL1TkTau_NHF;
  TH1D* hL1TkTau_NHFAbs;
  TH1D* hL1TkTau_NSigTks;
  TH1D* hL1TkTau_SigTksEt;
  TH1D* hL1TkTau_SigTksEta;
  TH1D* hL1TkTau_NIsoTks;
  TH1D* hL1TkTau_IsoTksEt;
  TH1D* hL1TkTau_IsoTksEta;
  TH1D* hL1TkTau_InvMass;
  TH1D* hL1TkTau_InvMassIncl;
  TH1D* hL1TkTau_SigConeRMin;
  TH1D* hL1TkTau_SigConeRMax;
  TH1D* hL1TkTau_IsoConeRMin;
  TH1D* hL1TkTau_IsoConeRMax;
  TH1D* hL1TkTau_Charge;
  TH1D* hL1TkTau_RelIso;
  TH1D* hL1TkTau_VtxIso;
  TH1D* hL1TkTau_VtxIsoAbs;
  TH1D* hL1TkTau_DeltaRGenP;

  // L1Taus: Resolutions
  TH1D* hL1Tau_ResolutionCaloEt;
  TH1D* hL1Tau_ResolutionCaloEta;
  TH1D* hL1Tau_ResolutionCaloPhi;

  // L1TkTaus: Resolutions
  TH1D* hL1TkTau_ResolutionCaloEt;
  TH1D* hL1TkTau_ResolutionCaloEta;
  TH1D* hL1TkTau_ResolutionCaloPhi;
  TH1D* hL1TkTau_ResolutionCaloEt_C;
  TH1D* hL1TkTau_ResolutionCaloEta_C;
  TH1D* hL1TkTau_ResolutionCaloPhi_C;
  TH1D* hL1TkTau_ResolutionCaloEt_I;
  TH1D* hL1TkTau_ResolutionCaloEta_I;
  TH1D* hL1TkTau_ResolutionCaloPhi_I;
  TH1D* hL1TkTau_ResolutionCaloEt_F;
  TH1D* hL1TkTau_ResolutionCaloEta_F;
  TH1D* hL1TkTau_ResolutionCaloPhi_F;
  
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
  // SingleTau: Efficiencies
  TH1D* hCalo_Eff;  // Inclusive = C+I+F
  TH1D* hCalo_Eff_C;
  TH1D* hCalo_Eff_I;
  TH1D* hCalo_Eff_F;
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

  // DiTau: Rates
  TH1D* hDiTau_Rate_Calo; // Inclusive = C+I+F
  TH1D* hDiTau_Rate_Calo_C;
  TH1D* hDiTau_Rate_Calo_I;
  TH1D* hDiTau_Rate_Calo_F;
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

  // DiTau: Efficiencies
  TH1D* hDiTau_Eff_Calo; // Inclusive = C+I+F
  TH1D* hDiTau_Eff_Calo_C;
  TH1D* hDiTau_Eff_Calo_I;
  TH1D* hDiTau_Eff_Calo_F;
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

  // DiTau: (Calo-Other)
  TH2D* hDiTau_Rate_Calo_Tk;
  TH2D* hDiTau_Rate_Calo_VtxIso;
  TH2D* hDiTau_Rate_Calo_RelIso;
  TH2D* hDiTau_Eff_Calo_Tk;
  TH2D* hDiTau_Eff_Calo_VtxIso;
  TH2D* hDiTau_Eff_Calo_RelIso;

  // DiTau (Tk-Other)
  TH2D* hDiTau_Rate_Tk_VtxIso;
  TH2D* hDiTau_Rate_Tk_RelIso;
  TH2D* hDiTau_Eff_Tk_VtxIso;
  TH2D* hDiTau_Eff_Tk_RelIso;

  // Turn-Ons
  // TEfficiency* pEff;
  TH1D* hMcHadronicTau_VisEt;
  TH1D* hCalo_TurnOn50;
  TH1D* hTk_TurnOn50;
  TH1D* hVtxIso_TurnOn50;
  TH1D* hRelIso_TurnOn50;
  //
  TH1D* hCalo_TurnOn25;
  TH1D* hTk_TurnOn25;
  TH1D* hVtxIso_TurnOn25;
  TH1D* hRelIso_TurnOn25;
  //
  TH1D* hCalo_TurnOn_SingleTau50KHz;
  TH1D* hTk_TurnOn_SingleTau50KHz;
  TH1D* hVtxIso_TurnOn_SingleTau50KHz;
  TH1D* hRelIso_TurnOn_SingleTau50KHz;
  //
  TH1D* hCalo_TurnOn_DiTau50KHz;
  TH1D* hTk_TurnOn_DiTau50KHz;
  TH1D* hVtxIso_TurnOn_DiTau50KHz;
  TH1D* hRelIso_TurnOn_DiTau50KHz;

};

#endif

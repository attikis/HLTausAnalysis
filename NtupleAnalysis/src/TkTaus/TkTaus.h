#ifndef TkTaus_h
#define TkTaus_h

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
#include "../DataFormat/interface/TTPixelTrack.h"
#include "../DataFormat/src/L1EG.C"
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

class TkTaus : public TreeAnalyserMC{
 public:
  // Constructors/Destructors
  ~TkTaus(){};
 TkTaus(const string SamplePath,
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

  void GetShrinkingConeSizes(double tk_pt,
			     double sigCone_Constant,
			     double isoCone_Constant,
			     const double sigCone_dRCutoff,
			     double &sigCone_dRMin,
			     double &sigCone_dRMax,
			     double &isoCone_dRMin,
			     double &isoCone_dRMax);

  double GetDonutRatio(L1TkTauParticle &L1TkTau, vector<TTTrack> isoTTTracks);


  double GetJetWidth(vector<TTTrack> sigTks, vector<TTTrack> isoTks,
		     TLorentzVector sigTks_p4, TLorentzVector isoTks_p4);

  void GetSigConeTracks(L1TkTauParticle &L1TkTau,
			vector<TTTrack> TTTracks,
			double sigConeTks_dPOCAz,
			double sigConeTks_invMass);
  
  void GetIsoConeTracks(L1TkTauParticle &L1TkTau,
			vector<TTTrack> TTTracks,
			double isoConeTks_dPOCAz);
    
  void GetIsolationValues(L1TkTauParticle &L1TkTau);
  
  void GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
			      vector<GenParticle> hadGenTaus);			    

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
  double vtxIso_WP;
  double relIso_WP;
  double relIso_dZ0;
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
  TH2D* hGenP_VisEt_Vs_dRMaxLdgPion;
  TH2D* hGenP_PtLdg_Vs_dRMaxLdgPion;

  // Counters
  TH1D* hCounters;

  // L1TkTaus
  TH1D* hL1TkTau_SeedTk_DeltaR;
  TH1D* hL1TkTau_SeedTk_PtRel;
  TH1D* hL1TkTau_SeedTk_Pt;
  TH1D* hL1TkTau_SeedTk_Eta;
  TH1D* hL1TkTau_SeedTk_POCAz;
  TH1D* hL1TkTau_SeedTk_NStubs;
  TH1D* hL1TkTau_SeedTk_NPsStubs;
  TH1D* hL1TkTau_SeedTk_NBarrelStubs;
  TH1D* hL1TkTau_SeedTk_NEndcapStubs;
  TH1D* hL1TkTau_SeedTk_ChiSquared;
  TH1D* hL1TkTau_SeedTk_RedChiSquared;
  TH1D* hL1TkTau_SeedTk_IsGenuine;
  TH1D* hL1TkTau_SeedTk_IsUnknown;
  TH1D* hL1TkTau_SeedTk_IsCombinatoric;

  TH1D* hL1TkTau_SigTks_Pt;
  TH1D* hL1TkTau_SigTks_PtRel;
  TH1D* hL1TkTau_SigTks_Eta;
  TH1D* hL1TkTau_SigTks_POCAz;
  TH1D* hL1TkTau_SigTks_DeltaPOCAz;
  TH1D* hL1TkTau_SigTks_DeltaR;
  TH1D* hL1TkTau_SigTks_NStubs;
  TH1D* hL1TkTau_SigTks_NPsStubs;
  TH1D* hL1TkTau_SigTks_NBarrelStubs;
  TH1D* hL1TkTau_SigTks_NEndcapStubs;
  TH1D* hL1TkTau_SigTks_ChiSquared;
  TH1D* hL1TkTau_SigTks_RedChiSquared;

  TH1D* hL1TkTau_IsoTks_Pt;
  TH1D* hL1TkTau_IsoTks_PtRel;
  TH1D* hL1TkTau_IsoTks_Eta;
  TH1D* hL1TkTau_IsoTks_POCAz;
  TH1D* hL1TkTau_IsoTks_DeltaPOCAz;
  TH1D* hL1TkTau_IsoTks_DeltaR;
  TH1D* hL1TkTau_IsoTks_NStubs;
  TH1D* hL1TkTau_IsoTks_NPsStubs;
  TH1D* hL1TkTau_IsoTks_NBarrelStubs;
  TH1D* hL1TkTau_IsoTks_NEndcapStubs;
  TH1D* hL1TkTau_IsoTks_ChiSquared;
  TH1D* hL1TkTau_IsoTks_RedChiSquared;

  TH1D* hL1TkTau_Multiplicity;
  TH1D* hL1TkTau_Multiplicity_MC;
  TH1D* hL1TkTau_JetWidth;
  TH1D* hL1TkTau_DonutRatio;
  TH2D* hL1TkTau_DonutRatio_Vs_JetWidth;
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
  TH2D* hL1TkTau_VtxIso_Vs_RelIso;
  TH1D* hL1TkTau_DeltaRGenP;

  // L1TkIsoTaus
  TH1D* hL1TkIsoTau_SeedTk_DeltaR;
  TH1D* hL1TkIsoTau_SeedTk_PtRel;
  TH1D* hL1TkIsoTau_SeedTk_Pt;
  TH1D* hL1TkIsoTau_SeedTk_Eta;
  TH1D* hL1TkIsoTau_SeedTk_POCAz;
  TH1D* hL1TkIsoTau_SeedTk_NStubs;
  TH1D* hL1TkIsoTau_SeedTk_NPsStubs;
  TH1D* hL1TkIsoTau_SeedTk_NBarrelStubs;
  TH1D* hL1TkIsoTau_SeedTk_NEndcapStubs;
  TH1D* hL1TkIsoTau_SeedTk_ChiSquared;
  TH1D* hL1TkIsoTau_SeedTk_RedChiSquared;
  TH1D* hL1TkIsoTau_SeedTk_IsGenuine;
  TH1D* hL1TkIsoTau_SeedTk_IsUnknown;
  TH1D* hL1TkIsoTau_SeedTk_IsCombinatoric;
  TH1D* hL1TkIsoTau_SigTks_Pt;
  TH1D* hL1TkIsoTau_SigTks_PtRel;
  TH1D* hL1TkIsoTau_SigTks_Eta;
  TH1D* hL1TkIsoTau_SigTks_POCAz;
  TH1D* hL1TkIsoTau_SigTks_DeltaPOCAz;
  TH1D* hL1TkIsoTau_SigTks_DeltaR;
  TH1D* hL1TkIsoTau_SigTks_NStubs;
  TH1D* hL1TkIsoTau_SigTks_NPsStubs;
  TH1D* hL1TkIsoTau_SigTks_NBarrelStubs;
  TH1D* hL1TkIsoTau_SigTks_NEndcapStubs;
  TH1D* hL1TkIsoTau_SigTks_ChiSquared;
  TH1D* hL1TkIsoTau_SigTks_RedChiSquared;
  TH1D* hL1TkIsoTau_IsoTks_Pt;
  TH1D* hL1TkIsoTau_IsoTks_PtRel;
  TH1D* hL1TkIsoTau_IsoTks_Eta;
  TH1D* hL1TkIsoTau_IsoTks_POCAz;
  TH1D* hL1TkIsoTau_IsoTks_DeltaPOCAz;
  TH1D* hL1TkIsoTau_IsoTks_DeltaR;
  TH1D* hL1TkIsoTau_IsoTks_NStubs;
  TH1D* hL1TkIsoTau_IsoTks_NPsStubs;
  TH1D* hL1TkIsoTau_IsoTks_NBarrelStubs;
  TH1D* hL1TkIsoTau_IsoTks_NEndcapStubs;
  TH1D* hL1TkIsoTau_IsoTks_ChiSquared;
  TH1D* hL1TkIsoTau_IsoTks_RedChiSquared;

  TH1D* hL1TkIsoTau_Multiplicity;
  TH1D* hL1TkIsoTau_Multiplicity_MC;
  TH1D* hL1TkIsoTau_JetWidth;
  TH1D* hL1TkIsoTau_DonutRatio;
  TH2D* hL1TkIsoTau_DonutRatio_Vs_JetWidth;
  TH1D* hL1TkIsoTau_NSigTks;
  TH1D* hL1TkIsoTau_SigTksEt;
  TH1D* hL1TkIsoTau_SigTksEta;
  TH1D* hL1TkIsoTau_NIsoTks;
  TH1D* hL1TkIsoTau_IsoTksEt;
  TH1D* hL1TkIsoTau_IsoTksEta;
  TH1D* hL1TkIsoTau_InvMass;
  TH1D* hL1TkIsoTau_InvMassIncl;
  TH1D* hL1TkIsoTau_SigConeRMin;
  TH1D* hL1TkIsoTau_SigConeRMax;
  TH1D* hL1TkIsoTau_IsoConeRMin;
  TH1D* hL1TkIsoTau_IsoConeRMax;
  TH1D* hL1TkIsoTau_Charge;
  TH1D* hL1TkIsoTau_RelIso;
  TH1D* hL1TkIsoTau_VtxIso;
  TH2D* hL1TkIsoTau_VtxIso_Vs_RelIso;
  TH1D* hL1TkIsoTau_DeltaRGenP;

  // Resolutions
  TH1D* hL1TkIsoTau_ResolutionEt;
  TH1D* hL1TkIsoTau_ResolutionEta;
  TH1D* hL1TkIsoTau_ResolutionPhi;
  TH1D* hL1TkIsoTau_ResolutionEt_withNeutrals;
  TH1D* hL1TkIsoTau_ResolutionEta_withNeutrals;
  TH1D* hL1TkIsoTau_ResolutionPhi_withNeutrals;
  TH1D* hL1TkIsoTau_ResolutionEt_noNeutrals;
  TH1D* hL1TkIsoTau_ResolutionEta_noNeutrals;
  TH1D* hL1TkIsoTau_ResolutionPhi_noNeutrals;
  TH1D* hL1TkIsoTau_ResolutionEt_C;
  TH1D* hL1TkIsoTau_ResolutionEta_C;
  TH1D* hL1TkIsoTau_ResolutionPhi_C;
  TH1D* hL1TkIsoTau_ResolutionEt_I;
  TH1D* hL1TkIsoTau_ResolutionEta_I;
  TH1D* hL1TkIsoTau_ResolutionPhi_I;
  TH1D* hL1TkIsoTau_ResolutionEt_F;
  TH1D* hL1TkIsoTau_ResolutionEta_F;
  TH1D* hL1TkIsoTau_ResolutionPhi_F;
  
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

  TH1D* hCalo_TurnOn50;
  TH1D* hTk_TurnOn50;
  TH1D* hVtxIso_TurnOn50;
  TH1D* hRelIso_TurnOn50;
  TH1D* hVtxIsoLoose_TurnOn50;
  TH1D* hVtxIsoTight_TurnOn50;
  TH1D* hRelIsoLoose_TurnOn50;
  TH1D* hRelIsoTight_TurnOn50;

  TH1D* hCalo_TurnOn25;
  TH1D* hTk_TurnOn25;
  TH1D* hVtxIso_TurnOn25;
  TH1D* hRelIso_TurnOn25;
  TH1D* hVtxIsoLoose_TurnOn25;
  TH1D* hVtxIsoTight_TurnOn25;
  TH1D* hRelIsoLoose_TurnOn25;
  TH1D* hRelIsoTight_TurnOn25;

};

#endif

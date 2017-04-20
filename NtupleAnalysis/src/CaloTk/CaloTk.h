#ifndef CaloTk_h
#define CaloTk_h

// System 
#include <iostream>
#include <stdlib.h> 

// User
#include "../Framework/src/TreeAnalyserMC.C"

#include "../Auxiliary/src/AuxTools.C"
#include "../Auxiliary/src/Table.C"
#include "../Auxiliary/src/MCTools.C"
#include "../Auxiliary/src/HistoTools.C"
#include "../Auxiliary/src/L1Tracks.C" // needed?
#include "../Auxiliary/src/Datasets.C" 

#include "../DataFormat/src/L1TkTauParticle.C"
#include "../DataFormat/src/GenParticle.C"
#include "../DataFormat/src/TrackingParticle.C"
#include "../DataFormat/interface/TTTrack.h"
#include "../DataFormat/interface/TTPixelTrack.h"
#include "../DataFormat/src/L1JetParticle.C"

#include "../Plugins/src/L1TkPrimaryVertex.C"
#include "../Plugins/src/L1PixelTrackFit.C"

// ROOT
#include "TEfficiency.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"


using namespace std;

class CaloTk : public TreeAnalyserMC{
 public:
  // Constructors/Destructors
  ~CaloTk(){};
 CaloTk(const string SamplePath,
		const string SampleName,
		const string text_, 
		const int maxEvents_ = -1, 
		TTree* tree=0) : 
  
  TreeAnalyserMC("CaloTk", SamplePath, SampleName, text_, maxEvents_, tree) 
    { 
      auxTools_.StopwatchStart();
      mcSample = SampleName;
      InitObjects();
    };
  
  // Public Variables
  virtual void Loop();

  void PrintSettings(void);

  void PrintGenParticleCollection(vector<GenParticle> collection);

  void PrintTrackingParticleCollection(vector<TrackingParticle> collection);

  void PrintTTTrackCollection(vector<TTTrack> collection);

  void PrintTTPixelTrackCollection(vector<TTPixelTrack> collection);

  void PrintL1JetParticleCollection(vector<L1JetParticle> collection);

  void PrintL1TkTauParticleCollection(vector<L1TkTauParticle> collection);

  void ApplyDiTauZMatching(string tkCollectionType,
			   vector<L1TkTauParticle> &L1TkTaus);

  vector<GenParticle> GetHadronicGenTaus(vector<GenParticle> GenTaus,
					 double visEt=20.0,
					 double visEta=2.3);

  void GetShrinkingConeSizes(double calo_et,
			     double sigCone_Constant,
			     double isoCone_Constant,
			     const double sigCone_dRCutoff,
			     double &sigCone_dRMin,
			     double &sigCone_dRMax,
			     double &isoCone_dRMin,
			     double &isoCone_dRMax);

  void GetMatchingTrack(L1TkTauParticle &L1TkTau,
			L1JetParticle L1CaloTau,
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
  double matchTk_minStubsPS;
  double matchTk_maxStubsPS;
  double matchTk_caloDeltaR;
  // Signal Cone Tracks
  string sigConeTks_Collection;
  int sigConeTks_nFitParams;
  double sigConeTks_minPt;
  double sigConeTks_minEta;
  double sigConeTks_maxEta;
  double sigConeTks_maxChiSqRed;
  unsigned int sigConeTks_minStubs;
  unsigned int sigConeTks_minStubsPS;
  unsigned int sigConeTks_maxStubsPS;
  // Isolation Cone Tracks
  string isoConeTks_Collection;
  int isoConeTks_nFitParams;
  double isoConeTks_minPt;
  double isoConeTks_minEta;
  double isoConeTks_maxEta;
  double isoConeTks_maxChiSqRed;
  unsigned int isoConeTks_minStubs;
  unsigned int isoConeTks_minStubsPS;
  unsigned int isoConeTks_maxStubsPS;

  double mcMatching_dRMax;
  double pv_deltaZMax;
  double pv_z;

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

  TrackingParticle GetTrackingParticle(unsigned int Index);
  vector<TrackingParticle> GetTrackingParticles(bool bPrintList=false);

  TTTrack GetTTTrack(unsigned int Index,
		     const unsigned int nFitParams = 5);
  
  vector<TTTrack> GetTTTracks(const double minPt = 0.0,
			      const double minEta = 0.0,
			      const double maxEta = 9999.9,
			      const double maxChiSqRed = 9999.9,
			      const unsigned int minStubs = 0,
			      const unsigned int minStubsPS = 0,
			      const unsigned int maxStubsPS = 999,
			      const unsigned nFitParams = 5,
			      bool bPrintList = false);

  double GetPVTTTracks(vector<TTTrack> &pvTTTracks,
		       bool bPrintList = false);

  TTPixelTrack GetTTPixelTrack(unsigned int Index);

  vector<TTPixelTrack> GetTTPixelTracks(const double minPt = 0.0,
					const double maxEta = 9999.9,
					const double maxChiSqRed = 9999.9,
					const int minHits = 0.0,
					bool bPrintList=false);

  L1JetParticle GetL1CaloTau(unsigned int Index);
  
  vector<L1JetParticle> GetL1CaloTaus(bool bPrintList=false);

  GenParticle GetGenParticle(unsigned int Index);

  void SetGenParticleMomsAndDaus(GenParticle &p);

  void SetGenParticleFinalDaughters(GenParticle &p);

  vector<GenParticle> GetGenParticles(bool bPrintList=false);

  vector<GenParticle> GetGenParticles(int pdgId,
				      bool isLastCopy=false);
  
  vector<L1TkTauParticle> GetMcMatchedL1TkTaus(vector<L1TkTauParticle> L1TkTaus);

  void GetHadronicTauFinalDaughters(GenParticle hTau,
				    vector<unsigned short> &Daug);

  
  // Variable declaration
  L1TkPrimaryVertex *pvProducer;
  AuxTools auxTools_;
  MCTools mcTools_;
  Datasets datasets_;
  HistoTools histoTools_;
  bool bFoundAllTaus_;

  // Event-Type Histograms
  TH1D* hHepMCEvt_VtxZ;
  TH2D* hHepMCEvt_VtxX_VtxY;
  TH1D* hL1TkPV_VtxZ;
  TH1D* hPrimaryVertex_DeltaZ;
  TH1D* hPrimaryVertex_AbsDeltaZ;  
  
  // L1TkTaus: Matching track
  TH1D* hL1TkTau_MatchTk_DeltaR;
  TH1D* hL1TkTau_MatchTk_PtRel;
  TH1D* hL1TkTau_MatchTk_Pt;
  TH1D* hL1TkTau_MatchTk_Eta;
  TH1D* hL1TkTau_MatchTk_POCAz;
  TH1D* hL1TkTau_MatchTk_d0;
  TH1D* hL1TkTau_MatchTk_d0Abs;
  TH1D* hL1TkTau_MatchTk_NStubs;
  TH1D* hL1TkTau_MatchTk_NPsStubs;
  TH1D* hL1TkTau_MatchTk_NBarrelStubs;
  TH1D* hL1TkTau_MatchTk_NEndcapStubs;
  TH1D* hL1TkTau_MatchTk_StubPtCons;
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
  TH1D* hL1TkTau_SigTks_d0;
  TH1D* hL1TkTau_SigTks_d0Abs;
  TH1D* hL1TkTau_SigTks_d0Sig;
  TH1D* hL1TkTau_SigTks_d0SigAbs;
  TH1D* hL1TkTau_SigTks_StubPtCons;
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
  TH1D* hL1TkTau_IsoTks_d0;
  TH1D* hL1TkTau_IsoTks_d0Abs;
  TH1D* hL1TkTau_IsoTks_d0Sig;
  TH1D* hL1TkTau_IsoTks_d0SigAbs;
  TH1D* hL1TkTau_IsoTks_StubPtCons;
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
  TH1D* hL1TkTau_Rtau;
  TH1D* hL1TkTau_CHF;
  TH1D* hL1TkTau_NHF;
  TH1D* hL1TkTau_NHFAbs;
  TH1D* hL1TkTau_NSigTks;
  TH1D* hL1TkTau_NIsoTks;
  TH1D* hL1TkTau_InvMass;
  TH1D* hL1TkTau_InvMassIncl;
  TH1D* hL1TkTau_SigConeRMin;
  TH1D* hL1TkTau_SigConeRMax;
  TH1D* hL1TkTau_IsoConeRMin;
  TH1D* hL1TkTau_IsoConeRMax;
  TH1D* hL1TkTau_Charge;
  TH1D* hL1TkTau_ChargeAbs;
  TH1D* hL1TkTau_RelIso;
  TH1D* hL1TkTau_VtxIso;
  TH1D* hL1TkTau_VtxIsoAbs;
  TH1D* hL1TkTau_DeltaRGenP;

  // L1CaloTaus: Resolutions
  TH1D* hL1CaloTau_ResolutionCaloEt;
  TH1D* hL1CaloTau_ResolutionCaloEta;
  TH1D* hL1CaloTau_ResolutionCaloPhi;

  // L1TkTaus: Resolutions
  TH1D* hL1TkTau_ResolutionCaloEt;
  TH1D* hL1TkTau_ResolutionCaloEta;
  TH1D* hL1TkTau_ResolutionCaloPhi;
  
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

  // DiTau: (Calo-Other)
  TH2D* hDiTau_Rate_Calo_Tk;
  TH2D* hDiTau_Rate_Calo_VtxIso;
  TH2D* hDiTau_Eff_Calo_Tk;
  TH2D* hDiTau_Eff_Calo_VtxIso;

  // DiTau (Tk-Other)
  TH2D* hDiTau_Rate_Tk_VtxIso;
  TH2D* hDiTau_Eff_Tk_VtxIso;

  // Turn-Ons
  // TEfficiency* pEff;
  TH1D* hMcHadronicTau_VisEt;
  TH1D* hCalo_TurnOn50;
  TH1D* hTk_TurnOn50;
  TH1D* hVtxIso_TurnOn50;
  //
  TH1D* hCalo_TurnOn25;
  TH1D* hTk_TurnOn25;
  TH1D* hVtxIso_TurnOn25;
  //
  TH1D* hCalo_TurnOn_SingleTau50KHz;
  TH1D* hTk_TurnOn_SingleTau50KHz;
  TH1D* hVtxIso_TurnOn_SingleTau50KHz;
  //
  TH1D* hCalo_TurnOn_DiTau50KHz;
  TH1D* hTk_TurnOn_DiTau50KHz;
  TH1D* hVtxIso_TurnOn_DiTau50KHz;

};

#endif

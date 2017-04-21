#ifndef CaloPlusTracks_h
#define CaloPlusTracks_h

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

class CaloPlusTracks : public TreeAnalyserMC{
 public:
  // Constructors/Destructors
  ~CaloPlusTracks(){};
 CaloPlusTracks(const string SamplePath,
		const string SampleName,
		const string text_, 
		const int maxEvents_ = -1, 
		TTree* tree=0) : 
  
  TreeAnalyserMC("CaloPlusTracks", SamplePath, SampleName, text_, maxEvents_, tree) 
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
  double matchTk_caloDeltaR_;
  double matchTk_minPt_;
  double mcMatching_dRMax;
  double pv_deltaZMax;
  double pv_z;
  //
  string selTks_Collection;
  double selTks_minPt;
  double selTks_maxEta;
  double selTks_maxChiSq;
  unsigned int selTks_minStubs;
  unsigned int selTks_minStubsPS;
  unsigned int selTks_nFitParams;
  int selTksPix_minHits;
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
		      TH1D *hEfficiency);

  void FillDiTau_(vector<L1TkTauParticle> L1TkTaus, 
		  TH1D *hRate,
		  TH1D *hEfficiency);

  void FillDiTau_(vector<L1TkTauParticle> L1TkTaus1,
		  vector<L1TkTauParticle> L1TkTaus2,
		  TH2D *hRate,
		  TH2D *hEfficiency);

  TrackingParticle GetTrackingParticle(unsigned int Index);
  vector<TrackingParticle> GetTrackingParticles(bool bPrintList=false);

  TTTrack GetTTTrack(unsigned int Index,
		     const unsigned int nFitParams = 5);
  
  vector<TTTrack> GetTTTracks(const double minPt = 0.0,
			      const double maxEta = 9999.9,
			      const double maxChiSq = 9999.9,
			      const unsigned int minStubs = 0,
			      const unsigned int minStubsPS = 0,
			      const unsigned nFitParams = 5,
			      bool bPrintList = false);

  double GetPVTTTracks(vector<TTTrack> &pvTTTracks,
		       bool bPrintList = false);

  TTPixelTrack GetTTPixelTrack(unsigned int Index);

  vector<TTPixelTrack> GetTTPixelTracks(const double minPt = 0.0,
					const double maxEta = 9999.9,
					const double maxChiSq = 9999.9,
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
  TH1D* hCounters;
  TH1D* hHepMCEvt_VtxZ;
  TH2D* hHepMCEvt_VtxX_VtxY;
  TH1D* hL1TkPV_VtxZ;
  TH1D* hPrimaryVertex_DeltaZ;
  TH1D* hPrimaryVertex_AbsDeltaZ;  

  // VtxIsolated L1TkTaus
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
  TH1D* hL1TkTau_RelIso;
  TH1D* hL1TkTau_VtxIso;
  TH1D* hL1TkTau_VtxIsoAbs;
  TH1D* hL1TkTau_DeltaRGenP;
  
  TH1D* hL1TkTau_SigTks_Pt;
  TH1D* hL1TkTau_SigTks_Eta;
  TH1D* hL1TkTau_SigTks_POCAz;
  TH1D* hL1TkTau_SigTks_DeltaPOCAz;
  TH1D* hL1TkTau_SigTks_d0;
  TH1D* hL1TkTau_SigTks_d0Abs;
  TH1D* hL1TkTau_SigTks_d0Sig;
  TH1D* hL1TkTau_SigTks_d0SigAbs;
  TH1D* hL1TkTau_SigTks_PtRel;
  TH1D* hL1TkTau_SigTks_StubPtCons;

  TH1D* hL1TkTau_IsoTks_Pt;
  TH1D* hL1TkTau_IsoTks_Eta;
  TH1D* hL1TkTau_IsoTks_POCAz;
  TH1D* hL1TkTau_IsoTks_DeltaPOCAz;
  TH1D* hL1TkTau_IsoTks_d0;
  TH1D* hL1TkTau_IsoTks_d0Abs;
  TH1D* hL1TkTau_IsoTks_d0Sig;
  TH1D* hL1TkTau_IsoTks_d0SigAbs;
  TH1D* hL1TkTau_IsoTks_PtRel;
  TH1D* hL1TkTau_IsoTks_StubPtCons;

  // VtxIsolated L1TkTaus, Signal TTPixelTracks
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

  // Resolution
  TH1D* hL1CaloTau_ResolutionCaloEt;
  TH1D* hL1CaloTau_ResolutionCaloEta;
  TH1D* hL1CaloTau_ResolutionCaloPhi;

  TH1D* hL1TkTau_ResolutionCaloEt;
  TH1D* hL1TkTau_ResolutionCaloEta;
  TH1D* hL1TkTau_ResolutionCaloPhi;
  
  // SingleTau: Rates
  TH1D* hCalo_Rate;
  TH1D* hTk_Rate;
  TH1D* hVtxIso_Rate;
  TH1D* hCalo_Eff;
  TH1D* hTk_Eff;
  TH1D* hVtxIso_Eff;

  // DiTau
  TH1D* hDiTau_Rate_Calo;
  TH1D* hDiTau_Rate_Tk;
  TH1D* hDiTau_Rate_VtxIso;
  TH1D* hDiTau_Eff_Calo;
  TH1D* hDiTau_Eff_Tk;
  TH1D* hDiTau_Eff_VtxIso;

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

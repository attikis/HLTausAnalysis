#ifndef PixelTauTrigger_h
#define PixelTauTrigger_h

// System 
#include <iostream>
#include <stdlib.h> 

// User
#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/MCTools.C"
#include "../utilities/AuxTools.C"
#include "../utilities/HistoTools.C"
#include "../utilities/L1Tracks.C"
#include "../utilities/L1TkTauParticle.C"
#include "../utilities/L1TkPrimaryVertex.C"
#include "../utilities/TTPixelTrack.h"
#include "../utilities/L1PixelTrackFit.C"
#include "../utilities/Table.C"
#include "../utilities/Datasets.C"

using namespace std;

class PixelTauTrigger : public TreeAnalyserMC{
 public:
  // Constructors/Destructors
 PixelTauTrigger(const string SamplePath,
	    const string SampleName,
	    const string text_, 
	    const int maxEvents_ = -1, 
	    TTree* tree=0) : 
  
  TreeAnalyserMC("PixelTauTrigger", SamplePath, SampleName, text_, maxEvents_, tree) 
    { 
      auxTools_.StopwatchStart();
      mcSample = SampleName;
      InitObjects();
    };

  ~PixelTauTrigger(){};

  // Function declaration
  virtual void Loop();
  bool EventWithForbiddenTauDecayMode(string allowedTauDecayMode, 
				      vector<int> hTaus_index);

  void GetL1TkTauMatchTrack(int &matchTk_index,
			    const int calo_index,
			    const string matchTk_CollectionType,
			    const double matchTk_minPt,
			    const double matchTk_maxDeltaR,
			    double &matchTk_deltaR);

  void GetL1TkTauMatchPixTrack(int &matchTk_index,
			       const int calo_index,
			       const string matchTk_CollectionType,
			       const double matchTk_minPt,
			       const double matchTk_maxDeltaR,
			       double &matchTk_deltaR);

  void GetL1TkTauMatchPixTrackRefit(int &matchTk_index,
				    const int calo_index,
				    const string matchTk_CollectionType,
				    const double matchTk_minPt,
				    const double matchTk_maxDeltaR,
				    double &matchTk_deltaR,
				    L1PixelTrackFit f);

  void GetL1TkTauIsolation(const int matchTk_index,
			   const vector<int> isoTks_index,
			   double &relIso,
			   double &vtxIso);

  void GetL1TkTauPixIsolation(const int matchTk_index,
			      const vector<int> isoTks_index,
			      double &relIso,
			      double &vtxIso);

  void GetL1TkTauPixIsolationRefit(const int matchTk_index,
				   const vector<int> isoTks_index,
				   double &relIso,
				   double &vtxIso,
				   L1PixelTrackFit f); 

  void GetL1TkTauAlgoConeSizes(const int calo_index,
			       const string algoType,
			       const double sigCone_dRMax, 
			       double &sigCone_minDeltaR_tmp,
			       double &sigCone_maxDeltaR_tmp,
			       double &isoCone_minDeltaR_tmp);

  void GetL1TkTauSigConeTracks(const int ldgTk_index, 
			       const double deltaR_min,
			       const double deltaR_max,
			       const string selectionType,
			       vector<int> &sigTks_index);
  
  void GetL1TkTauSigConePixTracks(const int matchTk_index,
				  const double deltaR_min,
				  const double deltaR_max,
				  const string selectionType,
				  vector<int> &sigTks_index);

  void GetL1TkTauSigConePixTracksRefit(int matchTk_index,
				       const double deltaR_min,
				       const double deltaR_max,
				       const string selectionType,
				       vector<int> &sigTks_index,
				       L1PixelTrackFit f);

  void GetL1TkTauIsoConeTracks(const int ldgTk_index, 
			       const double deltaR_min,
			       const double deltaR_max,
			       const string selectionType,
			       vector<int> &isoTks_index);
  
  void GetL1TkTauIsoConePixTracks(const int ldgTk_index, 
				  const double deltaR_min,
				  const double deltaR_max,
				  const string selectionType,
				  vector<int> &isoTks_index);

    void GetL1TkTauIsoConePixTracksRefit(const int matchTk_index, 
					 const double deltaR_min,
					 const double deltaR_max,
					 const string selectionType,
					 vector<int> &isoTks_index,
					 L1PixelTrackFit f);

  TLorentzVector GetL1TkTauSigTksP4(const string tkCollectionType, 
				    const L1TkTauParticle L1TkTau);

  TLorentzVector GetL1TkTauIsoTksP4(const string tkCollectionType, 
				    const L1TkTauParticle L1TkTau);

  double GetL1TkTauLdgTkPt(const string tkCollectionType, 
			   const L1TkTauParticle L1TkTau);

  void GetL1CaloTauUniqueMatchGenp(vector<L1TkTauParticle> &L1TkTaus);

  void GetL1CaloTauMatchGenp(const int calo_index,
			     vector<int> hTaus_index,
			     int &matchGenp_index, 
			     double &matchGenp_deltaR);

  void GetMcMatchedL1TkTaus(const vector<L1TkTauParticle> L1TkTaus, 
			    vector<L1TkTauParticle> &L1TkTaus_McMatched);
  void PrintSettings(void);

  void PrintL1CaloTau(int Indx);

  void PrintTrackProperties(int tk_Index);

  void PrintL1TkTauProperties(bool bDebug, 
			      L1TkTauParticle L1TkTau);

  void ApplyDiTauZMatching(const string tkCollectionType, 
			   vector<L1TkTauParticle> &L1TkTaus);

  void RemoveDuplicates(const vector<L1TkTauParticle> L1TkTaus1, 
			vector<L1TkTauParticle> &L1TkTaus2);
  
  // Variable declaration
  struct SortAscendingAbs{ bool operator() (double a, double b) const { return abs(a) > abs(b); } }; 
  struct SortDescendingAbs{ bool operator() (double a, double b) const { return abs(a) < abs(b); } }; 
  bool DEBUG;
  bool mcMatching_unique;
  double isoCone_ShrinkConstant;
  double isoCone_VtxIsoWP;
  double isoCone_maxDeltaR;
  double isoCone_minDeltaR;
  double matchTk_caloDeltaR_;
  double matchTk_minPt_;
  double mcMatching_maxDeltaR;
  double pv_deltaZMax;
  double pv_z;
  double diTau_deltaPOCAz;
  double sigCone_ShrinkConstant;
  double sigCone_maxDeltaR;
  double sigCone_maxTkDeltaPOCAz;
  double sigCone_maxTkInvMass;
  int nMaxNumOfHTausPossible;
  int realTauMom;
  string caloCorrectionType;
  string tk_CollectionType;
  string isoCone_tkCollectionWP;
  string mcSample;
  string pv_producer;
  string sigCone_tkCollectionWP;
  string tauAlgorithmMode;
  string tauDecayMode;
  
 private:
  // Function declaration
  void BookHistos_(void);

  void InitObjects(void);

  void InitVars_(void);

  void FillEfficiency_(TH1D *hSignalEfficiency,
		       const Double_t ldgEt);

  void FillEfficiency_(TH2D *hSignalEfficiency,
		       const Double_t ldgEt1,
		       const Double_t ldgEt2);  

  void FillRate_(TH1D *hRate,
		 const Double_t ldgEt);

  void FillRate_(TH2D *hRate,
		 const Double_t ldgEt1,
		 const Double_t ldgEt2);  

  void FinaliseEffHisto_(TH1D *histo, 
			 const int nEvtsTotal);

  void FinaliseEffHisto_(TH2D *histo, 
			 const int nEvtsTotal);  

  void FillTurnOn_Numerator_(const vector<L1TkTauParticle> L1TkTaus, 
			     const double minEt,
			     TH1D *hTurnOn);

  void FillTurnOn_Denominator_(vector<int> hTaus_index,
			       TH1D *hVisEt);

  void FillSingleTau_(const vector<L1TkTauParticle> L1TkTaus,
		      TH1D *hRate,
		      TH1D *hEfficiency);

  void FillDiTau_(const vector<L1TkTauParticle> L1TkTaus, 
		  TH1D *hRate,
		  TH1D *hEfficiency);

  void FillDiTau_(const vector<L1TkTauParticle> L1TkTaus1,
		  const vector<L1TkTauParticle> L1TkTaus2,
		  TH2D *hRate,
		  TH2D *hEfficiency);

  TLorentzVector GetHadronicTauVisP4_(const int genP_index);

  vector<int> GetTriggerHadronicTaus_(const int tauMom,
				      const double hTauVisEtCut);
  

  // Variable declaration
  L1TkPrimaryVertex *pvProducer;
  AuxTools auxTools_;
  Datasets datasets_;
  HistoTools histoTools_;
  bool bFoundAllTaus_;

  // Event-Type Histograms
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

  // VtxIsolated L1TkTaus, Signal TTTracks
  TH1D* hL1TkTau_MatchPixTk_DeltaR;
  TH1D* hL1TkTau_MatchPixTk_PtRel;
  TH1D* hL1TkTau_MatchPixTk_StubPtCons;
  TH1D* hL1TkTau_MatchPixTk_NPixHits;
  TH1D* hL1TkTau_MatchPixTk_NStubs;
  TH1D* hL1TkTau_MatchPixTk_NPsStubs;
  TH1D* hL1TkTau_MatchPixTk_NBarrelStubs;
  TH1D* hL1TkTau_MatchPixTk_NEndcapStubs;
  TH1D* hL1TkTau_MatchPixTk_StubPtConsistency;
  TH1D* hL1TkTau_MatchPixTk_Pt;
  TH1D* hL1TkTau_MatchPixTk_Eta;
  TH1D* hL1TkTau_MatchPixTk_POCAz;
  TH1D* hL1TkTau_MatchPixTk_POCAzSig;
  TH1D* hL1TkTau_MatchPixTk_d0;
  TH1D* hL1TkTau_MatchPixTk_d0Abs;
  TH1D* hL1TkTau_MatchPixTk_d0Sig;
  TH1D* hL1TkTau_MatchPixTk_d0SigAbs;
  TH1D* hL1TkTau_MatchPixTk_ChiSquared;
  TH1D* hL1TkTau_MatchPixTk_RedChiSquared;
  TH1D* hL1TkTau_MatchPixTk_SigmaRInv;
  TH1D* hL1TkTau_MatchPixTk_SigmaPhi0;
  TH1D* hL1TkTau_MatchPixTk_SigmaD0;
  TH1D* hL1TkTau_MatchPixTk_SigmaT;
  TH1D* hL1TkTau_MatchPixTk_SigmaZ0;

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

#ifndef TkEG_h
#define TkEG_h

// System 
#include <iostream>
#include <stdlib.h> 

// User
#include "../Framework/src/TreeAnalyserMC.C"
#include "../Auxiliary/src/AuxTools.C"
#include "../DataFormat/interface/TTTrack.h"
//#include "../DataFormat/src/L1EG.C"
#include "../DataFormat/src/L1TKEM.C"
#include "../DataFormat/src/EG.C"
#include "../DataFormat/src//L1TkEGParticle.C"
#include "../DataFormat/src/TrackingParticle.C" // for GetTTTrack function

//#include "../Auxiliary/src/Table.C"
//#include "../Auxiliary/src/MCTools.C"
#include "../Auxiliary/src/HistoTools.C"
//#include "../Auxiliary/src/Datasets.C" 
//#include "../DataFormat/src/L1TkEGParticle.C"
#include "../DataFormat/src/GenParticle.C"
//#include "../DataFormat/src/TrackingParticle.C"
//#include "../DataFormat/interface/TTTrack.h"
//#include "../DataFormat/interface/TTPixelTrack.h"
//#include "../DataFormat/src/L1Tau.C"
//#include "../DataFormat/src/L1Jet.C"

// #include "../Plugins/src/L1TkPrimaryVertex.C"
//#include "../Plugins/src/L1PixelTrackFit.C"

// ROOT
#include "TEfficiency.h"
#include "Math/VectorUtil.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
//#include "Math/Point3D.h"
#include "Math/Vector4D.h"

using namespace std;

// Class definition
class TkEG : public TreeAnalyserMC{
 public:
  // Destructor
  ~TkEG(){};
  // Constructor
  TkEG(const string SamplePath,
		const string SampleName,
		const string text_, 
		const int maxEvents_ = -1, 
		TChain* tree=0) :  
  TreeAnalyserMC("", SamplePath, SampleName, text_, maxEvents_, tree) 
    { 
      auxTools_.StopwatchStart();
      mcSample = SampleName;
      InitObjects();
    };
    
  // Public function declarations
  virtual void Loop();
  void PrintSettings(void);
  
  // Public variables
  string mcSample;
  bool cfg_AddL1Tks;  
  bool cfg_AddEGs; 
  bool cfg_AddGenP;
  string cfg_tk_Collection;
  int cfg_tk_nFitParams;
  double cfg_tk_minPt;
  double cfg_tk_minEta;
  double cfg_tk_maxEta;
  double cfg_tk_maxChiSqRed;
  double cfg_tk_minStubs;
  double cfg_tk_minStubsPS;
  double cfg_tk_maxStubsPS;
  bool cfg_DEBUG;  
    
 private:
  // Private function declarations
  void BookHistos_(void);
  void WriteHistos_(void);
  void InitObjects(void);
  void InitVars_(void);
  float DeltaPhi(float phi1, float phi2);
  float deltaR(float eta1, float eta2, float phi1, float phi2);
  float CorrectedEta(float eta, float zTrack);
  bool IsWithinEtaRegion(string etaRegion, double eta);
    vector<L1TkEGParticle> GetMcMatchedL1TkEGs(vector<L1TkEGParticle> L1TkEGs);
  double GetMatchingGenParticle(TTTrack track, GenParticle *genParticlePtr);
  void FillTurnOn_Numerator_(vector<L1TkEGParticle> L1TkEGs, const double minEt, TH1D *hTurnOn);

  void FillSingleTau_(vector<L1TkEGParticle> L1TkEGs,
		      TH1D *hRate,
		      TH1D *hEfficiency,
		      double minEta=0.0,
		      double maxEta=999.9);

  void FillDiTau_(vector<L1TkEGParticle> L1TkEGs, 
		  TH1D *hRate,
		  TH1D *hEfficiency,
		  double minEta=0.0,
		  double maxEta=999.9);

  void FillDiTau_(vector<L1TkEGParticle> L1TkEGs1,
		  vector<L1TkEGParticle> L1TkEGs2,
		  TH2D *hRate,
		  TH2D *hEfficiency);

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
				   
  // Private variables
  
  // Old parameters
  float ETmin; // min ET in GeV of L1EG objects
  float ZMAX; // |z_track| < ZMAX in cm
  float CHI2MAX;		
  float DRmin;
  float DRmax;
  float minPt_leading_track;
  bool PrimaryVtxConstrain; // use the primary vertex
  float DeltaZMax; // | z_track - z_primaryvtx | < DeltaZMax in cm. 
  float IsoCut;
  bool RelativeIsolation;

  // New parameters
  unsigned int minStubs_trk;  
  float maxChi2_trk;     
  float maxChi2_trk_alt;
  float minPt_leadtrk;    
  float maxEta_leadtrk;   
  float minDeltaR_leadtrk; 
  float maxDeltaR_leadtrk; 
  float maxDeltaZ_trk; 
  float maxInvMass_trk;  
  float minEt_EG;     
  float minDeltaR_EG;   
  float maxDeltaR_EG;   
  float maxInvMass_EG;  
  float maxDeltaR_MCmatch; 
  float minDeltaR_iso;  
  float maxDeltaR_iso;   
  float maxDeltaZ_iso;   
  float useRelIso;  
  float maxRelIso;  
  float maxVtxIso;

  // Eta regions
  double _eta_C;
  double _eta_F;
  
  vector<TTTrack> TTTracks;
  vector<EG> L1EGs;
  vector< vector <TTTrack> > trackTauCandidates;
  vector<L1TkEGParticle> TauCandidates;
  vector<L1TkEGParticle> TauCandidatesRelIsolated;
  vector<L1TkEGParticle> TauCandidatesVtxIsolated;
  vector<GenParticle> GenTausHadronic;
  vector<GenParticle> GenTausTrigger;

  int nMaxNumOfHTausPossible;  
  bool bFoundAllTaus_;

  // Histogram declarations

  // Event-Type Histograms                                                                                                                                     
  TH1D* h_Counters;
  TH1D* h_Counters_events_EGs;

  TH1D* h_genTausAll_N;
  TH1D* h_genTausAll_Pt;
  TH1D* h_genTausAll_Eta;
  TH1D* h_genTausAll_Phi;
  TH2D* h_genTausAll_Eta1VsEta2;
  TH2D* h_genTausAll_Phi1VsPhi2;

  TH1D* h_genTausHad_Daughters_N;
  TH1D* h_genTausHad_chargedDaugh_N;
  TH1D* h_genTausHad_neutralDaugh_N;
  TH1D* h_genTausHad_N;
  TH1D* h_genTau_chargedDaugh_Pt;
  TH1D* h_genTau_chargedDaugh_totalMass;
  TH1D* h_genTau_neutralDaugh_totalMass;
  TH1D* h_genTau_neutralDaugh_Et;
  TH2D* h_genTauHad_chargedPtVsneutralET;
  TH1D* h_genTau_CHF;
  TH1D* h_genTau_NHF;
  TH2D* h_genTau_chargedDaugh_visPt_dRmax;
  TH2D* h_genTau_chargedDaugh_PtLead_dRmax;
  TH2D* h_genTau_neutralDaugh_PtLead_dRmax;

  TH1D* h_Pion0_Et;
  TH1D* h_Photons_Et;
  TH1D* h_Photons_dR;
  TH2D* h_Pion0Et_Vs_PhotonsDR;
  TH1D* h_Photons_EGs_Matching;

  TH1D* h_trk_Chi2_all_4stubs;
  TH1D* h_trk_Chi2_all_5stubs;
  TH1D* h_leadTrk4stubs_MCmatched_Chi2;
  TH1D* h_trk_Chi2_all;
  TH1D* h_trk_Chi2Red_all;
  TH1D* h_trk_NStubs_all;                                                                                                                                              
  TH2D* h_trk_NStubsVsChi2_all;
  

  TH1D* h_Counters_leadTrks;
  TH1D* h_leadTrks_Multiplicity;
  
  TH1D* h_leadTrks_Pt;
  TH1D* h_leadTrks_Eta;
  TH1D* h_leadTrks_Phi;
  TH2D* h_leadTrks_Phi_Eta;

  TH1D* h_leadTrkSelection;

  TH1D* h_leadTrk_clustTrks_dZ0;

  TH1D* h_Counters_clustTrks;
  TH1D* h_clustTrks_Pt;
  TH1D* h_clustTrks_Eta;
  TH1D* h_clustTrks_Phi;
  TH2D* h_clustTrks_Phi_Eta;
  TH2D* h_clustTrks_counter;

  TH1D* h_trkClusters_MultiplicityPerCluster;
  TH1D* h_MCmatch_chargedDaugh_N;
  TH1D* h_MCmatch_neutralDaugh_N;
  TH1D* h_trkClusters_Pt;
  TH1D* h_trkClusters_M;
  //TH1D* h_trkClusters_M_beforeCut;

  TH1D* h_EGs_N;
  TH1D* h_EGs_N_OneHadTau;
  TH1D* h_EGs_MCmatched_Et;
  TH1D* h_EGs_Et;
  TH1D* h_EGs_Eta;
  TH1D* h_EGs_Phi;
  TH1D* h_EGs_IEta;
  TH1D* h_EGs_IPhi;

  TH1D* h_leadTrk_EG_dR;
  TH1D* h_leadTrk_EG_dR_beforecorrection;
  TH1D* h_leadTrk_EG_dPhi;
  TH1D* h_leadTrk_EG_dEta;
  TH1D* h_leadTrk_EG_dRmin;
  TH1D* h_leadTrk_EG_dPhiMin;
  TH1D* h_leadTrk_EG_dEtaMin;

  TH1D* h_clustEGs_allEGs;
  TH1D* h_clustEGs_passEt;
  TH1D* h_clustEGs_passDRmax;
  TH1D* h_clustEGs_passDRmin;
  TH1D* h_clustEGs_passInvMass;

  TH1D* h_clustEGs_Et;
  TH1D* h_clustEGs_Eta;
  TH1D* h_clustEGs_Phi;
  TH2D* h_clustEGs_Et_Eta;
  TH2D* h_clustEGs_Phi_Eta;
  TH1D* h_clustEGs_M;

  TH1D* h_clustEGs_counter;
  TH1D* h_EGClusters_MultiplicityPerCluster;
  TH1D* h_EGClusters_Et;
  TH1D* h_EGClusters_M;

  TH1D* h_TkEG_relIso;
  TH1D* h_TkEG_vtxIso;

  TH1D* hTkEG_matched_Et;
  TH1D* h_ldgTkEG_ET;
  TH1D* hTkEG_DeltaRmatch;
  TH1D* hTkEG_genVisEt;
  TH1D* hTkEG_genVisEt_clustEG;
  TH1D* hTkEG_genVisPt_clustEG;
  TH1D* hMcHadronicTau_VisEt;

  TH1D* h_TkEG_Pt;
  TH1D* h_TkEG_ET;
  TH1D* h_TkEG_Eta;
  TH1D* h_TkEG_Phi;
  TH1D* h_TkEG_InvMass;
  TH1D* h_TkEG_CHF;
  TH1D* h_TkEG_NHF;

  TH1D* h_TkEG_PtResolution;
  TH1D* h_TkEG_PtResolution_C;
  TH1D* h_TkEG_PtResolution_I;
  TH1D* h_TkEG_PtResolution_F;
  TH1D* h_TkEG_PtResolution_NoNeuDaugh;
  TH1D* h_TkEG_PtResolution_WhenNeuDaugh;
  TH1D* h_TkEG_PtResolution_OneProng;
  TH1D* h_TkEG_PtResolution_ThreeProng;
  
  TH1D* h_TkEG_EtResolution;
  TH1D* h_TkEG_EtResolution_C;
  TH1D* h_TkEG_EtResolution_I;
  TH1D* h_TkEG_EtResolution_F;
  TH1D* h_TkEG_EtResolution_NoNeuDaugh;
  TH1D* h_TkEG_EtResolution_WhenNeuDaugh;
  TH1D* h_TkEG_EtResolution_OneProng;
  TH1D* h_TkEG_EtResolution_ThreeProng;
  
  TH1D* h_TkEG_EtaResolution;
  TH1D* h_TkEG_EtaResolution_C;
  TH1D* h_TkEG_EtaResolution_I;
  TH1D* h_TkEG_EtaResolution_F;
  TH1D* h_TkEG_EtaResolution_NoNeuDaugh;
  TH1D* h_TkEG_EtaResolution_WhenNeuDaugh;
  TH1D* h_TkEG_EtaResolution_OneProng;
  TH1D* h_TkEG_EtaResolution_ThreeProng;

  TH1D* h_TkEG_PhiResolution;
  TH1D* h_TkEG_PhiResolution_C;
  TH1D* h_TkEG_PhiResolution_I;
  TH1D* h_TkEG_PhiResolution_F;
  TH1D* h_TkEG_PhiResolution_NoNeuDaugh;
  TH1D* h_TkEG_PhiResolution_WhenNeuDaugh;
  TH1D* h_TkEG_PhiResolution_OneProng;
  TH1D* h_TkEG_PhiResolution_ThreeProng;

  TH1D* h_nonMCmatched_EGenergyOverTracksPt;
  TH1D* h_nonMCmatchedCandidates_decayMode;

  TH1D* h_MCmatch_counters;
  TH1D* h_MCmatch_dR;
  TH1D* h_leadTrk_MCmatch;
  TH1D* h_leadTrk4stubs_MCmatch;

  TH1D* hTurnOn25_TkEG;
  TH1D* hTurnOn25_relIso;
  TH1D* hTurnOn25_vtxIso;
  TH1D* hTurnOn50_TkEG;
  TH1D* hTurnOn50_relIso;
  TH1D* hTurnOn50_vtxIso;
  
  TH1D* hRateSingleTau_TkEG; // Inclusive = C+I+F
  TH1D* hRateSingleTau_C;
  TH1D* hRateSingleTau_I;
  TH1D* hRateSingleTau_F;

  TH1D* hRateSingleTau_relIso; // Inclusive = C+I+F
  TH1D* hRateSingleTau_vtxIso; // Inclusive = C+I+F

  TH1D* hRateDiTau_TkEG; // Inclusive = C+I+F
  TH1D* hRateDiTau_C;
  TH1D* hRateDiTau_I;
  TH1D* hRateDiTau_F;

  TH1D* hRateDiTau_relIso; // Inclusive = C+I+F
  TH1D* hRateDiTau_vtxIso; // Inclusive = C+I+F

  TH1D* hEffSingleTau_TkEG;  // Inclusive = C+I+F
  TH1D* hEffSingleTau_C;
  TH1D* hEffSingleTau_I;
  TH1D* hEffSingleTau_F;

  TH1D* hEffSingleTau_relIso;  // Inclusive = C+I+F
  TH1D* hEffSingleTau_vtxIso;  // Inclusive = C+I+F

  TH1D* hEffDiTau_TkEG;  // Inclusive = C+I+F
  TH1D* hEffDiTau_C;
  TH1D* hEffDiTau_I;
  TH1D* hEffDiTau_F;

  TH1D* hEffDiTau_relIso;  // Inclusive = C+I+F
  TH1D* hEffDiTau_vtxIso;  // Inclusive = C+I+F
  
};

#endif

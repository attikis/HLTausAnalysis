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
  
  void ApplyDiTauZMatching(vector<L1TkEGParticle> &L1TkEGs);

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
  double cfg_tk_maxChiSq;
  double cfg_tk_minStubs;
  double cfg_eg_minEt;
  double cfg_eg_minEta;
  double cfg_eg_maxEta;

  bool DEBUG;  
    
 private:
  // Private function declarations
  void BookHistos_(void);
  void WriteHistos_(void);
  void InitObjects(void);
  void InitVars_(void);
  //void SortL1TkEGs();
  float DeltaPhi(float phi1, float phi2);
  float deltaR(float eta1, float eta2, float phi1, float phi2);
  float CorrectedEta(float eta, float zTrack);
  bool IsWithinEtaRegion(string etaRegion, double eta);
  vector<L1TkEGParticle> GetMcMatchedL1TkEGs(vector<L1TkEGParticle> L1TkEGs);

  double GetDonutRatio(L1TkEGParticle &L1TkEG, 
		       vector<TTTrack> isoTTTracks,
		       bool bUseCone);
  
  double GetMatchingGenParticle(TTTrack track, GenParticle *genParticlePtr);

  void FillTurnOn_Numerator_(vector<L1TkEGParticle> L1TkEGs,
			     const double minEt,
			     TH1D *hTurnOn, 
			     TH1D *hTurnOn_1pr, 
			     TH1D *hTurnOn_3pr, 
			     TH1D *hTurnOn_withNeutrals, 
			     TH1D *hTurnOn_noNeutrals);

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
  float maxDeltaR_const;
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
  float relIso_WP;  
  float vtxIso_WP;

  double diTau_deltaPOCAz;

  // Eta regions
  double _eta_C;
  double _eta_F;
  
  vector<TTTrack> TTTracks;
  vector<EG> L1EGs;
  vector< vector <TTTrack> > trackTauCandidates;
  vector<L1TkEGParticle> L1TkEGTauCandidates;
  vector<L1TkEGParticle> L1TkEGTaus_RelIso;
  vector<L1TkEGParticle> L1TkEGTaus_VtxIso;
  vector<L1TkEGParticle> L1TkEGTaus_RelIsoLoose;
  vector<L1TkEGParticle> L1TkEGTaus_VtxIsoLoose;
  vector<L1TkEGParticle> L1TkEGTaus_RelIsoTight;
  vector<L1TkEGParticle> L1TkEGTaus_VtxIsoTight;

  vector<GenParticle> GenTausHadronic;
  vector<GenParticle> GenTausTrigger;

  int nMaxNumOfHTausPossible;  
  bool bFoundAllTaus_;

  // Histogram declarations

  // Event-Type Histograms                                                                                                                                     
  TH1D* h_Counters;

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
  TH1D* h_Photons_dEta;
  TH1D* h_Photons_dPhi;
  TH2D* h_Pion0Et_Vs_PhotonsDR;
  TH1D* h_Photons_EGs_Matching;

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

  TH1D *h_SigCone_DeltaR;

  TH1D* h_trkClusters_MultiplicityPerCluster;
  TH1D* h_MCmatch_chargedDaugh_N;
  TH1D* h_MCmatch_neutralDaugh_N;
  TH1D* h_trkClusters_Pt;
  TH1D* h_trkClusters_M;
  TH1D* h_trkClusters_M_beforeCut;

  TH1D* h_EGs_N;
  TH1D* h_EGs_MCmatched_Et;
  TH1D* h_EGs_Et;
  TH1D* h_EGs_Eta;
  TH1D* h_EGs_Phi;
  TH1D* h_EGs_IEta;
  TH1D* h_EGs_IPhi;
  TH1D* h_EGs_HwQual;
  TH2D* h_EGs_EtaVsEt;

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
  TH1D* h_TkEGClusters_M_beforeCut;
  TH1D* h_EGClusters_MultiplicityPerCluster;
  TH1D* h_EGClusters_Et;
  TH1D* h_EGClusters_M;

  TH1D* h_TkEG_relIso;
  TH1D* h_TkEG_vtxIso;

  TH1D* h_ldgTkEG_ET;
  TH1D* hTkEG_DeltaRmatch;
  TH1D* hTkEG_genVisEt;
  TH1D* hTkEG_genVisEt_clustEG;
  TH1D* hTkEG_genVisPt_clustEG;

  TH1D* h_TkEG_N;
  TH1D* h_TkEG_RelIso_N;
  TH1D* h_TkEG_VtxIso_N;
  TH1D* h_TkEG_RelIsoLoose_N;
  TH1D* h_TkEG_VtxIsoLoose_N;
  TH1D* h_TkEG_RelIsoTight_N;
  TH1D* h_TkEG_VtxIsoTight_N;

  TH1D* h_TkEG_Pt;
  TH1D* h_TkEG_ET;
  TH1D* h_TkEG_Eta;
  TH1D* h_TkEG_Phi;
  TH1D* h_TkEG_InvMass;
  TH1D* h_TkEG_NEGs;
  TH1D* h_TkEG_NEGs_withEGs;
  TH1D* h_TkEG_NTks;
  TH1D* h_TkEG_NEGs_C;
  TH1D* h_TkEG_NEGs_I;
  TH1D* h_TkEG_NEGs_F;
  TH1D* h_TkEG_CHF;
  TH1D* h_TkEG_NHF;
  TH1D* h_TkEG_CHF_withNeutrals;
  TH1D* h_TkEG_NHF_withNeutrals;
  TH1D* h_TkEG_clustEGs_MCMatch;
  TH1D* h_TkEG_clustEGs_Matched_HwQual;
  TH1D* h_TkEG_clustEGs_nonMatched_HwQual;
  TH1D* h_TkEG_clustEGs_dET_matchPion0;
  TH1D* h_TkEG_clustEGs_ETResolution;
  TH1D* h_TkEG_isoTracks_InvMass;
  TH1D* h_TkEG_isoTracks_Multiplicity;
  //TH1D* h_TkEG_isoTracks_PtSum;
  //TH1D* h_TkEG_isoTracks_NStubs;
  //TH1D* h_TkEG_isoTracks_Chi2;
  TH1D* h_TkEG_isoTracks_Et;
  TH1D* h_TkEG_isoTracks_Eta;
  TH1D* h_TkEG_DonutRatio;
  TH1D* h_TkEG_signalEGs_Multiplicity;
  TH1D* h_TkEG_isoEGs_Multiplicity;

  TH1D* h_TkEG_PtResolution;
  TH1D* h_TkEG_PtResolution_C;
  TH1D* h_TkEG_PtResolution_I;
  TH1D* h_TkEG_PtResolution_F;
  TH1D* h_TkEG_PtResolution_noNeutrals;
  TH1D* h_TkEG_PtResolution_withNeutrals;
  TH1D* h_TkEG_PtResolution_1pr;
  TH1D* h_TkEG_PtResolution_3pr;

  TH1D* h_TkEG_PtResolution_F_withEGs;
  TH1D* h_TkEG_PtResolution_F_withEGs_posEta;
  TH1D* h_TkEG_PtResolution_F_withEGs_negEta;
  TH1D* h_TkEG_PtResolution_F_noEGs;
  TH1D* h_TkEG_PtResolution_F_noEGs_posEta;
  TH1D* h_TkEG_PtResolution_F_noEGs_negEta;

  TH1D* h_TkEG_PtResolution_noNeutrals_F;
  TH1D* h_TkEG_PtResolution_withNeutrals_F;
  TH1D* h_TkEG_PtResolution_1pr_F;
  TH1D* h_TkEG_PtResolution_3pr_F;

  TH1D* h_TkEG_PtResolution_noNeutrals_F_withEGs;
  TH1D* h_TkEG_PtResolution_withNeutrals_F_withEGs;
  TH1D* h_TkEG_PtResolution_1pr_F_withEGs;
  TH1D* h_TkEG_PtResolution_3pr_F_withEGs;
  TH1D* h_TkEG_PtResolution_noNeutrals_F_noEGs;
  TH1D* h_TkEG_PtResolution_withNeutrals_F_noEGs;
  TH1D* h_TkEG_PtResolution_1pr_F_noEGs;
  TH1D* h_TkEG_PtResolution_3pr_F_noEGs;

  // eta > 0
  TH1D* h_TkEG_PtResolution_noNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_PtResolution_withNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_PtResolution_1pr_F_withEGs_posEta;
  TH1D* h_TkEG_PtResolution_3pr_F_withEGs_posEta;
  TH1D* h_TkEG_PtResolution_noNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_PtResolution_withNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_PtResolution_1pr_F_noEGs_posEta;
  TH1D* h_TkEG_PtResolution_3pr_F_noEGs_posEta;

  // eta < 0
  TH1D* h_TkEG_PtResolution_noNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_PtResolution_withNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_PtResolution_1pr_F_withEGs_negEta;
  TH1D* h_TkEG_PtResolution_3pr_F_withEGs_negEta;
  TH1D* h_TkEG_PtResolution_noNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_PtResolution_withNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_PtResolution_1pr_F_noEGs_negEta;
  TH1D* h_TkEG_PtResolution_3pr_F_noEGs_negEta;

  TH1D* h_TkEG_EtResolution;
  TH1D* h_TkEG_EtResolution_C;
  TH1D* h_TkEG_EtResolution_I;
  TH1D* h_TkEG_EtResolution_F;
  TH1D* h_TkEG_EtResolution_noNeutrals;
  TH1D* h_TkEG_EtResolution_withNeutrals;
  TH1D* h_TkEG_EtResolution_1pr;
  TH1D* h_TkEG_EtResolution_3pr;
  TH1D* h_TkEG_EtResolution_withEGs;
  TH1D* h_TkEG_EtResolution_noEGs;

  TH1D* h_TkEG_EtResolution_withNeutrals_1pr;
  TH1D* h_TkEG_EtResolution_withNeutrals_3pr;
  TH1D* h_TkEG_EtResolution_withNeutrals_1pion0;
  TH1D* h_TkEG_EtResolution_withNeutrals_2pion0;
  TH1D* h_TkEG_EtResolution_withNeutrals_3pion0;
  TH1D* h_TkEG_EtResolution_withNeutrals_4pion0;

  TH1D* h_TkEG_EtResolution_noNeutrals_withEGs;
  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs;
  TH1D* h_TkEG_EtResolution_1pr_withEGs;
  TH1D* h_TkEG_EtResolution_3pr_withEGs;
  TH1D* h_TkEG_EtResolution_noNeutrals_noEGs;
  TH1D* h_TkEG_EtResolution_withNeutrals_noEGs;
  TH1D* h_TkEG_EtResolution_1pr_noEGs;
  TH1D* h_TkEG_EtResolution_3pr_noEGs;

  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs_0to5GeV;
  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs_5to10GeV;
  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs_10to15GeV;
  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs_15to20GeV;
  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs_20to30GeV;
  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs_30to40GeV;
  TH1D* h_TkEG_EtResolution_withNeutrals_withEGs_40to50GeV;

  TH1D* h_TkEG_EtResolution_F_withEGs;
  TH1D* h_TkEG_EtResolution_F_withEGs_posEta;
  TH1D* h_TkEG_EtResolution_F_withEGs_negEta;
  TH1D* h_TkEG_EtResolution_F_noEGs;
  TH1D* h_TkEG_EtResolution_F_noEGs_posEta;
  TH1D* h_TkEG_EtResolution_F_noEGs_negEta;

  TH1D* h_TkEG_EtResolution_noNeutrals_F;
  TH1D* h_TkEG_EtResolution_withNeutrals_F;
  TH1D* h_TkEG_EtResolution_1pr_F;
  TH1D* h_TkEG_EtResolution_3pr_F;

  TH1D* h_TkEG_EtResolution_noNeutrals_F_withEGs;
  TH1D* h_TkEG_EtResolution_withNeutrals_F_withEGs;
  TH1D* h_TkEG_EtResolution_1pr_F_withEGs;
  TH1D* h_TkEG_EtResolution_3pr_F_withEGs;
  TH1D* h_TkEG_EtResolution_noNeutrals_F_noEGs;
  TH1D* h_TkEG_EtResolution_withNeutrals_F_noEGs;
  TH1D* h_TkEG_EtResolution_1pr_F_noEGs;
  TH1D* h_TkEG_EtResolution_3pr_F_noEGs;

  // eta > 0
  TH1D* h_TkEG_EtResolution_noNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_EtResolution_withNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_EtResolution_1pr_F_withEGs_posEta;
  TH1D* h_TkEG_EtResolution_3pr_F_withEGs_posEta;
  TH1D* h_TkEG_EtResolution_noNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_EtResolution_withNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_EtResolution_1pr_F_noEGs_posEta;
  TH1D* h_TkEG_EtResolution_3pr_F_noEGs_posEta;

  // eta < 0
  TH1D* h_TkEG_EtResolution_noNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_EtResolution_withNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_EtResolution_1pr_F_withEGs_negEta;
  TH1D* h_TkEG_EtResolution_3pr_F_withEGs_negEta;
  TH1D* h_TkEG_EtResolution_noNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_EtResolution_withNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_EtResolution_1pr_F_noEGs_negEta;
  TH1D* h_TkEG_EtResolution_3pr_F_noEGs_negEta;
  
  TH1D* h_TkEG_EtaResolution;
  TH1D* h_TkEG_EtaResolution_C;
  TH1D* h_TkEG_EtaResolution_I;
  TH1D* h_TkEG_EtaResolution_F;
  TH1D* h_TkEG_EtaResolution_noNeutrals;
  TH1D* h_TkEG_EtaResolution_withNeutrals;
  TH1D* h_TkEG_EtaResolution_1pr;
  TH1D* h_TkEG_EtaResolution_3pr;

  TH1D* h_TkEG_EtaResolution_F_withEGs;
  TH1D* h_TkEG_EtaResolution_F_withEGs_posEta;
  TH1D* h_TkEG_EtaResolution_F_withEGs_negEta;
  TH1D* h_TkEG_EtaResolution_F_noEGs;
  TH1D* h_TkEG_EtaResolution_F_noEGs_posEta;
  TH1D* h_TkEG_EtaResolution_F_noEGs_negEta;

  TH1D* h_TkEG_EtaResolution_noNeutrals_F;
  TH1D* h_TkEG_EtaResolution_withNeutrals_F;
  TH1D* h_TkEG_EtaResolution_1pr_F;
  TH1D* h_TkEG_EtaResolution_3pr_F;

  TH1D* h_TkEG_EtaResolution_noNeutrals_F_withEGs;
  TH1D* h_TkEG_EtaResolution_withNeutrals_F_withEGs;
  TH1D* h_TkEG_EtaResolution_1pr_F_withEGs;
  TH1D* h_TkEG_EtaResolution_3pr_F_withEGs;
  TH1D* h_TkEG_EtaResolution_noNeutrals_F_noEGs;
  TH1D* h_TkEG_EtaResolution_withNeutrals_F_noEGs;
  TH1D* h_TkEG_EtaResolution_1pr_F_noEGs;
  TH1D* h_TkEG_EtaResolution_3pr_F_noEGs;

  // eta > 0
  TH1D* h_TkEG_EtaResolution_noNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_EtaResolution_withNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_EtaResolution_1pr_F_withEGs_posEta;
  TH1D* h_TkEG_EtaResolution_3pr_F_withEGs_posEta;
  TH1D* h_TkEG_EtaResolution_noNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_EtaResolution_withNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_EtaResolution_1pr_F_noEGs_posEta;
  TH1D* h_TkEG_EtaResolution_3pr_F_noEGs_posEta;

  // eta < 0
  TH1D* h_TkEG_EtaResolution_noNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_EtaResolution_withNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_EtaResolution_1pr_F_withEGs_negEta;
  TH1D* h_TkEG_EtaResolution_3pr_F_withEGs_negEta;
  TH1D* h_TkEG_EtaResolution_noNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_EtaResolution_withNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_EtaResolution_1pr_F_noEGs_negEta;
  TH1D* h_TkEG_EtaResolution_3pr_F_noEGs_negEta;

  TH1D* h_TkEG_PhiResolution;
  TH1D* h_TkEG_PhiResolution_C;
  TH1D* h_TkEG_PhiResolution_I;
  TH1D* h_TkEG_PhiResolution_F;
  TH1D* h_TkEG_PhiResolution_noNeutrals;
  TH1D* h_TkEG_PhiResolution_withNeutrals;
  TH1D* h_TkEG_PhiResolution_1pr;
  TH1D* h_TkEG_PhiResolution_3pr;

  TH1D* h_TkEG_PhiResolution_F_withEGs;
  TH1D* h_TkEG_PhiResolution_F_withEGs_posEta;
  TH1D* h_TkEG_PhiResolution_F_withEGs_negEta;
  TH1D* h_TkEG_PhiResolution_F_noEGs;
  TH1D* h_TkEG_PhiResolution_F_noEGs_posEta;
  TH1D* h_TkEG_PhiResolution_F_noEGs_negEta;

  TH1D* h_TkEG_PhiResolution_noNeutrals_F;
  TH1D* h_TkEG_PhiResolution_withNeutrals_F;
  TH1D* h_TkEG_PhiResolution_1pr_F;
  TH1D* h_TkEG_PhiResolution_3pr_F;

  TH1D* h_TkEG_PhiResolution_noNeutrals_F_withEGs;
  TH1D* h_TkEG_PhiResolution_withNeutrals_F_withEGs;
  TH1D* h_TkEG_PhiResolution_1pr_F_withEGs;
  TH1D* h_TkEG_PhiResolution_3pr_F_withEGs;
  TH1D* h_TkEG_PhiResolution_noNeutrals_F_noEGs;
  TH1D* h_TkEG_PhiResolution_withNeutrals_F_noEGs;
  TH1D* h_TkEG_PhiResolution_1pr_F_noEGs;
  TH1D* h_TkEG_PhiResolution_3pr_F_noEGs;

  // eta > 0
  TH1D* h_TkEG_PhiResolution_noNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_PhiResolution_withNeutrals_F_withEGs_posEta;
  TH1D* h_TkEG_PhiResolution_1pr_F_withEGs_posEta;
  TH1D* h_TkEG_PhiResolution_3pr_F_withEGs_posEta;
  TH1D* h_TkEG_PhiResolution_noNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_PhiResolution_withNeutrals_F_noEGs_posEta;
  TH1D* h_TkEG_PhiResolution_1pr_F_noEGs_posEta;
  TH1D* h_TkEG_PhiResolution_3pr_F_noEGs_posEta;

  // eta < 0
  TH1D* h_TkEG_PhiResolution_noNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_PhiResolution_withNeutrals_F_withEGs_negEta;
  TH1D* h_TkEG_PhiResolution_1pr_F_withEGs_negEta;
  TH1D* h_TkEG_PhiResolution_3pr_F_withEGs_negEta;
  TH1D* h_TkEG_PhiResolution_noNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_PhiResolution_withNeutrals_F_noEGs_negEta;
  TH1D* h_TkEG_PhiResolution_1pr_F_noEGs_negEta;
  TH1D* h_TkEG_PhiResolution_3pr_F_noEGs_negEta;

  // Charged Resolution
  TH1D* h_TkEG_ChargedResolution;
  TH1D* h_TkEG_ChargedResolution_noNeutrals;
  TH1D* h_TkEG_ChargedResolution_withNeutrals;
  TH1D* h_TkEG_ChargedResolution_1pr;
  TH1D* h_TkEG_ChargedResolution_3pr;

  // Neutrals Resolution
  TH1D* h_TkEG_NeutralsResolution;
  TH1D* h_TkEG_NeutralsResolution_noNeutrals;
  TH1D* h_TkEG_NeutralsResolution_withNeutrals;
  TH1D* h_TkEG_NeutralsResolution_1pr;
  TH1D* h_TkEG_NeutralsResolution_3pr;

  // Poor Et resolution candidates 
  TH1D* h_PoorEtResolCand_InvMass;
  TH1D* h_PoorEtResolCand_RelIso;
  TH1D* h_PoorEtResolCand_VtxIso;
  TH1D* h_PoorEtResolCand_CHF;
  TH1D* h_PoorEtResolCand_IsoTracks_N;
  TH1D* h_PoorEtResolCand_dR_EG_Seed;
  // Good Et resolution candidates
  TH1D* h_GoodEtResolCand_InvMass;
  TH1D* h_GoodEtResolCand_RelIso;
  TH1D* h_GoodEtResolCand_VtxIso;
  TH1D* h_GoodEtResolCand_CHF;
  TH1D* h_GoodEtResolCand_IsoTracks_N;
  TH1D* h_GoodEtResolCand_dR_EG_Seed;

  // Poor Neutral Resolution candidates
  TH1D* h_TkEG_PoorNeuResol_NeuMultiplicity;
  TH1D* h_TkEG_PoorNeuResol_dR_Pi0_visTau;
  TH1D* h_TkEG_PoorNeuResol_dEta_Pi0_visTau;
  TH1D* h_TkEG_PoorNeuResol_dPhi_Pi0_visTau;
  TH1D* h_TkEG_PoorNeuResol_Pi0_ET;
  TH1D* h_TkEG_PoorNeuResol_dRmin_Pi0_EG;
  TH1D* h_TkEG_PoorNeuResol_Pi0_closestEG_ET;
  TH2D* h_TkEG_PoorNeuResol_Pi0_closestEG_ET_Vs_dRmin_Pi0_EG;
  TH2D* h_TkEG_PoorNeuResol_Pi0_ET_Vs_closestEG_ET;


  TH1D* h_nonMCmatchedCandidates_decayMode;

  TH1D* h_MCmatch_dR;
  TH1D* h_leadTrk_MCmatch;
  TH1D* h_leadTrk4stubs_MCmatch;

  // Turn-Ons
  TH1D* hMcHadronicTau_VisEt;
  TH1D* hMcHadronicTau_VisEt_1pr;
  TH1D* hMcHadronicTau_VisEt_3pr;
  TH1D* hMcHadronicTau_VisEt_withNeutrals;
  TH1D* hMcHadronicTau_VisEt_noNeutrals;

  TH1D* hTkEG_TurnOn25;
  TH1D* hTkEG_TurnOn25_1pr;
  TH1D* hTkEG_TurnOn25_3pr;
  TH1D* hTkEG_TurnOn25_withNeutrals;
  TH1D* hTkEG_TurnOn25_noNeutrals;

  TH1D* hRelIso_TurnOn25;
  TH1D* hRelIso_TurnOn25_1pr;
  TH1D* hRelIso_TurnOn25_3pr;
  TH1D* hRelIso_TurnOn25_withNeutrals;
  TH1D* hRelIso_TurnOn25_noNeutrals;

  TH1D* hVtxIso_TurnOn25;
  TH1D* hVtxIso_TurnOn25_1pr;
  TH1D* hVtxIso_TurnOn25_3pr;
  TH1D* hVtxIso_TurnOn25_withNeutrals;
  TH1D* hVtxIso_TurnOn25_noNeutrals;
  
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

  TH1D* hTkEG_TurnOn50;
  TH1D* hTkEG_TurnOn50_1pr;
  TH1D* hTkEG_TurnOn50_3pr;
  TH1D* hTkEG_TurnOn50_withNeutrals;
  TH1D* hTkEG_TurnOn50_noNeutrals;

  TH1D* hRelIso_TurnOn50;
  TH1D* hRelIso_TurnOn50_1pr;
  TH1D* hRelIso_TurnOn50_3pr;
  TH1D* hRelIso_TurnOn50_withNeutrals;
  TH1D* hRelIso_TurnOn50_noNeutrals;

  TH1D* hVtxIso_TurnOn50;
  TH1D* hVtxIso_TurnOn50_1pr;
  TH1D* hVtxIso_TurnOn50_3pr;
  TH1D* hVtxIso_TurnOn50_withNeutrals;
  TH1D* hVtxIso_TurnOn50_noNeutrals;

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

  // Turn-ons for the best performing WP
  TH1D* hL1Taus_TurnOn25;
  TH1D* hL1Taus_TurnOn25_1pr;
  TH1D* hL1Taus_TurnOn25_3pr;
  TH1D* hL1Taus_TurnOn25_withNeutrals;
  TH1D* hL1Taus_TurnOn25_noNeutrals;
  TH1D* hL1Taus_TurnOn50;
  TH1D* hL1Taus_TurnOn50_1pr;
  TH1D* hL1Taus_TurnOn50_3pr;
  TH1D* hL1Taus_TurnOn50_withNeutrals;
  TH1D* hL1Taus_TurnOn50_noNeutrals;
  
  // SingleTau: Rates
  TH1D* hTkEG_Rate;
  TH1D* hTkEG_Rate_C;
  TH1D* hTkEG_Rate_I;
  TH1D* hTkEG_Rate_F;
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
  
  TH1D* hL1Taus_SingleTau_Rate;

  // DiTau: Rates
  TH1D* hDiTau_Rate_TkEG;
  TH1D* hDiTau_Rate_TkEG_C;
  TH1D* hDiTau_Rate_TkEG_I;
  TH1D* hDiTau_Rate_TkEG_F;
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

  TH1D* hL1Taus_DiTau_Rate;

  // SingleTau: Efficiencies
  TH1D* hTkEG_Eff;
  TH1D* hTkEG_Eff_C;
  TH1D* hTkEG_Eff_I;
  TH1D* hTkEG_Eff_F;
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
  
  TH1D* hL1Taus_SingleTau_Eff;

  // DiTau: Efficiencies
  TH1D* hDiTau_Eff_TkEG;
  TH1D* hDiTau_Eff_TkEG_C;
  TH1D* hDiTau_Eff_TkEG_I;
  TH1D* hDiTau_Eff_TkEG_F;
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
  
  TH1D* hL1Taus_DiTau_Eff;

};

#endif

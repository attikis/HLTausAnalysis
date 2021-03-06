#ifndef Sanity_h
#define Sanity_h

// System 
#include <iostream>
#include <stdlib.h> 

// User
#include "../Framework/src/TreeAnalyserMC.C"

#include "../Auxiliary/src/AuxTools.C"
#include "../Auxiliary/src/Table.C"
#include "../Auxiliary/src/MCTools.C"
//#include "../Auxiliary/src/HistoTools.C"
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
#include "../Plugins/src/L1PixelTrackFit.C"

// ROOT
#include "TEfficiency.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"


using namespace std;

class Sanity : public TreeAnalyserMC{
 public:
  // Constructors/Destructors
  ~Sanity(){};
 Sanity(const string SamplePath,
	  const string SampleName,
	  const string text_, 
	  const int maxEvents_ = -1, 
	  TTree* tree=0) : 
  
  TreeAnalyserMC("", SamplePath, SampleName, text_, maxEvents_, tree) 
    { 
      auxTools_.StopwatchStart();
      mcSample = SampleName;
    };
  
  // Public Variables
  virtual void Loop();

  void PrintSettings(void);

  // Public Variable
  string mcSample;
  bool cfg_AddGenP;   // Flag to enable sanity of GenParticles
  bool cfg_AddL1Tks;  // Flag to enable sanity of TTTracks
  bool cfg_AddStubs;  // Flag to enable sanity all TTSTubs
  bool cfg_AddTPs;    // Flag to enable sanity of all Tracking Particles
  bool cfg_AddEGs;    // Flag to enable sanity of all L1 EG objects
  bool cfg_AddTaus;   // Flag to enable sanity of all L1 Taus
  bool cfg_AddJets;   // Flag to enable sanity of all L1 Jets
  bool cfg_AddMuons;  // Flag to enable sanity of all L1 Muons
  bool cfg_AddSums;   // Flag to enable sanity of all L1 Et sums
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
  
  struct SortAscendingAbs{ bool operator() (double a, double b) const { return abs(a) > abs(b); } }; 
  struct SortDescendingAbs{ bool operator() (double a, double b) const { return abs(a) < abs(b); } };
  struct PtComparatorTP{ bool operator() (TrackingParticle a, TrackingParticle b) const { return a.getPt() > b.getPt(); } };
  struct PtComparatorTTTrack{ bool operator() (TTTrack a, TTTrack b) const { return a.getPt() > b.getPt(); } };
  
 private:
  // Private Functions
  void BookHistos_(void);
  
  void InitVars_(void);  
  
  // Private variables
  AuxTools auxTools_;
  MCTools mcTools_;
  Datasets datasets_;
  HistoTools histoTools_;

  // Histograms
  TH1D* hCounters;
  TH1D* hHepMCEvt_VtxZ;
  TH2D* hHepMCEvt_VtxX_VtxY;

    
  // GenParticles
  TH1D* hGenP_Index;
  TH1D* hGenP_Pt;
  TH1D* hGenP_Eta;
  TH1D* hGenP_Phi;
  TH1D* hGenP_PtVis;
  TH1D* hGenP_EtaVis;
  TH1D* hGenP_PhiVis;
  TH1D* hGenP_Mass;
  TH1D* hGenP_PdgId;
  TH1D* hGenP_Charge;
  TH1D* hGenP_Status;
  TH1D* hGenP_VertexX;
  TH1D* hGenP_VertexY;
  TH1D* hGenP_VertexZ;
  TH1D* hGenP_Mothers;
  TH1D* hGenP_Daughters;
  TH1D* hGenP_FinalDaughters;
  TH1D* hGenP_FinalDaughtersCharged;
  TH1D* hGenP_FinalDaughtersNeutral;

  // TrackingParticles
  TH1D* hTP_Index;
  TH1D* hTP_Pt;
  TH1D* hTP_Eta;
  TH1D* hTP_Phi;
  TH1D* hTP_PdgId;
  TH1D* hTP_X0;
  TH1D* hTP_Y0;
  TH1D* hTP_Z0;
  TH1D* hTP_Charge;
  TH1D* hTP_Dxy;
  TH1D* hTP_D0propagated;
  TH1D* hTP_Z0propagated;
  TH1D* hTP_D0;
  TH1D* hTP_D0Sign;
  TH1D* hTP_D0Phi;
  TH1D* hTP_EventId;
  TH1D* hTP_NMatch;
  TH1D* hTP_TTTrackIndex;
  TH1D* hTP_TTClusters;
  TH1D* hTP_TTStubs;
  TH1D* hTP_TTTracks;

  // TTTracks
  TH1D* hL1Tks_Pt;
  TH1D* hL1Tks_Eta;
  TH1D* hL1Tks_Phi;
  TH1D* hL1Tks_Charge;
  TH1D* hL1Tks_POCAx;
  TH1D* hL1Tks_POCAy;
  TH1D* hL1Tks_POCAz;
  TH1D* hL1Tks_ChiSq;
  TH1D* hL1Tks_ChiSqRed;
  TH1D* hL1Tks_StubPtConsistency;
  TH1D* hL1Tks_RInv;
  TH1D* hL1Tks_IsGenuine;
  TH1D* hL1Tks_IsUnknown;
  TH1D* hL1Tks_IsCombinatoric;
  TH1D* hL1Tks_IsLoose;
  TH1D* hL1Tks_IsFake;
  TH1D* hL1Tks_NStubs;
  TH1D* hL1Tks_NStubsPS;
  TH1D* hL1Tks_NStubsBarrel;
  TH1D* hL1Tks_NStubsEndcap;
  TH1D* hL1Tks_TP_Index;

  // L1Taus
  TH1D* hL1Tau_Index;
  TH1D* hL1Tau_nTaus;
  TH1D* hL1Tau_Et;
  TH1D* hL1Tau_Eta;
  TH1D* hL1Tau_Phi;
  TH1D* hL1Tau_IET;
  TH1D* hL1Tau_IEta;
  TH1D* hL1Tau_IPhi;
  TH1D* hL1Tau_Iso;
  TH1D* hL1Tau_Bx;
  TH1D* hL1Tau_TowerIPhi;
  TH1D* hL1Tau_TowerIEta;
  TH1D* hL1Tau_RawEt;
  TH1D* hL1Tau_IsoEt;
  TH1D* hL1Tau_NTT;
  TH1D* hL1Tau_HasEM;
  TH1D* hL1Tau_IsMerged;
  TH1D* hL1Tau_HwQual;

  // L1Jets
  TH1D* hL1Jet_Index;
  TH1D* hL1Jet_nJets;
  TH1D* hL1Jet_Et;
  TH1D* hL1Jet_Eta;
  TH1D* hL1Jet_Phi;
  TH1D* hL1Jet_IET;
  TH1D* hL1Jet_IEta;
  TH1D* hL1Jet_IPhi;
  TH1D* hL1Jet_Bx;
  TH1D* hL1Jet_RawEt;
  TH1D* hL1Jet_SeedEt;
  TH1D* hL1Jet_TowerIPhi;
  TH1D* hL1Jet_TowerIEta;
  TH1D* hL1Jet_PUEt;
  TH1D* hL1Jet_PUDonutEt0;
  TH1D* hL1Jet_PUDonutEt1;
  TH1D* hL1Jet_PUDonutEt2;
  TH1D* hL1Jet_PUDonutEt3;

  // L1EGs
  TH1D* hL1EG_Index;
  TH1D* hL1EG_nEGs;
  TH1D* hL1EG_Et;
  TH1D* hL1EG_Eta;
  TH1D* hL1EG_Phi;
  TH1D* hL1EG_IEt;
  TH1D* hL1EG_IEta;
  TH1D* hL1EG_IPhi;
  TH1D* hL1EG_Iso;
  TH1D* hL1EG_Bx;
  TH1D* hL1EG_TowerIPhi;
  TH1D* hL1EG_TowerIEta;
  TH1D* hL1EG_RawEt;
  TH1D* hL1EG_IsoEt;
  TH1D* hL1EG_FootprintEt;
  TH1D* hL1EG_NTT;
  TH1D* hL1EG_Shape;
  TH1D* hL1EG_TowerHoE;

  // L1Sums
  TH1D* hL1Sum_Index;
  TH1D* hL1Sum_nSums;
  TH1D* hL1Sum_Et;
  TH1D* hL1Sum_Phi;
  TH1D* hL1Sum_IEt;
  TH1D* hL1Sum_IPhi;
  TH1D* hL1Sum_Bx;
  TH1D* hL1Sum_Type;
  
};

#endif

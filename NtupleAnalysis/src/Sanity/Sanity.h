#ifndef Sanity_h
#define Sanity_h

// System 
#include <iostream>
#include <stdlib.h> 
#include <iomanip>

// User
#include "../Framework/src/TreeAnalyserMC.C"

#include "../Auxiliary/src/AuxTools.C"
#include "../Auxiliary/src/Table.C"
#include "../Auxiliary/src/MCTools.C"
#include "../Auxiliary/src/HistoTools.C"
#include "../Auxiliary/src/Datasets.C" 
#include "../DataFormat/src/GenParticle.C"
#include "../DataFormat/src/TrackingParticle.C"
#include "../DataFormat/interface/TTTrack.h"

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

  void PrintGenParticleCollection(vector<GenParticle> collection);

  void PrintTrackingParticleCollection(vector<TrackingParticle> collection);

  void PrintTTTrackCollection(vector<TTTrack> collection);

  
  // Public Variables
  bool DEBUG;
  string mcSample;
  string tk_Collection;
  int tk_nFitParams;
  double tk_minPt;
  double tk_minEta;
  double tk_maxEta;
  double tk_maxChiSqRed;
  double tk_minStubs;
  double tk_minStubsPS;
  double tk_maxStubsPS;

  struct SortAscendingAbs{ bool operator() (double a, double b) const { return abs(a) > abs(b); } }; 
  struct SortDescendingAbs{ bool operator() (double a, double b) const { return abs(a) < abs(b); } };
  struct PtComparatorTP{ bool operator() (TrackingParticle a, TrackingParticle b) const { return a.getPt() > b.getPt(); } };
  struct PtComparatorTTTrack{ bool operator() (TTTrack a, TTTrack b) const { return a.getPt() > b.getPt(); } };
  
 private:
  // Private Functions
  void BookHistos_(void);
  
  void InitVars_(void);

  TrackingParticle GetTrackingParticle(unsigned int Index);

  vector<TrackingParticle> GetTrackingParticles(bool bPrintList=false);
 
  TTTrack GetTTTrack(unsigned int Index, const unsigned int nFitParams = 5);

  vector<TTTrack> GetTTTracks(const double minPt = 0.0,
			      const double minEta = 0.0,
			      const double maxEta = 9999.9,
			      const double maxChiSqRed = 9999.9,
			      const unsigned int minStubs = 0,
			      const unsigned int minStubsPS = 0,
			      const unsigned int maxStubsPS = 999,
			      const unsigned nFitParams = 5,
			      bool bPrintList = false);


  void SetGenParticleMomsAndDaus(GenParticle &p);

  void SetGenParticleFinalDaughters(GenParticle &p);

  GenParticle GetGenParticle(unsigned int Index);

  vector<GenParticle> GetGenParticles(bool bPrintList=false);

  vector<GenParticle> GetGenParticles(int pdgId, bool isLastCopy=false);

  void GetHadronicTauFinalDaughters(GenParticle hTau, vector<unsigned short> &Daug);

  
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


};

#endif
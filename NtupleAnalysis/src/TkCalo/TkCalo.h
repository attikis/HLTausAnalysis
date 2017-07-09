#ifndef TkCalo_h
#define TkCalo_h

// System 
#include <iostream>
#include <stdlib.h> 

// User
#include "../Framework/src/TreeAnalyserMC.C"
#include "../Auxiliary/src/AuxTools.C"
#include "../DataFormat/interface/TTTrack.h"
#include "../DataFormat/src/L1EG.C"
#include "../DataFormat/src//L1TkEGParticle.C"
#include "../DataFormat/src/TrackingParticle.C" // for GetTTTrack function

//#include "../Auxiliary/src/Table.C"
//#include "../Auxiliary/src/MCTools.C"
#include "../Auxiliary/src/HistoTools.C"
//#include "../Auxiliary/src/Datasets.C" 
//#include "../DataFormat/src/L1TkTauParticle.C"
//#include "../DataFormat/src/GenParticle.C"
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
class TkCalo : public TreeAnalyserMC{
 public:
  // Destructor
  ~TkCalo(){};
  // Constructor
  TkCalo(const string SamplePath,
		const string SampleName,
		const string text_, 
		const int maxEvents_ = -1, 
		TTree* tree=0) :  
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
  bool cfg_AddL1Tks;  // Flag to enable sanity of TTTracks
  bool cfg_AddEGs;    // Flag to enable sanity of all L1 EG objects
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
  void InitObjects(void);
  void InitVars_(void);
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
  TrackingParticle GetTrackingParticle(unsigned int Index);
  L1EG GetL1EG(unsigned int Index);  
  vector<L1EG> GetL1EGs(bool bPrintList=false);
  float DeltaPhi(float phi1, float phi2);
  float deltaR(float eta1, float eta2, float phi1, float phi2);
  //struct PtComparatorTTTrack{ bool operator() (TTTrack a, TTTrack b) const { return a.getPt() > b.getPt(); } };
  //struct EtComparatorL1EG{ bool operator() (L1EG a, L1EG b) const { return a.et() > b.et(); } };
  
  // Private variables
  AuxTools auxTools_;
  Datasets datasets_;
  HistoTools histoTools_;
  
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
  float minDeltaR_iso;  
  float maxDeltaR_iso;   
  float maxDeltaZ_iso;   
  float useRelIso;  
  float maxRelIso;  
  
  vector<TTTrack> TTTracks;
  vector<L1EG> L1EGs;
  vector< vector <TTTrack> > trackTauCandidates;
  vector<L1TkEGParticle> TauCandidates;
  vector<L1TkEGParticle> TauCandidatesIsolated;
  
  // Histogram declarations
  TH1D* h_leadTrks_Multiplicity;
  
  TH1D* h_leadTrks_Pt;
  TH1D* h_leadTrks_Eta;
  TH1D* h_leadTrks_Phi;
  TH2D* h_leadTrks_Phi_Eta;

  TH1D* h_leadTrkSelection;

  TH1D* h_clustTrks_Pt;
  TH1D* h_clustTrks_Eta;
  TH1D* h_clustTrks_Phi;
  TH2D* h_clustTrks_Phi_Eta;

  TH1D* h_trkClusters_MultiplicityPerCluster;
  TH1D* h_trkClusters_Pt;
  TH1D* h_trkClusters_M;

  TH1D* h_clustEGs_MultiplicityPerCluster;
  TH1D* h_clustEGs_Et;
  TH1D* h_clustEGs_Eta;
  TH1D* h_clustEGs_Phi;
  TH2D* h_clustEGs_Phi_Eta;

  TH1D* h_EGClusters_Pt;
  TH1D* h_EGClusters_M;

  TH1D* h_trkClusters_relIso;
};

#endif

#ifndef TreeReaderMC_h
#define TreeReaderMC_h

// System
#include <vector>

// User
#include "TreeDefinitionGenP.h"
#include "../src/FileOpener.C"
#include "../src/TreeReaderReco.C"
#include "../../Auxiliary/src/AuxTools.C"
#include "../../Auxiliary/src/MCTools.C"
#include "../../Auxiliary/src/HistoTools.C"
#include "../../Auxiliary/src/Datasets.C" 
#include "../../DataFormat/src/L1TkTauParticle.C"
#include "../../DataFormat/src/GenParticle.C"
#include "../../DataFormat/src/TrackingParticle.C"
#include "../../DataFormat/interface/TTTrack.h"
#include "../../DataFormat/interface/TTPixelTrack.h"
#include "../../DataFormat/src/L1EG.C"
#include "../../DataFormat/src/L1Jet.C"
#include "../../DataFormat/src/L1Tau.C"
#include "../../DataFormat/src/L1Sum.C"
#include "../../Plugins/src/L1PixelTrackFit.C"

// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TEfficiency.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

// System 
#include <iostream>
#include <stdlib.h> 

class TreeReaderMC : public TreeReaderReco, public virtual TREEDEFINITIONGENP
{
  public:
  TreeReaderMC() { 
    std::cout << " TreeReaderMC : The constructor needs some arguments \n"; 
  };
  TreeReaderMC(const std::string SamplePath, const std::string SampleName, TChain *chain=0) :
  TreeReaderReco(SamplePath, SampleName, chain) { InitGenP(fChain); };

  // Public functions
  void PrintGenParticleCollection(vector<GenParticle> collection);

  void PrintTrackingParticleCollection(vector<TrackingParticle> collection);

  void PrintTTTrackCollection(vector<TTTrack> collection);

  void PrintL1SumCollection(vector<L1Sum> collection);

  void PrintL1EGCollection(vector<L1EG> collection);

  void PrintL1JetCollection(vector<L1Jet> collection);

  void PrintL1TauCollection(vector<L1Tau> collection);

  void PrintL1TkTauParticleCollection(vector<L1TkTauParticle> collection);

  // Public Variables
  bool cfg_DEBUG;
  
  struct SortAscendingAbs{ bool operator() (double a, double b) const { return abs(a) > abs(b); } }; 
  struct SortDescendingAbs{ bool operator() (double a, double b) const { return abs(a) < abs(b); } };
  struct PtComparatorTP{ bool operator() (TrackingParticle a, TrackingParticle b) const { return a.getPt() > b.getPt(); } };
  struct PtComparatorTTTrack{ bool operator() (TTTrack a, TTTrack b) const { return a.getPt() > b.getPt(); } };
  
  protected:
  // Protected Functions
  void InitVars_(void);

  TrackingParticle GetTrackingParticle(unsigned int Index);

  vector<TrackingParticle> GetTrackingParticles(bool bPrintList=false);
 
  TTTrack GetTTTrack(unsigned int Index, const unsigned int nFitParams = 5);

  vector<TTTrack> GetTTTracks(const double minPt = 0.0,
			      const double minEta = 0.0,
			      const double maxEta = 9999.9,
			      const double maxChiSqRed = 9999.9,
			      const unsigned int minStubs = 0,
			      const unsigned nFitParams = 5,
			      bool bPrintList = false);

  void SetGenParticleMomsAndDads(GenParticle &p);

  void SetGenParticleFinalDaughters(GenParticle &p);

  GenParticle GetGenParticle(unsigned int Index);

  vector<GenParticle> GetGenParticles(bool bPrintList=false);

  vector<GenParticle> GetGenParticles(int pdgId, bool isLastCopy=false);

  L1EG GetL1EG(unsigned int Index);
  
  vector<L1EG> GetL1EGs(bool bPrintList=false);
 
  void GetL1Muon(unsigned int Index);
  
  void GetL1Muons(bool bPrintList=false);

  L1Jet GetL1Jet(unsigned int Index);

  vector<L1Jet> GetL1Jets(bool bPrintList=false);
  
  L1Tau GetL1Tau(unsigned int Index);

  vector<L1Tau> GetL1Taus(bool bPrintList=false);
  
  L1Sum GetL1Sum(unsigned int Index);

  vector<L1Sum> GetL1Sums(bool bPrintList=false);
  
  vector<GenParticle> GetHadronicGenTaus(vector<GenParticle> GenTaus,
					 double visEt=20.0,
					 double visEta=2.3);
  
  void GetHadronicTauFinalDaughters(GenParticle hTau, vector<unsigned short> &Daug);
  
  // Protected variables
  AuxTools auxTools_;
  MCTools mcTools_;
  Datasets datasets_;
  HistoTools histoTools_;
};

#endif

#ifdef TreeReaderMC_cxx

// Nothing goes here

#endif // #ifdef TreeReaderMC_cxx

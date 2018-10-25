#ifndef L1TkEGParticle_h
#define L1TkEGParticle_h

// System
#include <iostream>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../DataFormat/interface/TTTrack.h"
#include "../../DataFormat/src/L1EG.C"
#include "../../DataFormat/src/EG.C"
#include "../../DataFormat/interface/GenParticle.h"
//#include "../../DataFormat/interface/TTPixelTrack.h"

using namespace std;

class L1TkEGParticle{     
 public:
  // Constructors
  L1TkEGParticle();
  L1TkEGParticle(vector<TTTrack> tracks,
		 vector<EG> EGs,
		 GenParticle genTau,
		 bool match);

  L1TkEGParticle(double vtxIso,
		 double relIso,
		 double CHF,
		 double NHF, 
		 int isoTracks_N);

  // Destructor
  ~L1TkEGParticle() {};
  
  // Function declarations
  void InitVars_(void);
  void PrintTTTracks();
  void AddTrack(TTTrack trk) { theTracks.push_back(trk); }
  TTTrack GetLeadingTrack() const { return theTracks[0]; } 
  vector<TTTrack> GetTracks() const {return theTracks; }
  vector<EG> GetEGs() const{ return theEGs;}
  bool HasMatchingGenParticle(void) const{return theMatching;}
  
  GenParticle GetMatchingGenParticle() const {return theGenTau;}
  TLorentzVector GetTotalP4();

  double GetVtxIso()  const { return vtxIso_;}
  double GetRelIso()  const { return relIso_;}
  void SetVtxIso(double vtxIso) { vtxIso_ = vtxIso;}
  void SetRelIso(double relIso) { relIso_ = relIso;}
  double GetCHF() const { return CHF_;}
  void SetCHF(double CHF) { CHF_ = CHF;}
  double GetNHF() const { return NHF_;}
  void SetNHF(double NHF) { NHF_ = NHF;}
  int GetIsoTracksN() const { return isoTracks_N_;}
  void SetIsoTracksN(int isoTracks_N) { isoTracks_N_ = isoTracks_N;}

  double vtxIso_;
  double relIso_;
  double CHF_;
  double NHF_;
  double isoTracks_N_;

  double GetTrackBasedPt();  
  double GetTotalPt();
  double GetTrackInvMass();
  double GetEGInvMass();
  double GetGenTauPt();

  double GetTrackBasedEt();
  double GetEGBasedEt(); 
  double GetTotalEt();
  double GetGenTauEt();  
  
 private:
  AuxTools auxTools;
  vector<TTTrack> theTracks;
  vector<EG> theEGs;
  GenParticle theGenTau;
  bool theMatching;
};

#endif

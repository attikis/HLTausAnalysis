#ifndef L1TkEMParticle_h
#define L1TkEMParticle_h

// System
#include <iostream>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../DataFormat/interface/TTTrack.h"
#include "../../DataFormat/src/L1CaloTP.C"
#include "../../DataFormat/interface/GenParticle.h"
//#include "../../DataFormat/interface/TTPixelTrack.h"

using namespace std;

class L1TkEMParticle{     
 public:
  // Constructors
  L1TkEMParticle();
  L1TkEMParticle(vector<TTTrack> tracks, vector<L1CaloTP> EcalTPs, GenParticle genTau, bool match);
  // Destructor
  ~L1TkEMParticle() {};
  
  // Function declarations
  void InitVars_(void);
  void PrintTTTracks();
  void AddTrack(TTTrack trk) { theTracks.push_back(trk); }
  TTTrack GetLeadingTrack() const { return theTracks[0]; } 
  bool HasMatchingGenParticle(void) const{return theMatching;}
  
  double GetTrackBasedPt();  
  double GetTrackInvMass();
  double GetEMInvMass();
  double GetGenTauPt();

  double GetTrackBasedEt();
  double GetEMBasedEt(); 
  double GetGenTauEt();  
  
 private:
  AuxTools auxTools;
  vector<TTTrack> theTracks;
  vector<L1CaloTP> theEcalTPs;
  GenParticle theGenTau;
  bool theMatching;
};

#endif

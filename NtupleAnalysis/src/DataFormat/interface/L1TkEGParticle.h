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
  L1TkEGParticle(vector<TTTrack> tracks, vector<EG> EGs, GenParticle genTau, bool match);
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
  
  double GetTrackBasedPt();  
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

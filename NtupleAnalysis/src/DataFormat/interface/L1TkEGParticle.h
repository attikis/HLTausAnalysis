#ifndef L1TkEGParticle_h
#define L1TkEGParticle_h

// System
#include <iostream>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../DataFormat/interface/TTTrack.h"
#include "../../DataFormat/src/L1EG.C"
//#include "../../DataFormat/interface/TTPixelTrack.h"
//#include "../../DataFormat/interface/GenParticle.h"

using namespace std;

class L1TkEGParticle{     
 public:
  // Constructors
  L1TkEGParticle();
  L1TkEGParticle(vector<TTTrack> tracks, vector<L1EG> EGs);
  // Destructor
  ~L1TkEGParticle() {};
  
  // Function declarations
  void InitVars_(void);
  void PrintTTTracks();
  void AddTrack(TTTrack trk) { theTracks.push_back(trk); }
  TTTrack GetLeadingTrack() const { return theTracks[0]; } 

  double GetTrackPtSum();  
  double GetTrackInvMass();
  double GetEGInvMass();
  
 private:
  AuxTools auxTools;
  vector<TTTrack> theTracks;
  vector<L1EG> theEGs;
};

#endif

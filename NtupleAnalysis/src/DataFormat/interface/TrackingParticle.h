#ifndef TrackingParticle_h
#define TrackingParticle_h

// System
#include <iostream>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../Auxiliary/src/Table.C"
#include "../../Auxiliary/interface/constants.h"
#include "../src/TTTrack.C"

// ROOT
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

using namespace std;

class TrackingParticle{     
 public:
  // Constructors/Destructors
  TrackingParticle();
  TrackingParticle(unsigned short index,
		   TVector3 momentum,
		   ROOT::Math::XYZVector poca, // ROOT::Math::XYZVector (https://root.cern.ch/doc/master/Vector3DPage.html)
		   int charge,
		   unsigned short pdgId,
		   unsigned short nMatch,
		   unsigned short ttTrackIndex,
		   unsigned short ttClusters,
		   unsigned short ttStubs,
		   unsigned short ttTracks);
  ~TrackingParticle();
  
  unsigned short index(void) const {return theIndex;} 
  TLorentzVector p4(double mass=pionMass);
  TVector3 p3(void) const {return theMomentum;}
  TVector3 getMomentum(void) const {return theMomentum;}
  ROOT::Math::XYZVector getPOCA(void) const {return thePOCA;} //https://root.cern.ch/doc/master/namespaceROOT_1_1Math.html#a676e5bc512b53b6999fd0200df1e81ab
  double getPt(void) const {return theMomentum.Perp(); }
  double getEta(void) const {return theMomentum.Eta(); }
  double getPhi(void) const {return theMomentum.Phi(); }
  double getX0(void) const {return thePOCA.X();}
  double getY0(void) const {return thePOCA.Y();}
  double getZ0(void) const {return thePOCA.Z();}
  int getCharge(void) const {return theCharge;}
  string getQ(void) const {return theQ;}
  double getD0(void)const {return theD0;}
  double getD0Sign(void)const {return theD0Sign;}
  double getD0Phi(void)const {return theD0Phi;}
  unsigned short getPdgId(void) const {return thePdgId;}
  unsigned short getNMatch(void) const {return theNMatch;}
  unsigned short getTTTrackIndex(void) const {return theTTTrackIndex;}
  unsigned short getTTClusters(void) const {return theTTClusters;}
  unsigned short getTTStubs(void) const {return theTTStubs;}
  unsigned short getTTTracks(void) const {return theTTTracks;}
  void SetTTTrack(TTTrack track) {theTTTrack = track;}
  void PrintProperties(void);
  void PrintAllProperties(void);
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Functions
  string _getQ(void);
  double _getD0(void);
  double _getD0Sign(void);
  double _getD0Phi(void);
  
  // Variable declaration
  unsigned short theIndex;
  TVector3 theMomentum;
  ROOT::Math::XYZVector thePOCA;
  unsigned short thePdgId;
  int theCharge;
  string theQ;
  double theD0;
  double theD0Sign;
  double theD0Phi;
  unsigned int theNMatch;
  unsigned int theTTTrackIndex;
  unsigned int theTTClusters;
  unsigned int theTTStubs;
  unsigned int theTTTracks;
  TTTrack theTTTrack;
  
};

#endif

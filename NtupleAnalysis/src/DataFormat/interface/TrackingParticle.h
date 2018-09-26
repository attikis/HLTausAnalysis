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
		   float d0_propagated,
		   float z0_propagated,
		   float d0_produced, 
		   float z0_produced, 
		   float dxy,
		   int charge,
		   int pdgId,
		   int nMatch,
		   int nStubs,
		   float ttTrackPt, 
		   float ttTrackEta, 
		   float ttTrackPhi, 
		   float ttTrackZ0, 
		   float ttTrackD0, 
		   float ttTrackChi2, 
		   int ttTrackNstubs, 
  		   int eventId);
  ~TrackingParticle();
  
  unsigned short index(void) const {return theIndex;} 
  TLorentzVector p4(float mass=pionMass);
  TVector3 p3(void) const {return theMomentum;}
  TVector3 getMomentum(void) const {return theMomentum;}
  float getPt(void) const {return theMomentum.Perp(); }
  float getEta(void) const {return theMomentum.Eta(); }
  float getPhi(void) const {return theMomentum.Phi(); }
  int getCharge(void) const {return theCharge;}
  string getQ(void) const {return theQ;}
  float getDxy(void)const {return theDxy;}
  float getD0propagated(void)const {return theD0propagated;}
  float getZ0propagated(void)const {return theZ0propagated;}
  float getD0produced(void)const {return theD0produced;}
  float getZ0produced(void)const {return theZ0produced;}
  int getEventId(void)const {return theEventId;}
  int getPdgId(void) const {return thePdgId;}
  int getNMatch(void) const {return theNMatch;}
  int getNStubs(void) const {return theNStubs;}
  float getTTTrackPt(void) const {return theTTTrackPt;}
  float getTTTrackEta(void) const {return theTTTrackEta;}
  float getTTTrackPhi(void) const {return theTTTrackPhi;}
  float getTTTrackZ0(void) const {return theTTTrackZ0;}
  float getTTTrackD0(void) const {return theTTTrackD0;}
  float getTTTrackChi2(void) const {return theTTTrackChi2;}
  float getTTTrackNStubs(void) const {return theTTTrackNStubs;}
  void SetTTTrack(TTTrack track) {theTTTrack = track;}
  void PrintProperties(bool bPrintTitleRow=true);
  void PrintAllProperties(void);
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Functions
  string _getQ(void);
  //float _getD0Sign(void);
  //float _getD0Phi(void);
  
  // Variable declaration
  unsigned short theIndex;
  TVector3 theMomentum;
  int thePdgId;
  int theCharge;
  float theD0propagated;
  float theZ0propagated;
  string theQ;
  float theDxy;
  float theD0produced;
  float theZ0produced;
  float theD0Sign;
  float theD0Phi;
  int theNMatch;
  int theNStubs;
  float theTTTrackPt;
  float theTTTrackEta;
  float theTTTrackPhi;
  float theTTTrackZ0;
  float theTTTrackD0;
  float theTTTrackChi2;
  int theTTTrackNStubs;
  int theEventId;
  TTTrack theTTTrack;
  
};

#endif

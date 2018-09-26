#ifndef TTTrack_h
#define TTTrack_h

// System
#include <iostream>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../Auxiliary/src/Table.C"
// #include "../src/TrackingParticle.C"  cyclic/fwd-declaration errors

// ROOT
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

using namespace std;

class TrackingParticle;

class TTTrack{     
 public:
  // Constructors/Destructors
  TTTrack();
  TTTrack(unsigned short index,
	  TVector3 aMomentum,
	  float d0, 
	  float z0,
	  float aChi2,	    
	  bool isGenuine,
	  bool isUnknown,
	  bool isCombinatoric,
	  bool isLoose,
	  bool isFake,
	  unsigned int nStubs,
	  int matchTP_pdgId,
	  float matchTP_pt,
	  float matchTP_eta,
	  float matchTP_phi,
	  float matchTP_z0,
	  float matchTP_dxy,
	  unsigned int nFitParams=5);
  ~TTTrack();
  
  unsigned short index(void) const {return theIndex;}
  TVector3 p3(void) const {return theMomentum;}
  TLorentzVector p4(float mass=pionMass);
  TVector3 getMomentum(void) const {return theMomentum;}
  float getPt(void) const {return theMomentum.Perp(); }
  float getEta(void) const {return theMomentum.Eta(); }
  float getPhi(void) const {return theMomentum.Phi(); }
  float getZ0(void) const {return theZ0;}
  float getChi2(void) const {return theChi2;}
  float getChi2Red(void) const {return theChi2Red;}
  bool getIsGenuine(void) const {return theIsGenuine;}
  bool getIsUnknown(void) const {return theIsUnknown;}
  bool getIsCombinatoric(void) const {return theIsCombinatoric;}
  bool getIsLoose(void) const {return theIsLoose;}
  bool getIsFake(void) const {return theIsFake;}
  unsigned int getNumOfStubs(void) const {return theNStubs;}
  int getTPPdgId(void) const {return theTPPdgId;}
  float getTPPt(void) const {return theTPPt;}
  float getTPEta(void) const {return theTPEta;}
  float getTPPhi(void) const {return theTPPhi;}
  float getTPZ0(void) const {return theTPZ0;}
  float getTPdxy(void) const {return theTPdxy;}
  float getD0(void) const {return theD0;}
  int getDOF(void);
  void PrintProperties(void);
  void PrintAllProperties(void);

  bool operator<(const TTTrack& other) const {return this->getPt() > other.getPt();}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  unsigned short theIndex;
  TVector3 theMomentum;
  float theZ0;
  float theD0;
  float theChi2;
  float theChi2Red;
  bool theIsGenuine;
  bool theIsUnknown;
  bool theIsCombinatoric;
  bool theIsLoose;
  bool theIsFake;
  unsigned int theNStubs;
  int theTPPdgId;
  float theTPPt;
  float theTPEta;
  float theTPPhi;
  float theTPZ0;
  float theTPdxy;
  int theFitParameters;
  int theDOF;
  unsigned int theStubs;
  
};

#endif

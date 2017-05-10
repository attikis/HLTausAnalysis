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
	  ROOT::Math::XYZVector aPOCA, // https://root.cern.ch/doc/master/Vector3DPage.html
	  double aRInv,
	  double aChi2,	    
	  double aStubPtConsistency,
	  bool isGenuine,
	  bool isUnknown,
	  bool isCombinatoric,
	  bool isLoose,
	  bool isFake,
	  unsigned int nStubs,
	  unsigned int nStubsPS,
	  unsigned int nStubsBarrel,
	  unsigned int nStubsEndcap,
	  int matchTP_index,
	  vector<unsigned int> stubs_isPS,
	  vector<unsigned int> stubs_iDisk,
	  vector<unsigned int> stubs_iLayer,
	  vector<unsigned int> stubs_iPhi,
	  vector<unsigned int> stubs_iRing,
	  vector<unsigned int> stubs_iSide,
	  vector<unsigned int> stubs_iZ,
	  unsigned int nFitParams=5);
  ~TTTrack();
  
  unsigned short index(void) const {return theIndex;}
  TVector3 p3(void) const {return theMomentum;}
  TLorentzVector p4(double mass=pionMass);
  TVector3 getMomentum(void) const {return theMomentum;}
  ROOT::Math::XYZVector getPOCA(void) const {return thePOCA;} //https://root.cern.ch/doc/master/namespaceROOT_1_1Math.html#a676e5bc512b53b6999fd0200df1e81ab
  double getPt(void) const {return theMomentum.Perp(); }
  double getEta(void) const {return theMomentum.Eta(); }
  double getPhi(void) const {return theMomentum.Phi(); }
  double getX0(void) const {return thePOCA.X();}
  double getY0(void) const {return thePOCA.Y();}
  double getZ0(void) const {return thePOCA.Z();}
  double getRInv(void) const {return theRInv;}
  int getCharge(void);
  string getQ(void);  
  double getChi2(void) const {return theChi2;}
  double getChi2Red(void) const {return theChi2Red;}
  double getStubPtConsistency(void) const {return theStubPtConsistency;}
  bool getIsGenuine(void) const {return theIsGenuine;}
  bool getIsUnknown(void) const {return theIsUnknown;}
  bool getIsCombinatoric(void) const {return theIsCombinatoric;}
  bool getIsLoose(void) const {return theIsLoose;}
  bool getIsFake(void) const {return theIsFake;}
  unsigned int getNumOfStubs(void) const {return theNStubs;}
  unsigned int getNumOfStubsPS(void) const {return theNStubsPS;}
  unsigned int getNumOfBarrelStubs(void) const {return theNStubsBarrel;}
  unsigned int getNumOfEndcapStubs(void) const {return theNStubsEndcap;}
  bool getTPIndex(void) const {return theTPIndex;}
  vector<unsigned int> getStubsIsPS(void) const {return theStubs_isPS;}
  vector<unsigned int> getStubsDisk(void) const {return theStubs_iDisk;}
  vector<unsigned int> getStubsLayer(void) const {return theStubs_iLayer;}
  vector<unsigned int> getStubsPhi(void) const {return theStubs_iPhi;}
  vector<unsigned int> getStubsRing(void) const {return theStubs_iRing;}
  vector<unsigned int> getStubsSide(void) const {return theStubs_iSide;}
  vector<unsigned int> getStubsZ(void) const {return theStubs_iZ;}
  double getD0(void);
  int getDOF(void);
  unsigned int getNumOfStubs(void);
  unsigned int getNumOfStubsPS(void);
  unsigned int getNumOfBarrelStubs(void);
  unsigned int getNumOfEndcapStubs(void);


  void PrintProperties(void);
  void PrintAllProperties(void);
  // int getL1Track() const {return 0; } //TTPixelTrack
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  unsigned short theIndex;
  TVector3 theMomentum;
  ROOT::Math::XYZVector thePOCA;
  double theZ0;
  double theRInv;
  int theCharge;
  double theChi2;
  double theChi2Red;
  double theStubPtConsistency;
  bool theIsGenuine;
  bool theIsUnknown;
  bool theIsCombinatoric;
  bool theIsLoose;
  bool theIsFake;
  unsigned int theNStubs;
  unsigned int theNStubsPS;
  unsigned int theNStubsBarrel;
  unsigned int theNStubsEndcap;
  int theTPIndex;
  vector<unsigned int> theStubs_isPS;
  vector<unsigned int> theStubs_iDisk;
  vector<unsigned int> theStubs_iLayer;
  vector<unsigned int> theStubs_iPhi;
  vector<unsigned int> theStubs_iRing;
  vector<unsigned int> theStubs_iSide;
  vector<unsigned int> theStubs_iZ;
  int theFitParameters;
  double theD0;
  int theDOF;
  unsigned int theStubs;
  unsigned int theStubsPS;

  
};

#endif

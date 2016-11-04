#ifndef TTPixelTrack_h
#define TTPixelTrack_h

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

class TTPixelTrack{     
 public:
  // Constructors/Destructors
  TTPixelTrack();
  TTPixelTrack(unsigned short aIndex,
	       TTTrack aTTTrack,
	       TVector3 aMomentum,
	       ROOT::Math::XYZVector aPOCA,
	       double aRInv,
	       double aChi2,
	       double aSigmaRInv, 
	       double aSigmaPhi0,
	       double aSigmaD0,
	       double aSigmaT,
	       double aSigmaZ0,
	       vector<ROOT::Math::XYZVector> pixHits,
	       vector<ROOT::Math::XYZVector> candidatePixHits);	      
  ~TTPixelTrack();

  // Function declaration
  void init(TVector3 aMomentum,
	    ROOT::Math::XYZVector aPOCA,
	    double aRInv,
	    double aChi2,
	    int nhit,
	    double sigmarinv, 
	    double sigmaphi0,
	    double sigmad0,
	    double sigmat,
	    double sigmaz0,
	    vector<ROOT::Math::XYZVector> pixHits,
	    vector<ROOT::Math::XYZVector> candidatePixHits); // Ported from CMSSW. For Pixel Re-fitting  

  unsigned short index() const {return theIndex;}
  TLorentzVector p4(double mass=pionMass);
  TVector3 p3() const {return theMomentum;}
  TVector3 getMomentum() const {return theMomentum;}
  double getPt() const {return theMomentum.Perp();}
  double getEta() const {return theMomentum.Eta();}
  double getPhi() const {return theMomentum.Phi();}
  ROOT::Math::XYZVector getPOCA() const {return thePOCA;}
  double getChi2() const {return theChi2;}
  double getChi2Red() const {return theChi2Red;} 
  double getD0();
  double getRInv() const {return theRInv;}
  double getSigmaD0() const {return theSigmaD0;}
  double getSigmaPhi0() const {return theSigmaPhi0;}
  double getSigmaRInv() const {return theSigmaRInv;}
  double getSigmaT() const {return theSigmaT;}
  double getSigmaZ0() const {return theSigmaZ0;}
  double getZ0() const {return thePOCA.Z();}
  int getCandidatePixelHitsPattern(void);
  int getCharge();
  string getQ();
  int getL1Track() const {return 0;}
  int getNcandidatehit() const {return thencandidatehit;}
  int getNhit() const {return thenhit;}  
  int getPixelHitType(ROOT::Math::XYZVector pixHit);   
  int getPixelHitsPattern(void);
  std::vector<ROOT::Math::XYZVector> getCandidatePixelHits() const {return theCandidatePixelHits;}
  std::vector<ROOT::Math::XYZVector> getPixelHits() const {return thePixelHits;}
  void PrintAllProperties(void);
  void PrintProperties(void);

  std::vector<double> getCandPixHitsPhi(void) const {return candPixHits_Phi;}
  std::vector<double> getCandPixHitsR(void) const {return candPixHits_R;}
  std::vector<double> getCandPixHitsType(void) const {return candPixHits_Type;}  
  std::vector<double> getCandPixHitsZ(void) const {return candPixHits_Z;}

  std::vector<double> getPixHitsPhi(void) const {return pixHits_Phi;}
  std::vector<double> getPixHitsR(void) const {return pixHits_R;}
  std::vector<double> getPixHitsType(void) const {return candPixHits_Type;} 
  std::vector<double> getPixHitsZ(void) const {return pixHits_Z;}

  
 private:
  void _FillAuxPixelHitVariables(void);
  void _InitVars(void);
  AuxTools auxTools;
  
  // Variable declaration
  unsigned short theIndex;
  TTTrack theTTTrack;
  TVector3 theMomentum;
  ROOT::Math::XYZVector thePOCA;
  double theChi2;
  double theChi2Red;
  double theD0;
  double theRInv;
  double theSigmaD0;
  double theSigmaPhi0;
  double theSigmaRInv;
  double theSigmaT;
  double theSigmaZ0;
  double theZ0;
  int theCandidatePixelHitsPattern;
  int theCharge;
  int thePixelHitsPattern;
  int thencandidatehit;
  int thenhit;
  std::vector<ROOT::Math::XYZVector> theCandidatePixelHits;
  std::vector<ROOT::Math::XYZVector> thePixelHits;
  std::vector<double> candPixHits_Phi;
  std::vector<double> candPixHits_R;
  std::vector<double> candPixHits_Type;  
  std::vector<double> candPixHits_Z;
  std::vector<double> pixHits_Phi;
  std::vector<double> pixHits_R;
  std::vector<double> pixHits_Type;
  std::vector<double> pixHits_Z;

};

#endif

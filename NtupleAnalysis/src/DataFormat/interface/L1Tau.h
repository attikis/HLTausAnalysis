#ifndef L1Tau_h
#define L1Tau_h

// System
#include <iostream>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../Auxiliary/src/Table.C"

// ROOT
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

using namespace std;

class L1Tau{     
 public:
  // Constructors/Destructors
  L1Tau();
  L1Tau(unsigned short Index,
		double Et,
		double Eta,
		double Phi,
		short int IEt,
		short int IEta,
		short int IPhi,
		short int Iso,
		short int Bx,
		short int TowerIPhi,
		short int TowerIEta,
		short int RawEt,
		short int IsoEt,
		short int NTT,
		short int HasEM,
		short int IsMerged,
		short int HwQual);
  ~L1Tau();
  
  void PrintProperties(bool printHeader=true);
  unsigned short index(void) const {return theIndex;}
  TLorentzVector p4(void) const {return theP4;}
  double et(void) const {return theEt;}
  double eta(void) const {return theEta;}
  double phi(void) const {return thePhi;}

  void setEt(double Et) {theEt = Et;}

  unsigned short getIndex(void) const {return theIndex;}
  double getEt(void) const {return theEt;}
  double getEta(void) const {return theEta;}
  double getPhi(void) const {return thePhi;}
  short int getIEt(void) const {return theIEt;}
  short int getIEta(void) const {return theIEta;}
  short int getIPhi(void) const {return theIPhi;}
  short int getIso(void) const {return theIso;}
  short int getBx(void) const {return theBx;}
  short int getTowerIPhi(void) const {return theTowerIPhi;}
  short int getTowerIEta(void) const {return theTowerIEta;}
  short int getRawEt(void) const {return theRawEt;}
  short int getIsoEt(void) const {return theIsoEt;}
  short int getNTT(void) const {return theNTT;}
  short int getHasEM(void) const {return theHasEM;}
  short int getIsMerged(void) const {return theIsMerged;}
  short int getHwQual(void) const {return theHwQual;}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  double theEt;
  double theEta;
  double thePhi;
  short int theIEt;
  short int theIEta;
  short int theIPhi;
  short int theIso;
  short int theBx;
  short int theTowerIPhi;
  short int theTowerIEta;
  short int theRawEt;
  short int theIsoEt;
  short int theNTT;
  short int theHasEM;
  short int theIsMerged;
  short int theHwQual;
  
};

#endif

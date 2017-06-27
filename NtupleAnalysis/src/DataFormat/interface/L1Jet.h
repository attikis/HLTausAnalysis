#ifndef L1Jet_h
#define L1Jet_h

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

class L1Jet{     
 public:
  // Constructors/Destructors
  L1Jet();
  L1Jet(unsigned short Index,
	float Et,
	float Eta,
	float Phi,
	short IET,
	short IEta,
	short IPhi,
	short Bx,
	short RawEt,
	short SeedEt,
	short TowerIEta,
	short TowerIPhi,
	short PUEt,
	short PUDonutEt0,
	short PUDonutEt1,
	short PUDonutEt2,
	short PUDonutEt3);
  
  ~L1Jet();
  
  void PrintProperties(bool printHeader=true);
  unsigned short index(void) const {return theIndex;}
  TLorentzVector p4(void) const {return theP4;}
  double et(void) const {return theEt;}
  double eta(void) const {return theEta;}
  double phi(void) const {return thePhi;}

  unsigned short getIndex(void) const {return theIndex;}
  double getEt(void) const {return theEt;}
  double getEta(void) const {return theEta;}
  double getPhi(void) const {return thePhi;}
  short getIET(void) const {return theIET;}
  short getIEta(void) const {return theIEta;}
  short getIPhi(void) const {return theIPhi;}
  short getBx(void) const {return theBx;}
  short getRawEt(void) const {return theRawEt;}
  short getSeedEt(void) const {return theSeedEt;}
  short getTowerIEta(void) const {return theTowerIEta;}
  short getTowerIPhi(void) const {return theTowerIPhi;}
  short getPUEt(void) const {return thePUEt;}
  short getPUDonutEt0(void) const {return thePUDonutEt0;}
  short getPUDonutEt1(void) const {return thePUDonutEt1;}
  short getPUDonutEt2(void) const {return thePUDonutEt2;}
  short getPUDonutEt3(void) const {return thePUDonutEt3;}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  float theEt;
  float theEta;
  float thePhi;
  short theIET;
  short theIEta;
  short theIPhi;
  short theBx;
  short theRawEt;
  short theSeedEt;
  short theTowerIEta;
  short theTowerIPhi;
  short thePUEt;
  short thePUDonutEt0;
  short thePUDonutEt1;
  short thePUDonutEt2;
  short thePUDonutEt3;
  
};

#endif

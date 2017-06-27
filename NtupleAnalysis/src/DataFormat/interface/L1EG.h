#ifndef L1EG_h
#define L1EG_h

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

class L1EG{     
 public:
  // Constructors/Destructors
  L1EG();
  L1EG(unsigned short Index,
       double Et,
       double Eta,
       double Phi,
       int IEt,
       int IEta,
       int IPhi,
       int Iso,
       int Bx,
       int TowerIPhi,
       int TowerIEta,
       int RawEt,
       int IsoEt,
       int FootprintEt,
       int NTT,
       int Shape,
       int TowerHoE);
  
  ~L1EG();
  
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
  int getIET(void) const {return theIEt;}
  int getIEta(void) const {return theIEta;}
  int getIPhi(void) const {return theIPhi;}
  int getIso(void) const {return theIso;}
  int getBx(void) const {return theBx;}
  int getTowerIPhi(void) const {return theTowerIPhi;}
  int getTowerIEta(void) const {return theTowerIEta;}
  int getRawEt(void) const {return theRawEt;}
  int getIsoEt(void) const {return theIsoEt;}
  int getFootprintEt(void) const {return theFootprintEt;}
  int getNTT(void) const {return theNTT;}
  int getShape(void) const {return theShape;}
  int getTowerHoE(void) const {return theTowerHoE;}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  float theEt;
  float theEta;
  float thePhi;
  int theIEt;
  int theIEta;
  int theIPhi;
  int theIso;
  int theBx;
  int theTowerIPhi;
  int theTowerIEta;
  int theRawEt;
  int theIsoEt;
  int theFootprintEt;
  int theNTT;
  int theShape;
  int theTowerHoE;
  
};

#endif

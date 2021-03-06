#ifndef L1TKEM_h
#define L1TKEM_h

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

class L1TKEM{     
 public:
  // Constructors/Destructors
  L1TKEM();
  L1TKEM(unsigned short Index,
       double Et,
       double Eta,
       double Phi,
       double EGRefPt,
       double EGRefEta,
       double EGRefPhi,
       double TrkIso,
       int    Bx,
       double HwQual,
       double zVtx);
  
  ~L1TKEM();
  
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
  double getEGRefPt(void) const {return theEGRefPt;}
  double getEGRefEta(void) const {return theEGRefEta;}
  double getEGRefPhi(void) const {return theEGRefPhi;}
  double getTrkIso(void) const {return theTrkIso;}
  int getBx(void) const {return theBx;}
  double getHwQual(void) const {return theHwQual;}
  double getzVtx(void) const {return thezVtx;}

  bool operator<(const L1TKEM& other) const {return this->et() > other.et();}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  double theEt;
  double theEta;
  double thePhi;
  double theEGRefPt;
  double theEGRefEta;
  double theEGRefPhi;
  double theTrkIso;
  int theBx;
  double theHwQual;
  double thezVtx;
  
};

#endif

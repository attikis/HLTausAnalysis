#ifndef EG_h
#define EG_h

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

class EG{     
 public:
  // Constructors/Destructors
  EG();
  EG(unsigned short Index,
       double Et,
       double Eta,
       double Phi,
       double Iso,
       int    Bx,
       int HwQual,
       double zVtx);
  
  ~EG();
  
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
  double getIso(void) const {return theIso;}
  int getBx(void) const {return theBx;}
  int getHwQual(void) const {return theHwQual;}
  double getzVtx(void) const {return thezVtx;}

  bool operator<(const EG& other) const {return this->et() > other.et();}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  double theEt;
  double theEta;
  double thePhi;
  double theIso;
  int theBx;
  int theHwQual;
  double thezVtx;
  
};

#endif

#ifndef L1CaloTP_h
#define L1CaloTP_h

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

class L1CaloTP{     
 public:
  // Constructors/Destructors
  L1CaloTP();
  L1CaloTP(unsigned short Index,
       double Et,
       double CompEt,
       int IEta,
       int IPhi,
       int CalIPhi);
  
  ~L1CaloTP();
  
  void PrintProperties(bool printHeader=true);
  unsigned short index(void) const {return theIndex;}
  unsigned short getIndex(void) const {return theIndex;}
  TLorentzVector p4(void) const {return theP4;}
  double et(void) const {return theEt;}
  double getEt(void) const {return theEt;}
  double getCompEt(void) const {return theCompEt;}
  double eta(void) const {return theEta;}
  double getEta(void) const {return theEta;}
  double phi(void) const {return thePhi;}
  double getPhi(void) const {return thePhi;}
  int getIEta(void) const {return theIEta;}
  int getIPhi(void) const {return theIPhi;}
  int getCalIPhi(void) const {return theCalIPhi;}
  double getEtaFromIEta(int) const {return 0.0;} // FIXME
  double getPhiFromIPhi(int) const {return 0.0;} // FIXME

  bool operator<(const L1CaloTP& other) const {return this->et() > other.et();}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  double theEt;
  double theCompEt;
  int    theIEta;
  int    theIPhi;
  int    theCalIPhi;
  double theEta;
  double thePhi;
};

#endif

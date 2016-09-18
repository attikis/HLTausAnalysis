#ifndef L1JetParticle_h
#define L1JetParticle_h

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

class L1JetParticle{     
 public:
  // Constructors/Destructors
  L1JetParticle();
  L1JetParticle(unsigned short Index,
		double E,
		double Et,
		double Eta,
		double Phi,
		double Bx,
		double Type);
  ~L1JetParticle();
  
  void PrintProperties(void);
  unsigned short index(void) const {return theIndex;}
  TLorentzVector p4(void) const {return theP4;}
  double energy(void) const {return theEnergy;}
  double et(void) const {return theEt;}
  double eta(void) const {return theEta;}
  double phi(void) const {return thePhi;}
  double bx(void) const {return theBx;}
  double type(void) const {return theType;}
  
 private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  double theEnergy;
  double theEt;
  double theEta;
  double thePhi;
  double theBx;
  double theType;
  
};

#endif

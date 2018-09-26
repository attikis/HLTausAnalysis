#ifndef L1Sum_h
#define L1Sum_h

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

class L1Sum{     
 public:
  // Constructors/Destructors
  L1Sum();
  L1Sum(unsigned short Index,
	float Et,
	float Phi,
	short int IEt,
	short int IPhi,
	short int Type,
	short int Bx);

~L1Sum();

void PrintProperties(bool printHeader=true);
unsigned short getIndex(void) const {return theIndex;}
double getEt(void) const {return theEt;}
double getPhi(void) const {return thePhi;}
int getIEt(void) const {return theIEt;}
int getIPhi(void) const {return theIPhi;}
int getType(void) const {return theType;}
int getBx(void) const {return theBx;}

private:
  void _InitVars(void);
  AuxTools auxTools;

  // Variable declaration
  unsigned short theIndex;
  float theEt;
  float thePhi;
  short int theIEt;
  short int theIPhi;
  short int theType;
  short int theBx;
    
};

#endif

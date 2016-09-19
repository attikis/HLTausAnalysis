#ifndef GenParticle_h
#define GenParticle_h

// System
#include <iostream>
#include <vector>
#include <algorithm>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../Auxiliary/src/Table.C"

// ROOT
#include <TROOT.h>
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

using namespace std;

class GenParticle{     
 public:
  // Constructors/Destructors
  GenParticle();
  GenParticle(unsigned short index,
	      double Pt,
	      double Eta,
	      double Phi,
	      double Mass,
	      int  Charge,
	      int PdgId,
	      int Status,
	      double VertexX,
	      double VertexY,
	      double VertexZ,
	      vector<unsigned short> MothersIndex,
	      vector<unsigned short> DaughtersIndex);
  ~GenParticle();
  
  void PrintProperties(void);
  unsigned short index(void) const {return theIndex;}
  TLorentzVector p4(void) const {return theP4;}
  double energy(void) const {return theP4.Energy();}
  double et(void) const {return theP4.Et();}
  double pt(void) const {return theP4.Pt();}
  double eta(void) const {return theP4.Eta();}
  double phi(void) const {return theP4.Phi();}
  double mass(void) const {return theP4.M();}
  int charge(void) const {return theCharge;}
  int pdgId(void) const {return thePdgId;}
  int status(void) const {return theStatus;}
  double vx(void) const {return theVertex.X();}
  double vy(void) const {return theVertex.Y();}
  double vz(void) const {return theVertex.Z();}
  ROOT::Math::XYZVector vertex(void) const {return theVertex;}
  vector<unsigned short> mothersIndex(void) const {return theMothersIndex;}
  vector<unsigned short> daughtersIndex(void) const {return theDaughtersIndex;}
  vector<GenParticle> mothers(void) const {return theMothers;}
  vector<GenParticle> daughters(void) const {return theDaughters;}
  void SetMothers(vector<GenParticle> mothers) {theMothers = mothers;}
  void SetDaughters(vector<GenParticle> daughters){theDaughters = daughters;}
  
 private:
  AuxTools auxTools;

  // Variable declaration
  TLorentzVector theP4;
  unsigned short theIndex;
  double thePt;
  double theEta;
  double thePhi;
  double theMass;
  int theCharge;
  int thePdgId;
  int theStatus;
  double theVertexX;
  double theVertexY;
  double theVertexZ;
  ROOT::Math::XYZVector theVertex;
  vector<unsigned short> theMothersIndex;
  vector<unsigned short> theDaughtersIndex;
  vector<GenParticle> theMothers;
  vector<GenParticle> theDaughters;

  
};

#endif

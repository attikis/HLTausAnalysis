#ifndef GenParticle_cxx
#define GenParticle_cxx

// User
#include "../interface/GenParticle.h"

//****************************************************************************
GenParticle::GenParticle()
//****************************************************************************
{

}

//****************************************************************************
GenParticle::GenParticle(unsigned short Index,
			 double Pt,
			 double Eta,
			 double Phi,
			 double Mass,
			 int Charge,
			 int PdgId,
			 int Status,
			 double VertexX,
			 double VertexY,
			 double VertexZ,
			 vector<unsigned short> MothersIndex,
			 vector<unsigned short> DaughtersIndex)
//****************************************************************************
{

  theIndex     = Index;
  thePt        = Pt;
  theEta       = Eta;
  thePhi       = Phi;
  theMass      = Mass;
  theCharge    = Charge;
  thePdgId     = PdgId;
  theStatus    = Status;
  theVertexX   = VertexX;
  theVertexY   = VertexY;
  theVertexZ   = VertexZ;
  theMothersIndex   = MothersIndex;
  theDaughtersIndex = DaughtersIndex;
  theP4.SetPtEtaPhiM(thePt, theEta, thePhi, theMass);
  theVertex.SetXYZ(VertexX, VertexY, VertexZ);
  if (0) PrintProperties();
  
}


//****************************************************************************
GenParticle::~GenParticle()
//****************************************************************************
{

}



//****************************************************************************
void GenParticle::PrintProperties(void)
//****************************************************************************
{
  
  Table info("Index | Pt | Eta | Phi | Mass | Charge | PdgId | Status | vX | vY | vZ | mothers | daughters |", "Text");
  info.AddRowColumn(0, auxTools.ToString( index() , 1) );
  info.AddRowColumn(0, auxTools.ToString( pt()    , 3) );
  info.AddRowColumn(0, auxTools.ToString( eta()   , 3) );
  info.AddRowColumn(0, auxTools.ToString( phi()   , 3) );
  info.AddRowColumn(0, auxTools.ToString( mass()  , 3) );
  info.AddRowColumn(0, auxTools.ToString( charge(), 1) );
  info.AddRowColumn(0, auxTools.ToString( pdgId() , 3) );
  info.AddRowColumn(0, auxTools.ToString( status(), 1) );
  info.AddRowColumn(0, auxTools.ToString( vx()    , 3) );
  info.AddRowColumn(0, auxTools.ToString( vy()    , 3) );
  info.AddRowColumn(0, auxTools.ToString( vz()    , 3) );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(mothersIndex()) );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(daughtersIndex()) );
  info.Print();

  return;
}


#endif

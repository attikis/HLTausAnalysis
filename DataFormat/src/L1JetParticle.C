#ifndef L1JetParticle_cxx
#define L1JetParticle_cxx

// User
#include "../interface/L1JetParticle.h"

//****************************************************************************
L1JetParticle::L1JetParticle()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void L1JetParticle::_InitVars(void)
//****************************************************************************
{

  theIndex  = 0;  
  theEnergy = 0.0;
  theEt     = 0.0;
  theEta    = 0.0;
  thePhi    = 0.0;
  theBx     = 0.0;
  theType   = 0.0;

  return;
}

//****************************************************************************
L1JetParticle::~L1JetParticle()
//****************************************************************************
{

}


//****************************************************************************
L1JetParticle::L1JetParticle(unsigned short Index,
			     double E,
			     double Et,
			     double Eta,
			     double Phi,
			     double Bx,
			     double Type)
//****************************************************************************
{

  _InitVars();
  theIndex  = Index;
  theEnergy = E;
  theEt     = Et;
  theEta    = Eta;
  thePhi    = Phi;
  theBx     = Bx;
  theType   = Type;
  theP4.SetPtEtaPhiE(Et, Eta, Phi, E); // WARNING: Approximation since we need Pt and Not Et 

  
  if (0) PrintProperties();
}

//****************************************************************************
void L1JetParticle::PrintProperties(void)
//****************************************************************************
{
  
  Table info("Energy | Et | Eta | Phi | Bx | Type", "Text");

  info.AddRowColumn(0, auxTools.ToString( energy(), 4) );
  info.AddRowColumn(0, auxTools.ToString( et(), 4)     );
  info.AddRowColumn(0, auxTools.ToString( eta(), 4)    );
  info.AddRowColumn(0, auxTools.ToString( phi(), 4)    );
  info.AddRowColumn(0, auxTools.ToString( bx(), 4)     );
  info.AddRowColumn(0, auxTools.ToString( type(), 4)   );

  info.Print();

  return;
}


#endif

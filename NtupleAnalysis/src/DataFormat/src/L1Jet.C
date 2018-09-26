#ifndef L1Jet_cxx
#define L1Jet_cxx

// User
#include "../interface/L1Jet.h"

//****************************************************************************
L1Jet::L1Jet()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void L1Jet::_InitVars(void)
//****************************************************************************
{

  theIndex      = 0;  
  theEt         = 0.0;
  theEta        = 0.0;
  thePhi        = 0.0;
  theIEt        = 0;
  theIEta       = 0;
  theIPhi       = 0;
  theBx         = 0;
  theRawEt      = 0;
  theSeedEt     = 0;
  theTowerIEta  = 0;
  theTowerIPhi  = 0;
  thePUEt       = 0;
  thePUDonutEt0 = 0;
  thePUDonutEt1 = 0;
  thePUDonutEt2 = 0;
  thePUDonutEt3 = 0;

  return;
}

//****************************************************************************
L1Jet::~L1Jet()
//****************************************************************************
{

}


//****************************************************************************
L1Jet::L1Jet(unsigned short Index,
	     float Et,
	     float Eta,
	     float Phi,
	     short IEt,
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
	     short PUDonutEt3)
//****************************************************************************
{

  _InitVars();
  theIndex      = Index;
  theEt         = Et;
  theEta        = Eta;
  thePhi        = Phi;
  theIEt        = IEt;
  theIEta       = IEta;
  theIPhi       = IPhi;
  theBx         = Bx;
  theRawEt      = RawEt;
  theSeedEt     = SeedEt;
  theTowerIEta  = TowerIEta;
  theTowerIPhi  = TowerIPhi;
  thePUEt       = PUEt;
  thePUDonutEt0 = PUDonutEt0;
  thePUDonutEt1 = PUDonutEt1;
  thePUDonutEt2 = PUDonutEt2;
  thePUDonutEt3 = PUDonutEt3;
  
  theP4.SetPtEtaPhiE(Et, Eta, Phi, Et); // WARNING: Approximation since we need Pt and Not Et 
  
  if (0) PrintProperties();
}

//****************************************************************************
void L1Jet::PrintProperties(bool printHeader)
//****************************************************************************
{
  
  Table info("Index | Et | Eta | Phi | IET | IEta | IPhi | Bx | Raw Et | Seed Et | TowerIPhi | TowerIEta | PU Et | PUDonutEt0 | PUDonutEt1 | PUDonutEt2 | PUDonutEt3", "Text");

  info.AddRowColumn(0, auxTools.ToString( getIndex() )        );
  info.AddRowColumn(0, auxTools.ToString( getEt(), 3)         );
  info.AddRowColumn(0, auxTools.ToString( getEta(), 3)        );
  info.AddRowColumn(0, auxTools.ToString( getPhi(), 3)        );
  info.AddRowColumn(0, auxTools.ToString( getIEt() , 3)       );
  info.AddRowColumn(0, auxTools.ToString( getIEta(), 3)       );
  info.AddRowColumn(0, auxTools.ToString( getIPhi(), 3)       );
  info.AddRowColumn(0, auxTools.ToString( getBx() , 3)        );
  info.AddRowColumn(0, auxTools.ToString( getRawEt(), 3)      );
  info.AddRowColumn(0, auxTools.ToString( getSeedEt(), 3)     );
  info.AddRowColumn(0, auxTools.ToString( getTowerIPhi(), 3)  );
  info.AddRowColumn(0, auxTools.ToString( getTowerIEta(), 3)  );
  info.AddRowColumn(0, auxTools.ToString( getPUEt(), 3)       );
  info.AddRowColumn(0, auxTools.ToString( getPUDonutEt0(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getPUDonutEt1(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getPUDonutEt2(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getPUDonutEt3(), 3) );

  info.Print(printHeader);

  return;
}


#endif

#ifndef L1Tau_cxx
#define L1Tau_cxx

// User
#include "../interface/L1Tau.h"

//****************************************************************************
L1Tau::L1Tau()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void L1Tau::_InitVars(void)
//****************************************************************************
{

  theIndex     = 0;  
  theEt        = 0.0;
  theEta       = 0.0;
  thePhi       = 0.0;
  theIEt       = 0;
  theIEta      = 0;
  theIPhi      = 0;
  theIso       = 0;
  theBx        = 0;
  theTowerIPhi = 0;
  theTowerIEta = 0;
  theRawEt     = 0;
  theIsoEt     = 0;
  theNTT       = 0;
  theHasEM     = 0;
  theIsMerged  = 0;
  theHwQual    = 0;

  return;
}

//****************************************************************************
L1Tau::~L1Tau()
//****************************************************************************
{

}


//****************************************************************************
L1Tau::L1Tau(unsigned short Index,
	     double Et,
	     double Eta,
	     double Phi,
	     short int IEt,
	     short int IEta,
	     short int IPhi,
	     short int Iso,
	     short int Bx,
	     short int TowerIPhi,
	     short int TowerIEta,
	     short int RawEt,
	     short int IsoEt,
	     short int NTT,
	     short int HasEM,
	     short int IsMerged,
	     short int HwQual)
//****************************************************************************
{

  _InitVars();
  theIndex     = Index;
  theEt        = Et;
  theEta       = Eta;
  thePhi       = Phi;
  theIEt       = IEt;
  theIEta      = IEta;
  theIPhi      = IPhi;
  theIso       = Iso;
  theBx        = Bx;
  theTowerIPhi = TowerIPhi;
  theTowerIEta = TowerIEta;
  theRawEt     = RawEt;
  theIsoEt     = IsoEt;
  theNTT       = NTT;
  theHasEM     = HasEM;
  theIsMerged  = IsMerged;
  theHwQual    = HwQual;

  theP4.SetPtEtaPhiE(Et, Eta, Phi, Et); // WARNING: Approximation since we need Pt and Not Et 

  
  if (0) PrintProperties();
}

//****************************************************************************
void L1Tau::PrintProperties(bool printHeader)
//****************************************************************************
{
  Table info("Index | Et | Eta | Phi | IEt | IEta | IPhi | Iso | Bx | TowerIPhi | TowerIEta | RawEt | IsoEt | NTT | HasEM | IsMerged | HwQual", "Text");
    
  info.AddRowColumn(0, auxTools.ToString( getIndex() )       );
  info.AddRowColumn(0, auxTools.ToString( getEt(), 2)        );
  info.AddRowColumn(0, auxTools.ToString( getEta(), 3)       );
  info.AddRowColumn(0, auxTools.ToString( getPhi(), 3)       );
  info.AddRowColumn(0, auxTools.ToString( getIEt() , 3)      );
  info.AddRowColumn(0, auxTools.ToString( getIEta(), 3)      );
  info.AddRowColumn(0, auxTools.ToString( getIPhi(), 3)      );
  info.AddRowColumn(0, auxTools.ToString( getIso() , 3)      );
  info.AddRowColumn(0, auxTools.ToString( getBx()  , 3)      );
  info.AddRowColumn(0, auxTools.ToString( getTowerIPhi(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getTowerIEta(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getRawEt(), 3)     );
  info.AddRowColumn(0, auxTools.ToString( getIsoEt(), 3)     );
  info.AddRowColumn(0, auxTools.ToString( getNTT(), 3)       );
  info.AddRowColumn(0, auxTools.ToString( getHasEM(), 3)     );
  info.AddRowColumn(0, auxTools.ToString( getIsMerged(), 3)  );
  info.AddRowColumn(0, auxTools.ToString( getHwQual(), 3)    );

  info.Print(printHeader);

  return;
}


#endif

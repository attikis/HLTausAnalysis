#ifndef L1EG_cxx
#define L1EG_cxx

// User
#include "../interface/L1EG.h"

//****************************************************************************
L1EG::L1EG()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void L1EG::_InitVars(void)
//****************************************************************************
{
  theP4.SetPtEtaPhiE(0, 0, 0, 0);
  theIndex       = 0;  
  theEt          = 0.0;
  theEta         = 0.0;
  thePhi         = 0.0;
  theIEt         = 0;
  theIEta        = 0;
  theIPhi        = 0;
  theIso         = 0;
  theBx          = 0;
  theTowerIPhi   = 0;
  theTowerIEta   = 0;
  theRawEt       = 0;
  theIsoEt       = 0;
  theFootprintEt = 0;
  theNTT         = 0;
  theShape       = 0;
  theTowerHoE    = 0;

  return;
}

//****************************************************************************
L1EG::~L1EG()
//****************************************************************************
{

}


//****************************************************************************
L1EG::L1EG(unsigned short Index,
	   double Et,
	   double Eta,
	   double Phi,
	   int IEt,
	   int IEta,
	   int IPhi,
	   int Iso,
	   int Bx,
	   int TowerIPhi,
	   int TowerIEta,
	   int RawEt,
	   int IsoEt,
	   int FootprintEt,
	   int NTT,
	   int Shape,
	   int TowerHoE)
//****************************************************************************
{

  _InitVars();
  theIndex       = Index;
  theEt          = Et;
  theEta         = Eta;
  thePhi         = Phi;
  theIEt         = IEt;
  theIEta        = IEta;
  theIPhi        = IPhi;
  theIso         = Iso;
  theBx          = Bx;
  theTowerIPhi   = TowerIPhi;
  theTowerIEta   = TowerIEta;
  theRawEt       = RawEt;
  theIsoEt       = IsoEt;
  theFootprintEt = FootprintEt;
  theNTT         = NTT;
  theShape       = Shape;
  theTowerHoE    = TowerHoE;
  theP4.SetPtEtaPhiM( sqrt(Et*Et-pionMass*pionMass), Eta, Phi, pionMass); // WARNING: Assumes pion mass for EG clusters 
  
  if (0) PrintProperties();
}

//****************************************************************************
void L1EG::PrintProperties(bool printHeader)
//****************************************************************************
{
  
  Table info("Index | Et | Eta | Phi | IET | IEta | IPhi | Iso | Bx | TowerIPhi | TowerIEta | Raw Et | Iso Et| Footprint Et | NTT | Shape | TowerHoE", "Text");

  info.AddRowColumn(0, auxTools.ToString( getIndex() ) );
  info.AddRowColumn(0, auxTools.ToString( getEt()         , 3) );
  info.AddRowColumn(0, auxTools.ToString( getEta()        , 3) );
  info.AddRowColumn(0, auxTools.ToString( getPhi()        , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIEt()        , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIEta()       , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIPhi()       , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIso()        , 3) );
  info.AddRowColumn(0, auxTools.ToString( getBx()         , 3) );
  info.AddRowColumn(0, auxTools.ToString( getTowerIPhi()  , 3) );
  info.AddRowColumn(0, auxTools.ToString( getTowerIEta()  , 3) );
  info.AddRowColumn(0, auxTools.ToString( getRawEt()      , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIsoEt()      , 3) );
  info.AddRowColumn(0, auxTools.ToString( getFootprintEt(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getNTT()        , 3) );
  info.AddRowColumn(0, auxTools.ToString( getShape()      , 3) );
  info.AddRowColumn(0, auxTools.ToString( getTowerHoE()   , 3) );

  info.Print(printHeader);

  return;
}


#endif

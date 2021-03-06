#ifndef L1TKEM_cxx
#define L1TKEM_cxx

// User
#include "../interface/L1TKEM.h"

//****************************************************************************
L1TKEM::L1TKEM()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void L1TKEM::_InitVars(void)
//****************************************************************************
{
  theP4.SetPtEtaPhiE(0, 0, 0, 0);
  theIndex       = 0;  
  theEt          = 0.0;
  theEta         = 0.0;
  thePhi         = 0.0;
  theEGRefPt     = 0.0;
  theEGRefEta    = 0.0;
  theEGRefPhi    = 0.0;
  theTrkIso      = 0;
  theBx          = 0;
  theHwQual      = 0;
  thezVtx        = 0;

  return;
}

//****************************************************************************
L1TKEM::~L1TKEM()
//****************************************************************************
{

}


//****************************************************************************
L1TKEM::L1TKEM(unsigned short Index,
	   double Et,
	   double Eta,
	   double Phi,
           double EGRefPt,
           double EgREfEta,
           double EGRefPhi,
	   double TrkIso,
	   int Bx,
           double HwQual,
           double zVtx)
//****************************************************************************
{

  _InitVars();
  theIndex       = Index;
  theEt          = Et;
  theEta         = Eta;
  thePhi         = Phi;
  theTrkIso      = TrkIso;
  theBx          = Bx;
  thezVtx        = zVtx;
  theHwQual      = HwQual;

  theP4.SetPtEtaPhiM( sqrt(Et*Et-pionMass*pionMass), Eta, Phi, pionMass); // WARNING: Assumes pion mass for EG clusters 
  
  if (0) PrintProperties();
}

//****************************************************************************
void L1TKEM::PrintProperties(bool printHeader)
//****************************************************************************
{
  
  Table info("Index | Et | Eta | Phi | IET | IEta | IPhi | Iso | Bx | TowerIPhi | TowerIEta | Raw Et | Iso Et| Footprint Et | NTT | Shape | TowerHoE", "Text");

  info.AddRowColumn(0, auxTools.ToString( getIndex() ) );
  info.AddRowColumn(0, auxTools.ToString( getEt()         , 3) );
  info.AddRowColumn(0, auxTools.ToString( getEta()        , 3) );
  info.AddRowColumn(0, auxTools.ToString( getPhi()        , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getIEt()        , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getIEta()       , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getIPhi()       , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getIso()        , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getBx()         , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getTowerIPhi()  , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getTowerIEta()  , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getRawEt()      , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getIsoEt()      , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getFootprintEt(), 3) );
//  info.AddRowColumn(0, auxTools.ToString( getNTT()        , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getShape()      , 3) );
//  info.AddRowColumn(0, auxTools.ToString( getTowerHoE()   , 3) );

  info.Print(printHeader);

  return;
}


#endif

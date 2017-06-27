#ifndef L1Sum_cxx
#define L1Sum_cxx

// User
#include "../interface/L1Sum.h"

//****************************************************************************
L1Sum::L1Sum()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void L1Sum::_InitVars(void)
//****************************************************************************
{
  theIndex = 0;
  theEt    = 0.0;
  thePhi   = 0.0;
  theIEt   = 0;
  theIPhi  = 0;
  theType  = 0;
  theBx    = 0;
  return;
}

//****************************************************************************
L1Sum::~L1Sum()
//****************************************************************************
{

}


//****************************************************************************
L1Sum::L1Sum(unsigned short Index,
	     float Et,
	     float Phi,
	     short int IEt,
	     short int IPhi,
	     short int Type,
	     short int Bx)
//****************************************************************************
{

  _InitVars();
  theIndex = Index;
  theEt    = Et;
  thePhi   = Phi;
  theIEt   = IEt;
  theIPhi  = IPhi;
  theType  = Type;
  theBx    = Bx;
  if (0) PrintProperties();
}

//****************************************************************************
void L1Sum::PrintProperties(bool printHeader)
//****************************************************************************
{
  
  Table info("Index | Et | Phi | IEt | IPhi | Type | Bx", "Text");

  info.AddRowColumn(0, auxTools.ToString( getIndex() ) );
  info.AddRowColumn(0, auxTools.ToString( getEt()  , 3) );
  info.AddRowColumn(0, auxTools.ToString( getPhi() , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIEt() , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIPhi(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getType(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getBx()  , 3) );

  info.Print(printHeader);

  return;
}


#endif

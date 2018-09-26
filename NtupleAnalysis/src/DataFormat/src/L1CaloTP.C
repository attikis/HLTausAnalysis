#ifndef L1CaloTP_cxx
#define L1CaloTP_cxx

// User
#include "../interface/L1CaloTP.h"

//****************************************************************************
L1CaloTP::L1CaloTP()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void L1CaloTP::_InitVars(void)
//****************************************************************************
{
  theP4.SetPtEtaPhiE(0, 0, 0, 0);
  theIndex       = 0;  
  theEt          = 0.0;
  theCompEt      = 0.0;
  theIEta         = 0.0;
  theIPhi         = 0.0;
  theCalIPhi      = 0.0;
  return;
}

//****************************************************************************
L1CaloTP::~L1CaloTP()
//****************************************************************************
{

}


//****************************************************************************
L1CaloTP::L1CaloTP(unsigned short Index,
	   double Et,
	   double CompEt,
	   int    IEta,
	   int    IPhi,
	   int    CalIPhi)
//****************************************************************************
{

  _InitVars();
  theIndex       = Index;
  theEt          = Et;
  theCompEt      = CompEt;
  theIEta        = IEta;
  theIPhi        = IPhi;
  theCalIPhi     = CalIPhi;
  theEta         = getEtaFromIEta(IEta);
  thePhi         = getPhiFromIPhi(IPhi);
  // theP4.SetPtEtaPhiM( sqrt(Et*Et-pionMass*pionMass), Eta, Phi, pionMass); // WARNING: Assumes pion mass for EG clusters //FIXME
  
  if (0) PrintProperties();
}

//****************************************************************************
double L1CaloTP:: getEtaFromIEta(int iEta)
//*****************************************************************************
{
  double Eta=-1;
  
  if (iEta>0){
    Eta = -(0.0870/2) + iEta*0.0870; }
  else if (iEta<0){
    Eta = +(0.0870/2) + iEta*0.0870; }

  return Eta;
}


//****************************************************************************                          
double L1CaloTP:: getPhiFromIPhi(int iPhi)
//*****************************************************************************                                                                
{
  double Phi;

  if (iPhi<37){
    Phi = -(0.0873/2) + iPhi*0.0873; }
  else {
    Phi = +(0.0873/2) + (iPhi-73)*0.873; }
  
  return Phi;
}




//****************************************************************************
void L1CaloTP::PrintProperties(bool printHeader)
//****************************************************************************
{
  
  Table info("Index | Et | CompEt | IEta | IPhi | CalIPhi | Eta | Phi", "Text");

  info.AddRowColumn(0, auxTools.ToString( getIndex() ) );
  info.AddRowColumn(0, auxTools.ToString( getEt()         , 3) );
  info.AddRowColumn(0, auxTools.ToString( getCompEt()     , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIEta()       , 3) );
  info.AddRowColumn(0, auxTools.ToString( getIPhi()       , 3) );
  info.AddRowColumn(0, auxTools.ToString( getCalIPhi()    , 3) );
  info.AddRowColumn(0, auxTools.ToString( getEta()        , 3) );
  info.AddRowColumn(0, auxTools.ToString( getPhi()        , 3) );

  info.Print(printHeader);

  return;
}


#endif

#ifndef TTTrack_cxx
#define TTTrack_cxx

// User
#include "../interface/TTTrack.h"

//****************************************************************************
TTTrack::TTTrack()
//****************************************************************************
{

  _InitVars();

}

//****************************************************************************
void TTTrack::_InitVars(void)
//****************************************************************************
{

  theIndex = 0.0;
  theChi2 = 0.0;
  theChi2Red = 0.0;
  theD0 = 0.0;
  theFitParameters = 0;
  theIsCombinatoric = false;
  theIsGenuine = false;
  theIsUnknown = false;
  theIsLoose = false;
  theIsFake = false;
  theNStubs = 0.0;
  theNStubsPS = 0.0;
  theNStubsBarrel = 0.0;
  theNStubsEndcap = 0.0;
  theMomentum.SetXYZ(0.0, 0.0, 0.0);
  thePOCA.SetXYZ(0.0, 0.0, 0.0);
  theRInv = 0.0;
  theStubPtConsistency = 0.0;
  theZ0 = 0.0;
  theCharge = -9999;
  
  return;
}

//****************************************************************************
TTTrack::~TTTrack()
//****************************************************************************
{

}


//****************************************************************************
TTTrack::TTTrack(unsigned short aIndex,
		 TVector3 aMomentum,
		 ROOT::Math::XYZVector aPOCA, // https://root.cern.ch/doc/master/Vector3DPage.html
		 double aRInv,
		 double aChi2,	    
		 double aStubPtConsistency,
		 bool isGenuine,
		 bool isUnknown,
		 bool isCombinatoric,
		 bool isLoose,
		 bool isFake,
		 unsigned int nStubs,
		 unsigned int nStubsPS,
		 unsigned int nStubsBarrel,
		 unsigned int nStubsEndcap,
		 int matchTP_index,
		 unsigned int nFitParams)
//****************************************************************************
{

  // theL1Track=aL1Track;
  _InitVars();
  
  theIndex             = aIndex;
  theChi2              = aChi2;
  theIsCombinatoric    = isCombinatoric;
  theIsGenuine         = isGenuine;
  theIsUnknown         = isUnknown;
  theIsLoose           = isLoose;
  theIsFake            = isFake;
  theNStubs            = nStubs;
  theNStubsPS          = nStubsPS;
  theNStubsBarrel      = nStubsBarrel;
  theNStubsEndcap      = nStubsEndcap;
  theTPIndex           = matchTP_index;
  theMomentum          = aMomentum;
  thePOCA              = aPOCA;
  theFitParameters     = nFitParams;
  theD0                = getD0(); // after "thePOCA" and nFitParams are assigned
  theRInv              = aRInv;
  theCharge            = getCharge(); // after "theRInv" is assigned
  theStubPtConsistency = aStubPtConsistency;
  theZ0                = getZ0();
  theStubs             = getNumOfStubs();
  theStubsPS           = getNumOfStubsPS();
  theDOF               = getDOF(); // after "theStubs" and "theFitParameters" are assigned
  theChi2Red           = double(theChi2)/double(theDOF); // after "theDOF" is assigned
  if (0) PrintProperties();
  
  return;
}


//****************************************************************************
double TTTrack::getD0(void)
//****************************************************************************
{
  if (theFitParameters == 4) return 1e6;
  theD0 = -thePOCA.X() * sin( theMomentum.Phi() ) + thePOCA.Y() * cos(theMomentum.Phi() );
  return theD0;
}


//****************************************************************************
int TTTrack::getDOF(void)
//****************************************************************************
{

  unsigned int nStubs = getNumOfStubs();
  int dof = 2 * nStubs - theFitParameters;
  return dof;
}



//****************************************************************************
int TTTrack::getCharge(void)
//****************************************************************************
{
  if (theRInv < 0.0) theCharge = -1;
  else if (theRInv > 0.0) theCharge = +1;
  else{
    std::cout << "E R R O R ! TTTrack::getCharge(...) - Invalid value for theRInv \"" << theRInv << "\". EXIT" << std::endl;
    exit(1);
  }

  return theCharge;
}


//****************************************************************************
string TTTrack::getQ(void)
//****************************************************************************
{

  string theQ = "N/A";
  if (theCharge > 0) theQ = "+";
  else if (theCharge < 0) theQ = "-";
  else
    {
      std::cout << "E R R O R ! TTTrack::getQ(...) - Invalid value for theCharge \"" << theCharge << "\". EXIT" << std::endl;
      exit(1);
    }

  return theQ;
}
  

//****************************************************************************
void TTTrack::PrintProperties(void)
//****************************************************************************
{
  
  Table info("Index | Pt | Eta | Phi | x0 | y0 | z0 | d0 | Q | Chi2 | DOF | Chi2Red | Stubs (PS) | StubPtCons. | Genuine | Unknown | Combinatoric", "Text");
  info.AddRowColumn(0, auxTools.ToString( theIndex) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Perp(), 3) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Eta() , 3) );  
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Phi() , 3) );
  info.AddRowColumn(0, auxTools.ToString( getX0(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getY0(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getZ0(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getD0(), 3) );  
  info.AddRowColumn(0, getQ());
  info.AddRowColumn(0, auxTools.ToString( getChi2(),3 ));
  info.AddRowColumn(0, auxTools.ToString( getDOF()) );
  info.AddRowColumn(0, auxTools.ToString( getChi2Red(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getNumOfStubs()) + " (" + auxTools.ToString(getNumOfStubsPS()) + ")");
  info.AddRowColumn(0, auxTools.ToString( getStubPtConsistency(), 3) );
  info.AddRowColumn(0, auxTools.ToString( getIsGenuine()) );
  info.AddRowColumn(0, auxTools.ToString( getIsUnknown()) );
  info.AddRowColumn(0, auxTools.ToString( getIsCombinatoric()) );
  info.Print();

  return;
}



//****************************************************************************
void TTTrack::PrintAllProperties(void)
//****************************************************************************
{
  
  Table info("Pt | Eta | Phi | x0 | y0 | z0 | d0 | Q | Chi2 | DOF | Chi2Red | Stubs (PS) | StubPtCons. | Genuine | Unknown | Comb. | Loose | Fake", "Text");  
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Perp(), 4) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Eta() , 4) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Phi() , 4) );
  info.AddRowColumn(0, auxTools.ToString( getX0(), 4) );
  info.AddRowColumn(0, auxTools.ToString( getY0(), 4) );
  info.AddRowColumn(0, auxTools.ToString( getZ0(), 4) );
  info.AddRowColumn(0, auxTools.ToString( getD0(), 4) );
  string theQ = "-";
  if(theCharge > 0) theQ = "+";
  info.AddRowColumn(0, theQ);
  info.AddRowColumn(0, auxTools.ToString( getChi2())    );
  info.AddRowColumn(0, auxTools.ToString( getDOF())     );
  info.AddRowColumn(0, auxTools.ToString( getChi2Red()) );
  info.AddRowColumn(0, auxTools.ToString( getNumOfStubs()) + " (" + auxTools.ToString(getNumOfStubsPS()) + ")");
  info.AddRowColumn(0, auxTools.ToString( getStubPtConsistency()) );
  info.AddRowColumn(0, auxTools.ToString( getIsGenuine()) );
  info.AddRowColumn(0, auxTools.ToString( getIsUnknown()) );
  info.AddRowColumn(0, auxTools.ToString( getIsCombinatoric()) );
  info.AddRowColumn(0, auxTools.ToString( getIsLoose()) );
  info.AddRowColumn(0, auxTools.ToString( getIsFake()) );
  info.Print();

  return;
}


//****************************************************************************
TLorentzVector TTTrack::p4(double mass)
//****************************************************************************
{
  TLorentzVector p4;
  p4.SetPtEtaPhiM(getPt(), getEta(), getPhi(), mass);
  return p4;
}


#endif

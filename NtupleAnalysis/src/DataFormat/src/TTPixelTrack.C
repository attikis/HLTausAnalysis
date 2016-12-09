#ifndef TTPixelTrack_cxx
#define TTPixelTrack_cxx

// User 
#include "../interface/TTPixelTrack.h"

//****************************************************************************
TTPixelTrack::TTPixelTrack()
//****************************************************************************
{

  _InitVars();

}


//****************************************************************************
TTPixelTrack::~TTPixelTrack()
//****************************************************************************
{

}


//****************************************************************************
TTPixelTrack::TTPixelTrack(unsigned short aIndex,
			   TTTrack aTTTrack,
			   TVector3 aMomentum,
			   ROOT::Math::XYZVector aPOCA,
			   double aRInv,
			   double aChi2,
			   double aSigmaRInv, 
			   double aSigmaPhi0,
			   double aSigmaD0,
			   double aSigmaT,
			   double aSigmaZ0,
			   vector<ROOT::Math::XYZVector> pixHits,
			   vector<ROOT::Math::XYZVector> candidatePixHits)
//****************************************************************************
{

  _InitVars();
  theIndex              = aIndex;
  theTTTrack            = aTTTrack;
  theMomentum           = aMomentum;
  thePOCA               = aPOCA;
  theRInv               = aRInv;
  theChi2               = aChi2;
  theSigmaD0            = aSigmaPhi0;
  theSigmaPhi0          = aSigmaD0;
  theSigmaRInv          = aSigmaRInv;
  theSigmaT             = aSigmaT;
  theSigmaZ0            = aSigmaZ0; 
  thePixelHits          = pixHits;
  theCandidatePixelHits = candidatePixHits;
  //
  theCharge             = getCharge();
  theD0                 = getD0();
  theZ0                 = getZ0();
  thenhit               = (int) pixHits.size();
  thencandidatehit      = (int) candidatePixHits.size();
  _FillAuxPixelHitVariables();

  if (0) PrintProperties();
}


//****************************************************************************
void TTPixelTrack::_InitVars(void)
//****************************************************************************
{

  theIndex     = 0.0;
  theMomentum.SetXYZ(0.0, 0.0, 0.0);
  theRInv      = 0.0;
  thePOCA.SetXYZ(0.0, 0.0, 0.0);
  theZ0        = 0.0;
  theD0        = 0.0;
  theChi2      = 0.0;
  theChi2Red   = 0.0;
  theSigmaRInv = 0.0;
  theSigmaPhi0 = 0.0;
  theSigmaD0   = 0.0;
  theSigmaT    = 0.0;
  theSigmaZ0   = 0.0;
  thenhit      = 0;
  thencandidatehit = 0;
  thePixelHits.clear();
  theCandidatePixelHits.clear();
  pixHits_R.clear();
  pixHits_Z.clear();
  pixHits_Phi.clear();
  pixHits_Type.clear();
  candPixHits_R.clear();
  candPixHits_Z.clear();
  candPixHits_Phi.clear();
  candPixHits_Type.clear();

  return;
}


//****************************************************************************
void TTPixelTrack::init(TVector3 aMomentum,
			ROOT::Math::XYZVector aPOCA,
                        double aRInv,
                        double aChi2,
                        int    anhit,
                        double sigmarinv,
                        double sigmad0,
                        double sigmaphi0,
                        double sigmat,
                        double sigmaz0,
                        std::vector<ROOT::Math::XYZVector> pixHits,
                        std::vector<ROOT::Math::XYZVector> candidatePixHits)
//****************************************************************************
{

  // theL1Track=aL1Track;
  _InitVars();
  theCandidatePixelHits = candidatePixHits;
  theCharge             = getCharge();
  theChi2               = aChi2;
  theD0                 = getD0();
  theMomentum           = aMomentum;
  thePOCA               = aPOCA;
  thePixelHits          = pixHits;
  theRInv               = aRInv;
  theSigmaD0            = sigmaphi0;
  theSigmaPhi0          = sigmad0;
  theSigmaRInv          = sigmarinv;
  theSigmaT             = sigmat;
  theSigmaZ0            = sigmaz0;
  theZ0                 = getZ0();
  thencandidatehit      = (int) candidatePixHits.size();
  thenhit               = anhit;
  _FillAuxPixelHitVariables();
  
  return;
}


//****************************************************************************
double TTPixelTrack::getD0(void)
//****************************************************************************
{
  theD0 = -thePOCA.X() * sin( theMomentum.Phi() ) + thePOCA.Y() * cos(theMomentum.Phi() );
  return theD0;
}


//****************************************************************************
int TTPixelTrack::getCharge(void)
//****************************************************************************
{
  if (theRInv < 0.0) theCharge = -1;
  else if (theRInv > 0.0) theCharge = +1;
  else{
    std::cout << "E R R O R ! TTPixelTrack::getCharge(...) - Invalid value for theRInv \"" << theRInv << "\". EXIT" << std::endl;
    exit(1);
  }

  return theCharge;
}


//****************************************************************************
string TTPixelTrack::getQ(void)
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
int TTPixelTrack::getPixelHitType(ROOT::Math::XYZVector pixHit)
//****************************************************************************
{
  
  int type = 4;
  if ( pixHit.Rho() < 12.0 ) type=3;
  if ( pixHit.Rho() <  8.0 ) type=2;
  if ( pixHit.Rho() <  5.0 ) type=1;
  
  if ( fabs( pixHit.Z() ) > 28.0) type=-1;
  if ( fabs( pixHit.Z() ) > 35.0) type=-2;
  if ( fabs( pixHit.Z() ) > 45.0) type=-3;

  return type;
 
}


//****************************************************************************
int TTPixelTrack::getPixelHitsPattern(void)
//****************************************************************************
{
  
  int nHits = (int) thePixelHits.size();  
  thePixelHitsPattern=0;
  
  for(int i =0; i< nHits; i++){

    int hitType = getPixelHitType(thePixelHits.at(i));
    if (hitType > 0) thePixelHitsPattern += pow(2, abs(hitType)-1);
    else thePixelHitsPattern += pow(2, abs(hitType)-1+4);
  }

  return thePixelHitsPattern;

}


//****************************************************************************
int TTPixelTrack::getCandidatePixelHitsPattern(void)
//****************************************************************************
{
  
  int nHits = (int) theCandidatePixelHits.size();  
  theCandidatePixelHitsPattern=0;
  
  for(int i =0; i< nHits; i++)
    {
      int hitType = getPixelHitType(theCandidatePixelHits.at(i));
      if (hitType > 0) theCandidatePixelHitsPattern += pow(2, abs(hitType)-1);
      else theCandidatePixelHitsPattern += pow(2, abs(hitType)-1+4);
    }

  return theCandidatePixelHitsPattern;

}


//****************************************************************************
void TTPixelTrack::_FillAuxPixelHitVariables(void)
//****************************************************************************
{

  // Pixel Hits
  for(int i = 0; i < (int) thePixelHits.size(); i++)
    {
      pixHits_R.push_back( thePixelHits.at(i).Rho()   );
      pixHits_Z.push_back( thePixelHits.at(i).Z()     );
      pixHits_Phi.push_back( thePixelHits.at(i).Phi() );
      pixHits_Type.push_back( getPixelHitType(thePixelHits.at(i) ) );
    }

  // Candidate Pixel Hits
  for(int i = 0; i < (int) theCandidatePixelHits.size(); i++)
    {
      candPixHits_R.push_back(    theCandidatePixelHits.at(i).Rho()  );
      candPixHits_Z.push_back(    theCandidatePixelHits.at(i).Z()    );
      candPixHits_Phi.push_back(  theCandidatePixelHits.at(i).Phi()  );
      candPixHits_Type.push_back( getPixelHitType( theCandidatePixelHits.at(i) ) );
    }
  
  return;
}


//****************************************************************************
void TTPixelTrack::PrintProperties(void)
//****************************************************************************
{
  
  Table info("Index | Pt | Eta | Phi | z0 | d0 | Q | Chi2 | RedChi2 | Hits | Hit-Pattern | Hit-Type | Hit-R | Hit-Z", "Text");
  info.AddRowColumn(0, auxTools.ToString( index()) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Perp(), 3) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Eta() , 3) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Phi() , 3) );
  info.AddRowColumn(0, auxTools.ToString( theZ0, 3) );
  info.AddRowColumn(0, auxTools.ToString( theD0, 3) );
  info.AddRowColumn(0, getQ() );
  info.AddRowColumn(0, auxTools.ToString( theChi2, 3 ) );
  // info.AddRowColumn(0, auxTools.ToString( theChi2Red, 3) );
  info.AddRowColumn(0, auxTools.ToString( thenhit) + " (" + auxTools.ToString(candPixHits_Type.size()) + ")" );
  int pixHits_Pattern = getPixelHitsPattern();
  info.AddRowColumn(0, auxTools.ToString( pixHits_Pattern) );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixHits_Type) );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixHits_R)    );
  info.Print();

  return;
}


//****************************************************************************
void TTPixelTrack::PrintAllProperties(void)
//****************************************************************************
{
  
  std::vector<double> pixHits_R;
  std::vector<double> pixHits_Z;
  std::vector<double> pixHits_Phi;
  for(int i = 0; i < (int) thePixelHits.size(); i++)
    {
      pixHits_R.push_back( thePixelHits.at(i).Rho() ); 
      pixHits_Z.push_back( thePixelHits.at(i).Z() );
      pixHits_Phi.push_back( thePixelHits.at(i).Phi() );
    }

  Table info("Pt | Eta | Phi | RInv | x0 | y0 | z0 | chi2 | redChi2 | err(RInv) | err(Phi0) | err(D0) | err(T) | err(Z0) | NHits | Hit R | Hit Z | Hit Phi", "Text");
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Perp() ) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Eta() ) );
  info.AddRowColumn(0, auxTools.ToString( theMomentum.Phi() ) );
  info.AddRowColumn(0, auxTools.ToString( theRInv ) );
  info.AddRowColumn(0, auxTools.ToString( thePOCA.X() ) );
  info.AddRowColumn(0, auxTools.ToString( thePOCA.Y() ) );
  info.AddRowColumn(0, auxTools.ToString( thePOCA.Z() ) );
  info.AddRowColumn(0, auxTools.ToString( theChi2 ) );
  info.AddRowColumn(0, auxTools.ToString( theChi2Red ) );
  info.AddRowColumn(0, auxTools.ToString( theSigmaRInv ) );
  info.AddRowColumn(0, auxTools.ToString( theSigmaPhi0 ) );
  info.AddRowColumn(0, auxTools.ToString( theSigmaD0 ) );
  info.AddRowColumn(0, auxTools.ToString( theSigmaT) );
  info.AddRowColumn(0, auxTools.ToString( theSigmaZ0) );
  info.AddRowColumn(0, auxTools.ToString( thenhit) );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixHits_R)    );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixHits_Z)    );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixHits_Phi)  );
  info.Print();

  return;
}

//****************************************************************************
TLorentzVector TTPixelTrack::p4(double mass)
//****************************************************************************
{
  TLorentzVector p4;
  p4.SetPtEtaPhiM(getPt(), getEta(), getPhi(), mass);
  return p4;
}


#endif

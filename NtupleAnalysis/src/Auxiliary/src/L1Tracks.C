#ifndef L1Tracks_cxx
#define L1Tracks_cxx

// System
#include <iostream>
#include <algorithm>
#include <cmath>

// User
#include "../interface/constants.h"
#include "../interface/L1Tracks.h"

#ifdef USING_MC
#include "../../Framework/interface/TreeDefinitionGenP.h"
#endif

// ROOT
#include <TLorentzVector.h>


//#define DEBUG

//****************************************************************************
bool L1Tracks::SelectTracks(int iTrack, 
			    const string selectionType,
			    const int nFitParams)
//****************************************************************************
{  // Get track variables
  double tk_pt                 = t->L1Tks_Pt->at(iTrack);
  double tk_eta                = t->L1Tks_Eta->at(iTrack);
  double tk_pocaz              = t->L1Tks_POCAz->at(iTrack);
  double tk_chiSq              = t->L1Tks_ChiSquared->at(iTrack);
  double tk_redChiSq           = GetRedChiSquared(iTrack, nFitParams);
  // double tk_redChiSq           = t->L1Tks_RedChiSquared->at(iTrack); // (tk_chiSq / d.o.f.   tk_dof = 2*tk_nStubs - NFitParams)
  double tk_stubPtConsist      = GetTrackStubPtConsistency(iTrack);
  // double tk_d0Abs              = fabs(GetTrackD0(iTrack));
  unsigned int tk_nStubs       = GetNumOfStubs(iTrack);
  unsigned int tk_nPSStubs     = GetNumOfPSStubs(iTrack);
  unsigned int tk_nEndcapStubs = GetNumOfEndcapStubs(iTrack);
  double tk_chiSq_Cut          = 0.0;
  double tk_redChiSq_Cut       = 0.0;
  bool bPass_pt                = false;
  bool bPass_eta               = false;
  bool bPass_pocaz             = false;
  bool bPass_chiSq             = false;
  bool bPass_redChiSq          = false;
  bool bPass_stubPtConsist     = false;
  bool bPass_d0Abs             = false;
  bool bPass_nStubs            = false;
  bool bPass_nPSStubs          = false;
  
  // Track Selection Working-Points
  if (selectionType.compare("TightWP") == 0) 
    {
    
      tk_chiSq_Cut        = 40.0;
      tk_redChiSq_Cut     = 10.0;
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = (tk_chiSq <= tk_chiSq_Cut);
      bPass_redChiSq      = (tk_redChiSq <= tk_redChiSq_Cut);
      bPass_stubPtConsist = (tk_stubPtConsist <= 5.0);
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 5);
      bPass_nPSStubs      = (tk_nPSStubs >= 3 && tk_nPSStubs <= 7);

    }
  else if (selectionType.compare("MediumWP") == 0) 
    {
    
      tk_redChiSq_Cut     = 10.0;
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = true;
      bPass_redChiSq      = (tk_redChiSq <= tk_redChiSq_Cut);
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 5);
      bPass_nPSStubs      = (tk_nPSStubs >= 3 && tk_nPSStubs <= 7);

    }
  else if (selectionType.compare("LooseWP") == 0) 
    {
      
      if (tk_nEndcapStubs > 0) tk_redChiSq_Cut = 10.0;
      else tk_redChiSq_Cut = 8.0;
      
      bPass_pt            = (tk_pt         >=   2.0);
      bPass_eta           = (abs(tk_eta)   <=   2.5);
      bPass_pocaz         = (abs(tk_pocaz) <=  30.0);
      bPass_chiSq         = true;
      bPass_redChiSq      = (tk_redChiSq   <=  tk_redChiSq_Cut);
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs     >=   4.0);
      bPass_nPSStubs      = true;

    }
  else if (selectionType.compare("vLooseWP") == 0) 
    {
      
      tk_chiSq_Cut        = 200.0;
      bPass_pt            = (tk_pt         >=   2.0);
      bPass_eta           = (abs(tk_eta)   <=   3.0);
      bPass_pocaz         = (abs(tk_pocaz) <=  30.0);
      bPass_chiSq         = (tk_chiSq <=  tk_chiSq_Cut);
      bPass_redChiSq      = true;
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 4.0);
      bPass_nPSStubs      = true;

    }
  else if (selectionType.compare("NoWP") == 0) 
    {
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = true;
      bPass_redChiSq      = true;
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = true;
      bPass_nPSStubs      = true;
      
    }
  else 
  {
    cout << "W A R N I N G ! L1Tracks::SelectTracks(...) - Selection string \"" << selectionType  << "\" is not defined. EXIT"; 
    exit(1);
  }

  return bPass_pt * bPass_eta * bPass_pocaz * bPass_chiSq * bPass_redChiSq * bPass_stubPtConsist * bPass_d0Abs * bPass_nStubs * bPass_nPSStubs;
}



//****************************************************************************
bool L1Tracks::SelectPixTracks(int iTrack, 
			       const string selectionType)
//****************************************************************************
{

  // Get track variables
  const int iTTTrack           =  t->L1PixTks_TTTrackIndex->at(iTrack);
  if (iTTTrack < 0) return false;
  double tk_pt                 = t->L1PixTks_Pt->at(iTrack);
  double tk_eta                = t->L1PixTks_Eta->at(iTrack);
  double tk_pocaz              = t->L1PixTks_POCAz->at(iTrack);
  double tk_chiSq              = t->L1PixTks_ChiSquared->at(iTrack);
  double tk_redChiSq           = GetPixTrackRedChiSq(iTrack);
  double tk_stubPtConsist      = GetTrackStubPtConsistency(iTTTrack);
  // double tk_d0Abs              = fabs(GetPixTrackD0(iTrack));
  unsigned int tk_nStubs       = GetNumOfStubs(iTTTrack);
  unsigned int tk_nPSStubs     = GetNumOfPSStubs(iTTTrack);
  unsigned int tk_nEndcapStubs = GetNumOfEndcapStubs(iTTTrack);
  double tk_chiSq_Cut          = 0.0;
  double tk_redChiSq_Cut       = 0.0;
  bool bPass_pt                = true;
  bool bPass_eta               = true;
  bool bPass_pocaz             = true;
  bool bPass_chiSq             = true;
  bool bPass_redChiSq          = true;
  bool bPass_stubPtConsist     = true;
  bool bPass_d0Abs             = true;
  bool bPass_nStubs            = true;
  bool bPass_nPSStubs          = true;
  
  // Track Selection Working-Points
  if (selectionType.compare("TightWP") == 0) 
    {
    
      tk_chiSq_Cut        = 40.0;
      tk_redChiSq_Cut     = 10.0;
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = (tk_chiSq <= tk_chiSq_Cut);
      bPass_redChiSq      = (tk_redChiSq <= tk_redChiSq_Cut);
      bPass_stubPtConsist = (tk_stubPtConsist <= 5.0);
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 5);
      bPass_nPSStubs      = (tk_nPSStubs >= 3 && tk_nPSStubs <= 7);

    }
  else if (selectionType.compare("MediumWP") == 0) 
    {
    
      tk_redChiSq_Cut     = 10.0;
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = true;
      bPass_redChiSq      = (tk_redChiSq <= tk_redChiSq_Cut);
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 5);
      bPass_nPSStubs      = (tk_nPSStubs >= 3 && tk_nPSStubs <= 7);

    }
  else if (selectionType.compare("LooseWP") == 0) 
    {
      
      if (tk_nEndcapStubs > 0) tk_redChiSq_Cut = 10.0;
      else tk_redChiSq_Cut = 8.0;
      
      bPass_pt            = (tk_pt         >=   2.0);
      bPass_eta           = (abs(tk_eta)   <=   2.5);
      bPass_pocaz         = (abs(tk_pocaz) <=  30.0);
      bPass_chiSq         = true;
      bPass_redChiSq      = (tk_redChiSq   <=  tk_redChiSq_Cut);
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs     >=   4.0);
      bPass_nPSStubs      = true;

    }
  else if (selectionType.compare("vLooseWP") == 0) 
    {
      
      tk_chiSq_Cut        = 200.0;
      bPass_pt            = (tk_pt         >=   2.0);
      bPass_eta           = (abs(tk_eta)   <=   3.0);
      bPass_pocaz         = (abs(tk_pocaz) <=  30.0);
      bPass_chiSq         = (tk_chiSq <=  tk_chiSq_Cut);
      bPass_redChiSq      = true;
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 4.0);
      bPass_nPSStubs      = true;

    }
  else if (selectionType.compare("NoWP") == 0) 
    {
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = true;
      bPass_redChiSq      = true;
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = true;
      bPass_nPSStubs      = true;
      
    }
  else 
  {
    cout << "W A R N I N G ! L1Tracks::SelectPixTracks(...) - Selection string \"" << selectionType  << "\" is not defined. EXIT"; 
    exit(1);
  }

  return bPass_pt * bPass_eta * bPass_pocaz * bPass_chiSq * bPass_redChiSq * bPass_stubPtConsist * bPass_d0Abs * bPass_nStubs * bPass_nPSStubs;
}



//****************************************************************************
bool L1Tracks::SelectPixTracksRefit(TTPixelTrack t,
				    int iTrack,
				    const string selectionType)
//****************************************************************************
{

  // The Refitted TTPixelTrack
  double tk_pt                 = t.getMomentum().Perp();
  double tk_eta                = t.getMomentum().Eta();
  double tk_pocaz              = t.getPOCA().z();
  double tk_chiSq              = t.getChi2();
  double tk_redChiSq           = t.getChi2Red();

  // The corresponding TTTrack
  double tk_stubPtConsist      = GetTrackStubPtConsistency(iTrack);
  /// double tk_d0Abs              = fabs(GetTrackD0(iTrack));
  unsigned int tk_nStubs       = GetNumOfStubs(iTrack);
  unsigned int tk_nPSStubs     = GetNumOfPSStubs(iTrack);
  unsigned int tk_nEndcapStubs = GetNumOfEndcapStubs(iTrack);
  double tk_chiSq_Cut          = 0.0;
  double tk_redChiSq_Cut       = 0.0;
  bool bPass_pt                = true;
  bool bPass_eta               = true;
  bool bPass_pocaz             = true;
  bool bPass_chiSq             = true;
  bool bPass_redChiSq          = true;
  bool bPass_stubPtConsist     = true;
  bool bPass_d0Abs             = true;
  bool bPass_nStubs            = true;
  bool bPass_nPSStubs          = true;
  
  // Track Selection Working-Points
  if (selectionType.compare("TightWP") == 0) 
    {
    
      tk_chiSq_Cut        = 40.0;
      tk_redChiSq_Cut     = 10.0;
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = (tk_chiSq <= tk_chiSq_Cut);
      bPass_redChiSq      = (tk_redChiSq <= tk_redChiSq_Cut);
      bPass_stubPtConsist = (tk_stubPtConsist <= 5.0);
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 5);
      bPass_nPSStubs      = (tk_nPSStubs >= 3 && tk_nPSStubs <= 7);

    }
  else if (selectionType.compare("MediumWP") == 0) 
    {
    
      tk_redChiSq_Cut     = 10.0;
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = true;
      bPass_redChiSq      = (tk_redChiSq <= tk_redChiSq_Cut);
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 5);
      bPass_nPSStubs      = (tk_nPSStubs >= 3 && tk_nPSStubs <= 7);

    }
  else if (selectionType.compare("LooseWP") == 0) 
    {
      
      if (tk_nEndcapStubs > 0) tk_redChiSq_Cut = 10.0;
      else tk_redChiSq_Cut = 8.0;
      
      bPass_pt            = (tk_pt         >=   2.0);
      bPass_eta           = (abs(tk_eta)   <=   2.5);
      bPass_pocaz         = (abs(tk_pocaz) <=  30.0);
      bPass_chiSq         = true;
      bPass_redChiSq      = (tk_redChiSq   <=  tk_redChiSq_Cut);
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs     >=   4.0);
      bPass_nPSStubs      = true;

    }
  else if (selectionType.compare("vLooseWP") == 0) 
    {
      
      tk_chiSq_Cut        = 200.0;
      bPass_pt            = (tk_pt         >=   2.0);
      bPass_eta           = (abs(tk_eta)   <=   3.0);
      bPass_pocaz         = (abs(tk_pocaz) <=  30.0);
      bPass_chiSq         = (tk_chiSq <=  tk_chiSq_Cut);
      bPass_redChiSq      = true;
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = (tk_nStubs >= 4.0);
      bPass_nPSStubs      = true;

    }
  else if (selectionType.compare("NoWP") == 0) 
    {
      bPass_pt            = true;
      bPass_eta           = true;
      bPass_pocaz         = true;
      bPass_chiSq         = true;
      bPass_redChiSq      = true;
      bPass_stubPtConsist = true;
      bPass_d0Abs         = true;
      bPass_nStubs        = true;
      bPass_nPSStubs      = true;
      
    }
  else 
  {
    cout << "W A R N I N G ! L1Tracks::SelectPixTracksRefit(...) - Selection string \"" << selectionType  << "\" is not defined. EXIT"; 
    exit(1);
  }

  return bPass_pt * bPass_eta * bPass_pocaz * bPass_chiSq * bPass_redChiSq * bPass_stubPtConsist * bPass_d0Abs * bPass_nStubs * bPass_nPSStubs;
}


//****************************************************************************
bool L1Tracks::HasMustPixelHitInLayerAndDisk(int iTrack,
					 std::vector<int> pixTk_Barrel_Type,
					 std::vector<int> pixTk_Endcap_Type,
					 double endcapEtaBoundary)
//****************************************************************************
{

  bool bFailedPixelHitType = false;
  bool bSuccess            = false; 
  double pixTk_Eta         = t->L1PixTks_Eta->at(iTrack);
  vector< int > pixTk_Type = t->L1PixTks_PixHits_Type->at(iTrack);

  if ( fabs(pixTk_Eta) <= endcapEtaBoundary ){
    for(int i = 0; i < (int) pixTk_Barrel_Type.size(); i++){	 
      int pixTk_HitType = pixTk_Barrel_Type.at(i);
      if( std::find(pixTk_Type.begin(), pixTk_Type.end(), pixTk_HitType) == pixTk_Type.end() ) bFailedPixelHitType = true;
    }
  }
  else{
    for(int i = 0; i < (int) pixTk_Endcap_Type.size(); i++){	 
      int pixTk_HitType = pixTk_Endcap_Type.at(i);
      if( std::find(pixTk_Type.begin(), pixTk_Type.end(), pixTk_HitType) == pixTk_Type.end() ) bFailedPixelHitType = true;
    }
  }

  bSuccess = !bFailedPixelHitType;
  return bSuccess;
}


//****************************************************************************
bool L1Tracks::HasMustPixelHitPattern(int pixTk_PixHits_Pattern,
				      vector<int> pixTk_PixHits_AllowedPatterns)
//****************************************************************************
{
  bool bSuccess = false;
  int numOfAcceptablePatterns = int(pixTk_PixHits_AllowedPatterns.size());
  if (numOfAcceptablePatterns < 1) return true;
  
  // For-Loop: All allowed hit patterns
  for (int i=0; i < numOfAcceptablePatterns; i++){
    
    if( pixTk_PixHits_Pattern != pixTk_PixHits_AllowedPatterns.at(i) ) continue;
    else{
      bSuccess = true;
      break;
    }
  }

  return bSuccess;
}



//****************************************************************************
void L1Tracks::GetPixTrackSharedHits(int iTrack,
				     vector<int>    &sharedPixHits_PixTkRefIndex,
				     vector<double> &sharedPixHits_X,
				     vector<double> &sharedPixHits_Y,
				     vector<double> &sharedPixHits_Z,
				     vector<double> &sharedPixHits_R,
				     vector<double> &sharedPixHits_Phi,
				     vector<int>    &sharedPixHits_Type,
				     vector<TVector3> pixHits_XYZ,
				     vector<int> pixHits_TTPixelTrackIndex)
//****************************************************************************
{
  if (iTrack < 0) return;

  // Get pixel hits for the TTPixelTrack with index iTrack
  vector<double> pixHits_X    = t->L1PixTks_PixHits_X->at(iTrack);
  vector<double> pixHits_Y    = t->L1PixTks_PixHits_Y->at(iTrack);
  vector<double> pixHits_Z    = t->L1PixTks_PixHits_Z->at(iTrack);
  vector<double> pixHits_R    = t->L1PixTks_PixHits_R->at(iTrack);
  vector<double> pixHits_Phi  = t->L1PixTks_PixHits_Phi->at(iTrack);
  vector< int >  pixHits_Type = t->L1PixTks_PixHits_Type->at(iTrack);

  // For-loop: All Pixel Hits
  for (int i = 0; i < (int) pixHits_XYZ.size(); i++){

    // Create 3D space point for current pixel hit
    TVector3 globalPixHit_XYZ = pixHits_XYZ.at(i);
    int refTrack = pixHits_TTPixelTrackIndex.at(i);
    bool bIsSharedHit = false;
    
    // Skip self
    if (refTrack == iTrack) continue;
    
    // For-loop: Pixel Hits of TTPixelTrack
    for (int j = 0; j < (int) pixHits_X.size(); j++){

      TVector3 tkPixHit_XYZ(pixHits_X.at(j), pixHits_Y.at(j), pixHits_Z.at(j) );

      bool bSameX = ( globalPixHit_XYZ.X() == tkPixHit_XYZ.X() );
      bool bSameY = ( globalPixHit_XYZ.Y() == tkPixHit_XYZ.Y() );
      bool bSameZ = ( globalPixHit_XYZ.Z() == tkPixHit_XYZ.Z() );
      bIsSharedHit = (bSameX * bSameY * bSameZ);

      // If no match found skip this hit
      if (!bIsSharedHit) continue;

      // Else save shared hits into their containers 
      sharedPixHits_PixTkRefIndex.push_back( refTrack );
      sharedPixHits_X.push_back( pixHits_X.at(j) );
      sharedPixHits_Y.push_back( pixHits_Y.at(j) );
      sharedPixHits_Z.push_back( pixHits_Z.at(j) );
      sharedPixHits_R.push_back( pixHits_R.at(j) );
      sharedPixHits_Phi.push_back( pixHits_Phi.at(j) );
      sharedPixHits_Type.push_back( pixHits_Type.at(j) );      
      
    }//for (int i = 0; i <  pixHits_X.size(); i++){

    
  } //for (int i = 0; i <  (int) pixHits.size(); i++){

  return;
  
}


//****************************************************************************
void L1Tracks::GetPixTrackAllHits(vector<TVector3> &pixHits_all,
				  vector<int> &pixHits_TTPixelTrackIndex)
//****************************************************************************
{

  // For-loop: L1PixTks
  for (int j = 0; j < (int) t->L1PixTks_Pt->size(); j++) {

    
    // Get all hits for given pixel track    
    vector<double> pixHits_all_X = t->L1PixTks_PixHits_X->at(j);
    vector<double> pixHits_all_Y = t->L1PixTks_PixHits_Y->at(j);
    vector<double> pixHits_all_Z = t->L1PixTks_PixHits_Z->at(j);
    
    // For-loop: Pixel Hits
    for (int k = 0; k <  (int) pixHits_all_X.size(); k++){
    
      // Create an XYZ vector for the pixel hit with x=pixHits_all_X, y=pixHits_all_Y , z=pixHits_all_Z
      TVector3 xyz( pixHits_all_X.at(k), pixHits_all_Y.at(k), pixHits_all_Z.at(k) );
      pixHits_all.push_back(xyz);

      // Save the index of the pixel track that has these hits
      pixHits_TTPixelTrackIndex.push_back(j);
      
    }// For-loop: Pixel Hits
    
  }// For-loop: L1PixTks

  return;
}


//****************************************************************************
bool L1Tracks::HasMustPixelHitInLayeOrDisk(int iTrack,
					   std::vector<int> pixTk_HitTypes)
//****************************************************************************
{

  bool bSuccess            = false;
  double pixTk_Eta         = t->L1PixTks_Eta->at(iTrack);
  vector< int > pixTk_Type = t->L1PixTks_PixHits_Type->at(iTrack);

  if (pixTk_HitTypes.size() < 1) return true;
  
  for(int i = 0; i < (int) pixTk_HitTypes.size(); i++){	 
    int pixTk_HitType = pixTk_HitTypes.at(i);

    if( std::find(pixTk_Type.begin(), pixTk_Type.end(), pixTk_HitType) != pixTk_Type.end() ){
      bSuccess = true;
      break;
    }
  }

  return bSuccess;
}


//****************************************************************************
double L1Tracks::GetTrackCharge(const int iTrack)
//****************************************************************************
{

  double tk_RInv    = t->L1Tks_RInv->at(iTrack);
  double tk_Charge  = +1.0;
  if (tk_RInv < 0.0) tk_Charge = -1.0;

  return tk_Charge;
}


//****************************************************************************
double L1Tracks::GetPixTrackCharge(const int iTrack)
//****************************************************************************
{

  double tk_RInv    = t->L1PixTks_RInv->at(iTrack);
  double tk_Charge  = +1.0;
  if (tk_RInv < 0.0) tk_Charge = -1.0;

  return tk_Charge;
}



//****************************************************************************
double L1Tracks::GetTrackPt(const int iTrack)
//****************************************************************************
{

  // pT [GeV/c] = 0.3 * B [T] * rho [m]
  double tk_RInv    = t->L1Tks_RInv->at(iTrack);
  double tk_Rho     = (1.0/tk_RInv); //cm
  double cms_BField = 3.8; // Tesla
  double tk_Pt      = (0.3) * (cms_BField) * (tk_Rho * 0.01); // GeV/c
  // double tk_Pt      = t->L1Tks_Pt->at(iTrack);
  
  return tk_Pt;
}



//****************************************************************************
double L1Tracks::GetPixTrackPt(const int iTrack)
//****************************************************************************
{

  // pT [GeV/c] = 0.3 * B [T] * rho [m]
  double tk_RInv    = t->L1PixTks_RInv->at(iTrack);
  double tk_Rho     = (1.0/tk_RInv); //cm
  double cms_BField = 3.8; // Tesla
  double tk_Pt      = (0.3) * (cms_BField) * (tk_Rho * 0.01); // GeV/c
  // double tk_Pt      = t->L1Tks_Pt->at(iTrack);
  
  return tk_Pt;
}



//****************************************************************************
double L1Tracks::GetTrackD0(const int iTrack)
//****************************************************************************
{
  
  double tk_Phi   = t->L1Tks_Phi->at(iTrack);
  double tk_POCAx = t->L1Tks_POCAx->at(iTrack);
  double tk_POCAy = t->L1Tks_POCAy->at(iTrack);
  double tk_d0    = -tk_POCAx * sin(tk_Phi) + tk_POCAy * cos(tk_Phi);

  return tk_d0;
  
}



//****************************************************************************
double L1Tracks::GetTrackD0Sign(const int iTrack)
//****************************************************************************
{

  // The track d0 sign has the same sign as the curvature (rho)
  // [and hence charge (Q) since Sgn(Q) = Sgn(rho)] if the z-axis is OUTSIDE
  // the circle in the x-y plane formed by the track helix projection.
  // In other words, if the track is going clockwise in the x-y plane it has a
  // positive d0 sign, otherwise negative. The latter is irrespective of charge.
  // One way to determine the track direction (clockwise or anticlockwise) is to
  // first take the cross product of the POCA vector with the instantaneous
  // momentum at the POCA and then project along the unit vector of the viewing direction.
  // If you are looking AT the x-y plane the we are looking along the negative z-axis.
  // Therefore we need to project along (0, 0, -1)
  // double tk_RInv = t->L1Tks_RInv->at(iTrack);

  // Q. How to determine if a vector is moving clockwise or anti-clockwise to another
  // A. The question only makes sense in 3-D since the notion only makes sense
  // once you have established a viewing direction. If V is the viewing direction vector
  // and A and B are the individual vectors, then what you look at is the sign of
  // (A x B) o V 
  // where A x B is the vector cross product and "o" is the dot product.
  // If your direction vector is looking down on the Cartesian plane then V = (0, 0, -1) and
  // then a positive sign indicates a clockwise move. 
  // To show this, let A = (1, 0, 0) and let B = (0, 1, 0); that is, A is along the x-axis and
  // vector B is along the y-axis. Then A x B = (0, 0, 1) and it is pointing up along the Z axis.
  // The dot product (0, 0, 1) o (0, 0, -1) = -1 and, indeed, from the X axis to the Y axis is
  // counter-clockwise when viewed from the +Z side of the plane.
    

  // Get the POCA vector
  double tk_POCAx = t->L1Tks_POCAx->at(iTrack);
  double tk_POCAy = t->L1Tks_POCAy->at(iTrack);
  double tk_POCAz = t->L1Tks_POCAz->at(iTrack);
  TVector3 v_POCA(tk_POCAx, tk_POCAy, tk_POCAz);
  
  // Get the track momentum vector at the POCA
  double tk_px = t->L1Tks_Px->at(iTrack);
  double tk_py = t->L1Tks_Py->at(iTrack);
  double tk_pz = t->L1Tks_Pz->at(iTrack);
  TVector3 v_p0(tk_px, tk_py, tk_pz);
  TVector3 v_p0_unit = v_p0.Unit();
  
  // Check the dot product sign: vec{a} dot vec{b}  = ax*bx + ay*by + az*bz
  // double dot_product = v_POCA.Dot(v_p0_unit);
  TVector3 cross_product = v_POCA.Cross(v_p0_unit);

  // For looking down on the Cartesian plane the viewing direction vector (0, 0, -1) is (0, 0, -1)
  TVector3 v_z(0, 0, -1);
  double dot_product = cross_product.Dot(v_z);  
  
  // If dot product is positive then track is clockwise. Otherwise anticlockwise
  if (dot_product >= 0.0) return +1.0;
  else return -1.0;
      
}


//****************************************************************************
double L1Tracks::GetTrackD0Mag(const int iTrack)
//****************************************************************************
{

  double tk_POCAx   = t->L1Tks_POCAx->at(iTrack);
  double tk_POCAy   = t->L1Tks_POCAy->at(iTrack);
  double tk_d0_mag  = sqrt(tk_POCAx * tk_POCAx + tk_POCAy * tk_POCAy);

  return tk_d0_mag;
  
}


//****************************************************************************
double L1Tracks::GetTrackD0Phi(const int iTrack)
//****************************************************************************
{

  double tk_POCAx    = t->L1Tks_POCAx->at(iTrack);
  double tk_POCAy    = t->L1Tks_POCAy->at(iTrack);
  double d0_phi      = atan2(tk_POCAy, tk_POCAx); // range(-pi, +pi)
  // double d0_phi      = atan(tk_POCAy/tk_POCAx);   // range(-pi/2, +pi/2)
  return d0_phi;
}




//****************************************************************************
double L1Tracks::GetPixTrackD0(const int iTrack)
//****************************************************************************
{

  double pixTk_Phi   = t->L1PixTks_Phi->at(iTrack);
  double pixTk_POCAx = t->L1PixTks_POCAx->at(iTrack);
  double pixTk_POCAy = t->L1PixTks_POCAy->at(iTrack);
  double pixTk_d0    = -pixTk_POCAx * sin(pixTk_Phi) + pixTk_POCAy * cos(pixTk_Phi);

  return pixTk_d0;
}


//****************************************************************************
double L1Tracks::GetPixTrackD0Sig(const int iTrack)
//****************************************************************************
{

  double tk_d0Sig = GetPixTrackD0(iTrack)/t->L1PixTks_SigmaD0->at(iTrack);

  return tk_d0Sig;
}


//****************************************************************************
double L1Tracks::GetTrackPOCAzSig(const int iTrack)
//****************************************************************************
{

  cout << "E R R O R ! L1Tracks::GetTrackPOCAzSig(...) - This is a placeholder. It will only become available once Anders Ryd provides the L1 TTTracks fit error for z0. EXIT" << endl;
  exit(1);
  double tk_POCAzSig = t->L1Tks_POCAz->at(iTrack);

  return tk_POCAzSig;
}


//****************************************************************************
double L1Tracks::GetPixTrackPOCAzSig(const int iTrack)
//****************************************************************************
{

  double tk_POCAzSig = t->L1PixTks_POCAz->at(iTrack) / t->L1PixTks_SigmaZ0->at(iTrack);

  return tk_POCAzSig;
}


//****************************************************************************
double L1Tracks::GetTrackSigmaPt(const int iTrack)
//****************************************************************************
{

  cout << "E R R O R ! L1Tracks::GetTrackPtSig(...) - This is a placeholder. It will only become available once Anders Ryd provides the L1 TTTracks fit error for RInv. EXIT" << endl;
  exit(1);
  double tk_PtSig = t->L1Tks_Pt->at(iTrack);

  return tk_PtSig;
}


//****************************************************************************
double L1Tracks::GetPixTrackSigmaPt(const int iTrack)
//****************************************************************************
{

  double tk_Pt        = t->L1Tks_Pt->at(iTrack);
  double from_cm_to_m = 100;
  double tk_PtSigma   = (tk_Pt * tk_Pt * t->L1PixTks_SigmaRInv->at(iTrack))/(0.299792*3.8*from_cm_to_m);
  //double tk_PtSig     = tk_Pt/tk_PtSigma;

  return tk_PtSigma;
}


//****************************************************************************
double L1Tracks::GetTrackRedChiSq(const int iTrack, 
				  const int nFitParams)
//****************************************************************************
{

  double tk_dof      = 2*GetNumOfStubs(iTrack) - nFitParams;
  double tk_chiSq    = t->L1Tks_ChiSquared->at(iTrack);
  double tk_redChiSq = tk_chiSq/tk_dof;

  return tk_redChiSq;
}


//****************************************************************************
double L1Tracks::GetPixTrackRedChiSq(const int iTrack)
//****************************************************************************
{
  
  // Anders definition as of 20 Feb 2015:   chi^2(L1track fit)+chi^2(pixel fit)]/[2*nstubs+2*npixel-5]
  // Each hit (stub or pixel) provides a x and y measurement and we have 5
  // parameters. So this is reasonable for the number of degrees of freedom.
  // But for the chi^2 I added the L1track fit and the pixel fit, since the chi^2
  // for the pixel fit does not include the contributions for the L1track 
  // fit. (This could be changed if we want.)
  const int nFitParams = 5;
  const int iTTTrack   = t->L1PixTks_TTTrackIndex->at(iTrack);
  // double tk_dof        = 2.0*t->L1PixTks_NPixHits->at(iTrack) - nFitParams; // Original Anders Version (in L1PixelTrackFit.cc)
  // double tk_chiSq      = t->L1PixTks_ChiSquared->at(iTrack);
  double tk_dof        = 2*GetNumOfStubs( t->L1PixTks_TTTrackIndex->at(iTrack) ) + 2*t->L1PixTks_NPixHits->at(iTrack) - nFitParams; // Revised Anders Version
  double tk_chiSq      = t->L1Tks_ChiSquared->at(iTTTrack) + t->L1PixTks_ChiSquared->at(iTrack);
  double tk_redChiSq   = tk_chiSq/tk_dof;


  // Temporarily disable until 100% sure.
  tk_redChiSq   = 0.0;
  
  return tk_redChiSq;
}


//****************************************************************************
double L1Tracks::GetTrackStubPtConsistency(const int iTrack)
//****************************************************************************
{

  // pT consistency = sum_{stubs} ( 1/pT_{stub} - 1/pt_{track})^2 / RMS(1/pT_{stub})^2
  return t->L1Tks_StubPtConsistency->at(iTrack);
}


//****************************************************************************
void L1Tracks::PrintTrackProperties(const int iTrack,
				    const int nFitParams)
//****************************************************************************
{

  if (iTrack < 0) return;
  
  // Sanity check
  const int nTks = t->L1Tks_Pt->size();
  if ( iTrack >= nTks ){
    cout << "E R R O R ! L1Tracks::PrintTrackProperties(...) - "
	 << "the selected track index \"" << iTrack << "\" is invalid. The total number of tracks is \" " << nTks << "\". EXIT";
    exit(1);
  }
 
  // Inform user of settings
  // Table tkProperties("Evt | iTk | Pt | Eta | Phi | z0 | d0 | Q | ChiSq | RedChiSq | Stubs (PS) | B-Stubs (E) | iLayer | iDisk | Genuine | Unknown | Combinatoric", "Text");
  Table tkProperties("Evt | iTk | Pt | Eta | Phi | z0 | d0 | Q | ChiSq | RedChiSq | Stubs (PS) | B-Stubs (E) | iLayer | iDisk | Genuine", "Text");  

  // Get track variables
  int EvtNumber                         = t->EvtNumber;
  double tk_pt                          = t->L1Tks_Pt->at(iTrack);
  double tk_eta                         = t->L1Tks_Eta->at(iTrack);
  double tk_phi                         = t->L1Tks_Phi->at(iTrack);
  double tk_pocaz                       = t->L1Tks_POCAz->at(iTrack);
  int tk_nStubs                         = t->L1Tks_Stubs_iLayer->at(iTrack).size();
  unsigned int tk_nPSStubs              = GetNumOfPSStubs(iTrack);
  unsigned int tk_nBarrelStubs          = GetNumOfBarrelStubs(iTrack);
  unsigned int tk_nEndcapStubs          = GetNumOfEndcapStubs(iTrack);
  double tk_ChiSq                       = t->L1Tks_ChiSquared->at(iTrack);
  double tk_RedChiSq                    = GetRedChiSquared(iTrack, nFitParams); // t->L1Tks_RedChiSquared->at(iTrack);
  double tk_d0                          = GetTrackD0(iTrack);
  double tk_Charge                      = t->L1Tks_Charge->at(iTrack);
  vector< unsigned int > tkStubs_iLayer = t->L1Tks_Stubs_iLayer->at(iTrack); // layer number (for endcap: 999999)
  vector< unsigned int > tkStubs_iDisk  = t->L1Tks_Stubs_iDisk->at(iTrack);  // disk number (for barrel: 999999 )
  // vector< unsigned int > tkStubs_iZ     = t->L1Tks_Stubs_iZ->at(iTrack);     // module number (for endcap: disk number overload)
  // vector< unsigned int > tkStubs_iPhi   = t->L1Tks_Stubs_iPhi->at(iTrack);   // ladder/rod number (for endcap: module number)
  // vector< unsigned int > tkStubs_iSide  = t->L1Tks_Stubs_iSide->at(iTrack);  // 1 or 2 for FrontWheel BackWheel endcaps (for barrel: 999999 )
  // vector< unsigned int > tkStubs_iRing  = t->L1Tks_Stubs_iRing->at(iTrack);  // ring number (for barrel: 999999 )
  // double tk_Sector                      = t->L1Tks_Sector->at(iTrack);
  // double tk_Wedge                       = t->L1Tks_Wedge->at(iTrack);
  // double tk_RInv                        = t->L1Tks_RInv->at(iTrack);
  bool tk_IsGenuine                     = t->L1Tks_IsGenuine->at(iTrack);  
  bool tk_IsUnknown                     = t->L1Tks_IsUnknown->at(iTrack);
  bool tk_IsCombinatoric                = t->L1Tks_IsCombinatoric->at(iTrack);
  string tk_Q = "-";
  if(tk_Charge > 0) tk_Q = "+";

  // Add table rows
  tkProperties.AddRowColumn(0, auxTools.ToString(EvtNumber)       );
  tkProperties.AddRowColumn(0, auxTools.ToString(iTrack) + " (" + auxTools.ToString(nTks-1) + ")" );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_pt   , 4)     );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_eta  , 4)     );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_phi  , 4)     );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_pocaz, 4)     );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_d0   , 4)     );
  tkProperties.AddRowColumn(0, tk_Q                               );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_ChiSq)        );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_RedChiSq)     );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_nStubs) + " (" + auxTools.ToString(tk_nPSStubs) + ")");
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_nBarrelStubs) + " (" + auxTools.ToString(tk_nEndcapStubs) + ")");
  tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(tkStubs_iLayer) );
  // tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(tkStubs_iZ) );
  // tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(tkStubs_iPhi) );
  tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(tkStubs_iDisk ) );
  // tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(tkStubs_iRing ) );
  tkProperties.AddRowColumn(0, auxTools.ToString(tk_IsGenuine) );
  // tkProperties.AddRowColumn(0, auxTools.ToString(tk_IsUnknown) );
  // tkProperties.AddRowColumn(0, auxTools.ToString(tk_IsCombinatoric));
  tkProperties.Print();
  
  return;
}


//****************************************************************************
void L1Tracks::PrintPixTrackProperties(const int iPixTrack, 
				       bool bPrintRefTTTrack)
//****************************************************************************
{
  if (iPixTrack < 0) return;
    
  // Sanity check
  const int nPixTks = t->L1PixTks_Pt->size();
  if ( iPixTrack >= nPixTks ){
    cout << "E R R O R ! L1Tracks::PrintPixTrackProperties(...) - "
	 << "the selected track index \"" << iPixTrack << "\" is invalid. The total number of tracks is \" " << nPixTks << "\". EXIT";
    exit(1);
  }

  // Inform user of settings
  Table tkProperties("Evt | iPixTk | Pt | Eta | Phi | z0 | d0 | Q | ChiSq | RedChiSq | Hits | Type | Hit R | Hit Z | Hit Phi", "Text");
 
  // Get pixel track variables
  int EvtNumber          = t->EvtNumber;
  double pixTk_pt        = t->L1PixTks_Pt->at(iPixTrack);
  double pixTk_eta       = t->L1PixTks_Eta->at(iPixTrack);
  double pixTk_phi       = t->L1PixTks_Phi->at(iPixTrack);
  double pixTk_pocaz     = t->L1PixTks_POCAz->at(iPixTrack);
  double pixTk_ChiSq     = t->L1PixTks_ChiSquared->at(iPixTrack);
  double pixTk_RedChiSq  = GetPixTrackRedChiSq(iPixTrack);
  double pixTk_d0        = GetPixTrackD0(iPixTrack);
  // double pixTk_d0Sig     = GetPixTrackD0Sig(iPixTrack);
  double pixTk_Charge    = t->L1PixTks_Charge->at(iPixTrack);
  double pixTk_SigmaRInv = t->L1PixTks_SigmaRInv->at(iPixTrack);
  double pixTk_SigmaPhi0 = t->L1PixTks_SigmaPhi0->at(iPixTrack);
  double pixTk_SigmaD0   = t->L1PixTks_SigmaD0->at(iPixTrack);
  double pixTk_SigmaT    = t->L1PixTks_SigmaT->at(iPixTrack);
  double pixTk_SigmaZ0   = t->L1PixTks_SigmaZ0->at(iPixTrack);
  const int iTTTrack     = t->L1PixTks_TTTrackIndex->at(iPixTrack);
  // Pixel Hits (best-fit)
  int pixTk_PixHits                    = t->L1PixTks_NPixHits->at(iPixTrack);
  vector<double> pixTk_PixHits_X       = t->L1PixTks_PixHits_X->at(iPixTrack);
  vector<double> pixTk_PixHits_Y       = t->L1PixTks_PixHits_Y->at(iPixTrack);
  vector<double> pixTk_PixHits_Z       = t->L1PixTks_PixHits_Z->at(iPixTrack);
  vector<double> pixTk_PixHits_R       = t->L1PixTks_PixHits_R->at(iPixTrack);
  vector<double> pixTk_PixHits_Phi     = t->L1PixTks_PixHits_Phi->at(iPixTrack);
  vector< int > pixTk_PixHits_Type     = t->L1PixTks_PixHits_Type->at(iPixTrack);
  int pixTk_CandPixHits_N              = t->L1PixTks_CandPixHits_Z->at(iPixTrack).size();
  // vector<double> pixTk_CandPixHits_Z   = t->L1PixTks_CandPixHits_Z->at(iPixTrack);
  // vector<double> pixTk_CandPixHits_R   = t->L1PixTks_CandPixHits_R->at(iPixTrack);
  // vector<double> pixTk_CandPixHits_Phi = t->L1PixTks_CandPixHits_Phi->at(iPixTrack);
  // vector< int > pixTk_CandPixHits_Type = t->L1PixTks_CandPixHits_Type->at(iPixTrack);
  
  // Pixel Hits (candidates)
  int pixTk_CandPixHits = t->L1PixTks_CandPixHits_X->at(iPixTrack).size();
 
  string pixTk_Q = "-";
  if(pixTk_Charge > 0) pixTk_Q = "+";

  // Add table rows
  tkProperties.AddRowColumn(0, auxTools.ToString(EvtNumber)      );
  tkProperties.AddRowColumn(0, auxTools.ToString(iPixTrack) + " (" + auxTools.ToString(nPixTks-1) + ")" );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_pt   , 4) );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_eta  , 4) );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_phi  , 4) );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_pocaz, 4) );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_d0   , 4) );
  tkProperties.AddRowColumn(0, pixTk_Q                           );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_ChiSq)    );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_RedChiSq) );
  tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_PixHits) + " (" + auxTools.ToString(pixTk_CandPixHits_N) + ")" );
  tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_PixHits_Type) );
  tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_PixHits_R)    );
  tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_PixHits_Z)    );
  tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_PixHits_Phi)  );
  // tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_CandPixHits) );
  // tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_CandPixHits_Type) );
  // tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_CandPixHits_R) );
  // tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_CandPixHits_Z) );
  // tkProperties.AddRowColumn(0, auxTools.ConvertIntVectorToString(pixTk_CandPixHits_Phi) );
  // tkProperties.AddRowColumn(0, auxTools.ToString(iTTTrack) );
  
  
  tkProperties.Print();
  if(bPrintRefTTTrack) PrintTrackProperties(iTTTrack);

  return;
}



//****************************************************************************
void L1Tracks::PrintPixTrackRefitProperties(TTPixelTrack tk,
 					    int pixTk_index,
					    bool bPrintRefTTPixelTrack)
//****************************************************************************
{

  
   // Inform user of settings
   Table tkProperties("Evt | Pt | Eta | Phi | z0 | d0 | Q | ChiSq | RedChiSq | Hits | Hits Pattern | sigma(RInv) | sigma(Phi0) | sigma(d0) | sigma(T) | sigma(Z0)", "Text");
   
   // Get pixel track variables
   int EvtNumber          = t->EvtNumber;
   double pixTk_pt        = tk.getMomentum().Perp();
   // double pixTk_RInv      = tk.getRInv();
   double pixTk_eta       = tk.getMomentum().Eta();
   double pixTk_phi       = tk.getMomentum().Phi();
   double pixTk_pocaz     = tk.getPOCA().Z();
   // double pixTk_z0        = tk.getZ0();
   double pixTk_d0        = tk.getD0();
   double pixTk_charge    = tk.getCharge();
   string pixTk_Q         = "-";
   if(pixTk_charge > 0) pixTk_Q = "+";
   double pixTk_chiSq        = tk.getChi2();
   double pixTk_redChiSq     = tk.getChi2Red();
   int pixTk_Nhit            = tk.getNhit();
   int pixTk_HitsPattern     = tk.getPixelHitsPattern();
   int pixTk_Ncandhit        = tk.getNcandidatehit();
   int pixTk_CandHitsPattern = tk.getCandidatePixelHitsPattern();
   double pixTk_sigmaRInv    = tk.getSigmaRInv();
   double pixTk_sigmaPhi0    = tk.getSigmaPhi0();
   double pixTk_sigmaD0      = tk.getSigmaD0();
   double pixTk_sigmaT       = tk.getSigmaT();
   double pixTk_sigmaZ0      = tk.getSigmaZ0();
   // std::vector<TVector3> pixTk_hits     = tk.getPixelHits();
   // std::vector<TVector3> pixTk_candHits = tk.getCandidatePixelHits();
   // int pixTk_HitType                    = tk.getPixelHitType(TVector3 pixHit); // L1, L2, L3, L4, D1, D2, D3
   // TVector3 pixTk_p3                    = tk.getMomentum();
   
   // Add table rows
   tkProperties.AddRowColumn(0, auxTools.ToString(EvtNumber)         );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_pt   , 4)    );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_eta  , 4)    );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_phi  , 4)    );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_pocaz, 4)    );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_d0   , 4)    );
   tkProperties.AddRowColumn(0, pixTk_Q                              );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_chiSq)       );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_redChiSq)    );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_Nhit) + " (" + auxTools.ToString(pixTk_HitsPattern) + ")" );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_Ncandhit) + " (" + auxTools.ToString(pixTk_CandHitsPattern) + ")" );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_sigmaRInv,  3) );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_sigmaPhi0,  3) );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_sigmaD0  ,  3) );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_sigmaT   ,  3) );
   tkProperties.AddRowColumn(0, auxTools.ToString(pixTk_sigmaZ0  ,  3) );

   tkProperties.Print();
   if (bPrintRefTTPixelTrack) PrintPixTrackProperties(pixTk_index, true);
   
   return;
}


//****************************************************************************
unsigned int L1Tracks::GetNumOfStubs(const int iTrack)
//****************************************************************************
{

  unsigned int tk_nStubs = t->L1Tks_Stubs_iLayer->at(iTrack).size();

  return tk_nStubs;
}


//****************************************************************************
int L1Tracks::GetPixelIndexOfTrack(const int tk_Index)
//****************************************************************************
{
  
  int pixTk_Index = -1;
  
  // For-loop: L1PixTks
  for (int i = 0; i < (int) t->L1PixTks_Pt->size(); i++) {
    
    int iTTTrack = t->L1PixTks_TTTrackIndex->at(i);
    if (iTTTrack == tk_Index) {
      pixTk_Index = i;
      break;
    }
  }
  
  return pixTk_Index;
}


//****************************************************************************
unsigned int L1Tracks::GetNumOfPSStubs(const int iTrack)
//****************************************************************************
{

  unsigned int tk_nPSStubs = 0;
  const int tk_nStubs = (int) t->L1Tks_Stubs_iLayer->at(iTrack).size();
  for(int iStub=0; iStub < tk_nStubs; iStub++){ 
    if (t->L1Tks_Stubs_isPS->at(iTrack).at(iStub)) tk_nPSStubs++;
    else continue;
  }

  return tk_nPSStubs;
}

//****************************************************************************
unsigned int L1Tracks::GetNumOfBarrelStubs(const int iTrack)
//****************************************************************************
{

  unsigned int tk_nBarrelStubs = 0;
  vector< unsigned int > tkStubs_iLayer = t->L1Tks_Stubs_iLayer->at(iTrack);
  for(int i=0; i < (int) tkStubs_iLayer.size(); i++){ if(tkStubs_iLayer.at(i) != 999999) tk_nBarrelStubs++;};

  return tk_nBarrelStubs;
}


//****************************************************************************
unsigned int L1Tracks::GetNumOfEndcapStubs(const int iTrack)
//****************************************************************************
{

  unsigned int tk_nEndcapStubs = 0;
  vector< unsigned int > tkStubs_iDisk  = t->L1Tks_Stubs_iDisk->at(iTrack);
  for(int i=0; i < (int) tkStubs_iDisk.size(); i++){ if(tkStubs_iDisk.at(i) != 999999) tk_nEndcapStubs++;};

  return tk_nEndcapStubs;
}

//****************************************************************************
TLorentzVector L1Tracks::GetP4FromTracks(const vector<int> tks_Index)
//****************************************************************************
{
  
  TLorentzVector p4;
  if (tks_Index.size() < 1) return p4;

  // For-loop: Tks
  for (int i = 0; i < (int) tks_Index.size(); i++){

    const int index = tks_Index.at(i);
    double pt       = t->L1Tks_Pt  ->at(index);
    double eta      = t->L1Tks_Eta ->at(index);
    double phi      = t->L1Tks_Phi ->at(index);
    double mass     = pionMass; //assume charged pion mass

    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(pt, eta, phi, mass);
    p4 += tmp;
    
  } // For-loop: Tks

  return p4;
}


//****************************************************************************
TLorentzVector L1Tracks::GetP4FromPixTracks(const vector<int> tks_Index)
//****************************************************************************
{
  
  TLorentzVector p4;
  if (tks_Index.size() < 1) return p4;

  // For-loop: Tks
  for (int i = 0; i < (int) tks_Index.size(); i++){

    const int index = tks_Index.at(i);
    double pt       = t->L1PixTks_Pt  ->at(index);
    double eta      = t->L1PixTks_Eta ->at(index);
    double phi      = t->L1PixTks_Phi ->at(index);
    double mass     = pionMass; //assume charged pion mass

    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(pt, eta, phi, mass);
    p4 += tmp;
    
  } // For-loop: Tks

  return p4;
}


//****************************************************************************
TLorentzVector L1Tracks::GetP4FromPixTracksRefit(const vector<int> tks_Index)
//****************************************************************************
{
  
  L1PixelTrackFit f(3.8112); // Bz in Tesla
  TLorentzVector p4;
  if (tks_Index.size() < 1) return p4;

  // For-loop: Refitted tracks (index is for TTTrack)
  for (int i = 0; i < (int) tks_Index.size(); i++){
    
    const int tk_index = tks_Index.at(i);

    // Get TTrack properties
    double tk_RInv      = t->L1Tks_RInv->at(tk_index);
    double tk_Phi       = t->L1Tks_Phi->at(tk_index);
    double tk_d0        = GetTrackD0 (tk_index);
    double tk_Eta       = t->L1Tks_Eta->at(tk_index);
    double tk_t         = sinh(tk_Eta);
    double tk_z0        = t->L1Tks_POCAz->at(tk_index);

    // Get the corresponding TTPixelTrack
    int pixTk_index = GetPixelIndexOfTrack(tk_index);
    if (pixTk_index < 0) return p4;
    vector<double> candPixHits_X    = t->L1PixTks_CandPixHits_X->at(pixTk_index);
    vector<double> candPixHits_Y    = t->L1PixTks_CandPixHits_Y->at(pixTk_index);
    vector<double> candPixHits_Z    = t->L1PixTks_CandPixHits_Z->at(pixTk_index);
    vector< int >  candPixHits_Type = t->L1PixTks_CandPixHits_Type->at(pixTk_index);
   
    // Re-fitted TTPixelTrack
    TTPixelTrack pixFitTk  = f.FitPixelTrack(tk_RInv, tk_Phi, tk_d0, tk_t, tk_z0, candPixHits_X, candPixHits_Y, candPixHits_Z, candPixHits_Type);
    double pt   = pixFitTk.getMomentum().Perp();
    double eta  = pixFitTk.getMomentum().Eta();
    double phi  = pixFitTk.getMomentum().Phi();
    double mass = pionMass; //assume charged pion mass

    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(pt, eta, phi, mass);
    p4 += tmp;
    
  } // For-loop: Tks

  return p4;
}



//****************************************************************************
int L1Tracks::GetChargeFromTracks(const vector<int> tks_Index)
//****************************************************************************
{
  
  int charge = 0;
  if (tks_Index.size() < 1) return -9999;

  // For-loop: Tks
  for (int i = 0; i < (int) tks_Index.size(); i++){

    const int index = tks_Index.at(i);
    charge += t->L1Tks_Charge->at(index);
    
  } // For-loop: Tks

  return charge;
}


//****************************************************************************
int L1Tracks::GetChargeFromPixTracks(const vector<int> tks_Index)
//****************************************************************************
{
  
  int charge = 0;
  if (tks_Index.size() < 1) return -9999;

  // For-loop: PixTks
  for (int i = 0; i < (int) tks_Index.size(); i++){

    const int index = tks_Index.at(i);
    charge += t->L1PixTks_Charge->at(index);
    
  } // For-loop: PixTks

  return charge;
}


//****************************************************************************
int L1Tracks::GetDOF(int iTrack,
		     int nFitParams)
//****************************************************************************
{

  unsigned int tk_nStubs     = GetNumOfStubs(iTrack);
  unsigned int tk_nFitParams = nFitParams;
  int tk_dof = 2*tk_nStubs - tk_nFitParams;

  return tk_dof;
}


//****************************************************************************
double L1Tracks::GetRedChiSquared(int iTrack,
				  int nFitParams)
//****************************************************************************
{

  double tk_chiSq    = t->L1Tks_ChiSquared->at(iTrack);
  int tk_dof         = GetDOF(iTrack, nFitParams);
  double tk_redChiSq = tk_chiSq/double(tk_dof);

  return tk_redChiSq;  
}

#endif

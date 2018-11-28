#ifndef L1TkEGParticle_cxx
#define L1TkEGParticle_cxx

// User
#include "../../Auxiliary/interface/constants.h"
#include "../interface/L1TkEGParticle.h"

//****************************************************************************
L1TkEGParticle::L1TkEGParticle()
//****************************************************************************
{

  InitVars_();
  
}


//****************************************************************************
L1TkEGParticle::L1TkEGParticle(vector<TTTrack> tracks, vector<EG> EGs,
                               GenParticle genTau, bool matching)
//****************************************************************************
{

  theTracks = tracks;
  theEGs = EGs;
  theGenTau = genTau;
  theMatching = matching;
}

//****************************************************************************
L1TkEGParticle::L1TkEGParticle(double vtxIso,
			       double relIso, 
			       double CHF, 
			       double NHF,
			       double shrinkConeConst,
			       double sigConeMaxOpen,
			       vector<TTTrack> isoTracks, 
			       vector<EG> signalEGs,
			       vector<EG> isoEGs)
  
//****************************************************************************
{
  SetVtxIso(vtxIso);
  SetRelIso(relIso);
  SetCHF(CHF);
  SetNHF(NHF);
  SetShrinkingConeConst(shrinkConeConst);
  SetSigConeMaxOpen(sigConeMaxOpen);
  SetIsoConeTracks(isoTracks);
  SetIsoConeEGs(isoEGs);
  SetSignalConeEGs(signalEGs);
}



//****************************************************************************
void L1TkEGParticle::InitVars_(void)
//****************************************************************************
{
  // theMatchingGenParticle_dR = 999.9;
  // return;
}

//============================================================================ 
float L1TkEGParticle::CorrectedEta(float eta, float zTrack)  {
//============================================================================ 
  // Correct the eta of the L1EG object once we know the zTrack 
  // (normaly we use the zvertex but since we care about the dR of the EG and the track it is not needed)
  
  bool IsBarrel = ( fabs(eta) < 1.479 );
  float REcal = 129. ;
  float ZEcal = 315.4 ;

  float theta = 2. * TMath::ATan( TMath::Exp( - eta ) );
  if (theta < 0) theta = theta + TMath::Pi();
  float tantheta = TMath::Tan( theta );

  float delta;
  if (IsBarrel) {
    delta = REcal / tantheta ;
  }
  else {
    if (theta > 0) delta =  ZEcal;
    if (theta < 0) delta = -ZEcal;
  }

  float tanthetaprime = delta * tantheta / (delta - zTrack );

  float thetaprime = TMath::ATan( tanthetaprime );
  if (thetaprime < 0) thetaprime = thetaprime + TMath::Pi();

  float etaprime = -TMath::Log( TMath::Tan( thetaprime / 2.) );

  return etaprime;

}


//****************************************************************************
TLorentzVector L1TkEGParticle::GetTotalP4()
//****************************************************************************
{
  TLorentzVector sum; // initialized to (0,0,0,0)

  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    sum += tk->p4(); //NB! Assumes pion mass!
  
  for (auto eg = theEGs.begin(); eg != theEGs.end(); eg++)
    sum += eg->p4(); //NB! Approximation Pt ~ Et used in TLorentzVector!

  return sum;
}


//****************************************************************************
double L1TkEGParticle::GetSigConeMax()
//****************************************************************************
{
  double maxDeltaR = GetShrinkingConeConst()/(theTracks[0].getPt());
  if (maxDeltaR > GetSigConeMaxOpen()) maxDeltaR = GetSigConeMaxOpen();

  return maxDeltaR;
  
}

//****************************************************************************
double L1TkEGParticle::GetIsoConeMin()
//****************************************************************************
{
  return GetSigConeMax();
}

//****************************************************************************
void L1TkEGParticle::FindIsoConeTracks(vector<TTTrack> TTTracks, bool useIsoCone)
//****************************************************************************
{
  vector<TTTrack> isoTracks;
  TTTrack leadingTrack = theTracks[0];
  vector<TTTrack> clustTracks  = theTracks;
  double deltaR;

  // For-loop: All TTTracks
  for (auto tk = TTTracks.begin(); tk != TTTracks.end(); tk++) {
  
    if (useIsoCone) {
      // Check if the track is clustered in the tau candidate
      bool clustered = false;
      
      for (auto clusttk = clustTracks.begin(); clusttk != clustTracks.end(); clusttk++){
	if (clusttk->index() == tk->index()) clustered = true;
      }
      if (clustered) continue;	
    }
    
    // Check if the track is in the iso-cone
    deltaR = auxTools.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), tk->getEta(), tk->getPhi());
    //std::cout << "minDr_FUNC = "<<GetIsoConeMin()<<std::endl; 

    if (deltaR > GetIsoConeMin() && deltaR < GetIsoConeMax()) {
      isoTracks.push_back(*tk);
    }
    
  } // For-loop: All TTTracks
  
  SetIsoConeTracks(isoTracks);

  return;

}

//****************************************************************************
void L1TkEGParticle::FindSignalConeEGs(vector<EG> EGs)
//****************************************************************************
{
  vector<EG> signalEGs;
  TTTrack leadingTrack = theTracks[0];
  double deltaR;

  // For-loop: All TTTracks
  for (auto eg = EGs.begin(); eg != EGs.end(); eg++) {
  
    // Check if the eg is in the signal-cone
    deltaR = auxTools.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), CorrectedEta(eg->getEta(), leadingTrack.getZ0()), eg->getPhi());

    if (deltaR > GetSigConeMin() && deltaR < GetSigConeMax()) {
      signalEGs.push_back(*eg);
    }
    
  } // For-loop: All EGs

  SetSignalConeEGs(signalEGs);
  
  return;
}


//****************************************************************************
void L1TkEGParticle::FindIsoConeEGs(vector<EG> EGs, bool useIsoCone)
//****************************************************************************
{
  vector<EG> isoEGs;
  TTTrack leadingTrack = theTracks[0];
  vector<EG> clustEGs  = theEGs;
  double deltaR;

  // For-loop: All EGs
  for (auto eg = EGs.begin(); eg != EGs.end(); eg++) {
  
    if (useIsoCone) {
      // Check if the eg is clustered in the tau candidate
      bool clustered = false;
      
      for (auto clusteg = clustEGs.begin(); clusteg != clustEGs.end(); clusteg++){
	if (clusteg->index() == eg->index()) clustered = true;
      }
      if (clustered) continue;	
    }
    
    // Check if the eg is in the iso-cone
    deltaR = auxTools.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), CorrectedEta(eg->getEta(), leadingTrack.getZ0()), eg->getPhi());

    if (deltaR > GetIsoConeMin() && deltaR < GetIsoConeMax()) {
      isoEGs.push_back(*eg);
    }
    
  } // For-loop: All EGs
  
  SetIsoConeEGs(isoEGs);

  return;

}

//****************************************************************************
double L1TkEGParticle::CalculateVtxIso(vector<TTTrack> TTTracks, bool useIsoCone)
//****************************************************************************
{
  TTTrack leadingTrack = theTracks[0];
  vector<TTTrack> clustTracks  = theTracks;
  double mindZ0 = 999.9;
  double dz0, deltaR;
  
  // For-loop: All TTTracks
  for (auto tk = TTTracks.begin(); tk != TTTracks.end(); tk++) {
    
    if (useIsoCone) {
      // Check if the track is clustered in the tau candidate
      bool clustered = false;
      
      for (auto clusttk = clustTracks.begin(); clusttk != clustTracks.end(); clusttk++){
	if (clusttk->index() == tk->index()) clustered = true;
      }
      if (clustered) continue;	
    }
    
    // Check if the track is in the iso-cone
    deltaR = auxTools.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), tk->getEta(), tk->getPhi());
    //std::cout << "minDr_FUNC = "<<GetIsoConeMin()<<std::endl; 

    if (deltaR > GetIsoConeMin() && deltaR < GetIsoConeMax()) {
      
      // Calculate mindz0
      dz0 = abs (tk->getZ0() - leadingTrack.getZ0() );
      if (dz0 < mindZ0) mindZ0 = dz0;
      }
    
  } // For-loop: All TTTracks
  
  double vtxIso = mindZ0;
  SetVtxIso(vtxIso);
  
  return vtxIso;
}


//****************************************************************************
double L1TkEGParticle::CalculateRelIso(vector<TTTrack> TTTracks, vector<EG> EGs, 
				       double deltaZ0_max, bool useIsoCone)
//****************************************************************************
{
  TTTrack leadingTrack = theTracks[0];
  vector<TTTrack> clustTracks  = theTracks;
  vector<EG> clustEGs  = theEGs;
  double relIso = -1.0;
  double ptSum  = 0.0;
  double etSum  = 0.0;
  double deltaR;

  // For-loop: All TTTracks
  for (auto tk = TTTracks.begin(); tk != TTTracks.end(); tk++) {
    
    if (useIsoCone) {
      // Check if the track is clustered in the tau candidate
      bool clustered = false;
      
      for (auto clusttk = clustTracks.begin(); clusttk != clustTracks.end(); clusttk++){
	if (clusttk->index() == tk->index()) clustered = true;
      }
      if (clustered) continue;	
    }
    
    // Apply dz0 cut
    if (abs (tk->getZ0() - leadingTrack.getZ0() ) > deltaZ0_max) continue;
    
    // Check if the track is in the iso-cone
    deltaR = auxTools.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), tk->getEta(), tk->getPhi());
  
    if (deltaR > GetIsoConeMin() && deltaR < GetIsoConeMax()) ptSum += tk -> getPt();
    
  } // For-loop: All TTTracks
  

  // For-loop: All EGs
  for (auto eg = EGs.begin(); eg != EGs.end(); eg++) {
  
    if (useIsoCone) {
      // Check if the eg is clustered in the tau candidate
      bool clustered = false;
      
      for (auto clusteg = clustEGs.begin(); clusteg != clustEGs.end(); clusteg++){
	if (clusteg->index() == eg->index()) clustered = true;
      }
      if (clustered) continue;	
    }
    
    // Check if the eg is in the iso-cone
    deltaR = auxTools.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), CorrectedEta(eg->getEta(), leadingTrack.getZ0()), eg->getPhi());

    if (deltaR > GetIsoConeMin() && deltaR < GetIsoConeMax()) etSum += eg -> getEt();
    
  } // For-loop: All EGs
  
  // Calculate relative isolation
  relIso = (ptSum + etSum) / GetTotalEt();
  SetRelIso(relIso);
  
  return relIso;
}




//****************************************************************************
double L1TkEGParticle::GetTrackBasedPt()
//****************************************************************************
{
  TLorentzVector sum; // initialized to (0,0,0,0)
  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    sum += tk->p4();
  return sum.Pt();
}


//****************************************************************************
double L1TkEGParticle::GetTotalPt()
//**************************************************************************** 
{
  TLorentzVector sum; // initialized to (0,0,0,0)

  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    sum += tk->p4(); //NB! Assumes pion mass!
  
  for (auto eg = theEGs.begin(); eg != theEGs.end(); eg++)
    sum += eg->p4(); //NB! Approximation Pt ~ Et used in TLorentzVector!

  return sum.Pt();
  
}


//****************************************************************************
double L1TkEGParticle::GetTrackInvMass()
//****************************************************************************
{
  TLorentzVector p4sum; // initialized to (0,0,0,0)
  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    p4sum += tk->p4();
  return p4sum.M();
}


//****************************************************************************
double L1TkEGParticle::GetEGInvMass()
//****************************************************************************
{
  TLorentzVector p4sum; // initialized to (0,0,0,0)
  for (auto eg = theEGs.begin(); eg != theEGs.end(); eg++)
    p4sum += eg->p4();
  return p4sum.M();
}


//****************************************************************************
double L1TkEGParticle::GetGenTauPt()
//****************************************************************************
{
  if (theMatching)
    return theGenTau.p4vis().Pt();
  else
    return -1.0;
}
	
	
//****************************************************************************
double L1TkEGParticle::GetGenTauEt()
//****************************************************************************
{
  if (theMatching)
    return theGenTau.p4vis().Et();
  else
    return -1.0;
}				 


//****************************************************************************
double L1TkEGParticle::GetTrackBasedEt()
//****************************************************************************
{
  TLorentzVector sum; // initialized to (0,0,0,0)
  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    sum += tk->p4(); //NB! Assumes pion mass!
  return sum.Et();
}


//****************************************************************************
double L1TkEGParticle::GetEGBasedEt()
//****************************************************************************
{
  TLorentzVector sum; // initialized to (0,0,0,0)
  for (auto eg = theEGs.begin(); eg != theEGs.end(); eg++)
    sum += eg->p4(); //NB! Approximation Pt ~ Et used in TLorentzVector!
  return sum.Et();
  
}

//****************************************************************************
double L1TkEGParticle::GetTotalEt()
//**************************************************************************** 
{
  TLorentzVector sum; // initialized to (0,0,0,0)

  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    sum += tk->p4(); //NB! Assumes pion mass!
  
  for (auto eg = theEGs.begin(); eg != theEGs.end(); eg++)
    sum += eg->p4(); //NB! Approximation Pt ~ Et used in TLorentzVector!

  return sum.Et();
  
}

//****************************************************************************
void L1TkEGParticle::PrintTTTracks()
//****************************************************************************
{

Table table(" # | Pt | Eta | Phi | z0 (cm) | d0 (cm) | Q | Chi2 | DOF | Chi2Red | Stubs (PS)", "Text");

// For-loop: All Tracks
for (size_t i = 0; i < theTracks.size(); i++)
  {

TTTrack tk = theTracks.at(i);
      
// Fill table
table.AddRowColumn(i, auxTools.ToString(i+1) );
table.AddRowColumn(i, auxTools.ToString( tk.getPt() , 3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getEta(), 3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getPhi(), 3  ) );
// table.AddRowColumn(i, auxTools.ToString( tk.getX0() , 3  ) );
// table.AddRowColumn(i, auxTools.ToString( tk.getY0() , 3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getZ0() , 3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getD0() , 3  ) );
//table.AddRowColumn(i, auxTools.ToString( tk.getQ()  , 3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getChi2(),3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getDOF()     ) );
table.AddRowColumn(i, auxTools.ToString( tk.getChi2Red(), 3 ) );
table.AddRowColumn(i, auxTools.ToString( tk.getNumOfStubs()) ); //+ " (" + auxTools.ToString(tk.getNumOfStubsPS()) + ")");
}
  
if (theTracks.size() > 0) table.Print();
  
return;
}


#endif

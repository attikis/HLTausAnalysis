#ifndef L1TkEMParticle_cxx
#define L1TkEMParticle_cxx

// User
#include "../../Auxiliary/interface/constants.h"
#include "../interface/L1TkEMParticle.h"

//****************************************************************************
L1TkEMParticle::L1TkEMParticle()
//****************************************************************************
{

  InitVars_();
  
}


//****************************************************************************
L1TkEMParticle::L1TkEMParticle(vector<TTTrack> tracks, vector<L1CaloTP> ecalTPs,
                               GenParticle genTau, bool matching)
//****************************************************************************
{

  theTracks = tracks;
  theEcalTPs = ecalTPs;
  theGenTau = genTau;
  theMatching = matching;
}


//****************************************************************************
void L1TkEMParticle::InitVars_(void)
//****************************************************************************
{
  // theMatchingGenParticle_dR = 999.9;
  // return;
}


//****************************************************************************
double L1TkEMParticle::GetTrackBasedPt()
//****************************************************************************
{
  TLorentzVector sum; // initialized to (0,0,0,0)
  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    sum += tk->p4();
  return sum.Pt();
}


//****************************************************************************
double L1TkEMParticle::GetTrackInvMass()
//****************************************************************************
{
  TLorentzVector p4sum; // initialized to (0,0,0,0)
  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    p4sum += tk->p4();
  return p4sum.M();
}


//****************************************************************************
double L1TkEMParticle::GetEMInvMass()
//****************************************************************************
{
  TLorentzVector p4sum; // initialized to (0,0,0,0)
  for (auto tp = theEcalTPs.begin(); tp != theEcalTPs.end(); tp++)
    p4sum += tp->p4();
  return p4sum.M();
}


//****************************************************************************
double L1TkEMParticle::GetGenTauPt()
//****************************************************************************
{
  if (theMatching)
    return theGenTau.p4vis().Pt();
  else
    return -1.0;
}
	
	
//****************************************************************************
double L1TkEMParticle::GetGenTauEt()
//****************************************************************************
{
  if (theMatching)
    return theGenTau.p4vis().Et();
  else
    return -1.0;
}				 


//****************************************************************************
double L1TkEMParticle::GetTrackBasedEt()
//****************************************************************************
{
  TLorentzVector sum; // initialized to (0,0,0,0)
  for (auto tk = theTracks.begin(); tk != theTracks.end(); tk++)
    sum += tk->p4(); //NB! Assumes pion mass!
  return sum.Et();
}


//****************************************************************************
double L1TkEMParticle::GetEMBasedEt()
//****************************************************************************
{
  TLorentzVector sum; // initialized to (0,0,0,0)
  for (auto tp = theEcalTPs.begin(); tp != theEcalTPs.end(); tp++)
    sum += tp->p4(); //NB! Approximation Pt ~ Et used in TLorentzVector!
  return sum.Et();
}


//****************************************************************************
void L1TkEMParticle::PrintTTTracks()
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
table.AddRowColumn(i, auxTools.ToString( tk.getQ()  , 3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getChi2(),3  ) );
table.AddRowColumn(i, auxTools.ToString( tk.getDOF()     ) );
table.AddRowColumn(i, auxTools.ToString( tk.getChi2Red(), 3 ) );
table.AddRowColumn(i, auxTools.ToString( tk.getNumOfStubs()) + " (" + auxTools.ToString(tk.getNumOfStubsPS()) + ")");
}
  
if (theTracks.size() > 0) table.Print();
  
return;
}


#endif

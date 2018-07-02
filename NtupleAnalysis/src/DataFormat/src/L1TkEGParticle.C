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
L1TkEGParticle::L1TkEGParticle(vector<TTTrack> tracks, vector<L1EG> EGs,
                               GenParticle genTau, bool matching)
//****************************************************************************
{

  theTracks = tracks;
  theEGs = EGs;
  theGenTau = genTau;
  theMatching = matching;
}


//****************************************************************************
void L1TkEGParticle::InitVars_(void)
//****************************************************************************
{
  // theMatchingGenParticle_dR = 999.9;
  // return;
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

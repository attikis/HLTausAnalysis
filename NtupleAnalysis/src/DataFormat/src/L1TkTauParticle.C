#ifndef L1TkTauParticle_cxx
#define L1TkTauParticle_cxx

// User
#include "../../Auxiliary/interface/constants.h"
#include "../interface/L1TkTauParticle.h"

//****************************************************************************
L1TkTauParticle::L1TkTauParticle()
//****************************************************************************
{
  
}

//****************************************************************************
L1TkTauParticle::L1TkTauParticle(double matchCone_dRMin,
				 double matchCone_dRMax,
				 double sigCone_dRMin,
				 double sigCone_dRMax,
				 double isoCone_dRMin,
				 double isoCone_dRMax)
//****************************************************************************
{

  theMatchCone_dRMin  = matchCone_dRMin;
  theMatchCone_dRMax  = matchCone_dRMax;
  theSigCone_dRMin    = sigCone_dRMin;
  theSigCone_dRMax    = sigCone_dRMax;
  theIsoCone_dRMin    = isoCone_dRMin;
  theIsoCone_dRMax    = isoCone_dRMax;

  InitVars_();

}


//****************************************************************************
void L1TkTauParticle::InitVars_(void)
//****************************************************************************
{
  theMatchingGenParticle_dR = 999.9;
  return;
}

				 
//****************************************************************************
L1TkTauParticle::L1TkTauParticle(int caloTau_Index, 
				 int matchTk_Index, 
				 double matchTk_deltaR, 
				 vector<int> sigTks_Index, 
				 vector<int> isoTks_Index,
				 double vtxIso, 
				 double relIso,
				 double sigCone_minDeltaR, 
				 double sigCone_maxDeltaR, 
				 double isoCone_minDeltaR, 
				 double isoCone_maxDeltaR)
//****************************************************************************
{
  
  SetCaloTau(caloTau_Index);
  SetMatchTk(matchTk_Index);
  SetMatchTkDeltaR(matchTk_deltaR);
  SetSigConeTks(sigTks_Index);
  SetIsoConeTks(isoTks_Index);
  SetVtxIso(vtxIso);
  SetRelIso(relIso);
  SetMatchGenp(-1.0, 9999.9);
  SetSignalConeSize(sigCone_minDeltaR,  sigCone_maxDeltaR);
  SetIsolationConeSize(isoCone_minDeltaR, isoCone_maxDeltaR);
}



//****************************************************************************
void L1TkTauParticle::SetSignalConeSize(double deltaR_min, double deltaR_max)
//****************************************************************************
{
  sigCone_minDeltaR_ = deltaR_min;
  sigCone_maxDeltaR_ = deltaR_max;
  
  return;
}



//****************************************************************************
void L1TkTauParticle::SetIsolationConeSize(double deltaR_min, double deltaR_max)
//****************************************************************************
{
  isoCone_minDeltaR_ = deltaR_min;
  isoCone_maxDeltaR_ = deltaR_max;
  
  return;
}



//============================================================================
double L1TkTauParticle::CalculateRelIso(const double deltaZ0_max,
					bool bStoreValue,
					bool bInvert_deltaZ0)

//============================================================================
{

  // Store default values
  // SetRelIsolation(0.0);
  
  // Return not Tk-Confirmed
  if (!HasMatchingTk()) return 0.0; 

  // If no tracks found in the isoalation cone return
  vector<TTTrack> isoConeTks = GetIsoConeTTTracks();
  if ( (isoConeTks.size() < 1) )  return 0.0;

  // Initialise variables
  TTTrack matchTk = GetMatchingTk();
  double isoTks_scalarSumPt  = 0.0;
  double deltaZ0 = 999.9;
  double relIso  = 0.0;
  
  // For-loop: All Tracks in isolation cone 
  for (size_t i = 0; i < isoConeTks.size(); i++)
    {
      TTTrack isoConeTk = isoConeTks.at(i);
      
      
      // Find the track closest in Z0
      deltaZ0 = abs(matchTk.getZ0() - isoConeTk.getZ0());

      // Decide on type of calculation
      bool considerTk  = false;
      bool considerTk_default = (deltaZ0 < deltaZ0_max);
      bool considerTk_invert  = (deltaZ0 > deltaZ0_max);
      if (bInvert_deltaZ0) considerTk = considerTk_invert;
      else considerTk = considerTk_default;
      
      // Add-up the pT of alltracks in isolation cone/annulus
      if (considerTk) isoTks_scalarSumPt += isoConeTk.getPt();
    }

  // Calculated + Assign value of relative isolation
  relIso = isoTks_scalarSumPt/matchTk.getPt();
  if (bStoreValue) SetRelIsolation(relIso);
  
  return relIso;
}


//============================================================================
double L1TkTauParticle::CalculateVtxIso(bool bStoreValue)
//============================================================================
{

  // Store default values
  double deltaZ0 = 999.9;
  double deltaZ0_tmp = 999.9;
  if (bStoreValue) SetVtxIsolation(deltaZ0);
  
  // Return not Tk-Confirmed
  if (!HasMatchingTk()) return 999.9; 

  // If no tracks found in the isoalation cone return
  vector<TTTrack> isoConeTks = GetIsoConeTTTracks();
  if ( (isoConeTks.size() < 1) )  return 999.9;

  // Initialise variables
  TTTrack matchTk = GetMatchingTk();
  
  // For-loop: All Tracks in isolation cone 
  for (size_t i = 0; i < isoConeTks.size(); i++)
    {
      TTTrack isoConeTk = isoConeTks.at(i);
      
      // Find the track closest in Z0
      deltaZ0_tmp = abs(matchTk.getZ0() - isoConeTk.getZ0());
      
      if (deltaZ0_tmp < deltaZ0)
      	{
	  deltaZ0 = deltaZ0_tmp;
      	  if (bStoreValue) SetVtxIsolation(deltaZ0);
	  if (bStoreValue) SetVtxIsolationTrack(isoConeTk);
      	}
    }

  // Calculated + Assign value of relative isolation
  return deltaZ0;
}


//****************************************************************************
void L1TkTauParticle::SetMatchGenp(int matchGenp_Index, double matchGenp_deltaR) 
//****************************************************************************
{ 
  matchGenp_Index_  = matchGenp_Index;
  matchGenp_deltaR_ = matchGenp_deltaR;
  
  return;
}



// //****************************************************************************
// void L1TkTauParticle::PrintProperties(void)
// //****************************************************************************
// {
  
//   Table general("iCalo | Match Tk | Match GenP | dR (Match GenP) | VtxIso (cm) | RelIso | Sig. Tks | Iso. Tks. | Sig_R (min) | Sig_R (max) | Iso_R (min) | Iso_R (max)", "Text");
//   general.AddRowColumn(0, auxTools.ToString( caloTau_Index_    ) );
//   general.AddRowColumn(0, auxTools.ToString( matchTk_Index_    ) );
//   general.AddRowColumn(0, auxTools.ToString( matchGenp_Index_  ) );
//   general.AddRowColumn(0, auxTools.ToString( matchGenp_deltaR_ ) );
//   general.AddRowColumn(0, auxTools.ToString( vtxIso_           ) );
//   general.AddRowColumn(0, auxTools.ToString( relIso_           ) );
//   general.AddRowColumn(0, auxTools.ConvertIntVectorToString(sigTks_Index_) );
//   general.AddRowColumn(0, auxTools.ConvertIntVectorToString(isoTks_Index_) );
//   general.AddRowColumn(0, auxTools.ToString(sigCone_minDeltaR_ ) );
//   general.AddRowColumn(0, auxTools.ToString(sigCone_maxDeltaR_ ) );
//   general.AddRowColumn(0, auxTools.ToString(isoCone_minDeltaR_ ) );
//   general.AddRowColumn(0, auxTools.ToString(isoCone_maxDeltaR_ ) );

//   general.Print();

//   return;
//}


//****************************************************************************
TTTrack L1TkTauParticle::GetSigConeLdgTk(void)
//****************************************************************************
{
  
  // Temporarily set the leading track
  vector<TTTrack> theTracks = GetSigConeTTTracks();
  TTTrack ldgTk = theTracks.at(0);
  
  // For-loop: All Tracks
  for (size_t i = 0; i < theTracks.size(); i++)
    {
      
      TTTrack tk = theTracks.at(i);
      // ldgTk = tk; //marina : this shouldn't be here

      // Find new leading track
      if (tk.getPt() > ldgTk.getPt() ) ldgTk = tk;
      
    }

  return ldgTk;  
}



//****************************************************************************
TTTrack L1TkTauParticle::GetIsoConeLdgTk(void)
//****************************************************************************
{

  // Temporarily set the leading track
  vector<TTTrack> theTracks = GetIsoConeTTTracks();
  TTTrack ldgTk = theTracks.at(0);
  
  // For-loop: All Tracks
  for (size_t i = 0; i < theTracks.size(); i++)
    {

      TTTrack tk = theTracks.at(i);
      ldgTk = tk;

      // Find new leading track
      if (tk.getPt() > ldgTk.getPt() ) ldgTk = tk;
      
    }
  
  return ldgTk;  
}


//****************************************************************************
void L1TkTauParticle::PrintProperties(bool bPrintCaloTau,
				      bool bPrintMatchTk,
				      bool bPrintSigConeTks,
				      bool bPrintIsoConeTks,
				      bool bPrintMatchGenParticle)
//****************************************************************************
{
  
  Table info("Match-Cone | Sig-Cone | Iso-Cone | Calo-Et | Calo-Eta | Tk-Pt | Tk-Eta | Tk-dR | Gen-Pt | Gen-Eta | Gen-dR | Sig-Tks | Iso-Tks | VtxIso | RelIso", "Text");
  info.AddRowColumn(0, auxTools.ToString( GetMatchConeMin(), 2) + " < dR < " + auxTools.ToString( GetMatchConeMax(), 2) );
  info.AddRowColumn(0, auxTools.ToString( GetSigConeMin()  , 2) + " < dR < " + auxTools.ToString( GetSigConeMax()  , 2) );
  info.AddRowColumn(0, auxTools.ToString( GetIsoConeMin()  , 2) + " < dR < " + auxTools.ToString( GetIsoConeMax()  , 2) ); 
  info.AddRowColumn(0, auxTools.ToString( GetCaloTau().et()       ,3  ) );
  info.AddRowColumn(0, auxTools.ToString( GetCaloTau().eta()      ,3  ) );
  info.AddRowColumn(0, auxTools.ToString( GetMatchingTk().getPt() ,3  ) );
  info.AddRowColumn(0, auxTools.ToString( GetMatchingTk().getEta(),3  ) );
  info.AddRowColumn(0, auxTools.ToString( GetMatchingTkDeltaR()   ,3  ) );
  info.AddRowColumn(0, auxTools.ToString( GetMatchingGenParticle().pt() , 3) );
  info.AddRowColumn(0, auxTools.ToString( GetMatchingGenParticle().eta(), 3) );
  info.AddRowColumn(0, auxTools.ToString( GetMatchingGenParticleDeltaR(), 3) );  
  info.AddRowColumn(0, auxTools.ToString( GetSigConeTTTracks().size() ) );
  info.AddRowColumn(0, auxTools.ToString( GetIsoConeTTTracks().size() ) );
  info.AddRowColumn(0, auxTools.ToString( GetVtxIsolation(), 3) );
  info.AddRowColumn(0, auxTools.ToString( GetRelIsolation(), 3) );
  info.Print();
  
  if (bPrintCaloTau) GetCaloTau().PrintProperties();
  if (bPrintMatchTk && HasMatchingTk()) GetMatchingTk().PrintProperties();

  if (bPrintSigConeTks)
    {
      PrintTTTracks(GetSigConeTTTracks(), "Sig-Cone Tks");
      PrintTTPixelTracks(GetSigConeTTPixelTracks(), "Sig-Cone Tks");
    }
  
  if (bPrintIsoConeTks)
    {
      PrintTTTracks(GetIsoConeTTTracks(), "Iso-Cone Tks");
      PrintTTPixelTracks(GetIsoConeTTPixelTracks(), "Iso-Cone Tks");
    }

  if (bPrintMatchGenParticle) GetMatchingGenParticle().PrintProperties();

      
  return;
}


//****************************************************************************
void L1TkTauParticle::PrintTTTracks(vector<TTTrack> theTracks,
				      string theTrackType)
//****************************************************************************
{

// Table table(theTrackType + " # | Pt | Eta | Phi | x0 | y0 | z0 (cm) | d0 (cm) | Q | Chi2 | DOF | Chi2Red | Stubs (PS)", "Text");
Table table(theTrackType + " # | Pt | Eta | Phi | z0 (cm) | d0 (cm) | Q | Chi2 | DOF | Chi2Red | Stubs (PS)", "Text");

// For-loop: All Tracks
 for (size_t i = 0; i < theTracks.size(); i++)
   {
     
     TTTrack tk = theTracks.at(i);
     
     // Fill table
     table.AddRowColumn(i, auxTools.ToString(i+1) );
     table.AddRowColumn(i, auxTools.ToString( tk.getPt() , 3  ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getEta(), 3  ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getPhi(), 3  ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getZ0() , 3  ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getD0() , 3  ) );
     //table.AddRowColumn(i, auxTools.ToString( tk.getQ()  , 3  ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getChi2(),3  ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getDOF()     ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getChi2Red(), 3 ) );
     table.AddRowColumn(i, auxTools.ToString( tk.getNumOfStubs()) );// + " (" + auxTools.ToString(tk.getNumOfStubsPS()) + ")");
   }
  
if (theTracks.size() > 0) table.Print();
  
return;
}


//****************************************************************************
void L1TkTauParticle::PrintTTPixelTracks(vector<TTPixelTrack> theTracks,
					 string theTrackType)
//****************************************************************************
{

Table table(theTrackType + " # | Pt | Eta | Phi | z0 | d0 | Q | Chi2 | RedChi2 | Hits | Hit-Pattern | Hit-Type | Hit-R | Hit-Z", "Text");

// For-loop: All Tracks
for (size_t i = 0; i < theTracks.size(); i++)
  {

TTPixelTrack tk = theTracks.at(i);

// Fill table
table.AddRowColumn(i, auxTools.ToString(i+1) );
table.AddRowColumn(i, auxTools.ToString( tk.getPt()  , 3) );
table.AddRowColumn(i, auxTools.ToString( tk.getEta() , 3) );
table.AddRowColumn(i, auxTools.ToString( tk.getPhi() , 3) );
table.AddRowColumn(i, auxTools.ToString( tk.getZ0()  , 3) );
table.AddRowColumn(i, auxTools.ToString( tk.getD0()  , 3) );
table.AddRowColumn(i, tk.getQ() );
table.AddRowColumn(i, auxTools.ToString( tk.getChi2()   , 3 ) );
table.AddRowColumn(i, auxTools.ToString( tk.getChi2Red(), 3) );
table.AddRowColumn(i, auxTools.ToString( tk.getNhit()) + " (" + auxTools.ToString(tk.getCandidatePixelHits().size()) + ")" );
table.AddRowColumn(i, auxTools.ToString( tk.getPixelHitsPattern()) );
table.AddRowColumn(i, auxTools.ConvertIntVectorToString( tk.getPixHitsType()) );
table.AddRowColumn(i, auxTools.ConvertIntVectorToString( tk.getPixHitsR() )   );
}
  
if (theTracks.size() > 0) table.Print();
  
return;
}


//****************************************************************************
void L1TkTauParticle::SetSigConeTTTracksP4_(void)
//****************************************************************************
{

  // Variables
  vector<TTTrack> TTTracks = GetSigConeTTTracks();
  TLorentzVector sigTks_p4(0, 0, 0, 0);
  for (vector<TTTrack>::iterator t = TTTracks.begin(); t != TTTracks.end(); t++)
    {
      TLorentzVector p4;
      p4.SetPtEtaPhiM(t->getPt(), t->getEta(), t->getPhi(), pionMass); // assume track is a pion (most likely true!)
      sigTks_p4 += p4;
    }

  theSigConeTTTracksP4 = sigTks_p4;
  return;
}


//****************************************************************************
void L1TkTauParticle::SetIsoConeTTTracksP4_(void)
//****************************************************************************
{

  // Variables
  vector<TTTrack> TTTracks = GetIsoConeTTTracks();
  TLorentzVector isoTks_p4(0, 0, 0, 0);
  for (vector<TTTrack>::iterator t = TTTracks.begin(); t != TTTracks.end(); t++)
    {
      TLorentzVector p4;
      p4.SetPtEtaPhiM(t->getPt(), t->getEta(), t->getPhi(), pionMass); // assume track is a pion (most likely true!)
      isoTks_p4 += p4;
    }

  theIsoConeTTTracksP4 = isoTks_p4;
  return;
}


//****************************************************************************
TLorentzVector L1TkTauParticle::GetSigConeTTTracksP4(void)
//****************************************************************************
{
  SetSigConeTTTracksP4_();
  return theSigConeTTTracksP4;
 
}

//****************************************************************************
TLorentzVector L1TkTauParticle::GetIsoConeTTTracksP4(void)
//****************************************************************************
{

  SetIsoConeTTTracksP4_();
  return theIsoConeTTTracksP4;
}

#endif

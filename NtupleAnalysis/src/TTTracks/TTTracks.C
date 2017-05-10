#ifndef TTTracks_cxx
#define TTTracks_cxx

// User
#include "TTTracks.h"
// C++
#include <stdexcept>

//============================================================================
void TTTracks::InitVars_()
//============================================================================
{

  DEBUG     = false;
  datasets_ = datasets_.GetDataset(mcSample);

  tk_Collection  =  "TTTracks"; // Default: Only "TTTracks" available (not "TTPixelTracks")
  tk_nFitParams  = 4;           // Default: 4
  tk_minPt       = 2.00;        // Default: 2.0
  tk_minEta      = 0.0;         // Default: 0.0
  tk_maxEta      = 1e6;         // Default: 1e6
  tk_maxChiSqRed = 1e6;         // Default: 1e6
  tk_minStubs    =   0;         // Default: 0
  tk_minStubsPS  =   0;         // Default: 0
  tk_maxStubsPS  = 1e6;         // Default: 1e6

  PrintSettings();

  return;
}


//============================================================================
void TTTracks::PrintSettings(void)
//============================================================================
{

  // Inform user of settings
  Table settings("Variable | Cut | Value | Default | Units", "Text");  // Table settingsTable("Variable & Value & Units", "LaTeX", "l l l");
  settings.AddRowColumn(0, "MC Sample");
  settings.AddRowColumn(0, "==" );
  settings.AddRowColumn(0, mcSample );
  settings.AddRowColumn(0, "");

  settings.AddRowColumn(1, "Tracks: Collection");
  settings.AddRowColumn(1, "==");
  settings.AddRowColumn(1, tk_Collection);
  settings.AddRowColumn(1, "TTTracks");
  settings.AddRowColumn(1, "");
  
  settings.AddRowColumn(2, "Tracks: Fit Parameters");
  settings.AddRowColumn(2, "==");
  settings.AddRowColumn(2, auxTools_.ToString( tk_nFitParams) );
  settings.AddRowColumn(2, "4");
  settings.AddRowColumn(2, "");

  settings.AddRowColumn(3, "Tracks: Pt");
  settings.AddRowColumn(3, ">=");
  settings.AddRowColumn(3, auxTools_.ToString( tk_minPt) );
  settings.AddRowColumn(3, "2" );
  settings.AddRowColumn(3, "GeV/c" );
  
  settings.AddRowColumn(4, "Tracks: |Eta|");
  settings.AddRowColumn(4, ">=");
  settings.AddRowColumn(4, auxTools_.ToString( tk_minEta) );
  settings.AddRowColumn(4, "0.0" );
  settings.AddRowColumn(4, "" );

  settings.AddRowColumn(5, "Tracks: |Eta|");
  settings.AddRowColumn(5, "<=");
  settings.AddRowColumn(5, auxTools_.ToString(tk_maxEta) );
  settings.AddRowColumn(5, auxTools_.ToString(2.5) );
  settings.AddRowColumn(5, "" );
  
  settings.AddRowColumn(6, "Tracks: ChiSqRed");
  settings.AddRowColumn(6, "<=");
  settings.AddRowColumn(6, auxTools_.ToString( tk_maxChiSqRed) );
  settings.AddRowColumn(6, "200/DOF");
  settings.AddRowColumn(6, "");

  settings.AddRowColumn(7, "Tracks: Stubs");
  settings.AddRowColumn(7, ">=");
  settings.AddRowColumn(7, auxTools_.ToString(tk_minStubs) );
  settings.AddRowColumn(7, auxTools_.ToString(4) );
  settings.AddRowColumn(7, "" );
  
  settings.AddRowColumn(8, "Tracks: PS-Stubs (min)");
  settings.AddRowColumn(8, ">=");
  settings.AddRowColumn(8, auxTools_.ToString(tk_minStubsPS) );
  settings.AddRowColumn(8, auxTools_.ToString(3) );
  settings.AddRowColumn(8, "" );

  settings.AddRowColumn(9, "Tracks: PS-Stubs (max)");
  settings.AddRowColumn(9, "<=");
  settings.AddRowColumn(9, auxTools_.ToString(tk_maxStubsPS) );
  settings.AddRowColumn(9, auxTools_.ToString(6) );
  settings.AddRowColumn(9, "" );

  settings.Print();
  
  return;
}

//============================================================================
void TTTracks::Loop()
//============================================================================
{

  // Sanity check
  if (fChain == 0) return;
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  cout << "=== TTTracks:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;

  // Initialisations
  InitVars_();
  BookHistos_();
  Long64_t nbytes       = 0;
  Long64_t nb           = 0;
  unsigned int nAllEvts = fChain->GetEntries();
  unsigned int nEvts    = 0;


  // For-loop: All TTree Entries
  for (int jentry = 0; jentry < nEntries; jentry++, nEvts++)
    {

      if(DEBUG) cout << "\tEntry = " << jentry << endl;

      // Load the tree && Get the entry
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      // GenParticles Collection
      vector<GenParticle> GenParticles = GetGenParticles(false);
      if (DEBUG) PrintGenParticleCollection(GenParticles);

      
      // Tracking Particle Collections
      vector<TrackingParticle> TPs  = GetTrackingParticles(false);
      sort( TPs.begin(), TPs.end(), PtComparatorTP() ); // not sorted by default
      if (DEBUG) PrintTrackingParticleCollection(TPs);

      // TTTracks Collections
      vector<TTTrack> TTTracks = GetTTTracks(tk_minPt, tk_minEta, tk_maxEta, tk_maxChiSqRed, tk_minStubs, tk_minStubsPS, tk_maxStubsPS, tk_nFitParams, false);
      sort( TTTracks.begin(), TTTracks.end(), PtComparatorTTTrack() ); // not sorted by default
      if (DEBUG) PrintTTTrackCollection(TTTracks);

    }//eof: Entries

  
  ////////////////////////////////////////////////
  // Fill histograms
  ////////////////////////////////////////////////
  hHepMCEvt_VtxX_VtxY->Fill(HepMCEvt_VtxX, HepMCEvt_VtxY);
  hHepMCEvt_VtxZ->Fill(HepMCEvt_VtxZ);

  
  ////////////////////////////////////////////////
  // Fill counters
  ////////////////////////////////////////////////
  hCounters->SetBinContent(1, nAllEvts);
  hCounters->SetBinContent(2, nEvts);
    
    
  ////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////
  outFile->cd();
  outFile->Write();
  auxTools_.StopwatchStop(5, "minutes", "Total Time");
    
}


//============================================================================
void TTTracks::BookHistos_(void)
//============================================================================
{
  
  // Event-Type Histograms
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 2, 0.0, +2.0);
  histoTools_.BookHisto_1D(hHepMCEvt_VtxZ      , "HepMCEvt_VtxZ"     ,  "", 600, -30.0,  +30.0);
  histoTools_.BookHisto_2D(hHepMCEvt_VtxX_VtxY , "HepMCEvt_VtxX_VtxY",  "", 400,  -0.01,  +0.01, 400,  -0.01,  +0.01);
  
  return;
}



//============================================================================
void TTTracks::GetHadronicTauFinalDaughters(GenParticle p,
					    vector<unsigned short> &Daug)
//============================================================================
{
  //
  // Description:
  // Keep only the detector-related particles (tracker/calo)
  //  pi+/-, pi0, K+/-, K0, K0L, KOS.
  // Avoid leptonic decays.
  //
  // Additional info:
  // We keep the K+/-(pdgId=321) and K0L(pdgId=130) and we don't get their decay
  // products (pi+/-) because typically they live long enought to reach the ECAL/HCAL
  // c*tau (K+/-) =  3.713 m
  // c*tau (K0L)  = 15.33  m
  // c*tau (K0S)  =  0.026842 m
  //
  
  // Sanity check
  if (p.daughtersIndex().size() < 1) return;

  // Define the allowed final daugther pdgIds
  std::vector<unsigned int> pdgIds;
  pdgIds.push_back(111); // pi0
  pdgIds.push_back(211); // pi+
  pdgIds.push_back(16);  // nu_tau
  pdgIds.push_back(130); // K0_L
  pdgIds.push_back(310); // K0_S
  pdgIds.push_back(311); // K0
  pdgIds.push_back(321); // K+
  // pdgIds.push_back(223);   // omega0
  // pdgIds.push_back(221);   // eta0
  // pdgIds.push_back(213);   // rho(770)+
  // pdgIds.push_back(20213); // alpha1(1260)+ 
    

  // For-loop: All daughters
  for (unsigned int i = 0; i < p.daughtersIndex().size(); i++) 
    {

      unsigned int index = p.daughtersIndex().at(i);
      GenParticle d      = GetGenParticle(index);
      int genD_pdgId     = abs(d.pdgId());


      // Check if particle is of interest
      if ( std::find(pdgIds.begin(), pdgIds.end(), genD_pdgId) == pdgIds.end() )
	{
	  GetHadronicTauFinalDaughters(d, Daug);
      	}
      else
	{
	  Daug.push_back( d.index() );
	}
    }
  
  return;

}


//============================================================================
vector<GenParticle> TTTracks::GetGenParticles(bool bPrintList)
//============================================================================
{
  vector<GenParticle> theGenParticles;
  Table info("Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | Mothers | Daughters |", "Text"); //LaTeX or Text    
  	
  for (Size_t genP_index = 0; genP_index < GenP_Pt->size(); genP_index++)
    {
      GenParticle p = GetGenParticle(genP_index);
      SetGenParticleMomsAndDaus(p);
      theGenParticles.push_back(p);

      info.AddRowColumn(genP_index, auxTools_.ToString(p.index())  );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.pdgId())  );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.status()) );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.charge()) );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.pt())     );
      if (p.pdgId() != 2212) info.AddRowColumn(genP_index, auxTools_.ToString(p.eta())    );
      else info.AddRowColumn(genP_index, "infty");
      info.AddRowColumn(genP_index, auxTools_.ToString(p.phi())    );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.energy()) );
      info.AddRowColumn(genP_index, "(" + auxTools_.ToString(p.vx())  + ", " + auxTools_.ToString(p.vy())  + ", " + auxTools_.ToString(p.vz())  + ", " + ")");
      if (p.mothersIndex().size() < 4)
	{
	  info.AddRowColumn(genP_index, auxTools_.ConvertIntVectorToString(p.mothersIndex() ) );
	}
      else info.AddRowColumn(genP_index, ".. Too many .." );
      if (p.daughtersIndex().size() < 6)
	{
	  info.AddRowColumn(genP_index, auxTools_.ConvertIntVectorToString(p.daughtersIndex()) );
	}
      else info.AddRowColumn(genP_index, ".. Too many .." );

      // p.PrintProperties();
      // p.PrintDaughters();
      // p.PrintFinalDaughters();
      // p.PrintMothers();
    }

 
  if (bPrintList) info.Print();
  
  return theGenParticles;
}


//============================================================================
vector<GenParticle> TTTracks::GetGenParticles(int pdgId, bool isLastCopy)
//============================================================================
{

  // Initialise variables
  vector<GenParticle> myGenParticles;

  // For-lool: All GenParticles
  for (unsigned int genP_index = 0; genP_index < GenP_Pt->size(); genP_index++)
    {

      // Only examine particles of specific pdgId
      if ( abs(GenP_PdgId->at(genP_index) ) != pdgId) continue;

      // Get the genParticle
      GenParticle p = GetGenParticle(genP_index);
      SetGenParticleMomsAndDaus(p);
		
      if (!isLastCopy) myGenParticles.push_back(p);
      else
	{
	  
	  // Determine if it's a last copy
	  bool save = true;
	  bool bDecaysToSelf = false;
	  for (unsigned int i = 0; i < p.daughtersIndex().size(); i++) 
	    {

	      GenParticle d = GetGenParticle(p.daughtersIndex().at(i) );
	      if (abs(d.pdgId()) == pdgId)
		{
		  bDecaysToSelf = true;
		  break;
		}
	    }

	  if (bDecaysToSelf) continue;	  
	  myGenParticles.push_back(p);
	  
	  // Print properties?
	  if (0) p.PrintDaughters();
	  
	}
    }

  return myGenParticles;
}

      
//============================================================================
GenParticle TTTracks::GetGenParticle(unsigned int Index)
//============================================================================
{

  // Get the GenParticle properties
  double Pt      = GenP_Pt->at(Index);
  double Eta     = GenP_Eta->at(Index);
  double Phi     = GenP_Phi->at(Index);
  double Mass    = GenP_Mass->at(Index);
  int Charge     = GenP_Charge->at(Index);
  int PdgId      = GenP_PdgId->at(Index);
  int Status     = GenP_Status->at(Index);
  double VertexX = GenP_VertexX->at(Index);
  double VertexY = GenP_VertexY->at(Index);
  double VertexZ = GenP_VertexZ->at(Index);
  vector<unsigned short> MothersIndex   = GenP_Mothers->at(Index);
  vector<unsigned short> DaughtersIndex = GenP_Daughters->at(Index);

  // Construct the GenParticle
  GenParticle theGenParticle(Index, Pt, Eta, Phi, Mass, Charge, PdgId, Status, VertexX, VertexY, VertexZ, MothersIndex, DaughtersIndex);
  
  return theGenParticle;
}



//============================================================================
void TTTracks::SetGenParticleMomsAndDaus(GenParticle &p)
//============================================================================
{

  // Variable declaration
  vector<GenParticle> Mothers;
  vector<GenParticle> Daughters;

  // For-loops:
  for (unsigned int i = 0; i < p.mothersIndex().size();   i++) Mothers.push_back  ( GetGenParticle( p.mothersIndex().at(i)   ) );
  for (unsigned int j = 0; j < p.daughtersIndex().size(); j++) Daughters.push_back( GetGenParticle( p.daughtersIndex().at(j) ) );
  p.SetMothers(Mothers);
  p.SetDaughters(Daughters);
  SetGenParticleFinalDaughters(p);  
  
  return;
}


//============================================================================
void TTTracks::SetGenParticleFinalDaughters(GenParticle &p)
//============================================================================
{

  // Variable declaration
  vector<GenParticle> finalDaughters;
  vector<unsigned short> finalDaughtersIndex;

  // For taus, get hadronic decay products (pi+/-, pi0, K+/-, K0, K0L)
  GetHadronicTauFinalDaughters(p, finalDaughtersIndex);

  // For-loop: All daughter indices
  for (unsigned int i = 0; i < finalDaughtersIndex.size(); i++)
    {
      GenParticle d = GetGenParticle( finalDaughtersIndex.at(i) );
      SetGenParticleMomsAndDaus(d);
      finalDaughters.push_back(d);
    }
  
  // Need to manually set the final daughters (not the direct ones)
  p.SetFinalDaughters(finalDaughters);

  return;
}


//============================================================================
vector<TTTrack> TTTracks::GetTTTracks(const double minPt,
				      const double minEta,
				      const double maxEta,
				      const double maxChiSqRed,
				      const unsigned int minStubs,
				      const unsigned int minStubsPS,
				      const unsigned int maxStubsPS,
				      const unsigned nFitParams,
				      bool bPrintList)
//============================================================================
{
  vector<TTTrack> theTTTracks;

  // For-loop: All TTTracsk
  for (Size_t iTk = 0; iTk < L1Tks_Pt->size(); iTk++)
    { 
      TTTrack tk = GetTTTrack(iTk, nFitParams);

      // Apply selection criteria
      if (tk.getPt() < minPt) continue;
      if (abs(tk.getEta()) > maxEta) continue;
      if (abs(tk.getEta()) < minEta) continue;
      if (tk.getChi2Red() > maxChiSqRed) continue;
      if (tk.getNumOfStubs() < minStubs) continue;
      if (tk.getNumOfStubsPS() < minStubsPS) continue;
      if (tk.getNumOfStubsPS() > maxStubsPS) continue;
      // double z0 = tk.getZ0();
      // double d0 = tk.getD0();
      theTTTracks.push_back( tk );

    }

  if (bPrintList) PrintTTTrackCollection(theTTTracks);
  
  return theTTTracks;
}


//============================================================================
TTTrack TTTracks::GetTTTrack(unsigned int Index,
			     const unsigned int nFitParams)
//============================================================================
{
  
  // Initialise variables
  TVector3 p;  

  // Get the track properties
  double pt  = L1Tks_Pt->at(Index);
  double eta = L1Tks_Eta->at(Index);
  double phi = L1Tks_Phi->at(Index);
  double x0  = L1Tks_POCAx->at(Index); // 1e6 if nFitParams == 4 
  double y0  = L1Tks_POCAy->at(Index); // 1e6 if nFitParams == 4    
  double z0  = L1Tks_POCAz->at(Index);
  p.SetPtEtaPhi(pt, eta, phi);
  ROOT::Math::XYZVector aPOCA(x0, y0, z0);
  double aRInv        = L1Tks_RInv->at(Index);
  double aChi2        = L1Tks_ChiSquared->at(Index);
  double aStubPtCons  = L1Tks_StubPtConsistency->at(Index);
  bool isGenuine      = L1Tks_IsGenuine->at(Index);
  bool isUnknown      = L1Tks_IsUnknown->at(Index);
  bool isCombinatoric = L1Tks_IsCombinatoric->at(Index);
  bool isLoose        = L1Tks_IsLoose->at(Index);
  bool isFake         = L1Tks_IsFake->at(Index);
  int  nStubs         = L1Tks_NStubs->at(Index);
  int  nStubsPS       = L1Tks_NStubsPS->at(Index);
  int  nStubsBarrel   = L1Tks_NStubsBarrel->at(Index);
  int  nStubsEndcap   = L1Tks_NStubsEndcap->at(Index);
  vector<unsigned int> stubs_isPS   ; //= L1Tks_Stubs_isPS->at(Index);
  vector<unsigned int> stubs_iDisk  ; //= L1Tks_Stubs_iDisk->at(Index);
  vector<unsigned int> stubs_iLayer ; //= L1Tks_Stubs_iLayer->at(Index);
  vector<unsigned int> stubs_iPhi   ; //= L1Tks_Stubs_iPhi->at(Index);
  vector<unsigned int> stubs_iRing  ; //= L1Tks_Stubs_iRing->at(Index);
  vector<unsigned int> stubs_iSide  ; //= L1Tks_Stubs_iSide->at(Index);
  vector<unsigned int> stubs_iZ     ; //= L1Tks_Stubs_iZ->at(Index);
  int matchTP_index = L1Tks_TP_Index->at(Index);
  // auxTools_.PrintVector(stubs_isPS);
  
  // Construct the TTTrack
  TTTrack theTTTrack(Index, p, aPOCA,  aRInv, aChi2, aStubPtCons, isGenuine, isUnknown, isCombinatoric, isLoose, isFake,
   		     nStubs, nStubsPS, nStubsBarrel, nStubsEndcap, matchTP_index,
   		     stubs_isPS, stubs_iDisk, stubs_iLayer, stubs_iPhi, stubs_iRing, stubs_iSide, stubs_iZ, nFitParams);
  
  // // Get the uniquely matched TP
  // TrackingParticle theTP;
  // if (matchTP_index >= 0) theTP= GetTrackingParticle(matchTP_index);
  // theTTTrack.SetTP(theTP);
  
  // // if (DEBUG) theTTTrack.PrintProperties();
  // // if (DEBUG) theTP.PrintProperties();
  // theTTTrack.PrintProperties();
  // theTP.PrintProperties();
  
  return theTTTrack;
}



//============================================================================
vector<TrackingParticle> TTTracks::GetTrackingParticles(bool bPrintList)
//============================================================================
{
  vector<TrackingParticle> theTrackingParticles;

  // For-loop: Tracking Particles
  for (Size_t iTP = 0; iTP < TP_Pt->size(); iTP++)
    {
      TrackingParticle p = GetTrackingParticle(iTP);
      theTrackingParticles.push_back(p);
      
    }
  if (bPrintList) PrintTrackingParticleCollection(theTrackingParticles);
  
  return theTrackingParticles;
}


//============================================================================
TrackingParticle TTTracks::GetTrackingParticle(unsigned int Index)
//============================================================================
{

  // Get the track properties
  double pt            = TP_Pt->at(Index);
  double eta           = TP_Eta->at(Index);
  double phi           = TP_Phi->at(Index);
  int charge           = TP_Charge->at(Index);
  int pdgId            = TP_PdgId->at(Index);
  double d0_propagated = TP_d0_propagated->at(Index);
  double z0_propagated = TP_z0_propagated->at(Index);
  int nMatch           = TP_NMatch->at(Index);
  int ttTrackIndex     = TP_TTTrackIndex->at(Index);
  int ttClusters       = TP_NClusters->at(Index);
  int ttStubs          = TP_NStubs->at(Index);
  int ttTracks         = TP_NTracks->at(Index);
  double x0            = TP_x0_produced->at(Index);
  double y0            = TP_y0_produced->at(Index);
  double z0            = TP_z0_produced->at(Index);
  int eventId          = TP_EventId->at(Index);

  // Construct higher-level variables
  TVector3 p;
  p.SetPtEtaPhi(pt, eta, phi);
  ROOT::Math::XYZVector poca(x0, y0, z0);

  // Construct the TP
  TrackingParticle theTrackingParticle(Index, p, poca, d0_propagated, z0_propagated, charge, pdgId, nMatch, ttTrackIndex, ttClusters, ttStubs, ttTracks, eventId);

  // Get the uniquely matched TTTrack
  TTTrack theTrack;
  if (ttTrackIndex >= 0) theTrack= GetTTTrack(ttTrackIndex, 4);
  // theTrackingParticle.SetTTTTrack(theTrack); // cannot use! cyclic problems 

  if (DEBUG) theTrackingParticle.PrintProperties();
  // if (DEBUG) theTrack.PrintProperties(); // cannot use! cyclic problems 
  return theTrackingParticle;
}


//============================================================================
void TTTracks::PrintTrackingParticleCollection(vector<TrackingParticle> collection)
//============================================================================
{

  Table info("Index | Pt | Eta | Phi | PdgId | Q | x0 | y0 | z0 | d0 | d0-phi | NMatch | TTTrackIndex | TTClusters | TTStubs | TTTracks | Event-Id", "Text");

  // For-loop: Tracking Particles
  int row=0;
  for (vector<TrackingParticle>::iterator p = collection.begin(); p != collection.end(); p++)
    {
      // Construct table
      info.AddRowColumn(row, auxTools_.ToString( p->index()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Perp(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Eta() , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Phi() , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getPdgId()) );
      info.AddRowColumn(row, p->getQ());
      info.AddRowColumn(row, auxTools_.ToString( p->getX0()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getY0()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getZ0()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getD0()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getD0Phi(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getNMatch())       );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackIndex()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTClusters())   );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTStubs())      );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTracks())     );
      info.AddRowColumn(row, auxTools_.ToString( p->getEventId())     );
      row++;
    }
  

  info.Print();
  
  return;    
}


//============================================================================
void TTTracks::PrintTTTrackCollection(vector<TTTrack> collection)
//============================================================================
{

  Table info("Index | Pt | Eta | Phi | x0 | y0 | z0 | d0 | Q | Chi2 | DOF | Chi2Red | Stubs (PS) | StubPtCons. | Genuine | Unknown | Comb.", "Text");
      
  // For-loop: All TTTracsk
  int row=0;
  for (vector<TTTrack>::iterator p = collection.begin(); p != collection.end(); p++)    
    {
      // Construct table
      info.AddRowColumn(row, auxTools_.ToString( p->index()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Perp(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Eta() , 3) );  
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Phi() , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getX0(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getY0(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getZ0(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getD0(), 3) );  
      info.AddRowColumn(row, p->getQ());
      info.AddRowColumn(row, auxTools_.ToString( p->getChi2(),3 ));
      info.AddRowColumn(row, auxTools_.ToString( p->getDOF()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getChi2Red(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getNumOfStubs()) + " (" + auxTools_.ToString(p->getNumOfStubsPS()) + ")");
      info.AddRowColumn(row, auxTools_.ToString( p->getStubPtConsistency(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getIsGenuine()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getIsUnknown()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getIsCombinatoric()) );
      row++;
    }

  info.Print();

  return;
}

//============================================================================
void TTTracks::PrintGenParticleCollection(vector<GenParticle> collection)
//============================================================================
{

  for (vector<GenParticle>::iterator p = collection.begin(); p != collection.end(); p++)
    {
      std::cout << "\nGenParticle:" << endl;
      p->PrintProperties();
      std::cout << "\nDaughters:" << endl;
      p->PrintDaughters();
      std::cout << "\nFinal Daughters:" << endl;
      p->PrintFinalDaughters();
      std::cout << "\nFinal Daughters (Charged):" << endl;
      p->PrintFinalDaughtersCharged();
      std::cout << "\nFinal Daughters (Neutral):" << endl;
      p->PrintFinalDaughtersNeutral();
      std::cout << "\nMothers:" << endl;
      p->PrintMothers();
    }
  
  return;
}

#endif

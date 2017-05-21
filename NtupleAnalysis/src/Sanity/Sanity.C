#ifndef Sanity_cxx
#define Sanity_cxx

// User
#include "Sanity.h"
// C++
#include <stdexcept>

//============================================================================
void Sanity::InitVars_()
//============================================================================
{

  DEBUG = false;
  if (DEBUG) std::cout << "=== Sanity::InitVars_()" << std::endl;
  
  datasets_      = datasets_.GetDataset(mcSample);
  tk_Collection  =  "TTTracks"; // Default: Only "TTTracks" available (not "TTPixelTracks")
  tk_nFitParams  = 4;           // Default: 4
  tk_minPt       = 2.00;        // Default: 2.0
  tk_minEta      = 0.0;         // Default: 0.0
  tk_maxEta      = 1e6;         // Default: 1e6
  tk_maxChiSqRed = 1e6;         // Default: 1e6
  tk_minStubs    =   0;         // Default: 0
  tk_minStubsPS  =   0;         // Default: 0
  tk_maxStubsPS  = 1e6;         // Default: 1e6

  return;
}


//============================================================================
void Sanity::PrintSettings(void)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::PrintSettings()" << std::endl;
  
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

  settings.AddRowColumn(10, "" );
  
  settings.Print();
  
  return;
}

//============================================================================
void Sanity::Loop()
//============================================================================
{

  if (1) std::cout << "=== Sanity::Loop()" << std::endl;
  
  // Sanity check
  if (fChain == 0) return;
  cout << "\tGetting Tree Entries ..." << endl;
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  // Initialisations
  InitVars_();
  BookHistos_();
  Long64_t nbytes       = 0;
  Long64_t nb           = 0;
  unsigned int nAllEvts = fChain->GetEntries();
  unsigned int nEvts    = 0;

  // Print User Settings
  cout << "\tPrinting user settings:" << endl;
  PrintSettings();

  cout << "\tAnalyzing " << nEntries << "/" << fChain->GetEntries() << " events" << endl;
  // For-loop: All TTree Entries
  for (int jentry = 0; jentry < nEntries; jentry++, nEvts++)
    {

      if(DEBUG) cout << "\tEntry = " << jentry << endl;

      // Load the tree && Get the entry
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;


      // ======================================================================
      // GenParticles Collection
      // ======================================================================
      if (DEBUG) cout << "=== GenParticles (" << GenP_Pt->size() << ")" << endl;
      vector<GenParticle> GenParticles = GetGenParticles(false);
      if (DEBUG) PrintGenParticleCollection(GenParticles);
      
      // For-loop: GenParticles
      for (vector<GenParticle>::iterator p = GenParticles.begin(); p != GenParticles.end(); p++)
	  {
	    
	    // Fill Histograms
	    hGenP_Index->Fill(p->index());
	    hGenP_Pt->Fill(p->pt());
	    hGenP_Eta->Fill(p->eta());
	    hGenP_Phi->Fill(p->phi());
	    hGenP_PtVis->Fill(p->p4vis().Pt());
	    hGenP_EtaVis->Fill(p->p4vis().Eta());
	    hGenP_PhiVis->Fill(p->p4vis().Phi());
	    hGenP_Mass->Fill(p->mass());
	    hGenP_PdgId->Fill(p->pdgId());
	    hGenP_Charge->Fill(p->charge());
	    hGenP_Status->Fill(p->status());
	    hGenP_VertexX->Fill(p->vx());
	    hGenP_VertexY->Fill(p->vy());
	    hGenP_VertexZ->Fill(p->vz());
	    hGenP_Mothers->Fill(p->mothers().size());
	    hGenP_Daughters->Fill(p->daughters().size());
	    hGenP_FinalDaughters->Fill(p->finalDaughters().size());
	    hGenP_FinalDaughtersCharged->Fill(p->finalDaughtersCharged().size());
	    hGenP_FinalDaughtersNeutral->Fill(p->finalDaughtersNeutral().size());
	    
	  }

  
      
      // ======================================================================
      // Sanity Particle Collections
      // ======================================================================
      if (DEBUG) cout << "=== Tracking Particles (" << TP_Pt->size() << ")" << endl;
      vector<TrackingParticle> TPs  = GetTrackingParticles(false);
      sort( TPs.begin(), TPs.end(), PtComparatorTP() ); // not sorted by default
      if (DEBUG) PrintTrackingParticleCollection(TPs);

      // For-loop: GenParticles
      for (vector<TrackingParticle>::iterator p = TPs.begin(); p != TPs.end(); p++)
	  {

	    // Fill histograms
	    hTP_Index->Fill(p->index());
	    hTP_Pt->Fill(p->getPt());
	    hTP_Eta->Fill(p->getEta());
	    hTP_Phi->Fill(p->getPhi());
	    hTP_PdgId->Fill(p->getPdgId());
	    hTP_X0->Fill(p->getX0());
	    hTP_Y0->Fill(p->getY0());
	    hTP_Z0->Fill(p->getZ0());
	    hTP_Charge->Fill(p->getCharge());
	    hTP_Dxy->Fill(p->getDxy());
	    hTP_D0propagated->Fill(p->getD0propagated());
	    hTP_Z0propagated->Fill(p->getZ0propagated());
	    hTP_D0->Fill(p->getD0());
	    hTP_D0Sign->Fill(p->getD0Sign());
	    hTP_D0Phi->Fill(p->getD0Phi());
	    hTP_EventId->Fill(p->getEventId());
	    hTP_NMatch->Fill(p->getNMatch());
	    hTP_TTTrackIndex->Fill(p->getTTTrackIndex());
	    hTP_TTClusters->Fill(p->getTTClusters());
	    hTP_TTStubs->Fill(p->getTTStubs());
	    hTP_TTTracks->Fill(p->getTTTracks());

	    // p->PrintProperties(false);
	    // p->PrintAllProperties();
	    
	  }



      // ======================================================================
      // TTTracks Collections
      // ======================================================================
      if (DEBUG) cout << "=== TTracks (" << L1Tks_Pt->size() << ")" << endl;
      vector<TTTrack> TTTracks = GetTTTracks(tk_minPt, tk_minEta, tk_maxEta, tk_maxChiSqRed, tk_minStubs, tk_minStubsPS, tk_maxStubsPS, tk_nFitParams, false);
      sort( TTTracks.begin(), TTTracks.end(), PtComparatorTTTrack() ); // not sorted by default
      if (DEBUG) PrintTTTrackCollection(TTTracks);

      // For-loop: TTTracks
      for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
	  {
	    // TVector3 tk_p3   = tk->getMomentum();
	    double tk_pt     = tk->getPt();
	    double tk_eta    = tk->getEta();
	    double tk_phi    = tk->getPhi();
	    int tk_charge    = tk->getCharge();
	    double tk_Chi2   = tk->getChi2();
	    double tk_Chi2Red= tk->getChi2Red();
	    bool tk_IsLoose  = tk->getIsLoose();
	    bool tk_IsFake   = tk->getIsFake();
	    bool tk_TPIndex  = tk->getTPIndex();
	    double tk_getD0  = tk->getD0();
	    int tk_DOF       = tk->getDOF();
	    double tk_StubPtCons   = tk->getStubPtConsistency();
	    bool tk_IsGenuine      = tk->getIsGenuine();
	    bool tk_IsUnknown      = tk->getIsUnknown();
	    bool tk_IsCombinatoric = tk->getIsCombinatoric();

	    unsigned int tk_NStubs       = tk->getNumOfStubs();
	    unsigned int tk_NStubsPS     = tk->getNumOfStubsPS();
	    unsigned int tk_NStubsBarrel = tk->getNumOfBarrelStubs();
	    unsigned int tk_NStubsEndcap = tk->getNumOfEndcapStubs();

	    double tk_X0   = tk->getX0();
	    double tk_Y0   = tk->getY0();
	    double tk_Z0   = tk->getZ0();
	    double tk_RInv = tk->getRInv();
	    // ROOT::Math::XYZVector tk_POCA = tk->getPOCA();
	    
	    if (0) tk->PrintProperties();
	    if (0) tk->PrintAllProperties();

	    // Fill Histograms
	    hL1Tks_Pt ->Fill(tk_pt);
	    hL1Tks_Eta->Fill(tk_eta);
	    hL1Tks_Phi->Fill(tk_phi);
	    hL1Tks_Charge->Fill(tk_charge);
	    hL1Tks_POCAx->Fill(tk_X0);
	    hL1Tks_POCAy->Fill(tk_Y0);
	    hL1Tks_POCAz->Fill(tk_Z0);
	    hL1Tks_ChiSq->Fill(tk_Chi2);
	    hL1Tks_ChiSqRed->Fill(tk_Chi2Red);
	    hL1Tks_StubPtConsistency->Fill(tk_StubPtCons);
	    hL1Tks_RInv->Fill(tk_RInv);
	    hL1Tks_IsGenuine->Fill(tk_IsGenuine);
	    hL1Tks_IsUnknown->Fill(tk_IsUnknown);
	    hL1Tks_IsCombinatoric->Fill(tk_IsCombinatoric);
	    hL1Tks_IsLoose->Fill(tk_IsLoose);
	    hL1Tks_IsFake->Fill(tk_IsLoose);
	    hL1Tks_NStubs->Fill(tk_NStubs);
	    hL1Tks_NStubsPS->Fill(tk_NStubsPS);
	    hL1Tks_NStubsBarrel->Fill(tk_NStubsBarrel);
	    hL1Tks_NStubsEndcap->Fill(tk_NStubsEndcap);
	    hL1Tks_TP_Index->Fill(tk_TPIndex);

	  }


      // Tau Collection
      if (DEBUG) cout << "=== Taus (" << tauEt->size() << ")" << endl;
      // vector<GenParticle> GenParticles = GetGenParticles(false);
      // if (DEBUG) PrintGenParticleCollection(GenParticles);

      ////////////////////////////////////////////////
      // Fill histograms
      ////////////////////////////////////////////////
      hHepMCEvt_VtxX_VtxY->Fill(HepMCEvt_VtxX, HepMCEvt_VtxY);
      hHepMCEvt_VtxZ->Fill(HepMCEvt_VtxZ);
      
    }//eof: Entries

  
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
void Sanity::BookHistos_(void)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::BookHistos_()" << std::endl;
  
  // Event-Type Histograms
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 2, 0.0, +2.0);
  histoTools_.BookHisto_1D(hHepMCEvt_VtxZ      , "HepMCEvt_VtxZ"     ,  "", 600, -30.0,  +30.0);
  histoTools_.BookHisto_2D(hHepMCEvt_VtxX_VtxY , "HepMCEvt_VtxX_VtxY",  "", 400,  -0.01,  +0.01, 400,  -0.01,  +0.01);

  // GenParticles
  histoTools_.BookHisto_1D(hGenP_Index, "GenP_Index", "; index; Entries", 501,  -0.5,  +500.5);
  histoTools_.BookHisto_1D(hGenP_Pt, "GenP_Pt", "; p_{T} (GeV/c); Entries", 200,  +0.0,  +200.0);
  histoTools_.BookHisto_1D(hGenP_Eta, "GenP_Eta", "; #eta; Entries", 130,  -2.6,  +2.6);
  histoTools_.BookHisto_1D(hGenP_Phi, "GenP_Phi", "; #phi (rads); Entries", 128,  -3.2,  +3.2);
  histoTools_.BookHisto_1D(hGenP_PtVis, "GenP_PtVis", "; p_{T}^{vis} (GeV/c); Entries", 200,  +0.0,  +200.0);
  histoTools_.BookHisto_1D(hGenP_EtaVis, "GenP_EtaVis", "; #eta^{vis}; Entries", 130,  -2.6,  +2.6);
  histoTools_.BookHisto_1D(hGenP_PhiVis, "GenP_PhiVis", "; #phi^{vis} (rads); Entries", 128,  -3.2,  +3.2);
  histoTools_.BookHisto_1D(hGenP_Mass, "GenP_Mass", "; mass (GeV/c^{2}); Entries", 200,  0.0,  +200.0);
  histoTools_.BookHisto_1D(hGenP_PdgId, "GenP_PdgId", "; pdgId; Entries", 1000,  -500.0,  +500.0);
  histoTools_.BookHisto_1D(hGenP_Charge, "GenP_Charge", "; charge (e); Entries", 5,  -2.5,  +2.5);
  histoTools_.BookHisto_1D(hGenP_Status, "GenP_Status", "; status (e); Entries", 5,   -0.5,  +4.5);
  histoTools_.BookHisto_1D(hGenP_VertexX, "GenP_VertexX", "; vtx x (cm); Entries", 400,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hGenP_VertexY, "GenP_VertexY", "; vtx y (cm); Entries", 400,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hGenP_VertexZ, "GenP_VertexZ", "; vtx z (cm); Entries", 200,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hGenP_Mothers, "GenP_Mothers", "; mother multiplicity; Entries", 10,  -0.5,  +9.5);
  histoTools_.BookHisto_1D(hGenP_Daughters, "GenP_Daughters", "; daughter multiplicity; Entries", 20,  -0.5,  +19.5);
  histoTools_.BookHisto_1D(hGenP_FinalDaughters, "GenP_FinalDaughters", "; final daughter multiplicity; Entries", 20,  -0.5,  +19.5);
  histoTools_.BookHisto_1D(hGenP_FinalDaughtersCharged, "GenP_FinalDaughtersCharged", "; final charged daughter multiplicity; Entries", 20,  -0.5,  +19.5);
  histoTools_.BookHisto_1D(hGenP_FinalDaughtersNeutral, "GenP_FinalDaughtersNeutral", "; final neutral daughter multiplicity; Entries", 20,  -0.5,  +19.5);


  // TrackingParticles
  histoTools_.BookHisto_1D(hTP_Index, "TP_Index", "; index; Entries", 501,  -0.5,  +500.5);
  histoTools_.BookHisto_1D(hTP_Pt, "TP_Pt", "; p_{T} (GeV/c); Entries", 200,  +0.0,  +200.0);
  histoTools_.BookHisto_1D(hTP_Eta, "TP_Eta", "; #eta; Entries", 130,  -2.6,  +2.6);
  histoTools_.BookHisto_1D(hTP_Phi, "TP_Phi", "; #phi (rads); Entries", 128,  -3.2,  +3.2);
  histoTools_.BookHisto_1D(hTP_PdgId, "TP_PdgId", "; pdgId; Entries", 1000,  -500.0,  +500.0);
  histoTools_.BookHisto_1D(hTP_X0, "TP_X0", "; x_{0} produced (cm); Entries", 400,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hTP_Y0, "TP_Y0", "; y_{0} produced (cm); Entries", 400,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hTP_Z0, "TP_Z0", "; z_{0} produced (cm); Entries", 400,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hTP_Charge, "TP_Charge", "; charge (e); Entries", 5,  -2.5,  +2.5);
  histoTools_.BookHisto_1D(hTP_Dxy, "TP_Dxy", "; d_{xy} (cm); Entries", 40, 0.0,  +20.0);
  histoTools_.BookHisto_1D(hTP_D0propagated, "TP_D0propagated", "; d_{0} propagated (cm); Entries", 40,  -20.0,  +20.0);
  histoTools_.BookHisto_1D(hTP_Z0propagated, "TP_Z0propagated", "; z_{0} propagated (cm); Entries", 40,  -20.0,  +20.0);
  histoTools_.BookHisto_1D(hTP_D0, "TP_D0", "; d_{0} (cm); Entries", 80,  -20.0,  +20.0);
    histoTools_.BookHisto_1D(hTP_D0Sign, "TP_D0Sign", "; d_{0} sign; Entries", 3,  -1.5,  +1.5);
    histoTools_.BookHisto_1D(hTP_D0Phi, "TP_D0Phi", "; d_{0} #phi (rads); Entries", 128,  -3.2,  +3.2);
    histoTools_.BookHisto_1D(hTP_EventId, "TP_EventId", "; Event Id; Entries", 1000,  -0.5,  +1000.5);
    histoTools_.BookHisto_1D(hTP_NMatch, "TP_NMatch", "number of tk-matches; Entries", 100,  -0.5,  +100.5);
    histoTools_.BookHisto_1D(hTP_TTTrackIndex, "TP_TTTrackIndex", "; index; Entries", 501,  -0.5,  +500.5);
    histoTools_.BookHisto_1D(hTP_TTClusters, "TP_TTClusters", "; TTClusters matched; Entries", 51,  -0.5,  +50.5);
    histoTools_.BookHisto_1D(hTP_TTStubs, "TP_TTStubs", "; TTStubs matched; Entries", 51,  -0.5,  +50.5);
    histoTools_.BookHisto_1D(hTP_TTTracks, "TP_TTTracks", "; TTStubs matched; Entries", 11,  -0.5,  +10.5);
  
  // TTracks
  histoTools_.BookHisto_1D(hL1Tks_Pt, "L1Tks_Pt", "; p_{T} (GeV/c); Entries", 200,  +0.0,  +200.0);
  histoTools_.BookHisto_1D(hL1Tks_Eta, "L1Tks_Eta", "; #eta; Entries", 130,  -2.6,  +2.6);
  histoTools_.BookHisto_1D(hL1Tks_Phi, "L1Tks_Phi", "; #phi (rads); Entries", 128,  -3.2,  +3.2);
  histoTools_.BookHisto_1D(hL1Tks_Charge, "L1Tks_Charge", "; charge (e); Entries", 6,  -2.0,  +2.0);
  histoTools_.BookHisto_1D(hL1Tks_POCAx, "L1Tks_POCAx", "; x_{0} (cm); Entries", 400,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hL1Tks_POCAy, "L1Tks_POCAy", "; y_{0} (cm); Entries", 400,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hL1Tks_POCAz, "L1Tks_POCAz", "; z_{0} (cm); Entries", 200,  -20.0,  +20.0); //units?
  histoTools_.BookHisto_1D(hL1Tks_ChiSq, "L1Tks_ChiSq", "; #chi^{2}; Entries", 200, 0.0, 200.0);
  histoTools_.BookHisto_1D(hL1Tks_ChiSqRed, "L1Tks_ChiSqRed", "; #chi^{2}/dof; Entries", 100, 0.0, 100.0);
  histoTools_.BookHisto_1D(hL1Tks_StubPtConsistency, "L1Tks_StubPtConsistency", "; p_{T}^{stub} constistency; Entries", 50, 0.0, 50.0);  
  histoTools_.BookHisto_1D(hL1Tks_RInv, "L1Tks_RInv", "; #rho^{-1} (m^{-1}); Entries", 200, -0.01, 0.01);
  histoTools_.BookHisto_1D(hL1Tks_IsGenuine, "L1Tks_IsGenuine", "; is genuine; Entries", 2, -0.5, 1.5);
  histoTools_.BookHisto_1D(hL1Tks_IsUnknown, "L1Tks_IsUnknown", "; is unknown; Entries", 2, -0.5, 1.5); 
  histoTools_.BookHisto_1D(hL1Tks_IsCombinatoric, "L1Tks_Combinatoric", "; is combinatoric; Entries", 2, -0.5, 1.5);
  histoTools_.BookHisto_1D(hL1Tks_IsLoose, "L1Tks_IsLoose", "; is loose; Entries", 2, -0.5, 1.5);
  histoTools_.BookHisto_1D(hL1Tks_IsFake, "L1Tks_IsFake", "; is fake; Entries", 2, -0.5, 1.5);
  histoTools_.BookHisto_1D(hL1Tks_NStubs, "L1Tks_NStubs", "; stub multiplicity", 13, -0.5, 12.5);
  histoTools_.BookHisto_1D(hL1Tks_NStubsPS, "L1Tks_NStubsPS", "; stub multiplicity (PS)", 13, -0.5, 12.5);
  histoTools_.BookHisto_1D(hL1Tks_NStubsBarrel, "L1Tks_NStubsBarrel", "; stub multiplicity (B)", 13, -0.5, 12.5);
  histoTools_.BookHisto_1D(hL1Tks_NStubsEndcap, "L1Tks_NStubsEndcap", "; stub multiplicity (E)", 13, -0.5, 12.5);
  histoTools_.BookHisto_1D(hL1Tks_TP_Index, "L1Tks_TP_Index", "; index of matched TP", 301, -0.5, 300.5);

  return;
}



//============================================================================
void Sanity::GetHadronicTauFinalDaughters(GenParticle p,
					    vector<unsigned short> &Daug)
//============================================================================
{

  if (DEBUG) std::cout << "=== Sanity::GetHadronicTauFinalDaughters()" << std::endl;
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
vector<GenParticle> Sanity::GetGenParticles(bool bPrintList)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::GetGenParticles()" << std::endl;
    
  vector<GenParticle> theGenParticles;
  Table info("Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | Mothers | Daughters |", "Text"); //LaTeX or Text    

  // For-loop: All GenParticles
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
vector<GenParticle> Sanity::GetGenParticles(int pdgId, bool isLastCopy)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::GetGenParticles()" << std::endl;
  
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
GenParticle Sanity::GetGenParticle(unsigned int Index)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::GetGenParticle()" << std::endl;

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
void Sanity::SetGenParticleMomsAndDaus(GenParticle &p)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::SetGenParticleMomsAndDaus()" << std::endl;
  
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
void Sanity::SetGenParticleFinalDaughters(GenParticle &p)
//============================================================================
{

  if (DEBUG) std::cout << "=== Sanity::SetGenParticleFinalDaughters()" << std::endl;
  
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
vector<TTTrack> Sanity::GetTTTracks(const double minPt,
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
  if (DEBUG) std::cout << "=== Sanity::GetTTTracks()" << std::endl;

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
TTTrack Sanity::GetTTTrack(unsigned int Index,
			     const unsigned int nFitParams)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::GetTTTrack()" << std::endl;
  
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
  int matchTP_index = L1Tks_TP_Index->at(Index);
  // auxTools_.PrintVector(stubs_isPS);
  
  // Construct the TTTrack
  TTTrack theTTTrack(Index, p, aPOCA,  aRInv, aChi2, aStubPtCons, isGenuine, isUnknown, isCombinatoric,
		     isLoose, isFake, nStubs, nStubsPS, nStubsBarrel, nStubsEndcap, matchTP_index, nFitParams);
  
  // Get the uniquely matched TP
  TrackingParticle theTP;
  if (matchTP_index >= 0) theTP= GetTrackingParticle(matchTP_index);
  // theTTTrack.SetTP(theTP); // cannot implement. cyclic
  
  if (DEBUG) theTTTrack.PrintProperties();
  if (DEBUG && matchTP_index >= 0)
    {
      theTP.PrintProperties();
      cout << "" << endl;
    }
  
  return theTTTrack;
}



//============================================================================
vector<TrackingParticle> Sanity::GetTrackingParticles(bool bPrintList)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::GetTrackingParticles()" << std::endl;
  
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
TrackingParticle Sanity::GetTrackingParticle(unsigned int Index)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::GetTrackingParticle()" << std::endl;
  
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
  if (DEBUG) theTrackingParticle.PrintProperties();

  // Get the uniquely matched TTTrack
  // TTTrack theTrack;
  // if (ttTrackIndex >= 0) theTrack= GetTTTrack(ttTrackIndex, 4);
  // theTrackingParticle.SetTTTTrack(theTrack); // cannot use! cyclic problems 
  // if (DEBUG) theTrack.PrintProperties(); // cannot use! cyclic problems
  
  return theTrackingParticle;
}


//============================================================================
void Sanity::PrintTrackingParticleCollection(vector<TrackingParticle> collection)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::PrintTrackingParticleCollection()" << std::endl;
  
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
void Sanity::PrintTTTrackCollection(vector<TTTrack> collection)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::PrintTTTrackCollection()" << std::endl;

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
void Sanity::PrintGenParticleCollection(vector<GenParticle> collection)
//============================================================================
{
  if (DEBUG) std::cout << "=== Sanity::PrintGenParticleCollection()" << std::endl;
  
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

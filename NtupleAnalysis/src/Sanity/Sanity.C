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

  cfg_DEBUG = true;
  datasets_ = datasets_.GetDataset(mcSample);
  if (cfg_DEBUG) std::cout << "=== Sanity::InitVars_()" << std::endl;
  
  cfg_AddGenP    = true;
  cfg_AddL1Tks   = true;
  cfg_AddStubs   = false;
  cfg_AddTPs     = true;
  cfg_AddEGs     = true;
  cfg_AddTaus    = true;
  cfg_AddJets    = true;
  cfg_AddMuons   = false;
  cfg_AddSums    = true;

  cfg_tk_Collection  =  "TTTracks"; // Default: "TTTracks" (not "TTPixelTracks")
  cfg_tk_nFitParams  = 4;           // Default: 4
  cfg_tk_minPt       = 2.00;        // Default: 2.0
  cfg_tk_minEta      = 0.0;         // Default: 0.0
  cfg_tk_maxEta      = 1e6;         // Default: 1e6
  cfg_tk_maxChiSqRed = 1e6;         // Default: 1e6
  cfg_tk_minStubs    =   0;         // Default: 0
  cfg_tk_minStubsPS  =   0;         // Default: 0
  cfg_tk_maxStubsPS  = 1e6;         // Default: 1e6
  
  PrintSettings();

  return;
}


//============================================================================
void Sanity::PrintSettings(void)
//============================================================================
{
  if (cfg_DEBUG) std::cout << "=== Sanity::PrintSettings()" << std::endl;
  
  // Inform user of settings
  Table settings("Variable | Cut | Value | Default | Units", "Text");
  settings.AddRowColumn(0, "MC Sample");
  settings.AddRowColumn(0, "==" );
  settings.AddRowColumn(0, mcSample );
  settings.AddRowColumn(0, "");
 
  settings.AddRowColumn(1, "Add: Gen Particles");
  settings.AddRowColumn(1, "==" );
  settings.AddRowColumn(1, auxTools_.ToString(cfg_AddGenP ) );
  settings.AddRowColumn(1, "true");
  
  settings.AddRowColumn(2, "Add: TTTracks");
  settings.AddRowColumn(2, "==" );
  settings.AddRowColumn(2, auxTools_.ToString(cfg_AddL1Tks) );
  settings.AddRowColumn(2, "true");
  
  settings.AddRowColumn(3, "Add: TTStubs");
  settings.AddRowColumn(3, "==" );
  settings.AddRowColumn(3, auxTools_.ToString(cfg_AddStubs) );
  settings.AddRowColumn(3, "false");

  settings.AddRowColumn(4, "Add: Tracking Particles");
  settings.AddRowColumn(4, "==" );
  settings.AddRowColumn(4, auxTools_.ToString(cfg_AddTPs) );
  settings.AddRowColumn(4, "true");

  settings.AddRowColumn(5, "Add: EGs");
  settings.AddRowColumn(5, "==" );
  settings.AddRowColumn(5, auxTools_.ToString(cfg_AddEGs) );
  settings.AddRowColumn(5, "true");

  settings.AddRowColumn(6, "Add: L1 Taus");
  settings.AddRowColumn(6, "==" );
  settings.AddRowColumn(6, auxTools_.ToString(cfg_AddTaus) );
  settings.AddRowColumn(6, "true");

  settings.AddRowColumn(7, "Add: L1 Jets");
  settings.AddRowColumn(7, "==" );
  settings.AddRowColumn(7, auxTools_.ToString(cfg_AddJets) );
  settings.AddRowColumn(7, "true");
  
  settings.AddRowColumn(8, "Add: L1 Muons");
  settings.AddRowColumn(8, "==" );
  settings.AddRowColumn(8, auxTools_.ToString(cfg_AddMuons) );
  settings.AddRowColumn(8, "false");
  
  settings.AddRowColumn(9, "Add: L1 Sums");
  settings.AddRowColumn(9, "==" );
  settings.AddRowColumn(9, auxTools_.ToString(cfg_AddSums) );
  settings.AddRowColumn(9, "true");
  
  settings.AddRowColumn(11, "Tracks: Collection");
  settings.AddRowColumn(11, "==");
  settings.AddRowColumn(11, cfg_tk_Collection);
  settings.AddRowColumn(11, "TTTracks");
  settings.AddRowColumn(11, "");
  
  settings.AddRowColumn(12, "Tracks: Fit Parameters");
  settings.AddRowColumn(12, "==");
  settings.AddRowColumn(12, auxTools_.ToString( cfg_tk_nFitParams) );
  settings.AddRowColumn(12, "4");
  settings.AddRowColumn(12, "");

  settings.AddRowColumn(13, "Tracks: Pt");
  settings.AddRowColumn(13, ">=");
  settings.AddRowColumn(13, auxTools_.ToString( cfg_tk_minPt) );
  settings.AddRowColumn(13, "2" );
  settings.AddRowColumn(13, "GeV/c" );
  
  settings.AddRowColumn(14, "Tracks: |Eta|");
  settings.AddRowColumn(14, ">=");
  settings.AddRowColumn(14, auxTools_.ToString( cfg_tk_minEta) );
  settings.AddRowColumn(14, "0.0" );
  settings.AddRowColumn(14, "" );

  settings.AddRowColumn(15, "Tracks: |Eta|");
  settings.AddRowColumn(15, "<=");
  settings.AddRowColumn(15, auxTools_.ToString(cfg_tk_maxEta) );
  settings.AddRowColumn(15, auxTools_.ToString(2.5) );
  settings.AddRowColumn(15, "" );
  
  settings.AddRowColumn(16, "Tracks: ChiSqRed");
  settings.AddRowColumn(16, "<=");
  settings.AddRowColumn(16, auxTools_.ToString( cfg_tk_maxChiSqRed) );
  settings.AddRowColumn(16, "200/DOF");
  settings.AddRowColumn(16, "");

  settings.AddRowColumn(17, "Tracks: Stubs");
  settings.AddRowColumn(17, ">=");
  settings.AddRowColumn(17, auxTools_.ToString(cfg_tk_minStubs) );
  settings.AddRowColumn(17, auxTools_.ToString(4) );
  settings.AddRowColumn(17, "" );
  
  settings.AddRowColumn(18, "Tracks: PS-Stubs (min)");
  settings.AddRowColumn(18, ">=");
  settings.AddRowColumn(18, auxTools_.ToString(cfg_tk_minStubsPS) );
  settings.AddRowColumn(18, auxTools_.ToString(3) );
  settings.AddRowColumn(18, "" );

  settings.AddRowColumn(19, "Tracks: PS-Stubs (max)");
  settings.AddRowColumn(19, "<=");
  settings.AddRowColumn(19, auxTools_.ToString(cfg_tk_maxStubsPS) );
  settings.AddRowColumn(19, auxTools_.ToString(6) );
  settings.AddRowColumn(19, "" );
 
  settings.AddRowColumn(20, "" );

  settings.AddRowColumn(21, "cfg_DEBUG");
  settings.AddRowColumn(21, "==");
  settings.AddRowColumn(21, auxTools_.ToString(cfg_DEBUG) );
  settings.AddRowColumn(21, "false" );
  settings.AddRowColumn(21, "" );

  
  settings.Print();
  
  return;
}

//============================================================================
void Sanity::Loop()
//============================================================================
{

  // Sanity check
  if (fChain == 0) return;
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);

  cout << "=== Sanity:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;
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

      if(cfg_DEBUG) cout << "\tEntry = " << jentry << endl;

      // Load the tree && Get the entry
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;


      // ======================================================================
      // GenParticles Collection
      // ======================================================================
      if (cfg_AddGenP)
	{
	  if (cfg_DEBUG) cout << "\n=== GenParticles (" << GenP_Pt->size() << ")" << endl;
	  vector<GenParticle> GenParticles = GetGenParticles(false);
	  if (cfg_DEBUG) PrintGenParticleCollection(GenParticles);
      
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
	}
  
      
      // ======================================================================
      // Tracking Particle Collections
      // ======================================================================
      if (cfg_AddTPs)
	{
	  if (cfg_DEBUG) cout << "\n=== Tracking Particles (" << TP_Pt->size() << ")" << endl;
	  vector<TrackingParticle> TPs  = GetTrackingParticles(false);
	  sort( TPs.begin(), TPs.end(), PtComparatorTP() ); // not sorted by default
	  if (cfg_DEBUG) PrintTrackingParticleCollection(TPs);
	  
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
	}


      // ======================================================================
      // TTTracks Collections
      // ======================================================================
      if (cfg_AddL1Tks)
	{

	  if (cfg_DEBUG) cout << "\n=== TTracks (" << L1Tks_Pt->size() << ")" << endl;
	  vector<TTTrack> TTTracks = GetTTTracks(cfg_tk_minPt, cfg_tk_minEta, cfg_tk_maxEta, cfg_tk_maxChiSqRed, cfg_tk_minStubs, cfg_tk_minStubsPS, cfg_tk_maxStubsPS, cfg_tk_nFitParams, false);
	  sort( TTTracks.begin(), TTTracks.end(), PtComparatorTTTrack() ); // not sorted by default
	  if (cfg_DEBUG) PrintTTTrackCollection(TTTracks);
	  
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
	}


      // ======================================================================
      // L1Tau Collection
      // ======================================================================
      if (cfg_AddTaus)
	{

	  if (cfg_DEBUG) cout << "\n=== L1Taus (" << L1Tau_Et->size() << ")" << endl;
	  vector<L1Tau> L1Taus = GetL1Taus(cfg_DEBUG);
	  
	  hL1Tau_nTaus->Fill(L1Tau_Et->size()); //L1Tau_nTaus);
	  // For-loop: L1 Taus
	  int index = 0;
	  for (auto tau = L1Taus.begin(); tau != L1Taus.end(); tau++, index++)
	    {
	      
	      double tau_et           = tau->getEt();
	      double tau_eta          = tau->getEta();
	      double tau_phi          = tau->getPhi();
	      short int tau_IET       = tau->getIET();
	      short int tau_IEta      = tau->getIEta();
	      short int tau_IPhi      = tau->getIPhi();
	      short int tau_Iso       = tau->getIso();
	      short int tau_Bx        = tau->getBx();
	      short int tau_TowerIPhi = 0; //tau->getTowerIPhi();
	      short int tau_TowerIEta = 0; //tau->getTowerIEta();
	      short int tau_RawEt     = tau->getRawEt();
	      short int tau_IsoEt     = tau->getIsoEt();
	      short int tau_NTT       = tau->getNTT();
	      short int tau_HasEM     = tau->getHasEM();
	      short int tau_IsMerged  = tau->getIsMerged();
	      short int tau_HwQual    = tau->getHwQual();
	      
	      // Fill Histograms
	      hL1Tau_Index->Fill(index);
	      hL1Tau_Et->Fill(tau_et);
	      hL1Tau_Eta->Fill(tau_eta);
	      hL1Tau_Phi->Fill(tau_phi);
	      hL1Tau_IET->Fill(tau_IET);
	      hL1Tau_IEta->Fill(tau_IEta);
	      hL1Tau_IPhi->Fill(tau_IPhi);
	      hL1Tau_Iso->Fill(tau_Iso);
	      hL1Tau_Bx->Fill(tau_Bx);
	      hL1Tau_TowerIPhi->Fill(tau_TowerIPhi);
	      hL1Tau_TowerIEta->Fill(tau_TowerIEta);
	      hL1Tau_RawEt->Fill(tau_RawEt);
	      hL1Tau_IsoEt->Fill(tau_IsoEt);
	      hL1Tau_NTT->Fill(tau_NTT);
	      hL1Tau_HasEM->Fill(tau_HasEM);
	      hL1Tau_IsMerged->Fill(tau_IsMerged);
	      hL1Tau_HwQual->Fill(tau_HwQual);
	  }
	}

      
      // ======================================================================
      // L1Jet Collection
      // ======================================================================
      if (cfg_AddJets)
	{
	  if (cfg_DEBUG) cout << "\n=== L1Jets (" << L1Jet_Et->size() << ")" << endl;
	  vector<L1Jet> L1Jets = GetL1Jets(cfg_DEBUG);
	  
	  hL1Jet_nJets->Fill(L1Jet_Et->size()); //L1Jet_nJets);      
	  // For-loop: L1 Jets
	  for (auto jet = L1Jets.begin(); jet != L1Jets.end(); jet++)
	    {
	      short Index      = jet->getIndex();
	      float Et         = jet->getEt();
	      float Eta        = jet->getEta();
	      float Phi        = jet->getPhi();
	      short IET        = jet->getIET();
	      short IEta       = jet->getIEta();
	      short IPhi       = jet->getIPhi();
	      short Bx         = jet->getBx();
	      short RawEt      = jet->getRawEt();
	      short SeedEt     = jet->getSeedEt();
	      short TowerIEta  = 0; //jet->getTowerIEta();
	      short TowerIPhi  = 0; //jet->getTowerIPhi();
	      short PUEt       = jet->getPUEt();
	      short PUDonutEt0 = jet->getPUDonutEt0();
	      short PUDonutEt1 = jet->getPUDonutEt1();
	      short PUDonutEt2 = jet->getPUDonutEt2();
	      short PUDonutEt3 = jet->getPUDonutEt3();
	      
	      // Fill Histograms
	      hL1Jet_Index->Fill(Index);
	      hL1Jet_Et->Fill(Et);
	      hL1Jet_Eta->Fill(Eta);
	      hL1Jet_Phi->Fill(Phi);
	      hL1Jet_IET->Fill(IET);
	      hL1Jet_IEta->Fill(IEta);
	      hL1Jet_IPhi->Fill(IPhi);
	      hL1Jet_Bx->Fill(Bx);
	      hL1Jet_RawEt->Fill(RawEt);
	      hL1Jet_SeedEt->Fill(SeedEt);
	      hL1Jet_TowerIPhi->Fill(TowerIPhi);
	      hL1Jet_TowerIEta->Fill(TowerIEta);
	      hL1Jet_PUEt->Fill(PUEt);
	      hL1Jet_PUDonutEt0->Fill(PUDonutEt0);
	      hL1Jet_PUDonutEt1->Fill(PUDonutEt1);
	      hL1Jet_PUDonutEt2->Fill(PUDonutEt2);
	      hL1Jet_PUDonutEt3->Fill(PUDonutEt3);
	    }
	}

      // ======================================================================
      // L1EG Collection
      // ======================================================================
      if (cfg_AddEGs)
	{
	  if (cfg_DEBUG) cout << "\n=== L1EG (" << L1EG_Et->size() << ")" << endl;
	  vector<L1EG> L1EGs = GetL1EGs(cfg_DEBUG);
	  
	  hL1EG_nEGs->Fill(L1EG_Et->size());
	  // For-loop: L1 Jets
	  for (auto eg = L1EGs.begin(); eg != L1EGs.end(); eg++)
	    {
	      int index      = eg->getIndex();
	      double Et      = eg->getEt();
	      double Eta     = eg->getEta();
	      double Phi     = eg->getPhi();
	      int IEt        = eg->getIET();
	      int IEta       = eg->getIEta();
	      int IPhi       = eg->getIPhi();
	      int Iso        = eg->getIso();
	      int Bx         = eg->getBx();
	      int TowerIPhi  = 0; // eg->getTowerIPhi();
	      int TowerIEta  = 0; // eg->getTowerIEta();
	      int RawEt      = eg->getRawEt();
	      int IsoEt      = eg->getIsoEt();
	      int FootprintEt= eg->getFootprintEt();
	      int NTT        = eg->getNTT();
	      int Shape      = eg->getShape();
	      int TowerHoE   = eg->getTowerHoE();
	      
	      // Fill Histograms
	      hL1EG_Index->Fill(index);
	      hL1EG_Et->Fill(Et);
	      hL1EG_Eta->Fill(Eta);
	      hL1EG_Phi->Fill(Phi);
	      hL1EG_IEt->Fill(IEt);
	      hL1EG_IEta->Fill(IEta);
	      hL1EG_IPhi->Fill(IPhi);
	      hL1EG_Iso->Fill(Iso);
	      hL1EG_Bx->Fill(Bx);
	      hL1EG_TowerIPhi->Fill(TowerIPhi);
	      hL1EG_TowerIEta->Fill(TowerIEta);
	      hL1EG_RawEt->Fill(RawEt);
	      hL1EG_IsoEt->Fill(IsoEt);
	      hL1EG_FootprintEt->Fill(FootprintEt);
	      hL1EG_NTT->Fill(NTT);
	      hL1EG_Shape->Fill(Shape);
	      hL1EG_TowerHoE->Fill(TowerHoE);
	  }
	}

      // ======================================================================
      // L1Sum Collection
      // ======================================================================
      if (cfg_AddSums)
	{
	  if (cfg_DEBUG) cout << "\n=== L1Sum (" << L1Sum_Et->size() << ")" << endl;
	  vector<L1Sum> L1Sums = GetL1Sums(cfg_DEBUG);
	  
	  hL1Sum_nSums->Fill(L1Sum_Et->size());
	  // For-loop: L1 Jets
	  for (auto sum = L1Sums.begin(); sum != L1Sums.end(); sum++)
	    {
	      unsigned short theIndex = sum->getIndex();
	      float theEt       = sum->getEt();
	      float thePhi      = sum->getPhi();
	      short int theIEt  = sum->getIEt();
	      short int theIPhi = sum->getIPhi();
	      short int theType = sum->getType();
	      short int theBx   = sum->getBx();
	      
	      // Fill Histograms
	      hL1Sum_Index->Fill(theIndex);
	      hL1Sum_Et   ->Fill(theEt);
	      hL1Sum_Phi  ->Fill(thePhi);
	      hL1Sum_IEt  ->Fill(theIEt);
	      hL1Sum_IPhi ->Fill(theIPhi);
	      hL1Sum_Bx   ->Fill(theBx);
	      hL1Sum_Type ->Fill(theType);
	    }	      
	}
	  
	  
      // ======================================================================
      // L1Muons Collection
      // ======================================================================
      if (cfg_AddMuons)
	{
	  
	}


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
  if (cfg_DEBUG) std::cout << "=== Sanity::BookHistos_()" << std::endl;
  
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

  // L1Taus 
  histoTools_.BookHisto_1D(hL1Tau_Index    , "L1Tau_Index"    , "; index of L1 Tau; Entries"    ,  20, +0.0,  20.0);
  histoTools_.BookHisto_1D(hL1Tau_nTaus    , "L1Tau_nTaus"    , "; L1Tau multiplicity; Entries" ,  20, +0.0,  20.0);
  histoTools_.BookHisto_1D(hL1Tau_Et       , "L1Tau_Et"       , "; E_{T} (GeV/; Entries"        ,  60, +0.0, 300.0);
  histoTools_.BookHisto_1D(hL1Tau_Eta      , "L1Tau_Eta"      , "; #eta; Entries"               , 130, -2.6,   2.6);
  histoTools_.BookHisto_1D(hL1Tau_Phi      , "L1Tau_Phi"      , "; #phi (rads); Entries"        , 128, -3.2,   3.2);
  histoTools_.BookHisto_1D(hL1Tau_IET      , "L1Tau_IET"      , "; I-E_{T} (GeV/c); Entries"    , 200,  0.0, 200.0);
  histoTools_.BookHisto_1D(hL1Tau_IEta     , "L1Tau_IEta"     , "; I-#eta; Entries"             ,  66,-33.0, +33.0);
  histoTools_.BookHisto_1D(hL1Tau_IPhi     , "L1Tau_IPhi"     , "; I-#phi (rads); Entries"      , 145,  0.0, 145.0);
  histoTools_.BookHisto_1D(hL1Tau_Iso      , "L1Tau_Iso"      , "; Iso; Entries"                , 100,  0.0, 100.0);
  histoTools_.BookHisto_1D(hL1Tau_Bx       , "L1Tau_Bx"       , "; BX; Entries"                 , 200,  0.0, 200.0);
  histoTools_.BookHisto_1D(hL1Tau_TowerIPhi, "L1Tau_TowerIPhi", "; Tower I-#phi (rads); Entries", 128, -3.2,   3.2); 
  histoTools_.BookHisto_1D(hL1Tau_TowerIEta, "L1Tau_TowerIEta", "; Tower I-#eta; Entries"       , 130, -2.6,   2.6);
  histoTools_.BookHisto_1D(hL1Tau_RawEt    , "L1Tau_RawEt"    , "; Raw E_{T} (GeV); Entries"    ,  60,  0.0, 300.0);
  histoTools_.BookHisto_1D(hL1Tau_IsoEt    , "L1Tau_IsoEt"    , "; Iso E_{T} (GeV); Entries"    ,  60,  0.0, 300.0);
  histoTools_.BookHisto_1D(hL1Tau_NTT      , "L1Tau_NTT"      , "; NTT; Entries"                , 200,  0.0, 200.0);
  histoTools_.BookHisto_1D(hL1Tau_HasEM    , "L1Tau_HasEM"    , "; has EM; Entries"             ,   2, -0.5,   1.5);
  histoTools_.BookHisto_1D(hL1Tau_IsMerged , "L1Tau_IsMerged" , "; is Merged; Entries"          ,   2, -0.5,   1.5);
  histoTools_.BookHisto_1D(hL1Tau_HwQual   , "L1Tau_HwQual"   , "; Hw Qual; Entries"            ,   2, -0.5,   1.5);

  // L1Jets
  histoTools_.BookHisto_1D(hL1Jet_Index     , "L1Jet_Index"     , "; index of L1 Tau; Entries"       ,  20,  +0.0,  20.0);
  histoTools_.BookHisto_1D(hL1Jet_nJets     , "L1Jet_nJets"     , "; L1Jet multiplicity; Entries"    ,  20,  +0.0,  20.0);
  histoTools_.BookHisto_1D(hL1Jet_Et        , "L1Jet_Et"        , "; E_{T} (GeV); Entries"           ,  60,  +0.0, 300.0);
  histoTools_.BookHisto_1D(hL1Jet_Eta       , "L1Jet_Eta"       , "; #eta; Entries"                  , 300,  -6.0,   6.0);
  histoTools_.BookHisto_1D(hL1Jet_Phi       , "L1Jet_Phi"       , "; #phi (rads); Entries"           , 128,  -3.2,   3.2);
  histoTools_.BookHisto_1D(hL1Jet_IET       , "L1Jet_IET"       , "; I-E_{T} (GeV/c); Entries"       , 200,   0.0, 200.0);
  histoTools_.BookHisto_1D(hL1Jet_IEta      , "L1Jet_IEta"      , "; I-#eta; Entries"                , 250,-125.0, 125.0);
  histoTools_.BookHisto_1D(hL1Jet_IPhi      , "L1Jet_IPhi"      , "; I-#phi (rads); Entries"         , 145,   0.0, 145.0);
  histoTools_.BookHisto_1D(hL1Jet_Bx        , "L1Jet_Bx"        , "; BX; Entries"                    , 200,   0.0, 200.0);
  histoTools_.BookHisto_1D(hL1Jet_RawEt     , "L1Jet_RawEt"     , "; Raw E_{T} (GeV); Entries"       ,  60,   0.0, 300.0);
  histoTools_.BookHisto_1D(hL1Jet_SeedEt    , "L1Jet_SeedEt"    , "; Seed E_{T} (GeV); Entries"      ,  60,   0.0, 300.0);
  histoTools_.BookHisto_1D(hL1Jet_TowerIPhi , "L1Jet_TowerIPhi" , "; Tower I-#phi (rads); Entries"   , 128,  -3.2,   3.2); 
  histoTools_.BookHisto_1D(hL1Jet_TowerIEta , "L1Jet_TowerIEta" , "; Tower I-#eta; Entries"          , 130,  -2.6,   2.6);
  histoTools_.BookHisto_1D(hL1Jet_PUEt      , "L1Jet_PUEt"      , "; PU E_{T} (GeV); Entries"        ,  60,   0.0, 300.0);
  histoTools_.BookHisto_1D(hL1Jet_PUDonutEt0, "L1Jet_PUDonutEt0", "; PU Donut E_{T,0} (GeV); Entries",  50,   0.0,  50.0);
  histoTools_.BookHisto_1D(hL1Jet_PUDonutEt1, "L1Jet_PUDonutEt1", "; PU Donut E_{T,1} (GeV); Entries", 200,   0.0, 200.0);
  histoTools_.BookHisto_1D(hL1Jet_PUDonutEt2, "L1Jet_PUDonutEt2", "; PU Donut E_{T,2} (GeV); Entries", 250,   0.0, 250.0);
  histoTools_.BookHisto_1D(hL1Jet_PUDonutEt3, "L1Jet_PUDonutEt3", "; PU Donut E_{T,3} (GeV); Entries", 500,   0.0, 500.0);

  // L1EGs
  histoTools_.BookHisto_1D(hL1EG_Index      , "L1EG_Index"       , "; index of L1 EG; Entries"       ,  20,  +0.0,   20.0);
  histoTools_.BookHisto_1D(hL1EG_nEGs       , "L1EG_nEGs"        , "; L1EG multiplicity; Entries"    ,  20,  +0.0,   20.0);
  histoTools_.BookHisto_1D(hL1EG_Et         , "L1EG_Et"          , "; E_{T} (GeV); Entries"          ,  60,  +0.0,  300.0);
  histoTools_.BookHisto_1D(hL1EG_Eta        , "L1EG_Eta"         , "; #eta; Entries"                 , 300,  -6.0,    6.0);
  histoTools_.BookHisto_1D(hL1EG_Phi        , "L1EG_Phi"         , "; #phi (rads); Entries"          , 128,  -3.2,    3.2);
  histoTools_.BookHisto_1D(hL1EG_IEt        , "L1EG_IEt"         , "; I E_{T} (GeV); Entries"        , 200,   0.0,  200.0);
  histoTools_.BookHisto_1D(hL1EG_IEta       , "L1EG_IEta"        , "; I-#eta; Entries"               ,  66, -33.0,   33.0);
  histoTools_.BookHisto_1D(hL1EG_IPhi       , "L1EG_IPhi"        , "; I-#phi (rads); Entries"        , 145,   0.0,  145.0);
  histoTools_.BookHisto_1D(hL1EG_Iso        , "L1EG_Iso"         , "; Iso; Entries"                  ,   2,  -0.5,    1.5);
  histoTools_.BookHisto_1D(hL1EG_Bx         , "L1EG_Bx"          , "; Bx; Entries"                   , 200,   0.0,  200.0);
  histoTools_.BookHisto_1D(hL1EG_TowerIPhi  , "L1EG_TowerIPhi"   , "; Tower I-#phi (rads); Entries"  , 200,   0.0,  200.0);
  histoTools_.BookHisto_1D(hL1EG_TowerIEta  , "L1EG_TowerIEta"   , "; Tower I-#eta; Entries"         , 200,   0.0,  200.0);
  histoTools_.BookHisto_1D(hL1EG_RawEt      , "L1EG_RawEt"       , "; Raw E_{T} (GeV); Entries"      , 400,   0.0, 2000.0);
  histoTools_.BookHisto_1D(hL1EG_IsoEt      , "L1EG_IsoEt"       , "; Iso E_{T} (GeV); Entries"      , 400,   0.0, 2000.0);
  histoTools_.BookHisto_1D(hL1EG_FootprintEt, "L1EG_FootprintEt" , "; Footprint E_{T} (GeV); Entries", 400,   0.0, 2000.0);
  histoTools_.BookHisto_1D(hL1EG_NTT        , "L1EG_NTT"         , "; NTT; Entries"                  , 200,   0.0,  200.0);
  histoTools_.BookHisto_1D(hL1EG_Shape      , "L1EG_Shape"       , "; Shape; Entries"                , 200,   0.0,  200.0);
  histoTools_.BookHisto_1D(hL1EG_TowerHoE   , "L1EG_TowerHoE"    , "; Tower H/E; Entries"            ,  30, -10.0,   20.0);

  // L1Sums
  histoTools_.BookHisto_1D(hL1Sum_Index, "L1Sum_Index", "; index of L1 EG; Entries"       ,  20,  +0.0,   20.0);
  histoTools_.BookHisto_1D(hL1Sum_nSums, "L1Sum_nSums", "; L1EG multiplicity; Entries"    ,  20,  +0.0,   20.0);
  histoTools_.BookHisto_1D(hL1Sum_Et   , "L1Sum_Et"   , "; E_{T} (GeV); Entries"          , 600,  +0.0, 3000.0);
  histoTools_.BookHisto_1D(hL1Sum_Phi  , "L1Sum_Phi"  , "; #phi (rads); Entries"          , 128,  -3.2,    3.2);
  histoTools_.BookHisto_1D(hL1Sum_IEt  , "L1Sum_IEt"  , "; I E_{T} (GeV); Entries"        , 600,  +0.0, 3000.0);
  histoTools_.BookHisto_1D(hL1Sum_IPhi , "L1Sum_IPhi" , "; I-#phi (rads); Entries"        , 145,   0.0,  145.0);
  histoTools_.BookHisto_1D(hL1Sum_Bx   , "L1Sum_Bx"   , "; Bx; Entries"                   , 200,   0.0,  200.0);
  histoTools_.BookHisto_1D(hL1Sum_Type , "L1Sum_Type" , "; Type; Entries"                 ,  30,   0.0,   30.0);

  return;
}

#endif

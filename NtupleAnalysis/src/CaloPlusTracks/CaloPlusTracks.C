#ifndef CaloPlusTracks_cxx
#define CaloPlusTracks_cxx

// User
#include "../Auxiliary/interface/constants.h"
#include "CaloPlusTracks.h"

// ROOT
#include "TFitResult.h"
#include "TF1.h"
#include "Math/VectorUtil.h"

//****************************************************************************
void CaloPlusTracks::InitObjects(void)
//****************************************************************************
{

  pvProducer = new L1TkPrimaryVertex(this->s);

  return; 
}


//****************************************************************************
void CaloPlusTracks::InitVars_()
//****************************************************************************
{
  
  DEBUG = false;

  // Dataset-related
  datasets_  = datasets_.GetDataset(mcSample);
  realTauMom = datasets_.McTauMomPdgId_;
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;

  // Tracks
  selTks_Collection = "TTTracks"; // "TTTracks"; // "TTPixelTracks";
  selTks_nFitParams =   5;   // TP:   5
  selTks_minPt      =   2.0; // TP:   2.0
  selTks_maxEta     = 999.9; // TP: 999.9
  selTks_maxChiSq   = 200.0; // TP: 200.0
  selTks_minStubs   =   0;   // TP:   0
  selTks_minStubsPS =   0;   // TP:   0
  selTksPix_minHits =   0;   // TP:   N/A
  
  // L1TkTau - Signal cone
  matchTk_minPt_          = +5.00; // TP: 5.00
  matchTk_caloDeltaR_     = +0.12; // TP: 0.10
  sigCone_Constant        = +0.00; // TP: 0.00
  sigCone_dRMin           = +0.00; // WARNING! If > 0 the matching Track will NOT be added in sigCone_TTTracks.
  sigCone_dRMax           = +0.15; // TP: 0.15
  sigCone_cutoffDeltaR    = +0.15; // TP: 0.15
  sigCone_maxTkInvMass    = +1.77; // TP: Unused (3-pr)
  sigCone_maxTkDeltaPOCAz = +0.20; // TP: Unused (3-pr)

  // Isolation cone
  isoCone_Constant = +3.50;         // TP: 3.50 GeV
  isoCone_VtxIsoWP = +0.40;         // TP: 1.0cm
  isoCone_dRMin    = sigCone_dRMax; // TP: 0.4cm
  isoCone_dRMax    = +0.30;         // TP: 0.4cm
  diTau_deltaPOCAz = +0.50;         // TP: 1.0cm

  // MC matching
  mcMatching_dRMax  = +0.05;
  mcMatching_unique = true;

  PrintSettings();

  return;
}


//****************************************************************************
void CaloPlusTracks::PrintSettings(void)
//****************************************************************************
{

  // Inform user of settings
  Table settings("Variable | Cut | Value | TP 2015 | Units", "Text");  // Table settingsTable("Variable & Value & Units", "LaTeX", "l l l");
  settings.AddRowColumn(0, "MC Sample");
  settings.AddRowColumn(0, "==" );
  settings.AddRowColumn(0, mcSample );

  settings.AddRowColumn(1, "");
  
  settings.AddRowColumn(2, "Tracks: Collection");
  settings.AddRowColumn(2, "==");
  settings.AddRowColumn(2, selTks_Collection);
  settings.AddRowColumn(2, "TTTracks");
  settings.AddRowColumn(2, "");
  
  settings.AddRowColumn(3, "Tracks: Fit Parameters");
  settings.AddRowColumn(3, "==");
  settings.AddRowColumn(3, auxTools_.ToString( selTks_nFitParams) );
  settings.AddRowColumn(3, "5");
  settings.AddRowColumn(3, "");

  settings.AddRowColumn(4, "Tracks: Pt");
  settings.AddRowColumn(4, ">=");
  settings.AddRowColumn(4, auxTools_.ToString( selTks_minPt) );
  settings.AddRowColumn(4, "2" );
  settings.AddRowColumn(4, "GeV/c" );
  
  settings.AddRowColumn(5, "Tracks: |Eta|");
  settings.AddRowColumn(5, "<=");
  settings.AddRowColumn(5, auxTools_.ToString( selTks_maxEta) );
  settings.AddRowColumn(5, "1e+03" );
  settings.AddRowColumn(5, "" );
  
  settings.AddRowColumn(6, "Tracks: ChiSq");
  settings.AddRowColumn(6, "<=");
  settings.AddRowColumn(6, auxTools_.ToString( selTks_maxChiSq) );
  settings.AddRowColumn(6, "10 (8) (chi2/DOF)");
  settings.AddRowColumn(6, "");

  settings.AddRowColumn(7, "Tracks: Stubs");
  settings.AddRowColumn(7, ">=");
  settings.AddRowColumn(7, auxTools_.ToString( selTks_minStubs) );
  settings.AddRowColumn(7, "4" );
  settings.AddRowColumn(7, "" );
  
  settings.AddRowColumn(8, "Tracks: PS-Stubs");
  settings.AddRowColumn(8, ">=");
  settings.AddRowColumn(8, auxTools_.ToString( selTks_minStubsPS) );
  settings.AddRowColumn(8, "0" );
  settings.AddRowColumn(8, "" );

  settings.AddRowColumn(9, "Tracks: Pixel Hits");
  settings.AddRowColumn(9, ">=");
  settings.AddRowColumn(9, auxTools_.ToString( selTksPix_minHits) );
  settings.AddRowColumn(9, "N/A" );
  settings.AddRowColumn(9, "" );
  			
  settings.AddRowColumn(10, "");

  settings.AddRowColumn(11, "Matching Track: Pt");
  settings.AddRowColumn(11, ">=");
  settings.AddRowColumn(11, auxTools_.ToString(matchTk_minPt_) );
  settings.AddRowColumn(11, "15" );
  settings.AddRowColumn(11, "GeV/c");

  settings.AddRowColumn(12, "Matching Track: DeltaR");
  settings.AddRowColumn(12, "<=");
  settings.AddRowColumn(12, auxTools_.ToString(matchTk_caloDeltaR_) );
  settings.AddRowColumn(12, "0.10" );
  settings.AddRowColumn(12, "");

  settings.AddRowColumn(13, "");

  settings.AddRowColumn(14, "Signal Cone: Shrink Constant");
  settings.AddRowColumn(14, "==");
  settings.AddRowColumn(14, auxTools_.ToString(sigCone_Constant) );
  settings.AddRowColumn(14, "0" );
  settings.AddRowColumn(14, "GeV");

  settings.AddRowColumn(15, "Signal Cone: DeltaR");
  settings.AddRowColumn(15, ">=");
  settings.AddRowColumn(15, auxTools_.ToString(sigCone_dRMin) );
  settings.AddRowColumn(15, "0.0" );
  settings.AddRowColumn(15, "" );

  settings.AddRowColumn(16, "Signal Cone: DeltaR");
  settings.AddRowColumn(16, "<=");
  settings.AddRowColumn(16, auxTools_.ToString(sigCone_dRMax) );
  settings.AddRowColumn(16, "0.15" );
  settings.AddRowColumn(16, "" );
  
  settings.AddRowColumn(17, "Signal Cone:-3pr InvMass");
  settings.AddRowColumn(17, "<=");
  settings.AddRowColumn(17, auxTools_.ToString(sigCone_maxTkInvMass) );
  settings.AddRowColumn(17, "N/A" );
  settings.AddRowColumn(17, "GeV/c^{-2}");

  settings.AddRowColumn(18, "Signal Cone:-3pr maxTkDeltaPOCAz");
  settings.AddRowColumn(18, "<=");
  settings.AddRowColumn(18, auxTools_.ToString(sigCone_maxTkDeltaPOCAz) );
  settings.AddRowColumn(18, "N/A" );
  settings.AddRowColumn(18, "cm");

  settings.AddRowColumn(19, "");

  settings.AddRowColumn(20, "Isolation Cone: Shrink Constant");
  settings.AddRowColumn(20, "==");
  settings.AddRowColumn(20, auxTools_.ToString(isoCone_Constant) );
  settings.AddRowColumn(20, "3.5");
  settings.AddRowColumn(20, "GeV");

  settings.AddRowColumn(21, "Isolation Cone: DeltaR");
  settings.AddRowColumn(21, ">=");
  settings.AddRowColumn(21, auxTools_.ToString(isoCone_dRMin) );
  settings.AddRowColumn(21, "0.15" );
  settings.AddRowColumn(21, "" );

  settings.AddRowColumn(22, "Isolation Cone: DeltaR");
  settings.AddRowColumn(22, "=<");
  settings.AddRowColumn(22, auxTools_.ToString(isoCone_dRMax) );
  settings.AddRowColumn(22, "0.30");
  settings.AddRowColumn(22, "");

  settings.AddRowColumn(23, "Isolation Cone: VtxIso" );
  settings.AddRowColumn(23, "<=" );
  settings.AddRowColumn(23, auxTools_.ToString(isoCone_VtxIsoWP) );
  settings.AddRowColumn(23, "1.0");
  settings.AddRowColumn(23, "cm");
  settings.AddRowColumn(23, "");

  settings.AddRowColumn(24, "Di-Tau |Delta z0|");
  settings.AddRowColumn(24, "<");
  settings.AddRowColumn(24, auxTools_.ToString(diTau_deltaPOCAz) );
  settings.AddRowColumn(24, "1.0" );
  settings.AddRowColumn(24, "cm");

  settings.AddRowColumn(25, "");

  settings.AddRowColumn(26, "MC-Matching DeltaR");
  settings.AddRowColumn(26, "<=");
  settings.AddRowColumn(26, auxTools_.ToString(mcMatching_dRMax) );
  settings.AddRowColumn(26, "0.05" );
  settings.AddRowColumn(26, "" );

  settings.AddRowColumn(27, "MC-Matching IsUnique");
  settings.AddRowColumn(27, "==");
  settings.AddRowColumn(27, auxTools_.ToString(mcMatching_unique) );
  settings.AddRowColumn(27, "1" );
  settings.AddRowColumn(27, "" );
  
  settings.AddRowColumn(28, "MC-Taus: Mom PdgId");
  settings.AddRowColumn(28, "==");
  settings.AddRowColumn(28, auxTools_.ToString(realTauMom));
  settings.AddRowColumn(28, "N/A" );
  settings.AddRowColumn(28, "" );

  settings.AddRowColumn(29, "MC-Taus: Number Expected");
  settings.AddRowColumn(29, ">=");
  settings.AddRowColumn(29, auxTools_.ToString(nMaxNumOfHTausPossible));
  settings.AddRowColumn(29, "N/A" );
  settings.AddRowColumn(29, "" );

  settings.Print();
  
  return;
}


//****************************************************************************
void CaloPlusTracks::Loop()
//****************************************************************************
{
  
  InitVars_(); 
  if (fChain == 0) return;  
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  cout << "Analyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;

  // Initialisations
  BookHistos_();
  Long64_t nbytes       = 0;
  Long64_t nb           = 0;
  int nEvtsWithMaxHTaus = 0; 
  L1PixelTrackFit f(3.8112); // Bz in Tesla

 
  ////////////////////////////////////////////////
  // For-loop: Entries
  ////////////////////////////////////////////////
  for (int jentry = 0; jentry < nEntries; jentry++){
  
    if(DEBUG) cout << "Entry = " << jentry << endl;
    
    // Init variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    
    // Gen-Collections
    vector<GenParticle> GenParticles    = GetGenParticles();
    vector<GenParticle> GenTaus         = GetGenParticles(15, true);
    vector<GenParticle> GenTausHadronic = GetHadronicGenTaus(GenTaus, 00.0, 999.9);
    vector<GenParticle> GenTausTrigger  = GetHadronicGenTaus(GenTaus, 20.0, 2.3);
    // Print Collections?
    if (0) for (vector<GenParticle>::iterator p = GenParticles.begin(); p != GenParticles.end(); p++) p->PrintProperties();
    if (0) for (vector<GenParticle>::iterator p = GenTaus.begin(); p != GenTaus.end(); p++) p->PrintProperties();
    if (0) for (vector<GenParticle>::iterator p = GenTausHadronic.begin(); p != GenTausHadronic.end(); p++) p->PrintProperties();
    if (0) for (vector<GenParticle>::iterator p = GenTausTrigger.begin(); p != GenTausTrigger.end(); p++) p->PrintProperties();

    
    // Track Collections
    vector<TrackingParticle> TPs       = GetTrackingParticles();
    vector<TTTrack> matchTTTracks      = GetTTTracks(matchTk_minPt_, selTks_maxEta, selTks_maxChiSq, selTks_minStubs, selTks_minStubsPS, selTks_nFitParams);
    vector<TTTrack> selTTTracks        = GetTTTracks(selTks_minPt  , selTks_maxEta, selTks_maxChiSq, selTks_minStubs, selTks_minStubsPS, selTks_nFitParams);
    vector<TTPixelTrack> TTPixelTracks = GetTTPixelTracks(selTks_minPt, selTks_maxEta, selTks_maxChiSq, selTksPix_minHits);
    vector<TTTrack> pvTTTracks;
    double pv_z = GetPVTTTracks(pvTTTracks);
    // Print Collections?
    if (0) for (vector<TrackingParticle>::iterator p = TPs.begin(); p != TPs.end(); p++) p->PrintProperties();
    if (0) for (vector<TTTrack>::iterator t = selTTTracks.begin(); t != selTTTracks.end(); t++) t->PrintProperties();
    if (0) for (vector<TTTrack>::iterator t = matchTTTracks.begin(); t != matchTTTracks.end(); t++) t->PrintProperties();
    if (0) for (vector<TTTrack>::iterator t = pvTTTracks.begin(); t != pvTTTracks.end(); t++) t->PrintProperties();
    if (0) for (vector<TTPixelTrack>::iterator t = TTPixelTracks.begin(); t != TTPixelTracks.end(); t++) t->PrintProperties();

    
    // Tau Collections
    vector<L1JetParticle> L1CaloTaus = GetL1CaloTaus();
    vector<L1TkTauParticle> L1TkTauCandidates;
    vector<L1TkTauParticle> L1TkTaus_Calo;
    vector<L1TkTauParticle> L1TkTaus_Tk;
    vector<L1TkTauParticle> L1TkTaus_VtxIso;
    // Print Collections?
    if (0) for (vector<L1JetParticle>::iterator j = L1CaloTaus.begin(); j != L1CaloTaus.end(); j++) j->PrintProperties();


    // Sanity checks
    auxTools_.EnsureVectorIsSorted(*L1CaloTau_Et, true);
    // auxTools_.EnsureVectorIsSorted(*TP_Pt, true); // is not sorted. does not have to be either
    auxTools_.EnsureVectorIsSorted(*L1Tks_Pt    , true);
    auxTools_.EnsureVectorIsSorted(*L1PixTks_Pt , true);


    // Ensure that all taus are found
    bool bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;
    
    
    ////////////////////////////////////////////////
    // For-loop: L1CaloTaus
    ////////////////////////////////////////////////
    for (vector<L1JetParticle>::iterator calo = L1CaloTaus.begin(); calo != L1CaloTaus.end(); calo++)
      {

	// Calculate the Et-dependent signal & isolation cone sizes
	GetShrinkingConeSizes(calo->et(), sigCone_Constant, isoCone_Constant, sigCone_cutoffDeltaR,
			      sigCone_dRMin, sigCone_dRMax, isoCone_dRMin, isoCone_dRMax);

	// Construct the CaloPlusTracks candidate
	L1TkTauParticle	L1TkTauCandidate(0.0, matchTk_caloDeltaR_, sigCone_dRMin, sigCone_dRMax, isoCone_dRMin, isoCone_dRMax);

	//
	if ( selTks_Collection.compare("TTPixelTracks") == 0 )
	  {
	    // FIXME 
	  }      
	else if ( selTks_Collection.compare("TTTracks") == 0 )
	  {
	    
	    GetMatchingTrack(L1TkTauCandidate, *calo, matchTTTracks);
	    GetSigConeTracks(L1TkTauCandidate, selTTTracks);
	    GetIsoConeTracks(L1TkTauCandidate, selTTTracks);
	    GetIsolationValues(L1TkTauCandidate);
	    GetMatchingGenParticle(L1TkTauCandidate, GenTaus);
	    if (0) L1TkTauCandidate.PrintProperties(false, false, true, true);
	  }
	else
	  {
	    cout << "E R R O R ! CaloPlusTracks::Loop(...) - Invalid Track Collection Type \"" << selTks_Collection << "\". EXIT" << endl;
	    exit(1);
	  }
	
	// Save L1TkTau Candidate
	L1TkTauCandidates.push_back(L1TkTauCandidate);
      }
    
    
    ////////////////////////////////////////////////
    /// Create L1TkTaus Collections
    ////////////////////////////////////////////////
    for (size_t i = 0; i < L1TkTauCandidates.size(); i++)
      {
	
	L1TkTauParticle L1TkTau = L1TkTauCandidates.at(i);

	// Calo
	if (!L1TkTau.HasMatchingTk() ) L1TkTaus_Calo.push_back(L1TkTau);
	else
	  {
	    // +Tk
	    L1TkTaus_Tk.push_back(L1TkTau);
	    
	    // +VtxIso
	    if ( L1TkTau.GetVtxIsolation() > isoCone_VtxIsoWP) L1TkTaus_VtxIso.push_back(L1TkTau);
	  }
	
      }// L1TkTauCandidates

    // Print the L1TkTauCollections
    if (0) for (vector<L1TkTauParticle>::iterator j = L1TkTauCandidates.begin(); j != L1TkTauCandidates.end(); j++) j->PrintProperties();
    if (0) for (vector<L1TkTauParticle>::iterator j = L1TkTaus_Calo.begin(); j != L1TkTaus_Calo.end(); j++) j->PrintProperties();
    if (0) for (vector<L1TkTauParticle>::iterator j = L1TkTaus_Tk.begin(); j != L1TkTaus_Tk.end(); j++) j->PrintProperties();
    if (0) for (vector<L1TkTauParticle>::iterator j = L1TkTaus_VtxIso.begin(); j != L1TkTaus_VtxIso.end(); j++) j->PrintProperties();

    
    ////////////////////////////////////////////////
    // Event-Type Histograms
    ////////////////////////////////////////////////
    hHepMCEvt_VtxX_VtxY     ->Fill(HepMCEvt_VtxX, HepMCEvt_VtxY);
    hHepMCEvt_VtxZ          ->Fill(HepMCEvt_VtxZ);
    hL1TkPV_VtxZ            ->Fill(pv_z);
    hPrimaryVertex_DeltaZ   ->Fill(pv_z - HepMCEvt_VtxZ);
    hPrimaryVertex_AbsDeltaZ->Fill( abs(pv_z - HepMCEvt_VtxZ) );

    
    ////////////////////////////////////////////////
    /// L1Tktau Properties 
    ////////////////////////////////////////////////
    hL1TkTau_Multiplicity ->Fill( L1TkTaus_VtxIso.size() );
    
    // For-loop: L1TkTaus_VtxIso
    for (vector<L1TkTauParticle>::iterator tau = L1TkTaus_VtxIso.begin(); tau != L1TkTaus_VtxIso.end(); tau++)
      {
	
	if (0) tau->PrintProperties(true, true, true, true);
	
	// Variables
	TLorentzVector sigTks_p4 = tau->GetSigConeTTTracksP4();
	TLorentzVector isoTks_p4 = tau->GetIsoConeTTTracksP4();
	// Do not Skip if using MinBias sample as no real taus exist!
	if (!tau->HasMatchingGenParticle() && ( mcSample.compare("MinBias") !=0 ) ) continue;
	
	// Fill Histos
	hL1TkTau_Rtau        ->Fill( tau->GetSigConeLdgTk().getPt() / tau->GetCaloTau().et() );
	hL1TkTau_CHF         ->Fill( tau->GetCaloTau().et()/sigTks_p4.Et() );
	hL1TkTau_NHF         ->Fill( (tau->GetCaloTau().et() - sigTks_p4.Et())/tau->GetCaloTau().et() );
	hL1TkTau_NHFAbs      ->Fill( std::abs( (tau->GetCaloTau().et() - sigTks_p4.Et())/tau->GetCaloTau().et() ) );
	hL1TkTau_NSigTks     ->Fill( tau->GetSigConeTTTracks().size() );
	hL1TkTau_NIsoTks     ->Fill( tau->GetIsoConeTTTracks().size() );
	hL1TkTau_InvMass     ->Fill( sigTks_p4.M() ); 
	hL1TkTau_InvMassIncl ->Fill( sigTks_p4.M() + isoTks_p4.M() );
	hL1TkTau_SigConeRMin ->Fill( tau->GetSigConeMin() );
	hL1TkTau_IsoConeRMin ->Fill( tau->GetIsoConeMin() );
	hL1TkTau_SigConeRMax ->Fill( tau->GetSigConeMax() );
	hL1TkTau_IsoConeRMax ->Fill( tau->GetIsoConeMax() );
	hL1TkTau_DeltaRGenP  ->Fill( tau->GetMatchingGenParticleDeltaR() );
	hL1TkTau_RelIso      ->Fill( tau->GetRelIsolation() );
	hL1TkTau_VtxIso      ->Fill( tau->GetVtxIsolation() );
	hL1TkTau_VtxIsoAbs   ->Fill( abs(tau->GetVtxIsolation()) );

	TTTrack matchTk   = tau->GetMatchingTk();
	double matchTk_dR = auxTools_.DeltaR(matchTk.getEta(), matchTk.getPhi(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi() );
	TLorentzVector caloTau_p4;
	caloTau_p4.SetPtEtaPhiE(tau->GetCaloTau().et(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi(), tau->GetCaloTau().energy() );
	hL1TkTau_MatchTk_DeltaR        ->Fill( matchTk_dR );
	hL1TkTau_MatchTk_PtRel         ->Fill( matchTk.p3().Perp(caloTau_p4.Vect()) );
	hL1TkTau_MatchTk_Pt            ->Fill( matchTk.getPt() );
	hL1TkTau_MatchTk_Eta           ->Fill( matchTk.getEta() );
	hL1TkTau_MatchTk_POCAz         ->Fill( matchTk.getZ0() );
	hL1TkTau_MatchTk_d0            ->Fill( matchTk.getD0() );
	hL1TkTau_MatchTk_d0Abs         ->Fill( std::abs(matchTk.getD0()) );
	hL1TkTau_MatchTk_NStubs        ->Fill( matchTk.getNumOfStubs() );
	hL1TkTau_MatchTk_NPsStubs      ->Fill( matchTk.getNumOfStubsPS() );
	hL1TkTau_MatchTk_NBarrelStubs  ->Fill( matchTk.getNumOfBarrelStubs() );
	hL1TkTau_MatchTk_NEndcapStubs  ->Fill( matchTk.getNumOfEndcapStubs() );
	hL1TkTau_MatchTk_StubPtCons    ->Fill( matchTk.getStubPtConsistency() );
	hL1TkTau_MatchTk_ChiSquared    ->Fill( matchTk.getChi2() );
	hL1TkTau_MatchTk_RedChiSquared ->Fill( matchTk.getChi2Red() );
	hL1TkTau_MatchTk_IsGenuine     ->Fill( matchTk.getIsGenuine() );
	hL1TkTau_MatchTk_IsUnknown     ->Fill( matchTk.getIsUnknown() );
	hL1TkTau_MatchTk_IsCombinatoric->Fill( matchTk.getIsCombinatoric() );
	

	// Matching TTTrack
	TTTrack match_tk = tau->GetMatchingTk();
	int sigTks_sumCharge = 0;
	  
	// For-loop: SigCone TTTracks
	vector<TTTrack> sigTks = tau->GetSigConeTTTracks();
	for (vector<TTTrack>::iterator sigTk = sigTks.begin(); sigTk != sigTks.end(); sigTk++)
	  {
	    if (0) sigTk->PrintProperties();

	    // Skip if tk is matching track
	    if ( sigTk->index() == match_tk.index() ) continue;
	    
	    // Get the transverse component of this track with respect to the matching track
	    TVector3 sigTk_p3  = sigTk->getMomentum();
	    double sigTk_PtRel = sigTk_p3.Perp( match_tk.getMomentum() );
	 
	    // Fill Histograms
	    hL1TkTau_SigTks_Pt        ->Fill( sigTk->getPt()  );
	    hL1TkTau_SigTks_Eta       ->Fill( sigTk->getEta() );
	    hL1TkTau_SigTks_POCAz     ->Fill( sigTk->getZ0()  );
	    hL1TkTau_SigTks_DeltaPOCAz->Fill( std::abs( sigTk->getZ0() - match_tk.getZ0() ) );
	    hL1TkTau_SigTks_d0        ->Fill( sigTk->getD0() );
	    hL1TkTau_SigTks_d0Abs     ->Fill( abs( sigTk->getD0()) );
	    // hL1TkTau_SigTks_d0Sig     ->Fill( sigTk->getD0()/sigTk->getSigmaD0() );             // TTPixelTracks
	    // hL1TkTau_SigTks_d0SigAbs  ->Fill( std::abs(sigTk->getD0()/sigTk->getSigmaD0() ) );  // TTPixelTracks
	    hL1TkTau_SigTks_d0Sig     ->Fill( -1.0 );
	    hL1TkTau_SigTks_d0SigAbs  ->Fill( -1.0 );
	    hL1TkTau_SigTks_PtRel     ->Fill( sigTk_PtRel );
	    hL1TkTau_SigTks_StubPtCons->Fill( sigTk->getStubPtConsistency() );

	    // Other variables
	    sigTks_sumCharge += sigTk->getCharge();
	    
	  }// SigCone_TTTracks

	// Fill histos for other variables
	hL1TkTau_Charge->Fill( sigTks_sumCharge);

	
	// For-loop: IsoCone TTTracks
	vector<TTTrack> isoTks = tau->GetIsoConeTTTracks();
	for (vector<TTTrack>::iterator isoTk = isoTks.begin(); isoTk != isoTks.end(); isoTk++)
	  {

	    if (0) isoTk->PrintProperties();
	    
	    // Skip if tk is matching track
	    if ( isoTk->index() == match_tk.index() ) continue;
	
	    // Get the transverse component of this track with respect to the matching track
	    TVector3 isoTk_p3  = isoTk->getMomentum();
	    double isoTk_PtRel = isoTk_p3.Perp( match_tk.getMomentum() );
	
	    // Fill Histograms
	    hL1TkTau_IsoTks_Pt        ->Fill( isoTk->getPt()  );
	    hL1TkTau_IsoTks_Eta       ->Fill( isoTk->getEta() );
	    hL1TkTau_IsoTks_POCAz     ->Fill( isoTk->getZ0()  );
	    hL1TkTau_IsoTks_DeltaPOCAz->Fill( std::abs( isoTk->getZ0() - match_tk.getZ0() ) );
	    hL1TkTau_IsoTks_d0        ->Fill( isoTk->getD0() );							     
	    hL1TkTau_IsoTks_d0Abs     ->Fill( abs( isoTk->getD0()) );
	    // hL1TkTau_IsoTks_d0Sig     ->Fill( isoTk->getD0()/isoTk->getSigmaD0() );             // TTPixelTracks
	    // hL1TkTau_IsoTks_d0SigAbs  ->Fill( std::abs(isoTk->getD0()/isoTk->getSigmaD0() ) );  // TTPixelTracks
	    hL1TkTau_IsoTks_d0Sig     ->Fill( -1.0 );
	    hL1TkTau_IsoTks_d0SigAbs  ->Fill( -1.0 );
	    hL1TkTau_IsoTks_PtRel     ->Fill( isoTk_PtRel );
	    hL1TkTau_IsoTks_StubPtCons->Fill( isoTk->getStubPtConsistency() );
	  }// IsoCone_TTTracks

      } // L1TkTaus_VtxIso
  

    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////
    FillTurnOn_Denominator_(GenTausHadronic, hMcHadronicTau_VisEt);
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 25.0, hCalo_TurnOn25  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 25.0, hTk_TurnOn25    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn25);
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 50.0, hCalo_TurnOn50  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 50.0, hTk_TurnOn50    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn50);
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 66.0, hCalo_TurnOn_SingleTau50KHz  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 65.0, hTk_TurnOn_SingleTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn_SingleTau50KHz);
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 42.0, hCalo_TurnOn_DiTau50KHz  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 40.0, hTk_TurnOn_DiTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn_DiTau50KHz);

    
    ////////////////////////////////////////////////
    // SingleTau
    ////////////////////////////////////////////////
    FillSingleTau_(L1TkTaus_Calo  , hCalo_Rate  , hCalo_Eff  );
    FillSingleTau_(L1TkTaus_Tk    , hTk_Rate    , hTk_Eff    );
    FillSingleTau_(L1TkTaus_VtxIso, hVtxIso_Rate, hVtxIso_Eff);


    ////////////////////////////////////////////////
    // DiTau
    ////////////////////////////////////////////////
    // FillDiTau_(L1TkTaus_Calo  , L1TkTaus_Tk       , hDiTau_Rate_Calo_Tk    , hDiTau_Eff_Calo_Tk    );
    // FillDiTau_(L1TkTaus_Calo  , L1TkTaus_VtxIso   , hDiTau_Rate_Calo_VtxIso, hDiTau_Eff_Calo_VtxIso);
    // FillDiTau_(L1TkTaus_Tk    , L1TkTaus_VtxIso   , hDiTau_Rate_Tk_VtxIso  , hDiTau_Eff_Tk_VtxIso  ); // FIXME
    

    ////////////////////////////////////////////////
    // WARNING: Removal of non Z-matching should be just before I need it!
    ////////////////////////////////////////////////
    ApplyDiTauZMatching(selTks_Collection, L1TkTaus_Tk);      // fixme - erases L1TkTaus from vector!
    ApplyDiTauZMatching(selTks_Collection, L1TkTaus_VtxIso);  // fixme - erases L1TkTaus from vector!
    FillDiTau_(L1TkTaus_Calo  , hDiTau_Rate_Calo  , hDiTau_Eff_Calo  );
    FillDiTau_(L1TkTaus_Tk    , hDiTau_Rate_Tk    , hDiTau_Eff_Tk    );
    FillDiTau_(L1TkTaus_VtxIso, hDiTau_Rate_VtxIso, hDiTau_Eff_VtxIso);
    
    ////////////////////////////////////////////////
    // Progress bar
    ////////////////////////////////////////////////
    auxTools_.ProgressBar(jentry, nEntries, 100, 100);
    
  }// For-loop: Entries


  ////////////////////////////////////////////////
  // Convert or Finalise Histos
  ////////////////////////////////////////////////
  // SingleTau
  histoTools_.ConvertToRateHisto_1D(hCalo_Rate  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hTk_Rate    , nEntries);
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate, nEntries);
  FinaliseEffHisto_( hCalo_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hTk_Eff    , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIso_Eff, nEvtsWithMaxHTaus);

  // DiTau
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Calo  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk    , nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso, nEntries);
  FinaliseEffHisto_( hDiTau_Eff_Calo  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Tk    , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIso, nEvtsWithMaxHTaus);

  // DiTau (Calo-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_Tk    , nEntries);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_VtxIso, nEntries);
  FinaliseEffHisto_( hDiTau_Eff_Calo_Tk    , nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Calo_VtxIso, nEvtsWithMaxHTaus);;

  // DiTau (Tk-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_VtxIso, nEntries);
  FinaliseEffHisto_( hDiTau_Eff_Tk_VtxIso, nEvtsWithMaxHTaus);;

  // Turn-Ons
  histoTools_.DivideHistos_1D(hCalo_TurnOn50   , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn50     , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50 , hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hCalo_TurnOn25   , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn25     , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25 , hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hCalo_TurnOn_SingleTau50KHz  , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn_SingleTau50KHz    , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn_SingleTau50KHz, hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hCalo_TurnOn_DiTau50KHz  , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn_DiTau50KHz    , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn_DiTau50KHz, hMcHadronicTau_VisEt);


  ////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////
  outFile->cd();
  outFile->Write();
  auxTools_.StopwatchStop(5, "minutes");

}


//****************************************************************************
void CaloPlusTracks::BookHistos_(void)
//****************************************************************************
{
  
  // Event-Type Histograms
  histoTools_.BookHisto_1D(hL1TkPV_VtxZ            , "L1TkPV_VtxZ"            ,  600, -30.0,  +30.0);
  histoTools_.BookHisto_1D(hHepMCEvt_VtxZ          , "HepMCEvt_VtxZ"          ,  600, -30.0,  +30.0);
  histoTools_.BookHisto_2D(hHepMCEvt_VtxX_VtxY     , "HepMCEvt_VtxX_VtxY"     ,  400,  -0.01,  +0.01, 400,  -0.01,  +0.01);
  histoTools_.BookHisto_1D(hPrimaryVertex_DeltaZ   , "PrimaryVertex_DeltaZ"   ,  600, -30.0,  +30.0);
  histoTools_.BookHisto_1D(hPrimaryVertex_AbsDeltaZ, "PrimaryVertex_AbsDeltaZ", 1200, + 0.0,  +60.0);

  // VtxIsolated L1TkTaus
  histoTools_.BookHisto_1D(hL1TkTau_Multiplicity   , "L1TkTau_Multiplicity",     30, -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_Rtau           , "L1TkTau_Rtau"        ,   1000,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_CHF            , "L1TkTau_CHF"         ,  10000,  0.0, +100.0);
  histoTools_.BookHisto_1D(hL1TkTau_NHF            , "L1TkTau_NHF"         ,   2000,-10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_NHFAbs         , "L1TkTau_NHFAbs"      ,   1000,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_NSigTks        , "L1TkTau_NSigTks"     ,     15, -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_NIsoTks        , "L1TkTau_NIsoTks"     ,     15, -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_InvMass        , "L1TkTau_InvMass"     ,   1000,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_InvMassIncl    , "L1TkTau_InvMassIncl" ,   1000,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_SigConeRMin    , "L1TkTau_SigConeRMin" ,   2000,  0.0,   +2.0);
  histoTools_.BookHisto_1D(hL1TkTau_SigConeRMax    , "L1TkTau_SigConeRMax" ,   2000,  0.0,   +2.0);
  histoTools_.BookHisto_1D(hL1TkTau_IsoConeRMin    , "L1TkTau_IsoConeRMin" ,   2000,  0.0,   +2.0);
  histoTools_.BookHisto_1D(hL1TkTau_IsoConeRMax    , "L1TkTau_IsoConeRMax" ,   2000,  0.0,   +2.0);
  histoTools_.BookHisto_1D(hL1TkTau_Charge         , "L1TkTau_Charge"      ,     19, -9.5,   +9.5);
  histoTools_.BookHisto_1D(hL1TkTau_RelIso         , "L1TkTau_RelIso"      ,   1000,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_VtxIso         , "L1TkTau_VtxIso"      ,   1200,-30.0,  +30.0);
  histoTools_.BookHisto_1D(hL1TkTau_VtxIsoAbs      , "L1TkTau_VtxIsoAbs"   ,    600,  0.0,  +30.0);
  histoTools_.BookHisto_1D(hL1TkTau_DeltaRGenP     , "L1TkTau_DeltaRGenP"  ,   2000,  0.0,   +2.0);
  //
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Pt               , "L1TkTau_SigTks_Pt"           ,  300,    +0.0,  +300.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Eta              , "L1TkTau_SigTks_Eta"          ,  600,    -3.0,    +3.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_POCAz            , "L1TkTau_SigTks_POCAz"        ,  600,   -30.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_DeltaPOCAz       , "L1TkTau_SigTks_DeltaPOCAz"   ,  600,    +0.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0               , "L1TkTau_SigTks_d0"           , 2000,   -10.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0Abs            , "L1TkTau_SigTks_d0Abs"        , 1000,     0.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0Sig            , "L1TkTau_SigTks_d0Sig"        , 4000,   -20.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0SigAbs         , "L1TkTau_SigTks_d0SigAbs"     , 2000,     0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_PtRel            , "L1TkTau_SigTks_PtRel"        , 2000,    +0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_StubPtCons       , "L1TkTau_SigTks_StubPtCons"   ,  200,    +0.0,  +200.0  );
  //
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Pt               , "L1TkTau_IsoTks_Pt"           ,  300,    +0.0,  +300.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Eta              , "L1TkTau_IsoTks_Eta"          ,  600,    -3.0,    +3.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_POCAz            , "L1TkTau_IsoTks_POCAz"        ,  600,   -30.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_DeltaPOCAz       , "L1TkTau_IsoTks_DeltaPOCAz"   ,  600,    +0.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0               , "L1TkTau_IsoTks_d0"           , 2000,   -10.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0Abs            , "L1TkTau_IsoTks_d0Abs"        , 1000,     0.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0Sig            , "L1TkTau_IsoTks_d0Sig"        , 4000,   -20.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0SigAbs         , "L1TkTau_IsoTks_d0SigAbs"     , 2000,     0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_PtRel            , "L1TkTau_IsoTks_PtRel"        , 2000,    +0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_StubPtCons       , "L1TkTau_IsoTks_StubPtCons"   ,  200,    +0.0,  +200.0  );

  // VtxIsolated L1TkTaus, Signal TTTracks
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_DeltaR          , "L1TkTau_MatchTk_DeltaR"        , 2000,  +0.0,   +2.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_PtRel           , "L1TkTau_MatchTk_PtRel"         , 2000,  +0.0,  +20.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_Pt              , "L1TkTau_MatchTk_Pt"            ,  300,  +0.0, +300.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_Eta             , "L1TkTau_MatchTk_Eta"           ,  600,  -3.0,   +3.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_POCAz           , "L1TkTau_MatchTk_POCAz"         ,  600, -30.0,  +30.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_d0              , "L1TkTau_MatchTk_d0"            , 2000, -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_d0Abs           , "L1TkTau_MatchTk_d0Abs"         , 1000,   0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NStubs          , "L1TkTau_MatchTk_NStubs"        ,   30,  -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NPsStubs        , "L1TkTau_MatchTk_NPsStubs"      ,   30,  -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NBarrelStubs    , "L1TkTau_MatchTk_NBarrelStubs"  ,   30,  -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NEndcapStubs    , "L1TkTau_MatchTk_NEndcapStubs"  ,   30,  -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_StubPtCons      , "L1TkTau_MatchTk_StubPtCons"    ,  200,  +0.0, +200.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_ChiSquared      , "L1TkTau_MatchTk_ChiSquared"    ,  200,  +0.0, +200.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_RedChiSquared   , "L1TkTau_MatchTk_RedChiSquared" , 2000,  +0.0, +200.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsGenuine       , "L1TkTau_MatchTk_IsGenuine"     ,    2,  -0.5,   +1.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsUnknown       , "L1TkTau_MatchTk_IsUnknown"     ,    2,  -0.5,   +1.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsCombinatoric  , "L1TkTau_MatchTk_IsCombinatoric",    2,  -0.5,   +1.5);


  // SingleTau
  histoTools_.BookHisto_1D(hCalo_Rate  , "Calo_Rate"   , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Rate    , "Tk_Rate"     , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Rate, "VtxIso_Rate" , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_Eff   , "Calo_Eff"    , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Eff     , "Tk_Eff"      , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Eff , "VtxIso_Eff"  , 200, 0.0,  +200.0);

  // DiTau
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo  , "DiTau_Rate_Calo"  , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk    , "DiTau_Rate_Tk"    , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso, "DiTau_Rate_VtxIso", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo   , "DiTau_Eff_Calo"   , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk     , "DiTau_Eff_Tk"     , 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso , "DiTau_Eff_VtxIso" , 200, 0.0,  +200.0);

  // Turn-Ons
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt, "McHadronicTau_VisEt" , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_TurnOn50      , "Calo_TurnOn50"       , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn50        , "Tk_TurnOn50"         , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50    , "VtxIso_TurnOn50"     , 40, 0.0,  +200.0);

  histoTools_.BookHisto_1D(hCalo_TurnOn25  , "Calo_TurnOn25"   , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn25    , "Tk_TurnOn25"     , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25, "VtxIso_TurnOn25" , 40, 0.0,  +200.0);
  
  histoTools_.BookHisto_1D(hCalo_TurnOn_SingleTau50KHz  , "Calo_TurnOn_SingleTau50KHz"  , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn_SingleTau50KHz    , "Tk_TurnOn_SingleTau50KHz"    , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_SingleTau50KHz, "VtxIso_TurnOn_SingleTau50KHz", 40, 0.0,  +200.0);

  histoTools_.BookHisto_1D(hCalo_TurnOn_DiTau50KHz  , "Calo_TurnOn_DiTau50KHz"  , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn_DiTau50KHz    , "Tk_TurnOn_DiTau50KHz"    , 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_DiTau50KHz, "VtxIso_TurnOn_DiTau50KHz", 40, 0.0,  +200.0);

  // DiTau (Calo-Other)
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_Tk    , "DiTau_Rate_Calo_Tk"    , 200, 0.0,  +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_VtxIso, "DiTau_Rate_Calo_VtxIso", 200, 0.0,  +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_Tk     , "DiTau_Eff_Calo_Tk"     , 200, 0.0,  +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_VtxIso , "DiTau_Eff_Calo_VtxIso" , 200, 0.0,  +200.0, 200, 0.0,  +200.0);

  // DiTau (Tk-Other)
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_VtxIso, "DiTau_Rate_Tk_VtxIso", 200, 0.0, +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_VtxIso , "DiTau_Eff_Tk_VtxIso" , 200, 0.0, +200.0, 200, 0.0,  +200.0);

  return;
}


//****************************************************************************
void CaloPlusTracks::FinaliseEffHisto_(TH1D *histo, 
				   const int nEvtsTotal)
//****************************************************************************
{

  const int nBins = histo->GetNbinsX()+1;
  double eff, err;

  // For-loop: Histogram bins
  for (int i = 0; i<= nBins; i++){
    
    const int nPass = histo->GetBinContent(i);
    auxTools_.Efficiency(nPass, nEvtsTotal, "binomial", eff, err );

    // Update current histo bin to true eff value and error
    histo->SetBinContent(i, eff);
    histo->SetBinError  (i, err);
  }

  return;
}


//****************************************************************************
void CaloPlusTracks::FinaliseEffHisto_(TH2D *histo, 
				   const int nEvtsTotal)
//****************************************************************************
{

  const int nBinsX  = histo->GetNbinsX()+1;
  const int nBinsY  = histo->GetNbinsY()+1;
  double eff, err;
  
  // For-loop: x-axis bins
  for (int bx=0; bx <= nBinsX; bx++){

    // For-loop: y-axis bins
    for (int by=0; by <= nBinsY; by++){

      const int nPass = histo->GetBinContent(bx, by);
      auxTools_.Efficiency(nPass, nEvtsTotal, "binomial", eff, err );

      // Update current histo bin to true eff value and error
      histo->SetBinContent(bx, by, eff);
      histo->SetBinError  (bx, by, err);

    } // For-loop: y-axis bins

  }// For-loop: x-axis bins

  return;
}


//****************************************************************************
void CaloPlusTracks::ApplyDiTauZMatching(string tkCollectionType, 
					 vector<L1TkTauParticle> &L1TkTaus)
//****************************************************************************
{
  
  // Sanity check
  if (L1TkTaus.size() < 2) return;


  // Initialise variables
  double deltaPOCAz = 9999.9;
  TTTrack match_tk0 = L1TkTaus.at(0).GetMatchingTk();
      
  // For-loop: L1TkTaus
  for (size_t i = 1; i < L1TkTaus.size(); i++)
    {

      TTTrack match_tk = L1TkTaus.at(i).GetMatchingTk();
      


    if ( tkCollectionType.compare("TTPixelTracks") == 0 )
      {
	cout << "E R R O R ! CaloPlusTracks::ApplyDiTauZMatching(...) - FIXME. EXIT" << endl;
	exit(1);
	
    }
    else if ( tkCollectionType.compare("TTTracks") == 0 ) {
      deltaPOCAz = abs( match_tk0.getZ0() - match_tk.getZ0() );
    }
    else{
      cout << "E R R O R ! CaloPlusTracks::ApplyDiTauZMatching(...) - Unknown sample \"" << mcSample << "\". EXIT" << endl;
      exit(1);
    }
    
    // If the Trigger objects is not within x-cm reject it
    if (deltaPOCAz > diTau_deltaPOCAz) L1TkTaus.erase ( L1TkTaus.begin()+i );
    
    }  // For-loop: L1TkTaus
  
  return;
}


//****************************************************************************
void CaloPlusTracks::FillSingleTau_(vector<L1TkTauParticle> L1TkTaus, 
				    TH1D *hRate,
				    TH1D *hEfficiency)
//****************************************************************************
{

  // Sanity check
  if( L1TkTaus.size() < 1 ) return;
  
  // Fill rate
  double ldgEt = L1TkTaus.at(0).GetCaloTau().et();
  FillRate_(hRate, ldgEt);
  
  // Get MC-matched trigger objects
  vector<L1TkTauParticle> L1TkTaus_mcMatched = GetMcMatchedL1TkTaus(L1TkTaus);
  if (L1TkTaus_mcMatched.size() < 1) return;
  
  // Check that all taus were found
  if(!bFoundAllTaus_) return;
  
  // Fill efficiency
  double ldgEt_mcMatched = L1TkTaus_mcMatched.at(0).GetCaloTau().et();
  FillEfficiency_(hEfficiency, ldgEt_mcMatched);

  return;
}


//****************************************************************************
void CaloPlusTracks::FillDiTau_(vector<L1TkTauParticle> L1TkTaus, 
				TH1D *hRate,
				TH1D *hEfficiency)
//****************************************************************************
{

  // Sanity check
  if( L1TkTaus.size() < 2 ) return;  

  // Fill rate
  L1TkTauParticle L1TkTau = L1TkTaus.at(1);
  double subLdgEt = L1TkTau.GetCaloTau().et();
  FillRate_(hRate, subLdgEt);

  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1TkTaus_mcMatched = GetMcMatchedL1TkTaus(L1TkTaus);
  if (L1TkTaus_mcMatched.size() < 2) return;
    
  // Check that all taus were found
  if(!bFoundAllTaus_) return;

  // Fill efficiency
  double subLdgEt_mcMatched = L1TkTaus_mcMatched.at(1).GetCaloTau().et();
  FillEfficiency_(hEfficiency, subLdgEt_mcMatched);

  return;
}


//****************************************************************************
void CaloPlusTracks::FillDiTau_(vector<L1TkTauParticle> L1TkTaus1,
				vector<L1TkTauParticle> L1TkTaus2, 
				TH2D *hRate,
				TH2D *hEfficiency)
//****************************************************************************
{

  // Sanity check
  if( L1TkTaus1.size() < 1 ) return;
  if( L1TkTaus2.size() < 1 ) return;
  
  vector<L1TkTauParticle> L1TkTaus1_ = L1TkTaus1;
  
  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1TkTaus_mcMatched1 = GetMcMatchedL1TkTaus(L1TkTaus1);
  vector<L1TkTauParticle> L1TkTaus_mcMatched2 = GetMcMatchedL1TkTaus(L1TkTaus1);

  // Fill rate 
  // RemoveDuplicates(L1TkTaus2, L1TkTaus1_); // FIXME - OBSOLETE
  if( L1TkTaus1_.size() < 1 ) return;
  
  double ldgEt1 = L1TkTaus1.at(0).GetCaloTau().et();
  double ldgEt2 = L1TkTaus2.at(0).GetCaloTau().et();

  // Make x-axis the ldgEt axis
  if (ldgEt1 > ldgEt2)      FillRate_(hRate, ldgEt1, ldgEt2); 
  else if (ldgEt1 < ldgEt2) FillRate_(hRate, ldgEt2, ldgEt1);
  else{
    if (L1TkTaus_mcMatched1.at(0).GetCaloTau().index() == L1TkTaus_mcMatched2.at(0).GetCaloTau().index() ) {
    cout << "E R R O R ! CaloPlusTracks::FillDiTau_(...) - Leading and Sub-Leading Trigger objects have the same eT (" << ldgEt1 << " , " << ldgEt2 << ")! "
	 << "This cannot happen. Exit." << endl;
    // exit(1); // FIXME 
    } 
  }

  // Get MC-matched trigger objects
  if (L1TkTaus_mcMatched1.size() < 1) return;
  if (L1TkTaus_mcMatched2.size() < 1) return;
  // RemoveDuplicates(L1TkTaus_mcMatched2, L1TkTaus_mcMatched1); // FIXME - OBSOLETE
  // if (L1TkTaus_mcMatched1.size() < 1) return; // FIXME - OBSOLETE

  double ldgEt_mcMatched1 = L1TkTaus_mcMatched1.at(1).GetCaloTau().et();
  double ldgEt_mcMatched2 = L1TkTaus_mcMatched2.at(1).GetCaloTau().et();
    
  // Check that all taus were found
  if(!bFoundAllTaus_) return;
  
  if (ldgEt_mcMatched1 > ldgEt_mcMatched2)      histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt_mcMatched1, ldgEt_mcMatched2);
  else if (ldgEt_mcMatched1 < ldgEt_mcMatched2) histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt_mcMatched2, ldgEt_mcMatched1);
  else{
    if (L1TkTaus_mcMatched1.at(0).GetCaloTau().index() == L1TkTaus_mcMatched2.at(0).GetCaloTau().index() ) {
      cout << "E R R O R ! CaloPlusTracks::FillDiTau_(...) - Leading and Sub-Leading MC-Matched Trigger objects have the same eT (" << ldgEt_mcMatched1 << " , " << ldgEt_mcMatched2 << ")! "
	   << "This cannot happen. Exit." << endl;
      // exit(1); // FIXME 
    }
  }

  return;
}


//****************************************************************************
void CaloPlusTracks::RemoveDuplicates(vector<L1TkTauParticle> L1TkTaus2, 
				      vector<L1TkTauParticle> &L1TkTaus1)
//****************************************************************************
{

  // Sanity check
  if (L1TkTaus1.size() < 1) return;
  if (L1TkTaus2.size() < 1) return;

  if ( L1TkTaus1.size() < L1TkTaus2.size() )
    {
      cout << "E R R O R ! CaloPlusTracks::RemoveDuplicates(...) - Size of L1TkTaus1 (" << L1TkTaus1.size() << ") is smaller than the size of L1TkTaus2 (" << L1TkTaus2.size() << "). "
	   <<  "Perhaps you should reverse their order when passing them as arguments to this function? EXIT" << endl;
      exit(1);
    }

  int calo_index1 = L1TkTaus1.at(0).caloTau_Index_;
  int calo_index2 = L1TkTaus2.at(0).caloTau_Index_;
  if( calo_index1 == calo_index2) L1TkTaus1.erase ( L1TkTaus1.begin()+0 );

  return;
}


//****************************************************************************
void CaloPlusTracks::FillRate_(TH1D *hRate,
			   const double ldgEt)
//****************************************************************************
{
  
  if (ldgEt < 0) return;
  hRate ->Fill( ldgEt );
  
  return;
}


//****************************************************************************
void CaloPlusTracks::FillRate_(TH2D *hRate,
			   const double ldgEt1,
			   const double ldgEt2)
//****************************************************************************
{
  
  if (ldgEt1 < 0) return;
  if (ldgEt2 < 0) return;

  hRate ->Fill( ldgEt1, ldgEt2 );
  
  return;
}


//****************************************************************************
void CaloPlusTracks::FillEfficiency_(TH1D *hEfficiency,
				 const double ldgEt)
//****************************************************************************
{
  
  histoTools_.FillAllBinsUpToValue_1D(hEfficiency, ldgEt);

  return;
}



//****************************************************************************
vector<GenParticle> CaloPlusTracks::GetHadronicGenTaus(vector<GenParticle> GenTaus,
						       double visEt,
						       double visEta)
//****************************************************************************
{

  // Sanity check
  vector<GenParticle> hadGenTaus;
  if (GenTaus.size() < 1 ) return hadGenTaus;

  
  // For-loop: GenTaus
  for (vector<GenParticle>::iterator tau = GenTaus.begin(); tau != GenTaus.end(); tau++)
    {

      if (0) tau->PrintProperties();

      // Get hadronic decay products (pi+/-,pi0, K+/-, K0, K0L, KOS, eta, omegas, gammas from tau->tau+gamma transition)
      std::vector<unsigned short> hadronicDaughters;
      std::cout << "1 ================================================" << std::endl;
      vector<unsigned short> test = tau->daughtersIndex();
      auxTools_.PrintVector(test);
      std::cout << "1 ================================================" << std::endl;
      mcTools_.GetHadronicTauFinalDaughters(tau->index(), hadronicDaughters); //fixme: iro

      // Sanity check
      if (hadronicDaughters.size() < 1) continue;

      // Acceptance cuts using visible 4-momenta
      TLorentzVector tau_visP4 = GetVisibleP4(hadronicDaughters);
      bool bPassVisEt  = ( tau_visP4.Et() >= visEt );
      bool bPassVisEta = ( std::abs(tau_visP4.Eta()) <= visEta );
      if (!(bPassVisEt * bPassVisEta)) continue;

      // Save this hadronic generator tau
      if (0) tau->PrintProperties();
      hadGenTaus.push_back(*tau);
    }
  
  return hadGenTaus;
}      



//****************************************************************************
void CaloPlusTracks::FillTurnOn_Numerator_(const vector<L1TkTauParticle> L1TkTaus, 
					   const double minEt,
					   TH1D *hTurnOn)
//****************************************************************************
{

  // Sanity check
  if (L1TkTaus.size() < 1) return;

  // For-loop: L1TkTaus
  for(int iTau = 0; iTau < (int) L1TkTaus.size(); iTau++){

    // Skip if trigger object is not MC matched
    L1TkTauParticle L1TkTau = L1TkTaus.at(iTau);
    if (!L1TkTau.HasMatchingGenParticle()) continue;    

    // Skip if trigger object has eT < minEt
    double calo_et = L1TkTau.GetCaloTau().et();
    if ( calo_et < minEt) continue;

    // Fill histo with visible eT of mathing MC Tau
    GenParticle match_GenP = L1TkTau.GetMatchingGenParticle();
    vector<unsigned short> tau_daughters;
    std::cout << "2 ================================================" << std::endl;
    vector<unsigned short> test = match_GenP.daughtersIndex();
    auxTools_.PrintVector(test);
    std::cout << "2 ================================================" << std::endl;
    mcTools_.GetHadronicTauFinalDaughters(match_GenP.index(), tau_daughters); //fixme: iro
    TLorentzVector tau_visP4 = GetVisibleP4(tau_daughters);
    auxTools_.PrintVector(tau_daughters);
    
    std::cout << "Filling Turn-On Histo with Et = " << tau_visP4.Et() << std::endl;
    hTurnOn->Fill( tau_visP4.Et() );
    
  } // For-loop: L1TkTaus
  
  return;
   
}


//****************************************************************************
void CaloPlusTracks::FillTurnOn_Denominator_(vector<GenParticle> GenTausHadronic,
					     TH1D *hVisEt)
//****************************************************************************
{

  if (GenTausHadronic.size() < 1) return;

  // For-loop: GenTausHadronic
  for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++)
    {
      if (0) tau->PrintProperties();

      // Get hadronic decay products (pi+/-,pi0, K+/-, K0, K0L, KOS, eta, omegas, gammas from tau->tau+gamma transition)
      std::vector<unsigned short> hadronicDaughters;
      std::cout << "3 ================================================" << std::endl;
      vector<unsigned short> test = tau->daughtersIndex();
      auxTools_.PrintVector(test);
      std::cout << "3 ================================================" << std::endl;
      mcTools_.GetHadronicTauFinalDaughters(tau->index(), hadronicDaughters); //fixme: iro

      // Sanity check
      if (hadronicDaughters.size() < 1) cout << "E R R O R ! CaloPlusTracks::FillTurnOn_Denominator_() - This should never be reached. EXIT" << endl;

      // Acceptance cuts using visible 4-momenta
      TLorentzVector tau_visP4 = GetVisibleP4(hadronicDaughters);
      hVisEt->Fill( tau_visP4.Et() );

    }// For-loop: GenTausHadronic

  return;

}


//****************************************************************************
void CaloPlusTracks::GetMatchingTrack(L1TkTauParticle &L1TkTau,
				      L1JetParticle L1CaloTau,
				      vector<TTTrack> TTTracks)

//****************************************************************************
{

  // Initialise variables
  TTTrack matchTk;
  double matchTk_dR = 999.9;
    
  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
    {
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), L1CaloTau.eta(), L1CaloTau.phi());

      // Only consider tracks within matching cone (0 <= DeltaR < matchTk_dR_max)
      if (dR > L1TkTau.GetMatchConeMax()) continue;

      // Find the closest matching tracks
      if (dR < matchTk_dR) 
	{
	  matchTk = *tk;
	  matchTk_dR = dR;
	}
    }

  // Assign values to the L1TkTau
  L1TkTau.SetCaloTau(L1CaloTau);
  L1TkTau.SetMatchingTk(matchTk);
  L1TkTau.SetMatchTkDeltaRNew(matchTk_dR);
  if (0) L1TkTau.PrintProperties(false, false, true, false);
  
  return;
}


//****************************************************************************
void CaloPlusTracks::GetSigConeTracks(L1TkTauParticle &L1TkTau,
				      vector<TTTrack> TTTracks)
//****************************************************************************
{
  if (!L1TkTau.HasMatchingTk()) return; 
  vector<TTTrack> sigConeTks;
  
  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
    {
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), L1TkTau.GetMatchingTk().getEta(), L1TkTau.GetMatchingTk().getPhi());

      // Only consider tracks within singal cone
      if (dR < L1TkTau.GetSigConeMin()) continue; // FIXME - Matching Track might NOT be added
      if (dR > L1TkTau.GetSigConeMax()) continue;
      
      sigConeTks.push_back(*tk);
    }
  
  L1TkTau.SetSigConeTracks(sigConeTks);
  if (0) L1TkTau.PrintProperties(false, false, true, false);
  
  return;
}


//****************************************************************************
void CaloPlusTracks::GetIsoConeTracks(L1TkTauParticle &L1TkTau,
				      vector<TTTrack> TTTracks)
//****************************************************************************
{
  if (!L1TkTau.HasMatchingTk()) return; 
  vector<TTTrack> isoConeTks;
  
  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
    {
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), L1TkTau.GetMatchingTk().getEta(), L1TkTau.GetMatchingTk().getPhi());

      // Only consider tracks within singal cone
      if (dR < L1TkTau.GetIsoConeMin()) continue;
      if (dR > L1TkTau.GetIsoConeMax()) continue;

      isoConeTks.push_back(*tk);
    }
  
  L1TkTau.SetIsoConeTracks(isoConeTks);

  return;
}


//****************************************************************************
vector<GenParticle> CaloPlusTracks::GetGenParticles(void)
//****************************************************************************
{
  vector<GenParticle> theGenParticles;
  for (Size_t genP_index = 0; genP_index < GenP_Pt->size(); genP_index++) theGenParticles.push_back( GetGenParticle(genP_index) );
  return theGenParticles;
}



//****************************************************************************
vector<GenParticle> CaloPlusTracks::GetGenParticles(int pdgId, bool isLastCopy)
//****************************************************************************
{

  // First get all the genParticles
  vector<GenParticle> allGenParticles = GetGenParticles();
  vector<GenParticle> myGenParticles;

  // For-loop: GenParticles
  for (vector<GenParticle>::iterator p = allGenParticles.begin(); p != allGenParticles.end(); p++)
    {

      // Apply criteria
      int genP_pdgId  = p->pdgId();
      if ( abs(genP_pdgId) != pdgId) continue;
      if (!isLastCopy) myGenParticles.push_back(*p);
      else
	{

	  // Determine if it's a last copy
	  bool save = true;
	  vector<unsigned short> genP_Daughters = p->daughtersIndex();
	  
	  // For-loop: Daughters
	  for (size_t i = 0; i < genP_Daughters.size(); i++)
	    {

	      unsigned short index = genP_Daughters.at(i);
	      GenParticle d        = GetGenParticle(index);
	      
	      if ( abs(d.pdgId()) == abs(pdgId) )
		{
		  save = false;
		  break;
		}
	    }
	  
	  if (!save) continue;	  
	  myGenParticles.push_back(*p);
	}
    }

  return myGenParticles;
}


//****************************************************************************
GenParticle CaloPlusTracks::GetGenParticle(unsigned int Index)
//****************************************************************************
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
  // theGenParticle.SetMothers(Mothers);
  // theGenParticle.SetDaughters(Daughters);
  
  return theGenParticle;
}



//****************************************************************************
double CaloPlusTracks::GetPVTTTracks(vector<TTTrack> &pvTTTracks)
//****************************************************************************
{
  // First get PV_z and pv-tracks indices
  vector<int> pvTks_index; 
  pv_z = 0.0;
  pv_z = pvProducer->GetPrimaryVertexZ("default", pvTks_index);
  
  // Uset tk-indices to get the actual pv-tracks
  for (Size_t iTk = 0; iTk < pvTks_index.size(); iTk++) pvTTTracks.push_back( GetTTTrack(iTk) );
  return pv_z;
}


//****************************************************************************
vector<TTTrack> CaloPlusTracks::GetTTTracks(const double minPt,
					    const double maxEta,
					    const double maxChiSq,
					    const unsigned int minStubs,
					    const unsigned int minStubsPS,
					    const unsigned nFitParams)
//****************************************************************************
{
  vector<TTTrack> theTTTracks; 
  for (Size_t iTk = 0; iTk < L1Tks_Pt->size(); iTk++)
    {
      TTTrack tk = GetTTTrack(iTk, nFitParams);
      
      if (tk.getPt() < minPt) continue;
      if (std::abs(tk.getEta()) > maxEta) continue;
      if (tk.getChi2() > maxChiSq) continue;
      if (tk.getNumOfStubs() < minStubs) continue;
      if (tk.getNumOfStubsPS() < minStubsPS) continue;
      // double z0 = tk.getZ0();
      // double d0 = tk.getD0();

      theTTTracks.push_back( tk );
    }
  return theTTTracks;
}


//****************************************************************************
TTTrack CaloPlusTracks::GetTTTrack(unsigned int Index,
				   const unsigned int nFitParams)
//****************************************************************************
{
  
  // Initialise variables
  TVector3 p;  

  // Get the track properties
  double pt  = L1Tks_Pt->at(Index);
  double eta = L1Tks_Eta->at(Index);
  double phi = L1Tks_Phi->at(Index);
  double x0  = L1Tks_POCAx->at(Index);
  double y0  = L1Tks_POCAy->at(Index);
  double z0  = L1Tks_POCAz->at(Index);
  p.SetPtEtaPhi(pt, eta, phi);
  ROOT::Math::XYZVector aPOCA(x0, y0, z0);
  double aRInv                      = L1Tks_RInv->at(Index);
  double aChi2                      = L1Tks_ChiSquared->at(Index);
  double aStubPtCons                = L1Tks_StubPtConsistency->at(Index);
  double aSector                    = L1Tks_Sector->at(Index);
  double aWedge                     = L1Tks_Wedge->at(Index);
  bool isGenuine                    = L1Tks_IsGenuine->at(Index);
  bool isUnknown                    = L1Tks_IsUnknown->at(Index);
  bool isCombinatoric               = L1Tks_IsCombinatoric->at(Index);
  vector<unsigned int> stubs_isPS   = L1Tks_Stubs_isPS->at(Index);
  vector<unsigned int> stubs_iDisk  = L1Tks_Stubs_iDisk->at(Index);
  vector<unsigned int> stubs_iLayer = L1Tks_Stubs_iLayer->at(Index);
  vector<unsigned int> stubs_iPhi   = L1Tks_Stubs_iPhi->at(Index);
  vector<unsigned int> stubs_iRing  = L1Tks_Stubs_iRing->at(Index);
  vector<unsigned int> stubs_iSide  = L1Tks_Stubs_iSide->at(Index);
  vector<unsigned int> stubs_iZ     = L1Tks_Stubs_iZ->at(Index);
  // auxTools_.PrintVector(stubs_isPS);
  
  // Construct the TTTrack
  TTTrack theTTTrack(Index, p, aPOCA,  aRInv, aChi2, aStubPtCons,
		     aSector, aWedge, isGenuine, isUnknown, isCombinatoric,
		     stubs_isPS, stubs_iDisk, stubs_iLayer, stubs_iPhi, stubs_iRing, stubs_iSide, stubs_iZ, nFitParams);

  return theTTTrack;
}



//****************************************************************************
vector<TrackingParticle> CaloPlusTracks::GetTrackingParticles(void)
//****************************************************************************
{
  vector<TrackingParticle> theTrackingParticles; 
  for (Size_t iTP = 0; iTP < TP_Pt->size(); iTP++) theTrackingParticles.push_back( GetTrackingParticle(iTP) );
  return theTrackingParticles;
}


//****************************************************************************
TrackingParticle CaloPlusTracks::GetTrackingParticle(unsigned int Index)
//****************************************************************************
{

  // Initialise variables
  TVector3 p;  

  // Get the track properties
  double pt  = TP_Pt->at(Index);
  double eta = TP_Eta->at(Index);
  double phi = TP_Phi->at(Index);
  double x0  = TP_POCAx->at(Index);
  double y0  = TP_POCAy->at(Index);
  double z0  = TP_POCAz->at(Index);
  int pdgId  = TP_PdgId->at(Index);
  int charge = TP_Charge->at(Index);
  p.SetPtEtaPhi(pt, eta, phi);
  ROOT::Math::XYZVector poca(x0, y0, z0);
  unsigned short nMatch       = TP_NMatch->at(Index);
  unsigned short ttTrackIndex = TP_TTTrackIndex->at(Index);
  unsigned short ttClusters   = TP_TTClusters->at(Index);
  unsigned short ttStubs      = TP_TTStubs->at(Index);
  unsigned short ttTracks     = TP_TTTracks->at(Index);

  // Construct the TTTrack
  TrackingParticle theTrackingParticle(Index, p, poca,  charge, pdgId, nMatch, ttTrackIndex, ttClusters, ttStubs, ttTracks);

  return theTrackingParticle;
}



//****************************************************************************
vector<TTPixelTrack> CaloPlusTracks::GetTTPixelTracks(const double minPt,
						      const double maxEta,
						      const double maxChiSq,
						      const int minHits)
//****************************************************************************
{
  vector<TTPixelTrack> theTTPixelTracks; 
  for (Size_t iTk = 0; iTk < L1PixTks_Pt->size(); iTk++)
    {
      TTPixelTrack pixTk = GetTTPixelTrack(iTk);
      // double z0 = tk.getZ0();
      // double d0 = tk.getD0();
      if (pixTk.getPt() < minPt) continue;
      if (std::abs(pixTk.getEta()) > maxEta) continue;
      if (pixTk.getChi2() > maxChiSq) continue;
      if (pixTk.getNcandidatehit() < minHits) continue;

      theTTPixelTracks.push_back( pixTk );
    }
  return theTTPixelTracks;
}


//****************************************************************************
TTPixelTrack CaloPlusTracks::GetTTPixelTrack(unsigned int Index)
//****************************************************************************
{

  // Initialise variables
  TVector3 theMomentum;
  ROOT::Math::XYZVector thePOCA;
  
  // Get the track properties
  double pt           = L1PixTks_Pt->at(Index);
  double eta          = L1PixTks_Eta->at(Index);
  double phi          = L1PixTks_Phi->at(Index);
  double x0           = L1PixTks_POCAx->at(Index);
  double y0           = L1PixTks_POCAy->at(Index);
  double z0           = L1PixTks_POCAz->at(Index);
  double theRInv      = L1PixTks_RInv->at(Index);
  double theChi2      = L1PixTks_ChiSquared->at(Index);
  double theSigmaRInv = L1PixTks_SigmaRInv->at(Index);
  double theSigmaPhi0 = L1PixTks_SigmaPhi0->at(Index);
  double theSigmaD0   = L1PixTks_SigmaD0->at(Index);
  double theSigmaT    = L1PixTks_SigmaT->at(Index);
  double theSigmaZ0   = L1PixTks_SigmaZ0->at(Index);
  int iTTTrack        = L1PixTks_TTTrackIndex->at(Index);
  TTTrack theTTTrack  = GetTTTrack(iTTTrack, 5); // TTPixelTracks can only have 5 fit parameters
  theMomentum.SetPtEtaPhi(pt, eta, phi);
  thePOCA.SetXYZ(x0, y0, z0);

  // Pixel Hits  
  int nPixHits = L1PixTks_NPixHits->at(Index);
  vector<ROOT::Math::XYZVector> thePixHits;
  for(int iHit = 0; iHit < nPixHits; iHit++)
    {
      double hitX = L1PixTks_PixHits_X->at(Index).at(iHit);
      double hitY = L1PixTks_PixHits_Y->at(Index).at(iHit);
      double hitZ = L1PixTks_PixHits_Z->at(Index).at(iHit);
      // double hitR    = L1PixTks_PixHits_R->at(Index).at(iHit);
      // double hitPhi  = L1PixTks_PixHits_Phi->at(Index).at(iHit);
      // double hitType = L1PixTks_PixHits_Type->at(Index).at(iHit);      

      ROOT::Math::XYZVector pixHit_XYZ( hitX, hitY, hitZ);
      thePixHits.push_back( pixHit_XYZ );
    }

  // Candidate Pixel Hits  
  int nCandPixHits = L1PixTks_CandPixHits_X->at(Index).size();
  vector<ROOT::Math::XYZVector> theCandPixHits;
  for(int iHit = 0; iHit < nCandPixHits; iHit++)
    {
      double hitX    = L1PixTks_CandPixHits_X->at(Index).at(iHit);
      double hitY    = L1PixTks_CandPixHits_Y->at(Index).at(iHit);
      double hitZ    = L1PixTks_CandPixHits_Z->at(Index).at(iHit);
      // double hitR    = L1PixTks_CandPixHits_R->at(Index).at(iHit);
      // double hitPhi  = L1PixTks_CandPixHits_Phi->at(Index).at(iHit);
      // double hitType = L1PixTks_CandPixHits_Type->at(Index).at(iHit);      

      ROOT::Math::XYZVector candPixHit_XYZ( hitX, hitY, hitZ);
      theCandPixHits.push_back( candPixHit_XYZ );
    }

  // Construct the TTPixelTrack
  TTPixelTrack theTTPixelTrack(Index, theTTTrack, theMomentum, thePOCA, theRInv, theChi2, theSigmaRInv, theSigmaPhi0, theSigmaD0, theSigmaT, theSigmaZ0, thePixHits, theCandPixHits);

  return theTTPixelTrack;
}


//****************************************************************************
vector<L1JetParticle> CaloPlusTracks::GetL1CaloTaus(void)
//****************************************************************************
{
  vector<L1JetParticle> theL1CaloTaus;
  for (Size_t iCalo = 0; iCalo < L1CaloTau_E->size(); iCalo++) theL1CaloTaus.push_back( GetL1CaloTau(iCalo) );
  return theL1CaloTaus;
}


//****************************************************************************
L1JetParticle CaloPlusTracks::GetL1CaloTau(unsigned int Index)
//****************************************************************************
{

  double E    = L1CaloTau_E->at(Index);
  double Et   = L1CaloTau_Et->at(Index);
  double Eta  = L1CaloTau_Eta->at(Index);
  double Phi  = L1CaloTau_Phi->at(Index);
  double Bx   = L1CaloTau_Bx->at(Index);
  double Type = L1CaloTau_Type->at(Index);
  L1JetParticle theL1CaloTau(Index, E, Et, Eta, Phi, Bx, Type); 
  return theL1CaloTau;
}

//****************************************************************************
void CaloPlusTracks::GetShrinkingConeSizes(double calo_et,
					   double sigCone_Constant,
					   double isoCone_Constant,
					   const double sigCone_dRCutoff,
					   double &sigCone_dRMin,
					   double &sigCone_dRMax,
					   double &isoCone_dRMin,
					   double &isoCone_dRMax)
//****************************************************************************
{
  
  double signalCone_min = (sigCone_Constant)/(calo_et);
  double signalCone_max = (isoCone_Constant)/(calo_et);
  if (signalCone_max > sigCone_dRCutoff) signalCone_max = sigCone_dRCutoff;
  else{}
  double isoCone_min    = signalCone_max;
  double isoCone_max    = (isoCone_dRMax/1.0);
      
  // Assign signal and isolation cone sizes
  sigCone_dRMin = signalCone_min;
  sigCone_dRMax = signalCone_max;  
  isoCone_dRMin = isoCone_min;
  isoCone_dRMax = isoCone_max;

  return;
}

//****************************************************************************
void CaloPlusTracks::GetIsolationValues(L1TkTauParticle &L1TkTau)
//****************************************************************************
{

  // Store default values
  L1TkTau.SetVtxIsolation(999.9);
  L1TkTau.SetRelIsolation(0.0);

  // Return if CaloTau is not Tk-Confirmed
  if (!L1TkTau.HasMatchingTk()) return; 

  // If no tracks found in the isoalation cone return
  vector<TTTrack> isoConeTks = L1TkTau.GetIsoConeTTTracks();
  if ( (isoConeTks.size() < 1) )  return;

  // Initialise variables
  TTTrack matchTk            = L1TkTau.GetMatchingTk();
  double isoTks_scalarSumPt  = 0.0;
  double deltaZ0             = 999.9;
  double relIso              = 0.0;
  
  // For-loop: All Tracks
  for (size_t i = 0; i < isoConeTks.size(); i++)
    {
      TTTrack isoConeTk = isoConeTks.at(i);
      
      // Add-up the pT of alltracks in isolation cone/annulus
      isoTks_scalarSumPt += isoConeTk.getPt();
      
      // Find the track closest in Z0
      deltaZ0 = std::abs(matchTk.getZ0() - isoConeTk.getZ0());
      if (deltaZ0 < L1TkTau.GetVtxIsolation())
	{
	  L1TkTau.SetVtxIsolation(deltaZ0);
	  L1TkTau.SetVtxIsolationTrack(isoConeTk);
	}
      
    }

  // Calculated + Assign value of relative isolation
  relIso = isoTks_scalarSumPt/matchTk.getPt();
  L1TkTau.SetRelIsolation(relIso);
  
  return;
}


//****************************************************************************
void CaloPlusTracks::GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
					    vector<GenParticle> GenTaus)					    
//****************************************************************************
{

  // Sanity check
  if (!L1TkTau.HasMatchingTk() ) return;
  if (GenTaus.size() < 1 ) return;
  
  // Initialise
  TTTrack matchTk  = L1TkTau.GetMatchingTk();
  double match_dR  = 9999.9;
  GenParticle match_GenParticle;

  // For-loop: GenTaus
  for (vector<GenParticle>::iterator tau = GenTaus.begin(); tau != GenTaus.end(); tau++)
    {

      if (0) tau->PrintProperties();

      // Get hadronic decay products (pi+/-,pi0, K+/-, K0, K0L, KOS, eta, omegas, gammas from tau->tau+gamma transition)
      std::vector<unsigned short> hadronicDaughters;
      std::cout << "4 ================================================" << std::endl;
      vector<unsigned short> test = tau->daughtersIndex();
      auxTools_.PrintVector(test);
      std::cout << "4 ================================================" << std::endl;
      mcTools_.GetHadronicTauFinalDaughters(tau->index(), hadronicDaughters); //fixme: iro

      // Sanity check
      if (hadronicDaughters.size() < 1) continue;
      
      // For-loop: Daughters
      vector<GenParticle> Daughters;
      for (size_t j = 0; j < hadronicDaughters.size(); j++)
       	{
	  GenParticle d = GetGenParticle( hadronicDaughters.at(j) );
	  if (d.charge() == 0) continue;
	  if (0) d.PrintProperties();

	  // Calculate deltaR from matching tracks
	  double deltaR = auxTools_.DeltaR( matchTk.getEta(), matchTk.getPhi(), d.eta(), d.phi() );
	  if (deltaR > mcMatching_dRMax) continue;

	  if (deltaR < match_dR)
	    {
	      match_dR = deltaR;
	      match_GenParticle = d;
	    }
	}            
    }
  
  // Save the matching
  L1TkTau.SetMatchingGenParticle(match_GenParticle);
  L1TkTau.SetMatchingGenParticleDeltaR(match_dR);

  if (0) match_GenParticle.PrintProperties();
  if (0) L1TkTau.PrintProperties(false, true, false, false);
  return;
}      

//****************************************************************************
vector<L1TkTauParticle> CaloPlusTracks::GetMcMatchedL1TkTaus(vector<L1TkTauParticle> L1TkTaus)
//****************************************************************************
{

  // Get MC-matched trigger objects above given Et Threshold
  vector<L1TkTauParticle> matchedL1TkTaus;
  for (vector<L1TkTauParticle>::iterator tau = L1TkTaus.begin(); tau != L1TkTaus.end(); tau++)
    {
      if (tau->HasMatchingGenParticle()) continue;
      matchedL1TkTaus.push_back(*tau);
    }
  
  return matchedL1TkTaus;
}

#endif

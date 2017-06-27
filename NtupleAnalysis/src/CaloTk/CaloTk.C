#ifndef CaloTk_cxx
#define CaloTk_cxx

// User
#include "../Auxiliary/interface/constants.h"
#include "CaloTk.h"

// ROOT
#include "TFitResult.h"
#include "TF1.h"
#include "Math/VectorUtil.h"

// C++
#include <stdexcept>

//============================================================================
void CaloTk::InitObjects(void)
//============================================================================
{

  // pvProducer = new L1TkPrimaryVertex(this->s);

  return; 
}


//============================================================================
void CaloTk::InitVars_()
//============================================================================
{

  DEBUG = false;

  // Dataset-related
  datasets_  = datasets_.GetDataset(mcSample);
  realTauMom = datasets_.McTauMomPdgId_;
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;

  // Matching tracks
  matchTk_Collection  =  "TTTracks"; // TP: "TTTracks" (not "TTPixelTracks")
  matchTk_nFitParams  =   5;         // TP:   5
  matchTk_minPt       =   5.00;      // TP:   5.0
  matchTk_minEta      =   0.0;       // TP:   0.0
  matchTk_maxEta      = 999.9;       // TP: 999.9  
  matchTk_maxChiSqRed = 200.0;       // TP: 200.0
  matchTk_minStubs    =   0;         // TP:   0
  matchTk_minStubsPS  =   0;         // TP:   0
  matchTk_maxStubsPS  = 999;         // TP: 999
  matchTk_caloDeltaR  =   0.10;      // TP:   0.10

  // Signal cone tracks
  sigConeTks_Collection  = matchTk_Collection; // TP: "TTTracks" (not "TTPixelTracks")
  sigConeTks_nFitParams  = matchTk_nFitParams; // TP:   5
  sigConeTks_minPt       =   2.0;              // TP:   2.0
  sigConeTks_minEta      =   0.0;              // TP:   0.0
  sigConeTks_maxEta      = 999.9;              // TP: 999.9
  sigConeTks_maxChiSqRed = 200.0;              // TP: 200.0
  sigConeTks_minStubs    =   0;                // TP:   0
  sigConeTks_minStubsPS  =   0;                // TP:   0
  sigConeTks_maxStubsPS  = 999;                // TP: 999

  // Isolation cone tracks
  isoConeTks_Collection  = matchTk_Collection;    // TP: "TTTracks" (not "TTPixelTracks")
  isoConeTks_nFitParams  = matchTk_nFitParams;    // TP:   5
  isoConeTks_minPt       =   2.0;                 // TP:   2.0
  isoConeTks_minEta      =   0.0;                 // TP:   0.0
  isoConeTks_maxEta      = 999.9;                 // TP: 999.9
  isoConeTks_maxChiSqRed = 200.0;                 // TP: 200.0
  isoConeTks_minStubs    =   0;                   // TP:   0
  isoConeTks_minStubsPS  =   0;                   // TP:   0
  isoConeTks_maxStubsPS  = 999;                   // TP: 999
  
  // Signal cone parameters
  sigCone_Constant        = +0.00; // TP: 0.00
  sigCone_dRMin           = +0.00; // WARNING! If > 0 the matching Track will NOT be added in sigCone_TTTracks
  sigCone_dRMax           = +0.15; // TP: 0.15
  sigCone_cutoffDeltaR    = +0.15; // TP: 0.15
  sigCone_maxTkInvMass    = +1.77; // TP: Unused (3-pr)
  sigCone_maxTkDeltaPOCAz = +0.20; // TP: Unused (3-pr)

  // Isolation cone
  isoCone_Constant = +3.50;         // TP: 3.50 GeV
  isoCone_VtxIsoWP = +0.50;         // TP: 1.0 cm
  isoCone_dRMin    = sigCone_dRMax; // TP: 0.4
  isoCone_dRMax    = +0.30;         // TP: 0.4
  diTau_deltaPOCAz = +1.00;         // TP: 1.0 cm

  // MC matching
  mcMatching_dRMax  = +0.10;        // TP: 0.05
  mcMatching_unique = true;

  PrintSettings();

  return;
}


//============================================================================
void CaloTk::PrintSettings(void)
//============================================================================
{

  // Inform user of settings
  Table settings("Variable | Cut | Value | TP 2015 | Units", "Text");  // Table settingsTable("Variable & Value & Units", "LaTeX", "l l l");
  settings.AddRowColumn(0, "MC Sample");
  settings.AddRowColumn(0, "==" );
  settings.AddRowColumn(0, mcSample );

  settings.AddRowColumn(1, "Matching Tracks: Collection");
  settings.AddRowColumn(1, "==");
  settings.AddRowColumn(1, matchTk_Collection);
  settings.AddRowColumn(1, "TTTracks");
  settings.AddRowColumn(1, "");
  
  settings.AddRowColumn(2, "Matching Tracks: Fit Parameters");
  settings.AddRowColumn(2, "==");
  settings.AddRowColumn(2, auxTools_.ToString( matchTk_nFitParams) );
  settings.AddRowColumn(2, "5");
  settings.AddRowColumn(2, "");

  settings.AddRowColumn(3, "Matching Tracks: Pt");
  settings.AddRowColumn(3, ">=");
  settings.AddRowColumn(3, auxTools_.ToString( matchTk_minPt) );
  settings.AddRowColumn(3, "2" );
  settings.AddRowColumn(3, "GeV/c" );
  
  settings.AddRowColumn(4, "Matching Tracks: |Eta|");
  settings.AddRowColumn(4, ">=");
  settings.AddRowColumn(4, auxTools_.ToString( matchTk_minEta) );
  settings.AddRowColumn(4, "0.0" );
  settings.AddRowColumn(4, "" );

  settings.AddRowColumn(5, "Matching Tracks: |Eta|");
  settings.AddRowColumn(5, "<=");
  settings.AddRowColumn(5, auxTools_.ToString( matchTk_maxEta) );
  settings.AddRowColumn(5, "1e+03" );
  settings.AddRowColumn(5, "" );
  
  settings.AddRowColumn(6, "Matching Tracks: ChiSqRed");
  settings.AddRowColumn(6, "<=");
  settings.AddRowColumn(6, auxTools_.ToString( matchTk_maxChiSqRed) );
  settings.AddRowColumn(6, "200/DOF"); // Cut was on ChiSq, not ChiSqRed
  settings.AddRowColumn(6, "");

  settings.AddRowColumn(7, "Matching Tracks: Stubs");
  settings.AddRowColumn(7, ">=");
  settings.AddRowColumn(7, auxTools_.ToString( matchTk_minStubs) );
  settings.AddRowColumn(7, "0" );
  settings.AddRowColumn(7, "" );
  
  settings.AddRowColumn(8, "Matching Tracks: PS-Stubs (min)");
  settings.AddRowColumn(8, ">=");
  settings.AddRowColumn(8, auxTools_.ToString( matchTk_minStubsPS) );
  settings.AddRowColumn(8, "" );
  settings.AddRowColumn(8, "" );

  settings.AddRowColumn(9, "Matching Tracks: PS-Stubs (max)");
  settings.AddRowColumn(9, "<=");
  settings.AddRowColumn(9, auxTools_.ToString( matchTk_maxStubsPS) );
  settings.AddRowColumn(9, "N/A" );
  settings.AddRowColumn(9, "" );

  settings.AddRowColumn(10, "Matching Tracks: DeltaR");
  settings.AddRowColumn(10, "<=");
  settings.AddRowColumn(10, auxTools_.ToString(matchTk_caloDeltaR) );
  settings.AddRowColumn(10, "0.10" );
  settings.AddRowColumn(10, "");

  settings.AddRowColumn(11, "Signal Cone Tks: Collection");
  settings.AddRowColumn(11, "==");
  settings.AddRowColumn(11, sigConeTks_Collection);
  settings.AddRowColumn(11, "TTTracks");
  settings.AddRowColumn(11, "");
  
  settings.AddRowColumn(12, "Signal Cone Tks: Fit Parameters");
  settings.AddRowColumn(12, "==");
  settings.AddRowColumn(12, auxTools_.ToString( sigConeTks_nFitParams) );
  settings.AddRowColumn(12, "5");
  settings.AddRowColumn(12, "");

  settings.AddRowColumn(13, "Signal Cone Tks: Pt");
  settings.AddRowColumn(13, ">=");
  settings.AddRowColumn(13, auxTools_.ToString( sigConeTks_minPt) );
  settings.AddRowColumn(13, "2" );
  settings.AddRowColumn(13, "GeV/c" );
  
  settings.AddRowColumn(14, "Signal Cone Tks: |Eta|");
  settings.AddRowColumn(14, ">=");
  settings.AddRowColumn(14, auxTools_.ToString( sigConeTks_minEta) );
  settings.AddRowColumn(14, "0.0" );
  settings.AddRowColumn(14, "" );

  settings.AddRowColumn(15, "Signal Cone Tks: |Eta|");
  settings.AddRowColumn(15, "<=");
  settings.AddRowColumn(15, auxTools_.ToString( sigConeTks_maxEta) );
  settings.AddRowColumn(15, "1e+03" );
  settings.AddRowColumn(15, "" );
  
  settings.AddRowColumn(16, "Signal Cone Tks: ChiSqRed");
  settings.AddRowColumn(16, "<=");
  settings.AddRowColumn(16, auxTools_.ToString( sigConeTks_maxChiSqRed) );
  settings.AddRowColumn(16, "200 (but on ChiSq, not ChiSqRed)");
  settings.AddRowColumn(16, "");

  settings.AddRowColumn(17, "Signal Cone Tks: Stubs");
  settings.AddRowColumn(17, ">=");
  settings.AddRowColumn(17, auxTools_.ToString( sigConeTks_minStubs) );
  settings.AddRowColumn(17, "" );
  settings.AddRowColumn(17, "" );
  
  settings.AddRowColumn(18, "Signal Cone Tks: PS-Stubs (min)");
  settings.AddRowColumn(18, ">=");
  settings.AddRowColumn(18, auxTools_.ToString( sigConeTks_minStubsPS) );
  settings.AddRowColumn(18, "0" );
  settings.AddRowColumn(18, "" );

  settings.AddRowColumn(19, "Signal Cone Tks: PS-Stubs (max)");
  settings.AddRowColumn(19, "<=");
  settings.AddRowColumn(19, auxTools_.ToString( sigConeTks_maxStubsPS) );
  settings.AddRowColumn(19, "N/A" );
  settings.AddRowColumn(19, "" );

  settings.AddRowColumn(20, "Isolation Cone Tks: Collection");
  settings.AddRowColumn(20, "==");
  settings.AddRowColumn(20, isoConeTks_Collection);
  settings.AddRowColumn(20, "TTTracks");
  settings.AddRowColumn(20, "");
  
  settings.AddRowColumn(21, "Isolation Cone Tks: Fit Parameters");
  settings.AddRowColumn(21, "==");
  settings.AddRowColumn(21, auxTools_.ToString( isoConeTks_nFitParams) );
  settings.AddRowColumn(21, "5");
  settings.AddRowColumn(21, "");

  settings.AddRowColumn(22, "Isolation Cone Tks: Pt");
  settings.AddRowColumn(22, ">=");
  settings.AddRowColumn(22, auxTools_.ToString( isoConeTks_minPt) );
  settings.AddRowColumn(22, "2" );
  settings.AddRowColumn(22, "GeV/c" );
  
  settings.AddRowColumn(23, "Isolation Cone Tks: |Eta|");
  settings.AddRowColumn(23, ">=");
  settings.AddRowColumn(23, auxTools_.ToString( isoConeTks_minEta) );
  settings.AddRowColumn(23, "0.0" );
  settings.AddRowColumn(23, "" );

  settings.AddRowColumn(24, "Isolation Cone Tks: |Eta|");
  settings.AddRowColumn(24, "<=");
  settings.AddRowColumn(24, auxTools_.ToString( isoConeTks_maxEta) );
  settings.AddRowColumn(24, "1e+03" );
  settings.AddRowColumn(24, "" );

  settings.AddRowColumn(25, "Isolation Cone Tks: ChiSqRed");
  settings.AddRowColumn(25, "<=");
  settings.AddRowColumn(25, auxTools_.ToString( isoConeTks_maxChiSqRed) );
  settings.AddRowColumn(25, "200 (but on ChiSq, not ChiSqRed)");
  settings.AddRowColumn(25, "");

  settings.AddRowColumn(26, "Isolation Cone Tks: Stubs");
  settings.AddRowColumn(26, ">=");
  settings.AddRowColumn(26, auxTools_.ToString( isoConeTks_minStubs) );
  settings.AddRowColumn(26, "" );
  settings.AddRowColumn(26, "" );
  
  settings.AddRowColumn(27, "Isolation Cone Tks: PS-Stubs (min)");
  settings.AddRowColumn(27, ">=");
  settings.AddRowColumn(27, auxTools_.ToString( isoConeTks_minStubsPS) );
  settings.AddRowColumn(27, "0" );
  settings.AddRowColumn(27, "" );

  settings.AddRowColumn(28, "Isolation Cone Tks: PS-Stubs (max)");
  settings.AddRowColumn(28, "<=");
  settings.AddRowColumn(28, auxTools_.ToString( isoConeTks_maxStubsPS) );
  settings.AddRowColumn(28, "N/A" );
  settings.AddRowColumn(28, "" );

  settings.AddRowColumn(29, "Signal Cone: Shrink Constant");
  settings.AddRowColumn(29, "==");
  settings.AddRowColumn(29, auxTools_.ToString(sigCone_Constant) );
  settings.AddRowColumn(29, "0" );
  settings.AddRowColumn(29, "GeV");

  settings.AddRowColumn(30, "Signal Cone: DeltaR");
  settings.AddRowColumn(30, ">=");
  settings.AddRowColumn(30, auxTools_.ToString(sigCone_dRMin) );
  settings.AddRowColumn(30, "0.0" );
  settings.AddRowColumn(30, "" );

  settings.AddRowColumn(31, "Signal Cone: DeltaR");
  settings.AddRowColumn(31, "<=");
  settings.AddRowColumn(31, auxTools_.ToString(sigCone_dRMax) );
  settings.AddRowColumn(31, "0.15" );
  settings.AddRowColumn(31, "" );
  
  settings.AddRowColumn(32, "Signal Cone:-3pr InvMass");
  settings.AddRowColumn(32, "<=");
  settings.AddRowColumn(32, auxTools_.ToString(sigCone_maxTkInvMass) );
  settings.AddRowColumn(32, "N/A" );
  settings.AddRowColumn(32, "GeV/c^{-2}");

  settings.AddRowColumn(33, "Signal Cone:-3pr maxTkDeltaPOCAz");
  settings.AddRowColumn(33, "<=");
  settings.AddRowColumn(33, auxTools_.ToString(sigCone_maxTkDeltaPOCAz) );
  settings.AddRowColumn(33, "N/A" );
  settings.AddRowColumn(33, "cm");

  // settings.AddRowColumn(19, "");

  settings.AddRowColumn(34, "Isolation Cone: Shrink Constant");
  settings.AddRowColumn(34, "==");
  settings.AddRowColumn(34, auxTools_.ToString(isoCone_Constant) );
  settings.AddRowColumn(34, "3.5");
  settings.AddRowColumn(34, "GeV");

  settings.AddRowColumn(35, "Isolation Cone: DeltaR");
  settings.AddRowColumn(35, ">=");
  settings.AddRowColumn(35, auxTools_.ToString(isoCone_dRMin) );
  settings.AddRowColumn(35, "0.15" );
  settings.AddRowColumn(35, "" );

  settings.AddRowColumn(36, "Isolation Cone: DeltaR");
  settings.AddRowColumn(36, "=<");
  settings.AddRowColumn(36, auxTools_.ToString(isoCone_dRMax) );
  settings.AddRowColumn(36, "0.30");
  settings.AddRowColumn(36, "");

  settings.AddRowColumn(37, "Isolation Cone: VtxIso" );
  settings.AddRowColumn(37, "<=" );
  settings.AddRowColumn(37, auxTools_.ToString(isoCone_VtxIsoWP) );
  settings.AddRowColumn(37, "1.0");
  settings.AddRowColumn(37, "cm");
  settings.AddRowColumn(37, "");

  settings.AddRowColumn(38, "Di-Tau |Delta z0|");
  settings.AddRowColumn(38, "<");
  settings.AddRowColumn(38, auxTools_.ToString(diTau_deltaPOCAz) );
  settings.AddRowColumn(38, "1.0" );
  settings.AddRowColumn(38, "cm");

  settings.AddRowColumn(39, "MC-Matching DeltaR");
  settings.AddRowColumn(39, "<=");
  settings.AddRowColumn(39, auxTools_.ToString(mcMatching_dRMax) );
  settings.AddRowColumn(39, "0.05" );
  settings.AddRowColumn(39, "" );

  settings.AddRowColumn(40, "MC-Matching IsUnique");
  settings.AddRowColumn(40, "==");
  settings.AddRowColumn(40, auxTools_.ToString(mcMatching_unique) );
  settings.AddRowColumn(40, "1" );
  settings.AddRowColumn(40, "" );
  
  settings.AddRowColumn(41, "MC-Taus: Mom PdgId");
  settings.AddRowColumn(41, "==");
  settings.AddRowColumn(41, auxTools_.ToString(realTauMom));
  settings.AddRowColumn(41, "N/A" );
  settings.AddRowColumn(41, "" );

  settings.AddRowColumn(42, "MC-Taus: Number Expected");
  settings.AddRowColumn(42, ">=");
  settings.AddRowColumn(42, auxTools_.ToString(nMaxNumOfHTausPossible));
  settings.AddRowColumn(42, "N/A" );
  settings.AddRowColumn(42, "" );

  settings.AddRowColumn(43, "" );
  settings.Print();
  
  return;
}


//============================================================================
void CaloTk::Loop()
//============================================================================
{

  // Sanity check
  if (fChain == 0) return;
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  cout << "=== CaloTk:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;
  // Initialisations
  InitVars_();
  BookHistos_();
  Long64_t nbytes       = 0;
  Long64_t nb           = 0;
  int nEvtsWithMaxHTaus = 0; 
  unsigned int nEvts    = 0;
  unsigned int nAllEvts = fChain->GetEntries();
  bool isMinBias        = false;  
  // L1PixelTrackFit f(3.8112); // Bz in Tesla (for pixel re-fitting)

  // Determine what sample this is
  std::size_t found = mcSample.find("SingleNeutrino");
    if (found!=std::string::npos)
      {
	isMinBias = true;
	std::cout << "Minimum Bias sample" << std::endl;
      }
    else
      {
	std::cout << "Not a Minimum Bias sample." << std::endl;
      }
    
  
  ////////////////////////////////////////////////
  // For-loop: Entries
  ////////////////////////////////////////////////
  for (int jentry = 0; jentry < nEntries; jentry++, nEvts++){
  
    if(DEBUG) cout << "\tEntry = " << jentry << endl;
    
    // Init variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    // Gen-Collections
    vector<GenParticle> GenParticles;
    vector<GenParticle> GenTaus;
    vector<GenParticle> GenTausHadronic;
    vector<GenParticle> GenTausTrigger;

    // Get the GenParticles
    if (!isMinBias)
      {
	if(DEBUG) cout << "\tGetting the GenParticles" << endl;
	if (0) GenParticles = GetGenParticles(false); // time-consuming
	GenTaus             = GetGenParticles(15, true);
	GenTausHadronic     = GetHadronicGenTaus(GenTaus, 00.0, 999.9);
	GenTausTrigger      = GetHadronicGenTaus(GenTaus, 20.0, 2.4);
      }
    
    if (DEBUG)
       {
	cout << "\tPrinting all GenParticle Collections" << endl;
	if (0) PrintGenParticleCollection(GenParticles);
	PrintGenParticleCollection(GenTaus);
	PrintGenParticleCollection(GenTausHadronic);    
	PrintGenParticleCollection(GenTausTrigger);
      }

    // Track Collections
    vector<TrackingParticle> TPs = GetTrackingParticles(false);
    vector<TTTrack> matchTTTracks = GetTTTracks(matchTk_minPt, matchTk_minEta, matchTk_maxEta, matchTk_maxChiSqRed,
						matchTk_minStubs, matchTk_minStubsPS, matchTk_maxStubsPS, matchTk_nFitParams, false);
    
    vector<TTTrack> sigTTTracks = GetTTTracks(sigConeTks_minPt , sigConeTks_minEta, sigConeTks_maxEta, sigConeTks_maxChiSqRed,
					      sigConeTks_minStubs, sigConeTks_minStubsPS, sigConeTks_maxStubsPS, sigConeTks_nFitParams, false);
    
    vector<TTTrack> isoTTTracks = GetTTTracks(isoConeTks_minPt , isoConeTks_minEta, isoConeTks_maxEta, isoConeTks_maxChiSqRed,
					      isoConeTks_minStubs, isoConeTks_minStubsPS, isoConeTks_maxStubsPS, isoConeTks_nFitParams, false);

    if (DEBUG)
      {
	cout << "\tPrinting all TTrack Collections" << endl;
	// cout << "sigTTTracks.size() = " << sigTTTracks.size() << ", isoTTTracks.size() = " << isoTTTracks.size() << endl;
	PrintTrackingParticleCollection(TPs);
	PrintTTTrackCollection(matchTTTracks);
	PrintTTTrackCollection(sigTTTracks);
	PrintTTTrackCollection(isoTTTracks);
      }

    
    // Tau Collections
    vector<L1Tau> L1Taus = GetL1Taus(false);
    vector<L1TkTauParticle> L1TkTauCandidates;
    vector<L1TkTauParticle> L1TkTaus_Calo;
    vector<L1TkTauParticle> L1TkTaus_Tk;
    vector<L1TkTauParticle> L1TkTaus_VtxIso;    


    // alex: new-start
    vector<L1Tau> L1Jets = GetL1Jets(true);
    GetL1EGs(true);
    GetL1Sums(true);
    // alex: new-end

    // Ensure that all taus are found
    bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;

    // ======================================================================================
    // For-loop: L1Taus
    // ======================================================================================
    for (vector<L1Tau>::iterator calo = L1Taus.begin(); calo != L1Taus.end(); calo++)
      {
	// Calculate the Et-dependent signal & isolation cone sizes
	GetShrinkingConeSizes(calo->et(), sigCone_Constant, isoCone_Constant, sigCone_cutoffDeltaR,
			      sigCone_dRMin, sigCone_dRMax, isoCone_dRMin, isoCone_dRMax);

	// Construct the CaloTk candidate
	L1TkTauParticle	L1TkTauCandidate(0.0, matchTk_caloDeltaR, sigCone_dRMin, sigCone_dRMax, isoCone_dRMin, isoCone_dRMax);


	if ( matchTk_Collection.compare("TTTracks") != 0 )
	  {
	    cout << "=== ERROR: Invalid Track Collection Type \"" << matchTk_Collection << "\". EXIT" << endl;
	    exit(1);
	  }
	    
	GetMatchingTrack(L1TkTauCandidate, *calo, matchTTTracks);
	GetSigConeTracks(L1TkTauCandidate, sigTTTracks);
	GetIsoConeTracks(L1TkTauCandidate, isoTTTracks);
	GetIsolationValues(L1TkTauCandidate);
	GetMatchingGenParticle(L1TkTauCandidate, GenTausHadronic);
	if (DEBUG) L1TkTauCandidate.PrintProperties(false, false, true, true);
	
	// Save L1TkTau Candidate
	L1TkTauCandidates.push_back(L1TkTauCandidate);
      }
    
    ////////////////////////////////////////////////
    /// Create L1TkTaus Collections
    ////////////////////////////////////////////////
    for(vector<L1TkTauParticle>::iterator L1TkTau = L1TkTauCandidates.begin(); L1TkTau != L1TkTauCandidates.end(); L1TkTau++)
      {
	// Calo
	L1TkTaus_Calo.push_back(*L1TkTau);

	// +Tk and +VtxIso
	if (L1TkTau->HasMatchingTk() )
	  {
	    // std::cout << "\n=== L1TkTau->GetMatchingTkDeltaR() = " << L1TkTau->GetMatchingTkDeltaR() << std::endl;
	    L1TkTaus_Tk.push_back(*L1TkTau);
	    if (L1TkTau->GetVtxIsolation() > isoCone_VtxIsoWP) L1TkTaus_VtxIso.push_back(*L1TkTau);
	  }
      }// L1TkTauCandidates

    if (DEBUG)
      {
	PrintL1TauCollection(L1Taus);
	PrintL1TkTauParticleCollection(L1TkTauCandidates);
	PrintL1TkTauParticleCollection(L1TkTaus_Calo);
	PrintL1TkTauParticleCollection(L1TkTaus_Tk);
	PrintL1TkTauParticleCollection(L1TkTaus_VtxIso);
      }

    ////////////////////////////////////////////////
    // Event-Type Histograms
    ////////////////////////////////////////////////
    hHepMCEvt_VtxX_VtxY->Fill(HepMCEvt_VtxX, HepMCEvt_VtxY);
    hHepMCEvt_VtxZ->Fill(HepMCEvt_VtxZ);
    hL1TkPV_VtxZ->Fill(pv_z);
    hPrimaryVertex_DeltaZ->Fill(pv_z - HepMCEvt_VtxZ);
    hPrimaryVertex_AbsDeltaZ->Fill( abs(pv_z - HepMCEvt_VtxZ) );

    
    ////////////////////////////////////////////////
    /// L1TkTaus_Calo Properties 
    ////////////////////////////////////////////////
    for (vector<L1TkTauParticle>::iterator tau = L1TkTaus_Calo.begin(); tau != L1TkTaus_Calo.end(); tau++)
      {

	// L1TkTau Resolution
	GenParticle p = tau->GetMatchingGenParticle();	
	hL1Tau_ResolutionCaloEt ->Fill( (tau->GetCaloTau().eta() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	hL1Tau_ResolutionCaloEta->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	hL1Tau_ResolutionCaloPhi->Fill( (tau->GetCaloTau().eta() - p.p4vis().Phi())/p.p4vis().Phi() );

      }
    
    ////////////////////////////////////////////////
    /// L1Tktau Properties 
    ////////////////////////////////////////////////
    vector<L1TkTauParticle> myL1TkTaus = L1TkTaus_Tk; // L1TkTaus_Calo, L1TkTaus_Tk, L1TkTaus_VtxIso
    hL1TkTau_Multiplicity ->Fill( myL1TkTaus.size() );
    for (vector<L1TkTauParticle>::iterator tau = myL1TkTaus.begin(); tau != myL1TkTaus.end(); tau++)
      {
	
	if (DEBUG) tau->PrintProperties(true, true, true, true);
	
	// Variables
	TLorentzVector sigTks_p4 = tau->GetSigConeTTTracksP4();
	TLorentzVector isoTks_p4 = tau->GetIsoConeTTTracksP4();

	// Do not skip if using MinBias sample as no real taus exist!
	if (!tau->HasMatchingGenParticle() && (isMinBias == false) ) continue;
	
	// L1TkTau Resolution
	GenParticle p = tau->GetMatchingGenParticle(); //fixme: more plots from MC info

	// Resolution
	hL1TkTau_ResolutionCaloEt ->Fill( (tau->GetCaloTau().eta() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	hL1TkTau_ResolutionCaloEta->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	hL1TkTau_ResolutionCaloPhi->Fill( (tau->GetCaloTau().eta() - p.p4vis().Phi())/p.p4vis().Phi() );


	// Apply isolation?
	if (0) if (tau->GetVtxIsolation() > isoCone_VtxIsoWP) continue;

	
	// Matching Track Variables
	TTTrack matchTk   = tau->GetMatchingTk();
	double matchTk_dR = auxTools_.DeltaR(matchTk.getEta(), matchTk.getPhi(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi() );
	TLorentzVector caloTau_p4;
      	caloTau_p4.SetPtEtaPhiE(tau->GetCaloTau().et(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi(), tau->GetCaloTau().et() );
	hL1TkTau_MatchTk_DeltaR        ->Fill( matchTk_dR );
	hL1TkTau_MatchTk_PtRel         ->Fill( matchTk.p3().Perp(caloTau_p4.Vect()) );
	hL1TkTau_MatchTk_Pt            ->Fill( matchTk.getPt() );
	hL1TkTau_MatchTk_Eta           ->Fill( matchTk.getEta() );
	hL1TkTau_MatchTk_POCAz         ->Fill( matchTk.getZ0() );
	hL1TkTau_MatchTk_d0            ->Fill( matchTk.getD0() );
	hL1TkTau_MatchTk_d0Abs         ->Fill( abs(matchTk.getD0()) );
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
	hL1TkTau_MatchTk_PtMinusCaloEt ->Fill( matchTk.getPt() - tau->GetCaloTau().et() );

	// Signal/Isolation Cone Variables
	hL1TkTau_Rtau        ->Fill( tau->GetSigConeLdgTk().getPt() / tau->GetCaloTau().et() );
	hL1TkTau_CHF         ->Fill( tau->GetCaloTau().et()/sigTks_p4.Et() );
	hL1TkTau_NHF         ->Fill( (tau->GetCaloTau().et() - sigTks_p4.Et())/tau->GetCaloTau().et() );
	hL1TkTau_NHFAbs      ->Fill( abs( (tau->GetCaloTau().et() - sigTks_p4.Et())/tau->GetCaloTau().et() ) );
	hL1TkTau_NSigTks     ->Fill( tau->GetSigConeTTTracks().size() );
	hL1TkTau_NIsoTks     ->Fill( tau->GetIsoConeTTTracks().size() );
	if (tau->GetSigConeTTTracks().size() > 1) hL1TkTau_InvMass->Fill( sigTks_p4.M() );
	hL1TkTau_InvMassIncl ->Fill( sigTks_p4.M() );
	hL1TkTau_SigConeRMin ->Fill( tau->GetSigConeMin() );
	hL1TkTau_IsoConeRMin ->Fill( tau->GetIsoConeMin() );
	hL1TkTau_SigConeRMax ->Fill( tau->GetSigConeMax() );
	hL1TkTau_IsoConeRMax ->Fill( tau->GetIsoConeMax() );
	hL1TkTau_DeltaRGenP  ->Fill( tau->GetMatchingGenParticleDeltaR() );
	hL1TkTau_RelIso      ->Fill( tau->GetRelIsolation() );
	hL1TkTau_VtxIso      ->Fill( tau->GetVtxIsolation() );
	// hL1TkTau_VtxIsoAbs   ->Fill( abs(tau->GetVtxIsolation()) );
	  

	// SigCone TTTracks
	int sigTks_sumCharge   = 0;
	vector<TTTrack> sigTks = tau->GetSigConeTTTracks();
	for (vector<TTTrack>::iterator sigTk = sigTks.begin(); sigTk != sigTks.end(); sigTk++)
	  {
	    // Print properties?
	    if (DEBUG) sigTk->PrintProperties();
	    
	    // Get the transverse component of this track with respect to the matching track
	    TVector3 sigTk_p3  = sigTk->getMomentum();
	    double sigTk_PtRel = sigTk_p3.Perp( matchTk.getMomentum() );
	    double sigTk_dR    = auxTools_.DeltaR(tau->GetMatchingTk().getEta(), tau->GetMatchingTk().getPhi(), sigTk->getEta(), sigTk->getPhi());
	    
	    // Fill Histograms
	    hL1TkTau_SigTks_Pt        ->Fill( sigTk->getPt()  );
	    hL1TkTau_SigTks_PtRel     ->Fill( sigTk_PtRel );
	    hL1TkTau_SigTks_Eta       ->Fill( sigTk->getEta() );
	    hL1TkTau_SigTks_POCAz     ->Fill( sigTk->getZ0()  );
	    if (sigTks.size() > 1)
	      {
		hL1TkTau_SigTks_DeltaPOCAz->Fill( abs( sigTk->getZ0() - matchTk.getZ0() ) );
	      }
	    hL1TkTau_SigTks_d0        ->Fill( sigTk->getD0() );
	    hL1TkTau_SigTks_d0Abs     ->Fill( abs( sigTk->getD0()) );
	    // hL1TkTau_SigTks_d0Sig     ->Fill( sigTk->getD0()/sigTk->getSigmaD0() );        // TTPixelTracks
	    // hL1TkTau_SigTks_d0SigAbs  ->Fill( abs(sigTk->getD0()/sigTk->getSigmaD0() ) );  // TTPixelTracks
	    hL1TkTau_SigTks_d0Sig     ->Fill( -1.0 );
	    hL1TkTau_SigTks_d0SigAbs  ->Fill( -1.0 );
	    hL1TkTau_SigTks_StubPtCons->Fill( sigTk->getStubPtConsistency() );
	    hL1TkTau_SigTks_DeltaR    ->Fill( sigTk_dR );
	    hL1TkTau_SigTks_NStubs    ->Fill( sigTk->getNumOfStubs() );
	    hL1TkTau_SigTks_NPsStubs  ->Fill( sigTk->getNumOfStubsPS() );
	    hL1TkTau_SigTks_NBarrelStubs->Fill( sigTk->getNumOfBarrelStubs() );
	    hL1TkTau_SigTks_NEndcapStubs->Fill( sigTk->getNumOfEndcapStubs() );
	    hL1TkTau_SigTks_ChiSquared->Fill( sigTk->getChi2() );
	    hL1TkTau_SigTks_RedChiSquared->Fill( sigTk->getChi2Red() );
	    hL1TkTau_SigTks_PtMinusCaloEt->Fill( sigTk->getPt() - tau->GetCaloTau().et() );

	    // Other variables
	    sigTks_sumCharge += sigTk->getCharge();
	    
	  }// SigCone_TTTracks

	// Fill histos for other variables
	hL1TkTau_Charge->Fill( sigTks_sumCharge);
	hL1TkTau_ChargeAbs->Fill( abs(sigTks_sumCharge) );


	// IsoCone TTTracks
	vector<TTTrack> isoTks = tau->GetIsoConeTTTracks();	
	for (vector<TTTrack>::iterator isoTk = isoTks.begin(); isoTk != isoTks.end(); isoTk++)
	  {

	    // Print properties?
	    if (0) isoTk->PrintProperties();
	    
	    // Get the transverse component of this track with respect to the matching track
	    TVector3 isoTk_p3  = isoTk->getMomentum();
	    double isoTk_PtRel = isoTk_p3.Perp( matchTk.getMomentum() );
	    double isoTk_dR    = auxTools_.DeltaR(tau->GetMatchingTk().getEta(), tau->GetMatchingTk().getPhi(), isoTk->getEta(), isoTk->getPhi());
	    
	    // Fill Histograms
	    hL1TkTau_IsoTks_Pt        ->Fill( isoTk->getPt()  );
	    hL1TkTau_IsoTks_PtRel     ->Fill( isoTk_PtRel );
	    hL1TkTau_IsoTks_Eta       ->Fill( isoTk->getEta() );
	    hL1TkTau_IsoTks_POCAz     ->Fill( isoTk->getZ0()  );
	    hL1TkTau_IsoTks_DeltaPOCAz->Fill( abs( isoTk->getZ0() - matchTk.getZ0() ) );
	    hL1TkTau_IsoTks_d0        ->Fill( isoTk->getD0() );							     
	    hL1TkTau_IsoTks_d0Abs     ->Fill( abs( isoTk->getD0()) );
	    // hL1TkTau_IsoTks_d0Sig     ->Fill( isoTk->getD0()/isoTk->getSigmaD0() );        // TTPixelTracks
	    // hL1TkTau_IsoTks_d0SigAbs  ->Fill( abs(isoTk->getD0()/isoTk->getSigmaD0() ) );  // TTPixelTracks
	    hL1TkTau_IsoTks_d0Sig     ->Fill( -1.0 );
	    hL1TkTau_IsoTks_d0SigAbs  ->Fill( -1.0 );
	    hL1TkTau_IsoTks_StubPtCons->Fill( isoTk->getStubPtConsistency() );
	    hL1TkTau_IsoTks_DeltaR    ->Fill( isoTk_dR );
	    hL1TkTau_IsoTks_NStubs    ->Fill( isoTk->getNumOfStubs() );
	    hL1TkTau_IsoTks_NPsStubs  ->Fill( isoTk->getNumOfStubsPS() );
	    hL1TkTau_IsoTks_NBarrelStubs->Fill( isoTk->getNumOfBarrelStubs() );
	    hL1TkTau_IsoTks_NEndcapStubs->Fill( isoTk->getNumOfEndcapStubs() );
	    hL1TkTau_IsoTks_ChiSquared->Fill( isoTk->getChi2() );
	    hL1TkTau_IsoTks_RedChiSquared->Fill( isoTk->getChi2Red() );
	    hL1TkTau_IsoTks_PtMinusCaloEt->Fill( isoTk->getPt() - tau->GetCaloTau().et() );
  
	  }// IsoCone_TTTracks

      } // L1TkTaus_VtxIso

    
    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////
    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++) hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() );
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
    FillSingleTau_(L1TkTaus_Calo, hCalo_Rate  , hCalo_Eff);
    FillSingleTau_(L1TkTaus_Calo, hCalo_Rate_C, hCalo_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_Calo, hCalo_Rate_I, hCalo_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_Calo, hCalo_Rate_F, hCalo_Eff_F, 1.6, 3.0); // 2.5 is max
    
    FillSingleTau_(L1TkTaus_Tk    , hTk_Rate  , hTk_Eff  );
    FillSingleTau_(L1TkTaus_Tk    , hTk_Rate_C, hTk_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_Tk    , hTk_Rate_I, hTk_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_Tk    , hTk_Rate_F, hTk_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkTaus_VtxIso, hVtxIso_Rate  , hVtxIso_Eff);
    FillSingleTau_(L1TkTaus_VtxIso, hVtxIso_Rate_C, hVtxIso_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_VtxIso, hVtxIso_Rate_I, hVtxIso_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_VtxIso, hVtxIso_Rate_F, hVtxIso_Eff_F, 1.6, 3.0); // 2.5 is max

    ////////////////////////////////////////////////
    // DiTau
    ////////////////////////////////////////////////
    FillDiTau_(L1TkTaus_Calo, L1TkTaus_Tk    , hDiTau_Rate_Calo_Tk    , hDiTau_Eff_Calo_Tk    );
    FillDiTau_(L1TkTaus_Calo, L1TkTaus_VtxIso, hDiTau_Rate_Calo_VtxIso, hDiTau_Eff_Calo_VtxIso);
    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_VtxIso, hDiTau_Rate_Tk_VtxIso  , hDiTau_Eff_Tk_VtxIso  );
    

    ////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////
    ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_Tk);
    ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_VtxIso); 

    FillDiTau_(L1TkTaus_Calo, hDiTau_Rate_Calo  , hDiTau_Eff_Calo);
    FillDiTau_(L1TkTaus_Calo, hDiTau_Rate_Calo_C, hDiTau_Eff_Calo_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_Calo, hDiTau_Rate_Calo_I, hDiTau_Eff_Calo_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_Calo, hDiTau_Rate_Calo_F, hDiTau_Eff_Calo_F, 1.6, 3.0);
	
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk  , hDiTau_Eff_Tk);
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_C, hDiTau_Eff_Tk_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_I, hDiTau_Eff_Tk_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_F, hDiTau_Eff_Tk_F, 1.6, 3.0);
    
    FillDiTau_(L1TkTaus_VtxIso, hDiTau_Rate_VtxIso  , hDiTau_Eff_VtxIso);
    FillDiTau_(L1TkTaus_VtxIso, hDiTau_Rate_VtxIso_C, hDiTau_Eff_VtxIso_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_VtxIso, hDiTau_Rate_VtxIso_I, hDiTau_Eff_VtxIso_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_VtxIso, hDiTau_Rate_VtxIso_F, hDiTau_Eff_VtxIso_F, 1.6, 3.0);
    
    // Progress bar
    if (!DEBUG) auxTools_.ProgressBar(jentry, nEntries, 100, 100);
    
  }// For-loop: Entries

  // Fill counters
  hCounters->SetBinContent(1, nAllEvts);
  hCounters->SetBinContent(2, nEvts);


  ////////////////////////////////////////////////
  // Convert/Finalise Histos
  ////////////////////////////////////////////////
  // SingleTau
  histoTools_.ConvertToRateHisto_1D(hCalo_Rate  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hCalo_Rate_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hCalo_Rate_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hCalo_Rate_F, nEntries);
  
  histoTools_.ConvertToRateHisto_1D(hTk_Rate  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hTk_Rate_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hTk_Rate_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hTk_Rate_F, nEntries);
      
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate_F, nEntries);
  
  FinaliseEffHisto_( hCalo_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hCalo_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hCalo_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hCalo_Eff_F, nEvtsWithMaxHTaus);
      
  FinaliseEffHisto_( hTk_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hTk_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hTk_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hTk_Eff_F, nEvtsWithMaxHTaus);
  
  FinaliseEffHisto_( hVtxIso_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIso_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIso_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIso_Eff_F, nEvtsWithMaxHTaus);

  // DiTau
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Calo  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Calo_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Calo_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Calo_F, nEntries);

  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk_F, nEntries);
  
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso_F, nEntries);
  
  FinaliseEffHisto_( hDiTau_Eff_Calo  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Calo_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Calo_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Calo_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hDiTau_Eff_Tk  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Tk_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Tk_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Tk_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hDiTau_Eff_VtxIso  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIso_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIso_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIso_F, nEvtsWithMaxHTaus);

  // DiTau (Calo-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_Tk    , nEntries);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_VtxIso, nEntries);
  FinaliseEffHisto_( hDiTau_Eff_Calo_Tk    , nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Calo_VtxIso, nEvtsWithMaxHTaus);;

  // DiTau (Tk-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_VtxIso, nEntries);
  FinaliseEffHisto_( hDiTau_Eff_Tk_VtxIso, nEvtsWithMaxHTaus);;

  // Turn-Ons: fixme. iro. alex
  // TEfficiency *pEff = 0;
  // pEff = new TEfficiency(*hCalo_TurnOn50_passed, *hMcHadronicTau_VisEt);
  // hCalo_TurnOn50 = (TH1D*) pEff->Clone();
  
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
  auxTools_.StopwatchStop(5, "minutes", "Total Time");

}


//============================================================================
void CaloTk::BookHistos_(void)
//============================================================================
{
  
  // Event-Type Histograms
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 2, 0.0, +2.0);
  histoTools_.BookHisto_1D(hL1TkPV_VtxZ            , "L1TkPV_VtxZ"            ,  "", 600, -30.0,  +30.0);
  histoTools_.BookHisto_1D(hHepMCEvt_VtxZ          , "HepMCEvt_VtxZ"          ,  "", 600, -30.0,  +30.0);
  histoTools_.BookHisto_2D(hHepMCEvt_VtxX_VtxY     , "HepMCEvt_VtxX_VtxY"     ,  "", 400,  -0.01,  +0.01, 400,  -0.01,  +0.01);
  histoTools_.BookHisto_1D(hPrimaryVertex_DeltaZ   , "PrimaryVertex_DeltaZ"   ,  "", 600, -30.0,  +30.0);
  histoTools_.BookHisto_1D(hPrimaryVertex_AbsDeltaZ, "PrimaryVertex_AbsDeltaZ",  "", 300, + 0.0,  +30.0);

  // VtxIsolated L1TkTaus
  histoTools_.BookHisto_1D(hL1TkTau_Multiplicity, "L1TkTau_Multiplicity", "",    30, -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_Rtau        , "L1TkTau_Rtau"        , "",   100,  0.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_CHF         , "L1TkTau_CHF"         , "",   500,  0.0,  +50.0);
  histoTools_.BookHisto_1D(hL1TkTau_NHF         , "L1TkTau_NHF"         , "",   200, -5.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_NHFAbs      , "L1TkTau_NHFAbs"      , "",   100,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_NSigTks     , "L1TkTau_NSigTks"     , "",    15, -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_NIsoTks     , "L1TkTau_NIsoTks"     , "",    15, -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_InvMass     , "L1TkTau_InvMass"     , "",   100,  0.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_InvMassIncl , "L1TkTau_InvMassIncl" , "",   100,  0.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_SigConeRMin , "L1TkTau_SigConeRMin" , "",   200,  0.0,   +2.0);
  histoTools_.BookHisto_1D(hL1TkTau_SigConeRMax , "L1TkTau_SigConeRMax" , "",   100,  0.0,   +1.0);
  histoTools_.BookHisto_1D(hL1TkTau_IsoConeRMin , "L1TkTau_IsoConeRMin" , "",   100,  0.0,   +1.0);
  histoTools_.BookHisto_1D(hL1TkTau_IsoConeRMax , "L1TkTau_IsoConeRMax" , "",   100,  0.0,   +1.0);
  histoTools_.BookHisto_1D(hL1TkTau_Charge      , "L1TkTau_Charge"      , "",    19, -9.5,   +9.5);
  histoTools_.BookHisto_1D(hL1TkTau_ChargeAbs   , "L1TkTau_ChargeAbs"   , "",    19, -9.5,   +9.5);
  histoTools_.BookHisto_1D(hL1TkTau_RelIso      , "L1TkTau_RelIso"      , "",   100,  0.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_VtxIso      , "L1TkTau_VtxIso"      , "",   300,  0.0,  +30.0);
  // histoTools_.BookHisto_1D(hL1TkTau_VtxIso      , "L1TkTau_VtxIso"      ,   1200,-30.0,  +30.0);
  // histoTools_.BookHisto_1D(hL1TkTau_VtxIsoAbs   , "L1TkTau_VtxIsoAbs"   ,    600,  0.0,  +30.0);
  histoTools_.BookHisto_1D(hL1TkTau_DeltaRGenP  , "L1TkTau_DeltaRGenP"  , "",   100,  0.0,   +1.0);
  //
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Pt           , "L1TkTau_SigTks_Pt"           , "",  300,   +0.0,  +300.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Eta          , "L1TkTau_SigTks_Eta"          , "",  600,   -3.0,    +3.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_POCAz        , "L1TkTau_SigTks_POCAz"        , "",  600,  -30.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_DeltaPOCAz   , "L1TkTau_SigTks_DeltaPOCAz"   , "",  600,   +0.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0           , "L1TkTau_SigTks_d0"           , "", 2000,  -10.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0Abs        , "L1TkTau_SigTks_d0Abs"        , "", 1000,    0.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0Sig        , "L1TkTau_SigTks_d0Sig"        , "",  400,  -20.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_d0SigAbs     , "L1TkTau_SigTks_d0SigAbs"     , "",  200,    0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_PtRel        , "L1TkTau_SigTks_PtRel"        , "",  200,   +0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_StubPtCons   , "L1TkTau_SigTks_StubPtCons"   , "",  200,   +0.0,  +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_DeltaR       , "L1TkTau_SigTks_DeltaR"       , "",  200,   +0.0,    +2.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_NStubs       , "L1TkTau_SigTks_NStubs"       , "",   30,   -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_NPsStubs     , "L1TkTau_SigTks_NPsStubs"     , "",   30,   -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_NBarrelStubs , "L1TkTau_SigTks_NBarrelStubs" , "",   30,   -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_NEndcapStubs , "L1TkTau_SigTks_NEndcapStubs" , "",   30,   -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_ChiSquared   , "L1TkTau_SigTks_ChiSquared"   , "",  200,   +0.0,  +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_RedChiSquared, "L1TkTau_SigTks_RedChiSquared", "", 2000,   +0.0,  +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_PtMinusCaloEt, "L1TkTau_SigTks_PtMinusCaloEt", "",  100, -100.0,  +100.0  );			   
  //
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Pt           , "L1TkTau_IsoTks_Pt"           , "",  300,    +0.0,  +300.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Eta          , "L1TkTau_IsoTks_Eta"          , "",  600,    -3.0,    +3.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_POCAz        , "L1TkTau_IsoTks_POCAz"        , "",  600,   -30.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_DeltaPOCAz   , "L1TkTau_IsoTks_DeltaPOCAz"   , "",  600,    +0.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0           , "L1TkTau_IsoTks_d0"           , "", 2000,   -10.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0Abs        , "L1TkTau_IsoTks_d0Abs"        , "", 1000,     0.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0Sig        , "L1TkTau_IsoTks_d0Sig"        , "", 4000,   -20.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_d0SigAbs     , "L1TkTau_IsoTks_d0SigAbs"     , "", 2000,     0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_PtRel        , "L1TkTau_IsoTks_PtRel"        , "", 2000,    +0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_StubPtCons   , "L1TkTau_IsoTks_StubPtCons"   , "",  200,    +0.0,  +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_DeltaR       , "L1TkTau_IsoTks_DeltaR"       , "",  200,   +0.0,     +2.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_NStubs       , "L1TkTau_IsoTks_NStubs"       , "",   30,   -0.5,    +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_NPsStubs     , "L1TkTau_IsoTks_NPsStubs"     , "",   30,   -0.5,    +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_NBarrelStubs , "L1TkTau_IsoTks_NBarrelStubs" , "",   30,   -0.5,    +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_NEndcapStubs , "L1TkTau_IsoTks_NEndcapStubs" , "",   30,   -0.5,    +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_ChiSquared   , "L1TkTau_IsoTks_ChiSquared"   , "",  200,   +0.0,   +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_RedChiSquared, "L1TkTau_IsoTks_RedChiSquared", "", 2000,   +0.0,   +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_PtMinusCaloEt, "L1TkTau_IsoTks_PtMinusCaloEt", "",  100, -100.0,   +100.0  );


  // VtxIsolated L1TkTaus, Signal TTTracks
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_DeltaR        , "L1TkTau_MatchTk_DeltaR"        , "",  200,   +0.0,   +2.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_PtRel         , "L1TkTau_MatchTk_PtRel"         , "",  200,   +0.0,  +20.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_Pt            , "L1TkTau_MatchTk_Pt"            , "",  150,   +0.0, +300.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_Eta           , "L1TkTau_MatchTk_Eta"           , "",  600,   -3.0,   +3.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_POCAz         , "L1TkTau_MatchTk_POCAz"         , "",  600,  -30.0,  +30.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_d0            , "L1TkTau_MatchTk_d0"            , "", 2000,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_d0Abs         , "L1TkTau_MatchTk_d0Abs"         , "", 1000,    0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NStubs        , "L1TkTau_MatchTk_NStubs"        , "",   30,   -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NPsStubs      , "L1TkTau_MatchTk_NPsStubs"      , "",   30,   -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NBarrelStubs  , "L1TkTau_MatchTk_NBarrelStubs"  , "",   30,   -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NEndcapStubs  , "L1TkTau_MatchTk_NEndcapStubs"  , "",   30,   -0.5,  +14.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_StubPtCons    , "L1TkTau_MatchTk_StubPtCons"    , "",  200,   +0.0, +200.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_ChiSquared    , "L1TkTau_MatchTk_ChiSquared"    , "",  200,   +0.0, +200.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_RedChiSquared , "L1TkTau_MatchTk_RedChiSquared" , "", 2000,   +0.0, +200.0);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsGenuine     , "L1TkTau_MatchTk_IsGenuine"     , "",    2,   -0.5,   +1.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsUnknown     , "L1TkTau_MatchTk_IsUnknown"     , "",    2,   -0.5,   +1.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsCombinatoric, "L1TkTau_MatchTk_IsCombinatoric", "",    2,   -0.5,   +1.5);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_PtMinusCaloEt , "L1TkTau_MatchTk_PtMinusCaloEt" , "",  100, -100.0, +100.0);

  // Resolutions 
  histoTools_.BookHisto_1D(hL1Tau_ResolutionCaloEt , "L1Tau_ResolutionCaloEt" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1Tau_ResolutionCaloEta, "L1Tau_ResolutionCaloEta", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1Tau_ResolutionCaloPhi, "L1Tau_ResolutionCaloPhi", ";#phi (rads);Events / %.2f rads", 200,  -10.0,  +10.0);

  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEt , "L1TkTau_ResolutionCaloEt" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEta, "L1TkTau_ResolutionCaloEta", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloPhi, "L1TkTau_ResolutionCaloPhi", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);

  // SingleTau
  histoTools_.BookHisto_1D(hCalo_Rate    , "Calo_Rate"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_Rate_C  , "Calo_Rate_C"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_Rate_I  , "Calo_Rate_I"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_Rate_F  , "Calo_Rate_F"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Rate      , "Tk_Rate"      , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Rate_C    , "Tk_Rate_C"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Rate_I    , "Tk_Rate_I"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Rate_F    , "Tk_Rate_F"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Rate  , "VtxIso_Rate"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Rate_C, "VtxIso_Rate_C", "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Rate_I, "VtxIso_Rate_I", "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Rate_F, "VtxIso_Rate_F", "", 200, 0.0,  +200.0);
  
  histoTools_.BookHisto_1D(hCalo_Eff     , "Calo_Eff"     , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_Eff_C   , "Calo_Eff_C"   , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_Eff_I   , "Calo_Eff_I"   , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_Eff_F   , "Calo_Eff_F"   , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Eff       , "Tk_Eff"       , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Eff_C     , "Tk_Eff_C"     , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Eff_I     , "Tk_Eff_I"     , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_Eff_F     , "Tk_Eff_F"     , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Eff   , "VtxIso_Eff"   , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Eff_C , "VtxIso_Eff_C" , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Eff_I , "VtxIso_Eff_I" , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_Eff_F , "VtxIso_Eff_F" , "", 200, 0.0,  +200.0);
  
  // DiTau
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo    , "DiTau_Rate_Calo"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo_C  , "DiTau_Rate_Calo_C"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo_I  , "DiTau_Rate_Calo_I"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo_F  , "DiTau_Rate_Calo_F"  , "", 200, 0.0,  +200.0);  
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk      , "DiTau_Rate_Tk"      , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk_C    , "DiTau_Rate_Tk_C"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk_I    , "DiTau_Rate_Tk_I"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk_F    , "DiTau_Rate_Tk_F"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso  , "DiTau_Rate_VtxIso"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso_C, "DiTau_Rate_VtxIso_C", "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso_I, "DiTau_Rate_VtxIso_I", "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso_F, "DiTau_Rate_VtxIso_F", "", 200, 0.0,  +200.0);
    
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo     , "DiTau_Eff_Calo"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo_C   , "DiTau_Eff_Calo_C"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo_I   , "DiTau_Eff_Calo_I"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo_F   , "DiTau_Eff_Calo_F"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk       , "DiTau_Eff_Tk"      , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk_C     , "DiTau_Eff_Tk_C"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk_I     , "DiTau_Eff_Tk_I"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk_F     , "DiTau_Eff_Tk_F"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso   , "DiTau_Eff_VtxIso"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso_C , "DiTau_Eff_VtxIso_C", "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso_I , "DiTau_Eff_VtxIso_I", "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso_F , "DiTau_Eff_VtxIso_F", "", 200, 0.0,  +200.0);
  
  // Turn-Ons
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt, "McHadronicTau_VisEt", "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hCalo_TurnOn50      , "Calo_TurnOn50"      , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn50        , "Tk_TurnOn50"        , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50    , "VtxIso_TurnOn50"    , "", 40, 0.0,  +200.0);

  histoTools_.BookHisto_1D(hCalo_TurnOn25  , "Calo_TurnOn25"   , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn25    , "Tk_TurnOn25"     , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25, "VtxIso_TurnOn25" , "", 40, 0.0,  +200.0);
  
  histoTools_.BookHisto_1D(hCalo_TurnOn_SingleTau50KHz  , "Calo_TurnOn_SingleTau50KHz"  , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn_SingleTau50KHz    , "Tk_TurnOn_SingleTau50KHz"    , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_SingleTau50KHz, "VtxIso_TurnOn_SingleTau50KHz", "", 40, 0.0,  +200.0);

  histoTools_.BookHisto_1D(hCalo_TurnOn_DiTau50KHz  , "Calo_TurnOn_DiTau50KHz"  , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTk_TurnOn_DiTau50KHz    , "Tk_TurnOn_DiTau50KHz"    , "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_DiTau50KHz, "VtxIso_TurnOn_DiTau50KHz", "", 40, 0.0,  +200.0);

  // DiTau (Calo-Other)
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_Tk    , "DiTau_Rate_Calo_Tk"    , "" ,200, 0.0,  +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_VtxIso, "DiTau_Rate_Calo_VtxIso", "" ,200, 0.0,  +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_Tk     , "DiTau_Eff_Calo_Tk"     , "" ,200, 0.0,  +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_VtxIso , "DiTau_Eff_Calo_VtxIso" , "" ,200, 0.0,  +200.0, 200, 0.0,  +200.0);

  // DiTau (Tk-Other)
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_VtxIso, "DiTau_Rate_Tk_VtxIso", "", 200, 0.0, +200.0, 200, 0.0,  +200.0);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_VtxIso , "DiTau_Eff_Tk_VtxIso" , "", 200, 0.0, +200.0, 200, 0.0,  +200.0);

  return;
}


//============================================================================
void CaloTk::FinaliseEffHisto_(TH1D *histo, 
				       const int nEvtsTotal)
//============================================================================
{

  const int nBins = histo->GetNbinsX()+1;
  double eff, err;

  // For-loop: Histogram bins
  for (int i = 0; i<= nBins; i++){
    
    const int nPass = histo->GetBinContent(i);
    auxTools_.Efficiency(nPass, nEvtsTotal, "binomial", eff, err ); //fixme: use TEfficiency?

    // Update current histo bin to true eff value and error
    histo->SetBinContent(i, eff);
    histo->SetBinError  (i, err);
  }

  return;
}


//============================================================================
void CaloTk::FinaliseEffHisto_(TH2D *histo, 
				       const int nEvtsTotal)
//============================================================================
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


//============================================================================
void CaloTk::ApplyDiTauZMatching(string tkCollectionType, 
					 vector<L1TkTauParticle> &L1TkTaus)
//============================================================================
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
	cout << "=== CaloTk::ApplyDiTauZMatching() - Unsupported track collection. Exit" << endl;
	exit(1);
      }
    else if ( tkCollectionType.compare("TTTracks") == 0 ) {
      deltaPOCAz = abs( match_tk0.getZ0() - match_tk.getZ0() );
    }
    else{
      cout << "=== CaloTk::ApplyDiTauZMatching() - Unknown sample \"" << mcSample << "\". EXIT" << endl;
      exit(1);
    }
    
    // If the Trigger objects is not within x-cm reject it
    if (deltaPOCAz > diTau_deltaPOCAz) L1TkTaus.erase ( L1TkTaus.begin()+i );
    
    }  // For-loop: L1TkTaus
  
  return;
}


//============================================================================
void CaloTk::FillSingleTau_(vector<L1TkTauParticle> L1TkTaus, 
			    TH1D *hRate,
			    TH1D *hEfficiency,
			    double minEta,
			    double maxEta)
//============================================================================
{

  // Sanity check
  if( L1TkTaus.size() < 1 ) return;
  
  // Fill rate
  double ldgEt = L1TkTaus.at(0).GetCaloTau().et();

  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(L1TkTaus.at(0).GetCaloTau().eta()) < minEta) return;
  if ( abs(L1TkTaus.at(0).GetCaloTau().eta()) > maxEta) return;
  // if ( abs(L1TkTaus.at(0).GetMatchingTk().getEta()) < minEta) return;
  // if ( abs(L1TkTaus.at(0).GetMatchingTk().getEta()) > maxEta) return;
    
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


//============================================================================
void CaloTk::FillDiTau_(vector<L1TkTauParticle> L1TkTaus, 
			TH1D *hRate,
			TH1D *hEfficiency,
			double minEta,
			double maxEta)
//============================================================================
{

  // Sanity check
  if( L1TkTaus.size() < 2 ) return;  

  // Fill rate
  L1TkTauParticle L1TkTau = L1TkTaus.at(1);
  double subLdgEt = L1TkTau.GetCaloTau().et();
  
  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(L1TkTaus.at(0).GetCaloTau().eta()) < minEta) return;
  if ( abs(L1TkTaus.at(0).GetCaloTau().eta()) > maxEta) return;
  if ( abs(L1TkTaus.at(1).GetCaloTau().eta()) < minEta) return;
  if ( abs(L1TkTaus.at(1).GetCaloTau().eta()) > maxEta) return;
  // if ( abs(L1TkTaus.at(0).GetMatchingTk().getEta()) < minEta) return;
  // if ( abs(L1TkTaus.at(0).GetMatchingTk().getEta()) > maxEta) return;
  // if ( abs(L1TkTaus.at(1).GetMatchingTk().getEta()) < minEta) return;
  // if ( abs(L1TkTaus.at(1).GetMatchingTk().getEta()) > maxEta) return;

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


//============================================================================
void CaloTk::FillDiTau_(vector<L1TkTauParticle> L1TkTaus1,
				vector<L1TkTauParticle> L1TkTaus2, 
				TH2D *hRate,
				TH2D *hEfficiency)
//============================================================================
{

  // Sanity check
  if( L1TkTaus1.size() < 1 ) return;
  if( L1TkTaus2.size() < 1 ) return;
  
  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1TkTaus1_mcMatched = GetMcMatchedL1TkTaus(L1TkTaus1);
  vector<L1TkTauParticle> L1TkTaus2_mcMatched = GetMcMatchedL1TkTaus(L1TkTaus2);

  // Fill rate 
  double ldgEt1 = L1TkTaus1.at(0).GetCaloTau().et();
  double ldgEt2 = L1TkTaus2.at(0).GetCaloTau().et();

  // Ensure that different calo objects are used
  unsigned int index1 = L1TkTaus1.at(0).GetCaloTau().index();
  unsigned int index2 = L1TkTaus2.at(0).GetCaloTau().index();
  if (index1==index2)
    {
      if (L1TkTaus2.size() < 2) return;
      index2 = L1TkTaus2.at(1).GetCaloTau().index();
    }

  // Make x-axis the ldgEt axis
  if (ldgEt1 > ldgEt2) FillRate_(hRate, ldgEt1, ldgEt2); 
  else FillRate_(hRate, ldgEt2, ldgEt1);

  
  // Get MC-matched trigger objects
  if (L1TkTaus1_mcMatched.size() < 1) return;
  if (L1TkTaus2_mcMatched.size() < 1) return;

  // Get MC-matched Et
  double ldgEt1_mcMatched = L1TkTaus1_mcMatched.at(0).GetCaloTau().et();
  double ldgEt2_mcMatched = L1TkTaus2_mcMatched.at(0).GetCaloTau().et();

  
  // Ensure that different calo objects are used
  index1 = L1TkTaus1_mcMatched.at(0).GetCaloTau().index();
  index2 = L1TkTaus2_mcMatched.at(0).GetCaloTau().index();
  if (index1==index2)
    {
      if (L1TkTaus2_mcMatched.size() < 2) return;
      index2 = L1TkTaus2_mcMatched.at(1).GetCaloTau().index();
    }

  // Check that all taus were found
  if(!bFoundAllTaus_) return;

  // Make x-axis the ldgEt axis
  if (ldgEt1_mcMatched > ldgEt2_mcMatched) histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt1_mcMatched, ldgEt2_mcMatched);
  else histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt2_mcMatched, ldgEt1_mcMatched);

  return;
}



//============================================================================
void CaloTk::FillRate_(TH1D *hRate,
			   const double ldgEt)
//============================================================================
{
  
  if (ldgEt < 0) return;
  hRate ->Fill( ldgEt );
  
  return;
}


//============================================================================
void CaloTk::FillRate_(TH2D *hRate,
			   const double ldgEt1,
			   const double ldgEt2)
//============================================================================
{
  
  if (ldgEt1 < 0) return;
  if (ldgEt2 < 0) return;

  hRate ->Fill( ldgEt1, ldgEt2 );
  
  return;
}


//============================================================================
void CaloTk::FillEfficiency_(TH1D *hEfficiency,
			     const double ldgEt)
//============================================================================
{
  
  histoTools_.FillAllBinsUpToValue_1D(hEfficiency, ldgEt);

  return;
}



//============================================================================
vector<GenParticle> CaloTk::GetHadronicGenTaus(vector<GenParticle> GenTaus,
						       double visEt,
						       double visEta)
//============================================================================
{
  // Sanity check
  vector<GenParticle> genP_hadGenTaus;
  if (GenTaus.size() < 1 ) return genP_hadGenTaus;
    
  // For-loop: GenTaus
  for (vector<GenParticle>::iterator tau = GenTaus.begin(); tau != GenTaus.end(); tau++)
    {

      // If no hadronic decay products found (pi+/-, pi0, K+/-, K0, K0L), skip this tau
      bool bIsHadronicTau = (tau->finalDaughtersCharged().size() > 0);

      if (0) tau->PrintFinalDaughters();
      
      // Apply acceptance cuts
      bool bPassVisEt  = ( tau->p4vis().Et() >= visEt);
      bool bPassVisEta = ( abs(tau->p4vis().Eta()) <= visEta );


      if (!(bPassVisEt)) continue;
      if (!(bPassVisEta)) continue;
      if (!bIsHadronicTau) continue;
      
      // Save this hadronic generator tau
      genP_hadGenTaus.push_back(*tau);
    }
  
  return genP_hadGenTaus;
}      


//============================================================================
void CaloTk::GetHadronicTauFinalDaughters(GenParticle p,
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
void CaloTk::FillTurnOn_Numerator_(vector<L1TkTauParticle> L1TkTaus, 
				   const double minEt,
				   TH1D *hTurnOn)
//============================================================================
{
  
  // For-loop: L1TkTaus
  for (vector<L1TkTauParticle>::iterator L1TkTau = L1TkTaus.begin(); L1TkTau != L1TkTaus.end(); L1TkTau++)
    {
      
      // Skip if trigger object is not MC matched
      if (!L1TkTau->HasMatchingGenParticle()) continue;	 
      
      // Skip if trigger object has eT < minEt
      if ( L1TkTau->GetCaloTau().et() < minEt) continue;
      
      // Get MC-match
      GenParticle p = L1TkTau->GetMatchingGenParticle();
      
      // Fill the turn-on
      hTurnOn->Fill( p.p4vis().Et() );

    } // For-loop: L1TkTaus

  return;
   
}


//============================================================================
void CaloTk::GetMatchingTrack(L1TkTauParticle &L1TkTau,
			      L1Tau myL1Tau,
			      vector<TTTrack> TTTracks)

//============================================================================
{

  // Initialise variables
  TTTrack matchTk;
  double matchTk_dR = 999.9;
    
  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
    {
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), myL1Tau.eta(), myL1Tau.phi());

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
  L1TkTau.SetCaloTau(myL1Tau.index());
  L1TkTau.SetMatchingTk(matchTk);
  L1TkTau.SetMatchTkDeltaRNew(matchTk_dR);
  if (0) L1TkTau.PrintProperties(false, false, true, false);
  
  return;
}


//============================================================================
void CaloTk::GetSigConeTracks(L1TkTauParticle &L1TkTau,
			      vector<TTTrack> TTTracks)
//============================================================================
{
  if (!L1TkTau.HasMatchingTk()) return; 
  vector<TTTrack> sigConeTks;
  
  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
    {

      // Skip the matching track
      if( tk->index() == L1TkTau.GetMatchingTk().index() ) continue;
      
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), L1TkTau.GetMatchingTk().getEta(), L1TkTau.GetMatchingTk().getPhi());

      // Only consider tracks within singal cone/annulus
      if (dR < L1TkTau.GetSigConeMin()) continue;
      if (dR > L1TkTau.GetSigConeMax()) continue;
      
      sigConeTks.push_back(*tk);
    }

  // Add the matching track
  sigConeTks.push_back(L1TkTau.GetMatchingTk());

  // Set the signal-cone tracks
  L1TkTau.SetSigConeTracks(sigConeTks);

  // Print properties
  if (0) L1TkTau.PrintProperties(false, false, true, false);
  
  return;
}


//============================================================================
void CaloTk::GetIsoConeTracks(L1TkTauParticle &L1TkTau,
				      vector<TTTrack> TTTracks)
//============================================================================
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


//============================================================================
vector<GenParticle> CaloTk::GetGenParticles(bool bPrintList)
//============================================================================
{
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
vector<GenParticle> CaloTk::GetGenParticles(int pdgId, bool isLastCopy)
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
GenParticle CaloTk::GetGenParticle(unsigned int Index)
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
void CaloTk::SetGenParticleMomsAndDaus(GenParticle &p)
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
void CaloTk::SetGenParticleFinalDaughters(GenParticle &p)
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
vector<TTTrack> CaloTk::GetTTTracks(const double minPt,
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
TTTrack CaloTk::GetTTTrack(unsigned int Index,
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
  if (DEBUG && matchTP_index >= 0) theTP.PrintProperties();

  return theTTTrack;
}



//============================================================================
vector<TrackingParticle> CaloTk::GetTrackingParticles(bool bPrintList)
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
TrackingParticle CaloTk::GetTrackingParticle(unsigned int Index)
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
  if (DEBUG) theTrackingParticle.PrintProperties();

  // Get the uniquely matched TTTrack
  // TTTrack theTrack;
  // if (ttTrackIndex >= 0) theTrack= GetTTTrack(ttTrackIndex, 4);
  // theTrackingParticle.SetTTTTrack(theTrack); // cannot use! cyclic problems 
  // if (DEBUG) theTrack.PrintProperties(); // cannot use! cyclic problems

  return theTrackingParticle;
}


//============================================================================
vector<L1Tau> CaloTk::GetL1Taus(bool bPrintList)
//============================================================================
{
  vector<L1Tau> theL1Taus;
  L1Tau theL1Tau;
  // For-loop: All L1Taus
  for (Size_t iCalo = 0; iCalo < L1Tau_Et->size(); iCalo++)
    {
      // cout << "iCalo = " << iCalo << endl;-
      theL1Tau = GetL1Tau(iCalo); 
      theL1Taus.push_back( theL1Tau );
    }
  
  if (bPrintList) PrintL1TauCollection(theL1Taus); 
  return theL1Taus;
}


//============================================================================
L1Tau CaloTk::GetL1Tau(unsigned int Index)
//============================================================================
{

  L1Tau theL1Tau(Index,
			 L1Tau_Et->at(Index),
			 L1Tau_Eta->at(Index),
			 L1Tau_Phi->at(Index),
			 L1Tau_IET->at(Index),
			 L1Tau_IEta->at(Index),
			 L1Tau_IPhi->at(Index),
			 L1Tau_Iso->at(Index),
			 L1Tau_Bx->at(Index),
			 0, //L1Tau_TowerIPhi->at(Index),
			 0, //L1Tau_TowerIEta->at(Index),
			 L1Tau_RawEt->at(Index),
			 L1Tau_IsoEt->at(Index),
			 L1Tau_NTT->at(Index),
			 L1Tau_HasEM->at(Index),
			 L1Tau_IsMerged->at(Index),
			 L1Tau_HwQual->at(Index)
			 );

  return theL1Tau;
}



//============================================================================
vector<L1Tau> CaloTk::GetL1Jets(bool bPrintList)
//============================================================================
{

  vector<L1Tau> theL1Jets;
  L1Tau theL1Jet;
  
  // For-loop: All L1 Jets
  for (Size_t i = 0; i < L1Jet_Et->size(); i++)
    {
      theL1Jet = GetL1Jet(i); 
      theL1Jets.push_back(theL1Jet);
    }
  
  if (bPrintList) PrintL1TauCollection(theL1Jets); 
  return theL1Jets;
}


//============================================================================
L1Tau CaloTk::GetL1Jet(unsigned int Index)
//============================================================================
{

  // int nJets = nJets->at(Index);
  double Et   = L1Jet_Et->at(Index);
  double Eta  = L1Jet_Eta->at(Index);
  double Phi  = L1Jet_Phi->at(Index);
  double Bx   = L1Jet_Bx->at(Index);
  double Type = -1.0;
  // L1Jet_IEt       
  // L1Jet_IEta      
  // L1Jet_IPhi      
  // L1Jet_Bx        
  // L1Jet_RawEt     
  // L1Jet_SeedEt    
  // L1Jet_TowerIEta 
  // L1Jet_TowerIPhi 
  // L1Jet_PUEt      
  // L1Jet_PUDonutEt0
  // L1Jet_PUDonutEt1
  // L1Jet_PUDonutEt2
  // L1Jet_PUDonutEt3

  // cout << "GeL1Jet returns theL1Jet(" << Index << ", " << Et << " , " << Eta << ", " << Phi << ", " << Bx << "," << Type << ")" << endl;
  L1Tau theL1Jet;//(Index, Et, Et, Eta, Phi, Bx, Type);

  return theL1Jet;
}


//============================================================================
void CaloTk::GetL1EGs(bool bPrintList)
//============================================================================
{

  // For-loop: All L1 EGs
  for (Size_t i = 0; i < L1EG_Et->size(); i++)
    {
      // cout << "L1EG " << i << "/" << L1EG_Et->size() << endl;
      if (bPrintList) GetL1EG(i);
    }
  
  return;
}


//============================================================================
void CaloTk::GetL1EG(unsigned int Index)
//============================================================================
{
  
  double Et      = L1EG_Et->at(Index);
  double Eta     = L1EG_Eta->at(Index);
  double Phi     = L1EG_Phi->at(Index);
  int IEt        = L1EG_IET->at(Index);
  int IEta       = L1EG_IEta->at(Index);
  int IPhi       = L1EG_IPhi->at(Index);
  int Iso        = L1EG_Iso->at(Index);
  int Bx         = L1EG_Bx->at(Index);
  // int TowerIPhi  = LEG_TowerIPhi->at(Index); 
  // int TowerIEta  = LEG_TowerIEta->at(Index); 
  int RawEt      = L1EG_RawEt->at(Index);
  int IsoEt      = L1EG_IsoEt->at(Index); 
  int FootprintEt= L1EG_FootprintEt->at(Index);
  int NTT        = L1EG_NTT->at(Index);
  int Shape      = L1EG_Shape->at(Index);
  int TowerHoE   = L1EG_TowerHoE->at(Index);
  
  cout << "nEGs        = " << L1EG_Et->size()
       << "\nEt          = " << Et
       << "\nEta         = " << Eta
       << "\nPhi         = " << Phi
       << "\nIEt         = " << IEt
       << "\nIEta        = " << IEta
       << "\nIPhi        = " << IPhi
       << "\nIso         = " << Iso
       << "\nBx          = " << Bx
    // << "\nTowerIPhi   = " << TowerIPhi
    // << "\nTowerIEta   = " << TowerIEta
       << "\nRawEt       = " << RawEt
       << "\nIsoEt       = " << IsoEt
       << "\nFootprintEt = " << FootprintEt
       << "\nNTT         = " << NTT
       << "\nShape       = " << Shape
       << "\nTowerHoE    = " << TowerHoE
       << "\n" << endl;
  
  return;
}


//============================================================================
void CaloTk::GetL1Sums(bool bPrintList)
//============================================================================
{

  // For-loop: All L1 Sums
  for (Size_t i = 0; i < L1Sum_Et->size(); i++)
    {
      // cout << "L1Sum(" << i << ")" << endl;
      if (bPrintList) GetL1Sum(i);
    }
  
  return;
}

//============================================================================
void CaloTk::GetL1Sum(unsigned int Index)
//============================================================================
{

  cout << "nSums        = " << L1Sum_Et->size()
       << "\nL1Sum_umType = " << L1Sum_Type->at(Index)
       << "\nL1Sum_umEt   = " << L1Sum_Et->at(Index)
       << "\nL1Sum_umPhi  = " << L1Sum_Phi->at(Index)
       << "\nL1Sum_umIET  = " << L1Sum_IET->at(Index)
       << "\nL1Sum_umIPhi = " << L1Sum_IPhi->at(Index)
       << "\nL1Sum_umBx   = " << L1Sum_Bx->at(Index)
       << "\n" << endl;
  
  return;
}

//============================================================================
void CaloTk::GetShrinkingConeSizes(double calo_et,
				   double sigCone_Constant,
				   double isoCone_Constant,
				   const double sigCone_dRCutoff,
				   double &sigCone_dRMin,
				   double &sigCone_dRMax,
				   double &isoCone_dRMin,
				   double &isoCone_dRMax)
//============================================================================
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


//============================================================================
void CaloTk::GetIsolationValues(L1TkTauParticle &L1TkTau)
//============================================================================
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
      deltaZ0 = abs(matchTk.getZ0() - isoConeTk.getZ0());
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


//============================================================================
void CaloTk::GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
					    vector<GenParticle> hadGenTaus)					    
//============================================================================
{

  //
  // Description: (to do)
  // Match the L1TkTau with a genParticle. At the moment try to match a genParticle (tau)
  // final decay products (pions, Kaons) with the matching track. If a match is found
  // the assign the tau (not the pion) to the L1TkTau as the matching genParticle.
  //

  // Sanity check
  if (hadGenTaus.size() < 1 ) return;

  // Initialise
  TTTrack matchTk  = L1TkTau.GetMatchingTk();
  double match_dR  = 9999.9;
  GenParticle match_GenParticle;

  // For-loop: All hadronic GenTaus
  for (vector<GenParticle>::iterator tau = hadGenTaus.begin(); tau != hadGenTaus.end(); tau++)
    {
      // If no hadronic decay products found (pi+/-, pi0, K+/-, K0, K0L), skip this tau
      if (tau->finalDaughtersCharged().size() < 1) continue;
      
      double deltaR = auxTools_.DeltaR( L1TkTau.GetCaloTau().eta(), L1TkTau.GetCaloTau().phi(), tau->eta(), tau->phi() );
      if (deltaR > mcMatching_dRMax) continue;
      if (deltaR < match_dR)
	{
	  match_dR = deltaR;
	  match_GenParticle = *tau;
	}
      
    }  // For-loop: All hadronic GenTaus

  // Save the matching
  L1TkTau.SetMatchingGenParticle(match_GenParticle);
  L1TkTau.SetMatchingGenParticleDeltaR(match_dR);

  // For debugging
  if (DEBUG)
    {
      L1TkTau.PrintProperties(false, true, false, false);
      // if (L1TkTau.HasMatchingGenParticle()) match_GenParticle.PrintProperties();
    }
  
  return;
}      


//============================================================================
vector<L1TkTauParticle> CaloTk::GetMcMatchedL1TkTaus(vector<L1TkTauParticle> L1TkTaus)
//============================================================================
{

  // Get all MC-matched trigger objects
  vector<L1TkTauParticle> matchedL1TkTaus;
  for (vector<L1TkTauParticle>::iterator tau = L1TkTaus.begin(); tau != L1TkTaus.end(); tau++)
    {
      if (!tau->HasMatchingGenParticle()) continue;
      matchedL1TkTaus.push_back(*tau);
    }
  
  return matchedL1TkTaus;
}


//============================================================================
void CaloTk::PrintGenParticleCollection(vector<GenParticle> collection)
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


//============================================================================
void CaloTk::PrintTrackingParticleCollection(vector<TrackingParticle> collection)
//============================================================================
{

  Table info("Index | Pt | Eta | Phi | PdgId | Q | x0 | y0 | z0 | d0 | d0-phi | NMatch | TTTrackIndex | TTClusters | TTStubs | TTTracks", "Text");

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
      row++;
    }

  info.Print();
  
  return;    
}


//============================================================================
void CaloTk::PrintTTTrackCollection(vector<TTTrack> collection)
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
void CaloTk::PrintL1TauCollection(vector<L1Tau> collection)
//============================================================================
{
  
  
  Table info("Index | Et | Eta | Phi | IET | IPhi | Iso | Bx | TowerIPhi | TowerIEta | RawEt | IsoEt | NTT | HasEM | IsMerged | HwQual | Type", "Text");

  // For-loop: All L1Taus
  int row=0;
  for (auto p = collection.begin(); p != collection.end(); p++)
  {
    // Construct table
    info.AddRowColumn(row, auxTools_.ToString( p->getIndex(), 1)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getEt(), 4)        );
    info.AddRowColumn(row, auxTools_.ToString( p->getEta(), 4)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getPhi(), 4)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getIET() , 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getIEta(), 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getIPhi(), 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getIso() , 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getBx()  , 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getTowerIPhi(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getTowerIEta(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getRawEt(), 3)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getIsoEt(), 3)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getNTT(), 3)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getHasEM(), 3)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getIsMerged(), 3)  );
    info.AddRowColumn(row, auxTools_.ToString( p->getHwQual(), 3)    );
    row++;
    
  }

  info.Print();
  return;
}


//============================================================================
void CaloTk::PrintL1TkTauParticleCollection(vector<L1TkTauParticle> collection)
//============================================================================
{
  
  Table info("Match-Cone | Sig-Cone | Iso-Cone | Calo-Et | Calo-Eta | Tk-Pt | Tk-Eta | Tk-dR | Gen-Pt | Gen-Eta | Gen-dR | Sig-Tks | Iso-Tks | VtxIso | RelIso", "Text");
  
  // For-loop: All L1TkTauParticles
  int row=0;
  for (vector<L1TkTauParticle>::iterator p = collection.begin(); p != collection.end(); p++)
  {
    // Construct table
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchConeMin(), 2) + " < dR < " + auxTools_.ToString( p->GetMatchConeMax(), 2) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetSigConeMin()  , 2) + " < dR < " + auxTools_.ToString( p->GetSigConeMax()  , 2) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetIsoConeMin()  , 2) + " < dR < " + auxTools_.ToString( p->GetIsoConeMax()  , 2) ); 
    info.AddRowColumn(0, auxTools_.ToString( p->GetCaloTau().et()       ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetCaloTau().eta()      ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingTk().getPt() ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingTk().getEta(),3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingTkDeltaR()   ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingGenParticle().pt() , 3) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingGenParticle().eta(), 3) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingGenParticleDeltaR(), 3) );  
    info.AddRowColumn(0, auxTools_.ToString( p->GetSigConeTTTracks().size() ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetIsoConeTTTracks().size() ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetVtxIsolation(), 3) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetRelIsolation(), 3) );
    row++;
  }

  info.Print();
  return;
}

#endif

#ifndef CaloTk_cxx
#define CaloTk_cxx

// User
#include "../Auxiliary/interface/constants.h"
#include "CaloTk.h"
//#include "../Framework/interface/TreeReaderReco.h"

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
  matchTk_nFitParams  =   4;         // TP:   5
  matchTk_minPt       =   5.00;      // TP:   5.0
  matchTk_minEta      =   0.0;       // TP:   0.0
  matchTk_maxEta      =  999.9;      // TP: 999.9  
  matchTk_maxChiSqRed =  200.0;      // TP: 200.0
  matchTk_minStubs    =   0;         // TP:   0
  matchTk_caloDeltaR  =   0.10;      // TP:   0.10

  // Signal cone tracks
  sigConeTks_Collection  = matchTk_Collection; // TP: "TTTracks" (not "TTPixelTracks")
  sigConeTks_nFitParams  = matchTk_nFitParams; // TP:   5
  sigConeTks_minPt       =   2.0;              // TP:   2.0
  sigConeTks_minEta      =   0.0;              // TP:   0.0
  sigConeTks_maxEta      = 999.9;              // TP: 999.9
  sigConeTks_maxChiSqRed = 200.0;              // TP: 200.0
  sigConeTks_minStubs    =   0;                // TP:   0

  // Isolation cone tracks
  isoConeTks_Collection  = matchTk_Collection;    // TP: "TTTracks" (not "TTPixelTracks")
  isoConeTks_nFitParams  = matchTk_nFitParams;    // TP:   5
  isoConeTks_minPt       =   2.0;                 // TP:   2.0
  isoConeTks_minEta      =   0.0;                 // TP:   0.0
  isoConeTks_maxEta      = 999.9;                 // TP: 999.9
  isoConeTks_maxChiSqRed = 200.0;                 // TP: 200.0
  isoConeTks_minStubs    =   0;                   // TP:   0
  
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
  isoCone_RelIsoWP = +0.20;         // TP: N/A
  isoCone_dRMin    = sigCone_dRMax; // TP: 0.4
  isoCone_dRMax    = +0.30;         // TP: 0.4
  diTau_deltaPOCAz = +1.00;         // TP: 1.0 cm

  // MC matching
  mcMatching_dRMax  = +0.1;        // TP: 0.05
  mcMatching_unique = true;

  // Eta Regions
  _eta_C = 0.8;
  _eta_F = 1.6;

  PrintSettings();

  return;
}


//============================================================================
void CaloTk::PrintSettings(void)
//============================================================================
{

  if (!DEBUG) return;

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

  settings.AddRowColumn(8, "Matching Tracks: DeltaR");
  settings.AddRowColumn(8, "<=");
  settings.AddRowColumn(8, auxTools_.ToString(matchTk_caloDeltaR) );
  settings.AddRowColumn(8, "0.10" );
  settings.AddRowColumn(8, "");

  settings.AddRowColumn(9, "Signal Cone Tks: Collection");
  settings.AddRowColumn(9, "==");
  settings.AddRowColumn(9, sigConeTks_Collection);
  settings.AddRowColumn(9, "TTTracks");
  settings.AddRowColumn(9, "");
  
  settings.AddRowColumn(10, "Signal Cone Tks: Fit Parameters");
  settings.AddRowColumn(10, "==");
  settings.AddRowColumn(10, auxTools_.ToString( sigConeTks_nFitParams) );
  settings.AddRowColumn(10, "5");
  settings.AddRowColumn(10, "");

  settings.AddRowColumn(11, "Signal Cone Tks: Pt");
  settings.AddRowColumn(11, ">=");
  settings.AddRowColumn(11, auxTools_.ToString( sigConeTks_minPt) );
  settings.AddRowColumn(11, "2" );
  settings.AddRowColumn(11, "GeV/c" );
  
  settings.AddRowColumn(12, "Signal Cone Tks: |Eta|");
  settings.AddRowColumn(12, ">=");
  settings.AddRowColumn(12, auxTools_.ToString( sigConeTks_minEta) );
  settings.AddRowColumn(12, "0.0" );
  settings.AddRowColumn(12, "" );

  settings.AddRowColumn(13, "Signal Cone Tks: |Eta|");
  settings.AddRowColumn(13, "<=");
  settings.AddRowColumn(13, auxTools_.ToString( sigConeTks_maxEta) );
  settings.AddRowColumn(13, "1e+03" );
  settings.AddRowColumn(13, "" );
  
  settings.AddRowColumn(14, "Signal Cone Tks: ChiSqRed");
  settings.AddRowColumn(14, "<=");
  settings.AddRowColumn(14, auxTools_.ToString( sigConeTks_maxChiSqRed) );
  settings.AddRowColumn(14, "200 (but on ChiSq, not ChiSqRed)");
  settings.AddRowColumn(14, "");

  settings.AddRowColumn(15, "Signal Cone Tks: Stubs");
  settings.AddRowColumn(15, ">=");
  settings.AddRowColumn(15, auxTools_.ToString( sigConeTks_minStubs) );
  settings.AddRowColumn(15, "" );
  settings.AddRowColumn(15, "" );

  settings.AddRowColumn(16, "Isolation Cone Tks: Collection");
  settings.AddRowColumn(16, "==");
  settings.AddRowColumn(16, isoConeTks_Collection);
  settings.AddRowColumn(16, "TTTracks");
  settings.AddRowColumn(16, "");
  
  settings.AddRowColumn(17, "Isolation Cone Tks: Fit Parameters");
  settings.AddRowColumn(17, "==");
  settings.AddRowColumn(17, auxTools_.ToString( isoConeTks_nFitParams) );
  settings.AddRowColumn(17, "5");
  settings.AddRowColumn(17, "");

  settings.AddRowColumn(18, "Isolation Cone Tks: Pt");
  settings.AddRowColumn(18, ">=");
  settings.AddRowColumn(18, auxTools_.ToString( isoConeTks_minPt) );
  settings.AddRowColumn(18, "2" );
  settings.AddRowColumn(18, "GeV/c" );
  
  settings.AddRowColumn(19, "Isolation Cone Tks: |Eta|");
  settings.AddRowColumn(19, ">=");
  settings.AddRowColumn(19, auxTools_.ToString( isoConeTks_minEta) );
  settings.AddRowColumn(19, "0.0" );
  settings.AddRowColumn(19, "" );

  settings.AddRowColumn(20, "Isolation Cone Tks: |Eta|");
  settings.AddRowColumn(20, "<=");
  settings.AddRowColumn(20, auxTools_.ToString( isoConeTks_maxEta) );
  settings.AddRowColumn(20, "1e+03" );
  settings.AddRowColumn(20, "" );

  settings.AddRowColumn(21, "Isolation Cone Tks: ChiSqRed");
  settings.AddRowColumn(21, "<=");
  settings.AddRowColumn(21, auxTools_.ToString( isoConeTks_maxChiSqRed) );
  settings.AddRowColumn(21, "200 (but on ChiSq, not ChiSqRed)");
  settings.AddRowColumn(21, "");

  settings.AddRowColumn(22, "Isolation Cone Tks: Stubs");
  settings.AddRowColumn(22, ">=");
  settings.AddRowColumn(22, auxTools_.ToString( isoConeTks_minStubs) );
  settings.AddRowColumn(22, "" );
  settings.AddRowColumn(22, "" );

  settings.AddRowColumn(23, "Signal Cone: Shrink Constant");
  settings.AddRowColumn(23, "==");
  settings.AddRowColumn(23, auxTools_.ToString(sigCone_Constant) );
  settings.AddRowColumn(23, "0" );
  settings.AddRowColumn(23, "GeV");

  settings.AddRowColumn(24, "Signal Cone: DeltaR");
  settings.AddRowColumn(24, ">=");
  settings.AddRowColumn(24, auxTools_.ToString(sigCone_dRMin) );
  settings.AddRowColumn(24, "0.0" );
  settings.AddRowColumn(24, "" );

  settings.AddRowColumn(25, "Signal Cone: DeltaR");
  settings.AddRowColumn(25, "<=");
  settings.AddRowColumn(25, auxTools_.ToString(sigCone_dRMax) );
  settings.AddRowColumn(25, "0.15" );
  settings.AddRowColumn(25, "" );
  
  settings.AddRowColumn(26, "Signal Cone:-3pr InvMass");
  settings.AddRowColumn(26, "<=");
  settings.AddRowColumn(26, auxTools_.ToString(sigCone_maxTkInvMass) );
  settings.AddRowColumn(26, "N/A" );
  settings.AddRowColumn(26, "GeV/c^{-2}");

  settings.AddRowColumn(27, "Signal Cone:-3pr maxTkDeltaPOCAz");
  settings.AddRowColumn(27, "<=");
  settings.AddRowColumn(27, auxTools_.ToString(sigCone_maxTkDeltaPOCAz) );
  settings.AddRowColumn(27, "N/A" );
  settings.AddRowColumn(27, "cm");

  // settings.AddRowColumn(19, "");

  settings.AddRowColumn(28, "Isolation Cone: Shrink Constant");
  settings.AddRowColumn(28, "==");
  settings.AddRowColumn(28, auxTools_.ToString(isoCone_Constant) );
  settings.AddRowColumn(28, "3.5");
  settings.AddRowColumn(28, "GeV");

  settings.AddRowColumn(29, "Isolation Cone: DeltaR");
  settings.AddRowColumn(29, ">=");
  settings.AddRowColumn(29, auxTools_.ToString(isoCone_dRMin) );
  settings.AddRowColumn(29, "0.15" );
  settings.AddRowColumn(29, "" );

  settings.AddRowColumn(30, "Isolation Cone: DeltaR");
  settings.AddRowColumn(30, "=<");
  settings.AddRowColumn(30, auxTools_.ToString(isoCone_dRMax) );
  settings.AddRowColumn(30, "0.30");
  settings.AddRowColumn(30, "");

  settings.AddRowColumn(31, "Isolation Cone: VtxIso" );
  settings.AddRowColumn(31, "<=" );
  settings.AddRowColumn(31, auxTools_.ToString(isoCone_VtxIsoWP) );
  settings.AddRowColumn(31, "1.0");
  settings.AddRowColumn(31, "cm");
  settings.AddRowColumn(31, "");

  settings.AddRowColumn(32, "Isolation Cone: RelIso" );
  settings.AddRowColumn(32, "<=" );
  settings.AddRowColumn(32, auxTools_.ToString(isoCone_RelIsoWP) );
  settings.AddRowColumn(32, "--");
  settings.AddRowColumn(32, "cm");
  settings.AddRowColumn(32, "");

  settings.AddRowColumn(33, "Di-Tau |Delta z0|");
  settings.AddRowColumn(33, "<");
  settings.AddRowColumn(33, auxTools_.ToString(diTau_deltaPOCAz) );
  settings.AddRowColumn(33, "1.0" );
  settings.AddRowColumn(33, "cm");

  settings.AddRowColumn(34, "MC-Matching DeltaR");
  settings.AddRowColumn(34, "<=");
  settings.AddRowColumn(34, auxTools_.ToString(mcMatching_dRMax) );
  settings.AddRowColumn(34, "0.05" );
  settings.AddRowColumn(34, "" );

  settings.AddRowColumn(35, "MC-Matching IsUnique");
  settings.AddRowColumn(35, "==");
  settings.AddRowColumn(35, auxTools_.ToString(mcMatching_unique) );
  settings.AddRowColumn(35, "1" );
  settings.AddRowColumn(35, "" );
  
  settings.AddRowColumn(36, "MC-Taus: Mom PdgId");
  settings.AddRowColumn(36, "==");
  settings.AddRowColumn(36, auxTools_.ToString(realTauMom));
  settings.AddRowColumn(36, "N/A" );
  settings.AddRowColumn(36, "" );

  settings.AddRowColumn(37, "MC-Taus: Number Expected");
  settings.AddRowColumn(37, ">=");
  settings.AddRowColumn(37, auxTools_.ToString(nMaxNumOfHTausPossible));
  settings.AddRowColumn(37, "N/A" );
  settings.AddRowColumn(37, "" );

  settings.AddRowColumn(38, "" );
  settings.Print();
  
  return;
}


//============================================================================
void CaloTk::Loop()
//============================================================================
{
  // Sanity check
  if (fChain == 0) return;
  
  const Long64_t nEntries   = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  if (DEBUG) cout << "=== CaloTk:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;

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
      if (DEBUG) std::cout << "Minimum Bias sample" << std::endl;
    }
  else
    {
      if (DEBUG) std::cout << "Not a Minimum Bias sample." << std::endl;
    }
  
  
  ////////////////////////////////////////////////
  // For-loop: Entries
  ////////////////////////////////////////////////
  for (int jentry = 0; jentry < nEntries; jentry++, nEvts++){
    
    if (DEBUG) cout << "\tEntry = (" << jentry << endl;
    
    // Init variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = fChain->GetEntry(jentry);
    nbytes += nb;
        
    //======================================================================================================
    
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
	// GenTausHadronic     = GetHadronicGenTaus(GenTaus, 00.0, 999.9); // tmp
	// GenTausTrigger      = GetHadronicGenTaus(GenTaus, 20.0, 2.4);   // tmp 
	GenTausHadronic     = GetHadronicGenTaus(GenTaus, 00.0, 1.479); // Calos restricted to Central Region
	GenTausTrigger      = GetHadronicGenTaus(GenTaus, 20.0, 1.479); // Calos restricted to Central Region
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
    if(DEBUG) cout << "\tGetting the Tracks and Track Particles Collections" << endl;
    vector<TrackingParticle> TPs = GetTrackingParticles(false);
    
    vector<TTTrack> matchTTTracks = GetTTTracks(matchTk_minPt, matchTk_minEta, matchTk_maxEta, matchTk_maxChiSqRed,
						matchTk_minStubs, matchTk_nFitParams, false);
    
    vector<TTTrack> sigTTTracks = GetTTTracks(sigConeTks_minPt , sigConeTks_minEta, sigConeTks_maxEta, sigConeTks_maxChiSqRed,
					      sigConeTks_minStubs, sigConeTks_nFitParams, false);
    
    vector<TTTrack> isoTTTracks = GetTTTracks(isoConeTks_minPt , isoConeTks_minEta, isoConeTks_maxEta, isoConeTks_maxChiSqRed,
					      isoConeTks_minStubs, isoConeTks_nFitParams, false);
    

    if (DEBUG)
      {
	cout << "\tPrinting all TTrack Collections" << endl;
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
    vector<L1TkTauParticle> L1TkTaus_RelIso;
    vector<L1TkTauParticle> L1TkTaus_Iso;
    // Jet Collections
    vector<L1Jet> L1Jets = GetL1Jets(false);
    vector<L1EG> L1EGs   = GetL1EGs(false);
    vector<L1Sum> L1Sums = GetL1Sums(false);

    // Ensure that all taus are found
    bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;

    // ======================================================================================
    // For-loop: GenTausHadronic
    // ======================================================================================
    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++)
      {
	
	// Get the charged pions from the tau decaying hadronically
	vector<unsigned short> chargedPionsIndices;
	vector<GenParticle>    chargedPions;

	//if (tau->HasMatchingGenParticle()) 
	GetHadronicTauChargedPions(tau->index(), chargedPionsIndices);

	// Ask for 3-prong or 5-prong decay
	if (chargedPionsIndices.size() >= 3)
	  {
	    double pions_dRMax = -1000.0;
	    double ldgPionPt  = -1000.0;
	    int ldgPionIndx   = -1;
	    GenParticle ldgPion;

	    // For-loop: All charged pions
	    for (unsigned int i=0; i< chargedPionsIndices.size(); i++)
	      {
		GenParticle pion = GetGenParticle(chargedPionsIndices.at(i));
		chargedPions.push_back(pion);
		
		if (0) PrintGenp(pion.index(), true);
		
		if (pion.pt() > ldgPionPt)
		  {
		    ldgPionIndx = chargedPionsIndices.at(i);
		    ldgPionPt   = pion.pt();
		    ldgPion     = pion;
		  }
		
	      }// For-loop: All charged pions
	    
	    for (vector<GenParticle>::iterator chPion = chargedPions.begin(); chPion != chargedPions.end(); chPion++)
	      {
		if (chPion->index() == ldgPionIndx) continue;
		double pions_dR = auxTools_.DeltaR(ldgPion.eta(), ldgPion.phi(), chPion->eta(), chPion->phi() );

		if ( pions_dR > pions_dRMax) pions_dRMax = pions_dR;
	      }
	    
	    // Fill histo
	    h_GenP_VisET_dRMaxLdgPion -> Fill(pions_dRMax, tau->p4vis().Et());

	  }// Ask for 3-prong or 5-prong decay

      }


    // ======================================================================================
    // For-loop: L1Taus
    // ======================================================================================
    unsigned int iCalo = -1;
    for (vector<L1Tau>::iterator calo = L1Taus.begin(); calo != L1Taus.end(); calo++)
      {
	iCalo += 1;
	if (0) calo->PrintProperties(iCalo==0);
	// std::cout << "NTT = " << calo->getNTT() << std::endl;

	// Fill histograms
	hL1CaloTau_Et        ->Fill( calo->getEt()    );
	hL1CaloTau_Eta       ->Fill( calo->getEta()   );
	hL1CaloTau_Phi       ->Fill( calo->getPhi()   );
	hL1CaloTau_IEt       ->Fill( calo->getIEt()   );
	hL1CaloTau_IEta      ->Fill( calo->getIEta()  );
	hL1CaloTau_IPhi      ->Fill( calo->getIPhi()  );
	hL1CaloTau_Iso       ->Fill( calo->getIso()   );
	hL1CaloTau_TowerIPhi ->Fill( calo->getIPhi()  );
	hL1CaloTau_TowerIEta ->Fill( calo->getIEta()  );
	hL1CaloTau_RawEt     ->Fill( calo->getRawEt() );
	hL1CaloTau_IsoEt     ->Fill( calo->getIsoEt() );
	hL1CaloTau_NTT       ->Fill( calo->getNTT()   );
	hL1CaloTau_HasEM     ->Fill( calo->getHasEM() );
	hL1CaloTau_IsMerged  ->Fill( calo->getIsMerged() );

	// Calculate the Et-dependent signal & isolation cone sizes
	GetShrinkingConeSizes(calo->et(), sigCone_Constant, isoCone_Constant, 
			      sigCone_cutoffDeltaR, sigCone_dRMin, sigCone_dRMax, 
			      isoCone_dRMin, isoCone_dRMax);
	
	// Construct the CaloTk candidate
	L1TkTauParticle	L1TkTauCandidate(0.0, matchTk_caloDeltaR, 
					 sigCone_dRMin, sigCone_dRMax, 
					 isoCone_dRMin, isoCone_dRMax);


	if ( matchTk_Collection.compare("TTTracks") != 0 )
	  {
	    cout << "=== ERROR: Invalid Track Collection Type \"" << matchTk_Collection << "\". EXIT" << endl;
	    exit(1);
	  }
	    
	GetMatchingTrack(L1TkTauCandidate, *calo, matchTTTracks);
	GetSigConeTracks(L1TkTauCandidate, sigTTTracks);
	GetIsoConeTracks(L1TkTauCandidate, isoTTTracks);
	GetIsolationValues(L1TkTauCandidate);
	//GetMatchingGenParticle(L1TkTauCandidate, GenTausHadronic); //marina
	GetMatchingGenParticle(L1TkTauCandidate, GenTausTrigger);
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

	// +Tk
	if (L1TkTau->HasMatchingTk() )
	  {
	    // std::cout << "\n=== L1TkTau->GetMatchingTkDeltaR() = " << L1TkTau->GetMatchingTkDeltaR() << std::endl;
	    L1TkTaus_Tk.push_back(*L1TkTau);

	    if (L1TkTau->GetVtxIsolation() > isoCone_VtxIsoWP) L1TkTaus_VtxIso.push_back(*L1TkTau); // +VtxIso
	    if (L1TkTau->GetRelIsolation() < isoCone_RelIsoWP) L1TkTaus_RelIso.push_back(*L1TkTau); // +RelIso
	    if (L1TkTau->GetVtxIsolation() > isoCone_VtxIsoWP && L1TkTau->GetRelIsolation() < isoCone_RelIsoWP) L1TkTaus_Iso.push_back(*L1TkTau); // +TkIso
	  }
      }// L1TkTauCandidates

    if (DEBUG)
      {
	PrintL1TauCollection(L1Taus);
	PrintL1TkTauParticleCollection(L1TkTauCandidates);
	PrintL1TkTauParticleCollection(L1TkTaus_Calo);
	PrintL1TkTauParticleCollection(L1TkTaus_Tk);
	PrintL1TkTauParticleCollection(L1TkTaus_VtxIso);
	PrintL1TkTauParticleCollection(L1TkTaus_RelIso);
	PrintL1TkTauParticleCollection(L1TkTaus_Iso);
      }
    
    ////////////////////////////////////////////////
    /// L1TkTaus_Calo Properties 
    ////////////////////////////////////////////////
    for (vector<L1TkTauParticle>::iterator tau = L1TkTaus_Calo.begin(); tau != L1TkTaus_Calo.end(); tau++)
      {
	// L1TkTau Resolution
	GenParticle p = tau->GetMatchingGenParticle();	
	hL1Tau_ResolutionCaloEt ->Fill( (tau->GetCaloTau().et() - p.p4vis().Et() )/p.p4vis().Et()  );
	hL1Tau_ResolutionCaloEta->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	hL1Tau_ResolutionCaloPhi->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() );
      }
    
    ////////////////////////////////////////////////
    /// L1Tktau Properties 
    ////////////////////////////////////////////////
    vector<L1TkTauParticle> myL1TkTaus = L1TkTaus_Tk;
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
	hL1TkTau_ResolutionCaloEt ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	hL1TkTau_ResolutionCaloEta->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	hL1TkTau_ResolutionCaloPhi->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() );
	
	double caloTau_eta =tau->GetCaloTau().eta();
	if ( IsWithinEtaRegion("Central", caloTau_eta) )
	  {
	    hL1TkTau_ResolutionCaloEt_C ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	    hL1TkTau_ResolutionCaloEta_C->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkTau_ResolutionCaloPhi_C->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	else if ( IsWithinEtaRegion("Intermediate", caloTau_eta) )
	  {
	    hL1TkTau_ResolutionCaloEt_I ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	    hL1TkTau_ResolutionCaloEta_I->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkTau_ResolutionCaloPhi_I->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	// currently no L1Taus in forward eta region
	else if ( IsWithinEtaRegion("Forward", caloTau_eta) )
	  {
	    hL1TkTau_ResolutionCaloEt_F ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	    hL1TkTau_ResolutionCaloEta_F->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkTau_ResolutionCaloPhi_F->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	else{                                                                                                                                                           
	  cout << "=== Tracking::Loop() - Unexpected Eta value of \"" << caloTau_eta << "\". EXIT" << endl;
	  exit(1);                                                                                                                                                      
	}                    
	
	// Apply isolation? 
	if (0) if (tau->GetVtxIsolation() <= isoCone_VtxIsoWP) continue; // Vertex Isolation
	if (0) if (tau->GetRelIsolation() >= isoCone_RelIsoWP) continue; // Relative Isolation

	// Matching Track Variables
	TTTrack matchTk   = tau->GetMatchingTk();
	double matchTk_dR = auxTools_.DeltaR(matchTk.getEta(), matchTk.getPhi(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi() ); // marina: can't we get the tau->GetMatchingTkDeltaR
	TLorentzVector caloTau_p4;
      	caloTau_p4.SetPtEtaPhiE(tau->GetCaloTau().et(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi(), tau->GetCaloTau().et() );
	hL1TkTau_MatchTk_DeltaR        ->Fill( matchTk_dR );
	hL1TkTau_MatchTk_PtRel         ->Fill( matchTk.p3().Perp(caloTau_p4.Vect()) );
	hL1TkTau_MatchTk_Pt            ->Fill( matchTk.getPt() );
	hL1TkTau_MatchTk_Eta           ->Fill( matchTk.getEta() );
	hL1TkTau_MatchTk_POCAz         ->Fill( matchTk.getZ0() );
	hL1TkTau_MatchTk_NStubs        ->Fill( matchTk.getNumOfStubs() );
	hL1TkTau_MatchTk_ChiSquared    ->Fill( matchTk.getChi2() );
	hL1TkTau_MatchTk_RedChiSquared ->Fill( matchTk.getChi2Red() );
	hL1TkTau_MatchTk_IsGenuine     ->Fill( matchTk.getIsGenuine() );
	hL1TkTau_MatchTk_IsUnknown     ->Fill( matchTk.getIsUnknown() );
	hL1TkTau_MatchTk_IsCombinatoric->Fill( matchTk.getIsCombinatoric() );
	hL1TkTau_MatchTk_PtMinusCaloEt ->Fill( matchTk.getPt() - tau->GetCaloTau().et() );

	// Signal/Isolation Cone Variables
	hL1TkTau_Rtau         ->Fill( tau->GetSigConeLdgTk().getPt() / tau->GetCaloTau().et() );
	hL1TkTau_CaloEt       ->Fill( tau->GetCaloTau().et() );
	hL1TkTau_CaloEta      ->Fill( tau->GetCaloTau().eta() );
	hL1TkTau_CaloEt       ->Fill( tau->GetCaloTau().getEt()    );
	hL1TkTau_CaloEta      ->Fill( tau->GetCaloTau().getEta()   );
	hL1TkTau_CaloPhi      ->Fill( tau->GetCaloTau().getPhi()   );
	hL1TkTau_CaloIEt      ->Fill( tau->GetCaloTau().getIEt()   );
	hL1TkTau_CaloIEta     ->Fill( tau->GetCaloTau().getIEta()  );
	hL1TkTau_CaloIPhi     ->Fill( tau->GetCaloTau().getIPhi()  );
	hL1TkTau_CaloIso      ->Fill( tau->GetCaloTau().getIso()   );
	hL1TkTau_CaloTowerIPhi->Fill( tau->GetCaloTau().getIPhi()  );
	hL1TkTau_CaloTowerIEta->Fill( tau->GetCaloTau().getIEta()  );
	hL1TkTau_CaloRawEt    ->Fill( tau->GetCaloTau().getRawEt() );
	hL1TkTau_CaloIsoEt    ->Fill( tau->GetCaloTau().getIsoEt() );
	hL1TkTau_CaloNTT      ->Fill( tau->GetCaloTau().getNTT()   );
	hL1TkTau_CaloHasEM    ->Fill( tau->GetCaloTau().getHasEM() );
	hL1TkTau_CaloIsMerged ->Fill( tau->GetCaloTau().getIsMerged() );
	hL1TkTau_CHF          ->Fill( tau->GetSigConeTTTracksP4().Et()/tau->GetCaloTau().et() );
	hL1TkTau_NHF          ->Fill( (tau->GetCaloTau().et() - tau->GetSigConeTTTracksP4().Et())/tau->GetCaloTau().et() );
	hL1TkTau_NHFAbs       ->Fill( abs( (tau->GetCaloTau().et() - tau->GetSigConeTTTracksP4().Et())/tau->GetCaloTau().et() ) );
	hL1TkTau_NSigTks      ->Fill( tau->GetSigConeTTTracks().size() );
	hL1TkTau_SigTksEt     ->Fill( tau->GetSigConeTTTracksP4().Et() );
	hL1TkTau_SigTksEta    ->Fill( tau->GetSigConeTTTracksP4().Eta() );
	hL1TkTau_NIsoTks      ->Fill( tau->GetIsoConeTTTracks().size() );
	if (tau->GetIsoConeTTTracks().size() > 0) 
	  {
	    hL1TkTau_IsoTksEt    ->Fill( tau->GetIsoConeTTTracksP4().Et() );
	    hL1TkTau_IsoTksEta   ->Fill( tau->GetIsoConeTTTracksP4().Eta() );
	  }
	if (tau->GetSigConeTTTracks().size() > 1) hL1TkTau_InvMass->Fill( tau->GetSigConeTTTracksP4().M() );
	hL1TkTau_InvMassIncl ->Fill( tau->GetSigConeTTTracksP4().M() ); 
	hL1TkTau_SigConeRMin ->Fill( tau->GetSigConeMin() );
	hL1TkTau_IsoConeRMin ->Fill( tau->GetIsoConeMin() );
	hL1TkTau_SigConeRMax ->Fill( tau->GetSigConeMax() );
	hL1TkTau_IsoConeRMax ->Fill( tau->GetIsoConeMax() );
	hL1TkTau_DeltaRGenP  ->Fill( tau->GetMatchingGenParticleDeltaR() );
	hL1TkTau_RelIso      ->Fill( tau->GetRelIsolation() );
	hL1TkTau_VtxIso      ->Fill( tau->GetVtxIsolation() );
	  

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
	    hL1TkTau_SigTks_DeltaR    ->Fill( sigTk_dR );
	    hL1TkTau_SigTks_NStubs    ->Fill( sigTk->getNumOfStubs() );
	    hL1TkTau_SigTks_ChiSquared->Fill( sigTk->getChi2() );
	    hL1TkTau_SigTks_RedChiSquared->Fill( sigTk->getChi2Red() );
	    hL1TkTau_SigTks_PtMinusCaloEt->Fill( sigTk->getPt() - tau->GetCaloTau().et() );

	    // Other variables
	    // sigTks_sumCharge += sigTk->getCharge(); //
	    
	  }// SigCone_TTTracks

	// Fill histos for other variables
	hL1TkTau_Charge->Fill( sigTks_sumCharge);

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
	    hL1TkTau_IsoTks_DeltaR    ->Fill( isoTk_dR );
	    hL1TkTau_IsoTks_NStubs    ->Fill( isoTk->getNumOfStubs() );
	    hL1TkTau_IsoTks_ChiSquared->Fill( isoTk->getChi2() );
	    hL1TkTau_IsoTks_RedChiSquared->Fill( isoTk->getChi2Red() );
	    hL1TkTau_IsoTks_PtMinusCaloEt->Fill( isoTk->getPt() - tau->GetCaloTau().et() );
  
	  }// IsoCone_TTTracks
      } // L1TkTaus_Tk

    ////////////////////////////////////////////////
    /// L1TkIsoTau Properties 
    ////////////////////////////////////////////////
    vector<L1TkTauParticle> myL1TkIsoTaus = L1TkTaus_Iso; // L1TkTaus_VtxIso;
    hL1TkIsoTau_Multiplicity ->Fill( myL1TkIsoTaus.size() );
    for (vector<L1TkTauParticle>::iterator tau = myL1TkIsoTaus.begin(); tau != myL1TkIsoTaus.end(); tau++)
      {
	
	if (DEBUG) tau->PrintProperties(true, true, true, true);

	// Variables
	TLorentzVector sigTks_p4 = tau->GetSigConeTTTracksP4();
	TLorentzVector isoTks_p4 = tau->GetIsoConeTTTracksP4();

	// Do not skip if using MinBias sample as no real taus exist!
	if (!tau->HasMatchingGenParticle() && (isMinBias == false) ) continue;
	
	// L1TkIsoTau Resolution
	GenParticle p = tau->GetMatchingGenParticle(); //fixme: more plots from MC info
      
	// Resolution
	hL1TkIsoTau_ResolutionCaloEt ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	hL1TkIsoTau_ResolutionCaloEta->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	hL1TkIsoTau_ResolutionCaloPhi->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() );
	
	double caloTau_eta =tau->GetCaloTau().eta();
	if ( IsWithinEtaRegion("Central", caloTau_eta) )
	  {
	    hL1TkIsoTau_ResolutionCaloEt_C ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	    hL1TkIsoTau_ResolutionCaloEta_C->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkIsoTau_ResolutionCaloPhi_C->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	else if ( IsWithinEtaRegion("Intermediate", caloTau_eta) )
	  {
	    hL1TkIsoTau_ResolutionCaloEt_I ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	    hL1TkIsoTau_ResolutionCaloEta_I->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkIsoTau_ResolutionCaloPhi_I->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	// currently no L1Taus in forward eta region
	else if ( IsWithinEtaRegion("Forward", caloTau_eta) )
	  {
	    hL1TkIsoTau_ResolutionCaloEt_F ->Fill( (tau->GetCaloTau().et() - p.p4vis().Pt() )/p.p4vis().Pt()  );
	    hL1TkIsoTau_ResolutionCaloEta_F->Fill( (tau->GetCaloTau().eta() - p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkIsoTau_ResolutionCaloPhi_F->Fill( (tau->GetCaloTau().phi() - p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	else{                                                                                                                                                           
	  cout << "=== Tracking::Loop() - Unexpected Eta value of \"" << caloTau_eta << "\". EXIT" << endl;
	  exit(1);                                                                                                                                                      
	}                    
	
	// Matching Track Variables
	TTTrack matchTk   = tau->GetMatchingTk();
	double matchTk_dR = auxTools_.DeltaR(matchTk.getEta(), matchTk.getPhi(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi() ); // marina: can't we get the tau->GetMatchingTkDeltaR
	TLorentzVector caloTau_p4;
      	caloTau_p4.SetPtEtaPhiE(tau->GetCaloTau().et(), tau->GetCaloTau().eta(), tau->GetCaloTau().phi(), tau->GetCaloTau().et() );
	hL1TkIsoTau_MatchTk_DeltaR        ->Fill( matchTk_dR );
	hL1TkIsoTau_MatchTk_PtRel         ->Fill( matchTk.p3().Perp(caloTau_p4.Vect()) );
	hL1TkIsoTau_MatchTk_Pt            ->Fill( matchTk.getPt() );
	hL1TkIsoTau_MatchTk_Eta           ->Fill( matchTk.getEta() );
	hL1TkIsoTau_MatchTk_POCAz         ->Fill( matchTk.getZ0() );
	hL1TkIsoTau_MatchTk_NStubs        ->Fill( matchTk.getNumOfStubs() );
	hL1TkIsoTau_MatchTk_ChiSquared    ->Fill( matchTk.getChi2() );
	hL1TkIsoTau_MatchTk_RedChiSquared ->Fill( matchTk.getChi2Red() );
	hL1TkIsoTau_MatchTk_IsGenuine     ->Fill( matchTk.getIsGenuine() );
	hL1TkIsoTau_MatchTk_IsUnknown     ->Fill( matchTk.getIsUnknown() );
	hL1TkIsoTau_MatchTk_IsCombinatoric->Fill( matchTk.getIsCombinatoric() );
	hL1TkIsoTau_MatchTk_PtMinusCaloEt ->Fill( matchTk.getPt() - tau->GetCaloTau().et() );

	// Signal/Isolation Cone Variables
	hL1TkIsoTau_Rtau         ->Fill( tau->GetSigConeLdgTk().getPt() / tau->GetCaloTau().et() );
	hL1TkIsoTau_CaloEt       ->Fill( tau->GetCaloTau().et() );
	hL1TkIsoTau_CaloEta      ->Fill( tau->GetCaloTau().eta() );
	hL1TkIsoTau_CaloEt       ->Fill( tau->GetCaloTau().getEt()    );
	hL1TkIsoTau_CaloEta      ->Fill( tau->GetCaloTau().getEta()   );
	hL1TkIsoTau_CaloPhi      ->Fill( tau->GetCaloTau().getPhi()   );
	hL1TkIsoTau_CaloIEt      ->Fill( tau->GetCaloTau().getIEt()   );
	hL1TkIsoTau_CaloIEta     ->Fill( tau->GetCaloTau().getIEta()  );
	hL1TkIsoTau_CaloIPhi     ->Fill( tau->GetCaloTau().getIPhi()  );
	hL1TkIsoTau_CaloIso      ->Fill( tau->GetCaloTau().getIso()   );
	hL1TkIsoTau_CaloTowerIPhi->Fill( tau->GetCaloTau().getIPhi()  );
	hL1TkIsoTau_CaloTowerIEta->Fill( tau->GetCaloTau().getIEta()  );
	hL1TkIsoTau_CaloRawEt    ->Fill( tau->GetCaloTau().getRawEt() );
	hL1TkIsoTau_CaloIsoEt    ->Fill( tau->GetCaloTau().getIsoEt() );
	hL1TkIsoTau_CaloNTT      ->Fill( tau->GetCaloTau().getNTT()   );
	hL1TkIsoTau_CaloHasEM    ->Fill( tau->GetCaloTau().getHasEM() );
	hL1TkIsoTau_CaloIsMerged ->Fill( tau->GetCaloTau().getIsMerged() );
	hL1TkIsoTau_CHF          ->Fill( tau->GetSigConeTTTracksP4().Et()/tau->GetCaloTau().et() );
	hL1TkIsoTau_NHF          ->Fill( (tau->GetCaloTau().et() - tau->GetSigConeTTTracksP4().Et())/tau->GetCaloTau().et() );
	hL1TkIsoTau_NHFAbs       ->Fill( abs( (tau->GetCaloTau().et() - tau->GetSigConeTTTracksP4().Et())/tau->GetCaloTau().et() ) );
	hL1TkIsoTau_NSigTks      ->Fill( tau->GetSigConeTTTracks().size() );
	hL1TkIsoTau_SigTksEt     ->Fill( tau->GetSigConeTTTracksP4().Et() );
	hL1TkIsoTau_SigTksEta    ->Fill( tau->GetSigConeTTTracksP4().Eta() );
	hL1TkIsoTau_NIsoTks      ->Fill( tau->GetIsoConeTTTracks().size() );
	if (tau->GetIsoConeTTTracks().size() > 0) 
	  {
	    hL1TkIsoTau_IsoTksEt    ->Fill( tau->GetIsoConeTTTracksP4().Et() );
	    hL1TkIsoTau_IsoTksEta   ->Fill( tau->GetIsoConeTTTracksP4().Eta() );
	  }
	if (tau->GetSigConeTTTracks().size() > 1) hL1TkIsoTau_InvMass->Fill( tau->GetSigConeTTTracksP4().M() );
	hL1TkIsoTau_InvMassIncl ->Fill( tau->GetSigConeTTTracksP4().M() ); 
	hL1TkIsoTau_SigConeRMin ->Fill( tau->GetSigConeMin() );
	hL1TkIsoTau_IsoConeRMin ->Fill( tau->GetIsoConeMin() );
	hL1TkIsoTau_SigConeRMax ->Fill( tau->GetSigConeMax() );
	hL1TkIsoTau_IsoConeRMax ->Fill( tau->GetIsoConeMax() );
	hL1TkIsoTau_DeltaRGenP  ->Fill( tau->GetMatchingGenParticleDeltaR() );
	hL1TkIsoTau_RelIso      ->Fill( tau->GetRelIsolation() );
	hL1TkIsoTau_VtxIso      ->Fill( tau->GetVtxIsolation() );

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
	    hL1TkIsoTau_SigTks_Pt        ->Fill( sigTk->getPt()  );
	    hL1TkIsoTau_SigTks_PtRel     ->Fill( sigTk_PtRel );
	    hL1TkIsoTau_SigTks_Eta       ->Fill( sigTk->getEta() );
	    hL1TkIsoTau_SigTks_POCAz     ->Fill( sigTk->getZ0()  );
	    if (sigTks.size() > 1)
	      {
		hL1TkIsoTau_SigTks_DeltaPOCAz->Fill( abs( sigTk->getZ0() - matchTk.getZ0() ) );
	      }
	    hL1TkIsoTau_SigTks_DeltaR    ->Fill( sigTk_dR );
	    hL1TkIsoTau_SigTks_NStubs    ->Fill( sigTk->getNumOfStubs() );
	    hL1TkIsoTau_SigTks_ChiSquared->Fill( sigTk->getChi2() );
	    hL1TkIsoTau_SigTks_RedChiSquared->Fill( sigTk->getChi2Red() );
	    hL1TkIsoTau_SigTks_PtMinusCaloEt->Fill( sigTk->getPt() - tau->GetCaloTau().et() );

	    // Other variables
	    // sigTks_sumCharge += sigTk->getCharge(); //
	    
	  }// SigCone_TTTracks

	// Fill histos for other variables
	hL1TkIsoTau_Charge->Fill( sigTks_sumCharge);

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
	    hL1TkIsoTau_IsoTks_Pt        ->Fill( isoTk->getPt()  );
	    hL1TkIsoTau_IsoTks_PtRel     ->Fill( isoTk_PtRel );
	    hL1TkIsoTau_IsoTks_Eta       ->Fill( isoTk->getEta() );
	    hL1TkIsoTau_IsoTks_POCAz     ->Fill( isoTk->getZ0()  );
	    hL1TkIsoTau_IsoTks_DeltaPOCAz->Fill( abs( isoTk->getZ0() - matchTk.getZ0() ) );
	    hL1TkIsoTau_IsoTks_DeltaR    ->Fill( isoTk_dR );
	    hL1TkIsoTau_IsoTks_NStubs    ->Fill( isoTk->getNumOfStubs() );
	    hL1TkIsoTau_IsoTks_ChiSquared->Fill( isoTk->getChi2() );
	    hL1TkIsoTau_IsoTks_RedChiSquared->Fill( isoTk->getChi2Red() );
	    hL1TkIsoTau_IsoTks_PtMinusCaloEt->Fill( isoTk->getPt() - tau->GetCaloTau().et() );
  
	  }// IsoCone_TTTracks
      } // L1TkTaus_VtxIso

    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////
    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++) hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() );
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 25.0, hCalo_TurnOn25  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 25.0, hTk_TurnOn25    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn25);
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 25.0, hRelIso_TurnOn25);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 25.0, hIso_TurnOn25);

    FillTurnOn_Numerator_(L1TkTaus_Calo   , 50.0, hCalo_TurnOn50  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 50.0, hTk_TurnOn50    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn50);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 50.0, hIso_TurnOn50);

    FillTurnOn_Numerator_(L1TkTaus_Calo   , 66.0, hCalo_TurnOn_SingleTau50KHz  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 65.0, hTk_TurnOn_SingleTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn_SingleTau50KHz);    
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 50.0, hRelIso_TurnOn_SingleTau50KHz);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 50.0, hIso_TurnOn_SingleTau50KHz);

    FillTurnOn_Numerator_(L1TkTaus_Calo   , 42.0, hCalo_TurnOn_DiTau50KHz  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 40.0, hTk_TurnOn_DiTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn_DiTau50KHz);
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 25.0, hRelIso_TurnOn_DiTau50KHz);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 25.0, hIso_TurnOn_DiTau50KHz);
    
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

    FillSingleTau_(L1TkTaus_RelIso, hRelIso_Rate  , hRelIso_Eff);
    FillSingleTau_(L1TkTaus_RelIso, hRelIso_Rate_C, hRelIso_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_RelIso, hRelIso_Rate_I, hRelIso_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_RelIso, hRelIso_Rate_F, hRelIso_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkTaus_Iso, hIso_Rate  , hIso_Eff);
    FillSingleTau_(L1TkTaus_Iso, hIso_Rate_C, hIso_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_Iso, hIso_Rate_I, hIso_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_Iso, hIso_Rate_F, hIso_Eff_F, 1.6, 3.0); // 2.5 is max
    
    ////////////////////////////////////////////////
    // DiTau
    ////////////////////////////////////////////////
    FillDiTau_(L1TkTaus_Calo, L1TkTaus_Tk    , hDiTau_Rate_Calo_Tk    , hDiTau_Eff_Calo_Tk     );
    FillDiTau_(L1TkTaus_Calo, L1TkTaus_VtxIso, hDiTau_Rate_Calo_VtxIso, hDiTau_Eff_Calo_VtxIso );
    FillDiTau_(L1TkTaus_Calo, L1TkTaus_RelIso, hDiTau_Rate_Calo_RelIso, hDiTau_Eff_Calo_RelIso );
    FillDiTau_(L1TkTaus_Calo, L1TkTaus_Iso   , hDiTau_Rate_Calo_Iso   , hDiTau_Eff_Calo_Iso    );

    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_VtxIso, hDiTau_Rate_Tk_VtxIso  , hDiTau_Eff_Tk_VtxIso   );
    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_RelIso, hDiTau_Rate_Tk_RelIso  , hDiTau_Eff_Tk_RelIso   );
    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_Iso   , hDiTau_Rate_Tk_Iso     , hDiTau_Eff_Tk_Iso      );

    ////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////
    ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_Tk);
    ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_VtxIso); 
    ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_RelIso); //fixme. correct?
    ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_Iso); //fixme. correct?

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

    FillDiTau_(L1TkTaus_RelIso, hDiTau_Rate_RelIso  , hDiTau_Eff_RelIso);
    FillDiTau_(L1TkTaus_RelIso, hDiTau_Rate_RelIso_C, hDiTau_Eff_RelIso_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_RelIso, hDiTau_Rate_RelIso_I, hDiTau_Eff_RelIso_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_RelIso, hDiTau_Rate_RelIso_F, hDiTau_Eff_RelIso_F, 1.6, 3.0);

    FillDiTau_(L1TkTaus_Iso, hDiTau_Rate_Iso  , hDiTau_Eff_Iso);
    FillDiTau_(L1TkTaus_Iso, hDiTau_Rate_Iso_C, hDiTau_Eff_Iso_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_Iso, hDiTau_Rate_Iso_I, hDiTau_Eff_Iso_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_Iso, hDiTau_Rate_Iso_F, hDiTau_Eff_Iso_F, 1.6, 3.0);

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

  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate_F, nEntries);

  histoTools_.ConvertToRateHisto_1D(hIso_Rate  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hIso_Rate_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hIso_Rate_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hIso_Rate_F, nEntries);

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

  FinaliseEffHisto_( hRelIso_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIso_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIso_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIso_Eff_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hIso_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hIso_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hIso_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hIso_Eff_F, nEvtsWithMaxHTaus);

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
  
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso_F, nEntries);

  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso  , nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso_F, nEntries);
  
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

  FinaliseEffHisto_( hDiTau_Eff_RelIso  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIso_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIso_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIso_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hDiTau_Eff_Iso  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Iso_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Iso_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_Iso_F, nEvtsWithMaxHTaus);

  // DiTau (Calo-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_Tk    , nEntries);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_VtxIso, nEntries);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_RelIso, nEntries);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Calo_Iso   , nEntries);

  FinaliseEffHisto_( hDiTau_Eff_Calo_Tk    , nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Calo_VtxIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Calo_RelIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Calo_Iso   , nEvtsWithMaxHTaus);;

  // DiTau (Tk-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_VtxIso, nEntries);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_RelIso, nEntries);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_Iso   , nEntries);

  FinaliseEffHisto_( hDiTau_Eff_Tk_VtxIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_RelIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_Iso   , nEvtsWithMaxHTaus);;

  // Turn-Ons: fixme. iro. alex
  // TEfficiency *pEff = 0;
  // pEff = new TEfficiency(*hCalo_TurnOn50_passed, *hMcHadronicTau_VisEt);
  // hCalo_TurnOn50 = (TH1D*) pEff->Clone();
  histoTools_.DivideHistos_1D(hCalo_TurnOn50   , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn50     , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn50 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hIso_TurnOn50    , hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hCalo_TurnOn25   , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn25     , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn25 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hIso_TurnOn25    , hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hCalo_TurnOn_SingleTau50KHz  , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn_SingleTau50KHz    , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn_SingleTau50KHz, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn_SingleTau50KHz, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hIso_TurnOn_SingleTau50KHz   , hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hCalo_TurnOn_DiTau50KHz  , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn_DiTau50KHz    , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn_DiTau50KHz, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn_DiTau50KHz, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hIso_TurnOn_DiTau50KHz   , hMcHadronicTau_VisEt);


  ////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////
  WriteHistos_();
  if (DEBUG) auxTools_.StopwatchStop(5, "minutes", "Total Time");
  
}


//============================================================================
void CaloTk::BookHistos_(void)
//============================================================================
{
  
  // Binning
  const unsigned int nEt = 300;
  const float minEt = 0.0;
  const float maxEt = +300.0;

  const unsigned int nPt = 200;
  const float minPt = 0.0;
  const float maxPt = +200.0;

  const unsigned int nPtR = 400;
  const float minPtR = 0.0;
  const float maxPtR = +20.0;

  const unsigned int nIEta = 120;
  const float minIEta = -60.0;
  const float maxIEta = +60.0;

  const unsigned int nIPhi = 150;
  const float minIPhi =   0.0;
  const float maxIPhi = 150.0;

  const unsigned int nEta = 60;
  const float minEta = -3.0;
  const float maxEta = +3.0;

  const unsigned int nPhi = 64;
  const float minPhi = -3.2;
  const float maxPhi = +3.2;

  const unsigned int nN = 15;
  const float minN =  -0.5;
  const float maxN = +14.5;

  const unsigned int nBool = 2;
  const float minBool = -0.5;
  const float maxBool = +1.5;

  const unsigned int nM = 100;
  const float minM =  0.0;
  const float maxM = 10.0;

  const unsigned int nR = 100;
  const float minR = 0.0;
  const float maxR = 1.0;

  const unsigned int nDR = 100;
  const float minDR = 0.0;
  const float maxDR = 1.0;

  const unsigned int nChi = 500;
  const float minChi =   0.0;
  const float maxChi = 500.0;

  const unsigned int nRChi = 200;
  const float minRChi =   0.0;
  const float maxRChi = 200.0;

  const unsigned int nZ0 = 600;
  const float minZ0 = -30.0;
  const float maxZ0 = +30.0;

  const unsigned int nRIso = 1000;
  const float minRIso =   0.0;
  const float maxRIso = +50.0;

  const unsigned int nVIso = 500;
  const float minVIso =   0.0;
  const float maxVIso = +50.0;
  
  const char* tEt   = ";E_{T} (GeV); Entries / %.0f GeV";
  const char* tPt   = ";p_{T} (GeV/c);Entries / %.0f GeV/c";
  const char* tPtR  = ";p_{T}^{rel} (GeV/c);Entries / %.1f GeV/c";
  const char* tEta  = ";#eta;Entries / %.2f";
  const char* tPhi  = ";#phi (rads);Entries / %.2f (rads)";
  const char* tM    = ";M (GeV/c^{2});Entries / %.2f (GeV/c^{2});";
  const char* tN    = ";multiplicity;Entries / %.0f";
  const char* tQ    = ";charge (e);Entries / %.0f e";
  const char* tR    = ";R;Entries / %.2f";
  const char* tDR   = ";#DeltaR;Entries / %.2f";
  const char* tZ0   = ";z_{0} (cm);Entries / %.2f cm";
  const char* tDZ0  = ";#Deltaz_{0} (cm);Entries / %.2f cm";
  const char* tRIso = ";relative isolation;Entries / %.2f";
  const char* tVIso = ";min(z_{0}^{m} - z_{0}^{iso} (cm);Entries / %.2f cm";
  const char* tChi  = ";#chi^{2};Entries / %.2f";
  const char* tRChi = ";#chi^{2}_{#nu};Entries / %.2f";
  const char* tRate = ";#E_{T} (GeV);Rate (kHz) / %.0f GeV";

  // GenParticles Histograms
  histoTools_.BookHisto_2D(h_GenP_VisET_dRMaxLdgPion, "GenP_VisET_dRMaxLdgPion", ";#DeltaR_{max}(#pi_{ldg}^{#pm},#pi^{#pm});E_{T}^{vis}",  50,  0.0, +0.25, 100, 0.0, +200.0);

  // Counters
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 2, 0.0, +2.0);

  // L1CaloTaus
  histoTools_.BookHisto_1D(hL1CaloTau_Et       , "L1CaloTau_Et"       , tEt ,  nEt ,  minEt , maxEt   );
  histoTools_.BookHisto_1D(hL1CaloTau_Eta      , "L1CaloTau_Eta"      , tEta,  nEta,  minEta, maxEta  );
  histoTools_.BookHisto_1D(hL1CaloTau_Phi      , "L1CaloTau_Phi"      , tPhi,  nPhi,  minPhi, maxPhi  );
  histoTools_.BookHisto_1D(hL1CaloTau_IEt      , "L1CaloTau_IEt"      , tEt ,  nEt ,  minEt , maxEt   );
  histoTools_.BookHisto_1D(hL1CaloTau_IEta     , "L1CaloTau_IEta"     , tEta, nIEta, minIEta, maxIEta );
  histoTools_.BookHisto_1D(hL1CaloTau_IPhi     , "L1CaloTau_IPhi"     , tPhi, nIPhi, minIPhi, maxIPhi );
  histoTools_.BookHisto_1D(hL1CaloTau_Iso      , "L1CaloTau_Iso"      ,   "", nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hL1CaloTau_TowerIEta, "L1CaloTau_TowerIEta", tEta, nIEta, minIEta, maxIEta );
  histoTools_.BookHisto_1D(hL1CaloTau_TowerIPhi, "L1CaloTau_TowerIPhi", tPhi, nIPhi, minIPhi, maxIPhi );
  histoTools_.BookHisto_1D(hL1CaloTau_RawEt    , "L1CaloTau_RawEt"    , tEt , nEt*3, minEt  , maxEt*3 );
  histoTools_.BookHisto_1D(hL1CaloTau_IsoEt    , "L1CaloTau_IsoEt"    , tEt , nEt  , minEt  , maxEt   );
  histoTools_.BookHisto_1D(hL1CaloTau_NTT      , "L1CaloTau_NTT"      , tN  ,nN*100, minN   , maxN*100);
  histoTools_.BookHisto_1D(hL1CaloTau_HasEM    , "L1CaloTau_HasEM"    , ""  , nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hL1CaloTau_IsMerged , "L1CaloTau_IsMerged" , ""  , nBool, minBool, maxBool );

  // L1TkTaus
  histoTools_.BookHisto_1D(hL1TkTau_Multiplicity , "L1TkTau_Multiplicity" , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_CaloEt       , "L1TkTau_CaloEt"       , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkTau_CaloEta      , "L1TkTau_CaloEta"      , tEta , nEta , minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkTau_CaloPhi      , "L1TkTau_CaloPhi"      , tPhi , nPhi , minPhi , maxPhi );
  histoTools_.BookHisto_1D(hL1TkTau_CaloIEt      , "L1TkTau_CaloIEt"      , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkTau_CaloIEta     , "L1TkTau_CaloIEta"     , tEta , nIEta, minIEta, maxIEta);
  histoTools_.BookHisto_1D(hL1TkTau_CaloIPhi     , "L1TkTau_CaloIPhi"     , tPhi , nIPhi, minIPhi, maxIPhi);
  histoTools_.BookHisto_1D(hL1TkTau_CaloIso      , "L1TkTau_CaloIso"      ,   "" , nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_CaloTowerIEta, "L1TkTau_CaloTowerIEta", tEta , nIEta, minIEta, maxIEta);
  histoTools_.BookHisto_1D(hL1TkTau_CaloTowerIPhi, "L1TkTau_CaloTowerIPhi", tPhi , nIPhi, minIPhi, maxIPhi);
  histoTools_.BookHisto_1D(hL1TkTau_CaloRawEt    , "L1TkTau_CaloRawEt"    , tEt  , nEt*3, minEt  , maxEt*3);
  histoTools_.BookHisto_1D(hL1TkTau_CaloIsoEt    , "L1TkTau_CaloIsoEt"    , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkTau_CaloNTT      , "L1TkTau_CaloNTT"      , tN   ,nN*100, minN   , maxN*100);
  histoTools_.BookHisto_1D(hL1TkTau_CaloHasEM    , "L1TkTau_CaloHasEM"    , ""   , nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_CaloIsMerged , "L1TkTau_CaloIsMerged" , ""   , nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_NSigTks      , "L1TkTau_NSigTks"      , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTksEt     , "L1TkTau_SigTksEt"     , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTksEta    , "L1TkTau_SigTksEta"    , tEta , nEta , minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkTau_NIsoTks      , "L1TkTau_NIsoTks"      , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTksEt     , "L1TkTau_IsoTksEt"     , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTksEta    , "L1TkTau_IsoTksEta"    , tEta , nEta , minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkTau_InvMass      , "L1TkTau_InvMass"      , tM   , nM   , minM   , maxM   );
  histoTools_.BookHisto_1D(hL1TkTau_InvMassIncl  , "L1TkTau_InvMassIncl"  , tM   , nM   , minM   , maxM   );
  histoTools_.BookHisto_1D(hL1TkTau_SigConeRMin  , "L1TkTau_SigConeRMin"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_SigConeRMax  , "L1TkTau_SigConeRMax"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_IsoConeRMin  , "L1TkTau_IsoConeRMin"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_IsoConeRMax  , "L1TkTau_IsoConeRMax"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_Charge       , "L1TkTau_Charge"       , tQ   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_RelIso       , "L1TkTau_RelIso"       , tRIso, nRIso, minRIso, maxRIso);
  histoTools_.BookHisto_1D(hL1TkTau_VtxIso       , "L1TkTau_VtxIso"       , tVIso, nVIso, minVIso, maxVIso);
  histoTools_.BookHisto_1D(hL1TkTau_DeltaRGenP   , "L1TkTau_DeltaRGenP"   , tDR  , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_Rtau         , "L1TkTau_Rtau"  , ";Entries / 0.1f; R_{#tau}", 100,  0.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_CHF          , "L1TkTau_CHF"   , ";Entries / 0.1f; CHF"     , 500,  0.0,  +50.0);
  histoTools_.BookHisto_1D(hL1TkTau_NHF          , "L1TkTau_NHF"   , ";Entries / 0.1f; NHF"     , 200, -5.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_NHFAbs       , "L1TkTau_NHFAbs", ";Entries / 0.1f; |NHF|"   , 100,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Pt           , "L1TkTau_SigTks_Pt"           , tPt  ,  nPt  , minPt  , maxPt   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Eta          , "L1TkTau_SigTks_Eta"          , tEta ,  nEta , minEta , maxEta  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_POCAz        , "L1TkTau_SigTks_POCAz"        , tZ0  ,  nZ0  , minZ0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_DeltaPOCAz   , "L1TkTau_SigTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,     0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_PtRel        , "L1TkTau_SigTks_PtRel"        , tPtR ,  nPtR , minPtR , maxPtR  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_DeltaR       , "L1TkTau_SigTks_DeltaR"       , tDR  ,  nR   , minR   , maxR    );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_NStubs       , "L1TkTau_SigTks_NStubs"       , tN   ,  nN   , minN   , maxN    );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_ChiSquared   , "L1TkTau_SigTks_ChiSquared"   , tChi ,  nChi , minChi , maxChi  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_RedChiSquared, "L1TkTau_SigTks_RedChiSquared", tRChi,  nRChi, minRChi, maxRChi );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_PtMinusCaloEt, "L1TkTau_SigTks_PtMinusCaloEt",    "",   300 , -250.0 , +250.0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Pt           , "L1TkTau_IsoTks_Pt"           , tPt  ,  nPt  ,  minPt ,  maxPt  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Eta          , "L1TkTau_IsoTks_Eta"          , tEta ,  nEta ,  minEta,  maxEta );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_POCAz        , "L1TkTau_IsoTks_POCAz"        , tZ0  ,  nZ0  ,  minZ0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_DeltaPOCAz   , "L1TkTau_IsoTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,      0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_PtRel        , "L1TkTau_IsoTks_PtRel"        , tPtR ,  nPtR ,  minPtR,  maxPtR );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_DeltaR       , "L1TkTau_IsoTks_DeltaR"       , tDR  ,  nR   ,  minR  ,  maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_NStubs       , "L1TkTau_IsoTks_NStubs"       , tN   ,   nN  ,  minN  ,  maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_ChiSquared   , "L1TkTau_IsoTks_ChiSquared"   , tChi ,  nChi ,  minChi,  maxChi );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_RedChiSquared, "L1TkTau_IsoTks_RedChiSquared", tRChi,  nRChi, minRChi,  maxRChi);
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_PtMinusCaloEt, "L1TkTau_IsoTks_PtMinusCaloEt",    "",   500 ,  -250.0,  +250.0 );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_DeltaR        , "L1TkTau_MatchTk_DeltaR"        , tDR  ,  nR  , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_PtRel         , "L1TkTau_MatchTk_PtRel"         , tPtR ,  nPtR, minPtR , maxPtR );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_Pt            , "L1TkTau_MatchTk_Pt"            , tPt  ,  nPt , minPt  , maxPt  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_Eta           , "L1TkTau_MatchTk_Eta"           , tEta ,  nEta, minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_POCAz         , "L1TkTau_MatchTk_POCAz"         , tZ0  ,  nZ0 , minZ0  , maxZ0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_NStubs        , "L1TkTau_MatchTk_NStubs"        , tN   ,   nN , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_ChiSquared    , "L1TkTau_MatchTk_ChiSquared"    , tChi ,  nChi, minChi , maxChi );
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_RedChiSquared , "L1TkTau_MatchTk_RedChiSquared" , tRChi, nRChi, minRChi, maxRChi);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsGenuine     , "L1TkTau_MatchTk_IsGenuine"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsUnknown     , "L1TkTau_MatchTk_IsUnknown"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_IsCombinatoric, "L1TkTau_MatchTk_IsCombinatoric",    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_MatchTk_PtMinusCaloEt , "L1TkTau_MatchTk_PtMinusCaloEt" ,    "",   500,  -250.0, +250.0 );

  // L1TkIsoTaus
  histoTools_.BookHisto_1D(hL1TkIsoTau_Multiplicity , "L1TkIsoTau_Multiplicity" , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloEt       , "L1TkIsoTau_CaloEt"       , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloEta      , "L1TkIsoTau_CaloEta"      , tEta , nEta , minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloPhi      , "L1TkIsoTau_CaloPhi"      , tPhi , nPhi , minPhi , maxPhi );
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloIEt      , "L1TkIsoTau_CaloIEt"      , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloIEta     , "L1TkIsoTau_CaloIEta"     , tEta , nIEta, minIEta, maxIEta);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloIPhi     , "L1TkIsoTau_CaloIPhi"     , tPhi , nIPhi, minIPhi, maxIPhi);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloIso      , "L1TkIsoTau_CaloIso"      ,   "" , nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloTowerIEta, "L1TkIsoTau_CaloTowerIEta", tEta , nIEta, minIEta, maxIEta);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloTowerIPhi, "L1TkIsoTau_CaloTowerIPhi", tPhi , nIPhi, minIPhi, maxIPhi);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloRawEt    , "L1TkIsoTau_CaloRawEt"    , tEt  , nEt*3, minEt  , maxEt*3);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloIsoEt    , "L1TkIsoTau_CaloIsoEt"    , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloNTT      , "L1TkIsoTau_CaloNTT"      , tN   ,nN*100, minN   , maxN*100);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloHasEM    , "L1TkIsoTau_CaloHasEM"    , ""   , nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CaloIsMerged , "L1TkIsoTau_CaloIsMerged" , ""   , nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_NSigTks      , "L1TkIsoTau_NSigTks"      , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTksEt     , "L1TkIsoTau_SigTksEt"     , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTksEta    , "L1TkIsoTau_SigTksEta"    , tEta , nEta , minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkIsoTau_NIsoTks      , "L1TkIsoTau_NIsoTks"      , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTksEt     , "L1TkIsoTau_IsoTksEt"     , tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTksEta    , "L1TkIsoTau_IsoTksEta"    , tEta , nEta , minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkIsoTau_InvMass      , "L1TkIsoTau_InvMass"      , tM   , nM   , minM   , maxM   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_InvMassIncl  , "L1TkIsoTau_InvMassIncl"  , tM   , nM   , minM   , maxM   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigConeRMin  , "L1TkIsoTau_SigConeRMin"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigConeRMax  , "L1TkIsoTau_SigConeRMax"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoConeRMin  , "L1TkIsoTau_IsoConeRMin"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoConeRMax  , "L1TkIsoTau_IsoConeRMax"  , tR   , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_Charge       , "L1TkIsoTau_Charge"       , tQ   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_RelIso       , "L1TkIsoTau_RelIso"       , tRIso, nRIso, minRIso, maxRIso);
  histoTools_.BookHisto_1D(hL1TkIsoTau_VtxIso       , "L1TkIsoTau_VtxIso"       , tVIso, nVIso, minVIso, maxVIso);
  histoTools_.BookHisto_1D(hL1TkIsoTau_DeltaRGenP   , "L1TkIsoTau_DeltaRGenP"   , tDR  , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_Rtau         , "L1TkIsoTau_Rtau"  , ";Entries / 0.1f; R_{#tau}", 100,  0.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_CHF          , "L1TkIsoTau_CHF"   , ";Entries / 0.1f; CHF"     , 500,  0.0,  +50.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_NHF          , "L1TkIsoTau_NHF"   , ";Entries / 0.1f; NHF"     , 200, -5.0,   +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_NHFAbs       , "L1TkIsoTau_NHFAbs", ";Entries / 0.1f; |NHF|"   , 100,  0.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_Pt           , "L1TkIsoTau_SigTks_Pt"           , tPt  ,  nPt  , minPt  , maxPt   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_Eta          , "L1TkIsoTau_SigTks_Eta"          , tEta ,  nEta , minEta , maxEta  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_POCAz        , "L1TkIsoTau_SigTks_POCAz"        , tZ0  ,  nZ0  , minZ0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_DeltaPOCAz   , "L1TkIsoTau_SigTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,     0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_PtRel        , "L1TkIsoTau_SigTks_PtRel"        , tPtR ,  nPtR , minPtR , maxPtR  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_DeltaR       , "L1TkIsoTau_SigTks_DeltaR"       , tDR  ,  nR   , minR   , maxR    );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_NStubs       , "L1TkIsoTau_SigTks_NStubs"       , tN   ,  nN   , minN   , maxN    );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_ChiSquared   , "L1TkIsoTau_SigTks_ChiSquared"   , tChi ,  nChi , minChi , maxChi  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_RedChiSquared, "L1TkIsoTau_SigTks_RedChiSquared", tRChi,  nRChi, minRChi, maxRChi );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_PtMinusCaloEt, "L1TkIsoTau_SigTks_PtMinusCaloEt",    "",   300 , -250.0 , +250.0  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_Pt           , "L1TkIsoTau_IsoTks_Pt"           , tPt  ,  nPt  ,  minPt ,  maxPt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_Eta          , "L1TkIsoTau_IsoTks_Eta"          , tEta ,  nEta ,  minEta,  maxEta );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_POCAz        , "L1TkIsoTau_IsoTks_POCAz"        , tZ0  ,  nZ0  ,  minZ0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_DeltaPOCAz   , "L1TkIsoTau_IsoTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,      0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_PtRel        , "L1TkIsoTau_IsoTks_PtRel"        , tPtR ,  nPtR ,  minPtR,  maxPtR );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_DeltaR       , "L1TkIsoTau_IsoTks_DeltaR"       , tDR  ,  nR   ,  minR  ,  maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_NStubs       , "L1TkIsoTau_IsoTks_NStubs"       , tN   ,   nN  ,  minN  ,  maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_ChiSquared   , "L1TkIsoTau_IsoTks_ChiSquared"   , tChi ,  nChi ,  minChi,  maxChi );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_RedChiSquared, "L1TkIsoTau_IsoTks_RedChiSquared", tRChi,  nRChi, minRChi,  maxRChi);
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_PtMinusCaloEt, "L1TkIsoTau_IsoTks_PtMinusCaloEt",    "",   500 ,  -250.0,  +250.0 );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_DeltaR        , "L1TkIsoTau_MatchTk_DeltaR"        , tDR  ,  nR  , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_PtRel         , "L1TkIsoTau_MatchTk_PtRel"         , tPtR ,  nPtR, minPtR , maxPtR );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_Pt            , "L1TkIsoTau_MatchTk_Pt"            , tPt  ,  nPt , minPt  , maxPt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_Eta           , "L1TkIsoTau_MatchTk_Eta"           , tEta ,  nEta, minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_POCAz         , "L1TkIsoTau_MatchTk_POCAz"         , tZ0  ,  nZ0 , minZ0  , maxZ0  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_NStubs        , "L1TkIsoTau_MatchTk_NStubs"        , tN   ,   nN , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_ChiSquared    , "L1TkIsoTau_MatchTk_ChiSquared"    , tChi ,  nChi, minChi , maxChi );
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_RedChiSquared , "L1TkIsoTau_MatchTk_RedChiSquared" , tRChi, nRChi, minRChi, maxRChi);
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_IsGenuine     , "L1TkIsoTau_MatchTk_IsGenuine"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_IsUnknown     , "L1TkIsoTau_MatchTk_IsUnknown"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_IsCombinatoric, "L1TkIsoTau_MatchTk_IsCombinatoric",    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_MatchTk_PtMinusCaloEt , "L1TkIsoTau_MatchTk_PtMinusCaloEt" ,    "",   500,  -250.0, +250.0 );

  // Resolutions 
  histoTools_.BookHisto_1D(hL1Tau_ResolutionCaloEt , "L1Tau_ResolutionCaloEt" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1Tau_ResolutionCaloEta, "L1Tau_ResolutionCaloEta", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1Tau_ResolutionCaloPhi, "L1Tau_ResolutionCaloPhi", ";#phi (rads);Events / %.2f rads", 200,  -10.0,  +10.0);
  // L1TkTau
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEt , "L1TkTau_ResolutionCaloEt" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEta, "L1TkTau_ResolutionCaloEta", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloPhi, "L1TkTau_ResolutionCaloPhi", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEt_C, "L1TkTau_ResolutionCaloEt_C" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEta_C, "L1TkTau_ResolutionCaloEta_C", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloPhi_C, "L1TkTau_ResolutionCaloPhi_C", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEt_I , "L1TkTau_ResolutionCaloEt_I" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEta_I, "L1TkTau_ResolutionCaloEta_I", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloPhi_I, "L1TkTau_ResolutionCaloPhi_I", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEt_F , "L1TkTau_ResolutionCaloEt_F" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloEta_F, "L1TkTau_ResolutionCaloEta_F", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkTau_ResolutionCaloPhi_F, "L1TkTau_ResolutionCaloPhi_F", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  // L1TkIsoTau
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEt , "L1TkIsoTau_ResolutionCaloEt" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEta, "L1TkIsoTau_ResolutionCaloEta", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloPhi, "L1TkIsoTau_ResolutionCaloPhi", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEt_C, "L1TkIsoTau_ResolutionCaloEt_C" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEta_C, "L1TkIsoTau_ResolutionCaloEta_C", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloPhi_C, "L1TkIsoTau_ResolutionCaloPhi_C", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEt_I , "L1TkIsoTau_ResolutionCaloEt_I" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEta_I, "L1TkIsoTau_ResolutionCaloEta_I", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloPhi_I, "L1TkIsoTau_ResolutionCaloPhi_I", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEt_F , "L1TkIsoTau_ResolutionCaloEt_F" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloEta_F, "L1TkIsoTau_ResolutionCaloEta_F", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionCaloPhi_F, "L1TkIsoTau_ResolutionCaloPhi_F", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);

  // SingleTau
  histoTools_.BookHisto_1D(hCalo_Rate    , "Calo_Rate"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_C  , "Calo_Rate_C"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_I  , "Calo_Rate_I"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_F  , "Calo_Rate_F"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Rate      , "Tk_Rate"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Rate_C    , "Tk_Rate_C"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Rate_I    , "Tk_Rate_I"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Rate_F    , "Tk_Rate_F"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Rate  , "VtxIso_Rate"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Rate_C, "VtxIso_Rate_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Rate_I, "VtxIso_Rate_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Rate_F, "VtxIso_Rate_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Rate  , "RelIso_Rate"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Rate_C, "RelIso_Rate_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Rate_I, "RelIso_Rate_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Rate_F, "RelIso_Rate_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Rate     , "Iso_Rate"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Rate_C   , "Iso_Rate_C"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Rate_I   , "Iso_Rate_I"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Rate_F   , "Iso_Rate_F"   , "", nEt , minEt , maxEt );
  
  histoTools_.BookHisto_1D(hCalo_Eff     , "Calo_Eff"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_C   , "Calo_Eff_C"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_I   , "Calo_Eff_I"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_F   , "Calo_Eff_F"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Eff       , "Tk_Eff"       , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Eff_C     , "Tk_Eff_C"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Eff_I     , "Tk_Eff_I"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_Eff_F     , "Tk_Eff_F"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Eff   , "VtxIso_Eff"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Eff_C , "VtxIso_Eff_C" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Eff_I , "VtxIso_Eff_I" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_Eff_F , "VtxIso_Eff_F" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Eff   , "RelIso_Eff"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Eff_C , "RelIso_Eff_C" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Eff_I , "RelIso_Eff_I" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_Eff_F , "RelIso_Eff_F" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Eff      , "Iso_Eff"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Eff_C    , "Iso_Eff_C"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Eff_I    , "Iso_Eff_I"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_Eff_F    , "Iso_Eff_F"    , "", nEt , minEt , maxEt );

  // DiTau
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo    , "DiTau_Rate_Calo"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo_C  , "DiTau_Rate_Calo_C"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo_I  , "DiTau_Rate_Calo_I"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Calo_F  , "DiTau_Rate_Calo_F"  , "", nEt , minEt , maxEt );  
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk      , "DiTau_Rate_Tk"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk_C    , "DiTau_Rate_Tk_C"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk_I    , "DiTau_Rate_Tk_I"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Tk_F    , "DiTau_Rate_Tk_F"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso  , "DiTau_Rate_VtxIso"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso_C, "DiTau_Rate_VtxIso_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso_I, "DiTau_Rate_VtxIso_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIso_F, "DiTau_Rate_VtxIso_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIso  , "DiTau_Rate_RelIso"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIso_C, "DiTau_Rate_RelIso_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIso_I, "DiTau_Rate_RelIso_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIso_F, "DiTau_Rate_RelIso_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Iso     , "DiTau_Rate_Iso"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Iso_C   , "DiTau_Rate_Iso_C"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Iso_I   , "DiTau_Rate_Iso_I"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_Iso_F   , "DiTau_Rate_Iso_F"   , "", nEt , minEt , maxEt );

  histoTools_.BookHisto_1D(hDiTau_Eff_Calo     , "DiTau_Eff_Calo"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo_C   , "DiTau_Eff_Calo_C"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo_I   , "DiTau_Eff_Calo_I"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Calo_F   , "DiTau_Eff_Calo_F"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk       , "DiTau_Eff_Tk"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk_C     , "DiTau_Eff_Tk_C"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk_I     , "DiTau_Eff_Tk_I"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Tk_F     , "DiTau_Eff_Tk_F"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso   , "DiTau_Eff_VtxIso"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso_C , "DiTau_Eff_VtxIso_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso_I , "DiTau_Eff_VtxIso_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIso_F , "DiTau_Eff_VtxIso_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIso   , "DiTau_Eff_RelIso"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIso_C , "DiTau_Eff_RelIso_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIso_I , "DiTau_Eff_RelIso_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIso_F , "DiTau_Eff_RelIso_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Iso      , "DiTau_Eff_Iso"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Iso_C    , "DiTau_Eff_Iso_C"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Iso_I    , "DiTau_Eff_Iso_I"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_Iso_F    , "DiTau_Eff_Iso_F"   , "", nEt , minEt , maxEt );
  
  // Turn-Ons
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt, "McHadronicTau_VisEt", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn50      , "Calo_TurnOn50"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn50        , "Tk_TurnOn50"        , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50    , "VtxIso_TurnOn50"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn50    , "RelIso_TurnOn50"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn50       , "Iso_TurnOn50"    , "", nEt , minEt , maxEt );

  histoTools_.BookHisto_1D(hCalo_TurnOn25  , "Calo_TurnOn25"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn25    , "Tk_TurnOn25"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25, "VtxIso_TurnOn25" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn25, "RelIso_TurnOn25" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn25   , "Iso_TurnOn25"    , "", nEt , minEt , maxEt );
  
  histoTools_.BookHisto_1D(hCalo_TurnOn_SingleTau50KHz  , "Calo_TurnOn_SingleTau50KHz"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn_SingleTau50KHz    , "Tk_TurnOn_SingleTau50KHz"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_SingleTau50KHz, "VtxIso_TurnOn_SingleTau50KHz", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn_SingleTau50KHz, "RelIso_TurnOn_SingleTau50KHz", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn_SingleTau50KHz   , "Iso_TurnOn_SingleTau50KHz"   , "", nEt , minEt , maxEt );

  histoTools_.BookHisto_1D(hCalo_TurnOn_DiTau50KHz  , "Calo_TurnOn_DiTau50KHz"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn_DiTau50KHz    , "Tk_TurnOn_DiTau50KHz"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_DiTau50KHz, "VtxIso_TurnOn_DiTau50KHz", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn_DiTau50KHz, "RelIso_TurnOn_DiTau50KHz", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn_DiTau50KHz   , "Iso_TurnOn_DiTau50KHz"   , "", nEt , minEt , maxEt );

  // DiTau (Calo-Other)
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_Tk    , "DiTau_Rate_Calo_Tk"    , "" , nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_VtxIso, "DiTau_Rate_Calo_VtxIso", "" , nEt, minEt, maxEt, nEt, minEt, maxEt); 
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_RelIso, "DiTau_Rate_Calo_RelIso", "" , nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Calo_Iso   , "DiTau_Rate_Calo_Iso"   , "" , nEt, minEt, maxEt, nEt, minEt, maxEt);

  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_Tk     , "DiTau_Eff_Calo_Tk"     , "" , nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_VtxIso , "DiTau_Eff_Calo_VtxIso" , "" , nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_RelIso , "DiTau_Eff_Calo_RelIso" , "" , nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Calo_Iso    , "DiTau_Eff_Calo_Iso"    , "" , nEt, minEt, maxEt, nEt, minEt, maxEt);

  // DiTau (Tk-Other)
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_VtxIso, "DiTau_Rate_Tk_VtxIso", "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_RelIso, "DiTau_Rate_Tk_RelIso", "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_Iso   , "DiTau_Rate_Tk_Iso"   , "", nEt, minEt, maxEt, nEt, minEt, maxEt);

  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_VtxIso , "DiTau_Eff_Tk_VtxIso" , "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_RelIso , "DiTau_Eff_Tk_RelIso" , "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_Iso    , "DiTau_Eff_Tk_Iso"    , "", nEt, minEt, maxEt, nEt, minEt, maxEt);

  return;
}

//============================================================================
void CaloTk::WriteHistos_(void)
//============================================================================
{
  // Location -> outFile 
  outFile->cd();
  
  // GenParticles Histograms
  h_GenP_VisET_dRMaxLdgPion->Write();

  // Counters
  hCounters->Write();

  // L1CaloTaus
  hL1CaloTau_Et->Write();
  hL1CaloTau_Eta->Write();
  hL1CaloTau_Phi->Write();
  hL1CaloTau_IEt->Write();
  hL1CaloTau_IEta->Write();
  hL1CaloTau_IPhi->Write();
  hL1CaloTau_Iso->Write();
  hL1CaloTau_TowerIEta->Write();
  hL1CaloTau_TowerIPhi->Write();
  hL1CaloTau_RawEt->Write();
  hL1CaloTau_IsoEt->Write();
  hL1CaloTau_NTT->Write();
  hL1CaloTau_HasEM->Write();
  hL1CaloTau_IsMerged->Write();

  // L1TkTaus: Matching track
  hL1TkTau_MatchTk_DeltaR->Write();
  hL1TkTau_MatchTk_PtRel->Write();
  hL1TkTau_MatchTk_Pt->Write();
  hL1TkTau_MatchTk_Eta->Write();
  hL1TkTau_MatchTk_NStubs->Write();
  hL1TkTau_MatchTk_ChiSquared->Write();
  hL1TkTau_MatchTk_RedChiSquared->Write();
  hL1TkTau_MatchTk_IsGenuine->Write();
  hL1TkTau_MatchTk_IsUnknown->Write();
  hL1TkTau_MatchTk_IsCombinatoric->Write();
  hL1TkTau_MatchTk_PtMinusCaloEt->Write();
  hL1TkTau_SigTks_Pt->Write();
  hL1TkTau_SigTks_PtRel->Write();
  hL1TkTau_SigTks_Eta->Write();
  hL1TkTau_SigTks_POCAz->Write();
  hL1TkTau_SigTks_DeltaPOCAz->Write();
  hL1TkTau_SigTks_DeltaR->Write();
  hL1TkTau_SigTks_NStubs->Write();
  hL1TkTau_SigTks_ChiSquared->Write();
  hL1TkTau_SigTks_RedChiSquared->Write();
  hL1TkTau_SigTks_PtMinusCaloEt->Write();
  hL1TkTau_IsoTks_Pt->Write();
  hL1TkTau_IsoTks_PtRel->Write();
  hL1TkTau_IsoTks_Eta->Write();
  hL1TkTau_IsoTks_POCAz->Write();
  hL1TkTau_IsoTks_DeltaPOCAz->Write();
  hL1TkTau_IsoTks_DeltaR->Write();
  hL1TkTau_IsoTks_NStubs->Write();
  hL1TkTau_IsoTks_ChiSquared->Write();
  hL1TkTau_IsoTks_RedChiSquared->Write();
  hL1TkTau_IsoTks_PtMinusCaloEt->Write();
  hL1TkTau_Multiplicity->Write();
  hL1TkTau_CaloEt->Write();
  hL1TkTau_CaloEta->Write();
  hL1TkTau_CaloPhi->Write();
  hL1TkTau_CaloIEt->Write();
  hL1TkTau_CaloIEta->Write();
  hL1TkTau_CaloIPhi->Write();
  hL1TkTau_CaloIso->Write();
  hL1TkTau_CaloTowerIEta->Write();
  hL1TkTau_CaloTowerIPhi->Write();
  hL1TkTau_CaloRawEt->Write();
  hL1TkTau_CaloIsoEt->Write();
  hL1TkTau_CaloNTT->Write();
  hL1TkTau_CaloHasEM->Write();
  hL1TkTau_CaloIsMerged->Write();
  hL1TkTau_Rtau->Write();
  hL1TkTau_CHF->Write();
  hL1TkTau_NHF->Write();
  hL1TkTau_NHFAbs->Write();
  hL1TkTau_NSigTks->Write();
  hL1TkTau_SigTksEt->Write();
  hL1TkTau_SigTksEta->Write();
  hL1TkTau_NIsoTks->Write();
  hL1TkTau_IsoTksEt->Write();
  hL1TkTau_IsoTksEta->Write();
  hL1TkTau_InvMass->Write();
  hL1TkTau_InvMassIncl->Write();
  hL1TkTau_SigConeRMin->Write();
  hL1TkTau_SigConeRMax->Write();
  hL1TkTau_IsoConeRMin->Write();
  hL1TkTau_IsoConeRMax->Write();
  hL1TkTau_Charge->Write();
  hL1TkTau_RelIso->Write();
  hL1TkTau_VtxIso->Write();
  hL1TkTau_DeltaRGenP->Write();

  // L1TkIsoTaus: Matching track
  hL1TkIsoTau_MatchTk_DeltaR->Write();
  hL1TkIsoTau_MatchTk_PtRel->Write();
  hL1TkIsoTau_MatchTk_Pt->Write();
  hL1TkIsoTau_MatchTk_Eta->Write();
  hL1TkIsoTau_MatchTk_NStubs->Write();
  hL1TkIsoTau_MatchTk_ChiSquared->Write();
  hL1TkIsoTau_MatchTk_RedChiSquared->Write();
  hL1TkIsoTau_MatchTk_IsGenuine->Write();
  hL1TkIsoTau_MatchTk_IsUnknown->Write();
  hL1TkIsoTau_MatchTk_IsCombinatoric->Write();
  hL1TkIsoTau_MatchTk_PtMinusCaloEt->Write();
  hL1TkIsoTau_SigTks_Pt->Write();
  hL1TkIsoTau_SigTks_PtRel->Write();
  hL1TkIsoTau_SigTks_Eta->Write();
  hL1TkIsoTau_SigTks_POCAz->Write();
  hL1TkIsoTau_SigTks_DeltaPOCAz->Write();
  hL1TkIsoTau_SigTks_DeltaR->Write();
  hL1TkIsoTau_SigTks_NStubs->Write();
  hL1TkIsoTau_SigTks_ChiSquared->Write();
  hL1TkIsoTau_SigTks_RedChiSquared->Write();
  hL1TkIsoTau_SigTks_PtMinusCaloEt->Write();
  hL1TkIsoTau_IsoTks_Pt->Write();
  hL1TkIsoTau_IsoTks_PtRel->Write();
  hL1TkIsoTau_IsoTks_Eta->Write();
  hL1TkIsoTau_IsoTks_POCAz->Write();
  hL1TkIsoTau_IsoTks_DeltaPOCAz->Write();
  hL1TkIsoTau_IsoTks_DeltaR->Write();
  hL1TkIsoTau_IsoTks_NStubs->Write();
  hL1TkIsoTau_IsoTks_ChiSquared->Write();
  hL1TkIsoTau_IsoTks_RedChiSquared->Write();
  hL1TkIsoTau_IsoTks_PtMinusCaloEt->Write();
  hL1TkIsoTau_Multiplicity->Write();
  hL1TkIsoTau_CaloEt->Write();
  hL1TkIsoTau_CaloEta->Write();
  hL1TkIsoTau_CaloPhi->Write();
  hL1TkIsoTau_CaloIEt->Write();
  hL1TkIsoTau_CaloIEta->Write();
  hL1TkIsoTau_CaloIPhi->Write();
  hL1TkIsoTau_CaloIso->Write();
  hL1TkIsoTau_CaloTowerIEta->Write();
  hL1TkIsoTau_CaloTowerIPhi->Write();
  hL1TkIsoTau_CaloRawEt->Write();
  hL1TkIsoTau_CaloIsoEt->Write();
  hL1TkIsoTau_CaloNTT->Write();
  hL1TkIsoTau_CaloHasEM->Write();
  hL1TkIsoTau_CaloIsMerged->Write();
  hL1TkIsoTau_Rtau->Write();
  hL1TkIsoTau_CHF->Write();
  hL1TkIsoTau_NHF->Write();
  hL1TkIsoTau_NHFAbs->Write();
  hL1TkIsoTau_NSigTks->Write();
  hL1TkIsoTau_SigTksEt->Write();
  hL1TkIsoTau_SigTksEta->Write();
  hL1TkIsoTau_NIsoTks->Write();
  hL1TkIsoTau_IsoTksEt->Write();
  hL1TkIsoTau_IsoTksEta->Write();
  hL1TkIsoTau_InvMass->Write();
  hL1TkIsoTau_InvMassIncl->Write();
  hL1TkIsoTau_SigConeRMin->Write();
  hL1TkIsoTau_SigConeRMax->Write();
  hL1TkIsoTau_IsoConeRMin->Write();
  hL1TkIsoTau_IsoConeRMax->Write();
  hL1TkIsoTau_Charge->Write();
  hL1TkIsoTau_RelIso->Write();
  hL1TkIsoTau_VtxIso->Write();
  hL1TkIsoTau_DeltaRGenP->Write();

  // Resolutions
  hL1Tau_ResolutionCaloEt->Write();
  hL1Tau_ResolutionCaloEta->Write();
  hL1Tau_ResolutionCaloPhi->Write();
  // L1TkTau
  hL1TkTau_ResolutionCaloEt->Write();
  hL1TkTau_ResolutionCaloEta->Write();
  hL1TkTau_ResolutionCaloPhi->Write();
  hL1TkTau_ResolutionCaloEt_C->Write();
  hL1TkTau_ResolutionCaloEta_C->Write();
  hL1TkTau_ResolutionCaloPhi_C->Write();
  hL1TkTau_ResolutionCaloEt_I->Write();
  hL1TkTau_ResolutionCaloEta_I->Write();
  hL1TkTau_ResolutionCaloPhi_I->Write();
  hL1TkTau_ResolutionCaloEt_F->Write();
  hL1TkTau_ResolutionCaloEta_F->Write();
  hL1TkTau_ResolutionCaloPhi_F->Write();
  // L1TkIsoTau
  hL1TkIsoTau_ResolutionCaloEt->Write();
  hL1TkIsoTau_ResolutionCaloEta->Write();
  hL1TkIsoTau_ResolutionCaloPhi->Write();
  hL1TkIsoTau_ResolutionCaloEt_C->Write();
  hL1TkIsoTau_ResolutionCaloEta_C->Write();
  hL1TkIsoTau_ResolutionCaloPhi_C->Write();
  hL1TkIsoTau_ResolutionCaloEt_I->Write();
  hL1TkIsoTau_ResolutionCaloEta_I->Write();
  hL1TkIsoTau_ResolutionCaloPhi_I->Write();
  hL1TkIsoTau_ResolutionCaloEt_F->Write();
  hL1TkIsoTau_ResolutionCaloEta_F->Write();
  hL1TkIsoTau_ResolutionCaloPhi_F->Write();

  // SingleTau: Rates
  hCalo_Rate->Write(); // Inclusive = C+I+F
  hCalo_Rate_C->Write();
  hCalo_Rate_I->Write();
  hCalo_Rate_F->Write();
  hTk_Rate->Write();
  hTk_Rate_C->Write();
  hTk_Rate_I->Write();
  hTk_Rate_F->Write();
  hVtxIso_Rate->Write();
  hVtxIso_Rate_C->Write();
  hVtxIso_Rate_I->Write();
  hVtxIso_Rate_F->Write();
  hRelIso_Rate->Write();
  hRelIso_Rate_C->Write();
  hRelIso_Rate_I->Write();
  hRelIso_Rate_F->Write();
  hIso_Rate->Write();
  hIso_Rate_C->Write();
  hIso_Rate_I->Write();
  hIso_Rate_F->Write();

  // SingleTau: Efficiencies
  hCalo_Eff->Write();  // Inclusive = C+I+F
  hCalo_Eff_C->Write();
  hCalo_Eff_I->Write();
  hCalo_Eff_F->Write();
  hTk_Eff->Write();
  hTk_Eff_C->Write();
  hTk_Eff_I->Write();
  hTk_Eff_F->Write();
  hVtxIso_Eff->Write();
  hVtxIso_Eff_C->Write();
  hVtxIso_Eff_I->Write();
  hVtxIso_Eff_F->Write();
  hRelIso_Eff->Write();
  hRelIso_Eff_C->Write();
  hRelIso_Eff_I->Write();
  hRelIso_Eff_F->Write();
  hIso_Eff->Write();
  hIso_Eff_C->Write();
  hIso_Eff_I->Write();
  hIso_Eff_F->Write();

  // DiTau: Rates
  hDiTau_Rate_Calo->Write(); // Inclusive = C+I+F
  hDiTau_Rate_Calo_C->Write();
  hDiTau_Rate_Calo_I->Write();
  hDiTau_Rate_Calo_F->Write();
  hDiTau_Rate_Tk->Write();
  hDiTau_Rate_Tk_C->Write();
  hDiTau_Rate_Tk_I->Write();
  hDiTau_Rate_Tk_F->Write();
  hDiTau_Rate_VtxIso->Write();
  hDiTau_Rate_VtxIso_C->Write();
  hDiTau_Rate_VtxIso_I->Write();
  hDiTau_Rate_VtxIso_F->Write();
  hDiTau_Rate_RelIso->Write();
  hDiTau_Rate_RelIso_C->Write();
  hDiTau_Rate_RelIso_I->Write();
  hDiTau_Rate_RelIso_F->Write();
  hDiTau_Rate_Iso->Write();
  hDiTau_Rate_Iso_C->Write();
  hDiTau_Rate_Iso_I->Write();
  hDiTau_Rate_Iso_F->Write();

  // DiTau: Efficiencies
  hDiTau_Eff_Calo->Write(); // Inclusive = C+I+F
  hDiTau_Eff_Calo_C->Write();
  hDiTau_Eff_Calo_I->Write();
  hDiTau_Eff_Calo_F->Write();
  hDiTau_Eff_Tk->Write();
  hDiTau_Eff_Tk_C->Write();
  hDiTau_Eff_Tk_I->Write();
  hDiTau_Eff_Tk_F->Write();
  hDiTau_Eff_VtxIso->Write();
  hDiTau_Eff_VtxIso_C->Write();
  hDiTau_Eff_VtxIso_I->Write();
  hDiTau_Eff_VtxIso_F->Write();
  hDiTau_Eff_RelIso->Write();
  hDiTau_Eff_RelIso_C->Write();
  hDiTau_Eff_RelIso_I->Write();
  hDiTau_Eff_RelIso_F->Write();
  hDiTau_Eff_Iso->Write();
  hDiTau_Eff_Iso_C->Write();
  hDiTau_Eff_Iso_I->Write();
  hDiTau_Eff_Iso_F->Write();

  // DiTau: (Calo-Other)
  hDiTau_Rate_Calo_Tk->Write();
  hDiTau_Rate_Calo_VtxIso->Write();
  hDiTau_Rate_Calo_RelIso->Write();
  hDiTau_Rate_Calo_Iso->Write();

  hDiTau_Eff_Calo_Tk->Write();
  hDiTau_Eff_Calo_VtxIso->Write();
  hDiTau_Eff_Calo_RelIso->Write();
  hDiTau_Eff_Calo_Iso->Write();

  // DiTau (Tk-Other)
  hDiTau_Rate_Tk_VtxIso->Write();
  hDiTau_Rate_Tk_RelIso->Write();
  hDiTau_Rate_Tk_Iso->Write();

  hDiTau_Eff_Tk_VtxIso->Write();
  hDiTau_Eff_Tk_RelIso->Write();
  hDiTau_Eff_Tk_Iso->Write();

  // Turn-Ons
  hMcHadronicTau_VisEt->Write();
  hCalo_TurnOn50->Write();
  hTk_TurnOn50->Write();
  hVtxIso_TurnOn50->Write();
  hRelIso_TurnOn50->Write();
  hIso_TurnOn50->Write();

  hCalo_TurnOn25->Write();
  hTk_TurnOn25->Write();
  hVtxIso_TurnOn25->Write();
  hRelIso_TurnOn25->Write();
  hIso_TurnOn25->Write();

  hCalo_TurnOn_SingleTau50KHz->Write();
  hTk_TurnOn_SingleTau50KHz->Write();
  hVtxIso_TurnOn_SingleTau50KHz->Write();
  hRelIso_TurnOn_SingleTau50KHz->Write();
  hIso_TurnOn_SingleTau50KHz->Write();

  hCalo_TurnOn_DiTau50KHz->Write();
  hTk_TurnOn_DiTau50KHz->Write();
  hVtxIso_TurnOn_DiTau50KHz->Write();
  hRelIso_TurnOn_DiTau50KHz->Write();
  hIso_TurnOn_DiTau50KHz->Write();

  // Write the outfile
  outFile->Write();

  return;
}


//****************************************************************************
bool CaloTk::IsWithinEtaRegion(string etaRegion,
                                 double eta)
//****************************************************************************
{

  bool bWithinEtaRegion = false;
  if ( etaRegion.compare("Central") == 0 )           bWithinEtaRegion = (fabs(eta) <= _eta_C);
  else if ( etaRegion.compare("Intermediate") == 0 ) bWithinEtaRegion = (fabs(eta) <= _eta_F && fabs(eta) > _eta_C);
  else if ( etaRegion.compare("Forward") == 0 )      bWithinEtaRegion = (fabs(eta) > _eta_F);
  else{
    cout << "=== Tracking::IsWithinEtaRegion() - Invalid eta region type \"" << etaRegion << "\". EXIT" << endl;
    exit(1);
  }

  return bWithinEtaRegion;
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
  L1TkTau.SetCaloTau(myL1Tau);
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
  
  // For-loop: All Tracks in isolation cone 
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
      if (tau->finalDaughtersCharged().size() < 1) continue; //iro

      double deltaR = auxTools_.DeltaR( L1TkTau.GetCaloTau().eta(), L1TkTau.GetCaloTau().phi(), tau->eta(), tau->phi() );
      if (deltaR > mcMatching_dRMax) continue;
      if (deltaR < match_dR)
	{
	  match_dR = deltaR;
	  match_GenParticle = *tau;
	}
      
    }  // For-loop: All hadronic GenTaus

  // Save the matching
  //cout << "index of match genp ="<<match_GenParticle.index()<<endl;
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

#endif

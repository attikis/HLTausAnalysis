#ifndef TkTaus_cxx
#define TkTaus_cxx

// User
#include "../Auxiliary/interface/constants.h"
#include "TkTaus.h"
//#include "../Framework/interface/TreeReaderReco.h"

// ROOT
#include "TFitResult.h"
#include "TF1.h"
#include "Math/VectorUtil.h"

// C++
#include <stdexcept>

//============================================================================
void TkTaus::InitObjects(void)
//============================================================================
{
  // pvProducer = new L1TkPrimaryVertex(this->s);
  return; 
}


//============================================================================
void TkTaus::InitVars_()
//============================================================================
{

  DEBUG = false;

  // Dataset-related
  datasets_  = datasets_.GetDataset(mcSample);
  realTauMom = datasets_.McTauMomPdgId_;
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;

  // Matching tracks
  seedTk_Collection  =  "TTTracks"; // "TTTracks"
  seedTk_nFitParams  =   4;         // 4
  seedTk_minPt       =   5.0;       //  10.0
  seedTk_minEta      =   0.0;       //   0.0
  seedTk_maxEta      =   2.5;       // 999.9
  seedTk_maxChiSq    =  50.0;       //  25.0
  seedTk_minStubs    =    5;        //   5

  // Signal cone tracks
  sigConeTks_Collection  = seedTk_Collection;
  sigConeTks_nFitParams  = seedTk_nFitParams;
  sigConeTks_minPt       =   2.0;  //   2.0
  sigConeTks_minEta      =   0.0;  //   0.0
  sigConeTks_maxEta      = 999.9;  // 999.9
  sigConeTks_maxChiSq    =  50.0;  //  25.0
  sigConeTks_minStubs    =   5;    //   4
  sigConeTks_dPOCAz      =   0.8;  // 0.80 (A. Ryd)
  sigConeTks_maxInvMass  =   1.5;  // 1.77 (A. Ryd)
 
  // Isolation cone tracks
  isoConeTks_Collection  = seedTk_Collection;
  isoConeTks_nFitParams  = seedTk_nFitParams;
  isoConeTks_minPt       =   2.0; //   2.0
  isoConeTks_minEta      =   0.0; //   0.0
  isoConeTks_maxEta      = 999.9; // 999.9
  isoConeTks_maxChiSq    = 100.0; // 100.00
  isoConeTks_minStubs    =   4;   //   4
  isoConeTks_dPOCAz      =   1.0; // 0.6 (A. Ryd) 

  // Signal cone parameters
  sigCone_Constant        = +0.00; // 0.0
  sigCone_dRMin           = +0.00; // WARNING! If > 0 the matching Track will NOT be added in sigCone_TTTracks
  sigCone_dRMax           = +0.15; // 0.15
  sigCone_cutoffDeltaR    = +0.15; // 0.15

  // Isolation cone
  isoCone_Constant = +0.65;         // 3.50 GeV (for ET_vis, 0
  isoCone_VtxIsoWP = +0.50;         // 1.0 cm
  isoCone_RelIsoWP = +0.50;         // 0.2
  isoCone_dRMin    = sigCone_dRMax; // 0.4
  isoCone_dRMax    = +0.35;         // 0.3
  diTau_deltaPOCAz = +1.00;         // 1.0 cm

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
void TkTaus::PrintSettings(void)
//============================================================================
{

  if (!DEBUG) return;

  // Inform user of settings
  Table settings("Variable | Cut | Value | Default | Units", "Text");
  settings.AddRowColumn(0, "MC Sample");
  settings.AddRowColumn(0, "==" );
  settings.AddRowColumn(0, mcSample );

  settings.AddRowColumn(1, "Matching Tracks: Collection");
  settings.AddRowColumn(1, "==");
  settings.AddRowColumn(1, seedTk_Collection);
  settings.AddRowColumn(1, "TTTracks");
  settings.AddRowColumn(1, "");
  
  settings.AddRowColumn(2, "Matching Tracks: Fit Parameters");
  settings.AddRowColumn(2, "==");
  settings.AddRowColumn(2, auxTools_.ToString( seedTk_nFitParams) );
  settings.AddRowColumn(2, "5");
  settings.AddRowColumn(2, "");

  settings.AddRowColumn(3, "Matching Tracks: Pt");
  settings.AddRowColumn(3, ">=");
  settings.AddRowColumn(3, auxTools_.ToString( seedTk_minPt) );
  settings.AddRowColumn(3, "2" );
  settings.AddRowColumn(3, "GeV/c" );
  
  settings.AddRowColumn(4, "Matching Tracks: |Eta|");
  settings.AddRowColumn(4, ">=");
  settings.AddRowColumn(4, auxTools_.ToString( seedTk_minEta) );
  settings.AddRowColumn(4, "0.0" );
  settings.AddRowColumn(4, "" );

  settings.AddRowColumn(5, "Matching Tracks: |Eta|");
  settings.AddRowColumn(5, "<=");
  settings.AddRowColumn(5, auxTools_.ToString( seedTk_maxEta) );
  settings.AddRowColumn(5, "1e+03" );
  settings.AddRowColumn(5, "" );
  
  settings.AddRowColumn(6, "Matching Tracks: ChiSqRed");
  settings.AddRowColumn(6, "<=");
  settings.AddRowColumn(6, auxTools_.ToString( seedTk_maxChiSq) );
  settings.AddRowColumn(6, "200/DOF"); // Cut was on ChiSq, not ChiSqRed
  settings.AddRowColumn(6, "");

  settings.AddRowColumn(7, "Matching Tracks: Stubs");
  settings.AddRowColumn(7, ">=");
  settings.AddRowColumn(7, auxTools_.ToString( seedTk_minStubs) );
  settings.AddRowColumn(7, "0" );
  settings.AddRowColumn(7, "" );

  settings.AddRowColumn(8, "Signal Cone Tks: Collection");
  settings.AddRowColumn(8, "==");
  settings.AddRowColumn(8, sigConeTks_Collection);
  settings.AddRowColumn(8, "TTTracks");
  settings.AddRowColumn(8, "");
  
  settings.AddRowColumn(9, "Signal Cone Tks: Fit Parameters");
  settings.AddRowColumn(9, "==");
  settings.AddRowColumn(9, auxTools_.ToString( sigConeTks_nFitParams) );
  settings.AddRowColumn(9, "5");
  settings.AddRowColumn(9, "");

  settings.AddRowColumn(10, "Signal Cone Tks: Pt");
  settings.AddRowColumn(10, ">=");
  settings.AddRowColumn(10, auxTools_.ToString( sigConeTks_minPt) );
  settings.AddRowColumn(10, "2" );
  settings.AddRowColumn(10, "GeV/c" );
  
  settings.AddRowColumn(11, "Signal Cone Tks: |Eta|");
  settings.AddRowColumn(11, ">=");
  settings.AddRowColumn(11, auxTools_.ToString( sigConeTks_minEta) );
  settings.AddRowColumn(11, "0.0" );
  settings.AddRowColumn(11, "" );

  settings.AddRowColumn(12, "Signal Cone Tks: |Eta|");
  settings.AddRowColumn(12, "<=");
  settings.AddRowColumn(12, auxTools_.ToString( sigConeTks_maxEta) );
  settings.AddRowColumn(12, "1e+03" );
  settings.AddRowColumn(12, "" );
  
  settings.AddRowColumn(13, "Signal Cone Tks: ChiSqRed");
  settings.AddRowColumn(13, "<=");
  settings.AddRowColumn(13, auxTools_.ToString( sigConeTks_maxChiSq) );
  settings.AddRowColumn(13, "200 (but on ChiSq, not ChiSqRed)");
  settings.AddRowColumn(13, "");

  settings.AddRowColumn(14, "Signal Cone Tks: Stubs");
  settings.AddRowColumn(14, ">=");
  settings.AddRowColumn(14, auxTools_.ToString( sigConeTks_minStubs) );
  settings.AddRowColumn(14, "" );
  settings.AddRowColumn(14, "" );

  settings.AddRowColumn(15, "Isolation Cone Tks: Collection");
  settings.AddRowColumn(15, "==");
  settings.AddRowColumn(15, isoConeTks_Collection);
  settings.AddRowColumn(15, "TTTracks");
  settings.AddRowColumn(15, "");
  
  settings.AddRowColumn(16, "Isolation Cone Tks: Fit Parameters");
  settings.AddRowColumn(16, "==");
  settings.AddRowColumn(16, auxTools_.ToString( isoConeTks_nFitParams) );
  settings.AddRowColumn(16, "5");
  settings.AddRowColumn(16, "");

  settings.AddRowColumn(17, "Isolation Cone Tks: Pt");
  settings.AddRowColumn(17, ">=");
  settings.AddRowColumn(17, auxTools_.ToString( isoConeTks_minPt) );
  settings.AddRowColumn(17, "2" );
  settings.AddRowColumn(17, "GeV/c" );
  
  settings.AddRowColumn(18, "Isolation Cone Tks: |Eta|");
  settings.AddRowColumn(18, ">=");
  settings.AddRowColumn(18, auxTools_.ToString( isoConeTks_minEta) );
  settings.AddRowColumn(18, "0.0" );
  settings.AddRowColumn(18, "" );

  settings.AddRowColumn(19, "Isolation Cone Tks: |Eta|");
  settings.AddRowColumn(19, "<=");
  settings.AddRowColumn(19, auxTools_.ToString( isoConeTks_maxEta) );
  settings.AddRowColumn(19, "1e+03" );
  settings.AddRowColumn(19, "" );

  settings.AddRowColumn(20, "Isolation Cone Tks: ChiSqRed");
  settings.AddRowColumn(20, "<=");
  settings.AddRowColumn(20, auxTools_.ToString( isoConeTks_maxChiSq) );
  settings.AddRowColumn(20, "200 (but on ChiSq, not ChiSqRed)");
  settings.AddRowColumn(20, "");

  settings.AddRowColumn(21, "Isolation Cone Tks: Stubs");
  settings.AddRowColumn(21, ">=");
  settings.AddRowColumn(21, auxTools_.ToString( isoConeTks_minStubs) );
  settings.AddRowColumn(21, "" );
  settings.AddRowColumn(21, "" );

  settings.AddRowColumn(22, "Signal Cone: Shrink Constant");
  settings.AddRowColumn(22, "==");
  settings.AddRowColumn(22, auxTools_.ToString(sigCone_Constant) );
  settings.AddRowColumn(22, "0" );
  settings.AddRowColumn(22, "GeV");

  settings.AddRowColumn(23, "Signal Cone: DeltaR");
  settings.AddRowColumn(23, ">=");
  settings.AddRowColumn(23, auxTools_.ToString(sigCone_dRMin) );
  settings.AddRowColumn(23, "0.0" );
  settings.AddRowColumn(23, "" );

  settings.AddRowColumn(24, "Signal Cone: DeltaR");
  settings.AddRowColumn(24, "<=");
  settings.AddRowColumn(24, auxTools_.ToString(sigCone_dRMax) );
  settings.AddRowColumn(24, "0.15" );
  settings.AddRowColumn(24, "" );
  
  settings.AddRowColumn(25, "Signal Cone:-3pr InvMass");
  settings.AddRowColumn(25, "<=");
  settings.AddRowColumn(25, auxTools_.ToString(sigConeTks_maxInvMass) );
  settings.AddRowColumn(25, "N/A" );
  settings.AddRowColumn(25, "GeV/c^{-2}");

  settings.AddRowColumn(26, "Signal Cone:-3pr dPOCAz");
  settings.AddRowColumn(26, "<=");
  settings.AddRowColumn(26, auxTools_.ToString(sigConeTks_dPOCAz) );
  settings.AddRowColumn(26, "N/A" );
  settings.AddRowColumn(26, "cm");

  settings.AddRowColumn(27, "Isolation Cone: Shrink Constant");
  settings.AddRowColumn(27, "==");
  settings.AddRowColumn(27, auxTools_.ToString(isoCone_Constant) );
  settings.AddRowColumn(27, "3.5");
  settings.AddRowColumn(27, "GeV");

  settings.AddRowColumn(28, "Isolation Cone: DeltaR");
  settings.AddRowColumn(28, ">=");
  settings.AddRowColumn(28, auxTools_.ToString(isoCone_dRMin) );
  settings.AddRowColumn(28, "0.15" );
  settings.AddRowColumn(28, "" );

  settings.AddRowColumn(29, "Isolation Cone: DeltaR");
  settings.AddRowColumn(29, "=<");
  settings.AddRowColumn(29, auxTools_.ToString(isoCone_dRMax) );
  settings.AddRowColumn(29, "0.30");
  settings.AddRowColumn(29, "");

  settings.AddRowColumn(30, "Isolation Cone: VtxIso" );
  settings.AddRowColumn(30, "<=" );
  settings.AddRowColumn(30, auxTools_.ToString(isoCone_VtxIsoWP) );
  settings.AddRowColumn(30, "1.0");
  settings.AddRowColumn(30, "cm");
  settings.AddRowColumn(30, "");

  settings.AddRowColumn(31, "Isolation Cone: RelIso" );
  settings.AddRowColumn(31, "<=" );
  settings.AddRowColumn(31, auxTools_.ToString(isoCone_RelIsoWP) );
  settings.AddRowColumn(31, "--");
  settings.AddRowColumn(31, "cm");
  settings.AddRowColumn(31, "");

  settings.AddRowColumn(32, "Di-Tau |Delta z0|");
  settings.AddRowColumn(32, "<");
  settings.AddRowColumn(32, auxTools_.ToString(diTau_deltaPOCAz) );
  settings.AddRowColumn(32, "1.0" );
  settings.AddRowColumn(32, "cm");

  settings.AddRowColumn(33, "MC-Matching DeltaR");
  settings.AddRowColumn(33, "<=");
  settings.AddRowColumn(33, auxTools_.ToString(mcMatching_dRMax) );
  settings.AddRowColumn(33, "0.05" );
  settings.AddRowColumn(33, "" );

  settings.AddRowColumn(34, "MC-Matching IsUnique");
  settings.AddRowColumn(34, "==");
  settings.AddRowColumn(34, auxTools_.ToString(mcMatching_unique) );
  settings.AddRowColumn(34, "1" );
  settings.AddRowColumn(34, "" );
  
  settings.AddRowColumn(35, "MC-Taus: Mom PdgId");
  settings.AddRowColumn(35, "==");
  settings.AddRowColumn(35, auxTools_.ToString(realTauMom));
  settings.AddRowColumn(35, "N/A" );
  settings.AddRowColumn(35, "" );

  settings.AddRowColumn(36, "MC-Taus: Number Expected");
  settings.AddRowColumn(36, ">=");
  settings.AddRowColumn(36, auxTools_.ToString(nMaxNumOfHTausPossible));
  settings.AddRowColumn(36, "N/A" );
  settings.AddRowColumn(36, "" );

  settings.AddRowColumn(37, "" );
  settings.Print();
  
  return;
}


//============================================================================
void TkTaus::Loop()
//============================================================================
{
  // Sanity check
  if (fChain == 0) return;
  
  const Long64_t nEntries   = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  if (DEBUG) cout << "=== TkTaus:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;

  // Initialisations
  InitVars_();
  BookHistos_();
  Long64_t nbytes       = 0;
  Long64_t nb           = 0;
  int nEvtsWithMaxHTaus = 0; 
  unsigned int nEvts    = 0;
  unsigned int nAllEvts = fChain->GetEntries();
  bool isMinBias        = false;  
  
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
    
    if (DEBUG) cout << "\tEntry = " << jentry << " / " << nEntries << endl;
    
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
	GenTausHadronic     = GetHadronicGenTaus(GenTaus, 00.0, 999.9); // for visEt and genP plots
	GenTausTrigger      = GetHadronicGenTaus(GenTaus, 20.0, 2.3);
      }

    if (DEBUG)
    {
	cout << "\tPrinting all GenParticle Collections" << endl;
	if (0) PrintGenParticleCollection(GenParticles);
	// PrintGenParticleCollection(GenTaus);
	// PrintGenParticleCollection(GenTausHadronic);    
	PrintGenParticleCollection(GenTausTrigger);
      }

    
    // Track Collections
    if(DEBUG) cout << "\tGetting the Tracks and Track Particles Collections" << endl;
    vector<TrackingParticle> TPs = GetTrackingParticles(false);
    
    vector<TTTrack> seedTracks = GetTTTracks(seedTk_minPt, seedTk_minEta, seedTk_maxEta, seedTk_maxChiSq,
					     seedTk_minStubs, seedTk_nFitParams, false);
    
    vector<TTTrack> sigTTTracks = GetTTTracks(sigConeTks_minPt, sigConeTks_minEta, sigConeTks_maxEta, 
					      sigConeTks_maxChiSq, sigConeTks_minStubs, 
					      sigConeTks_nFitParams, false);
    
    vector<TTTrack> isoTTTracks = GetTTTracks(isoConeTks_minPt , isoConeTks_minEta, isoConeTks_maxEta, 
					      isoConeTks_maxChiSq, isoConeTks_minStubs, 
					      isoConeTks_nFitParams, false);
    

    if (0) // DEBUG
      {
	cout << "\tPrinting all TTrack Collections" << endl;
	PrintTrackingParticleCollection(TPs);
	PrintTTTrackCollection(seedTracks);
	PrintTTTrackCollection(sigTTTracks);
	PrintTTTrackCollection(isoTTTracks);
      }

    // Tau Collections
    vector<L1TkTauParticle> L1TkTauCandidates;
    vector<L1TkTauParticle> L1TkTaus_Tk;
    vector<L1TkTauParticle> L1TkTaus_VtxIso;    
    vector<L1TkTauParticle> L1TkTaus_RelIso;
    vector<L1TkTauParticle> L1TkTaus_Iso;

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

	// Get the charged pions
	GetHadronicTauChargedPions(tau->index(), chargedPionsIndices);

	// Ask for 3-prong or 5-prong decay
	if (chargedPionsIndices.size() >= 3)
	  {
	    double pions_dRMax = -1000.0;
	    double ldgPionPt   = -1000.0;
	    int ldgPionIndx    = -1;
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
	    
	    // For-loop: All charged pions (3pr or 5pr)
	    for (vector<GenParticle>::iterator chPion = chargedPions.begin(); chPion != chargedPions.end(); chPion++)
	      {
		if (chPion->index() == ldgPionIndx) continue;
		double pions_dR = auxTools_.DeltaR(ldgPion.eta(), ldgPion.phi(), chPion->eta(), chPion->phi() );
		
		if ( pions_dR > pions_dRMax) pions_dRMax = pions_dR;
	      }
	    
	    // Fill histo
	    hGenP_VisEt_Vs_dRMaxLdgPion -> Fill(pions_dRMax, tau->p4vis().Et());
	    hGenP_PtLdg_Vs_dRMaxLdgPion -> Fill(pions_dRMax, ldgPionPt);
	    
	  }// Ask for 3-prong or 5-prong decay

      }


    // ======================================================================================
    // For-loop: Seed Tracks
    // ======================================================================================
    for (size_t i = 0; i < seedTracks.size(); i++)
      {
	TTTrack tk = seedTracks.at(i);

	// Calculate the Et-dependent signal & isolation cone sizes
	GetShrinkingConeSizes(tk.getPt(), sigCone_Constant, isoCone_Constant, 
			      sigCone_cutoffDeltaR, sigCone_dRMin, sigCone_dRMax, 
			      isoCone_dRMin, isoCone_dRMax);
	
	// Initialise the TkTau candidate (matching cone, signal cone edges, isolation cone edges)
	L1TkTauParticle	L1TkTauCandidate(0.0, 0.1, sigCone_dRMin, sigCone_dRMax, isoCone_dRMin, isoCone_dRMax);

	// The "matching track" is the seed track itself
	L1TkTauCandidate.SetMatchingTk(tk);
	L1TkTauCandidate.SetMatchTkDeltaRNew(0.0);
	if (0) L1TkTauCandidate.PrintProperties(false, false, true, false);

	//  Get signal-cone tracks
	GetSigConeTracks(L1TkTauCandidate, sigTTTracks, sigConeTks_dPOCAz, sigConeTks_maxInvMass);

	//  Get isolation-annulus tracks
	GetIsoConeTracks(L1TkTauCandidate, isoTTTracks, isoConeTks_dPOCAz);

	// Calculate isolation variables
	GetIsolationValues(L1TkTauCandidate);

	// Get the matching gen-particle
	GetMatchingGenParticle(L1TkTauCandidate, GenTausTrigger);

	// Print information on L1TkTauCandidate ??
	if (0) L1TkTauCandidate.PrintProperties(false, false, true, true);

	// Save L1TkTau Candidate
	L1TkTauCandidates.push_back(L1TkTauCandidate);
      }
    
    
    ////////////////////////////////////////////////
    /// Create Collections
    ////////////////////////////////////////////////
    for(vector<L1TkTauParticle>::iterator L1TkTau = L1TkTauCandidates.begin(); L1TkTau != L1TkTauCandidates.end(); L1TkTau++)
      {
	// +Tk
	if (L1TkTau->HasMatchingTk() )
	  {
	    
	    bool bIsLdgTrack = true;
	    vector<TTTrack> myTks;
	    vector<TTTrack> sigTks = L1TkTau->GetSigConeTTTracks();
	    vector<TTTrack> isoTks = L1TkTau->GetIsoConeTTTracks();
	    myTks.insert(myTks.end(), sigTks.begin(), sigTks.end());
	    myTks.insert(myTks.end(), isoTks.begin(), isoTks.end());

	    // For-loop: All signal tracks
	    for (vector<TTTrack>::iterator tk = myTks.begin(); tk != myTks.end(); tk++)
	      {
		double eta_seed = L1TkTau->GetMatchingTk().getEta(); // matchingTk = seeTk
		double phi_seed = L1TkTau->GetMatchingTk().getPhi();
		double eta_tk   = tk->getEta();
		double phi_tk   = tk->getPhi();
		double deltaPt  = L1TkTau->GetMatchingTk().getPt() - tk->getPt(); 

		// Skip identical tracks
		if ( (eta_seed == eta_tk) && (phi_seed == phi_tk) )
		  {
		    if (0) std::cout << "SAME TRACK! Continue ..." << std::endl;
		    continue;
		  }

		// Calculate dR
		double dR = auxTools_.DeltaR(eta_seed, phi_seed, eta_tk, phi_tk);

		// Consider only tracks within dR <= 0.3
		if (dR > 0.3) continue;
		
		// Compare pT of seed track with all tracks within dR = 0 (NEW)
		if (deltaPt < 0) 
		  {
		    if (0) std::cout << "Seed track not the leading track! Reject candidate! " << std::endl;
		    bIsLdgTrack = false;
		    break;
		  }
	      }

	    // No higher pT track within dR < 0.3
	    if (!bIsLdgTrack) continue;

	    // Save the tau candidates
	    L1TkTaus_Tk.push_back(*L1TkTau);
	    bool bPassVtxIso = (L1TkTau->GetVtxIsolation() > isoCone_VtxIsoWP);
	    bool bPassRelIso = (L1TkTau->GetRelIsolation() < isoCone_RelIsoWP);
	    bool bPassIso    = bPassVtxIso * bPassRelIso;
	    
	    // Fill containers with TkTaus
	    if (bPassVtxIso) L1TkTaus_VtxIso.push_back(*L1TkTau);
	    if (bPassRelIso) L1TkTaus_RelIso.push_back(*L1TkTau);
	    if (bPassIso)    L1TkTaus_Iso.push_back(*L1TkTau);
	  }
      }// L1TkTauCandidates

    if (DEBUG)
      {
	PrintL1TkTauParticleCollection(L1TkTaus_Tk);
	PrintL1TkTauParticleCollection(L1TkTaus_VtxIso);
	PrintL1TkTauParticleCollection(L1TkTaus_RelIso);
	PrintL1TkTauParticleCollection(L1TkTaus_Iso);
      }
          
    ////////////////////////////////////////////////
    /// TkTaus
    ////////////////////////////////////////////////
    vector<L1TkTauParticle> myL1TkTaus = L1TkTaus_Tk;
    hL1TkTau_Multiplicity ->Fill( myL1TkTaus.size() );
    unsigned int nMCTaus = 0;

    // For-loop: TkTaus
    for (vector<L1TkTauParticle>::iterator tau = myL1TkTaus.begin(); tau != myL1TkTaus.end(); tau++)
      {
	if (DEBUG) tau->PrintProperties(true, true, true, true);

	// Variables
	TLorentzVector sigTks_p4 = tau->GetSigConeTTTracksP4();
	TLorentzVector isoTks_p4 = tau->GetIsoConeTTTracksP4();

	// Do not skip if using MinBias sample as no real taus exist!
	if (DEBUG) std::cout << "=== Checking matching condition" << std::endl;
	if (!tau->HasMatchingGenParticle() && (isMinBias == false) ) continue;
	
	// Keep track of MC-matched taus
	nMCTaus++;

	// Get matching gen particle
	GenParticle p = tau->GetMatchingGenParticle(); //fixme: separate to inclusive, 1pr, and 3pr 
      
	// Seed Track Variables
	TTTrack matchTk   = tau->GetMatchingTk();
	double seedTk_dR  = auxTools_.DeltaR(matchTk.getEta(), matchTk.getPhi(), sigTks_p4.Eta(), sigTks_p4.Phi() );
	hL1TkTau_SeedTk_DeltaR        ->Fill( seedTk_dR );
	hL1TkTau_SeedTk_PtRel         ->Fill( matchTk.p3().Perp(sigTks_p4.Vect()) );
	hL1TkTau_SeedTk_Pt            ->Fill( matchTk.getPt() );
	hL1TkTau_SeedTk_Eta           ->Fill( matchTk.getEta() );
	hL1TkTau_SeedTk_POCAz         ->Fill( matchTk.getZ0() );
	hL1TkTau_SeedTk_NStubs        ->Fill( matchTk.getNumOfStubs() );
	hL1TkTau_SeedTk_ChiSquared    ->Fill( matchTk.getChi2() );
	hL1TkTau_SeedTk_RedChiSquared ->Fill( matchTk.getChi2Red() );
	hL1TkTau_SeedTk_IsGenuine     ->Fill( matchTk.getIsGenuine() );
	hL1TkTau_SeedTk_IsUnknown     ->Fill( matchTk.getIsUnknown() );
	hL1TkTau_SeedTk_IsCombinatoric->Fill( matchTk.getIsCombinatoric() );

	// Signal/Isolation Cone Variables
	double jetWidth = GetJetWidth(tau->GetSigConeTTTracks(), tau->GetIsoConeTTTracks(), sigTks_p4, isoTks_p4);
	hL1TkTau_JetWidth     ->Fill( jetWidth );
	hL1TkTau_DonutRatio   ->Fill( GetDonutRatio(*tau, isoTTTracks) );
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
	hL1TkTau_VtxIso_Vs_RelIso->Fill( tau->GetVtxIsolation(), tau->GetRelIsolation() );
	  
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

	    // Other variables
	    // sigTks_sumCharge += sigTk->getCharge(); // fixme
	    
	  }// SigCone_TTTracks

	// Fill histos for other variables
	hL1TkTau_Charge->Fill( sigTks_sumCharge); // fixme

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
  
	  }// IsoCone_TTTracks
      } // L1TkTaus_Tk
    
    // Fill MC-truth histos
    hL1TkTau_Multiplicity_MC ->Fill( nMCTaus );
 
    ////////////////////////////////////////////////
    /// L1TkIsoTau Properties 
    ////////////////////////////////////////////////
    vector<L1TkTauParticle> myL1TkIsoTaus = L1TkTaus_Iso; //L1TkTaus_VtxIso;
    unsigned int nMCIsoTaus = 0;
    hL1TkIsoTau_Multiplicity ->Fill( myL1TkIsoTaus.size() );

    // For-loop: All isolated tau candidates
    for (vector<L1TkTauParticle>::iterator tau = myL1TkIsoTaus.begin(); tau != myL1TkIsoTaus.end(); tau++)
      {
	
	if (DEBUG) tau->PrintProperties(true, true, true, true);

	// Variables
	TLorentzVector sigTks_p4 = tau->GetSigConeTTTracksP4();
	TLorentzVector isoTks_p4 = tau->GetIsoConeTTTracksP4();

	// Do not skip if using MinBias sample as no real taus exist!
	if (!tau->HasMatchingGenParticle() && (isMinBias == false) ) continue;
	
	// Keep track of MC-matched isolated taus
	nMCIsoTaus++;

	// Get matching gen particle
	GenParticle p = tau->GetMatchingGenParticle();  //fixme: separate to inclusive, 1pr, and 3pr
      
	// Resolution
	hL1TkIsoTau_ResolutionEt ->Fill( (tau->GetSigConeTTTracksP4().Et()-p.p4vis().Et() )/p.p4vis().Et()  );
	hL1TkIsoTau_ResolutionEta->Fill( (tau->GetSigConeTTTracksP4().Eta()-p.p4vis().Eta())/p.p4vis().Eta() );
	hL1TkIsoTau_ResolutionPhi->Fill( (tau->GetSigConeTTTracksP4().Phi()-p.p4vis().Phi())/p.p4vis().Phi() );
	
	double tauEta =tau->GetSigConeTTTracksP4().Eta();

	if ( IsWithinEtaRegion("Central", tauEta) )
	  {
	    hL1TkIsoTau_ResolutionEt_C ->Fill( (tau->GetSigConeTTTracksP4().Et()-p.p4vis().Et() )/p.p4vis().Et()  );
	    hL1TkIsoTau_ResolutionEta_C->Fill( (tau->GetSigConeTTTracksP4().Eta()-p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkIsoTau_ResolutionPhi_C->Fill( (tau->GetSigConeTTTracksP4().Phi()-p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	else if ( IsWithinEtaRegion("Intermediate", tauEta) )
	  {
	    hL1TkIsoTau_ResolutionEt_I ->Fill( (tau->GetSigConeTTTracksP4().Et()-p.p4vis().Et() )/p.p4vis().Et()  );
	    hL1TkIsoTau_ResolutionEta_I->Fill( (tau->GetSigConeTTTracksP4().Eta()-p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkIsoTau_ResolutionPhi_I->Fill( (tau->GetSigConeTTTracksP4().Phi()-p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	// currently no L1Taus in forward eta region
	else if ( IsWithinEtaRegion("Forward", tauEta) )
	  {
	    hL1TkIsoTau_ResolutionEt_F ->Fill( (tau->GetSigConeTTTracksP4().Et()-p.p4vis().Et() )/p.p4vis().Et()  );
	    hL1TkIsoTau_ResolutionEta_F->Fill( (tau->GetSigConeTTTracksP4().Eta()-p.p4vis().Eta())/p.p4vis().Eta() );
	    hL1TkIsoTau_ResolutionPhi_F->Fill( (tau->GetSigConeTTTracksP4().Phi()-p.p4vis().Phi())/p.p4vis().Phi() ); 
	  }
	else{                                                                                                                                                           
	  cout << "=== Tracking::Loop() - Unexpected Eta value of \"" << tauEta << "\". EXIT" << endl;
	  exit(1);
	}                    
	
	// Matching Track Variables
	TTTrack matchTk   = tau->GetMatchingTk();
	double seedTk_dR = auxTools_.DeltaR(matchTk.getEta(), matchTk.getPhi(), tau->GetSigConeTTTracksP4().Eta(), tau->GetSigConeTTTracksP4().Phi() ); // marina: can't we get the tau->GetMatchingTkDeltaR
	hL1TkIsoTau_SeedTk_DeltaR        ->Fill( seedTk_dR );
	hL1TkIsoTau_SeedTk_PtRel         ->Fill( matchTk.p3().Perp(sigTks_p4.Vect()) );
	hL1TkIsoTau_SeedTk_Pt            ->Fill( matchTk.getPt() );
	hL1TkIsoTau_SeedTk_Eta           ->Fill( matchTk.getEta() );
	hL1TkIsoTau_SeedTk_POCAz         ->Fill( matchTk.getZ0() );
	hL1TkIsoTau_SeedTk_NStubs        ->Fill( matchTk.getNumOfStubs() );
	hL1TkIsoTau_SeedTk_ChiSquared    ->Fill( matchTk.getChi2() );
	hL1TkIsoTau_SeedTk_RedChiSquared ->Fill( matchTk.getChi2Red() );
	hL1TkIsoTau_SeedTk_IsGenuine     ->Fill( matchTk.getIsGenuine() );
	hL1TkIsoTau_SeedTk_IsUnknown     ->Fill( matchTk.getIsUnknown() );
	hL1TkIsoTau_SeedTk_IsCombinatoric->Fill( matchTk.getIsCombinatoric() );

	// Signal/Isolation Cone Variables
	double jetWidth = GetJetWidth(tau->GetSigConeTTTracks(), tau->GetIsoConeTTTracks(), sigTks_p4, isoTks_p4);
	hL1TkIsoTau_JetWidth     ->Fill( jetWidth );
	hL1TkIsoTau_DonutRatio   ->Fill( GetDonutRatio(*tau, isoTTTracks) );
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
	hL1TkIsoTau_VtxIso_Vs_RelIso->Fill( tau->GetVtxIsolation(), tau->GetRelIsolation() );

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
  
	  }// IsoCone_TTTracks
      } // L1TkTaus_Iso

   // Fill MC-truth histos
    hL1TkIsoTau_Multiplicity_MC ->Fill( nMCIsoTaus );
 
    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////
    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++)
      {
	hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() );
      }

    FillTurnOn_Numerator_(L1TkTaus_Tk     , 25.0, hTk_TurnOn25    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn25);
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 25.0, hRelIso_TurnOn25);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 25.0, hIso_TurnOn25);

    FillTurnOn_Numerator_(L1TkTaus_Tk     , 50.0, hTk_TurnOn50    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn50);
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 50.0, hRelIso_TurnOn50);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 50.0, hIso_TurnOn50);

    FillTurnOn_Numerator_(L1TkTaus_Tk     , 65.0, hTk_TurnOn_SingleTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn_SingleTau50KHz);    
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 50.0, hRelIso_TurnOn_SingleTau50KHz);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 50.0, hIso_TurnOn_SingleTau50KHz);

    FillTurnOn_Numerator_(L1TkTaus_Tk     , 40.0, hTk_TurnOn_DiTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn_DiTau50KHz);
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 25.0, hRelIso_TurnOn_DiTau50KHz);
    FillTurnOn_Numerator_(L1TkTaus_Iso    , 25.0, hIso_TurnOn_DiTau50KHz);
    
    ////////////////////////////////////////////////
    // SingleTau
    ////////////////////////////////////////////////
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
    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_VtxIso, hDiTau_Rate_Tk_VtxIso, hDiTau_Eff_Tk_VtxIso   );
    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_RelIso, hDiTau_Rate_Tk_RelIso, hDiTau_Eff_Tk_RelIso   );
    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_Iso   , hDiTau_Rate_Tk_Iso   , hDiTau_Eff_Tk_Iso      );

    ////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_Tk);
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_VtxIso); 
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_RelIso);
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_Iso);

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
  double N = nEntries;
  histoTools_.ConvertToRateHisto_1D(hTk_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hTk_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hTk_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hTk_Rate_F, N);
      
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hVtxIso_Rate_F, N);

  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hRelIso_Rate_F, N);

  histoTools_.ConvertToRateHisto_1D(hIso_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hIso_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hIso_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hIso_Rate_F, N);

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
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Tk_F, N);
  
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIso_F, N);
  
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIso_F, N);

  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_Iso_F, N);
  
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

  // DiTau (Tk-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_VtxIso, N);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_RelIso, N);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_Iso   , N);

  FinaliseEffHisto_( hDiTau_Eff_Tk_VtxIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_RelIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_Iso   , nEvtsWithMaxHTaus);;

  // Turn-Ons: fixme. iro. alex
  // TEfficiency *pEff = 0;
  // pEff = new TEfficiency(*hCalo_TurnOn50_passed, *hMcHadronicTau_VisEt);
  // hCalo_TurnOn50 = (TH1D*) pEff->Clone();
  histoTools_.DivideHistos_1D(hTk_TurnOn50     , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn50 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hIso_TurnOn50    , hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hTk_TurnOn25     , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn25 , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hIso_TurnOn25    , hMcHadronicTau_VisEt);

  histoTools_.DivideHistos_1D(hTk_TurnOn_SingleTau50KHz    , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn_SingleTau50KHz, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn_SingleTau50KHz, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hIso_TurnOn_SingleTau50KHz   , hMcHadronicTau_VisEt);

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
void TkTaus::BookHistos_(void)
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

  const unsigned int nN = 20;
  const float minN =   0.0;
  const float maxN = +20.0;

  const unsigned int nBool = 2;
  const float minBool = -0.5;
  const float maxBool = +1.5;

  const unsigned int nM = 100;
  const float minM =  0.0;
  const float maxM = 10.0;

  const unsigned int nR = 100;
  const float minR = 0.0;
  const float maxR = 1.0;

  // const unsigned int nDR = 100;
  // const float minDR = 0.0;
  // const float maxDR = 1.0;

  const unsigned int nChi = 500;
  const float minChi =   0.0;
  const float maxChi = 500.0;

  const unsigned int nG = 1000;
  const float minG      =    0.0;
  const float maxG      =   10.0;

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
  const char* tIso = ";vertex isolation / %.2f cm;relativeisolation / %.2f";
  const char* tChi  = ";#chi^{2};Entries / %.2f";
  const char* tRChi = ";#chi^{2}_{#nu};Entries / %.2f";
  const char* tW    = ";w_{#tau};Entries / %.2f";
  const char* tG    = ";#gamma;Entries / %.2f";
  // const char* tRate = ";#E_{T} (GeV);Rate (kHz) / %.0f GeV";

  // GenParticles Histograms
  histoTools_.BookHisto_2D(hGenP_VisEt_Vs_dRMaxLdgPion, "GenP_VisEt_Vs_dRMaxLdgPion", ";#DeltaR_{max}(#pi_{ldg}^{#pm},#pi^{#pm});E_{T}^{vis}",  50,  0.0, +0.25, 100, 0.0, +200.0);
  histoTools_.BookHisto_2D(hGenP_PtLdg_Vs_dRMaxLdgPion, "GenP_PtLdg_Vs_dRMaxLdgPion", ";#DeltaR_{max}(#pi_{ldg}^{#pm},#pi^{#pm});p_{T}^{#pi_{ldg}^{#pm}}",  50,  0.0, +0.25, 100, 0.0, +200.0);

  // Counters
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 2, 0.0, +2.0);

  // L1TkTaus
  histoTools_.BookHisto_1D(hL1TkTau_Multiplicity , "L1TkTau_Multiplicity" , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_Multiplicity_MC, "L1TkTau_Multiplicity_MC", tN, nN, minN, maxN );
  histoTools_.BookHisto_1D(hL1TkTau_JetWidth     , "L1TkTau_JetWidth"     , tW   , nM   , minM   , maxM   );
  histoTools_.BookHisto_1D(hL1TkTau_DonutRatio   , "L1TkTau_DonutRatio"   , tG   , nG   , minG   , maxG   );
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
  histoTools_.BookHisto_2D(hL1TkTau_VtxIso_Vs_RelIso, "L1TkTau_VtxIso_Vs_RelIso", tIso, nVIso, minVIso, maxVIso, nRIso, minRIso, maxRIso);
  histoTools_.BookHisto_1D(hL1TkTau_DeltaRGenP   , "L1TkTau_DeltaRGenP"   , tDR  , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Pt           , "L1TkTau_SigTks_Pt"           , tPt  ,  nPt  , minPt  , maxPt   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_Eta          , "L1TkTau_SigTks_Eta"          , tEta ,  nEta , minEta , maxEta  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_POCAz        , "L1TkTau_SigTks_POCAz"        , tZ0  ,  nZ0  , minZ0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_DeltaPOCAz   , "L1TkTau_SigTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,     0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_PtRel        , "L1TkTau_SigTks_PtRel"        , tPtR ,  nPtR , minPtR , maxPtR  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_DeltaR       , "L1TkTau_SigTks_DeltaR"       , tDR  ,  nR   , minR   , maxR    );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_NStubs       , "L1TkTau_SigTks_NStubs"       , tN   ,  nN   , minN   , maxN    );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_ChiSquared   , "L1TkTau_SigTks_ChiSquared"   , tChi ,  nChi , minChi , maxChi  );
  histoTools_.BookHisto_1D(hL1TkTau_SigTks_RedChiSquared, "L1TkTau_SigTks_RedChiSquared", tRChi,  nRChi, minRChi, maxRChi );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Pt           , "L1TkTau_IsoTks_Pt"           , tPt  ,  nPt  ,  minPt ,  maxPt  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_Eta          , "L1TkTau_IsoTks_Eta"          , tEta ,  nEta ,  minEta,  maxEta );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_POCAz        , "L1TkTau_IsoTks_POCAz"        , tZ0  ,  nZ0  ,  minZ0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_DeltaPOCAz   , "L1TkTau_IsoTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,      0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_PtRel        , "L1TkTau_IsoTks_PtRel"        , tPtR ,  nPtR ,  minPtR,  maxPtR );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_DeltaR       , "L1TkTau_IsoTks_DeltaR"       , tDR  ,  nR   ,  minR  ,  maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_NStubs       , "L1TkTau_IsoTks_NStubs"       , tN   ,   nN  ,  minN  ,  maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_ChiSquared   , "L1TkTau_IsoTks_ChiSquared"   , tChi ,  nChi ,  minChi,  maxChi );
  histoTools_.BookHisto_1D(hL1TkTau_IsoTks_RedChiSquared, "L1TkTau_IsoTks_RedChiSquared", tRChi,  nRChi, minRChi,  maxRChi);
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_DeltaR        , "L1TkTau_SeedTk_DeltaR"        , tDR  ,  nR  , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_PtRel         , "L1TkTau_SeedTk_PtRel"         , tPtR ,  nPtR, minPtR , maxPtR );
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_Pt            , "L1TkTau_SeedTk_Pt"            , tPt  ,  nPt , minPt  , maxPt  );
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_Eta           , "L1TkTau_SeedTk_Eta"           , tEta ,  nEta, minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_POCAz         , "L1TkTau_SeedTk_POCAz"         , tZ0  ,  nZ0 , minZ0  , maxZ0  );
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_NStubs        , "L1TkTau_SeedTk_NStubs"        , tN   ,   nN , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_ChiSquared    , "L1TkTau_SeedTk_ChiSquared"    , tChi ,  nChi, minChi , maxChi );
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_RedChiSquared , "L1TkTau_SeedTk_RedChiSquared" , tRChi, nRChi, minRChi, maxRChi);
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_IsGenuine     , "L1TkTau_SeedTk_IsGenuine"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_IsUnknown     , "L1TkTau_SeedTk_IsUnknown"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkTau_SeedTk_IsCombinatoric, "L1TkTau_SeedTk_IsCombinatoric",    "", nBool, minBool, maxBool);

  // L1TkIsoTaus
  histoTools_.BookHisto_1D(hL1TkIsoTau_Multiplicity , "L1TkIsoTau_Multiplicity" , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_Multiplicity_MC , "L1TkIsoTau_Multiplicity_MC" , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_JetWidth     , "L1TkIsoTau_JetWidth"     , tW   , nM   , minM   , maxM   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_DonutRatio   , "L1TkIsoTau_DonutRatio"   , tG   , nG   , minG   , maxG   );
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
  histoTools_.BookHisto_2D(hL1TkIsoTau_VtxIso_Vs_RelIso, "L1TkIsoTau_VtxIso_Vs_RelIso", tIso, nVIso, minVIso, maxVIso, nRIso, minRIso, maxRIso);
  histoTools_.BookHisto_1D(hL1TkIsoTau_DeltaRGenP   , "L1TkIsoTau_DeltaRGenP"   , tDR  , nR   , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_Pt           , "L1TkIsoTau_SigTks_Pt"           , tPt  ,  nPt  , minPt  , maxPt   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_Eta          , "L1TkIsoTau_SigTks_Eta"          , tEta ,  nEta , minEta , maxEta  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_POCAz        , "L1TkIsoTau_SigTks_POCAz"        , tZ0  ,  nZ0  , minZ0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_DeltaPOCAz   , "L1TkIsoTau_SigTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,     0  , maxZ0   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_PtRel        , "L1TkIsoTau_SigTks_PtRel"        , tPtR ,  nPtR , minPtR , maxPtR  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_DeltaR       , "L1TkIsoTau_SigTks_DeltaR"       , tDR  ,  nR   , minR   , maxR    );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_NStubs       , "L1TkIsoTau_SigTks_NStubs"       , tN   ,  nN   , minN   , maxN    );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_ChiSquared   , "L1TkIsoTau_SigTks_ChiSquared"   , tChi ,  nChi , minChi , maxChi  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SigTks_RedChiSquared, "L1TkIsoTau_SigTks_RedChiSquared", tRChi,  nRChi, minRChi, maxRChi );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_Pt           , "L1TkIsoTau_IsoTks_Pt"           , tPt  ,  nPt  ,  minPt ,  maxPt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_Eta          , "L1TkIsoTau_IsoTks_Eta"          , tEta ,  nEta ,  minEta,  maxEta );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_POCAz        , "L1TkIsoTau_IsoTks_POCAz"        , tZ0  ,  nZ0  ,  minZ0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_DeltaPOCAz   , "L1TkIsoTau_IsoTks_DeltaPOCAz"   , tDZ0 ,  nZ0  ,      0 ,  maxZ0  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_PtRel        , "L1TkIsoTau_IsoTks_PtRel"        , tPtR ,  nPtR ,  minPtR,  maxPtR );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_DeltaR       , "L1TkIsoTau_IsoTks_DeltaR"       , tDR  ,  nR   ,  minR  ,  maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_NStubs       , "L1TkIsoTau_IsoTks_NStubs"       , tN   ,   nN  ,  minN  ,  maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_ChiSquared   , "L1TkIsoTau_IsoTks_ChiSquared"   , tChi ,  nChi ,  minChi,  maxChi );
  histoTools_.BookHisto_1D(hL1TkIsoTau_IsoTks_RedChiSquared, "L1TkIsoTau_IsoTks_RedChiSquared", tRChi,  nRChi, minRChi,  maxRChi);
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_DeltaR        , "L1TkIsoTau_SeedTk_DeltaR"        , tDR  ,  nR  , minR   , maxR   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_PtRel         , "L1TkIsoTau_SeedTk_PtRel"         , tPtR ,  nPtR, minPtR , maxPtR );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_Pt            , "L1TkIsoTau_SeedTk_Pt"            , tPt  ,  nPt , minPt  , maxPt  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_Eta           , "L1TkIsoTau_SeedTk_Eta"           , tEta ,  nEta, minEta , maxEta );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_POCAz         , "L1TkIsoTau_SeedTk_POCAz"         , tZ0  ,  nZ0 , minZ0  , maxZ0  );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_NStubs        , "L1TkIsoTau_SeedTk_NStubs"        , tN   ,   nN , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_ChiSquared    , "L1TkIsoTau_SeedTk_ChiSquared"    , tChi ,  nChi, minChi , maxChi );
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_RedChiSquared , "L1TkIsoTau_SeedTk_RedChiSquared" , tRChi, nRChi, minRChi, maxRChi);
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_IsGenuine     , "L1TkIsoTau_SeedTk_IsGenuine"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_IsUnknown     , "L1TkIsoTau_SeedTk_IsUnknown"     ,    "", nBool, minBool, maxBool);
  histoTools_.BookHisto_1D(hL1TkIsoTau_SeedTk_IsCombinatoric, "L1TkIsoTau_SeedTk_IsCombinatoric",    "", nBool, minBool, maxBool);

  // Resolutions 
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt , "L1TkIsoTau_ResolutionEt" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta, "L1TkIsoTau_ResolutionEta", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi, "L1TkIsoTau_ResolutionPhi", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt_C, "L1TkIsoTau_ResolutionEt_C" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta_C, "L1TkIsoTau_ResolutionEta_C", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi_C, "L1TkIsoTau_ResolutionPhi_C", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt_I , "L1TkIsoTau_ResolutionEt_I" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta_I, "L1TkIsoTau_ResolutionEta_I", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi_I, "L1TkIsoTau_ResolutionPhi_I", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt_F , "L1TkIsoTau_ResolutionEt_F" , ";E_{T} (GeV);Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta_F, "L1TkIsoTau_ResolutionEta_F", ";#eta;Events / %.2f", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi_F, "L1TkIsoTau_ResolutionPhi_F", ";#phi (rads);Events / %.2f rads", 100,  -10.0,  +10.0);

  // SingleTau
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
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt, "McHadronicTau_VisEt", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn50        , "Tk_TurnOn50"        , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50    , "VtxIso_TurnOn50"    , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn50    , "RelIso_TurnOn50"    , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn50       , "Iso_TurnOn50"       , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hTk_TurnOn25    , "Tk_TurnOn25"     , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25, "VtxIso_TurnOn25" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn25, "RelIso_TurnOn25" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn25   , "Iso_TurnOn25"    , "", 60 , minEt , maxEt );
  
  histoTools_.BookHisto_1D(hTk_TurnOn_SingleTau50KHz    , "Tk_TurnOn_SingleTau50KHz"    , "", 60, minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_SingleTau50KHz, "VtxIso_TurnOn_SingleTau50KHz", "", 60, minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn_SingleTau50KHz, "RelIso_TurnOn_SingleTau50KHz", "", 60, minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn_SingleTau50KHz   , "Iso_TurnOn_SingleTau50KHz"   , "", 60, minEt , maxEt );

  histoTools_.BookHisto_1D(hTk_TurnOn_DiTau50KHz    , "Tk_TurnOn_DiTau50KHz"    , "", 60, minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn_DiTau50KHz, "VtxIso_TurnOn_DiTau50KHz", "", 60, minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn_DiTau50KHz, "RelIso_TurnOn_DiTau50KHz", "", 60, minEt , maxEt );
  histoTools_.BookHisto_1D(hIso_TurnOn_DiTau50KHz   , "Iso_TurnOn_DiTau50KHz"   , "", 60, minEt , maxEt );

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
void TkTaus::WriteHistos_(void)
//============================================================================
{
  // Location -> outFile 
  outFile->cd();
  
  // GenParticles Histograms
  hGenP_VisEt_Vs_dRMaxLdgPion->Write();
  hGenP_PtLdg_Vs_dRMaxLdgPion->Write();

  // Counters
  hCounters->Write();

  // L1TkTaus: Matching track
  hL1TkTau_SeedTk_DeltaR->Write();
  hL1TkTau_SeedTk_PtRel->Write();
  hL1TkTau_SeedTk_Pt->Write();
  hL1TkTau_SeedTk_Eta->Write();
  hL1TkTau_SeedTk_NStubs->Write();
  hL1TkTau_SeedTk_POCAz->Write();
  hL1TkTau_SeedTk_ChiSquared->Write();
  hL1TkTau_SeedTk_RedChiSquared->Write();
  hL1TkTau_SeedTk_IsGenuine->Write();
  hL1TkTau_SeedTk_IsUnknown->Write();
  hL1TkTau_SeedTk_IsCombinatoric->Write();
  hL1TkTau_SigTks_Pt->Write();
  hL1TkTau_SigTks_PtRel->Write();
  hL1TkTau_SigTks_Eta->Write();
  hL1TkTau_SigTks_POCAz->Write();
  hL1TkTau_SigTks_DeltaPOCAz->Write();
  hL1TkTau_SigTks_DeltaR->Write();
  hL1TkTau_SigTks_NStubs->Write();
  hL1TkTau_SigTks_ChiSquared->Write();
  hL1TkTau_SigTks_RedChiSquared->Write();
  hL1TkTau_IsoTks_Pt->Write();
  hL1TkTau_IsoTks_PtRel->Write();
  hL1TkTau_IsoTks_Eta->Write();
  hL1TkTau_IsoTks_POCAz->Write();
  hL1TkTau_IsoTks_DeltaPOCAz->Write();
  hL1TkTau_IsoTks_DeltaR->Write();
  hL1TkTau_IsoTks_NStubs->Write();
  hL1TkTau_IsoTks_ChiSquared->Write();
  hL1TkTau_IsoTks_RedChiSquared->Write();
  hL1TkTau_Multiplicity->Write();
  hL1TkTau_Multiplicity_MC->Write();
  hL1TkTau_JetWidth->Write();
  hL1TkTau_DonutRatio->Write();
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
  hL1TkTau_VtxIso_Vs_RelIso->Write();
  hL1TkTau_DeltaRGenP->Write();

  // L1TkIsoTaus: Matching track
  hL1TkIsoTau_SeedTk_DeltaR->Write();
  hL1TkIsoTau_SeedTk_PtRel->Write();
  hL1TkIsoTau_SeedTk_Pt->Write();
  hL1TkIsoTau_SeedTk_Eta->Write();
  hL1TkIsoTau_SeedTk_NStubs->Write();
  hL1TkIsoTau_SeedTk_POCAz->Write();
  hL1TkIsoTau_SeedTk_ChiSquared->Write();
  hL1TkIsoTau_SeedTk_RedChiSquared->Write();
  hL1TkIsoTau_SeedTk_IsGenuine->Write();
  hL1TkIsoTau_SeedTk_IsUnknown->Write();
  hL1TkIsoTau_SeedTk_IsCombinatoric->Write();
  hL1TkIsoTau_SigTks_Pt->Write();
  hL1TkIsoTau_SigTks_PtRel->Write();
  hL1TkIsoTau_SigTks_Eta->Write();
  hL1TkIsoTau_SigTks_POCAz->Write();
  hL1TkIsoTau_SigTks_DeltaPOCAz->Write();
  hL1TkIsoTau_SigTks_DeltaR->Write();
  hL1TkIsoTau_SigTks_NStubs->Write();
  hL1TkIsoTau_SigTks_ChiSquared->Write();
  hL1TkIsoTau_SigTks_RedChiSquared->Write();
  hL1TkIsoTau_IsoTks_Pt->Write();
  hL1TkIsoTau_IsoTks_PtRel->Write();
  hL1TkIsoTau_IsoTks_Eta->Write();
  hL1TkIsoTau_IsoTks_POCAz->Write();
  hL1TkIsoTau_IsoTks_DeltaPOCAz->Write();
  hL1TkIsoTau_IsoTks_DeltaR->Write();
  hL1TkIsoTau_IsoTks_NStubs->Write();
  hL1TkIsoTau_IsoTks_ChiSquared->Write();
  hL1TkIsoTau_IsoTks_RedChiSquared->Write();
  hL1TkIsoTau_Multiplicity->Write();
  hL1TkIsoTau_Multiplicity_MC->Write();
  hL1TkIsoTau_JetWidth->Write();
  hL1TkIsoTau_DonutRatio->Write();
  hL1TkIsoTau_SigTksEt->Write();
  hL1TkIsoTau_SigTksEta->Write();
  hL1TkIsoTau_NSigTks->Write();
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
  hL1TkIsoTau_VtxIso_Vs_RelIso->Write();
  hL1TkIsoTau_DeltaRGenP->Write();

  // Resolutions
  hL1TkIsoTau_ResolutionEt->Write();
  hL1TkIsoTau_ResolutionEta->Write();
  hL1TkIsoTau_ResolutionPhi->Write();
  hL1TkIsoTau_ResolutionEt_C->Write();
  hL1TkIsoTau_ResolutionEta_C->Write();
  hL1TkIsoTau_ResolutionPhi_C->Write();
  hL1TkIsoTau_ResolutionEt_I->Write();
  hL1TkIsoTau_ResolutionEta_I->Write();
  hL1TkIsoTau_ResolutionPhi_I->Write();
  hL1TkIsoTau_ResolutionEt_F->Write();
  hL1TkIsoTau_ResolutionEta_F->Write();
  hL1TkIsoTau_ResolutionPhi_F->Write();

  // SingleTau: Rates
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

  // DiTau (Tk-Other)
  hDiTau_Rate_Tk_VtxIso->Write();
  hDiTau_Rate_Tk_RelIso->Write();
  hDiTau_Rate_Tk_Iso->Write();

  hDiTau_Eff_Tk_VtxIso->Write();
  hDiTau_Eff_Tk_RelIso->Write();
  hDiTau_Eff_Tk_Iso->Write();

  // Turn-Ons
  hMcHadronicTau_VisEt->Write();
  hTk_TurnOn50->Write();
  hVtxIso_TurnOn50->Write();
  hRelIso_TurnOn50->Write();
  hIso_TurnOn50->Write();

  hTk_TurnOn25->Write();
  hVtxIso_TurnOn25->Write();
  hRelIso_TurnOn25->Write();
  hIso_TurnOn25->Write();

  hTk_TurnOn_SingleTau50KHz->Write();
  hVtxIso_TurnOn_SingleTau50KHz->Write();
  hRelIso_TurnOn_SingleTau50KHz->Write();
  hIso_TurnOn_SingleTau50KHz->Write();

  hTk_TurnOn_DiTau50KHz->Write();
  hVtxIso_TurnOn_DiTau50KHz->Write();
  hRelIso_TurnOn_DiTau50KHz->Write();
  hIso_TurnOn_DiTau50KHz->Write();

  // Write the outfile
  outFile->Write();

  return;
}


//****************************************************************************
bool TkTaus::IsWithinEtaRegion(string etaRegion,
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
void TkTaus::FinaliseEffHisto_(TH1D *histo, 
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
void TkTaus::FinaliseEffHisto_(TH2D *histo, 
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
void TkTaus::ApplyDiTauZMatching(string tkCollectionType, 
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
	cout << "=== TkTaus::ApplyDiTauZMatching() - Unsupported track collection. Exit" << endl;
	exit(1);
      }
    else if ( tkCollectionType.compare("TTTracks") == 0 ) {
      deltaPOCAz = abs( match_tk0.getZ0() - match_tk.getZ0() );
    }
    else{
      cout << "=== TkTaus::ApplyDiTauZMatching() - Unknown sample \"" << mcSample << "\". EXIT" << endl;
      exit(1);
    }
    
    // If the Trigger objects is not within x-cm reject it
    if (deltaPOCAz > diTau_deltaPOCAz) L1TkTaus.erase ( L1TkTaus.begin()+i );
    
    }  // For-loop: L1TkTaus
  
  return;
}


//============================================================================
void TkTaus::FillSingleTau_(vector<L1TkTauParticle> L1TkTaus, 
			    TH1D *hRate,
			    TH1D *hEfficiency,
			    double minEta,
			    double maxEta)
//============================================================================
{

  // Sanity check
  if( L1TkTaus.size() < 1 ) return;
  
  // Fill rate
  TLorentzVector sigTks_p4 = L1TkTaus.at(0).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  double ldgEt  = sigTks_p4.Et();  
  double ldgEta = sigTks_p4.Eta();
  
  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(ldgEta) < minEta) return;
  if ( abs(ldgEta) > maxEta) return;
    
  FillRate_(hRate, ldgEt);
  
  // Get MC-matched trigger objects
  vector<L1TkTauParticle> L1TkTaus_mcMatched = GetMcMatchedL1TkTaus(L1TkTaus);
  if (L1TkTaus_mcMatched.size() < 1) return;

  // Check that all taus were found
  if(!bFoundAllTaus_) return;
  
  // Fill efficiency
  TLorentzVector sigTks_p4_mc = L1TkTaus_mcMatched.at(0).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  double ldgEt_mcMatched = sigTks_p4_mc.Et();
  FillEfficiency_(hEfficiency, ldgEt_mcMatched);

  return;
}


//============================================================================
void TkTaus::FillDiTau_(vector<L1TkTauParticle> L1TkTaus, 
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
  TLorentzVector sigTks_p4 = L1TkTaus.at(1).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  double subLdgEt  = sigTks_p4.Et();  
  double subLdgEta = sigTks_p4.Eta();  

  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(subLdgEta) < minEta) return;
  if ( abs(subLdgEta) > maxEta) return;
  if ( abs(subLdgEta) < minEta) return;
  if ( abs(subLdgEta) > maxEta) return;

  FillRate_(hRate, subLdgEt);

  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1TkTaus_mcMatched = GetMcMatchedL1TkTaus(L1TkTaus);
  if (L1TkTaus_mcMatched.size() < 2) return;
    
  // Check that all taus were found
  if(!bFoundAllTaus_) return;

  // Fill efficiency
  TLorentzVector sigTks_p4_mc = L1TkTaus_mcMatched.at(1).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  double subLdgEt_mcMatched  = sigTks_p4_mc.Et();  
  FillEfficiency_(hEfficiency, subLdgEt_mcMatched);

  return;
}


//============================================================================
void TkTaus::FillDiTau_(vector<L1TkTauParticle> L1TkTaus1,
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
  TLorentzVector sigTks1_p4 = L1TkTaus1.at(0).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  TLorentzVector sigTks2_p4 = L1TkTaus2.at(0).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  double ldgEt1 = sigTks1_p4.Et();
  double ldgEt2 = sigTks2_p4.Et();

  // Ensure that different calo objects are used
  double eta1 = L1TkTaus1.at(0).GetMatchingTk().getEta();
  double phi1 = L1TkTaus1.at(0).GetMatchingTk().getPhi();
  double eta2 = L1TkTaus2.at(0).GetMatchingTk().getEta();
  double phi2 = L1TkTaus2.at(0).GetMatchingTk().getPhi();
  double dR = auxTools_.DeltaR(eta1, phi1, eta2, phi2);
  if (dR < 0.4)
    {
      if (L1TkTaus2.size() < 2) return;
      ldgEt2 = L1TkTaus2.at(1).GetCaloTau().et(); //fixme! better clean-up needed
    }

  // Make x-axis the ldgEt axis
  if (ldgEt1 > ldgEt2) FillRate_(hRate, ldgEt1, ldgEt2); 
  else FillRate_(hRate, ldgEt2, ldgEt1);

  
  // Get MC-matched trigger objects
  if (L1TkTaus1_mcMatched.size() < 1) return;
  if (L1TkTaus2_mcMatched.size() < 1) return;

  // Get MC-matched Et
  TLorentzVector sigTks1_p4_mc = L1TkTaus1_mcMatched.at(0).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  TLorentzVector sigTks2_p4_mc = L1TkTaus2_mcMatched.at(0).GetSigConeTTTracksP4(); // fixme (filter before or sort with pT)?
  double ldgEt1_mcMatched = sigTks1_p4_mc.Et();
  double ldgEt2_mcMatched = sigTks2_p4_mc.Et();

  // Ensure that different calo objects are used
  eta1 = L1TkTaus1_mcMatched.at(0).GetMatchingTk().getEta();
  phi1 = L1TkTaus1_mcMatched.at(0).GetMatchingTk().getPhi();
  eta2 = L1TkTaus2_mcMatched.at(0).GetMatchingTk().getEta();
  phi2 = L1TkTaus2_mcMatched.at(0).GetMatchingTk().getPhi();
  dR    = auxTools_.DeltaR(eta1, phi1, eta2, phi2);
  if (dR < 0.4)
    {
      if (L1TkTaus2_mcMatched.size() < 2) return;
      ldgEt2_mcMatched = L1TkTaus2_mcMatched.at(1).GetCaloTau().et(); //fixme! better clean-up needed
    }
  
  // Check that all taus were found
  if(!bFoundAllTaus_) return;

  // Make x-axis the ldgEt axis
  if (ldgEt1_mcMatched > ldgEt2_mcMatched) histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt1_mcMatched, ldgEt2_mcMatched);
  else histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt2_mcMatched, ldgEt1_mcMatched);

  return;
}



//============================================================================
void TkTaus::FillRate_(TH1D *hRate,
		       const double ldgEt)
//============================================================================
{
  
  if (ldgEt < 0) return;
  hRate ->Fill( ldgEt );
  
  return;
}


//============================================================================
void TkTaus::FillRate_(TH2D *hRate,
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
void TkTaus::FillEfficiency_(TH1D *hEfficiency,
			     const double ldgEt)
//============================================================================
{
  
  histoTools_.FillAllBinsUpToValue_1D(hEfficiency, ldgEt);

  return;
}


//============================================================================
void TkTaus::FillTurnOn_Numerator_(vector<L1TkTauParticle> L1TkTaus, 
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
      TLorentzVector sigTks2_p4 = L1TkTau->GetSigConeTTTracksP4();
      double tkTau_et = sigTks2_p4.Et();
      if ( tkTau_et < minEt) continue;
      
      // Get MC-match
      GenParticle p = L1TkTau->GetMatchingGenParticle();
      
      // Fill the turn-on
      hTurnOn->Fill( p.p4vis().Et() );

    } // For-loop: L1TkTaus

  return;
   
}


//============================================================================
void TkTaus::GetSigConeTracks(L1TkTauParticle &L1TkTau,
			      vector<TTTrack> TTTracks,
			      double sigConeTks_dPOCAz,
			      double sigConeTks_maxInvMass)
//============================================================================
{
  if (!L1TkTau.HasMatchingTk()) return; 
  vector<TTTrack> sigConeTks;
  TTTrack seedTk = L1TkTau.GetMatchingTk();
  TLorentzVector sigTks_p4 = seedTk.p4();

  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
    {

      // Skip the matching track
      if( tk->index() == L1TkTau.GetMatchingTk().index() ) continue;
      
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), L1TkTau.GetMatchingTk().getEta(), L1TkTau.GetMatchingTk().getPhi());
      
      double dPOCAz = abs(seedTk.getZ0() - tk->getZ0());

      // Only consider tracks within singal cone/annulus
      if (dR < L1TkTau.GetSigConeMin()) continue;
      if (dR > L1TkTau.GetSigConeMax()) continue;

      // Require signal tracks come from same vertex as the seed track (NEW)
      if (dPOCAz > sigConeTks_dPOCAz) continue;

      // Apply invariant mass cut-off
      sigTks_p4 += tk->p4();
      if (sigTks_p4.M() > sigConeTks_maxInvMass)
	{
	  sigTks_p4 -= tk->p4();
	  continue;
	}

      // Save the track as a signal-cone track
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
void TkTaus::GetIsoConeTracks(L1TkTauParticle &L1TkTau,
			      vector<TTTrack> TTTracks,
			      double isoConeTks_dPOCAz)
//============================================================================
{
  if (!L1TkTau.HasMatchingTk()) return; 
  vector<TTTrack> isoConeTks;
  TTTrack seedTk = L1TkTau.GetMatchingTk();

  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = TTTracks.begin(); tk != TTTracks.end(); tk++)
    {
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), L1TkTau.GetMatchingTk().getEta(), L1TkTau.GetMatchingTk().getPhi());
      double dPOCAz = abs(seedTk.getZ0() - tk->getZ0());

      // Only consider tracks within singal cone
      if (dR < L1TkTau.GetIsoConeMin()) continue;
      if (dR > L1TkTau.GetIsoConeMax()) continue;

      // Require isolation tracks come from same vertex as the seed track (NEW)
      if (dPOCAz > isoConeTks_dPOCAz) continue;
      
      // Save the track as an isolation-cone track
      isoConeTks.push_back(*tk);
    }
  
  L1TkTau.SetIsoConeTracks(isoConeTks);

  return;
}


//============================================================================
double TkTaus::GetJetWidth(vector<TTTrack> sigTks,
			   vector<TTTrack> isoTks, 
			   TLorentzVector sigTks_p4,
			   TLorentzVector isoTks_p4)
//============================================================================
{
  double jetWidth = 0.0;
  
  vector<TTTrack> allTks;
  allTks.insert(allTks.end(), sigTks.begin(), sigTks.end());
  allTks.insert(allTks.end(), isoTks.begin(), isoTks.end());

  // For-loop: Signal Tracks
  // std::cout << "allTks.size() = " << allTks.size() << std::endl;
  for (vector<TTTrack>::iterator tk = allTks.begin(); tk != allTks.end(); tk++)
    {
      if (0) tk->PrintProperties();

      double pT = tk->getPt();
      double dR = auxTools_.DeltaR(sigTks_p4.Eta(), sigTks_p4.Phi(), tk->getEta(), tk->getPhi());

      // Calculate the jet widths
      // std::cout << "=== jetWidth += (" << pT << " x " << dR << ")/" << pT << std::endl;
      jetWidth += (pT*dR)/pT;
    }
  
  // std::cout << "=== jetWidth = " << jetWidth;
  return jetWidth;
}

//============================================================================
double TkTaus::GetDonutRatio(L1TkTauParticle &L1TkTau,
			     vector<TTTrack> isoTTTracks)
//============================================================================
{
  /*
    Check ratio of sum of pT of tracks in isolation annulus 
    wrt the sum of pT of tracks in another outer annulus 
    which has the same area as the isolation annulus. The idea
    is that if both are only populated by PU then the ratio should be 
    close to 1 (assumption is that PU energy distribution is isotropic in phi
    for a given eta ring)
   */

  TTTrack seedTk = L1TkTau.GetMatchingTk();  
  vector<TTTrack> isoConeTks = L1TkTau.GetIsoConeTTTracks();
  double sumPt_smallRAnnulus = 0.0;
  double sumPt_largeRAnnulus = 0.0;
  double smallR = L1TkTau.GetIsoConeMax();
  double largeR = sqrt(2)*smallR;
  double sumPt_ratio = 0.0;

  // For-loop: Isolation annulus Tracks
  for (vector<TTTrack>::iterator tk = isoConeTks.begin(); tk != isoConeTks.end(); tk++)
    {
      double pT = tk->getPt();
      sumPt_smallRAnnulus += pT;
    }

  // For-loop: Isolation-annulus-quality tracks
  for (vector<TTTrack>::iterator tk = isoTTTracks.begin(); tk != isoTTTracks.end(); tk++)
    {
      double pT = tk->getPt();
      double dR = auxTools_.DeltaR(seedTk.getEta(), seedTk.getPhi(), tk->getEta(), tk->getPhi());

      // Ensure track is inside the large annulus
      if (dR < smallR) continue;
      if (dR > largeR) continue;

      sumPt_largeRAnnulus += pT;
    }
  
  if (sumPt_smallRAnnulus > 0.0) sumPt_ratio = (sumPt_largeRAnnulus/sumPt_smallRAnnulus);
  // std::cout << "=== sumPt_ratio = (" << sumPt_largeRAnnulus << "/" << sumPt_smallRAnnulus << ") = " << sumPt_ratio << std::endl;
  return sumPt_ratio;
}

//============================================================================
void TkTaus::GetShrinkingConeSizes(double tk_pt,
				   double sigCone_Constant,
				   double isoCone_Constant,
				   const double sigCone_dRCutoff,
				   double &sigCone_dRMin,
				   double &sigCone_dRMax,
				   double &isoCone_dRMin,
				   double &isoCone_dRMax)
//============================================================================
{
  

  double signalCone_min = (sigCone_Constant)/(tk_pt);
  double signalCone_max = (isoCone_Constant)/(tk_pt);
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
void TkTaus::GetIsolationValues(L1TkTauParticle &L1TkTau)
//============================================================================
{

  // Store default values
  L1TkTau.SetVtxIsolation(999.9);
  L1TkTau.SetRelIsolation(0.0);
  
  // Return not Tk-Confirmed
  if (!L1TkTau.HasMatchingTk()) return; 

  // If no tracks found in the isoalation cone return
  vector<TTTrack> isoConeTks = L1TkTau.GetIsoConeTTTracks();
  if ( (isoConeTks.size() < 1) )  return;

  // Initialise variables
  TTTrack seedTk = L1TkTau.GetMatchingTk();
  double isoTks_scalarSumPt  = 0.0;
  double deltaZ0 = 999.9;
  double relIso  = 0.0;
  
  // For-loop: All Tracks in isolation cone 
  for (size_t i = 0; i < isoConeTks.size(); i++)
    {
      TTTrack isoConeTk = isoConeTks.at(i);
      
      // Add-up the pT of alltracks in isolation cone/annulus
      isoTks_scalarSumPt += isoConeTk.getPt();
      
      // Find the track closest in Z0
      deltaZ0 = abs(seedTk.getZ0() - isoConeTk.getZ0());
      if (deltaZ0 < L1TkTau.GetVtxIsolation())
	{
	  L1TkTau.SetVtxIsolation(deltaZ0);
	  L1TkTau.SetVtxIsolationTrack(isoConeTk);
	}
      
    }

  // Calculated + Assign value of relative isolation
  relIso = isoTks_scalarSumPt/seedTk.getPt();
  L1TkTau.SetRelIsolation(relIso);
  
  return;
}


//============================================================================
void TkTaus::GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
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

      if (0) tau->PrintFinalDaughtersCharged();

      TLorentzVector p4charged = tau->p4charged(false);
      double deltaR = auxTools_.DeltaR( p4charged.Eta(), p4charged.Phi(), matchTk.getEta(), matchTk.getPhi() );
      if (deltaR > mcMatching_dRMax) continue;
      if (deltaR < match_dR)
	{
	  match_dR = deltaR;
	  match_GenParticle = *tau;
	}
    }  // For-loop: All hadronic GenTaus
      
  // Save the matching (DataFormat/src/L1TkTauParticle.C and DataFormat/interface/L1TkTauParticle.h)
  if (match_dR <= mcMatching_dRMax)
    {
      L1TkTau.SetMatchingGenParticle(match_GenParticle);
      L1TkTau.SetMatchingGenParticleDeltaR(match_dR);
    }

  // For debugging
  if (DEBUG)
    {
      if (L1TkTau.HasMatchingGenParticle()) 
	{
	  L1TkTau.PrintProperties(false, true, false, false);
	  match_GenParticle.PrintProperties();
	}
    }
  
  return;
}      


//============================================================================
vector<L1TkTauParticle> TkTaus::GetMcMatchedL1TkTaus(vector<L1TkTauParticle> L1TkTaus)
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

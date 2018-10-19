#ifndef TkTaus_cxx
#define TkTaus_cxx

// User
#include "../Auxiliary/interface/constants.h"
#include "TkTaus.h"

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
  seedTk_maxChiSq    =  50.0;       //  50.0
  seedTk_minStubs    =    5;        //   5

  // Signal cone tracks
  sigConeTks_Collection  = seedTk_Collection;
  sigConeTks_nFitParams  = seedTk_nFitParams;
  sigConeTks_minPt       =   2.0;  //   2.0
  sigConeTks_minEta      =   0.0;  //   0.0
  sigConeTks_maxEta      =   2.5;  // 999.9
  sigConeTks_maxChiSq    =  50.0;  //  50.0
  sigConeTks_minStubs    =   5;    //   4
  sigConeTks_dPOCAz      =   1.0;  // 0.80 (A. Ryd)
  sigConeTks_maxInvMass  =   1.6;  // 1.77 (A. Ryd)
 
  // Isolation cone tracks
  isoConeTks_Collection  = seedTk_Collection;
  isoConeTks_nFitParams  = seedTk_nFitParams;
  isoConeTks_minPt       =   2.0; //   2.0
  isoConeTks_minEta      =   0.0; //   0.0
  isoConeTks_maxEta      =   2.5; // 999.9
  isoConeTks_maxChiSq    =  50.0; // 100.00
  isoConeTks_minStubs    =   4;   //   4

  // Signal cone parameters
  sigCone_Constant        = +0.00; // 0.0
  sigCone_dRMin           = +0.00; // WARNING! If > 0 the matching Track will NOT be added in sigCone_TTTracks
  sigCone_dRMax           = +0.25; // 0.20
  sigCone_cutoffDeltaR    = sigCone_dRMax; // ??? do i need this? (0.15)

  // Isolation cone
  isoCone_Constant = +2.5;          // 2.3 by fit on fit on ldg pT (Fotis)
  isoCone_dRMin    = sigCone_dRMax; // 0.4
  isoCone_dRMax    = +0.30;         // 0.30
  isoCone_useCone  = false; // instead of annulus

  // Isolation variables
  vtxIso_WP  = +0.50;  // 0.5 cm
  relIso_WP  = +0.20;  // 0.2
  relIso_dZ0 = +0.50;  // 0.6 from A. Ryd

  // Double-tau
  diTau_deltaPOCAz = +1.00; // cm

  // MC matching
  mcMatching_dRMax  = +0.10; // if too big get problem with turn-on curves (numerator>denominator)
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
  settings.AddRowColumn(30, auxTools_.ToString(vtxIso_WP) );
  settings.AddRowColumn(30, "1.0");
  settings.AddRowColumn(30, "cm");
  settings.AddRowColumn(30, "");

  settings.AddRowColumn(31, "Isolation Cone: RelIso" );
  settings.AddRowColumn(31, "<=" );
  settings.AddRowColumn(31, auxTools_.ToString(relIso_WP) );
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
  Long64_t nbytes = 0;
  Long64_t nb     = 0;
  bool isMinBias  = false;  
  int nEvtsWithMaxHTaus = 0; 
  unsigned int nEvts           = 0;
  unsigned int nEvtsSeedPt     = 0;
  unsigned int nEvtsSeedEta    = 0;
  unsigned int nEvtsSeedChiSq  = 0;
  unsigned int nEvtsSeedStubs  = 0;
  unsigned int nEvtsMcMatch    = 0; 
  unsigned int nEvtsVtxIso     = 0;
  unsigned int nEvtsRelIso     = 0;
  unsigned int nEvtsVtxIsoLoose= 0;
  unsigned int nEvtsVtxIsoTight= 0;
  unsigned int nEvtsRelIsoLoose= 0;
  unsigned int nEvtsRelIsoTight= 0;
  unsigned int nEvtsIso        = 0;
  unsigned int nAllEvts        = fChain->GetEntries();
  
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

    // Smart Counter
    vector<TTTrack> tmp;    
    tmp = GetTTTracks(seedTk_minPt, 0.0, 999.9, 999.9, 0, seedTk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedPt++;
    tmp = GetTTTracks(seedTk_minPt, seedTk_minEta, seedTk_maxEta, 999.9, 0, seedTk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedEta++;
    tmp = GetTTTracks(seedTk_minPt, seedTk_minEta, seedTk_maxEta, seedTk_maxChiSq, 0, seedTk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedChiSq++;
    tmp = GetTTTracks(seedTk_minPt, seedTk_minEta, seedTk_maxEta, seedTk_maxChiSq, seedTk_minStubs, seedTk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedStubs++;
    tmp.clear();

    vector<TTTrack> seedTTTracks = GetTTTracks(seedTk_minPt, seedTk_minEta, seedTk_maxEta, seedTk_maxChiSq,
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
	PrintTTTrackCollection(seedTTTracks);
	PrintTTTrackCollection(sigTTTracks);
	PrintTTTrackCollection(isoTTTracks);
      }

    // Tau Collections
    vector<L1TkTauParticle> L1TkTauCandidates;
    vector<L1TkTauParticle> L1TkTaus_Tk;
    vector<L1TkTauParticle> L1TkTaus_VtxIso;    
    vector<L1TkTauParticle> L1TkTaus_RelIso;
    vector<L1TkTauParticle> L1TkTaus_VtxIsoLoose;
    vector<L1TkTauParticle> L1TkTaus_VtxIsoTight;
    vector<L1TkTauParticle> L1TkTaus_RelIsoLoose;
    vector<L1TkTauParticle> L1TkTaus_RelIsoTight;

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
	    
	    // Fill histos
	    hGenP_VisEt_Vs_dRMaxLdgPion -> Fill(pions_dRMax, tau->p4vis().Et());
	    hGenP_PtLdg_Vs_dRMaxLdgPion -> Fill(pions_dRMax, ldgPionPt);
	    
	  }// Ask for 3-prong or 5-prong decay

      }


    // ======================================================================================
    // Construct the TkTau candidates from the seed tracks
    // ======================================================================================
    bool bFoundMC = false;

    // For-loop: Seed tracks
    for (size_t i = 0; i < seedTTTracks.size(); i++)
      {
	TTTrack tk = seedTTTracks.at(i);

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
	GetIsoConeTracks(L1TkTauCandidate, isoTTTracks, 999.99);

	// Calculate isolation variables
	L1TkTauCandidate.CalculateRelIso(relIso_dZ0, true, false, isoCone_useCone);
	L1TkTauCandidate.CalculateVtxIso(true, isoCone_useCone);

	// Get the matching gen-particle
	GetMatchingGenParticle(L1TkTauCandidate, GenTausTrigger); // GenTausHadronic
	if ( L1TkTauCandidate.HasMatchingGenParticle() ) bFoundMC = true;
	      
	// Print information on L1TkTauCandidate ??
	if (0) L1TkTauCandidate.PrintProperties(false, false, true, true);

	// Save L1TkTau Candidate
	L1TkTauCandidates.push_back(L1TkTauCandidate);
      }
	
    if (bFoundMC) nEvtsMcMatch++;
    
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

		// Consider only tracks within enitre jet definition
		if (dR > L1TkTau->GetIsoConeMax()) continue; //redundant but keep for safety
		// if (dR > L1TkTau->GetSigConeMax()) continue;
		
		// Compare pT of seed track with all tracks within dR = 0 (NEW)
		if (deltaPt < 0) 
		  {
		    if (0) std::cout << "Seed track not the leading track! Reject candidate! " << std::endl;
		    bIsLdgTrack = false;
		    break;
		  }
	      }

	    // No higher pT track within the entire jet (signal & isolation cone)
	    if (!bIsLdgTrack) continue;
	    
	    // Save the tau candidates
	    L1TkTaus_Tk.push_back(*L1TkTau);

	    // Calculate isolation variables
	    const double vtxIso  = L1TkTau->GetVtxIsolation(); // L1TkTau->CalculateVtxIso(false);
	    const double relIso  = L1TkTau->GetRelIsolation(); // L1TkTau->CalculateRelIso(0.5, false);
	    bool bPassVtxIso      = (vtxIso > vtxIso_WP); // orthogona1 to RelIso
	    bool bPassVtxIsoLoose = (vtxIso > 0.2);
	    bool bPassVtxIsoTight = (vtxIso > 1.0);
	    // bool bPassVtxIsoTight = (vtxIso > 0.5) && (L1TkTau->CalculateRelIso(relIso_dZ0, false, true) < 0.30);
	    bool bPassRelIso      = (relIso < relIso_WP); // orthogonal to VtxIso
	    bool bPassRelIsoLoose = (relIso < 0.3);
	    bool bPassRelIsoTight = (relIso < 0.15);
	      
	    // Fill containers with TkTaus
	    if (bPassVtxIso) L1TkTaus_VtxIso.push_back(*L1TkTau);
	    if (bPassRelIso) L1TkTaus_RelIso.push_back(*L1TkTau);
	    if (bPassVtxIsoLoose) L1TkTaus_VtxIsoLoose.push_back(*L1TkTau);
	    if (bPassVtxIsoTight) L1TkTaus_VtxIsoTight.push_back(*L1TkTau);
	    if (bPassRelIsoLoose) L1TkTaus_RelIsoLoose.push_back(*L1TkTau);
	    if (bPassRelIsoTight) L1TkTaus_RelIsoTight.push_back(*L1TkTau);
	  }
      }// L1TkTauCandidates

    // Counters
    if (L1TkTaus_VtxIso.size() > 0) nEvtsVtxIso++;
    if (L1TkTaus_RelIso.size() > 0) nEvtsRelIso++;
    if (L1TkTaus_VtxIsoLoose.size() > 0) nEvtsVtxIsoLoose++;
    if (L1TkTaus_VtxIsoTight.size() > 0) nEvtsVtxIsoTight++;
    if (L1TkTaus_RelIsoLoose.size() > 0) nEvtsRelIsoLoose++;
    if (L1TkTaus_RelIsoTight.size() > 0) nEvtsRelIsoTight++;

    if (DEBUG)
      {
	PrintL1TkTauParticleCollection(L1TkTaus_Tk);
	PrintL1TkTauParticleCollection(L1TkTaus_VtxIso);
	PrintL1TkTauParticleCollection(L1TkTaus_RelIso);
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
	GenParticle p = tau->GetMatchingGenParticle();
      
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
	hL1TkTau_DonutRatio_Vs_JetWidth->Fill( GetDonutRatio(*tau, isoTTTracks) , jetWidth );
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
	//if (  tau->GetIsoConeTTTracks().size() > 0 ) 
	//{
	// hL1TkTau_VtxIso_Vs_RelIso->Fill( tau->GetVtxIsolation(), tau->GetRelIsolation() ); // relIso shown only up to vtxIso = relIso_dZ0 => not very informative
	hL1TkTau_VtxIso_Vs_RelIso->Fill( tau->GetVtxIsolation(), tau->CalculateRelIso(999.9, false) ); // shows entire range for relIso
	//}
	
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
  
	  }// IsoCone_TTTracks
      } // L1TkTaus_Tk
    
    // Fill MC-truth histos
    hL1TkTau_Multiplicity_MC ->Fill( nMCTaus );
 
    ////////////////////////////////////////////////
    /// L1TkIsoTau Properties 
    ////////////////////////////////////////////////
    vector<L1TkTauParticle> myL1TkIsoTaus = L1TkTaus_VtxIso;
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
	GenParticle p = tau->GetMatchingGenParticle();
	double etRes  = (tau->GetSigConeTTTracksP4().Et()-p.p4vis().Et() )/p.p4vis().Et();
	double etaRes = (tau->GetSigConeTTTracksP4().Eta()-p.p4vis().Eta())/p.p4vis().Eta();
	double phiRes = (tau->GetSigConeTTTracksP4().Phi()-p.p4vis().Phi())/p.p4vis().Phi();

	// Resolution
	hL1TkIsoTau_ResolutionEt ->Fill( etRes  );
	hL1TkIsoTau_ResolutionEta->Fill( etaRes );
	hL1TkIsoTau_ResolutionPhi->Fill( phiRes );

	// 
	if (p.finalDaughtersNeutral().size() > 0)
	  {
	    hL1TkIsoTau_ResolutionEt_withNeutrals ->Fill( etRes  );
	    hL1TkIsoTau_ResolutionEta_withNeutrals->Fill( etaRes );
	    hL1TkIsoTau_ResolutionPhi_withNeutrals->Fill( phiRes );
	  }
	else{
	  hL1TkIsoTau_ResolutionEt_noNeutrals ->Fill( etRes  );
	  hL1TkIsoTau_ResolutionEta_noNeutrals->Fill( etaRes );
	  hL1TkIsoTau_ResolutionPhi_noNeutrals->Fill( phiRes );
	}

	if (p.finalDaughtersCharged().size() == 1) 
	  {
	    hL1TkIsoTau_ResolutionEt_1pr ->Fill( etRes  );
	    hL1TkIsoTau_ResolutionEta_1pr->Fill( etaRes );
	    hL1TkIsoTau_ResolutionPhi_1pr->Fill( phiRes );
	  }
	else if (p.finalDaughtersCharged().size() == 3) 
	  {
	    hL1TkIsoTau_ResolutionEt_3pr ->Fill( etRes  );
	    hL1TkIsoTau_ResolutionEta_3pr->Fill( etaRes );
	    hL1TkIsoTau_ResolutionPhi_3pr->Fill( phiRes );
	  }

	double tauEta =tau->GetSigConeTTTracksP4().Eta();

	if ( IsWithinEtaRegion("Central", tauEta) )
	  {
	    hL1TkIsoTau_ResolutionEt_C ->Fill( etRes  );
	    hL1TkIsoTau_ResolutionEta_C->Fill( etaRes );
	    hL1TkIsoTau_ResolutionPhi_C->Fill( phiRes ); 
	  }
	else if ( IsWithinEtaRegion("Intermediate", tauEta) )
	  {
	    hL1TkIsoTau_ResolutionEt_I ->Fill( etRes  );
	    hL1TkIsoTau_ResolutionEta_I->Fill( etaRes );
	    hL1TkIsoTau_ResolutionPhi_I->Fill( phiRes ); 
	  }
	// currently no L1Taus in forward eta region
	else if ( IsWithinEtaRegion("Forward", tauEta) )
	  {
	    hL1TkIsoTau_ResolutionEt_F ->Fill( etRes  );
	    hL1TkIsoTau_ResolutionEta_F->Fill( etaRes );
	    hL1TkIsoTau_ResolutionPhi_F->Fill( phiRes ); 
	  }
	else{
	  cout << "=== TkTaus::Loop() - Unexpected Eta value of \"" << tauEta << "\". EXIT" << endl;
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
	hL1TkIsoTau_DonutRatio_Vs_JetWidth->Fill( GetDonutRatio(*tau, isoTTTracks) , jetWidth );
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
	if (  tau->GetIsoConeTTTracks().size() > 0 ) 
	  {
	    hL1TkIsoTau_VtxIso_Vs_RelIso->Fill( tau->GetVtxIsolation(), tau->GetRelIsolation() );
	  }

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
      } // myL1TkIsoTaus

   // Fill MC-truth histos
    hL1TkIsoTau_Multiplicity_MC ->Fill( nMCIsoTaus );
 
    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////
    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++)
      {
	hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() ); // turn-on fill
	
	if (tau->finalDaughtersNeutral().size() > 0)
	  {
	    hMcHadronicTau_VisEt_withNeutrals->Fill( tau->p4vis().Et() );
	  }
	else
	  {
	    hMcHadronicTau_VisEt_noNeutrals->Fill( tau->p4vis().Et() );
	  }
	if (tau->finalDaughtersCharged().size() == 1) 
	  {
	    hMcHadronicTau_VisEt_1pr->Fill( tau->p4vis().Et() );
	  }
	else if (tau->finalDaughtersCharged().size() == 3) 
	  {
	    hMcHadronicTau_VisEt_3pr->Fill( tau->p4vis().Et() );
	  }
      }

    FillTurnOn_Numerator_(L1TkTaus_Tk     , 25.0, hTk_TurnOn25, hTk_TurnOn25_1pr, hTk_TurnOn25_3pr, hTk_TurnOn25_withNeutrals, hTk_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn25, hVtxIso_TurnOn25_1pr, hVtxIso_TurnOn25_3pr, hVtxIso_TurnOn25_withNeutrals, hVtxIso_TurnOn25_noNeutrals); 
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 25.0, hRelIso_TurnOn25, hRelIso_TurnOn25_1pr, hRelIso_TurnOn25_3pr, hRelIso_TurnOn25_withNeutrals, hRelIso_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_VtxIsoLoose , 25.0, hVtxIsoLoose_TurnOn25, hVtxIsoLoose_TurnOn25_1pr, hVtxIsoLoose_TurnOn25_3pr, hVtxIsoLoose_TurnOn25_withNeutrals, hVtxIsoLoose_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_VtxIsoTight , 25.0, hVtxIsoTight_TurnOn25, hVtxIsoTight_TurnOn25_1pr, hVtxIsoTight_TurnOn25_3pr, hVtxIsoTight_TurnOn25_withNeutrals, hVtxIsoTight_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_RelIsoLoose , 25.0, hRelIsoLoose_TurnOn25, hRelIsoLoose_TurnOn25_1pr, hRelIsoLoose_TurnOn25_3pr, hRelIsoLoose_TurnOn25_withNeutrals, hRelIsoLoose_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_RelIsoTight , 25.0, hRelIsoTight_TurnOn25, hRelIsoTight_TurnOn25_1pr, hRelIsoTight_TurnOn25_3pr, hRelIsoTight_TurnOn25_withNeutrals, hRelIsoTight_TurnOn25_noNeutrals);

    FillTurnOn_Numerator_(L1TkTaus_Tk     , 50.0, hTk_TurnOn50, hTk_TurnOn50_1pr, hTk_TurnOn50_3pr, hTk_TurnOn50_withNeutrals, hTk_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn50, hVtxIso_TurnOn50_1pr, hVtxIso_TurnOn50_3pr, hVtxIso_TurnOn50_withNeutrals, hVtxIso_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_RelIso , 50.0, hRelIso_TurnOn50, hRelIso_TurnOn50_1pr, hRelIso_TurnOn50_3pr, hRelIso_TurnOn50_withNeutrals, hRelIso_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_VtxIsoLoose , 50.0, hVtxIsoLoose_TurnOn50, hVtxIsoLoose_TurnOn50_1pr, hVtxIsoLoose_TurnOn50_3pr, hVtxIsoLoose_TurnOn50_withNeutrals, hVtxIsoLoose_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_VtxIsoTight , 50.0, hVtxIsoTight_TurnOn50, hVtxIsoTight_TurnOn50_1pr, hVtxIsoTight_TurnOn50_3pr, hVtxIsoTight_TurnOn50_withNeutrals, hVtxIsoTight_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_RelIsoLoose , 50.0, hRelIsoLoose_TurnOn50, hRelIsoLoose_TurnOn50_1pr, hRelIsoLoose_TurnOn50_3pr, hRelIsoLoose_TurnOn50_withNeutrals, hRelIsoLoose_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkTaus_RelIsoTight , 50.0, hRelIsoTight_TurnOn50, hRelIsoTight_TurnOn50_1pr, hRelIsoTight_TurnOn50_3pr, hRelIsoTight_TurnOn50_withNeutrals, hRelIsoTight_TurnOn50_noNeutrals);
    
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

    FillSingleTau_(L1TkTaus_VtxIsoLoose, hVtxIsoLoose_Rate  , hVtxIsoLoose_Eff);
    FillSingleTau_(L1TkTaus_VtxIsoLoose, hVtxIsoLoose_Rate_C, hVtxIsoLoose_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_VtxIsoLoose, hVtxIsoLoose_Rate_I, hVtxIsoLoose_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_VtxIsoLoose, hVtxIsoLoose_Rate_F, hVtxIsoLoose_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkTaus_VtxIsoTight, hVtxIsoTight_Rate  , hVtxIsoTight_Eff);
    FillSingleTau_(L1TkTaus_VtxIsoTight, hVtxIsoTight_Rate_C, hVtxIsoTight_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_VtxIsoTight, hVtxIsoTight_Rate_I, hVtxIsoTight_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_VtxIsoTight, hVtxIsoTight_Rate_F, hVtxIsoTight_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkTaus_RelIsoLoose, hRelIsoLoose_Rate  , hRelIsoLoose_Eff);
    FillSingleTau_(L1TkTaus_RelIsoLoose, hRelIsoLoose_Rate_C, hRelIsoLoose_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_RelIsoLoose, hRelIsoLoose_Rate_I, hRelIsoLoose_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_RelIsoLoose, hRelIsoLoose_Rate_F, hRelIsoLoose_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkTaus_RelIsoTight, hRelIsoTight_Rate  , hRelIsoTight_Eff);
    FillSingleTau_(L1TkTaus_RelIsoTight, hRelIsoTight_Rate_C, hRelIsoTight_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkTaus_RelIsoTight, hRelIsoTight_Rate_I, hRelIsoTight_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkTaus_RelIsoTight, hRelIsoTight_Rate_F, hRelIsoTight_Eff_F, 1.6, 3.0); // 2.5 is max

    ////////////////////////////////////////////////
    // DiTau
    ////////////////////////////////////////////////
    FillDiTau_(L1TkTaus_Tk, L1TkTaus_VtxIso     , hDiTau_Rate_Tk_VtxIso     , hDiTau_Eff_Tk_VtxIso );
    FillDiTau_(L1TkTaus_Tk, L1TkTaus_RelIso     , hDiTau_Rate_Tk_RelIso     , hDiTau_Eff_Tk_RelIso );
    FillDiTau_(L1TkTaus_Tk, L1TkTaus_VtxIsoLoose, hDiTau_Rate_Tk_VtxIsoLoose, hDiTau_Eff_Tk_VtxIsoLoose );
    FillDiTau_(L1TkTaus_Tk, L1TkTaus_VtxIsoTight, hDiTau_Rate_Tk_VtxIsoTight, hDiTau_Eff_Tk_VtxIsoTight );
    FillDiTau_(L1TkTaus_Tk, L1TkTaus_RelIsoLoose, hDiTau_Rate_Tk_RelIsoLoose, hDiTau_Eff_Tk_RelIsoLoose );
    FillDiTau_(L1TkTaus_Tk, L1TkTaus_RelIsoTight, hDiTau_Rate_Tk_RelIsoTight, hDiTau_Eff_Tk_RelIsoTight );

    ////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_Tk);
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_VtxIso); 
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_RelIso);
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_VtxIsoLoose);
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_VtxIsoTight);
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_RelIsoLoose);
    ApplyDiTauZMatching(seedTk_Collection, L1TkTaus_RelIsoTight);

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

    FillDiTau_(L1TkTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose  , hDiTau_Eff_VtxIsoLoose);
    FillDiTau_(L1TkTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose_C, hDiTau_Eff_VtxIsoLoose_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose_I, hDiTau_Eff_VtxIsoLoose_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose_F, hDiTau_Eff_VtxIsoLoose_F, 1.6, 3.0);

    FillDiTau_(L1TkTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight  , hDiTau_Eff_VtxIsoTight);
    FillDiTau_(L1TkTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight_C, hDiTau_Eff_VtxIsoTight_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight_I, hDiTau_Eff_VtxIsoTight_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight_F, hDiTau_Eff_VtxIsoTight_F, 1.6, 3.0);

    FillDiTau_(L1TkTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose  , hDiTau_Eff_RelIsoLoose);
    FillDiTau_(L1TkTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose_C, hDiTau_Eff_RelIsoLoose_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose_I, hDiTau_Eff_RelIsoLoose_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose_F, hDiTau_Eff_RelIsoLoose_F, 1.6, 3.0);

    FillDiTau_(L1TkTaus_RelIsoTight, hDiTau_Rate_RelIsoTight  , hDiTau_Eff_RelIsoTight);
    FillDiTau_(L1TkTaus_RelIsoTight, hDiTau_Rate_RelIsoTight_C, hDiTau_Eff_RelIsoTight_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_RelIsoTight, hDiTau_Rate_RelIsoTight_I, hDiTau_Eff_RelIsoTight_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_RelIsoTight, hDiTau_Rate_RelIsoTight_F, hDiTau_Eff_RelIsoTight_F, 1.6, 3.0);

    // Progress bar
    if (!DEBUG) auxTools_.ProgressBar(jentry, nEntries, 100, 100);
    
  }// For-loop: Entries

  // Fill counters
  hCounters->SetBinContent( 1, nAllEvts);
  hCounters->SetBinContent( 2, nEvts);
  hCounters->SetBinContent( 3, nEvtsSeedPt);
  hCounters->SetBinContent( 4, nEvtsSeedEta);
  hCounters->SetBinContent( 5, nEvtsSeedChiSq);
  hCounters->SetBinContent( 6, nEvtsSeedStubs);
  hCounters->SetBinContent( 7, nEvtsVtxIso);
  hCounters->SetBinContent( 8, nEvtsRelIso);
  hCounters->SetBinContent( 9, nEvtsVtxIsoLoose);
  hCounters->SetBinContent(10, nEvtsVtxIsoTight);
  hCounters->SetBinContent(11, nEvtsRelIsoLoose);
  hCounters->SetBinContent(12, nEvtsRelIsoTight);
  hCounters->SetBinContent(13, 0);
  hCounters->SetBinContent(14, nEvtsMcMatch);
  hCounters->GetXaxis()->SetBinLabel( 1, "All Evts");
  hCounters->GetXaxis()->SetBinLabel( 2, "Evts");
  hCounters->GetXaxis()->SetBinLabel( 3, "Seed Pt");
  hCounters->GetXaxis()->SetBinLabel( 4, "Seed Eta");
  hCounters->GetXaxis()->SetBinLabel( 5, "Seed ChiSq");
  hCounters->GetXaxis()->SetBinLabel( 6, "Seed Stubs");
  hCounters->GetXaxis()->SetBinLabel( 7, "VtxIso");
  hCounters->GetXaxis()->SetBinLabel( 8, "RelIso");
  hCounters->GetXaxis()->SetBinLabel( 9, "VtxIso (L)");
  hCounters->GetXaxis()->SetBinLabel(10, "VtxIso (T)");
  hCounters->GetXaxis()->SetBinLabel(11, "RelIso (L)");
  hCounters->GetXaxis()->SetBinLabel(12, "RelIso (T)");
  hCounters->GetXaxis()->SetBinLabel(13, "");
  hCounters->GetXaxis()->SetBinLabel(14, "Matched");
  
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

  histoTools_.ConvertToRateHisto_1D(hVtxIsoLoose_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hVtxIsoLoose_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hVtxIsoLoose_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hVtxIsoLoose_Rate_F, N);

  histoTools_.ConvertToRateHisto_1D(hVtxIsoTight_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hVtxIsoTight_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hVtxIsoTight_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hVtxIsoTight_Rate_F, N);

  histoTools_.ConvertToRateHisto_1D(hRelIsoLoose_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hRelIsoLoose_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hRelIsoLoose_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hRelIsoLoose_Rate_F, N);

  histoTools_.ConvertToRateHisto_1D(hRelIsoTight_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hRelIsoTight_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hRelIsoTight_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hRelIsoTight_Rate_F, N);

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

  FinaliseEffHisto_( hVtxIsoLoose_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIsoLoose_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIsoLoose_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIsoLoose_Eff_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hVtxIsoTight_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIsoTight_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIsoTight_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hVtxIsoTight_Eff_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hRelIsoLoose_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIsoLoose_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIsoLoose_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIsoLoose_Eff_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hRelIsoTight_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIsoTight_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIsoTight_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hRelIsoTight_Eff_F, nEvtsWithMaxHTaus);

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

  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoLoose  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoLoose_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoLoose_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoLoose_F, N);

  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoTight  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoTight_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoTight_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_VtxIsoTight_F, N);

  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoLoose  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoLoose_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoLoose_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoLoose_F, N);

  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoTight  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoTight_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoTight_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_RelIsoTight_F, N);
  
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

  FinaliseEffHisto_( hDiTau_Eff_VtxIsoLoose  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIsoLoose_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIsoLoose_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIsoLoose_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hDiTau_Eff_VtxIsoTight  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIsoTight_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIsoTight_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_VtxIsoTight_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hDiTau_Eff_RelIsoLoose  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIsoLoose_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIsoLoose_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIsoLoose_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_( hDiTau_Eff_RelIsoTight  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIsoTight_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIsoTight_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_RelIsoTight_F, nEvtsWithMaxHTaus);

  // DiTau (Tk-Other)
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_VtxIso, N);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_RelIso, N);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_VtxIsoLoose, N);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_VtxIsoTight, N);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_RelIsoLoose, N);
  histoTools_.ConvertToRateHisto_2D(hDiTau_Rate_Tk_RelIsoTight, N);

  FinaliseEffHisto_( hDiTau_Eff_Tk_VtxIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_RelIso, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_VtxIsoLoose, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_VtxIsoTight, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_RelIsoLoose, nEvtsWithMaxHTaus);;
  FinaliseEffHisto_( hDiTau_Eff_Tk_RelIsoTight, nEvtsWithMaxHTaus);;

  // Turn-Ons 
  // TEfficiency *pEff = 0;
  // pEff = new TEfficiency(*hCalo_TurnOn50_passed, *hMcHadronicTau_VisEt);
  // hCalo_TurnOn50 = (TH1D*) pEff->Clone();
  histoTools_.DivideHistos_1D(hTk_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hTk_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hTk_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hTk_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hRelIso_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hTk_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTk_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hTk_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hTk_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hTk_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hVtxIso_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hRelIso_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hRelIso_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hVtxIsoLoose_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hVtxIsoTight_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hRelIsoLoose_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hRelIsoTight_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);


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
  const char* tIso  = ";vertex isolation;relative isolation";
  const char* tChi  = ";#chi^{2};Entries / %.2f";
  const char* tRChi = ";#chi^{2}_{#nu};Entries / %.2f";
  const char* tW    = ";w_{#tau};Entries / %.2f";
  const char* tG    = ";#gamma;Entries / %.2f";
  const char* tGW    = ";#gamma;w_{#tau}";
  // const char* tRate = ";#E_{T} (GeV);Rate (kHz) / %.0f GeV";

  // GenParticles Histograms
  histoTools_.BookHisto_2D(hGenP_VisEt_Vs_dRMaxLdgPion, "GenP_VisEt_Vs_dRMaxLdgPion", ";#DeltaR_{max}(#pi_{ldg}^{#pm},#pi^{#pm});E_{T}^{vis}",  50,  0.0, +0.25, 100, 0.0, +200.0);
  histoTools_.BookHisto_2D(hGenP_PtLdg_Vs_dRMaxLdgPion, "GenP_PtLdg_Vs_dRMaxLdgPion", ";#DeltaR_{max}(#pi_{ldg}^{#pm},#pi^{#pm});p_{T}^{#pi_{ldg}^{#pm}}",  50,  0.0, +0.25, 100, 0.0, +200.0);

  // Counters
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 15, 0.0, +15.0);

  // L1TkTaus
  histoTools_.BookHisto_1D(hL1TkTau_Multiplicity , "L1TkTau_Multiplicity" , tN   , nN   , minN   , maxN   );
  histoTools_.BookHisto_1D(hL1TkTau_Multiplicity_MC, "L1TkTau_Multiplicity_MC", tN, nN, minN, maxN );
  histoTools_.BookHisto_1D(hL1TkTau_JetWidth     , "L1TkTau_JetWidth"     , tW   , nM   , minM   , maxM   );
  histoTools_.BookHisto_1D(hL1TkTau_DonutRatio   , "L1TkTau_DonutRatio"   , tG   , nG   , minG   , maxG   );
  histoTools_.BookHisto_2D(hL1TkTau_DonutRatio_Vs_JetWidth, "L1TkTau_DonutRatio_Vs_JetWidth", tGW, nG, minG, maxG, nM, minM, maxM);
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
  histoTools_.BookHisto_2D(hL1TkTau_VtxIso_Vs_RelIso, "L1TkTau_VtxIso_Vs_RelIso", tIso, 200, 0.0, 10.0, 200, 0.0, 10.0);
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
  histoTools_.BookHisto_2D(hL1TkIsoTau_DonutRatio_Vs_JetWidth, "L1TkIsoTau_DonutRatio_Vs_JetWidth", tGW, nG, minG, maxG, nM, minM, maxM   );
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
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt    , "L1TkIsoTau_ResolutionEt"     , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt_1pr, "L1TkIsoTau_ResolutionEt_1pr" , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt_3pr, "L1TkIsoTau_ResolutionEt_3pr" , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt_withNeutrals, "L1TkIsoTau_ResolutionEt_withNeutrals", ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 100,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEt_noNeutrals  , "L1TkIsoTau_ResolutionEt_noNeutrals"  , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 100,  -5.0,  +5.0);

  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta    , "L1TkIsoTau_ResolutionEta"    , ";#delta #eta / #eta^{vis};Events / %.2f", 200,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta_1pr, "L1TkIsoTau_ResolutionEta_1pr", ";#delta #eta / #eta^{vis};Events / %.2f", 200,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta_3pr, "L1TkIsoTau_ResolutionEta_3pr", ";#delta #eta / #eta^{vis};Events / %.2f", 200,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta_withNeutrals, "L1TkIsoTau_ResolutionEta_withNeutrals", ";#delta #eta / #eta^{vis};Events / %.2f", 200,  -5.0,  +5.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionEta_noNeutrals  , "L1TkIsoTau_ResolutionEta_noNeutrals"  , ";#delta #eta / #eta^{vis};Events / %.2f", 200,  -5.0,  +5.0);

  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi    , "L1TkIsoTau_ResolutionPhi"    , ";#delta #phi / #phi^{vis};Events / %.2f rads", 200,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi_1pr, "L1TkIsoTau_ResolutionPhi_1pr", ";#delta #phi / #phi^{vis};Events / %.2f rads", 200,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi_3pr, "L1TkIsoTau_ResolutionPhi_3pr", ";#delta #phi / #phi^{vis};Events / %.2f rads", 200,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi_withNeutrals, "L1TkIsoTau_ResolutionPhi_withNeutrals", ";#phi (rads);Events / %.2f rads", 200,  -10.0,  +10.0);
  histoTools_.BookHisto_1D(hL1TkIsoTau_ResolutionPhi_noNeutrals  , "L1TkIsoTau_ResolutionPhi_noNeutrals"  , ";#phi (rads);Events / %.2f rads", 200,  -10.0,  +10.0);

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
  histoTools_.BookHisto_1D(hVtxIsoLoose_Rate  , "VtxIsoLoose_Rate"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_Rate_C, "VtxIsoLoose_Rate_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_Rate_I, "VtxIsoLoose_Rate_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_Rate_F, "VtxIsoLoose_Rate_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Rate  , "VtxIsoTight_Rate"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Rate_C, "VtxIsoTight_Rate_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Rate_I, "VtxIsoTight_Rate_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Rate_F, "VtxIsoTight_Rate_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Rate  , "RelIsoLoose_Rate"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Rate_C, "RelIsoLoose_Rate_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Rate_I, "RelIsoLoose_Rate_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Rate_F, "RelIsoLoose_Rate_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Rate  , "RelIsoTight_Rate"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Rate_C, "RelIsoTight_Rate_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Rate_I, "RelIsoTight_Rate_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Rate_F, "RelIsoTight_Rate_F", "", nEt , minEt , maxEt );
  
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
  histoTools_.BookHisto_1D(hVtxIsoLoose_Eff  , "VtxIsoLoose_Eff"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_Eff_C, "VtxIsoLoose_Eff_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_Eff_I, "VtxIsoLoose_Eff_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_Eff_F, "VtxIsoLoose_Eff_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Eff  , "VtxIsoTight_Eff"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Eff_C, "VtxIsoTight_Eff_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Eff_I, "VtxIsoTight_Eff_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_Eff_F, "VtxIsoTight_Eff_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Eff  , "RelIsoLoose_Eff"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Eff_C, "RelIsoLoose_Eff_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Eff_I, "RelIsoLoose_Eff_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_Eff_F, "RelIsoLoose_Eff_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Eff  , "RelIsoTight_Eff"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Eff_C, "RelIsoTight_Eff_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Eff_I, "RelIsoTight_Eff_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_Eff_F, "RelIsoTight_Eff_F", "", nEt , minEt , maxEt );
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
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoLoose  , "DiTau_Rate_VtxIsoLoose"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoLoose_C, "DiTau_Rate_VtxIsoLoose_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoLoose_I, "DiTau_Rate_VtxIsoLoose_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoLoose_F, "DiTau_Rate_VtxIsoLoose_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoTight  , "DiTau_Rate_VtxIsoTight"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoTight_C, "DiTau_Rate_VtxIsoTight_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoTight_I, "DiTau_Rate_VtxIsoTight_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_VtxIsoTight_F, "DiTau_Rate_VtxIsoTight_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoLoose  , "DiTau_Rate_RelIsoLoose"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoLoose_C, "DiTau_Rate_RelIsoLoose_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoLoose_I, "DiTau_Rate_RelIsoLoose_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoLoose_F, "DiTau_Rate_RelIsoLoose_F", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoTight  , "DiTau_Rate_RelIsoTight"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoTight_C, "DiTau_Rate_RelIsoTight_C", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoTight_I, "DiTau_Rate_RelIsoTight_I", "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_RelIsoTight_F, "DiTau_Rate_RelIsoTight_F", "", nEt , minEt , maxEt );

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
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoLoose  , "DiTau_Eff_VtxIsoLoose"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoLoose_C, "DiTau_Eff_VtxIsoLoose_C" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoLoose_I, "DiTau_Eff_VtxIsoLoose_I" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoLoose_F, "DiTau_Eff_VtxIsoLoose_F" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoTight  , "DiTau_Eff_VtxIsoTight"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoTight_C, "DiTau_Eff_VtxIsoTight_C" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoTight_I, "DiTau_Eff_VtxIsoTight_I" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_VtxIsoTight_F, "DiTau_Eff_VtxIsoTight_F" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoLoose  , "DiTau_Eff_RelIsoLoose"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoLoose_C, "DiTau_Eff_RelIsoLoose_C" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoLoose_I, "DiTau_Eff_RelIsoLoose_I" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoLoose_F, "DiTau_Eff_RelIsoLoose_F" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoTight  , "DiTau_Eff_RelIsoTight"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoTight_C, "DiTau_Eff_RelIsoTight_C" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoTight_I, "DiTau_Eff_RelIsoTight_I" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_RelIsoTight_F, "DiTau_Eff_RelIsoTight_F" , "", nEt , minEt , maxEt );
  
  // Turn-Ons
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt    , "McHadronicTau_VisEt"    , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_1pr, "McHadronicTau_VisEt_1pr", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_3pr, "McHadronicTau_VisEt_3pr", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_withNeutrals, "McHadronicTau_VisEt_withNeutrals", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_noNeutrals, "McHadronicTau_VisEt_noNeutrals", "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hTk_TurnOn25, "Tk_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn25_1pr, "Tk_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn25_3pr, "Tk_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn25_withNeutrals, "Tk_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn25_noNeutrals, "Tk_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hVtxIso_TurnOn25, "VtxIso_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25_1pr, "VtxIso_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25_3pr, "VtxIso_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25_withNeutrals, "VtxIso_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn25_noNeutrals, "VtxIso_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hRelIso_TurnOn25, "RelIso_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn25_1pr, "RelIso_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn25_3pr, "RelIso_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn25_withNeutrals, "RelIso_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn25_noNeutrals, "RelIso_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn25, "VtxIsoLoose_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn25_1pr, "VtxIsoLoose_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn25_3pr, "VtxIsoLoose_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn25_withNeutrals, "VtxIsoLoose_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn25_noNeutrals, "VtxIsoLoose_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn25, "VtxIsoTight_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn25_1pr, "VtxIsoTight_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn25_3pr, "VtxIsoTight_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn25_withNeutrals, "VtxIsoTight_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn25_noNeutrals, "VtxIsoTight_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn25, "RelIsoLoose_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn25_1pr, "RelIsoLoose_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn25_3pr, "RelIsoLoose_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn25_withNeutrals, "RelIsoLoose_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn25_noNeutrals, "RelIsoLoose_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn25, "RelIsoTight_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn25_1pr, "RelIsoTight_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn25_3pr, "RelIsoTight_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn25_withNeutrals, "RelIsoTight_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn25_noNeutrals, "RelIsoTight_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hTk_TurnOn50, "Tk_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn50_1pr, "Tk_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn50_3pr, "Tk_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn50_withNeutrals, "Tk_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTk_TurnOn50_noNeutrals, "Tk_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hVtxIso_TurnOn50, "VtxIso_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50_1pr, "VtxIso_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50_3pr, "VtxIso_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50_withNeutrals, "VtxIso_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIso_TurnOn50_noNeutrals, "VtxIso_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hRelIso_TurnOn50, "RelIso_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn50_1pr, "RelIso_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn50_3pr, "RelIso_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn50_withNeutrals, "RelIso_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIso_TurnOn50_noNeutrals, "RelIso_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn50, "VtxIsoLoose_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn50_1pr, "VtxIsoLoose_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn50_3pr, "VtxIsoLoose_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn50_withNeutrals, "VtxIsoLoose_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoLoose_TurnOn50_noNeutrals, "VtxIsoLoose_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn50, "VtxIsoTight_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn50_1pr, "VtxIsoTight_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn50_3pr, "VtxIsoTight_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn50_withNeutrals, "VtxIsoTight_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hVtxIsoTight_TurnOn50_noNeutrals, "VtxIsoTight_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn50, "RelIsoLoose_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn50_1pr, "RelIsoLoose_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn50_3pr, "RelIsoLoose_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn50_withNeutrals, "RelIsoLoose_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoLoose_TurnOn50_noNeutrals, "RelIsoLoose_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn50, "RelIsoTight_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn50_1pr, "RelIsoTight_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn50_3pr, "RelIsoTight_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn50_withNeutrals, "RelIsoTight_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hRelIsoTight_TurnOn50_noNeutrals, "RelIsoTight_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

  // DiTau (Tk-Other)
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_VtxIso, "DiTau_Rate_Tk_VtxIso", "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_RelIso, "DiTau_Rate_Tk_RelIso", "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_VtxIsoLoose, "DiTau_Rate_Tk_VtxIsoLoose"   , "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_VtxIsoTight, "DiTau_Rate_Tk_VtxIsoTight"   , "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_RelIsoLoose, "DiTau_Rate_Tk_RelIsoLoose"   , "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Rate_Tk_RelIsoTight, "DiTau_Rate_Tk_RelIsoTight"   , "", nEt, minEt, maxEt, nEt, minEt, maxEt);

  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_VtxIso , "DiTau_Eff_Tk_VtxIso" , "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_RelIso , "DiTau_Eff_Tk_RelIso" , "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_VtxIsoLoose, "DiTau_Eff_Tk_VtxIsoLoose", "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_VtxIsoTight, "DiTau_Eff_Tk_VtxIsoTight", "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_RelIsoLoose, "DiTau_Eff_Tk_RelIsoLoose", "", nEt, minEt, maxEt, nEt, minEt, maxEt);
  histoTools_.BookHisto_2D(hDiTau_Eff_Tk_RelIsoTight, "DiTau_Eff_Tk_RelIsoTight", "", nEt, minEt, maxEt, nEt, minEt, maxEt);

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
  hL1TkTau_DonutRatio_Vs_JetWidth->Write();
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
  hL1TkIsoTau_DonutRatio_Vs_JetWidth->Write();
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
  hL1TkIsoTau_ResolutionEt_1pr->Write();
  hL1TkIsoTau_ResolutionEt_3pr->Write();
  hL1TkIsoTau_ResolutionEt_withNeutrals->Write();
  hL1TkIsoTau_ResolutionEt_noNeutrals->Write();

  hL1TkIsoTau_ResolutionEta->Write();
  hL1TkIsoTau_ResolutionEta_1pr->Write();
  hL1TkIsoTau_ResolutionEta_3pr->Write();
  hL1TkIsoTau_ResolutionEta_withNeutrals->Write();
  hL1TkIsoTau_ResolutionEta_noNeutrals->Write();

  hL1TkIsoTau_ResolutionPhi->Write();
  hL1TkIsoTau_ResolutionPhi_1pr->Write();
  hL1TkIsoTau_ResolutionPhi_3pr->Write();
  hL1TkIsoTau_ResolutionPhi_withNeutrals->Write();
  hL1TkIsoTau_ResolutionPhi_noNeutrals->Write();

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
  hVtxIsoLoose_Rate->Write();
  hVtxIsoLoose_Rate_C->Write();
  hVtxIsoLoose_Rate_I->Write();
  hVtxIsoLoose_Rate_F->Write();
  hVtxIsoTight_Rate->Write();
  hVtxIsoTight_Rate_C->Write();
  hVtxIsoTight_Rate_I->Write();
  hVtxIsoTight_Rate_F->Write();
  hRelIsoLoose_Rate->Write();
  hRelIsoLoose_Rate_C->Write();
  hRelIsoLoose_Rate_I->Write();
  hRelIsoLoose_Rate_F->Write();
  hRelIsoTight_Rate->Write();
  hRelIsoTight_Rate_C->Write();
  hRelIsoTight_Rate_I->Write();
  hRelIsoTight_Rate_F->Write();

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
  hVtxIsoLoose_Eff->Write();
  hVtxIsoLoose_Eff_C->Write();
  hVtxIsoLoose_Eff_I->Write();
  hVtxIsoLoose_Eff_F->Write();
  hVtxIsoTight_Eff->Write();
  hVtxIsoTight_Eff_C->Write();
  hVtxIsoTight_Eff_I->Write();
  hVtxIsoTight_Eff_F->Write();
  hRelIsoLoose_Eff->Write();
  hRelIsoLoose_Eff_C->Write();
  hRelIsoLoose_Eff_I->Write();
  hRelIsoLoose_Eff_F->Write();
  hRelIsoTight_Eff->Write();
  hRelIsoTight_Eff_C->Write();
  hRelIsoTight_Eff_I->Write();
  hRelIsoTight_Eff_F->Write();

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
  hDiTau_Rate_VtxIsoLoose->Write();
  hDiTau_Rate_VtxIsoLoose_C->Write();
  hDiTau_Rate_VtxIsoLoose_I->Write();
  hDiTau_Rate_VtxIsoLoose_F->Write();
  hDiTau_Rate_VtxIsoTight->Write();
  hDiTau_Rate_VtxIsoTight_C->Write();
  hDiTau_Rate_VtxIsoTight_I->Write();
  hDiTau_Rate_VtxIsoTight_F->Write();
  hDiTau_Rate_RelIsoLoose->Write();
  hDiTau_Rate_RelIsoLoose_C->Write();
  hDiTau_Rate_RelIsoLoose_I->Write();
  hDiTau_Rate_RelIsoLoose_F->Write();
  hDiTau_Rate_RelIsoTight->Write();
  hDiTau_Rate_RelIsoTight_C->Write();
  hDiTau_Rate_RelIsoTight_I->Write();
  hDiTau_Rate_RelIsoTight_F->Write();

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
  hDiTau_Eff_VtxIsoLoose->Write();
  hDiTau_Eff_VtxIsoLoose_C->Write();
  hDiTau_Eff_VtxIsoLoose_I->Write();
  hDiTau_Eff_VtxIsoLoose_F->Write();
  hDiTau_Eff_VtxIsoTight->Write();
  hDiTau_Eff_VtxIsoTight_C->Write();
  hDiTau_Eff_VtxIsoTight_I->Write();
  hDiTau_Eff_VtxIsoTight_F->Write();
  hDiTau_Eff_RelIsoLoose->Write();
  hDiTau_Eff_RelIsoLoose_C->Write();
  hDiTau_Eff_RelIsoLoose_I->Write();
  hDiTau_Eff_RelIsoLoose_F->Write();
  hDiTau_Eff_RelIsoTight->Write();
  hDiTau_Eff_RelIsoTight_C->Write();
  hDiTau_Eff_RelIsoTight_I->Write();
  hDiTau_Eff_RelIsoTight_F->Write();

  // DiTau (Tk-Other)
  hDiTau_Rate_Tk_VtxIso->Write();
  hDiTau_Rate_Tk_RelIso->Write();
  hDiTau_Rate_Tk_VtxIsoLoose->Write();
  hDiTau_Rate_Tk_VtxIsoTight->Write();
  hDiTau_Rate_Tk_RelIsoLoose->Write();
  hDiTau_Rate_Tk_RelIsoTight->Write();

  hDiTau_Eff_Tk_VtxIso->Write();
  hDiTau_Eff_Tk_RelIso->Write();
  hDiTau_Eff_Tk_VtxIsoLoose->Write();
  hDiTau_Eff_Tk_VtxIsoTight->Write();
  hDiTau_Eff_Tk_RelIsoLoose->Write();
  hDiTau_Eff_Tk_RelIsoTight->Write();

  // Turn-Ons
  hMcHadronicTau_VisEt->Write();
  hMcHadronicTau_VisEt_1pr->Write();
  hMcHadronicTau_VisEt_3pr->Write();
  hMcHadronicTau_VisEt_withNeutrals->Write();
  hMcHadronicTau_VisEt_noNeutrals->Write();

  hTk_TurnOn25->Write();
  hTk_TurnOn25_1pr->Write();
  hTk_TurnOn25_3pr->Write();
  hTk_TurnOn25_withNeutrals->Write();
  hTk_TurnOn25_noNeutrals->Write();
  hVtxIso_TurnOn25->Write();
  hVtxIso_TurnOn25_1pr->Write();
  hVtxIso_TurnOn25_3pr->Write();
  hVtxIso_TurnOn25_withNeutrals->Write();
  hVtxIso_TurnOn25_noNeutrals->Write();
  hRelIso_TurnOn25->Write();
  hRelIso_TurnOn25_1pr->Write();
  hRelIso_TurnOn25_3pr->Write();
  hRelIso_TurnOn25_withNeutrals->Write();
  hRelIso_TurnOn25_noNeutrals->Write();
  hVtxIsoLoose_TurnOn25->Write();
  hVtxIsoLoose_TurnOn25_1pr->Write();
  hVtxIsoLoose_TurnOn25_3pr->Write();
  hVtxIsoLoose_TurnOn25_withNeutrals->Write();
  hVtxIsoLoose_TurnOn25_noNeutrals->Write();
  hVtxIsoTight_TurnOn25->Write();
  hVtxIsoTight_TurnOn25_1pr->Write();
  hVtxIsoTight_TurnOn25_3pr->Write();
  hVtxIsoTight_TurnOn25_withNeutrals->Write();
  hVtxIsoTight_TurnOn25_noNeutrals->Write();
  hRelIsoLoose_TurnOn25->Write();
  hRelIsoLoose_TurnOn25_1pr->Write();
  hRelIsoLoose_TurnOn25_3pr->Write();
  hRelIsoLoose_TurnOn25_withNeutrals->Write();
  hRelIsoLoose_TurnOn25_noNeutrals->Write();
  hRelIsoTight_TurnOn25->Write();
  hRelIsoTight_TurnOn25_1pr->Write();
  hRelIsoTight_TurnOn25_3pr->Write();
  hRelIsoTight_TurnOn25_withNeutrals->Write();
  hRelIsoTight_TurnOn25_noNeutrals->Write();

  hTk_TurnOn50->Write();
  hTk_TurnOn50_1pr->Write();
  hTk_TurnOn50_3pr->Write();
  hTk_TurnOn50_withNeutrals->Write();
  hTk_TurnOn50_noNeutrals->Write();
  hVtxIso_TurnOn50->Write();
  hVtxIso_TurnOn50_1pr->Write();
  hVtxIso_TurnOn50_3pr->Write();
  hVtxIso_TurnOn50_withNeutrals->Write();
  hVtxIso_TurnOn50_noNeutrals->Write();
  hRelIso_TurnOn50->Write();
  hRelIso_TurnOn50_1pr->Write();
  hRelIso_TurnOn50_3pr->Write();
  hRelIso_TurnOn50_withNeutrals->Write();
  hRelIso_TurnOn50_noNeutrals->Write();
  hVtxIsoLoose_TurnOn50->Write();
  hVtxIsoLoose_TurnOn50_1pr->Write();
  hVtxIsoLoose_TurnOn50_3pr->Write();
  hVtxIsoLoose_TurnOn50_withNeutrals->Write();
  hVtxIsoLoose_TurnOn50_noNeutrals->Write();
  hVtxIsoTight_TurnOn50->Write();
  hVtxIsoTight_TurnOn50_1pr->Write();
  hVtxIsoTight_TurnOn50_3pr->Write();
  hVtxIsoTight_TurnOn50_withNeutrals->Write();
  hVtxIsoTight_TurnOn50_noNeutrals->Write();
  hRelIsoLoose_TurnOn50->Write();
  hRelIsoLoose_TurnOn50_1pr->Write();
  hRelIsoLoose_TurnOn50_3pr->Write();
  hRelIsoLoose_TurnOn50_withNeutrals->Write();
  hRelIsoLoose_TurnOn50_noNeutrals->Write();
  hRelIsoTight_TurnOn50->Write();
  hRelIsoTight_TurnOn50_1pr->Write();
  hRelIsoTight_TurnOn50_3pr->Write();
  hRelIsoTight_TurnOn50_withNeutrals->Write();
  hRelIsoTight_TurnOn50_noNeutrals->Write();


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
    auxTools_.Efficiency(nPass, nEvtsTotal, "binomial", eff, err );

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
  TLorentzVector sigTks_p4 = L1TkTaus.at(0).GetSigConeTTTracksP4(); // fixme (sort with pT)?
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
  TLorentzVector sigTks_p4_mc = L1TkTaus_mcMatched.at(0).GetSigConeTTTracksP4(); // fixme (sort with pT)?
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
  TLorentzVector sigTks_p4 = L1TkTaus.at(1).GetSigConeTTTracksP4(); // fixme (sort with pT)?
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
  TLorentzVector sigTks_p4_mc = L1TkTaus_mcMatched.at(1).GetSigConeTTTracksP4(); // fixme (sort with pT)?
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
  TLorentzVector sigTks1_p4 = L1TkTaus1.at(0).GetSigConeTTTracksP4(); // fixme (sort with pT)?
  TLorentzVector sigTks2_p4 = L1TkTaus2.at(0).GetSigConeTTTracksP4(); // fixme (sort with pT)?
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
  TLorentzVector sigTks1_p4_mc = L1TkTaus1_mcMatched.at(0).GetSigConeTTTracksP4(); // fixme (sort with pT)?
  TLorentzVector sigTks2_p4_mc = L1TkTaus2_mcMatched.at(0).GetSigConeTTTracksP4(); // fixme (sort with pT)?
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
				   TH1D *hTurnOn,
				   TH1D *hTurnOn_1pr,
				   TH1D *hTurnOn_3pr,
				   TH1D *hTurnOn_withNeutrals,
				   TH1D *hTurnOn_noNeutrals)
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
      hTurnOn->Fill( p.p4vis().Et() ); // turn-on fill: warning if MC matchind dR is too lose this can have more entries than denominator

	
      if (p.finalDaughtersNeutral().size() > 0)
	  {
	    hTurnOn_withNeutrals->Fill( p.p4vis().Et() );
	  }
      else
	  {
	    hTurnOn_noNeutrals->Fill( p.p4vis().Et() );
	  }
      if (p.finalDaughtersCharged().size() == 1) 
	  {
	    hTurnOn_1pr->Fill( p.p4vis().Et() );	    
	  }
      else if (p.finalDaughtersCharged().size() == 3) 
	{
	  hTurnOn_3pr->Fill( p.p4vis().Et() );
	}
      
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
			      vector<TTTrack> isoTks,
			      double isoConeTks_dPOCAz)
//============================================================================
{
  if (!L1TkTau.HasMatchingTk()) return; 
  vector<TTTrack> isoConeTks_tmp;
  vector<TTTrack> isoConeTks;
  vector<TTTrack> isoAnnulusTks;
  vector<TTTrack> sigConeTks = L1TkTau.GetSigConeTTTracks();
  TTTrack seedTk = L1TkTau.GetMatchingTk();

  // For-loop: All Tracks
  for (vector<TTTrack>::iterator tk = isoTks.begin(); tk != isoTks.end(); tk++)
    {
      double dR = auxTools_.DeltaR(tk->getEta(), tk->getPhi(), L1TkTau.GetMatchingTk().getEta(), L1TkTau.GetMatchingTk().getPhi());
      double dPOCAz = abs(seedTk.getZ0() - tk->getZ0());

      // Only consider tracks within singal cone
      bool bInSigCone    = (dR < L1TkTau.GetIsoConeMin());
      bool bInIsoCone    = (dR < L1TkTau.GetIsoConeMax());
      bool bInIsoAnnulus = bInIsoCone*(!bInSigCone); 

      // Require isolation tracks come from same vertex as the seed track (NEW)
      if (dPOCAz > isoConeTks_dPOCAz) continue;
            
      // Save the track as an isolation-cone track
      if (bInIsoCone) isoConeTks_tmp.push_back(*tk);
      if (bInIsoAnnulus) isoAnnulusTks.push_back(*tk);
    }
  
  // For-loop: All tracks in isolation cone
  for (vector<TTTrack>::iterator iTk = isoConeTks_tmp.begin(); iTk != isoConeTks_tmp.end(); iTk++)
    {
      bool bSaveTk = true;

      // For-loop: All tracks in signal cone
      for (vector<TTTrack>::iterator sTk = sigConeTks.begin(); sTk != sigConeTks.end(); sTk++)
	{
      
	  double dPt = abs(sTk->getPt() - iTk->getPt());
	  double dR  = auxTools_.DeltaR(sTk->getEta(), sTk->getPhi(), iTk->getEta(), iTk->getPhi());
	  double dZ  = abs(sTk->getZ0() - iTk->getZ0());
	  
	  // Isolation tracks matches a signal cone track that was clustered
	  if ((dPt == 0.0) && (dR == 0.0) && (dZ == 0.0) )
	    {
	      bSaveTk = false;
	      break;
	    }
	}
      // Now only add isolation tracks which are not signal clustered tracks
      if (bSaveTk) isoConeTks.push_back(*iTk);
    }

  if (0)
    {
      std::cout << "\nsigConeTks.size() = " << sigConeTks.size() << std::endl;
      std::cout << "isoAnnulusTks.size() = " << isoAnnulusTks.size() << std::endl;
      std::cout << "isoConeTks.size() = " << isoConeTks.size() << std::endl;
    }

  // Sanity check
  if (isoConeTks.size() < isoAnnulusTks.size())
    {
      cout << "=== TkTaus::GetIsoConeTracks() - Unexpected number of tracks in isolation cone. Cannot be less than those in isolation annulus! EXIT" << endl;
      exit(1);
    }

  // Save the tracks
  L1TkTau.SetIsoConeTracks(isoConeTks);
  L1TkTau.SetIsoAnnulusTracks(isoAnnulusTks);
  
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
  if (signalCone_max > sigCone_dRCutoff) signalCone_max = sigCone_dRCutoff; // fixme
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

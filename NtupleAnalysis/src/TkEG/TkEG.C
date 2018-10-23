#ifndef TkEG_cxx
#define TkEG_cxx

// User
//#include "../Auxiliary/interface/constants.h"
#include "TkEG.h"

// ROOT
//#include "TFitResult.h"
//#include "TF1.h"

// C++
#include <stdexcept>

//============================================================================
void TkEG::InitObjects(void)
//============================================================================
{
  return; 
}


//============================================================================
void TkEG::InitVars_()
//============================================================================
{

  datasets_  = datasets_.GetDataset(mcSample);
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;
  
  cfg_DEBUG = false;
  if (cfg_DEBUG) std::cout << "=== TkEG::InitVars_()" << std::endl;
  
  cfg_AddL1Tks   = true;
  cfg_AddEGs     = true;
  cfg_AddGenP    = true;
  
  // Track parameters
  cfg_tk_Collection  =  "TTTracks"; // Default: "TTTracks" (not "TTPixelTracks")
  cfg_tk_nFitParams  = 4;           // Default: 4
  cfg_tk_minPt       = 2.00;        // Default: 2.0
  cfg_tk_minEta      = 0.0;         // Default: 0.0
  cfg_tk_maxEta      = 2.4;        // Default: 1e6
  cfg_tk_maxChiSq    = 50.0;         // Default: 1e6
  cfg_tk_minStubs    =   5;         // Default: 0

  // TkEG algorithm parameters
  minStubs_trk      = 5; 
  maxChi2_trk       = 50.0; // GeV
  minPt_leadtrk     = 5.0; // GeV
  maxEta_leadtrk    = 2.4;
  minDeltaR_leadtrk = 0.0;
  maxDeltaR_leadtrk = 0.15;//0.3;
  maxDeltaZ_trk     = 0.8;  // cm
  maxInvMass_trk    = 1.4; // GeV 
  minEt_EG          = 2.0;  // GeV
  minDeltaR_EG      = 0.0;
  maxDeltaR_EG      = 0.15;//0.3;
  maxInvMass_EG     = 1.77; // GeV
  maxDeltaR_MCmatch = 0.2;    
  minDeltaR_iso     = 0.15;
  maxDeltaR_iso     = 0.3;
  maxDeltaZ_iso     = 0.6;  // cm
  useRelIso         = true;
  maxRelIso         = 0.10; //0.2
  maxVtxIso         = 0.50;

  // Eta regions
  _eta_C             =   0.8; 
  _eta_F             =   1.6;


  return;

}


//============================================================================
void TkEG::PrintSettings(void)
//============================================================================
{

  // TODO
  
  return;
}


//============================================================================
void TkEG::Loop()
//============================================================================
{
  
  // Sanity check
  if (fChain == 0) return;
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  // Initialisations
  InitVars_();
  BookHistos_();
  Long64_t nbytes              = 0;
  Long64_t nb                  = 0;
  int nEvtsWithMaxHTaus        = 0; 
  unsigned int nEvts           = 0;
  unsigned int nEvtsSeedPt     = 0;
  unsigned int nEvtsSeedEta    = 0;
  unsigned int nEvtsSeedChiSq  = 0;
  unsigned int nEvtsSeedStubs  = 0;
  unsigned int nEvtsNoHigherPt = 0;
  unsigned int nEvtsMcMatch    = 0;
  unsigned int nEvtsVtxIso     = 0;
  unsigned int nEvtsRelIso     = 0;
  unsigned int nAllEvts = fChain->GetEntries();
  bool isMinBias        = false;  
  // L1PixelTrackFit f(3.8112); // Bz in Tesla (for pixel re-fitting)
  
  if (cfg_DEBUG)  cout << "=== TkEG:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;
  
  // Determine what sample this is
  std::size_t found = mcSample.find("SingleNeutrino");
  std::size_t found_electron = mcSample.find("SingleE");
  if ((found!=std::string::npos) || (found_electron!=std::string::npos))
    {
      isMinBias = true;
      if (cfg_DEBUG) std::cout << "Minimum Bias sample" << std::endl;
    }
  else
    {
      if (cfg_DEBUG) std::cout << "Not a Minimum Bias sample." << std::endl;
    }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize Counters
  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Event Counters 
  unsigned int counters_events_EGs = 0;
  unsigned int counters_events_1GenTauHad_EGs = 0;
  unsigned int counters_events_2GenTauHad_EGs = 0;
  unsigned int counters_events_noEGs = 0;
  
  // Gen-Level Counters 
  unsigned int counter_allHadGenTaus = 0;
  unsigned int counter_ldgPtCut = 0; 
  unsigned int counter_sumETcut =0;

  // Leading Tracks Counters
  unsigned int counter_allTracks           = 0;
  unsigned int counter_passChi2100         = 0;
  unsigned int counter_passPtCut           = 0;
  unsigned int counter_passEtaCut          = 0;
  unsigned int counter_passChi2Nstubs      = 0;
  unsigned int counter_passNoHigherPtNeigh = 0;

  // Track clustering counters
  unsigned int trkcounter_allTracks = 0;
  unsigned int trkcounter_allNonLeading = 0;
  unsigned int trkcounter_passZ = 0;
  unsigned int trkcounter_passDRmax = 0;
  unsigned int trkcounter_passDRmin = 0;
  unsigned int trkcounter_passInvMass = 0;    

  // Clustering Counters 
  unsigned int counter_allTrkTauCand = 0;
  unsigned int counter_1trk =0;
  unsigned int counter_2trk =0;
  unsigned int counter_3trk =0;
  unsigned int counter_4trk =0;
  unsigned int counter_5trk =0;
  unsigned int counter_1eg =0;
  unsigned int counter_2eg =0;
  unsigned int counter_3eg =0;
  unsigned int counter_4eg =0;
  unsigned int counter_5eg =0;

  // MC Match Counters 
  unsigned int counter_MCmatched_Tk   = 0;
  unsigned int counter_MCmatched_TkEG = 0;
  unsigned int counter_nonMCmatched_Tk = 0;
  unsigned int counter_nonMCmatched_TkEG = 0;
  unsigned int counter_allEGspassET = 0;
  unsigned int counter_MCmatched_EG = 0;

  
  // Tau Candidate Composition Counters 
  unsigned int counter_hasGenTau = 0;
  unsigned int counter_hasEG = 0;



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // For-loop: Entries
  ////////////////////////////////////////////////////////////////////////////////////////////////

  for (int jentry = 0; jentry < nEntries; jentry++, nEvts++){
    
    if (cfg_DEBUG) cout << "\t------------ Entry = " << jentry << endl;
    
    // Init variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Build collections
    ////////////////////////////////////////////////////////////////////////////////////////////////       

    // Tracks
    if (cfg_AddL1Tks) {
      TTTracks = GetTTTracks(cfg_tk_minPt, cfg_tk_minEta, cfg_tk_maxEta, cfg_tk_maxChiSq, cfg_tk_minStubs, cfg_tk_nFitParams, false);
      sort( TTTracks.begin(), TTTracks.end() ); // Sort from highest Pt to lowest (not done by default)
      if (cfg_DEBUG*0) cout << "\n=== TTracks (" << L1Tks_Pt->size() << ")" << endl;
      if (cfg_DEBUG*0) PrintTTTrackCollection(TTTracks);
    }    
    
    // EGs
    if (cfg_AddEGs) {
      L1EGs = GetEGs(false);
      sort( L1EGs.begin(), L1EGs.end() ); // Sort from highest Et to lowest Et (should be already done by default)
      if (cfg_DEBUG) cout << "\n=== L1EGs (" << L1EGs.size() << ")" << endl;
//      if (cfg_DEBUG*0) PrintL1EGCollection(L1EGs);
    }

    // GenParticles (skip for MinBias samples as no real taus exist)
    vector<GenParticle> GenTaus;
    if (!isMinBias) GenTaus = GetGenParticles(15, true);
    if (cfg_DEBUG*0) PrintGenParticleCollection(GenTaus);
    
    // Hadronic GenTaus (skip for MinBias samples)
    GenTausHadronic.clear();
    if (cfg_AddGenP) {
      if (cfg_DEBUG) cout << "\n=== GenParticles (" << GenP_Pt.size() << ")" << endl;
      if (!isMinBias) GenTausHadronic = GetHadronicGenTaus(GenTaus, 00.0, 999.9);
      if (cfg_DEBUG*0) PrintGenParticleCollection(GenTausHadronic);
    }
    
    // Triggred GenTaus (skip for MinBias samples)
    vector<GenParticle> GenTausTrigger;  
    if (!isMinBias) GenTausTrigger = GetHadronicGenTaus(GenTaus, 20.0, 2.4);  
    
    // Ensure that all taus are found, needed by the current efficiency definition 
    // E.g. for ttbar, only events with two taus within trigger acceptance are considered for efficiency calculation)
    bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;


    // Fill histos with genTau info
    if (cfg_DEBUG) cout << "\n=== All gen-taus info " << endl;

    if (!isMinBias){
      h_genTausAll_N->Fill(GenTaus.size());
      for (vector<GenParticle>::iterator genTau = GenTaus.begin(); genTau != GenTaus.end(); genTau++) {
	h_genTausAll_Pt  -> Fill(genTau->pt());
	h_genTausAll_Eta -> Fill(genTau->eta());
	h_genTausAll_Phi -> Fill(genTau->phi());
      }
      
      if (GenTaus.size() == 2) {
	h_genTausAll_Eta1VsEta2 -> Fill (GenTaus.at(0).eta(), GenTaus.at(1).eta());
	h_genTausAll_Phi1VsPhi2 -> Fill (GenTaus.at(0).phi(), GenTaus.at(1).phi());
      }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // INITIAL STUDIES 
    ////////////////////////////////////////////////////////////////////////////////////////////////    
    
    // Check whether EGs exist or not when asking for a specific number of gen hadronic taus
    if (cfg_DEBUG) cout << "\n=== Count EGs multiplicity" << endl;
    // Without asking for genuine hadronic taus 
    if (L1EGs.size() >=1 ) counters_events_EGs++;

    // Consider only events with at least one genuine hadronic tau (except for MinBias sample)    
    h_genTausHad_N->Fill( GenTausHadronic.size() );
    //if (GenTausHadronic.size() < 1 && !isMinBias) continue;

    // Ask for at least 1 genuine hadronic tau 
    if (L1EGs.size() >=1 ) counters_events_1GenTauHad_EGs++;

    /*
    // Consider only events with at least two genuine hadronic taus (except for MinBias sample)
    if (GenTausHadronic.size() < 2 && !isMinBias) continue;

    // Ask for at least 2 genuine hadronic tau 
    if (L1EGs.size() >=1 ) counters_events_2GenTauHad_EGs++;
    */

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Gen-Level Studies
    ////////////////////////////////////////////////////////////////////////////////////////////////    
    if (cfg_DEBUG) cout << "\n\t=== Gen-Level Studies" <<endl;


    // For-loop: All hadronic gen-taus
    for (vector<GenParticle>::iterator genTau = GenTausHadronic.begin(); genTau != GenTausHadronic.end(); genTau++) {
      counter_allHadGenTaus++;
      // Fill histo with the number of the decay products of hadronic taus
      h_genTausHad_Daughters_N-> Fill(genTau->finalDaughters().size());
      h_genTausHad_chargedDaugh_N-> Fill(genTau->finalDaughtersCharged().size());
      h_genTausHad_neutralDaugh_N-> Fill(genTau->finalDaughtersNeutral().size());

      // Consider only 1-prong decays with 1 neutral product (pion0)
      if ((genTau->finalDaughtersCharged().size() == 1) && (genTau->finalDaughtersNeutral().size() == 1)) {

	if (genTau->finalDaughtersNeutral().at(0).pdgId() == 111) {

	  // Get the photons from pion0
	  GenParticle photon_1 = genTau->finalDaughtersNeutral().at(0).daughters().at(0);
	  GenParticle photon_2 = genTau->finalDaughtersNeutral().at(0).daughters().at(1);

	  double photons_dR = auxTools_.DeltaR( photon_1.eta(), photon_1.phi(), photon_2.eta(), photon_2.phi() );

	  // Fill histos (photons properties)
	  h_Pion0_Et   -> Fill (genTau->finalDaughtersNeutral().at(0).et());
	  h_Photons_Et -> Fill (photon_1.et());
	  h_Photons_Et -> Fill (photon_2.et());
	  h_Photons_dR -> Fill (photons_dR);

	  h_Pion0Et_Vs_PhotonsDR -> Fill (genTau->finalDaughtersNeutral().at(0).et(), photons_dR);
	  
	  // Check if the photons from pion0 are matched with EGs
	  EG match_EG1;
	  EG match_EG2;
	  double deltaR1;
	  double deltaR2;
	  double match_dR1 = 999;
	  double match_dR2 = 999;
	  bool photon_1_matched = false;
	  bool photon_2_matched = false;

	  // Find the minimum dR between photons and EGs
	  
	  // --- will check the matching only for pion0s with ET>10GeV (because EGs ET is >10GeV)
	  if (genTau->finalDaughtersNeutral().at(0).et() < 10) continue;

	  // For-loop: All the EGs in the event                                                                                                                         
	  for (auto eg = L1EGs.begin(); eg != L1EGs.end() ; eg++) { 
	    // 1st photon
	    deltaR1 = auxTools_.DeltaR( eg->getEta(), eg->getPhi(), photon_1.eta(), photon_1.phi() );
	    
	    if (deltaR1 < match_dR1) {
	      match_dR1 = deltaR1;
	      match_EG1 = *eg;
	    }

	    // 2nd photon
	    deltaR2 = auxTools_.DeltaR( eg->getEta(), eg->getPhi(), photon_2.eta(), photon_2.phi() );
	    
	    if (deltaR2 < match_dR2) {
	      match_dR2 = deltaR2;
	      match_EG2 = *eg;
	    }
	    
	  } // For-loop: All the EGs in the event                                                                                                                    

	  // Matching of 1st photon
	  if (match_dR1 < 0.15) {
	    photon_1_matched = true;
	  }
	  // Matching of 2nd photon
	  if (match_dR2 < 0.15) {
	    photon_2_matched = true;
	  }
	  
	  // Fill histo for the photon-EG matching (splitted into cases)
	  // None of the photons is matched
	  if ((!photon_1_matched) && (!photon_2_matched)){
	    h_Photons_EGs_Matching -> Fill(0);
	  }
	  // Only the 1 photon is matched
	  else if ( ((photon_1_matched) && (!photon_2_matched)) || ((!photon_1_matched) && (photon_2_matched)) ){
	    h_Photons_EGs_Matching  -> Fill(1);
	  }
	  // Both photons are matched
	  else if ( photon_1_matched && photon_2_matched) {
	    
	    // --- with the same EG
	    if (match_EG1.index() == match_EG2.index()) {
	      h_Photons_EGs_Matching  -> Fill(2);
	    }
	    // --- with different EGs
	    else {
	      h_Photons_EGs_Matching  -> Fill(3);
	    }
	   
	  }
	    
	}
      } 
            
      // Angular Separation of charged and neutral daughters of hadronic gen tau
      TLorentzVector p4sum;
      float deltaRcharged;
      float deltaRneutral;
      float deltaRcharged_max = -1.0;
      float deltaRneutral_max = -1.0;
      float sumET;

      // For-loop: All charged daughters
      for (unsigned int i = 0; i < genTau->finalDaughtersCharged().size(); i++)	{

	GenParticle d      = genTau->finalDaughtersCharged().at(i);
	//int genD_pdgId     = d.pdgId();
	p4sum += d.p4();
	
	// Leading charged product should have pt > 10 GeV
	if (d.pt() < 10.0) continue;

	bool notLeading = false;

	// For-loop: All charged daughters
	for (unsigned int j = 0; j < genTau->finalDaughtersCharged().size(); j++) {
	
	  GenParticle d_2      = genTau->finalDaughtersCharged().at(j);
	  
	  // Skip the leading charged product 
	  if ( d.index() == d_2.index() ) continue;

	  // Other products should have pt > 2 GeV
	  //if ( d_2.pt() < 2.0 ) continue;

	  // Skip if it is not the leading product (if it has a higher pt neighbour)
	  if ( d_2.pt() > d.pt() ) {
	    notLeading = true;
	    continue;
	  }
	  
	  // Other products should have pt > 2 GeV
	  if ( d_2.pt() >= 2.0 ) {
	    // Find the maximum deltaR of the leading charged product and the other charged products
	    deltaRcharged = auxTools_.DeltaR(d.eta(), d.phi(), d_2.eta(), d_2.phi());
	    if (deltaRcharged > deltaRcharged_max) deltaRcharged_max = deltaRcharged;
	  }

	}
	if (notLeading) continue;
	counter_ldgPtCut++;

	// Each neutral product should have ET > 10GeV (most of the times it is absorbed by one EG - EG ET threshold is 10 GeV)
	//sumET =0;
	//for (unsigned int i = 0; i < genTau->finalDaughtersNeutral().size(); i++) sumET += genTau->finalDaughtersNeutral().at(i).et();
	//if (sumET >= 2.0) {
	//counter_sumETcut++;
	// For-loop: All neutral daughters
	for (unsigned int i = 0; i < genTau->finalDaughtersNeutral().size(); i++) {
	  GenParticle d_n      = genTau->finalDaughtersNeutral().at(i);
	  
	  if (d_n.et() < 10) continue;

	  // Find the maximum deltaR of the leading charged product and the other neutral products
	  deltaRneutral = auxTools_.DeltaR(d.eta(), d.phi(), d_n.eta(), d_n.phi());
	  if (deltaRneutral > deltaRneutral_max) deltaRneutral_max = deltaRneutral;
	  //}
	}
	
	// Fill histo if it is the leading charged product
	if (!notLeading) h_genTau_chargedDaugh_visPt_dRmax -> Fill( genTau->p4vis().Pt(), deltaRcharged_max);
	if (!notLeading) h_genTau_chargedDaugh_PtLead_dRmax -> Fill( d.pt(), deltaRcharged_max);
	if (!notLeading) h_genTau_neutralDaugh_PtLead_dRmax -> Fill (d.pt(), deltaRneutral_max);
      }                
           
      if (cfg_DEBUG) genTau -> PrintFinalDaughtersCharged();

      // Fill the invariant mass and Pt of all charged products
      float chargedPt = p4sum.Pt();
      h_genTau_chargedDaugh_totalMass -> Fill(p4sum.M());
      h_genTau_chargedDaugh_Pt -> Fill (chargedPt);


      p4sum.SetPtEtaPhiM(0., 0., 0., 0.);
      float neutralET = 0.0;
      // For-loop: All neutral daughters
      for (unsigned int i = 0; i < genTau->finalDaughtersNeutral().size(); i++)	{

	GenParticle d      = genTau->finalDaughtersNeutral().at(i);
	//int genD_pdgId     = d.pdgId();
	p4sum += d.p4();
	neutralET += d.et();
      }                           

      if (cfg_DEBUG) genTau -> PrintFinalDaughtersNeutral();
      
      // Fill the invariant mass and ET of all neutral products (if exist)
      if (genTau->finalDaughtersNeutral().size()!=0){
	h_genTau_neutralDaugh_totalMass -> Fill(p4sum.M());
	h_genTau_neutralDaugh_Et -> Fill (neutralET);
	
	// Fill scatter plot of the chargedPt vs neutralET
	h_genTauHad_chargedPtVsneutralET -> Fill (chargedPt, neutralET);
	
	h_genTau_CHF -> Fill(chargedPt/ genTau->p4vis().Pt());
	h_genTau_NHF -> Fill(neutralET/ genTau->p4vis().Et());
      }

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //  Tracks + EG Algorithm
    ////////////////////////////////////////////////////////////////////////////////////////////////    
    if (cfg_DEBUG) cout << "\n\t=== Tracks + EG Algorithm" <<endl;


    //====================================================
    // Smart-counter
    //====================================================
    vector<TTTrack> tmp, tmp1;
    bool highPtNeighbourFound;
	
    // Pt Cut
    tmp = GetTTTracks(minPt_leadtrk, 0.0, 999.9, 999.9, 0, cfg_tk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedPt++;
    
    // Eta cut
    tmp = GetTTTracks(minPt_leadtrk, 0.0, maxEta_leadtrk, 999.9, 0, cfg_tk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedEta++;

    // Chi^2 cut
    tmp = GetTTTracks(minPt_leadtrk, 0.0, maxEta_leadtrk, maxChi2_trk, 0, cfg_tk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedChiSq++;

    // Stubs Multiplicity cut
    tmp = GetTTTracks(minPt_leadtrk, 0.0, maxEta_leadtrk, maxChi2_trk, minStubs_trk, cfg_tk_nFitParams, false, false);
    if (tmp.size() > 0) nEvtsSeedStubs++;

    // No higher Pt neighbour in the signal cone
    double deltaR;
    for (vector<TTTrack>::iterator seedTrkIter = tmp.begin(); seedTrkIter != tmp.end(); seedTrkIter++){
      highPtNeighbourFound = false;
      for (vector<TTTrack>::iterator trackIter = TTTracks.begin(); !highPtNeighbourFound && trackIter != TTTracks.end(); trackIter++) {
	deltaR = auxTools_.DeltaR(seedTrkIter->getEta(), seedTrkIter->getPhi(), trackIter->getEta(), trackIter->getPhi());
	if (deltaR > minDeltaR_leadtrk && deltaR < maxDeltaR_leadtrk && trackIter->getPt() > seedTrkIter->getPt())
	  highPtNeighbourFound = true;
      }
      
      if (!highPtNeighbourFound) {
	tmp1.push_back(*seedTrkIter);
      }
    }
    if (tmp1.size() > 0) nEvtsNoHigherPt++;
    
    // MC Match cut
    bool hasGenTau = false;
    GenParticle genTau;
    for (vector<TTTrack>::iterator seedTrkIter = tmp1.begin(); seedTrkIter != tmp1.end(); seedTrkIter++){
      double deltaR_match = GetMatchingGenParticle(*seedTrkIter, &genTau);
      if (deltaR_match < maxDeltaR_MCmatch) hasGenTau = true;
    }
    if (hasGenTau) nEvtsMcMatch++;

    tmp.clear();
    tmp1.clear();
    

    //====================================================
    // Leading-Tracks (tau-seeds)
    //====================================================
    if (cfg_DEBUG) std::cout << "\n=== Leading-Tracks (tau-seeds)" << endl;


    // Select lead tracks  
    vector<unsigned short> leadTrackIndices;
    trackTauCandidates.clear();
    vector< TTTrack> newCandidateTracks; // temporary container
    for (vector<TTTrack>::iterator leadTrkIter = TTTracks.begin(); leadTrkIter != TTTracks.end(); leadTrkIter++) {
      counter_allTracks++;
      if (cfg_DEBUG*0) std:: cout << "Considering leading track candidate " << leadTrkIter->index() << std::endl;
      
      // Plot quality quantities for all the TTTracks 
      h_trk_NStubs_all       -> Fill (leadTrkIter->getNumOfStubs());
      h_trk_Chi2_all         -> Fill (leadTrkIter->getChi2());
      h_trk_Chi2Red_all      -> Fill (leadTrkIter->getChi2Red());
      h_trk_NStubsVsChi2_all -> Fill (leadTrkIter->getNumOfStubs(), leadTrkIter->getChi2());

      // Pt and eta cuts
      if (leadTrkIter->getPt() < minPt_leadtrk) continue;
      counter_passPtCut++;
      if (leadTrkIter->getEta() > maxEta_leadtrk) continue;
      counter_passEtaCut++;

      // Plot chi^2 for >=5 stub leading tracks
      if (leadTrkIter->getNumOfStubs() >= minStubs_trk) h_trk_Chi2_all_5stubs -> Fill (leadTrkIter->getChi2());

      // Only use high quality tracks
      if (leadTrkIter->getNumOfStubs() < minStubs_trk) continue;
      if (leadTrkIter->getChi2() > maxChi2_trk) continue;
      //if ((leadTrkIter->getNumOfStubs() > minStubs_trk) && (leadTrkIter->getChi2() > maxChi2_trk_alt)) continue;
      counter_passChi2Nstubs++;


      if (cfg_DEBUG) std:: cout << "Leading track candidate " << leadTrkIter->index() << "passed cuts!" << std::endl;
      // Check that there are no close tracks (in terms of deltaR) with higher Pt
      highPtNeighbourFound = false;
      for (vector<TTTrack>::iterator trackIter = TTTracks.begin(); !highPtNeighbourFound && trackIter != TTTracks.end(); trackIter++) {
         if (cfg_DEBUG*0) std::cout << "Considering neighbour track " << trackIter->index() << std::endl;
         deltaR = auxTools_.DeltaR(leadTrkIter->getEta(), leadTrkIter->getPhi(), trackIter->getEta(), trackIter->getPhi());
         if (cfg_DEBUG*0) std::cout << "DeltaR = " << deltaR << std::endl;
         if (deltaR > minDeltaR_leadtrk && deltaR < maxDeltaR_leadtrk && trackIter->getPt() > leadTrkIter->getPt())
           highPtNeighbourFound = true;
           if (cfg_DEBUG*0) std::cout << "High-pT neighbour found, leading track cadidate " << leadTrkIter->index() << " discarded" << std::endl;
      }
      // If not, save the lead track to trackTauCandidates vector
      if (!highPtNeighbourFound) {
	counter_passNoHigherPtNeigh++;
        leadTrackIndices.push_back(leadTrkIter->index());
        newCandidateTracks.clear();
        newCandidateTracks.push_back(*leadTrkIter);
        trackTauCandidates.push_back(newCandidateTracks);
        // Fill lead track histograms
        h_leadTrks_Pt->Fill( leadTrkIter->getPt() );
        h_leadTrks_Eta->Fill( leadTrkIter->getEta() );
        h_leadTrks_Phi->Fill( leadTrkIter->getPhi() );
        h_leadTrks_Phi_Eta->Fill( leadTrkIter->getPhi(),leadTrkIter->getEta() );
      }
    }
    // Fill number of lead tracks (in this event) to a histogram
    h_leadTrks_Multiplicity->Fill( trackTauCandidates.size() );
    // Debug prints
    if (cfg_DEBUG) std::cout << "Lead tracks:" << std::endl;
    if (cfg_DEBUG) std::cout << trackTauCandidates.size() << " lead tracks found, here are 3 first:" << std::endl;
    for (std::size_t i=0; i<3 && i<trackTauCandidates.size(); i++) {
      if (cfg_DEBUG*0) std::cout << "leading track index = " << trackTauCandidates[i][0].index() << ", Pt = " << trackTauCandidates[i][0].getPt() << std::endl;
    }  



    //====================================================
    // Tracks Clustering (around tau-seed)
    //====================================================
    if (cfg_DEBUG) std::cout << "\n=== Track Clustering (around tau-seed)" << endl;
	
    // Cluster surrounding tracks with lead tracks
    TTTrack* leadTrackPtr = NULL;
    float invMass;
    float pT;
    int Ntrks;
    bool stopClustering;
    TLorentzVector p4sum; // initialized to (0,0,0,0)
    //PrintTTTrackCollection(TTTracks);
    for (size_t i=0; i<trackTauCandidates.size(); i++) {
      leadTrackPtr = &(trackTauCandidates[i][0]);

      if (cfg_DEBUG) cout << "Starting to cluster lead track " << leadTrackPtr->index();
      // Loop over other tracks
      stopClustering = false;

      for (vector<TTTrack>::iterator trackIter = TTTracks.begin(); ( !stopClustering && (trackIter != TTTracks.end()) ); trackIter++) {

	trkcounter_allTracks += 1;

	// Do not double-counts the lead track
	if (trackIter->index() == leadTrackPtr->index()) continue;
	trkcounter_allNonLeading += 1;

	// Skip tracks that are not close in terms of deltaR and deltaZ
	h_leadTrk_clustTrks_dZ0 -> Fill(abs (trackIter->getZ0() - leadTrackPtr->getZ0() ));
	if (abs (trackIter->getZ0() - leadTrackPtr->getZ0() ) > maxDeltaZ_trk) continue;

	trkcounter_passZ += 1;
	deltaR = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), trackIter->getEta(), trackIter->getPhi());

	if (deltaR > maxDeltaR_leadtrk) continue;
	trkcounter_passDRmax += 1;
	if (deltaR < minDeltaR_leadtrk) continue;
	trkcounter_passDRmin += 1;

	// Calculate invariant mass and pT of already clustered tracks
	p4sum.SetPtEtaPhiM(0., 0., 0., 0.);
	for(size_t j=0; j < trackTauCandidates[i].size(); j++){
	  p4sum += trackTauCandidates[i][j].p4();
	}
	invMass = p4sum.M();
	pT = p4sum.Pt();

	// Add the p4 of the new track
	p4sum += trackIter->p4();

	// Cluster the track if the total invariant mass is below the maximum  
	if (p4sum.M() < maxInvMass_trk) {
	  trackTauCandidates[i].push_back(*trackIter);
	  trkcounter_passInvMass += 1;

	  // Fill histos for clustered tracks
	  h_clustTrks_Pt->Fill( trackIter->getPt() );
	  h_clustTrks_Eta->Fill( trackIter->getEta() );
	  h_clustTrks_Phi->Fill( trackIter->getPhi() );
	  h_clustTrks_Phi_Eta->Fill( trackIter->getPhi(),trackIter->getEta() );

	  // Remove clustered track from TTTrack's collction so that it will not be used for another tau candidate
	  trackIter = TTTracks.erase(trackIter);
	  trackIter--; //WARNING!
	} 
	else {
	  // Stop clustering
	  stopClustering = true;
	}
      }

      //h_trkClusters_M_beforeCut->Fill(p4sum.M());
      // Calculate final multiplicity, pt and invariant mass of the track clusters
      Ntrks = trackTauCandidates[i].size();
      p4sum.SetPtEtaPhiM(0., 0., 0., 0.);
      for(size_t j=0; j < trackTauCandidates[i].size(); j++){
	p4sum += trackTauCandidates[i][j].p4();
      }
      invMass = p4sum.M();
      pT = p4sum.Pt();

      if (cfg_DEBUG) cout << "After clustering trackTauCandidates[" << i << "] with leadTrack " << leadTrackPtr->index() << ", it has " 
			  << trackTauCandidates[i].size() << " tracks and a mass of " << invMass << endl;

      // Fill histos with Track Clusters info
      h_trkClusters_MultiplicityPerCluster->Fill( Ntrks );
      h_trkClusters_Pt->Fill( pT );
      h_trkClusters_M->Fill( invMass );
            
      //PrintTTTrackCollection(trackTauCandidates[i]);
      counter_allTrkTauCand ++; 
      if (trackTauCandidates[i].size()==1) counter_1trk++;
      if (trackTauCandidates[i].size()==2) counter_2trk++;
      if (trackTauCandidates[i].size()==3) counter_3trk++;
      if (trackTauCandidates[i].size()==4) counter_4trk++;
      if (trackTauCandidates[i].size()==5) counter_5trk++;

    }


    //====================================================
    //  EGs Properties
    //====================================================
    if (cfg_DEBUG) std::cout << "\n=== EGs Properties" << endl;

    h_EGs_N -> Fill( L1EGs.size() );

    // For-loop: All the EGs in the event                                                                                                                               
    for (auto eg = L1EGs.begin(); eg != L1EGs.end() ; eg++) {

      h_EGs_Et  -> Fill( eg->getEt() );
      h_EGs_Eta -> Fill( eg->getEta() );
      h_EGs_Phi -> Fill( eg->getPhi() );
      
    }// For-loop: All the EGs in the event
    

    //====================================================
    //  EGs Clustering
    //====================================================
    if (cfg_DEBUG) std::cout << "\n=== EGs Clustering" << endl;

    // EGs properties (when Leading Track exist)
    if (trackTauCandidates.size() != 0) {
      
      // Fill histo with number of EGs in each event
      //h_EGs_N -> Fill( L1EGs.size() );
      //if (GenTausHadronic.size() == 1) h_EGs_N_OneHadTau -> Fill(L1EGs.size());
      
      // For-loop: All the EGs in the event
      for (auto eg = L1EGs.begin(); eg != L1EGs.end() ; eg++) {
	
	// Fill histos with EGs properties 
	//h_EGs_Et  -> Fill( eg->getEt() );
	//h_EGs_Eta -> Fill( eg->getEta() );
	//h_EGs_Phi -> Fill( eg->getPhi() );
        //h_EGs_IEta -> Fill( eg->getIEta() );
        //h_EGs_IPhi -> Fill( eg->getIPhi() );
	
      }// For-loop: All the EGs in the event
   
    }
        
    // Build EG clusters and create tau candidates
    vector<EG> EGcluster;
    TauCandidates.clear();
    double ET;
    float deltaEta, deltaPhi;
    float deltaRmin, deltaEtaMin, deltaPhiMin;
    float deltaR_beforecorrection;
    for (size_t i=0; i<trackTauCandidates.size(); i++) {
      //cout << "Clustering EGs around trackTauCandidates[" << i << "]" << endl;
      //h_EGs_N -> Fill( L1EGs.size() );
      if (GenTausHadronic.size() == 1) h_EGs_N_OneHadTau -> Fill(L1EGs.size());

      EGcluster.clear();
      leadTrackPtr = &(trackTauCandidates[i][0]);
      stopClustering = false;
      TLorentzVector p4sum; // initialized to (0,0,0,0)  //marina
      // EG clustering counters
      unsigned int counter_allEG = 0;
      unsigned int counter_passEt = 0;
      unsigned int counter_passDRmax = 0;
      unsigned int counter_passDRmin = 0;
      unsigned int counter_passInvMass = 0;    
      
      // Calculate the p4sum of the clustered tracks for each trackTau Candidate 
      for(size_t j=0; j < trackTauCandidates[i].size(); j++){
	p4sum += trackTauCandidates[i][j].p4();
      }
      
      // Initialize the ET (scalar sum of the ET of EGs) and the InvMass as InvMass of the clustered tracks
      ET      = 0.0; //= p4sum.Et();
      invMass = p4sum.M();
      deltaRmin   = 1000000.0;
      deltaEtaMin = 1000000.0;
      deltaPhiMin = 1000000.0;
      // For-loop: All the EGs in the event
      for (auto eg = L1EGs.begin(); ( !stopClustering && (eg != L1EGs.end()) ); eg++) {
	counter_allEG++;

	// Skip small-Et EGs 
	if (eg->getEt() < minEt_EG) continue;
	counter_passEt++;

	// Skip EGs which do not match to lead track (in terms of DeltaR)
	// Usage of CorrectedEta() of EG based on vertex position of the leading track
	deltaR = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), CorrectedEta(eg->getEta(), leadTrackPtr->getZ0()), eg->getPhi());
	deltaR_beforecorrection = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), eg->getEta(), eg->getPhi());

	deltaPhi = auxTools_.DeltaPhi( leadTrackPtr->getPhi(), eg->getPhi());
	deltaEta = auxTools_.DeltaEta( leadTrackPtr->getEta(), CorrectedEta(eg->getEta(), leadTrackPtr->getZ0()));

	// Fill  dR, dEta, dPhi histograms
	h_leadTrk_EG_dR -> Fill (deltaR); 
	h_leadTrk_EG_dR_beforecorrection -> Fill(deltaR_beforecorrection);
	
	h_leadTrk_EG_dPhi -> Fill (deltaPhi);
	h_leadTrk_EG_dEta -> Fill (deltaEta);

	// Find minimum dR, dEta, dPhi
	if (deltaR < deltaRmin) deltaRmin=deltaR;
	if (abs(deltaPhi) < abs(deltaPhiMin)) deltaPhiMin = deltaPhi;
	if (abs(deltaEta) < abs(deltaEtaMin)) deltaEtaMin = deltaEta;

	// Aply dR cuts
	// Skip if EG is not within the signal cone 
	if (deltaR > maxDeltaR_EG) continue;
	counter_passDRmax++;
	if (deltaR < minDeltaR_EG) continue;
	counter_passDRmin++;

	// Add the p4 of the EG to the total p4 of the clustered tracks
	p4sum += eg->p4();
	
	// Cluster the EG if the total invariant mass is below the maximum (the mass of the tau)
	if (p4sum.M() < maxInvMass_EG){
	  EGcluster.push_back(*eg);
	  counter_passInvMass++;

	  // Fill EG histograms
	  h_clustEGs_Et->Fill( eg->getEt() );
	  h_clustEGs_Eta->Fill( eg->getEta() );
	  h_clustEGs_Phi->Fill( eg->getPhi() );
	  h_clustEGs_Et_Eta->Fill( eg->getEt(),eg->getEta() );
	  h_clustEGs_Phi_Eta->Fill( eg->getPhi(),eg->getEta() );
	  h_clustEGs_M->Fill( eg->p4().M() );

	  // Update EG cluster variables
	  ET     += eg->getEt();//= p4sum.Et();
	  invMass = p4sum.M();
	  //cout<< "invMass =  "<<invMass<< "     ET = "<<endl;
	  // Remove clustered EGs from L1EG's collction so that it will not be used for another tau candidate
	  L1EGs.erase(eg);
	  eg--;
	}
	else {
	  stopClustering = true;
	}
	
	if (cfg_DEBUG) cout << "Invariant mass of cluster is " << invMass << " and ET is " << ET << endl; 
	
      }// For-loop: All the EGs in the event
      
      // Fill EG cluster variables to histograms (if >=1 EG is clustered)
      if (EGcluster.size() >=1) {
	h_EGClusters_Et->Fill( ET );            
	h_EGClusters_M->Fill( invMass );
      }
      
      h_leadTrk_EG_dRmin -> Fill (deltaRmin);
      h_leadTrk_EG_dPhiMin -> Fill (abs(deltaPhiMin));      
      h_leadTrk_EG_dEtaMin -> Fill (abs(deltaEtaMin));

      if (EGcluster.size()==1) counter_1eg++;
      if (EGcluster.size()==2) counter_2eg++;
      if (EGcluster.size()==3) counter_3eg++;
      if (EGcluster.size()==4) counter_4eg++;
      if (EGcluster.size()==5) counter_5eg++;

      // Fill EG clustering counter histogram
      h_clustEGs_allEGs      -> Fill(counter_allEG);
      h_clustEGs_passEt      -> Fill(counter_passEt);
      h_clustEGs_passDRmax   -> Fill(counter_passDRmax);
      h_clustEGs_passDRmin   -> Fill(counter_passDRmin);
      h_clustEGs_passInvMass -> Fill(counter_passInvMass);
      
      // Fill the number of EGs in cluster
      h_EGClusters_MultiplicityPerCluster->Fill( EGcluster.size() );
      
      //====================================================
      // MC Matching of Tau Candidates
      //====================================================
      if (cfg_DEBUG) std::cout << "\n=== MC Matching of Tau Candidates" << endl;


      // Find the genTau matching to the lead track
      double deltaR_match = GetMatchingGenParticle(trackTauCandidates[i][0], &genTau);
      h_MCmatch_dR->Fill(deltaR_match);

      bool hasGenTau = false;
      if (deltaR_match < maxDeltaR_MCmatch) hasGenTau = true;
      h_leadTrk_MCmatch -> Fill (hasGenTau);

      if (trackTauCandidates[i][0].getNumOfStubs() == 4) h_leadTrk4stubs_MCmatch -> Fill(hasGenTau);
      if ((trackTauCandidates[i][0].getNumOfStubs() == 4) && (hasGenTau)) h_leadTrk4stubs_MCmatched_Chi2 -> Fill (trackTauCandidates[i][0].getChi2());
      
      if (hasGenTau) {
	//Increase counter for trackTauCandidates which have a genTau
	counter_hasGenTau++;

	if (EGcluster.size() >=1) counter_hasEG++;
	
	// ET of TkEG when it has a matched genTau
	hTkEG_matched_Et  -> Fill( ET );

	// Properties of the genTau
	hTkEG_genVisEt    -> Fill( genTau.p4vis().Et() );
	hTkEG_DeltaRmatch -> Fill( deltaR_match );


	//////////////////////////////////////////////////////////////////////////////////////////
	// Tk+EG algorithm representation in Gen-Level (Charged and Neutral daughters multiplicity) 
	//////////////////////////////////////////////////////////////////////////////////////////
	if (cfg_DEBUG) std::cout << "\n=== Tk+EG algorithm representation in Gen-Level" << endl;

	// Find the leading charged daughter
	GenParticle d_ldg = genTau.finalDaughtersCharged().at(0);
	// For-loop: All charged daughters
	for (unsigned int i = 0; i < genTau.finalDaughtersCharged().size(); i++)	{
	  GenParticle d      = genTau.finalDaughtersCharged().at(i);
	  if (d.pt() > d_ldg.pt() ) d_ldg = d;
	  //p4sum += d.p4();
	}
	float sumET;
	// Leading charged product should have pt > 10 GeV
	if (d_ldg.pt() > 10.0) {
	  // If leading daughter has pt > 10GeV find the number of charged daughters with pt > 2 GeV
	  unsigned int daugh_passPtCut = 0;
	  for (unsigned int i = 0; i < genTau.finalDaughtersCharged().size(); i++) {
	    GenParticle d      = genTau.finalDaughtersCharged().at(i);
	    if (d.pt() > 2.0) daugh_passPtCut++;
	  }
	  h_MCmatch_chargedDaugh_N -> Fill(daugh_passPtCut);
	  
	  // Sum of ET of neutral products should be greater than 2 GeV                                                                                                  
	  sumET =0;
	  for (unsigned int i = 0; i < genTau.finalDaughtersNeutral().size(); i++) sumET += genTau.finalDaughtersNeutral().at(i).et();
	  if (sumET > 2.0) h_MCmatch_neutralDaugh_N -> Fill(genTau.finalDaughtersNeutral().size());
	}

      }

      // Fill counters for number of mcmatched taus (only Tks, or TkEG particles) and non mcmatched Tk or TkEGs
      if (hasGenTau) {
       
	counter_MCmatched_Tk++;

	if (EGcluster.size()>=1) {
	  counter_MCmatched_TkEG++;
	}
      }
      else {	

	counter_nonMCmatched_Tk++;

	if (EGcluster.size()>=1) {
	    counter_nonMCmatched_TkEG++;
	    
	    // If the trkTau Candidate is not MCmatched but it has EG check the fraction of the tracks energy which is deposited on the EG (for electrons ~1)
	    p4sum.SetPtEtaPhiM(0., 0., 0., 0.);
	    for(size_t j=0; j < trackTauCandidates[i].size(); j++){
	      p4sum += trackTauCandidates[i][j].p4();
	    }
	    
	    float EGenergyOverTracksPt = ET/p4sum.Pt();
	    h_nonMCmatched_EGenergyOverTracksPt -> Fill(EGenergyOverTracksPt);
	}

	// Find how the non-hadronically decaying taus (the non MCmatched) are decaying (electron or muon)  
	
	// Initialise the matching GenParticle
	GenParticle match_GenTau;
	double deltaR;
	double match_dR = 999;
	
	// For-loop: All GenTaus
	for (vector<GenParticle>::iterator tau = GenTaus.begin(); tau != GenTaus.end(); tau++) {
	  // Skip gen-taus outside of our eta acceptance 
	  if (tau->eta() > 2.4) continue;

	  // Find the closest gen-tau to the leading track
	  deltaR = auxTools_.DeltaR( trackTauCandidates[i][0].getEta(), trackTauCandidates[i][0].getPhi(), tau->eta(), tau->phi() );
	  if (deltaR < match_dR) {
	    match_dR = deltaR;
	    match_GenTau = *tau;
	  }

	}// For-loop: All GenTaus

	// If it has a matching gen tau check to what it is decaying 
	if (match_dR < maxDeltaR_MCmatch) {
	  	  
	  int decayMode = 2 ; // others (not e not muon)
	  
	  for (unsigned int i = 0; i < match_GenTau.daughters().size(); i++) {
	    
	    GenParticle d      = match_GenTau.daughters().at(i);

	    if ((abs(d.pdgId()) == 11) || (abs(d.pdgId()) == 12)) decayMode = 0; // electron
	    if ((abs(d.pdgId()) == 13) || (abs(d.pdgId()) == 14)) decayMode = 1; // muon
	  
	  }
	  h_nonMCmatchedCandidates_decayMode -> Fill (decayMode);
	  
	}
	
      }


      //====================================================
      // Build Tau Candidate
      //====================================================
      if (cfg_DEBUG) std::cout << "\n=== Build Tau Candidate" << endl;


      // Build a tau candidate from tracks and EGs
      L1TkEGParticle newTauCandidate(trackTauCandidates[i], EGcluster, genTau, hasGenTau);
      if (cfg_DEBUG) { //marina
	cout << "Constructed a new tau candidate with a track invariant mass of " << newTauCandidate.GetTrackInvMass() 
	     << ", an EG invariant mass of " << newTauCandidate.GetEGInvMass();
	if (hasGenTau) cout  << " and a generator tau with visible pT " << newTauCandidate.GetGenTauPt() << endl;
	else cout << " and does not have a matching generator tau" << endl;
      }
      
      TauCandidates.push_back(newTauCandidate);
      
    }
    
    if (TauCandidates.size() >=1){
      h_ldgTkEG_ET -> Fill( TauCandidates.at(0).GetTotalEt());
    }
    
    

    for (auto tkeg = TauCandidates.begin(); tkeg != TauCandidates.end(); tkeg++) {

      h_TkEG_Pt      -> Fill (tkeg->GetTotalPt());
      h_TkEG_ET      -> Fill (tkeg->GetTotalEt());
      h_TkEG_Eta     -> Fill (tkeg->GetTotalP4().Eta());
      h_TkEG_Phi     -> Fill (tkeg->GetTotalP4().Phi());
      h_TkEG_InvMass -> Fill (tkeg->GetTotalP4().M());
      h_TkEG_CHF     -> Fill (tkeg->GetTrackBasedPt() / tkeg->GetTotalP4().Pt());
      h_TkEG_NHF     -> Fill ( (tkeg->GetTotalP4().Pt() - tkeg->GetTrackBasedPt() ) / tkeg->GetTotalP4().Pt());
      
      if (tkeg->GetEGs().size() > 0) {
	h_TkEG_CHF_withNeutrals     -> Fill (tkeg->GetTrackBasedPt() / tkeg->GetTotalP4().Pt());
	h_TkEG_NHF_withNeutrals     -> Fill ( (tkeg->GetTotalP4().Pt() - tkeg->GetTrackBasedPt() ) / tkeg->GetTotalP4().Pt());
      }

      if (tkeg->HasMatchingGenParticle()){
	
	// Pt resolution
	double pTresolution = ( tkeg->GetTotalPt()- tkeg->GetGenTauPt() ) / tkeg->GetGenTauPt();
	h_TkEG_PtResolution->Fill(pTresolution);

	// ET resolution
	double ETresolution = ( tkeg->GetTotalEt()- tkeg->GetGenTauEt() ) / tkeg->GetGenTauEt();
	h_TkEG_EtResolution->Fill(ETresolution);

	// Eta Resolution
	double EtaResolution = ( tkeg->GetTotalP4().Eta() - tkeg->GetMatchingGenParticle().eta() ) / tkeg->GetMatchingGenParticle().eta();
        h_TkEG_EtaResolution->Fill(EtaResolution);

	// Phi Resolution
	double PhiResolution = ( tkeg->GetTotalP4().Phi() - tkeg->GetMatchingGenParticle().phi() ) / tkeg->GetMatchingGenParticle().phi();
        h_TkEG_PhiResolution->Fill(PhiResolution);


	// Resolution in different eta regions
	float tkeg_eta = tkeg->GetTotalP4().Eta();

	if ( IsWithinEtaRegion("Central", tkeg_eta) ) {
	  h_TkEG_PtResolution_C  -> Fill(pTresolution);
	  h_TkEG_EtResolution_C  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_C -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_C -> Fill(PhiResolution);
	}
	else if ( IsWithinEtaRegion("Intermediate", tkeg_eta) ) {
	  h_TkEG_PtResolution_I  -> Fill(pTresolution);
	  h_TkEG_EtResolution_I  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_I -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_I -> Fill(PhiResolution);
	}
	else if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	  h_TkEG_PtResolution_F  -> Fill(pTresolution);
	  h_TkEG_EtResolution_F  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_F -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_F -> Fill(PhiResolution);
	}
	
	// Resolution in case of zero or more than one neutral products 
	int neuDaugh_N = tkeg->GetMatchingGenParticle().finalDaughtersNeutral().size();
	
	if (neuDaugh_N == 0) {
	  h_TkEG_PtResolution_noNeutrals  -> Fill(pTresolution);
	  h_TkEG_EtResolution_noNeutrals  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_noNeutrals -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_noNeutrals -> Fill(PhiResolution);
	  
	  if ( IsWithinEtaRegion("Forward", tkeg_eta) ){
	    h_TkEG_PtResolution_noNeutrals_F  -> Fill(pTresolution);
	    h_TkEG_EtResolution_noNeutrals_F  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_noNeutrals_F -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_noNeutrals_F -> Fill(PhiResolution);
	  }
	}
	else {
	  h_TkEG_PtResolution_withNeutrals  -> Fill(pTresolution);
	  h_TkEG_EtResolution_withNeutrals  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_withNeutrals -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_withNeutrals -> Fill(PhiResolution);
	  
	  if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	    h_TkEG_PtResolution_withNeutrals_F  -> Fill(pTresolution);
	    h_TkEG_EtResolution_withNeutrals_F  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_withNeutrals_F -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_withNeutrals_F -> Fill(PhiResolution);
	  }
	}
	
	// Resolution for 1- and 3-prong decays
	int chargedDaugh_N = tkeg->GetMatchingGenParticle().finalDaughtersCharged().size();
	
	if (chargedDaugh_N == 1) {
	  h_TkEG_PtResolution_1pr  -> Fill(pTresolution);
	  h_TkEG_EtResolution_1pr  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_1pr -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_1pr -> Fill(PhiResolution);
	  
	  if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	    h_TkEG_PtResolution_1pr_F  -> Fill(pTresolution);
	    h_TkEG_EtResolution_1pr_F  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_1pr_F -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_1pr_F -> Fill(PhiResolution);
	  }
	}
	else if (chargedDaugh_N == 3){
	  h_TkEG_PtResolution_3pr  -> Fill(pTresolution);
	  h_TkEG_EtResolution_3pr  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_3pr -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_3pr -> Fill(PhiResolution);
	  
	  if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	    h_TkEG_PtResolution_3pr_F  -> Fill(pTresolution);
	    h_TkEG_EtResolution_3pr_F  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_3pr_F -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_3pr_F -> Fill(PhiResolution);
	  }
	}	
      }
      
      /*
      // Pt resolution
      TLorentzVector p4tracks; // initialized to (0,0,0,0) 
      for (auto tk = trackTauCandidates[i].begin(); tk != trackTauCandidates[i].end(); tk++) {
      p4tracks += tk->p4();
      }
      double pTresolution = ( p4tracks.Pt()-genTau.p4vis().Pt() ) / genTau.p4vis().Pt();
      h_trkClusters_PtResolution->Fill(pTresolution);
      // FIX ME: ET resolution
      if (EGcluster.size() >=1) {
      hTkEG_genVisEt_clustEG -> Fill( genTau.p4vis().Et() );
      hTkEG_genVisPt_clustEG -> Fill( genTau.p4vis().Pt() );
      
      double ETresolution = ( ET -genTau.p4vis().Et() ) / genTau.p4vis().Et(); //marina
      h_TkEG_EtResolution->Fill(ETresolution);
      }
      */
      
      
      //cout<< tkeg->GetTotalEt()<<endl;
      
      //float sumET=0;
      //for (auto clusteg = tkeg->GetEGs().begin(); clusteg != tkeg->GetEGs().end(); clusteg++){ 
      //sumET += clusteg->getEt();
      //}
      //cout <<"--- "<< tkeg->GetTrackBasedPt()+sumET<<endl;
    }
    

    //====================================================
    // Apply Isolation (Vertex, Relative)
    //====================================================
    if (cfg_DEBUG) std::cout << "\n=== Apply Isolation (Vertex, Relative)" << endl;

    // Check relative isolation of tau candidates
    float relIso = -1.0;
    float ptSum;
    TTTrack leadingTrack;
    vector<TTTrack> clustTracks;
    TauCandidatesRelIsolated.clear();
    
    for (auto tkeg = TauCandidates.begin(); tkeg != TauCandidates.end(); tkeg++) {
    
      leadingTrack = tkeg->GetLeadingTrack();
      clustTracks  = tkeg->GetTracks();

      // Sum Pt's inside the isolation cone
      ptSum = 0.0;
      // For-loop: All TTTracks
      for (auto tk = TTTracks.begin(); tk != TTTracks.end(); tk++) {

	// Check if the track is clustered in the tau candidate
	bool clustered = false;
	
	for (auto clusttk = clustTracks.begin(); clusttk != clustTracks.end(); clusttk++){
	  if (clusttk->index() == tk->index()) clustered = true;
	}
	if (clustered) continue;	
	
	// Check if the track is from the same vertex and within an isolation cone of 0.5
	if (abs (tk->getZ0() - leadingTrack.getZ0() ) > maxDeltaZ_iso) continue;
	deltaR = auxTools_.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), tk->getEta(), tk->getPhi());

        if (deltaR > minDeltaR_iso && deltaR < maxDeltaR_iso)
	  
          ptSum += tk -> getPt();
	
      } // For-loop: All TTTracks

      // Calculate relative isolation
      relIso = ptSum / tkeg->GetTotalEt();
      // Fill relative isolation histogram
      h_TkEG_relIso->Fill( relIso );

      if(relIso < maxRelIso) {
        TauCandidatesRelIsolated.push_back(*tkeg);
      }
    }

    // Check vertex isolation of tau candidates
    float mindZ0 = 1000000.0;
    float dz0;
    TauCandidatesVtxIsolated.clear();

    for (auto tkeg = TauCandidates.begin(); tkeg != TauCandidates.end(); tkeg++) {

      leadingTrack = tkeg->GetLeadingTrack();
      clustTracks  = tkeg->GetTracks();

      // For-loop: All TTTracks
      for (auto tk = TTTracks.begin(); tk != TTTracks.end(); tk++) {

	// Check if the track is clustered in the tau candidate
	bool clustered = false;
	
	for (auto clusttk = clustTracks.begin(); clusttk != clustTracks.end(); clusttk++){
	  if (clusttk->index() == tk->index()) clustered = true;
	}
	if (clustered) continue;	
	
	// Check if the track is in the iso-cone
	deltaR = auxTools_.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), tk->getEta(), tk->getPhi());
	if (deltaR > minDeltaR_iso && deltaR < maxDeltaR_iso) {
	  dz0 = abs (tk->getZ0() - leadingTrack.getZ0() );

	  if (dz0 < mindZ0) mindZ0 = dz0;
	  
	}

      } // For-loop: All TTTracks

      // Fill relative isolation histogram
      h_TkEG_vtxIso->Fill( mindZ0 );

      if(mindZ0 > maxVtxIso) {
        TauCandidatesVtxIsolated.push_back(*tkeg);
      }
    }
    
    if (TauCandidatesRelIsolated.size() > 0) nEvtsRelIso++;
    if (TauCandidatesVtxIsolated.size() > 0) nEvtsVtxIso++;



    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (cfg_DEBUG) std::cout << "\n=== Fill Turn-On histograms" << endl;


    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++)
      {
	hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() );
	
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

    FillTurnOn_Numerator_(TauCandidates     , 25.0, hTkEG_TurnOn25, hTkEG_TurnOn25_1pr, hTkEG_TurnOn25_3pr, hTkEG_TurnOn25_withNeutrals, hTkEG_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(TauCandidatesRelIsolated , 25.0, hVtxIso_TurnOn25, hVtxIso_TurnOn25_1pr, hVtxIso_TurnOn25_3pr, hVtxIso_TurnOn25_withNeutrals, hVtxIso_TurnOn25_noNeutrals); 
    FillTurnOn_Numerator_(TauCandidatesVtxIsolated , 25.0, hRelIso_TurnOn25, hRelIso_TurnOn25_1pr, hRelIso_TurnOn25_3pr, hRelIso_TurnOn25_withNeutrals, hRelIso_TurnOn25_noNeutrals);

    FillTurnOn_Numerator_(TauCandidates     , 50.0, hTkEG_TurnOn50, hTkEG_TurnOn50_1pr, hTkEG_TurnOn50_3pr, hTkEG_TurnOn50_withNeutrals, hTkEG_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(TauCandidatesRelIsolated , 50.0, hVtxIso_TurnOn50, hVtxIso_TurnOn50_1pr, hVtxIso_TurnOn50_3pr, hVtxIso_TurnOn50_withNeutrals, hVtxIso_TurnOn50_noNeutrals); 
    FillTurnOn_Numerator_(TauCandidatesVtxIsolated , 50.0, hRelIso_TurnOn50, hRelIso_TurnOn50_1pr, hRelIso_TurnOn50_3pr, hRelIso_TurnOn50_withNeutrals, hRelIso_TurnOn50_noNeutrals);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Rates and efficiencies for single tau
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (cfg_DEBUG) std::cout << "\n=== Rates and efficiencies for single tau" << endl;


    FillSingleTau_(TauCandidates, hRateSingleTau_TkEG, hEffSingleTau_TkEG);
    FillSingleTau_(TauCandidates, hRateSingleTau_C, hEffSingleTau_C, 0.0, 1.0);
    FillSingleTau_(TauCandidates, hRateSingleTau_I, hEffSingleTau_I, 1.0, 1.6);
    FillSingleTau_(TauCandidates, hRateSingleTau_F, hEffSingleTau_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(TauCandidatesRelIsolated, hRateSingleTau_relIso, hEffSingleTau_relIso);
    FillSingleTau_(TauCandidatesVtxIsolated, hRateSingleTau_vtxIso, hEffSingleTau_vtxIso);
    

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Rates and efficiencies for ditau
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (cfg_DEBUG) std::cout << "\n=== Rates and efficiencies for ditau" << endl;


    FillDiTau_(TauCandidates, hRateDiTau_TkEG, hEffDiTau_TkEG);
    FillDiTau_(TauCandidates, hRateDiTau_C, hEffDiTau_C, 0.0, 1.0);
    FillDiTau_(TauCandidates, hRateDiTau_I, hEffDiTau_I, 1.0, 1.6);
    FillDiTau_(TauCandidates, hRateDiTau_F, hEffDiTau_F, 1.6, 3.0); // 2.5 is max

    FillDiTau_(TauCandidatesRelIsolated, hRateDiTau_relIso, hEffDiTau_relIso);
    FillDiTau_(TauCandidatesVtxIsolated, hRateDiTau_vtxIso, hEffDiTau_vtxIso);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_Tk);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk  , hDiTau_Eff_Tk);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_C, hDiTau_Eff_Tk_C, 0.0, 1.0);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_I, hDiTau_Eff_Tk_I, 1.0, 1.6);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_F, hDiTau_Eff_Tk_F, 1.6, 3.0); 
    
    // Progress bar
    if (!cfg_DEBUG) auxTools_.ProgressBar(jentry, nEntries, 100, 100);
   
  } // For-loop: Entries


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Fill counters
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (cfg_DEBUG) std::cout << "\n=== Fill counters" << endl;


  // Event counters (smart counters)
  h_Counters->SetBinContent( 1, nAllEvts);
  h_Counters->SetBinContent( 2, nEvts);
  h_Counters->SetBinContent( 3, nEvtsSeedPt);
  h_Counters->SetBinContent( 4, nEvtsSeedEta);
  h_Counters->SetBinContent( 5, nEvtsSeedChiSq);
  h_Counters->SetBinContent( 6, nEvtsSeedStubs);
  h_Counters->SetBinContent( 7, nEvtsNoHigherPt);
  h_Counters->SetBinContent( 8, nEvtsRelIso);
  h_Counters->SetBinContent( 9, nEvtsVtxIso);
  //h_Counters->SetBinContent(10, nEvtsMcMatch);
  
  // events with EGs counters
  h_Counters_events_EGs->SetBinContent(1, nAllEvts);
  h_Counters_events_EGs->SetBinContent(2, counters_events_EGs);
  h_Counters_events_EGs->SetBinContent(3, counters_events_1GenTauHad_EGs);
  h_Counters_events_EGs->SetBinContent(4, counters_events_2GenTauHad_EGs);

  //cout<< "counters_events_EGs = "<< counters_events_EGs <<endl;
  //  cout<< "counters_events_noEGs = "<< counters_events_noEGs <<endl;

  // Leading Tracks Counters 
  h_Counters_leadTrks -> SetBinContent(1, counter_allTracks);
  h_Counters_leadTrks -> SetBinContent(2, counter_passChi2100);
  h_Counters_leadTrks -> SetBinContent(3, counter_passPtCut);
  h_Counters_leadTrks -> SetBinContent(4, counter_passEtaCut);
  h_Counters_leadTrks -> SetBinContent(5, counter_passChi2Nstubs);
  h_Counters_leadTrks -> SetBinContent(6, counter_passNoHigherPtNeigh);
  

  // Clustered Tracks Counters 
  h_Counters_clustTrks -> SetBinContent(1,trkcounter_allTracks);
  h_Counters_clustTrks -> SetBinContent(2,trkcounter_allNonLeading);
  h_Counters_clustTrks -> SetBinContent(3,trkcounter_passZ);
  h_Counters_clustTrks -> SetBinContent(4,trkcounter_passDRmax);
  h_Counters_clustTrks -> SetBinContent(5,trkcounter_passDRmin);
  h_Counters_clustTrks -> SetBinContent(6,trkcounter_passInvMass);
  
  // cout<<"counter_allHadGenTaus = "<< counter_allHadGenTaus<<endl;
  // cout<<"counter_ldgPtCut = " << counter_ldgPtCut<<endl;
  // cout<<"counter_sumETcut = "<< counter_sumETcut<<endl;

  // cout <<"counter_allTrkTauCand = " << counter_allTrkTauCand<<endl;
  // cout <<"counter_1trk = "<< counter_1trk<<endl;
  // cout <<"counter_2trk = "<< counter_2trk<<endl;
  // cout <<"counter_3trk = "<< counter_3trk<<endl;
  // cout <<"counter_4trk = "<< counter_4trk<<endl;
  // cout <<"counter_5trk = "<< counter_5trk<<endl;

  // cout <<"counter_1eg = "<< counter_1eg<<endl;
  // cout <<"counter_2eg = "<< counter_2eg<<endl;
  // cout <<"counter_3eg = "<< counter_3eg<<endl;
  // cout <<"counter_4eg = "<< counter_4eg<<endl;
  // cout <<"counter_5eg = "<< counter_5eg<<endl;
  
  // Fill MCmatch counters 
  h_MCmatch_counters -> SetBinContent(1, counter_allTrkTauCand);
  h_MCmatch_counters -> SetBinContent(2, counter_MCmatched_Tk );
  h_MCmatch_counters -> SetBinContent(3, counter_MCmatched_TkEG);
  h_MCmatch_counters -> SetBinContent(4, counter_nonMCmatched_Tk);
  h_MCmatch_counters -> SetBinContent(5, counter_nonMCmatched_TkEG);
  h_MCmatch_counters -> SetBinContent(6, counter_allEGspassET);
  h_MCmatch_counters -> SetBinContent(7, counter_MCmatched_EG);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Convert/Finalise histograms
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (cfg_DEBUG) std::cout << "\n=== Convert/Finalise histograms" << endl;

  // Divide turn-on numerators by the denumerator
  histoTools_.DivideHistos_1D(hTkEG_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

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

  //

  histoTools_.DivideHistos_1D(hTkEG_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hTkEG_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

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


  // Convert rate histograms
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_TkEG, nEntries);

  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_F, nEntries);

  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_relIso, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_vtxIso, nEntries);

  histoTools_.ConvertToRateHisto_1D(hRateDiTau_TkEG, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateDiTau_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateDiTau_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateDiTau_F, nEntries);

  histoTools_.ConvertToRateHisto_1D(hRateDiTau_relIso, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateDiTau_vtxIso, nEntries);

  // Finalise efficiency histograms
  FinaliseEffHisto_(hEffSingleTau_TkEG  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_(hEffSingleTau_relIso  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_vtxIso  , nEvtsWithMaxHTaus);

  FinaliseEffHisto_(hRateDiTau_TkEG  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hRateDiTau_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hRateDiTau_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hRateDiTau_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_(hRateDiTau_relIso  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hRateDiTau_vtxIso  , nEvtsWithMaxHTaus);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (cfg_DEBUG) std::cout << "\n=== Write the histograms to the file" << endl;
  
  WriteHistos_();
  auxTools_.StopwatchStop(5, "minutes", "Total Time");

}


//============================================================================
float TkEG::DeltaPhi(float phi1, float phi2) 
//============================================================================
{  float dphi = phi1 - phi2;
   if (dphi < 0) dphi = dphi + 2.*TMath::Pi();
   if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
   return dphi;
}


//============================================================================
float TkEG::deltaR(float eta1, float eta2, float phi1, float phi2) 
//============================================================================
{   float deta = eta1 - eta2;
    float dphi = DeltaPhi(phi1, phi2);
    float DR= sqrt( deta*deta + dphi*dphi );
    return DR;
}

//============================================================================ 
float TkEG::CorrectedEta(float eta, float zTrack)  {
//============================================================================ 
  // Correct the eta of the L1EG object once we know the zTrack 
  // (normaly we use the zvertex but since we care about the dR of the EG and the track it is not needed)
  
  bool IsBarrel = ( fabs(eta) < 1.479 );
  float REcal = 129. ;
  float ZEcal = 315.4 ;

  float theta = 2. * TMath::ATan( TMath::Exp( - eta ) );
  if (theta < 0) theta = theta + TMath::Pi();
  float tantheta = TMath::Tan( theta );

  float delta;
  if (IsBarrel) {
    delta = REcal / tantheta ;
  }
  else {
    if (theta > 0) delta =  ZEcal;
    if (theta < 0) delta = -ZEcal;
  }

  float tanthetaprime = delta * tantheta / (delta - zTrack );

  float thetaprime = TMath::ATan( tanthetaprime );
  if (thetaprime < 0) thetaprime = thetaprime + TMath::Pi();

  float etaprime = -TMath::Log( TMath::Tan( thetaprime / 2.) );

  return etaprime;

}

//============================================================================
bool TkEG::IsWithinEtaRegion(string etaRegion, 
				 double eta)
//============================================================================
{

  bool bWithinEtaRegion = false;
  if ( etaRegion.compare("Central") == 0 )           bWithinEtaRegion = (fabs(eta) <= _eta_C);
  else if ( etaRegion.compare("Intermediate") == 0 ) bWithinEtaRegion = (fabs(eta) <= _eta_F && fabs(eta) > _eta_C);
  else if ( etaRegion.compare("Forward") == 0 )      bWithinEtaRegion = (fabs(eta) > _eta_F);
  else{
    cout << "E R R O R ! Tracking::IsWithinEtaRegion(...) - Invalid eta region type \"" << etaRegion << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinEtaRegion;  
}



//============================================================================
void TkEG::BookHistos_(void)
//============================================================================
{

  // Binning
  const unsigned int nEt = 300;
  const float minEt = 0.0;
  const float maxEt = +300.0;


  // Event-Type Histograms
  histoTools_.BookHisto_1D(h_Counters, "Counters",  "", 10, 0.0, +9.0);
  h_Counters->GetXaxis()->SetBinLabel( 1, "All Evts");
  h_Counters->GetXaxis()->SetBinLabel( 2, "Evts");
  h_Counters->GetXaxis()->SetBinLabel( 3, "Seed Pt");
  h_Counters->GetXaxis()->SetBinLabel( 4, "Seed Eta");
  h_Counters->GetXaxis()->SetBinLabel( 5, "Seed ChiSq");
  h_Counters->GetXaxis()->SetBinLabel( 6, "Seed Stubs");
  h_Counters->GetXaxis()->SetBinLabel( 7, "Seed ldg in SigCone");
  h_Counters->GetXaxis()->SetBinLabel( 8, "RelIso");
  h_Counters->GetXaxis()->SetBinLabel( 9, "VtxIso");
  //h_Counters->GetXaxis()->SetBinLabel(10, "Matched");


  histoTools_.BookHisto_1D(h_Counters_events_EGs, "Counters_events_EGs", ";;Entries", 4, 0, 4);
  h_Counters_events_EGs->GetXaxis()->SetBinLabel(1, "All events");
  h_Counters_events_EGs->GetXaxis()->SetBinLabel(2, "Events with EG");
  h_Counters_events_EGs->GetXaxis()->SetBinLabel(3, "Events with 1 #tau_{h} & EG");
  h_Counters_events_EGs->GetXaxis()->SetBinLabel(4, "Events with 2 #tau_{h} & EG");

  // Number of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_N, "genTausAll_N", ";Number of all genuine taus in event;Events / bin", 5, -0.5, +4.5);

  // Pt of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_Pt, "genTausAll_Pt", ";p_{T} (GeV);Events / bin",20, +0.0, +100.0);

  // Eta of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_Eta, "genTausAll_Eta", ";#eta ;Events / bin", 60, -3.0, +3.0);

  // Phi of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_Phi, "genTausAll_Phi", ";#phi (rads) ;Events / bin", 21,  -3.15,  +3.15);

  // Eta of first tau Vs Eta of the second tau
  histoTools_.BookHisto_2D(h_genTausAll_Eta1VsEta2, "genTausAll_Eta1VsEta2", ";#eta_{1}; #eta_{2};Events / bin", 60, -3.0, +3.0, 60, -3.0, +3.0);

  // Phi of first tau Vs Phi of the second tau
  histoTools_.BookHisto_2D(h_genTausAll_Phi1VsPhi2, "genTausAll_Phi1VsPhi2", ";#phi_{1}; #phi_{2};Events / bin", 21,  -3.15,  +3.15, 21,  -3.15,  +3.15);

  // Number of daughters of genuine taus
  histoTools_.BookHisto_1D(h_genTausHad_Daughters_N, "genTausHad_Daughters_N", ";N_{daughters}; Events / bin", 11, -0.5, +10.5);

  // Number of charged daughters of genuine taus
  histoTools_.BookHisto_1D(h_genTausHad_chargedDaugh_N, "genTausHad_chargedDaugh_N", ";N_{charged}; Events / bin", 11, -0.5, +10.5);

  // Number of neutral daughters of genuine taus
  histoTools_.BookHisto_1D(h_genTausHad_neutralDaugh_N, "genTausHad_neutralDaugh_N", ";N_{neutral}; Events / bin", 11, -0.5, +10.5);

  // Number of genuine hadronic taus (per event)
  histoTools_.BookHisto_1D(h_genTausHad_N, "genTausHadronic_N", ";Number of genuine hadronic taus in event;Events / bin", 5, -0.5, +4.5);

  // Pt of the charged products of the tau decay
  histoTools_.BookHisto_1D(h_genTau_chargedDaugh_Pt, "genTau_chargedDaugh_Pt", ";p_{T}^{charged} (GeV);Particles / bin", 300, +0.0, +300.0);

  // TotalMass of the charged products of the tau decay 
  histoTools_.BookHisto_1D(h_genTau_chargedDaugh_totalMass, "genTau_chargedDaugh_totalMass", ";M_{charged}; GenTaus/ bin", 40, 0.0, +4.0);

  // TotalMass of the neutral products of the tau decay 
  histoTools_.BookHisto_1D(h_genTau_neutralDaugh_totalMass, "genTau_neutralDaugh_totalMass", ";M_{neutral}; GenTaus/ bin", 40, 0.0, +4.0);

  // Et of the neutral products of the tau decay
  histoTools_.BookHisto_1D(h_genTau_neutralDaugh_Et, "genTau_neutralDaugh_Et", ";E_{T}^{neutral} (GeV);Particles / bin", 300, +0.0, +300.0);

  // Pt of charged products Vs Et of neutral products of hadronic tau
  histoTools_.BookHisto_2D(h_genTauHad_chargedPtVsneutralET, "genTauHad_chargedPtVsneutralET", ";p_{T}^{charged} (GeV); E_{T}^{neutral} (GeV);Events / bin", 120, +0.0, +120.0, 20, +0.0, +120.0);
  
  // Charged Hadron Fraction
  histoTools_.BookHisto_1D(h_genTau_CHF, "genTau_CHF", ";CHF; Entries /bin", 20, 0.0, 1.0);
  
  // Neutral Hadron Fraction
  histoTools_.BookHisto_1D(h_genTau_NHF, "genTau_NHF", ";NHF; Entries /bin", 20, 0.0, 1.0);
  
  // Visible pt of charged decay products (leading pt>10GeV) Vs dRmax
  histoTools_.BookHisto_2D(h_genTau_chargedDaugh_visPt_dRmax, "genTau_chargedDaugh_visPt_dRmax", ";p_{T}^{vis} (#tau)  (GeV); #DeltaR_{max} (charged^{ldg}, charged)", 100, 0.0, 100.0, 50, 0.0, 0.5);

  // Pt of leading charged product (leading pt>10GeV) Vs dRmax with other charged products
  histoTools_.BookHisto_2D(h_genTau_chargedDaugh_PtLead_dRmax, "genTau_chargedDaugh_PtLead_dRmax", ";p_{T} (charged^{ldg}) (GeV); #DeltaR_{max} (charged^{ldg}, charged)", 100, 0.0, 100.0, 50, 0.0, 0.5);

  // Pt of leading charged decay product  (leading pt>10GeV) Vs dRmax with  other neutral products
  histoTools_.BookHisto_2D(h_genTau_neutralDaugh_PtLead_dRmax, "genTau_neutralDaugh_PtLead_dRmax", ";p_{T} (charged^{ldg}) (GeV); #DeltaR_{max} (charged^{ldg}, neutral)", 100, 0.0, 100.0, 50, 0.0, 0.5);

  // Et of pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Pion0_Et, "Pion0_Et", ";E_{T} (GeV); Particles / bin", 100, +0.0, +100.0);

  // Et of the photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Photons_Et, "Photons_Et", ";E_{T} (GeV); Particles / bin", 100, +0.0, +100.0);

  // dR of the two photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Photons_dR  , "Photons_dR"  , ";#Delta R; Entries / bin",   200,  0.0,   +2.0);

  // Et of pi0 Vs dR of the two photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_2D(h_Pion0Et_Vs_PhotonsDR, "Pion0Et_Vs_PhotonsDR", ";E_{T} (#pi^{0})  (GeV); #DeltaR (photon_{1}, photon_{2})", 100, 0.0, 100.0, 100, 0.0, 1.0);

  // Matching od EGs and Gen-Photons
  histoTools_.BookHisto_1D(h_Photons_EGs_Matching, "Photons_EGs_Matching", "Entries / bin", 4, 0, 4);
  h_Photons_EGs_Matching->GetXaxis()->SetBinLabel(1, "No matching");
  h_Photons_EGs_Matching->GetXaxis()->SetBinLabel(2, "1 Photon Matched");
  h_Photons_EGs_Matching->GetXaxis()->SetBinLabel(3, "2 Photons Matched (same EG)");
  h_Photons_EGs_Matching->GetXaxis()->SetBinLabel(4, "2 Photons Matched (diff EG)");
  h_Photons_EGs_Matching->GetXaxis()->SetLabelSize(15);
  // Number of stubs of all tracks
  histoTools_.BookHisto_1D(h_trk_NStubs_all, "trk_NStubs_all", ";N_{stubs} (all tracks); Events / bin", 16, -0.5, 15.5);

  // Chi Squared of the tracks (>=5 stubs)
  histoTools_.BookHisto_1D(h_trk_Chi2_all_5stubs, "trk_Chi2_all_5stubs", ";#chi^{2} ( #ge 5-stubs tracks); Events / bin", 100, 0.0, 100.0);

  //Chi Squared of the leading tracks when it has 4 stubs and it is MC matched
  histoTools_.BookHisto_1D(h_leadTrk4stubs_MCmatched_Chi2, "leadTrk4stubs_MCmatched_Chi2", ";#chi^{2} (4-stubs, MC matched ldg track); Events / bin", 100, 0.0, 100.0);

  // Chi Squared of all tracks
  histoTools_.BookHisto_1D(h_trk_Chi2_all, "trk_Chi2_all", ";#chi^{2} (all tracks); Events / bin", 100, 0.0, 100.0);

  // Chi Squared/ dof of all tracks
  histoTools_.BookHisto_1D(h_trk_Chi2Red_all, "trk_Chi2Red_all", ";#chi^{2}/dof (all tracks); Events / bin", 50, 0.0, 50.0);

  // Number of stubs Vs  Chi Squared for all tracks
  histoTools_.BookHisto_2D(h_trk_NStubsVsChi2_all, "trk_NStubsVsChi2_all", ";N_{stubs}; #chi^{2}; Events / bin",11, -0.5, 11.5, 100, 0.0, 100.0);

  // Leading Tracks Counters 
  histoTools_.BookHisto_1D(h_Counters_leadTrks, "Counters_leadTrks", ";;Events / bin", 6, 0, 6);
  h_Counters_leadTrks->GetXaxis()->SetBinLabel(1, "All tracks");
  h_Counters_leadTrks->GetXaxis()->SetBinLabel(2, "Pass #chi^{2}<100");
  h_Counters_leadTrks->GetXaxis()->SetBinLabel(3, "Pass p_{T} cut");
  h_Counters_leadTrks->GetXaxis()->SetBinLabel(4, "Pass #eta cut");
  h_Counters_leadTrks->GetXaxis()->SetBinLabel(5, "Pass N_{stubs} & #chi^{2} cut");
  h_Counters_leadTrks->GetXaxis()->SetBinLabel(6, "No higher p_{T} neighbour");

  // Number of lead tracks (per event)
  histoTools_.BookHisto_1D(h_leadTrks_Multiplicity, "leadTrks_Multiplicity", ";Number of lead tracks in event;Events / bin", 11, -0.5, +10.5);

  // Lead track Pt
  histoTools_.BookHisto_1D(h_leadTrks_Pt, "leadTrks_Pt", ";p_{T} (GeV);Tracks / bin", 300, +0.0, +300.0);

  // Lead track Eta
  histoTools_.BookHisto_1D(h_leadTrks_Eta, "leadTrks_Eta", ";#eta;Tracks / bin", 60, -3.0, +3.0);

  // Lead track Phi
  histoTools_.BookHisto_1D(h_leadTrks_Phi, "leadTrks_Phi", ";#phi (rads);Tracks / bin", 21,  -3.15,  +3.15);
  
  // Lead tracks in Eta-Phi plane
  // (Syntax: BookHisto_2D(histogram, hName, hTitle, binsX, xMin, xMax, binsY, yMin, yMax)
  histoTools_.BookHisto_2D(h_leadTrks_Phi_Eta, "leadTrks_Phi_Eta",  ";#phi (rads);#eta", 210,  -3.15,  +3.15, 600,  -3.0,  +3.0);  

  // DeltaZ0 of leading track with other tracks
  histoTools_.BookHisto_1D(h_leadTrk_clustTrks_dZ0, "leadTrk_clustTrks_dZ0", ";#Deltaz_{0} (cm); entries / bin", 250, 0.0, 50.0); 

  // Clustered tracks Pt
  histoTools_.BookHisto_1D(h_clustTrks_Pt, "clustTrks_Pt", ";p_{T} (GeV);Tracks / bin", 300, +0.0, +300.0);

  // Clustered tracks Eta
  histoTools_.BookHisto_1D(h_clustTrks_Eta, "clustTrks_Eta", ";#eta;Tracks / bin", 60, -3.0, +3.0);

  // Clustered tracks Phi
  histoTools_.BookHisto_1D(h_clustTrks_Phi, "clustTrks_Phi", ";#phi (rads);Tracks / bin", 21,  -3.15,  +3.15);
  
  // Clustered tracks in Eta-Phi plane
  histoTools_.BookHisto_2D(h_clustTrks_Phi_Eta, "clustTrks_Phi_Eta",  ";#phi (rads);#eta", 210,  -3.15,  +3.15, 300,  -3.0,  +3.0);  			      

  // Track clustering counter
  histoTools_.BookHisto_1D(h_Counters_clustTrks, "Counters_clustTrks", ";;Entries / bin", 6, 0., 6.);
  h_Counters_clustTrks->GetXaxis()->SetBinLabel(1,"allTracks");
  h_Counters_clustTrks->GetXaxis()->SetBinLabel(2,"allNonLeading");
  h_Counters_clustTrks->GetXaxis()->SetBinLabel(3,"passZ");
  h_Counters_clustTrks->GetXaxis()->SetBinLabel(4,"passDRmax");
  h_Counters_clustTrks->GetXaxis()->SetBinLabel(5,"passDRmin");
  h_Counters_clustTrks->GetXaxis()->SetBinLabel(6,"passInvMass");
  //h_Counters_clustTrks->SetMarkerStyle(24);

  // Number of clustered tracks (per cluster)
  histoTools_.BookHisto_1D(h_trkClusters_MultiplicityPerCluster, "trkClusters_MultiplicityPerCluster", ";Number of tracks in cluster;Clusters / bin", 16, -0.5, +15.5);

  histoTools_.BookHisto_1D(h_MCmatch_chargedDaugh_N, "MCmatch_chargedDaugh_N", ";Number of daughters^{+} of matched gen-#tau; entries / bin", 16, -0.5, +15.5);

  histoTools_.BookHisto_1D(h_MCmatch_neutralDaugh_N, "MCmatch_neutralDaugh_N", ";Number of daughters^{0} of matched gen-#tau; entries / bin", 16, -0.5, +15.5);

  // Pt of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_Pt, "trkClusters_Pt", ";p_{T} (GeV);Clusters / bin", 150, +0.0, +300.0);

  
  // Invariant mass of track clusters (before cut)
  // histoTools_.BookHisto_1D(h_trkClusters_M_beforeCut, "trkClusters_M_beforeCut", ";Invariant mass;Clusters / bin", 40, 0.0, +4.0);

  // Invariant mass of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_M, "trkClusters_M", ";Invariant mass;Clusters / bin", 40, 0.0, +4.0);

  // Number of EGs                                                                                                                                                       
  histoTools_.BookHisto_1D(h_EGs_N, "EGs_Multiplicity", ";Number of EGs in the event ; entries / bin", 21, -0.5, +20.5);

  // Number of EGs when there is only 1 hadronic tau
  histoTools_.BookHisto_1D(h_EGs_N_OneHadTau, "EGs_Multiplicity_OneHadTau", ";Number of EGs in the event ; entries / bin", 21, -0.5, +20.5);

  // MCmatched EGs Et 
  histoTools_.BookHisto_1D(h_EGs_MCmatched_Et, "EGs_MCmatched_Et", ";E_{T} (GeV) ;EGs / bin", 50, 0.0, +100.0);

  // EGs Et                                                                                                                                                       
  histoTools_.BookHisto_1D(h_EGs_Et, "EGs_Et", ";E_{T} (GeV) ;EGs / bin", 50, 0.0, +100.0);

  // EGs Eta 
  histoTools_.BookHisto_1D(h_EGs_Eta, "EGs_Eta", ";#eta ;EGs / bin", 100, -2.5, +2.5);

  // EGs Phi
  histoTools_.BookHisto_1D(h_EGs_Phi, "EGs_Phi", ";#phi (rads);EGs / bin", 36,  -3.15,  +3.15);

  // EGs IEta
//  histoTools_.BookHisto_1D(h_EGs_IEta, "EGs_IEta", ";i#eta;EGs / bin", 70, -35, +35); 

  // EGs IPhi
//  histoTools_.BookHisto_1D(h_EGs_IPhi, "EGs_IPhi", ";i#phi;EGs / bin", 36,  0,  145);

  // DR of lead trk and EGs
  histoTools_.BookHisto_1D(h_leadTrk_EG_dR, "leadTrk_EG_dR", ";#DeltaR(trk_{ldg}, EG); entries / bin",   60,  0.0,   +6.0);

  // DR of lead trk and EGs (before the correction of EG's eta
  histoTools_.BookHisto_1D(h_leadTrk_EG_dR_beforecorrection, "leadTrk_EG_dR_beforecorrection", ";#DeltaR(trk_{ldg}, EG); entries / bin",   60,  0.0,   +6.0);

  // DPhi of lead trk and EGs                                                                      
  histoTools_.BookHisto_1D(h_leadTrk_EG_dPhi, "leadTrk_EG_dPhi", ";#Delta#phi(trk_{ldg}, EG); entries / bin",  21,  -3.15,  +3.15);

  // DEta of lead trk and EGs
  histoTools_.BookHisto_1D(h_leadTrk_EG_dEta, "leadTrk_EG_dEta", ";#Delta#eta(trk_{ldg}, EG); entries / bin",  100,  -10.0,  10.0);

  // minimum DR of lead trk and EGs
  histoTools_.BookHisto_1D(h_leadTrk_EG_dRmin, "leadTrk_EG_dRmin", ";#DeltaR_{min}(trk_{ldg}, EG); entries / bin",   60,  0.0,   +6.0);

  // minimum DPhi of lead trk and EGs
  histoTools_.BookHisto_1D(h_leadTrk_EG_dPhiMin, "leadTrk_EG_dPhiMin", ";|#Delta#phi_{min}(trk_{ldg}, EG)|; entries / bin",  16,  0,  +3.2);

  // minimum DEta of lead trk and EGs
  histoTools_.BookHisto_1D(h_leadTrk_EG_dEtaMin, "leadTrk_EG_dEtaMin", ";|#Delta#eta_{min}(trk_{ldg}, EG)|; entries / bin",   50,  0.0,  10.0);

  // Clustered EGs counters
  histoTools_.BookHisto_1D(h_clustEGs_allEGs, "clustEGs_allEGs", ";N_{EGs} (all) ; entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_clustEGs_passEt, "clustEGs_passEt", ";N_{EGs} (pass E_{T}) ; entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_clustEGs_passDRmax, "clustEGs_passDRmax", ";N_{EGs} (pass #DeltaR_{max}) ; entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_clustEGs_passDRmin, "clustEGs_passDRmin", ";N_{EGs} (pass #DeltaR_{min}) ; entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_clustEGs_passInvMass, "clustEGs_passInvMass", ";N_{EGs} (pass InvMass) ; entries / bin", 16, -0.5, +15.5);

  // Clustered EGs Pt
  histoTools_.BookHisto_1D(h_clustEGs_Et, "clustEGs_Et", ";E_{T} (GeV);Events / bin", 60, +0.0, +300.0);

  // Clustered EGs Eta
  histoTools_.BookHisto_1D(h_clustEGs_Eta, "clustEGs_Eta", ";#eta;EGs / bin", 60, -3.0, +3.0);

  // Clustered EGs Phi
  histoTools_.BookHisto_1D(h_clustEGs_Phi, "clustEGs_Phi", ";#phi (rads);EGs / bin", 21,  -3.15,  +3.15);

  // Clustered EGs in Pt-Eta plane                                                           
  histoTools_.BookHisto_2D(h_clustEGs_Et_Eta, "clustEGs_Et_Eta",  ";E_{T} (GeV);#eta", 60,  0.0,  +300.0, 60,  -3.0,  +3.0);
  
  // Clustered EGs in Eta-Phi plane
  histoTools_.BookHisto_2D(h_clustEGs_Phi_Eta, "clustEGs_Phi_Eta",  ";#phi (rads);#eta", 230,  -3.15,  +3.15, 300,  -3.0,  +3.0);		      

  // Mass of clustered EGs
  histoTools_.BookHisto_1D(h_clustEGs_M, "clustEGs_M", ";Mass (GeV);clustered EGs / bin", 40, 0.0, +4.0);

  // EG clustering counter
  //histoTools_.BookHisto_2D(h_clustEGs_counter, "clustEGs_counter", ";Selection steps;clustered EGs / selection step", 5, 0.0, +5.0, 15, 0., 15.);
  histoTools_.BookHisto_1D(h_clustEGs_counter, "clustEGs_counter", ";Selection steps;clustered EGs / selection step", 5, 0.0, +5.0);
  h_clustEGs_counter->GetXaxis()->SetBinLabel(1,"allEGs");
  h_clustEGs_counter->GetXaxis()->SetBinLabel(2,"passEt");
  h_clustEGs_counter->GetXaxis()->SetBinLabel(3,"passDRmax");
  h_clustEGs_counter->GetXaxis()->SetBinLabel(4,"passDRmin");
  h_clustEGs_counter->GetXaxis()->SetBinLabel(5,"passInvMass");
  h_clustEGs_counter->SetMarkerStyle(24);

  // Number of EGs (per cluster)
  histoTools_.BookHisto_1D(h_EGClusters_MultiplicityPerCluster, "EGClusters_MultiplicityPerCluster", ";Number of EGs in cluster;Clusters / bin", 16, -0.5, +15.5);

  // Et of EG clusters
  histoTools_.BookHisto_1D(h_EGClusters_Et, "EGClusters_Et", ";Et (GeV);Clusters / bin", 150, +0.0, +150.0);

  
  // Invariant mass of EG clusters
  histoTools_.BookHisto_1D(h_EGClusters_M, "EGClusters_M", ";Invariant mass (GeV);Clusters / bin", 40, +0.0, +4.0);

  // Track-based relative isolation of tau candidates
  histoTools_.BookHisto_1D(h_TkEG_relIso, "TkEG_relIso", ";Relative isolation;Clusters / bin", 100, 0.0, +5.0);

    // Track-based relative isolation of tau candidates
  histoTools_.BookHisto_1D(h_TkEG_vtxIso, "TkEG_vtxIso", ";Vertex isolation;Clusters / bin", 100, 0.0, +5.0);
  
  // Et of the tkeg when it has a matched genTau
  histoTools_.BookHisto_1D(hTkEG_matched_Et, "TkEG_matched_Et", ";E_{T} (GeV); TkEGs / bin", 40, 0.0,  +200.0);

  
  histoTools_.BookHisto_1D(h_ldgTkEG_ET, "TkEG_ldgET", ";E_{T} (GeV); TkEGs / bin", 300, 0.0,  +300.0); 

  // DeltaR in MC matching of the lead track of tau candidates
  histoTools_.BookHisto_1D(hTkEG_DeltaRmatch  , "TkEG_DeltaRmatch"  , ";#Delta R;Particles / bin",   100,  0.0,   +1.0);
  
  // Visible Et of the matched generator taus
  histoTools_.BookHisto_1D(hTkEG_genVisEt, "TkEG_genVisEt", ";E_{T}^{vis} (gen) (GeV);Particles / bin", 40, 0.0,  +200.0);

  // Visible Et of the matched generator taus (when there is at least one EG clustered)
  histoTools_.BookHisto_1D(hTkEG_genVisEt_clustEG, "TkEG_genVisEt_clustEG", ";E_{T}^{vis} (gen) (GeV);Particles / bin", 40, 0.0,  +200.0);
  
  // Visible Pt of the matched generator taus (when there is at least one EG clustered)
  histoTools_.BookHisto_1D(hTkEG_genVisPt_clustEG, "TkEG_genVisPt_clustEG", ";p_{T}^{vis} (gen) (GeV);Particles / bin", 40, 0.0,  +200.0);

  // Properties of tau candidate 
  histoTools_.BookHisto_1D(h_TkEG_Pt, "TkEG_Pt", ";p_{T} (GeV);Entries / bin", 300, +0.0, +300.0);
  histoTools_.BookHisto_1D(h_TkEG_ET, "TkEG_ET", ";E_{T} (GeV);Entries / bin", 300, +0.0, +300.0);
  histoTools_.BookHisto_1D(h_TkEG_Eta, "TkEG_Eta", ";#eta;Entries / bin", 60, -3.0, +3.0);
  histoTools_.BookHisto_1D(h_TkEG_Phi, "TkEG_Phi", ";#phi (rads); Entries / bin", 21,  -3.15,  +3.15);
  histoTools_.BookHisto_1D(h_TkEG_InvMass, "TkEG_InvMass", ";Mass (GeV); Entries / bin", 40, 0.0, +4.0);
  histoTools_.BookHisto_1D(h_TkEG_CHF, "TkEG_CHF", ";CHF; Entries / bin", 200,  0.0,   +2.0);
  histoTools_.BookHisto_1D(h_TkEG_NHF, "TkEG_NHF", ";NHF; Entries / bin", 200,  0.0,   +2.0);
  histoTools_.BookHisto_1D(h_TkEG_CHF_withNeutrals, "TkEG_CHF_withNeutrals", ";CHF; Entries / bin", 200,  0.0,   +2.0);
  histoTools_.BookHisto_1D(h_TkEG_NHF_withNeutrals, "TkEG_NHF_withNeutrals", ";NHF; Entries / bin", 200,  0.0,   +2.0);


  // Pt resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_PtResolution, "TkEG_PtResolution", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_C, "TkEG_PtResolution_C", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_I, "TkEG_PtResolution_I", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F, "TkEG_PtResolution_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals, "TkEG_PtResolution_noNeutrals", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals, "TkEG_PtResolution_withNeutrals", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr, "TkEG_PtResolution_1pr", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr, "TkEG_PtResolution_3pr", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F, "TkEG_PtResolution_noNeutrals_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F, "TkEG_PtResolution_withNeutrals_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F, "TkEG_PtResolution_1pr_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F, "TkEG_PtResolution_3pr_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  
  // Et resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_EtResolution, "TkEG_EtResolution", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_C, "TkEG_EtResolution_C", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_I, "TkEG_EtResolution_I", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F, "TkEG_EtResolution_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals, "TkEG_EtResolution_noNeutrals", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals, "TkEG_EtResolution_withNeutrals", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr, "TkEG_EtResolution_1pr", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr, "TkEG_EtResolution_3pr", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F, "TkEG_EtResolution_noNeutrals_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F, "TkEG_EtResolution_withNeutrals_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F, "TkEG_EtResolution_1pr_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F, "TkEG_EtResolution_3pr_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);


  // Eta resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution, "TkEG_EtaResolution", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_C, "TkEG_EtaResolution_C", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_I, "TkEG_EtaResolution_I", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F, "TkEG_EtaResolution_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals, "TkEG_EtaResolution_noNeutrals", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals, "TkEG_EtaResolution_withNeutrals", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr, "TkEG_EtaResolution_1pr", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr, "TkEG_EtaResolution_3pr", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F, "TkEG_EtaResolution_noNeutrals_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F, "TkEG_EtaResolution_withNeutrals_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F, "TkEG_EtaResolution_1pr_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F, "TkEG_EtaResolution_3pr_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  
  // Phi resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution, "TkEG_PhiResolution", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_C, "TkEG_PhiResolution_C", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_I, "TkEG_PhiResolution_I", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F, "TkEG_PhiResolution_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals, "TkEG_PhiResolution_noNeutrals", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals, "TkEG_PhiResolution_withNeutrals", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr, "TkEG_PhiResolution_1pr", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr, "TkEG_PhiResolution_3pr", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F, "TkEG_PhiResolution_noNeutrals_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F, "TkEG_PhiResolution_withNeutrals_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F, "TkEG_PhiResolution_1pr_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F, "TkEG_PhiResolution_3pr_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  
  // MC match counters 
  histoTools_.BookHisto_1D(h_MCmatch_counters, "MCmatch_counters", ";;Events", 7, 0, 7);
  h_MCmatch_counters->GetXaxis()->SetBinLabel(1, "All TrkTau Candidates");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(2, "MC Matched Tks");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(3, "MC Matched TkEGs");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(4, "not MC Matched Tks");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(5, "not MC Matched TkEGs");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(6, "All EG Candidates");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(7, "MC Matched EGs");

  // Decay Mode of non MCmatched Candidates 
  histoTools_.BookHisto_1D(h_nonMCmatchedCandidates_decayMode, "nonMCmatchedCandidates_decayMode", ";;Events", 3, -0.5 , 2.5);
  h_nonMCmatchedCandidates_decayMode->GetXaxis()->SetBinLabel(1, "e");
  h_nonMCmatchedCandidates_decayMode->GetXaxis()->SetBinLabel(2, "#mu");
  h_nonMCmatchedCandidates_decayMode->GetXaxis()->SetBinLabel(3, "other");
  // TkEG particle counters 
  //histoTools_.BookHisto_1D(h_Counters_TkEGParticle, "Counters_TkEGParticle", ";;Events", 11, 0, 11);

  histoTools_.BookHisto_1D(h_nonMCmatched_EGenergyOverTracksPt, "nonMCmatched_EGenergyOverTracksPt", ";EGCluster E_{T} / Track Cluster p_{T}; Entries / bin", 100, 0.0, 5.0);

  // MC matching 
  histoTools_.BookHisto_1D(h_leadTrk_MCmatch, "leadTrk_MCmatch", ";;Events", 2, 0, 2);
  h_leadTrk_MCmatch->GetXaxis()->SetBinLabel(1, "NOT MC Matched");
  h_leadTrk_MCmatch->GetXaxis()->SetBinLabel(2, "MC Matched");

  histoTools_.BookHisto_1D(h_leadTrk4stubs_MCmatch, "leadTrk4stubs_MCmatch", ";;Events", 2, 0, 2);
  h_leadTrk4stubs_MCmatch->GetXaxis()->SetBinLabel(1, "NOT MC Matched");
  h_leadTrk4stubs_MCmatch->GetXaxis()->SetBinLabel(2, "MC Matched");

  // minimum Delta R of leading track and gen particle 
  histoTools_.BookHisto_1D(h_MCmatch_dR, "MCmatch_dR", ";#DeltaR_{min}(trk_{ldg}, gen-#tau); entries / bin",   60,  0.0,   +6.0);

  // Turn-on histograms
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt    , "McHadronicTau_VisEt"    , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_1pr, "McHadronicTau_VisEt_1pr", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_3pr, "McHadronicTau_VisEt_3pr", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_withNeutrals, "McHadronicTau_VisEt_withNeutrals", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_noNeutrals, "McHadronicTau_VisEt_noNeutrals", "", 60 , minEt , maxEt );
  

  histoTools_.BookHisto_1D(hTkEG_TurnOn25, "TkEG_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn25_1pr, "TkEG_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn25_3pr, "TkEG_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn25_withNeutrals, "TkEG_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn25_noNeutrals, "TkEG_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

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

  histoTools_.BookHisto_1D(hTkEG_TurnOn50, "TkEG_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn50_1pr, "TkEG_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn50_3pr, "TkEG_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn50_withNeutrals, "TkEG_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_TurnOn50_noNeutrals, "TkEG_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );

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
  
  // Single-tau rates
  histoTools_.BookHisto_1D(hRateSingleTau_TkEG, "Rate_SingleTau_TkEG"    , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hRateSingleTau_C  , "Rate_SingleTau_C"  , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hRateSingleTau_I  , "Rate_SingleTau_I"  , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hRateSingleTau_F  , "Rate_SingleTau_F"  , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);  
  
  histoTools_.BookHisto_1D(hRateSingleTau_relIso  , "Rate_SingleTau_RelIso"    , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hRateSingleTau_vtxIso  , "Rate_SingleTau_VtxIso"    , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);

  // Di-tau rates
  histoTools_.BookHisto_1D(hRateDiTau_TkEG, "Rate_DiTau_TkEG"    , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hRateDiTau_C  , "Rate_DiTau_C"  , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hRateDiTau_I  , "Rate_DiTau_I"  , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hRateDiTau_F  , "Rate_DiTau_F"  , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);  

  histoTools_.BookHisto_1D(hRateDiTau_relIso    , "Rate_DiTau_RelIso"    , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0); 
  histoTools_.BookHisto_1D(hRateDiTau_vtxIso    , "Rate_DiTau_VtxIso"    , ";E_{T} threshold (GeV); Rate (kHz) / bin", 400, 0.0, +400.0);
 
  // Single-tau efficiencies
  histoTools_.BookHisto_1D(hEffSingleTau_TkEG , "Eff_SingleTau_TkEG"     , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffSingleTau_C   , "Eff_SingleTau_C"   , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffSingleTau_I   , "Eff_SingleTau_I"   , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffSingleTau_F   , "Eff_SingleTau_F"   , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);

  histoTools_.BookHisto_1D(hEffSingleTau_relIso     , "Eff_SingleTau_RelIso"     , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffSingleTau_vtxIso     , "Eff_SingleTau_VtxIso"     , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);

  // Di-tau efficiencies
  histoTools_.BookHisto_1D(hEffDiTau_TkEG , "Eff_DiTau_TkEG"     , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffDiTau_C   , "Eff_DiTau_C"   , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffDiTau_I   , "Eff_DiTau_I"   , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffDiTau_F   , "Eff_DiTau_F"   , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);

  histoTools_.BookHisto_1D(hEffDiTau_relIso     , "Eff_DiTau_RelIso"     , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  histoTools_.BookHisto_1D(hEffDiTau_vtxIso     , "Eff_DiTau_VtxIso"     , ";E_{T} threshold (GeV); Efficiency / bin", 400, 0.0, +400.0);
  
  return;
}


//============================================================================
void TkEG::WriteHistos_(void)
//============================================================================
{
  // Location -> outFile
  outFile->cd();
  
  // Write Counters
  h_Counters->Write();
  h_Counters_events_EGs->Write();

  // Write 1-D histograms
  
  //--Gen-Leve
  h_genTausAll_N->Write();
  h_genTausAll_Pt->Write();
  h_genTausAll_Eta->Write();
  h_genTausAll_Phi->Write();
  h_genTausAll_Eta1VsEta2->Write();
  h_genTausAll_Phi1VsPhi2->Write();

  h_genTausHad_N->Write();
  h_genTausHad_Daughters_N->Write();
  h_genTausHad_chargedDaugh_N->Write();
  h_genTausHad_neutralDaugh_N->Write();
  h_genTau_chargedDaugh_Pt->Write();
  h_genTau_chargedDaugh_totalMass->Write();
  h_genTau_neutralDaugh_totalMass->Write();
  h_genTau_neutralDaugh_Et->Write();
  h_genTauHad_chargedPtVsneutralET->Write();
  h_genTau_CHF->Write();
  h_genTau_NHF->Write();
  h_genTau_chargedDaugh_visPt_dRmax->Write();
  h_genTau_chargedDaugh_PtLead_dRmax->Write();
  h_genTau_neutralDaugh_PtLead_dRmax->Write();

  h_Pion0_Et->Write();
  h_Photons_Et->Write();
  h_Photons_dR->Write();
  h_Pion0Et_Vs_PhotonsDR->Write();
  h_Photons_EGs_Matching->Write();

  h_trk_NStubs_all->Write();
  h_trk_Chi2_all_5stubs->Write();
  h_leadTrk4stubs_MCmatched_Chi2->Write();
  h_trk_Chi2_all->Write();
  h_trk_Chi2Red_all->Write();
  h_trk_NStubsVsChi2_all->Write();

  h_Counters_leadTrks->Write();

  h_leadTrks_Multiplicity->Write();
  h_leadTrks_Pt->Write();
  h_leadTrks_Eta->Write();
  h_leadTrks_Phi->Write();

  h_leadTrk_clustTrks_dZ0->Write();

  h_clustTrks_Pt->Write();
  h_clustTrks_Eta->Write();
  h_clustTrks_Phi->Write();
  h_Counters_clustTrks->Write();

  h_trkClusters_MultiplicityPerCluster->Write();
  h_MCmatch_chargedDaugh_N->Write();
  h_MCmatch_neutralDaugh_N->Write();
  h_trkClusters_Pt->Write();
  h_trkClusters_M->Write();
  //h_trkClusters_M_beforeCut->Write();

  h_clustEGs_Et->Write();
  h_clustEGs_Eta->Write();
  h_clustEGs_Phi->Write();
  h_clustEGs_M->Write();

  h_EGs_N -> Write();
  h_EGs_N_OneHadTau->Write();
  h_EGs_MCmatched_Et->Write();
  h_EGs_Et->Write();
  h_EGs_Eta->Write();
  h_EGs_Phi->Write();
//  h_EGs_IEta->Write();
//  h_EGs_IPhi->Write();

  h_clustEGs_allEGs      -> Write();
  h_clustEGs_passEt      -> Write();
  h_clustEGs_passDRmax   -> Write();
  h_clustEGs_passDRmin   -> Write();
  h_clustEGs_passInvMass -> Write();

  h_leadTrk_EG_dR->Write();
  h_leadTrk_EG_dR_beforecorrection->Write();
  h_leadTrk_EG_dPhi->Write();
  h_leadTrk_EG_dEta->Write();
  h_leadTrk_EG_dRmin->Write();
  h_leadTrk_EG_dPhiMin->Write();
  h_leadTrk_EG_dEtaMin->Write();

  h_EGClusters_MultiplicityPerCluster->Write();
  h_EGClusters_Et->Write();
  h_EGClusters_M->Write();

  h_TkEG_relIso->Write();
  h_TkEG_vtxIso->Write();

  h_TkEG_Pt->Write();
  h_TkEG_ET->Write();
  h_TkEG_Eta->Write();
  h_TkEG_Phi->Write();
  h_TkEG_InvMass->Write();
  h_TkEG_CHF->Write();
  h_TkEG_NHF->Write();
  h_TkEG_CHF_withNeutrals->Write();
  h_TkEG_NHF_withNeutrals->Write();
  
  h_TkEG_PtResolution->Write();
  h_TkEG_PtResolution_C->Write();
  h_TkEG_PtResolution_I->Write();
  h_TkEG_PtResolution_F->Write();
  h_TkEG_PtResolution_noNeutrals->Write();
  h_TkEG_PtResolution_withNeutrals->Write();
  h_TkEG_PtResolution_1pr->Write();
  h_TkEG_PtResolution_3pr->Write();

  h_TkEG_PtResolution_noNeutrals_F->Write();
  h_TkEG_PtResolution_withNeutrals_F->Write();
  h_TkEG_PtResolution_1pr_F->Write();
  h_TkEG_PtResolution_3pr_F->Write();
  
  h_TkEG_EtResolution->Write();
  h_TkEG_EtResolution_C->Write();
  h_TkEG_EtResolution_I->Write();
  h_TkEG_EtResolution_F->Write();
  h_TkEG_EtResolution_noNeutrals->Write();
  h_TkEG_EtResolution_withNeutrals->Write();
  h_TkEG_EtResolution_1pr->Write();
  h_TkEG_EtResolution_3pr->Write();

  h_TkEG_EtResolution_noNeutrals_F->Write();
  h_TkEG_EtResolution_withNeutrals_F->Write();
  h_TkEG_EtResolution_1pr_F->Write();
  h_TkEG_EtResolution_3pr_F->Write();
  
  h_TkEG_EtaResolution->Write();
  h_TkEG_EtaResolution_C->Write();
  h_TkEG_EtaResolution_I->Write();
  h_TkEG_EtaResolution_F->Write();
  h_TkEG_EtaResolution_noNeutrals->Write();
  h_TkEG_EtaResolution_withNeutrals->Write();
  h_TkEG_EtaResolution_1pr->Write();
  h_TkEG_EtaResolution_3pr->Write();

  h_TkEG_EtaResolution_noNeutrals_F->Write();
  h_TkEG_EtaResolution_withNeutrals_F->Write();
  h_TkEG_EtaResolution_1pr_F->Write();
  h_TkEG_EtaResolution_3pr_F->Write();

  h_TkEG_PhiResolution->Write();
  h_TkEG_PhiResolution_C->Write();
  h_TkEG_PhiResolution_I->Write();
  h_TkEG_PhiResolution_F->Write();
  h_TkEG_PhiResolution_noNeutrals->Write();
  h_TkEG_PhiResolution_withNeutrals->Write();
  h_TkEG_PhiResolution_1pr->Write();
  h_TkEG_PhiResolution_3pr->Write();

  h_TkEG_PhiResolution_noNeutrals_F->Write();
  h_TkEG_PhiResolution_withNeutrals_F->Write();
  h_TkEG_PhiResolution_1pr_F->Write();
  h_TkEG_PhiResolution_3pr_F->Write();

  hTkEG_matched_Et->Write();
  h_ldgTkEG_ET->Write();
  hTkEG_DeltaRmatch->Write();
  hTkEG_genVisEt->Write();
  hTkEG_genVisEt_clustEG->Write();
  hTkEG_genVisPt_clustEG->Write();

  
  h_MCmatch_counters->Write(); 
  h_MCmatch_dR->Write(); 
  h_leadTrk_MCmatch->Write(); 
  h_leadTrk4stubs_MCmatch->Write();
  
  h_nonMCmatched_EGenergyOverTracksPt->Write();
  h_nonMCmatchedCandidates_decayMode->Write();
  
  // Turn-on histograms
  hMcHadronicTau_VisEt->Write();
  hMcHadronicTau_VisEt_1pr->Write();
  hMcHadronicTau_VisEt_3pr->Write();
  hMcHadronicTau_VisEt_withNeutrals->Write();
  hMcHadronicTau_VisEt_noNeutrals->Write();

  hTkEG_TurnOn25->Write();
  hTkEG_TurnOn25_1pr->Write();
  hTkEG_TurnOn25_3pr->Write();
  hTkEG_TurnOn25_withNeutrals->Write();
  hTkEG_TurnOn25_noNeutrals->Write();
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

  hTkEG_TurnOn50->Write();
  hTkEG_TurnOn50_1pr->Write();
  hTkEG_TurnOn50_3pr->Write();
  hTkEG_TurnOn50_withNeutrals->Write();
  hTkEG_TurnOn50_noNeutrals->Write();
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


  hRateSingleTau_TkEG->Write(); // Inclusive = C+I+F                                                                                                                   
  hRateSingleTau_C->Write();
  hRateSingleTau_I->Write();
  hRateSingleTau_F->Write();

  hRateSingleTau_relIso->Write(); // Inclusive = C+I+F                                                                                                                   
  hRateSingleTau_vtxIso->Write(); // Inclusive = C+I+F                                                                                                                   

  hRateDiTau_TkEG->Write(); // Inclusive = C+I+F                                                                                                                       
  hRateDiTau_C->Write();
  hRateDiTau_I->Write();
  hRateDiTau_F->Write();

  hRateDiTau_relIso->Write(); // Inclusive = C+I+F                                                                                                                       
  hRateDiTau_vtxIso->Write(); // Inclusive = C+I+F                                                                                                                       

  hEffSingleTau_TkEG->Write();  // Inclusive = C+I+F                                                                                                                   
  hEffSingleTau_C->Write();
  hEffSingleTau_I->Write();
  hEffSingleTau_F->Write();

  hEffSingleTau_relIso->Write();  // Inclusive = C+I+F                                                                                                                   
  hEffSingleTau_vtxIso->Write();  // Inclusive = C+I+F                                                                                                                   

  hEffDiTau_TkEG->Write();  // Inclusive = C+I+F                                                                                                                       
  hEffDiTau_C->Write();
  hEffDiTau_I->Write();
  hEffDiTau_F->Write();

  hEffDiTau_relIso->Write();  // Inclusive = C+I+F                                                                                                                       
  hEffDiTau_vtxIso->Write();  // Inclusive = C+I+F                                                                                                                       
  
  // Write 2-D histograms
  h_leadTrks_Phi_Eta->Write();
  h_clustTrks_Phi_Eta->Write();
  h_clustEGs_Et_Eta->Write();
  h_clustEGs_Phi_Eta->Write();
  h_clustEGs_counter->Write();
  

  // Write the outfile
  outFile->Write();

}


//============================================================================
vector<L1TkEGParticle> TkEG::GetMcMatchedL1TkEGs(vector<L1TkEGParticle> L1TkEGs)
//============================================================================
{

  // Get all MC-matched trigger objects
  vector<L1TkEGParticle> matchedL1TkEGs;
  for (vector<L1TkEGParticle>::iterator tau = L1TkEGs.begin(); tau != L1TkEGs.end(); tau++)
    {
      if (!tau->HasMatchingGenParticle()) continue;
      matchedL1TkEGs.push_back(*tau);
    }
  
  return matchedL1TkEGs;
}


//============================================================================
double TkEG::GetMatchingGenParticle(TTTrack track, GenParticle *genParticlePtr)					    
//============================================================================
{

  // Match the track with a genParticle. At the moment try to match a genParticle (tau)
  // final decay products (pions, Kaons) with the track. If a match is found
  // the assign the tau (not the pion) to the L1TkEG as the matching genParticle.
  // Returns the match dR.

  // Sanity check
  if (GenTausHadronic.size() < 1 ) return 999.9; // FIXME: move outside this function

  // Initialise the GenParticle (to be returned)
  GenParticle match_GenParticle;
  double deltaR;
  double match_dR = 9999.9;
  
  // For-loop: All hadronic GenTaus
  for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++) {
  //for (vector<GenParticle>::iterator tau = GenTausTrigger.begin(); tau != GenTausTrigger.end(); tau++) {
      // If no hadronic decay products found (pi+/-, pi0, K+/-, K0, K0L), skip this tau
      if (tau->finalDaughtersCharged().size() < 1) continue;
      deltaR = auxTools_.DeltaR( track.getEta(), track.getPhi(), tau->eta(), tau->phi() );
      if (deltaR > maxDeltaR_MCmatch) continue;
      if (deltaR < match_dR) {
        match_dR = deltaR;
        *genParticlePtr = *tau;
	}
      
  } 
  
  return match_dR;
  
}      


//============================================================================
void TkEG::FillSingleTau_(vector<L1TkEGParticle> L1TkEGs, 
			    TH1D *hRate,
			    TH1D *hEfficiency,
			    double minEta,
			    double maxEta)
//============================================================================
{

  // Sanity check
  if( L1TkEGs.size() < 1 ) return;
  
  // Fill rate
  double ldgEt = L1TkEGs.at(0).GetTotalEt();

  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) < minEta) return;
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) > maxEta) return;
    
  FillRate_(hRate, ldgEt);
  
  // Get MC-matched trigger objects
  vector<L1TkEGParticle> L1TkEGs_mcMatched = GetMcMatchedL1TkEGs(TauCandidates);
  if (L1TkEGs_mcMatched.size() < 1) return;
  
  // Check that all taus were found (needed due to the efficiency defintion)
  if(!bFoundAllTaus_) return;

  // Fill efficiency
  double ldgEt_mcMatched = L1TkEGs_mcMatched.at(0).GetTotalEt();
  FillEfficiency_(hEfficiency, ldgEt_mcMatched);

  return;
}


//============================================================================
void TkEG::FillDiTau_(vector<L1TkEGParticle> L1TkEGs, 
			TH1D *hRate,
			TH1D *hEfficiency,
			double minEta,
			double maxEta)
//============================================================================
{

  // Sanity check
  if( L1TkEGs.size() < 2 ) return;  

  // Fill rate
  L1TkEGParticle L1TkEG = L1TkEGs.at(1);
  double subLdgEt = L1TkEG.GetTotalEt();
  
  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) < minEta) return;
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) > maxEta) return;
  if ( abs(L1TkEGs.at(1).GetLeadingTrack().getEta()) < minEta) return;
  if ( abs(L1TkEGs.at(1).GetLeadingTrack().getEta()) > maxEta) return;

  FillRate_(hRate, subLdgEt);

  // Get MC-Matched trigger objects
  vector<L1TkEGParticle> L1TkEGs_mcMatched = GetMcMatchedL1TkEGs(L1TkEGs);
  if (L1TkEGs_mcMatched.size() < 2) return;
    
  // Check that all taus were found (needed due to the efficiency defintion)
  if(!bFoundAllTaus_) return;

  // Fill efficiency
  double subLdgEt_mcMatched = L1TkEGs_mcMatched.at(1).GetTotalEt();
  FillEfficiency_(hEfficiency, subLdgEt_mcMatched);

  return;
}


//============================================================================
void TkEG::FillDiTau_(vector<L1TkEGParticle> L1TkEGs1,
				vector<L1TkEGParticle> L1TkEGs2, 
				TH2D *hRate,
				TH2D *hEfficiency)
//============================================================================
{

  // Sanity check
  if( L1TkEGs1.size() < 1 ) return;
  if( L1TkEGs2.size() < 1 ) return;
  
  // Get MC-Matched trigger objects
  vector<L1TkEGParticle> L1TkEGs1_mcMatched = GetMcMatchedL1TkEGs(L1TkEGs1);
  vector<L1TkEGParticle> L1TkEGs2_mcMatched = GetMcMatchedL1TkEGs(L1TkEGs2);

  // Fill rate 
  double ldgEt1 = L1TkEGs1.at(0).GetTotalEt();
  double ldgEt2 = L1TkEGs2.at(0).GetTotalEt();

  // Ensure that different calo objects are used
  unsigned int index1 = L1TkEGs1.at(0).GetLeadingTrack().index();
  unsigned int index2 = L1TkEGs2.at(0).GetLeadingTrack().index();
  if (index1==index2)
    {
      if (L1TkEGs2.size() < 2) return;
      index2 = L1TkEGs2.at(1).GetLeadingTrack().index();
    }

  // Make x-axis the ldgEt axis
  if (ldgEt1 > ldgEt2) FillRate_(hRate, ldgEt1, ldgEt2); 
  else FillRate_(hRate, ldgEt2, ldgEt1);

  
  // Get MC-matched trigger objects
  if (L1TkEGs1_mcMatched.size() < 1) return;
  if (L1TkEGs2_mcMatched.size() < 1) return;

  // Get MC-matched Et
  double ldgEt1_mcMatched = L1TkEGs1_mcMatched.at(0).GetTotalEt();
  double ldgEt2_mcMatched = L1TkEGs2_mcMatched.at(0).GetTotalEt();

  
  // Ensure that different calo objects are used
  index1 = L1TkEGs1_mcMatched.at(0).GetLeadingTrack().index();
  index2 = L1TkEGs2_mcMatched.at(0).GetLeadingTrack().index();
  if (index1==index2)
    {
      if (L1TkEGs2_mcMatched.size() < 2) return;
      index2 = L1TkEGs2_mcMatched.at(1).GetLeadingTrack().index();
    }

  // Check that all taus were found
  //if(!bFoundAllTaus_) return; //TODO

  // Make x-axis the ldgEt axis
  if (ldgEt1_mcMatched > ldgEt2_mcMatched) histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt1_mcMatched, ldgEt2_mcMatched);
  else histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt2_mcMatched, ldgEt1_mcMatched);

  return;
}



//============================================================================
void TkEG::FillRate_(TH1D *hRate,
			   const double ldgEt)
//============================================================================
{
  
  if (ldgEt < 0) return;
  hRate ->Fill( ldgEt );
  return;
}


//============================================================================
void TkEG::FillRate_(TH2D *hRate,
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
void TkEG::FillEfficiency_(TH1D *hEfficiency,
			     const double ldgEt)
//============================================================================
{
  
  histoTools_.FillAllBinsUpToValue_1D(hEfficiency, ldgEt);

  return;
}


//============================================================================
void TkEG::FillTurnOn_Numerator_(vector<L1TkEGParticle> L1TkEGs, 
				 const double minEt,
				 TH1D *hTurnOn,
				 TH1D *hTurnOn_1pr,
				 TH1D *hTurnOn_3pr,
				 TH1D *hTurnOn_withNeutrals,
				 TH1D *hTurnOn_noNeutrals)
//============================================================================
{
  
  // For-loop: L1TkEGs
  for (vector<L1TkEGParticle>::iterator L1TkEG = L1TkEGs.begin(); L1TkEG != L1TkEGs.end(); L1TkEG++)
    {
      
      // Skip if trigger object is not MC matched
      if (!L1TkEG->HasMatchingGenParticle()) continue;	 
      
      // Skip if trigger object has eT < minEt
      if (L1TkEG->GetTotalEt() < minEt) continue;
            
      // Get MC-match
      GenParticle p = L1TkEG->GetMatchingGenParticle();
      
      // Fill the turn-on
      hTurnOn->Fill( p.p4vis().Et() );

	
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
      
    } // For-loop: L1TkEGs

  return;

}

//============================================================================
void TkEG::FinaliseEffHisto_(TH1D *histo, 
				       const int nEvtsTotal)
//============================================================================
{

  const int nBins = histo->GetNbinsX()+1;
  double eff, err;

  // For-loop: Histogram bins
  for (int i = 1; i<= nBins; i++){
    
    const int nPass = histo->GetBinContent(i);
    auxTools_.Efficiency(nPass, nEvtsTotal, "binomial", eff, err ); //fixme: use TEfficiency?

    // Update current histo bin to true eff value and error
    histo->SetBinContent(i, eff);
    histo->SetBinError  (i, err);
  }

  return;
}


//============================================================================
void TkEG::FinaliseEffHisto_(TH2D *histo, 
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


#endif

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
  DEBUG = false;

  // Dataset-related
  datasets_  = datasets_.GetDataset(mcSample);
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;
  
  if (DEBUG) std::cout << "=== TkEG::InitVars_()" << std::endl;
  
  cfg_AddL1Tks   = true;
  cfg_AddEGs     = true;
  cfg_AddGenP    = true;
  
  // Track parameters
  cfg_tk_Collection  =  "TTTracks"; // Default: "TTTracks" (not "TTPixelTracks")
  cfg_tk_nFitParams  = 4;           // Default: 4
  cfg_tk_minPt       = 2.00;        // Default: 2.0
  cfg_tk_minEta      = 0.0;         // Default: 0.0
  cfg_tk_maxEta      = 1.5;//2.5;        // Default: 1e6
  cfg_tk_maxChiSq    = 50.0;         // Default: 1e6
  cfg_tk_minStubs    =   5;         // Default: 0

  // EGs parameters
  cfg_eg_minEt       = 2.00;        // Default: 2.0
  cfg_eg_minEta      = 0.0;         // Default: 0.0
  cfg_eg_maxEta      = 1.5;//2.5;        // Default: 1e6   

  // TkEG algorithm parameters
  minStubs_trk      = 5; 
  maxChi2_trk       = 50.0; // GeV
  minPt_leadtrk     = 5.0; // GeV
  maxEta_leadtrk    = 1.5;//2.5;
  minDeltaR_leadtrk = 0.0;
  maxDeltaR_leadtrk = 0.15;//0.3;
  maxDeltaR_const   = 2.5;
  maxDeltaZ_trk     = 0.8;  // cm
  maxInvMass_trk    = 1.5; // GeV 
  minEt_EG          = 1.5;  // GeV
  minDeltaR_EG      = 0.0;
  maxDeltaR_EG      = 0.15;//0.3;
  maxInvMass_EG     = 1.77; // GeV
  maxDeltaR_MCmatch = 0.1;    

  
  //minDeltaR_iso     = 0.15;
  maxDeltaR_iso     = 0.3;
  maxDeltaZ_iso     = 0.6;  // cm
  useRelIso         = true;
  relIso_WP         = 0.10; //0.2
  vtxIso_WP         = 0.50;

  // Double-tau
  diTau_deltaPOCAz = +1.00; // cm

  // Eta regions
  _eta_C = 0.8; 
  _eta_F = 1.6;


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
  unsigned int nEvtsVtxIsoLoose= 0;
  unsigned int nEvtsVtxIsoTight= 0;
  unsigned int nEvtsRelIsoLoose= 0;
  unsigned int nEvtsRelIsoTight= 0;
  unsigned int nEvtsIso        = 0;
  unsigned int nAllEvts = fChain->GetEntries();
  bool isMinBias        = false;  
  // L1PixelTrackFit f(3.8112); // Bz in Tesla (for pixel re-fitting)
  
  if (DEBUG)  cout << "=== TkEG:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;
  
  // Determine what sample this is
  std::size_t found = mcSample.find("SingleNeutrino");
  std::size_t found_electron = mcSample.find("SingleE");
  if ((found!=std::string::npos) || (found_electron!=std::string::npos))
    {
      isMinBias = true;
      if (DEBUG) std::cout << "Minimum Bias sample" << std::endl;
    }
  else
    {
      if (DEBUG) std::cout << "Not a Minimum Bias sample." << std::endl;
    }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialize Counters
  ////////////////////////////////////////////////////////////////////////////////////////////////

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

  // Tau Candidate Composition Counters 
  unsigned int counter_hasGenTau = 0;
  unsigned int counter_hasEG = 0;



  ////////////////////////////////////////////////////////////////////////////////////////////////
  // For-loop: Entries
  ////////////////////////////////////////////////////////////////////////////////////////////////

  for (int jentry = 0; jentry < nEntries; jentry++, nEvts++){
    
    if (DEBUG) cout << "\t------------ Entry = " << jentry << endl;
    
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
      if (DEBUG*0) cout << "\n=== TTracks (" << L1Tks_Pt->size() << ")" << endl;
      if (DEBUG*0) PrintTTTrackCollection(TTTracks);
    }    
    
    // EGs
    if (cfg_AddEGs) {
      L1EGs = GetEGs(cfg_eg_minEt, cfg_eg_minEta, cfg_eg_maxEta, false);
      sort( L1EGs.begin(), L1EGs.end() ); // Sort from highest Et to lowest Et (should be already done by default)
      if (DEBUG) cout << "\n=== L1EGs (" << L1EGs.size() << ")" << endl;
//      if (DEBUG*0) PrintL1EGCollection(L1EGs);
    }

    // GenParticles (skip for MinBias samples as no real taus exist)
    vector<GenParticle> GenTaus;
    if (!isMinBias) GenTaus = GetGenParticles(15, true);
    if (DEBUG*0) PrintGenParticleCollection(GenTaus);
    
    // Hadronic GenTaus (skip for MinBias samples)
    GenTausHadronic.clear();
    if (cfg_AddGenP) {
      if (DEBUG) cout << "\n=== GenParticles (" << GenP_Pt.size() << ")" << endl;
      if (!isMinBias) GenTausHadronic = GetHadronicGenTaus(GenTaus, 00.0, 1.4);//999.9
      if (DEBUG*0) PrintGenParticleCollection(GenTausHadronic);
    }
    
    // Triggred GenTaus (skip for MinBias samples)
    vector<GenParticle> GenTausTrigger;  
    if (!isMinBias) GenTausTrigger = GetHadronicGenTaus(GenTaus, 20.0, 1.4);//2.3);  
    
    // Ensure that all taus are found, needed by the current efficiency definition 
    // E.g. for ttbar, only events with two taus within trigger acceptance are considered for efficiency calculation)
    bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;


    // Fill histos with genTau info
    if (DEBUG) cout << "\n=== All gen-taus info " << endl;

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
    
    // Consider only events with at least one genuine hadronic tau (except for MinBias sample)    
    h_genTausHad_N->Fill( GenTausHadronic.size() );
    //if (GenTausHadronic.size() < 1 && !isMinBias) continue;


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Gen-Level Studies
    ////////////////////////////////////////////////////////////////////////////////////////////////    
    if (DEBUG) cout << "\n\t=== Gen-Level Studies" <<endl;


    // For-loop: All hadronic gen-taus
    for (vector<GenParticle>::iterator genTau = GenTausHadronic.begin(); genTau != GenTausHadronic.end(); genTau++) {

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

	  double photons_dR   = auxTools_.DeltaR( photon_1.eta(), photon_1.phi(), photon_2.eta(), photon_2.phi() );
	  double photons_dEta = auxTools_.DeltaEta( photon_1.eta(), photon_2.eta());
	  double photons_dPhi = auxTools_.DeltaPhi( photon_1.phi(), photon_2.phi());
						   
	  // Fill histos (photons properties)
	  h_Pion0_Et     -> Fill (genTau->finalDaughtersNeutral().at(0).et());
	  h_Photons_Et   -> Fill (photon_1.et());
	  h_Photons_Et   -> Fill (photon_2.et());
	  h_Photons_dR   -> Fill (photons_dR);
	  h_Photons_dEta -> Fill (photons_dEta);
	  h_Photons_dPhi -> Fill (photons_dPhi);
	
	  h_Photons_dEtaVsdPhi   -> Fill (photons_dEta, photons_dPhi);
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
	  
	  // --- will check the matching only for pion0s with ET>1.5GeV (because EGs ET is >10GeV)
	  if (genTau->finalDaughtersNeutral().at(0).et() < 1.5) continue;

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
	
	// Leading charged product should have pt > 5 GeV
	if (d.pt() < 5.0) continue;

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

	// Each neutral product should have ET > 10GeV (most of the times it is absorbed by one EG - EG ET threshold is 10 GeV)

	// For-loop: All neutral daughters
	for (unsigned int i = 0; i < genTau->finalDaughtersNeutral().size(); i++) {
	  GenParticle d_n      = genTau->finalDaughtersNeutral().at(i);
	  
	  if (d_n.et() < 1.5) continue;

	  // Find the maximum deltaR of the leading charged product and the other neutral products
	  deltaRneutral = auxTools_.DeltaR(d.eta(), d.phi(), d_n.eta(), d_n.phi());
	  if (deltaRneutral > deltaRneutral_max) deltaRneutral_max = deltaRneutral;
	  //}
	}
	
	// Fill histo if it is the leading charged product
	if (!notLeading) h_genTau_chargedDaugh_visPt_dRmax -> Fill( genTau->p4vis().Pt(), deltaRcharged_max);
	if (!notLeading) h_genTau_chargedDaugh_PtLead_dRmax -> Fill( deltaRcharged_max, d.pt());
	if (!notLeading) h_genTau_neutralDaugh_PtLead_dRmax -> Fill( deltaRneutral_max, d.pt());
      }                
           
      if (DEBUG) genTau -> PrintFinalDaughtersCharged();

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

      if (DEBUG) genTau -> PrintFinalDaughtersNeutral();
      
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
    if (DEBUG) cout << "\n\t=== Tracks + EG Algorithm" <<endl;


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
    double deltaR, maxDeltaR;
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
    if (DEBUG) std::cout << "\n=== Leading-Tracks (tau-seeds)" << endl;


    // Select lead tracks  
    vector<unsigned short> leadTrackIndices;
    trackTauCandidates.clear();
    vector< TTTrack> newCandidateTracks; // temporary container
    for (vector<TTTrack>::iterator leadTrkIter = TTTracks.begin(); leadTrkIter != TTTracks.end(); leadTrkIter++) {
      counter_allTracks++;
      if (DEBUG*0) std:: cout << "Considering leading track candidate " << leadTrkIter->index() << std::endl;
      
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
      counter_passChi2Nstubs++;

      if (DEBUG) std:: cout << "Leading track candidate " << leadTrkIter->index() << "passed cuts!" << std::endl;
      // Check that there are no close tracks (in terms of deltaR) with higher Pt
      highPtNeighbourFound = false;
      for (vector<TTTrack>::iterator trackIter = TTTracks.begin(); !highPtNeighbourFound && trackIter != TTTracks.end(); trackIter++) {
         if (DEBUG*0) std::cout << "Considering neighbour track " << trackIter->index() << std::endl;
         deltaR = auxTools_.DeltaR(leadTrkIter->getEta(), leadTrkIter->getPhi(), trackIter->getEta(), trackIter->getPhi());
         if (DEBUG*0) std::cout << "DeltaR = " << deltaR << std::endl;
         if (deltaR > minDeltaR_leadtrk && deltaR < maxDeltaR_leadtrk && trackIter->getPt() > leadTrkIter->getPt())
           highPtNeighbourFound = true;
           if (DEBUG*0) std::cout << "High-pT neighbour found, leading track cadidate " << leadTrkIter->index() << " discarded" << std::endl;
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
    if (DEBUG) std::cout << "Lead tracks:" << std::endl;
    if (DEBUG) std::cout << trackTauCandidates.size() << " lead tracks found, here are 3 first:" << std::endl;
    for (std::size_t i=0; i<3 && i<trackTauCandidates.size(); i++) {
      if (DEBUG*0) std::cout << "leading track index = " << trackTauCandidates[i][0].index() << ", Pt = " << trackTauCandidates[i][0].getPt() << std::endl;
    }  


    //====================================================
    // Tracks Clustering (around tau-seed)
    //====================================================
    if (DEBUG) std::cout << "\n=== Track Clustering (around tau-seed)" << endl;
    
    // Cluster surrounding tracks with lead tracks
    TTTrack* leadTrackPtr = NULL;
    float invMass, pT;
    int Ntrks;
    bool stopClustering;
        
    for (size_t i=0; i<trackTauCandidates.size(); i++) {
      
      // Initialize the p4sum for each tau candidate
      TLorentzVector p4sum; // initialized to (0,0,0,0)

      leadTrackPtr = &(trackTauCandidates[i][0]);
      p4sum += leadTrackPtr->p4();
      
      // Shrinking cone size
      maxDeltaR = maxDeltaR_const/leadTrackPtr->getPt();
      if (maxDeltaR > maxDeltaR_leadtrk) maxDeltaR = maxDeltaR_leadtrk;

      h_SigCone_DeltaR -> Fill (maxDeltaR);

      if (DEBUG) cout << "Starting to cluster lead track " << leadTrackPtr->index();
      // Loop over other tracks
      //stopClustering = false;

      for (vector<TTTrack>::iterator trackIter = TTTracks.begin(); trackIter != TTTracks.end(); trackIter++) {

	trkcounter_allTracks += 1;
	// Do not double-counts the lead track
	if (trackIter->index() == leadTrackPtr->index()) continue;
	trkcounter_allNonLeading += 1;

	// Skip tracks that are not close in terms of deltaR and deltaZ
	h_leadTrk_clustTrks_dZ0 -> Fill(abs (trackIter->getZ0() - leadTrackPtr->getZ0() ));
	if (abs (trackIter->getZ0() - leadTrackPtr->getZ0() ) > maxDeltaZ_trk) continue;

	trkcounter_passZ += 1;
	deltaR = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), trackIter->getEta(), trackIter->getPhi());

	if (deltaR > maxDeltaR) continue;
	trkcounter_passDRmax += 1;
	if (deltaR < minDeltaR_leadtrk) continue;
	trkcounter_passDRmin += 1;

	// Add the p4 of the new track
	p4sum += trackIter->p4();
	trackTauCandidates[i].push_back(*trackIter);

	// Fill histos for clustered tracks
	// h_clustTrks_Pt->Fill( trackIter->getPt() );
	// h_clustTrks_Eta->Fill( trackIter->getEta() );
	// h_clustTrks_Phi->Fill( trackIter->getPhi() );
	// h_clustTrks_Phi_Eta->Fill( trackIter->getPhi(),trackIter->getEta() );

      }
      
      h_trkClusters_M_beforeCut->Fill(p4sum.M());
      
      // Keep the trackCandidate if the total invariant mass is below the maximum  
      if (p4sum.M() > maxInvMass_trk) {

	trackTauCandidates.erase(trackTauCandidates.begin()+i);
	i--;
      }

      else {

	// Calculate final multiplicity, pt and invariant mass of the track clusters
	Ntrks = trackTauCandidates[i].size();
	p4sum.SetPtEtaPhiM(0., 0., 0., 0.);
	for(size_t j=0; j < trackTauCandidates[i].size(); j++){
	  p4sum += trackTauCandidates[i][j].p4();
	}
	invMass = p4sum.M();
	pT = p4sum.Pt();
	
	if (DEBUG) cout << "After clustering trackTauCandidates[" << i << "] with leadTrack " << leadTrackPtr->index() << ", it has " 
			    << trackTauCandidates[i].size() << " tracks and a mass of " << invMass << endl;
	
	// Fill histos with Track Clusters info
	h_trkClusters_MultiplicityPerCluster->Fill( Ntrks );
	h_trkClusters_Pt->Fill( pT );
	h_trkClusters_M->Fill( invMass );
	
      }
    }

    //====================================================
    //  EGs Properties
    //====================================================
    if (DEBUG) std::cout << "\n=== EGs Properties" << endl;

    h_EGs_N -> Fill( L1EGs.size() );

    // For-loop: All the EGs in the event
    for (auto eg = L1EGs.begin(); eg != L1EGs.end() ; eg++) {

	h_EGs_Et      -> Fill( eg->getEt() );
	h_EGs_Eta     -> Fill( eg->getEta() );
	h_EGs_Phi     -> Fill( eg->getPhi() );
	h_EGs_EtaVsEt -> Fill( eg->getEta(), eg->getEt() ); 
	
	h_EGs_HwQual  -> Fill (eg->getHwQual());

    }// For-loop: All the EGs in the event
    

    //====================================================
    //  EGs Clustering
    //====================================================
    if (DEBUG) std::cout << "\n=== EGs Clustering" << endl;

    // Build EG clusters and create tau candidates
    vector<EG> EGcluster;
    L1TkEGTauCandidates.clear();
    double ET;
    float deltaEta, deltaPhi;
    float deltaRmin, deltaEtaMin, deltaPhiMin;
    float deltaR_beforecorrection;
    for (size_t i=0; i<trackTauCandidates.size(); i++) {

      EGcluster.clear();
      leadTrackPtr = &(trackTauCandidates[i][0]);
      stopClustering = false;

      // Shrinking cone size
      maxDeltaR = maxDeltaR_const/leadTrackPtr->getPt();
      if (maxDeltaR > maxDeltaR_EG) maxDeltaR = maxDeltaR_leadtrk;
      
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
	if (deltaR > maxDeltaR) continue;
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
	  //L1EGs.erase(eg);
	  //eg--;
	}
	else {
	  stopClustering = true;
	}
	
	if (DEBUG) cout << "Invariant mass of cluster is " << invMass << " and ET is " << ET << endl; 
	
      }// For-loop: All the EGs in the event
      
      // Fill EG cluster variables to histograms (if >=1 EG is clustered)
      if (EGcluster.size() >=1) {
	h_EGClusters_Et->Fill( ET );            
	h_EGClusters_M->Fill( invMass );
      }
      
      h_leadTrk_EG_dRmin -> Fill (deltaRmin);
      h_leadTrk_EG_dPhiMin -> Fill (abs(deltaPhiMin));      
      h_leadTrk_EG_dEtaMin -> Fill (abs(deltaEtaMin));

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
      if (DEBUG) std::cout << "\n=== MC Matching of Tau Candidates" << endl;


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
	
	// Properties of the genTau
	hTkEG_genVisEt    -> Fill( genTau.p4vis().Et() );
	hTkEG_DeltaRmatch -> Fill( deltaR_match );


	//////////////////////////////////////////////////////////////////////////////////////////
	// Tk+EG algorithm representation in Gen-Level (Charged and Neutral daughters multiplicity) 
	//////////////////////////////////////////////////////////////////////////////////////////
	if (DEBUG) std::cout << "\n=== Tk+EG algorithm representation in Gen-Level" << endl;

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
      if (!hasGenTau) {
	
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
      if (DEBUG) std::cout << "\n=== Build Tau Candidate" << endl;


      // Build a tau candidate from tracks and EGs
      L1TkEGParticle newTauCandidate(trackTauCandidates[i], EGcluster, genTau, hasGenTau);
      if (DEBUG) { //marina
	cout << "Constructed a new tau candidate with a track invariant mass of " << newTauCandidate.GetTrackInvMass() 
	     << ", an EG invariant mass of " << newTauCandidate.GetEGInvMass();
	if (hasGenTau) cout  << " and a generator tau with visible pT " << newTauCandidate.GetGenTauPt() << endl;
	else cout << " and does not have a matching generator tau" << endl;
      }

      newTauCandidate.SetShrinkingConeConst(maxDeltaR_const);
      newTauCandidate.SetSigConeMaxOpen(maxDeltaR_leadtrk);
      newTauCandidate.FindIsoConeTracks(TTTracks);
      newTauCandidate.FindSignalConeEGs(L1EGs);
      newTauCandidate.FindIsoConeEGs(L1EGs);
      
      L1TkEGTauCandidates.push_back(newTauCandidate);
      
    }
    
    // Sort L1TkEGCandidates according to ET
    //SortL1TkEGs();
    sort(L1TkEGTauCandidates.begin(), L1TkEGTauCandidates.end(), [](L1TkEGParticle& a, L1TkEGParticle& b) {return a.GetTotalEt() > b.GetTotalEt();});
    
    
    if (L1TkEGTauCandidates.size() >=1){
      h_ldgTkEG_ET -> Fill( L1TkEGTauCandidates.at(0).GetTotalEt());
    }
    
    
    //====================================================
    // Apply Isolation (Vertex, Relative)
    //====================================================
    if (DEBUG) std::cout << "\n=== Apply Isolation (Vertex, Relative)" << endl;

    // Clear containers
    L1TkEGTaus_RelIso.clear();
    L1TkEGTaus_RelIsoLoose.clear();
    L1TkEGTaus_RelIsoTight.clear();
    L1TkEGTaus_VtxIso.clear();
    L1TkEGTaus_VtxIsoLoose.clear();
    L1TkEGTaus_VtxIsoTight.clear();
    
    // For-loop: All L1TkEGTauCandidates
    for (auto tkeg = L1TkEGTauCandidates.begin(); tkeg != L1TkEGTauCandidates.end(); tkeg++) {

      // ------ RELATIVE ISOLATION ------
      double relIso = tkeg->CalculateRelIso(TTTracks, maxDeltaZ_iso);
      
      // Fill relative isolation histogram
      h_TkEG_relIso->Fill( relIso );
      
      // Check Relative Isolation WPs
      bool bPassRelIso      = (relIso < relIso_WP); // orthogonal to VtxIso
      bool bPassRelIsoLoose = (relIso < 0.3);
      bool bPassRelIsoTight = (relIso < 0.1);
      
      // Fill containers with TkEGTaus
      if (bPassRelIso) L1TkEGTaus_RelIso.push_back(*tkeg);
      if (bPassRelIsoLoose) L1TkEGTaus_RelIsoLoose.push_back(*tkeg);
      if (bPassRelIsoTight) L1TkEGTaus_RelIsoTight.push_back(*tkeg);


      // ------ VERTEX ISOLATION ------
      
      double vtxIso = tkeg->CalculateVtxIso(TTTracks);
      
      // Fill vertex isolation histogram
      h_TkEG_vtxIso                 -> Fill( vtxIso );
      
      // Check Vertex Isolation WPs
      bool bPassVtxIso      = (vtxIso > vtxIso_WP); // orthogona1 to RelIso
      bool bPassVtxIsoLoose = (vtxIso > 0.2);
      bool bPassVtxIsoTight = (vtxIso > 1.0);
      
      // Fill containers with TkEGTaus
      if (bPassVtxIso) L1TkEGTaus_VtxIso.push_back(*tkeg);
      if (bPassVtxIsoLoose) L1TkEGTaus_VtxIsoLoose.push_back(*tkeg);
      if (bPassVtxIsoTight) L1TkEGTaus_VtxIsoTight.push_back(*tkeg);

    }// For-loop: All L1TkEGTauCandidates

    if (L1TkEGTaus_RelIso.size() > 0) nEvtsRelIso++;
    if (L1TkEGTaus_VtxIso.size() > 0) nEvtsVtxIso++;
    if (L1TkEGTaus_RelIsoLoose.size() > 0) nEvtsRelIsoLoose++;
    if (L1TkEGTaus_VtxIsoLoose.size() > 0) nEvtsVtxIsoLoose++;
    if (L1TkEGTaus_RelIsoTight.size() > 0) nEvtsRelIsoTight++;
    if (L1TkEGTaus_VtxIsoTight.size() > 0) nEvtsVtxIsoTight++;

    // Fill histos with tau candidates multiplicity
    h_TkEG_N              -> Fill(L1TkEGTauCandidates.size());
    h_TkEG_RelIso_N       -> Fill(L1TkEGTaus_RelIso.size());
    h_TkEG_VtxIso_N       -> Fill(L1TkEGTaus_VtxIso.size());
    h_TkEG_RelIsoLoose_N  -> Fill(L1TkEGTaus_RelIsoLoose.size());
    h_TkEG_VtxIsoLoose_N  -> Fill(L1TkEGTaus_VtxIsoLoose.size());
    h_TkEG_RelIsoTight_N  -> Fill(L1TkEGTaus_RelIsoTight.size());
    h_TkEG_VtxIsoTight_N  -> Fill(L1TkEGTaus_VtxIsoTight.size());


    //====================================================
    //  Tau Candidates (Properties & Resolution)
    //====================================================
        
    if (DEBUG) std::cout << "\n=== Tau Candidates (Properties & Resolution)" << endl;

    // Tau Candidates Properties & Resolution
    //for (auto tkeg = L1TkEGTauCandidates.begin(); tkeg != L1TkEGTauCandidates.end(); tkeg++) {
    for (auto tkeg = L1TkEGTaus_VtxIsoTight.begin(); tkeg != L1TkEGTaus_VtxIsoTight.end(); tkeg++) {
      
      if (!tkeg->HasMatchingGenParticle() && (isMinBias == false) ) continue;
      
      float tkeg_eta       = tkeg->GetTotalP4().Eta();
      vector<EG> clustEGs  = tkeg->GetEGs();
      int NEGs             = tkeg->GetEGs().size();
      int NTracks          = tkeg->GetTracks().size();
      TTTrack leadingTrack = tkeg->GetLeadingTrack();
      GenParticle        matchGenP = tkeg->GetMatchingGenParticle();
      vector<GenParticle> neuDaugh = matchGenP.finalDaughtersNeutral();
      

      h_TkEG_Pt      -> Fill (tkeg->GetTotalPt());
      h_TkEG_ET      -> Fill (tkeg->GetTotalEt());
      h_TkEG_Eta     -> Fill (tkeg->GetTotalP4().Eta());
      h_TkEG_Phi     -> Fill (tkeg->GetTotalP4().Phi());
      h_TkEG_InvMass -> Fill (tkeg->GetTotalP4().M());
      h_TkEG_NEGs    -> Fill (NEGs);
      h_TkEG_NTks    -> Fill (NTracks);
      
      
      if ( IsWithinEtaRegion("Central", tkeg_eta) ) {
	h_TkEG_NEGs_C  -> Fill(NEGs);
      }
      else if ( IsWithinEtaRegion("Intermediate", tkeg_eta) ) {	
	h_TkEG_NEGs_I  -> Fill(NEGs);
	}
      else if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	  h_TkEG_NEGs_F  -> Fill(NEGs);
      }
      
    
      // h_TkEG_CHF     -> Fill (tkeg->GetTrackBasedPt() / tkeg->GetTotalP4().Pt());
      // h_TkEG_NHF     -> Fill ( (tkeg->GetTotalP4().Pt() - tkeg->GetTrackBasedPt() ) / tkeg->GetTotalP4().Pt());
      
      if (NEGs > 0) {
	h_TkEG_NEGs_withEGs -> Fill (NEGs);

	double CHF = tkeg->GetTrackBasedPt() / tkeg->GetTotalP4().Pt();
	double NHF = (tkeg->GetTotalP4().Pt() - tkeg->GetTrackBasedPt() ) / tkeg->GetTotalP4().Pt();
	
	tkeg->SetCHF(CHF);
	tkeg->SetNHF(NHF);
	
	h_TkEG_CHF_withNeutrals     -> Fill (CHF);
	h_TkEG_NHF_withNeutrals     -> Fill (NHF);


	for (vector<EG>::iterator eg = clustEGs.begin(); eg != clustEGs.end() ; eg++) {
	  GenParticle match_pion0;
	  double deltaR;
	  double match_dR = 999;
	  bool isMatched = false;

	  // For-loop: All pion0s in the event                                                                                                                         
	  for (vector<GenParticle>::iterator daugh = neuDaugh.begin() ; daugh != neuDaugh.end() ; daugh ++){

	    if (daugh->et() < 1.5) continue;

	    deltaR = auxTools_.DeltaR( eg->getEta(), eg->getPhi(), daugh->eta(), daugh->phi() );
	    
	    if (deltaR < match_dR) {
	      match_dR = deltaR;
	      match_pion0 = *daugh;
	    }
	    
	  } // For-loop: All pion0s in the event                                                                                                                    
	  
	  // Matching of clustered EG
	  if (match_dR < 0.1) {
	    isMatched = true;
	  }

	  h_TkEG_clustEGs_MCMatch -> Fill(isMatched);

	  // Plot the quality bit of the matched and non Matched candidates
	  if (isMatched)
	    h_TkEG_clustEGs_Matched_HwQual    -> Fill( eg->getHwQual() );
	  else
	    h_TkEG_clustEGs_nonMatched_HwQual -> Fill( eg->getHwQual() );
	  
	  if (!isMatched) continue;
	  
	  h_TkEG_clustEGs_dET_matchPion0 ->  Fill( eg->getEt() - match_pion0.et() );
	  h_TkEG_clustEGs_ETResolution   ->  Fill( (eg->getEt() - match_pion0.et())/match_pion0.et());
	  
	}
      }
          
      // Isolation annulus properties 
      vector<TTTrack> isoTracks = tkeg->GetIsoConeTracks();
      vector<EG>      signalEGs = tkeg->GetSignalConeEGs();
      vector<EG>      isoEGs    = tkeg->GetIsoConeEGs();
      
      h_TkEG_isoTracks_Multiplicity -> Fill( isoTracks.size() );
      h_TkEG_signalEGs_Multiplicity -> Fill( signalEGs.size() );
      h_TkEG_isoEGs_Multiplicity    -> Fill( isoEGs.size() );
      
      
      TLorentzVector p4sum;
      // For-loop: All isoTracks
      for (auto tk = isoTracks.begin(); tk != isoTracks.end(); tk++) {
	p4sum += tk->p4();
      }

      if (isoTracks.size() > 0) {
        h_TkEG_isoTracks_InvMass -> Fill( p4sum.M() );
	h_TkEG_isoTracks_Et      -> Fill( p4sum.Et());
	h_TkEG_isoTracks_Eta     -> Fill( p4sum.Eta());
	h_TkEG_DonutRatio        -> Fill( GetDonutRatio(*tkeg, isoTracks, false) );
      }

      //if (tkeg->HasMatchingGenParticle()){ //marina
      
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

      // Charged Resolution
      int ChargedResolution = (tkeg->GetTracks().size() - tkeg->GetMatchingGenParticle().finalDaughtersCharged().size());
      h_TkEG_ChargedResolution->Fill(ChargedResolution);

      // Neutrals Resolution (EGs-GenPi0)
      int NeutralsResolution = (tkeg->GetEGs().size() - tkeg->GetMatchingGenParticle().finalDaughtersNeutral().size()); 
      h_TkEG_NeutralsResolution->Fill(NeutralsResolution); 
      

      // Studies for underestimating neutrals multiplicity
      if (NeutralsResolution < 0) {

	// GenParticle        matchGenP = tkeg->GetMatchingGenParticle();
	// vector<GenParticle> neuDaugh = matchGenP.finalDaughtersNeutral();

	h_TkEG_PoorNeuResol_NeuMultiplicity -> Fill( neuDaugh.size() );


	for (vector<GenParticle>::iterator daugh = neuDaugh.begin() ; daugh != neuDaugh.end() ; daugh ++){

	  float dR = auxTools_.DeltaR(daugh->eta(), daugh->phi(), matchGenP.p4vis().Eta(), matchGenP.p4vis().Phi());
	  float dEta = auxTools_.DeltaEta( daugh->eta(), matchGenP.p4vis().Eta() );
	  float dPhi = auxTools_.DeltaPhi( daugh->phi(), matchGenP.p4vis().Phi() );
	  
	  h_TkEG_PoorNeuResol_dR_Pi0_visTau   -> Fill ( dR );
	  h_TkEG_PoorNeuResol_dEta_Pi0_visTau -> Fill ( dEta);
	  h_TkEG_PoorNeuResol_dPhi_Pi0_visTau -> Fill ( dPhi);
	  
	  h_TkEG_PoorNeuResol_Pi0_ET -> Fill ( daugh->et() );

	  // Check for the closest EG if there are EGs in the event
	  if (L1EGs.size() == 0) continue;

	  double deltaR;
	  double min_dR = 999;
	  EG closestEG;
	  
	  for (vector<EG>::iterator eg = L1EGs.begin(); eg != L1EGs.end() ; eg++) {
	    
	    deltaR = auxTools_.DeltaR( eg->getEta(), eg->getPhi(), daugh->eta(), daugh->phi() );
	    
	    if (deltaR < min_dR)  {
	      min_dR = deltaR; 
	      closestEG = *eg;
	    }
	  }
	  
	  float min_dR_seed = auxTools_.DeltaR( closestEG.getEta(), closestEG.getPhi(), leadingTrack.getEta(), leadingTrack.getPhi() );

	  // Fill histos
	  h_TkEG_PoorNeuResol_dRmin_Pi0_EG         -> Fill( min_dR );
	  h_TkEG_PoorNeuResol_dRmin_Seed_closestEG -> Fill( min_dR_seed);
	  h_TkEG_PoorNeuResol_Pi0_closestEG_ET     -> Fill( closestEG.getEt() );
	  h_TkEG_PoorNeuResol_Pi0_closestEG_ET_Vs_dRmin_Pi0_EG -> Fill( closestEG.getEt() , min_dR );
	  h_TkEG_PoorNeuResol_Pi0_ET_Vs_closestEG_ET           -> Fill( daugh->et(), closestEG.getEt() );
	  
	}
      }
      
      // Resolution with and without clustered EGs
      if (tkeg->GetEGs().size() > 0) {
	h_TkEG_EtResolution_withEGs  -> Fill(ETresolution);
      }
      else {
	h_TkEG_EtResolution_noEGs  -> Fill(ETresolution);
      }

      // Resolution in different eta regions
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
	
	if (tkeg->GetEGs().size() > 0) {
	  h_TkEG_PtResolution_F_withEGs  -> Fill(pTresolution);
	  h_TkEG_EtResolution_F_withEGs  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_F_withEGs -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_F_withEGs -> Fill(PhiResolution);
	  
	  if (tkeg_eta >= 0){
	    h_TkEG_PtResolution_F_withEGs_posEta -> Fill(pTresolution);
	    h_TkEG_EtResolution_F_withEGs_posEta  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_F_withEGs_posEta -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_F_withEGs_posEta -> Fill(PhiResolution);
	  }
	  else{
	    h_TkEG_PtResolution_F_withEGs_negEta -> Fill(pTresolution);
	    h_TkEG_EtResolution_F_withEGs_negEta  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_F_withEGs_negEta -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_F_withEGs_negEta -> Fill(PhiResolution);
	  }
	  
	}
	else {
	  h_TkEG_PtResolution_F_noEGs  -> Fill(pTresolution);
	  h_TkEG_EtResolution_F_noEGs  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_F_noEGs -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_F_noEGs -> Fill(PhiResolution);
	  
	  if (tkeg_eta >= 0){
	    h_TkEG_PtResolution_F_noEGs_posEta -> Fill(pTresolution);
	    h_TkEG_EtResolution_F_noEGs_posEta  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_F_noEGs_posEta -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_F_noEGs_posEta -> Fill(PhiResolution);
	  }
	  else{
	    h_TkEG_PtResolution_F_noEGs_negEta -> Fill(pTresolution);
	    h_TkEG_EtResolution_F_noEGs_negEta  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_F_noEGs_negEta -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_F_noEGs_negEta -> Fill(PhiResolution);
	  }
	}
	
	// Study the properties of poor ET resolution candidates in forward region
	if ((ETresolution > 0.4) && (ETresolution < 0.8)) {
	  h_PoorEtResolCand_InvMass -> Fill(tkeg->GetTotalP4().M());
	  h_PoorEtResolCand_RelIso  -> Fill(tkeg->GetRelIso());
	  h_PoorEtResolCand_VtxIso  -> Fill(tkeg->GetVtxIso());
	  h_PoorEtResolCand_CHF     -> Fill(tkeg->GetCHF());
	  h_PoorEtResolCand_IsoTracks_N -> Fill(tkeg->GetIsoConeTracks().size());
	  
	  vector <EG> EGs = tkeg->GetEGs();
	  for (auto eg = EGs.begin(); eg != EGs.end(); eg++) {
	    for (unsigned int i = 0; i<= tkeg->GetEGs().size(); i++){
	      
	      deltaR = auxTools_.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), CorrectedEta(eg->getEta(), leadingTrack.getZ0()), eg->getPhi());
	      
	      h_PoorEtResolCand_dR_EG_Seed -> Fill(deltaR);
	    }
	  }
	}
	else{
	  h_GoodEtResolCand_InvMass -> Fill(tkeg->GetTotalP4().M());
	  h_GoodEtResolCand_RelIso  -> Fill(tkeg->GetRelIso());
	  h_GoodEtResolCand_VtxIso  -> Fill(tkeg->GetVtxIso());
	  h_GoodEtResolCand_CHF     -> Fill(tkeg->GetCHF());
	  h_GoodEtResolCand_IsoTracks_N -> Fill(tkeg->GetIsoConeTracks().size());
	  
	  vector <EG> EGs = tkeg->GetEGs();
	  for (auto eg = EGs.begin(); eg != EGs.end(); eg++) {
	    for (unsigned int i = 0; i<= tkeg->GetEGs().size(); i++){
	      leadingTrack = tkeg->GetLeadingTrack();
	      
	      deltaR = auxTools_.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), CorrectedEta(eg->getEta(), leadingTrack.getZ0()), eg->getPhi());
	      
	      h_GoodEtResolCand_dR_EG_Seed -> Fill(deltaR);
	    }
	  }
	}
      }

      // Daughters' multiplicity of matched GenParticle 
      int neuDaugh_N = tkeg->GetMatchingGenParticle().finalDaughtersNeutral().size();
      int chargedDaugh_N = tkeg->GetMatchingGenParticle().finalDaughtersCharged().size();

      // Resolution in case of zero or more than one neutral products 
      
      if (neuDaugh_N == 0) {
	h_TkEG_PtResolution_noNeutrals  -> Fill(pTresolution);
	h_TkEG_EtResolution_noNeutrals  -> Fill(ETresolution);
	h_TkEG_EtaResolution_noNeutrals -> Fill(EtaResolution);
	h_TkEG_PhiResolution_noNeutrals -> Fill(PhiResolution);
	h_TkEG_ChargedResolution_noNeutrals -> Fill(ChargedResolution);
	h_TkEG_NeutralsResolution_noNeutrals -> Fill(NeutralsResolution);

	if (tkeg->GetEGs().size() > 0) {
	  h_TkEG_EtResolution_noNeutrals_withEGs -> Fill(ETresolution);
	}
	else {
	  h_TkEG_EtResolution_noNeutrals_noEGs   -> Fill(ETresolution);
	}

		
	if ( IsWithinEtaRegion("Forward", tkeg_eta) ){
	  h_TkEG_PtResolution_noNeutrals_F  -> Fill(pTresolution);
	  h_TkEG_EtResolution_noNeutrals_F  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_noNeutrals_F -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_noNeutrals_F -> Fill(PhiResolution);
	  
	  if (tkeg->GetEGs().size() > 0) {
	    h_TkEG_PtResolution_noNeutrals_F_withEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_noNeutrals_F_withEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_noNeutrals_F_withEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_noNeutrals_F_withEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_noNeutrals_F_withEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_noNeutrals_F_withEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_noNeutrals_F_withEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_noNeutrals_F_withEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_noNeutrals_F_withEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_noNeutrals_F_withEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_noNeutrals_F_withEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_noNeutrals_F_withEGs_negEta -> Fill(PhiResolution);
	    }
	  }
	  else {
	    h_TkEG_PtResolution_noNeutrals_F_noEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_noNeutrals_F_noEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_noNeutrals_F_noEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_noNeutrals_F_noEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_noNeutrals_F_noEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_noNeutrals_F_noEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_noNeutrals_F_noEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_noNeutrals_F_noEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_noNeutrals_F_noEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_noNeutrals_F_noEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_noNeutrals_F_noEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_noNeutrals_F_noEGs_negEta -> Fill(PhiResolution);
	    }
	  }
	}
      }
      else {
	h_TkEG_PtResolution_withNeutrals  -> Fill(pTresolution);
	h_TkEG_EtResolution_withNeutrals  -> Fill(ETresolution);
	h_TkEG_EtaResolution_withNeutrals -> Fill(EtaResolution);
	h_TkEG_PhiResolution_withNeutrals -> Fill(PhiResolution);
	h_TkEG_ChargedResolution_withNeutrals -> Fill(ChargedResolution);
	h_TkEG_NeutralsResolution_withNeutrals -> Fill(NeutralsResolution);

	if (chargedDaugh_N == 1) {
	  h_TkEG_EtResolution_withNeutrals_1pr  -> Fill(ETresolution);
	}
	else if (chargedDaugh_N == 3) {
	  h_TkEG_EtResolution_withNeutrals_3pr  -> Fill(ETresolution);
	}

	if (neuDaugh_N == 1) {
	  h_TkEG_EtResolution_withNeutrals_1pion0 -> Fill(ETresolution);
	}
	else if (neuDaugh_N == 2) {
          h_TkEG_EtResolution_withNeutrals_2pion0 -> Fill(ETresolution);
        }
	else if (neuDaugh_N == 3) {
          h_TkEG_EtResolution_withNeutrals_3pion0 -> Fill(ETresolution);
	}
	else if (neuDaugh_N == 4) {
          h_TkEG_EtResolution_withNeutrals_4pion0 -> Fill(ETresolution);
	}

	
	if (tkeg->GetEGs().size() > 0) {
	  h_TkEG_EtResolution_withNeutrals_withEGs -> Fill(ETresolution);
	  if ( (tkeg->GetEGs().at(0).getEt() >= 0) && (tkeg->GetEGs().at(0).getEt() < 5) ) {
	    h_TkEG_EtResolution_withNeutrals_withEGs_0to5GeV ->Fill(ETresolution);
	    }
	  else if ( (tkeg->GetEGs().at(0).getEt() >= 5) && (tkeg->GetEGs().at(0).getEt() < 10) ) {
	    h_TkEG_EtResolution_withNeutrals_withEGs_5to10GeV ->Fill(ETresolution);
	  }
	  else if ( (tkeg->GetEGs().at(0).getEt() >= 10) && (tkeg->GetEGs().at(0).getEt() < 15) ) {
	    h_TkEG_EtResolution_withNeutrals_withEGs_10to15GeV ->Fill(ETresolution);
	  }
	  else if ( (tkeg->GetEGs().at(0).getEt() >= 15) && (tkeg->GetEGs().at(0).getEt() < 20) ) {
	    h_TkEG_EtResolution_withNeutrals_withEGs_15to20GeV ->Fill(ETresolution);
	  }
	  else if ( (tkeg->GetEGs().at(0).getEt() >= 20) && (tkeg->GetEGs().at(0).getEt() < 30) ) {
	    h_TkEG_EtResolution_withNeutrals_withEGs_20to30GeV ->Fill(ETresolution);
	  }
	  else if ( (tkeg->GetEGs().at(0).getEt() >= 30) && (tkeg->GetEGs().at(0).getEt() < 40) ) {
	    h_TkEG_EtResolution_withNeutrals_withEGs_30to40GeV ->Fill(ETresolution);
	  }
	  else if ( (tkeg->GetEGs().at(0).getEt() >= 40) && (tkeg->GetEGs().at(0).getEt() < 50) ) {
	    h_TkEG_EtResolution_withNeutrals_withEGs_40to50GeV ->Fill(ETresolution);
	  }

	}
	else {
	  h_TkEG_EtResolution_withNeutrals_noEGs   -> Fill(ETresolution);
	}

	if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	  h_TkEG_PtResolution_withNeutrals_F  -> Fill(pTresolution);
	  h_TkEG_EtResolution_withNeutrals_F  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_withNeutrals_F -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_withNeutrals_F -> Fill(PhiResolution);
	  
	  if (tkeg->GetEGs().size() > 0) {
	    h_TkEG_PtResolution_withNeutrals_F_withEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_withNeutrals_F_withEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_withNeutrals_F_withEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_withNeutrals_F_withEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_withNeutrals_F_withEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_withNeutrals_F_withEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_withNeutrals_F_withEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_withNeutrals_F_withEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_withNeutrals_F_withEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_withNeutrals_F_withEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_withNeutrals_F_withEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_withNeutrals_F_withEGs_negEta -> Fill(PhiResolution);
	    }
	    
	  }
	  else {
	    h_TkEG_PtResolution_withNeutrals_F_noEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_withNeutrals_F_noEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_withNeutrals_F_noEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_withNeutrals_F_noEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_withNeutrals_F_noEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_withNeutrals_F_noEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_withNeutrals_F_noEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_withNeutrals_F_noEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_withNeutrals_F_noEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_withNeutrals_F_noEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_withNeutrals_F_noEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_withNeutrals_F_noEGs_negEta -> Fill(PhiResolution);
	    }
	    
	  }
	}
      }
	
      // Resolution for 1- and 3-prong decays
 
      if (chargedDaugh_N == 1) {
	h_TkEG_PtResolution_1pr  -> Fill(pTresolution);
	h_TkEG_EtResolution_1pr  -> Fill(ETresolution);
	h_TkEG_EtaResolution_1pr -> Fill(EtaResolution);
	h_TkEG_PhiResolution_1pr -> Fill(PhiResolution);
	h_TkEG_ChargedResolution_1pr -> Fill(ChargedResolution);
	h_TkEG_NeutralsResolution_1pr -> Fill(NeutralsResolution);
	
	if (tkeg->GetEGs().size() > 0) {
	  h_TkEG_EtResolution_1pr_withEGs -> Fill(ETresolution);
	}
	else {
	  h_TkEG_EtResolution_1pr_noEGs   -> Fill(ETresolution);
	}

	if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	  h_TkEG_PtResolution_1pr_F  -> Fill(pTresolution);
	  h_TkEG_EtResolution_1pr_F  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_1pr_F -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_1pr_F -> Fill(PhiResolution);
	  
	  if (tkeg->GetEGs().size() > 0) {
	    h_TkEG_PtResolution_1pr_F_withEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_1pr_F_withEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_1pr_F_withEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_1pr_F_withEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_1pr_F_withEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_1pr_F_withEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_1pr_F_withEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_1pr_F_withEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_1pr_F_withEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_1pr_F_withEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_1pr_F_withEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_1pr_F_withEGs_negEta -> Fill(PhiResolution);
	    }
	    
	  }
	  else {
	    h_TkEG_PtResolution_1pr_F_noEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_1pr_F_noEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_1pr_F_noEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_1pr_F_noEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_1pr_F_noEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_1pr_F_noEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_1pr_F_noEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_1pr_F_noEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_1pr_F_noEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_1pr_F_noEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_1pr_F_noEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_1pr_F_noEGs_negEta -> Fill(PhiResolution);
	    }
	    
	  }
	}
      }
      else if (chargedDaugh_N == 3){
	h_TkEG_PtResolution_3pr  -> Fill(pTresolution);
	h_TkEG_EtResolution_3pr  -> Fill(ETresolution);
	h_TkEG_EtaResolution_3pr -> Fill(EtaResolution);
	h_TkEG_PhiResolution_3pr -> Fill(PhiResolution);
	h_TkEG_ChargedResolution_3pr -> Fill(ChargedResolution);
	h_TkEG_NeutralsResolution_3pr -> Fill(NeutralsResolution);
	
	if (tkeg->GetEGs().size() > 0) {
	  h_TkEG_EtResolution_3pr_withEGs -> Fill(ETresolution);
	}
	else {
	  h_TkEG_EtResolution_3pr_noEGs   -> Fill(ETresolution);
	}

	if ( IsWithinEtaRegion("Forward", tkeg_eta) ) {
	  h_TkEG_PtResolution_3pr_F  -> Fill(pTresolution);
	  h_TkEG_EtResolution_3pr_F  -> Fill(ETresolution);
	  h_TkEG_EtaResolution_3pr_F -> Fill(EtaResolution);
	  h_TkEG_PhiResolution_3pr_F -> Fill(PhiResolution);
	  
	  if (tkeg->GetEGs().size() > 0) {
	    h_TkEG_PtResolution_3pr_F_withEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_3pr_F_withEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_3pr_F_withEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_3pr_F_withEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_3pr_F_withEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_3pr_F_withEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_3pr_F_withEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_3pr_F_withEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_3pr_F_withEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_3pr_F_withEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_3pr_F_withEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_3pr_F_withEGs_negEta -> Fill(PhiResolution);
	    }
	    
	  }
	  else {
	    h_TkEG_PtResolution_3pr_F_noEGs  -> Fill(pTresolution);
	    h_TkEG_EtResolution_3pr_F_noEGs  -> Fill(ETresolution);
	    h_TkEG_EtaResolution_3pr_F_noEGs -> Fill(EtaResolution);
	    h_TkEG_PhiResolution_3pr_F_noEGs -> Fill(PhiResolution);
	    
	    if (tkeg_eta >= 0){
	      h_TkEG_PtResolution_3pr_F_noEGs_posEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_3pr_F_noEGs_posEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_3pr_F_noEGs_posEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_3pr_F_noEGs_posEta -> Fill(PhiResolution);
	    }
	    else{
	      h_TkEG_PtResolution_3pr_F_noEGs_negEta -> Fill(pTresolution);
	      h_TkEG_EtResolution_3pr_F_noEGs_negEta  -> Fill(ETresolution);
	      h_TkEG_EtaResolution_3pr_F_noEGs_negEta -> Fill(EtaResolution);
	      h_TkEG_PhiResolution_3pr_F_noEGs_negEta -> Fill(PhiResolution);
	    }
	  }
	}
      }	
    }

        
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (DEBUG) std::cout << "\n=== Fill Turn-On histograms" << endl;

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
    
    FillTurnOn_Numerator_(L1TkEGTauCandidates     , 25.0, hTkEG_TurnOn25, hTkEG_TurnOn25_1pr, hTkEG_TurnOn25_3pr, hTkEG_TurnOn25_withNeutrals, hTkEG_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIso , 25.0, hVtxIso_TurnOn25, hVtxIso_TurnOn25_1pr, hVtxIso_TurnOn25_3pr, hVtxIso_TurnOn25_withNeutrals, hVtxIso_TurnOn25_noNeutrals); 
    FillTurnOn_Numerator_(L1TkEGTaus_RelIso , 25.0, hRelIso_TurnOn25, hRelIso_TurnOn25_1pr, hRelIso_TurnOn25_3pr, hRelIso_TurnOn25_withNeutrals, hRelIso_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIsoLoose , 25.0, hVtxIsoLoose_TurnOn25, hVtxIsoLoose_TurnOn25_1pr, hVtxIsoLoose_TurnOn25_3pr, hVtxIsoLoose_TurnOn25_withNeutrals, hVtxIsoLoose_TurnOn25_noNeutrals); 
    FillTurnOn_Numerator_(L1TkEGTaus_RelIsoLoose , 25.0, hRelIsoLoose_TurnOn25, hRelIsoLoose_TurnOn25_1pr, hRelIsoLoose_TurnOn25_3pr, hRelIsoLoose_TurnOn25_withNeutrals, hRelIsoLoose_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIsoTight , 25.0, hVtxIsoTight_TurnOn25, hVtxIsoTight_TurnOn25_1pr, hVtxIsoTight_TurnOn25_3pr, hVtxIsoTight_TurnOn25_withNeutrals, hVtxIsoTight_TurnOn25_noNeutrals); 
    FillTurnOn_Numerator_(L1TkEGTaus_RelIsoTight , 25.0, hRelIsoTight_TurnOn25, hRelIsoTight_TurnOn25_1pr, hRelIsoTight_TurnOn25_3pr, hRelIsoTight_TurnOn25_withNeutrals, hRelIsoTight_TurnOn25_noNeutrals);


    FillTurnOn_Numerator_(L1TkEGTauCandidates     , 50.0, hTkEG_TurnOn50, hTkEG_TurnOn50_1pr, hTkEG_TurnOn50_3pr, hTkEG_TurnOn50_withNeutrals, hTkEG_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIso , 50.0, hVtxIso_TurnOn50, hVtxIso_TurnOn50_1pr, hVtxIso_TurnOn50_3pr, hVtxIso_TurnOn50_withNeutrals, hVtxIso_TurnOn50_noNeutrals); 
    FillTurnOn_Numerator_(L1TkEGTaus_RelIso , 50.0, hRelIso_TurnOn50, hRelIso_TurnOn50_1pr, hRelIso_TurnOn50_3pr, hRelIso_TurnOn50_withNeutrals, hRelIso_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIsoLoose , 50.0, hVtxIsoLoose_TurnOn50, hVtxIsoLoose_TurnOn50_1pr, hVtxIsoLoose_TurnOn50_3pr, hVtxIsoLoose_TurnOn50_withNeutrals, hVtxIsoLoose_TurnOn50_noNeutrals); 
    FillTurnOn_Numerator_(L1TkEGTaus_RelIsoLoose , 50.0, hRelIsoLoose_TurnOn50, hRelIsoLoose_TurnOn50_1pr, hRelIsoLoose_TurnOn50_3pr, hRelIsoLoose_TurnOn50_withNeutrals, hRelIsoLoose_TurnOn50_noNeutrals);
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIsoTight , 50.0, hVtxIsoTight_TurnOn50, hVtxIsoTight_TurnOn50_1pr, hVtxIsoTight_TurnOn50_3pr, hVtxIsoTight_TurnOn50_withNeutrals, hVtxIsoTight_TurnOn50_noNeutrals); 
    FillTurnOn_Numerator_(L1TkEGTaus_RelIsoTight , 50.0, hRelIsoTight_TurnOn50, hRelIsoTight_TurnOn50_1pr, hRelIsoTight_TurnOn50_3pr, hRelIsoTight_TurnOn50_withNeutrals, hRelIsoTight_TurnOn50_noNeutrals);

    // Turn-ons for the best performing WP
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIsoTight , 25.0, hL1TkEGTaus_TurnOn25, hL1TkEGTaus_TurnOn25_1pr, hL1TkEGTaus_TurnOn25_3pr, hL1TkEGTaus_TurnOn25_withNeutrals, hL1TkEGTaus_TurnOn25_noNeutrals); 
    FillTurnOn_Numerator_(L1TkEGTaus_VtxIsoTight , 50.0, hL1TkEGTaus_TurnOn50, hL1TkEGTaus_TurnOn50_1pr, hL1TkEGTaus_TurnOn50_3pr, hL1TkEGTaus_TurnOn50_withNeutrals, hL1TkEGTaus_TurnOn50_noNeutrals); 


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Rates and efficiencies for single tau
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (DEBUG) std::cout << "\n=== Rates and efficiencies for single tau" << endl;

    FillSingleTau_(L1TkEGTauCandidates    , hTkEG_Rate  , hTkEG_Eff  );
    FillSingleTau_(L1TkEGTauCandidates    , hTkEG_Rate_C, hTkEG_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkEGTauCandidates    , hTkEG_Rate_I, hTkEG_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkEGTauCandidates    , hTkEG_Rate_F, hTkEG_Eff_F, 1.6, 3.0); // 2.5 is max
    
    FillSingleTau_(L1TkEGTaus_VtxIso, hVtxIso_Rate  , hVtxIso_Eff);
    FillSingleTau_(L1TkEGTaus_VtxIso, hVtxIso_Rate_C, hVtxIso_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkEGTaus_VtxIso, hVtxIso_Rate_I, hVtxIso_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkEGTaus_VtxIso, hVtxIso_Rate_F, hVtxIso_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkEGTaus_RelIso, hRelIso_Rate  , hRelIso_Eff);
    FillSingleTau_(L1TkEGTaus_RelIso, hRelIso_Rate_C, hRelIso_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkEGTaus_RelIso, hRelIso_Rate_I, hRelIso_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkEGTaus_RelIso, hRelIso_Rate_F, hRelIso_Eff_F, 1.6, 3.0); // 2.5 is max
    
    FillSingleTau_(L1TkEGTaus_VtxIsoLoose, hVtxIsoLoose_Rate  , hVtxIsoLoose_Eff);
    FillSingleTau_(L1TkEGTaus_VtxIsoLoose, hVtxIsoLoose_Rate_C, hVtxIsoLoose_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkEGTaus_VtxIsoLoose, hVtxIsoLoose_Rate_I, hVtxIsoLoose_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkEGTaus_VtxIsoLoose, hVtxIsoLoose_Rate_F, hVtxIsoLoose_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkEGTaus_VtxIsoTight, hVtxIsoTight_Rate  , hVtxIsoTight_Eff);
    FillSingleTau_(L1TkEGTaus_VtxIsoTight, hVtxIsoTight_Rate_C, hVtxIsoTight_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkEGTaus_VtxIsoTight, hVtxIsoTight_Rate_I, hVtxIsoTight_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkEGTaus_VtxIsoTight, hVtxIsoTight_Rate_F, hVtxIsoTight_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkEGTaus_RelIsoLoose, hRelIsoLoose_Rate  , hRelIsoLoose_Eff);
    FillSingleTau_(L1TkEGTaus_RelIsoLoose, hRelIsoLoose_Rate_C, hRelIsoLoose_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkEGTaus_RelIsoLoose, hRelIsoLoose_Rate_I, hRelIsoLoose_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkEGTaus_RelIsoLoose, hRelIsoLoose_Rate_F, hRelIsoLoose_Eff_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(L1TkEGTaus_RelIsoTight, hRelIsoTight_Rate  , hRelIsoTight_Eff);
    FillSingleTau_(L1TkEGTaus_RelIsoTight, hRelIsoTight_Rate_C, hRelIsoTight_Eff_C, 0.0, 1.0);
    FillSingleTau_(L1TkEGTaus_RelIsoTight, hRelIsoTight_Rate_I, hRelIsoTight_Eff_I, 1.0, 1.6);
    FillSingleTau_(L1TkEGTaus_RelIsoTight, hRelIsoTight_Rate_F, hRelIsoTight_Eff_F, 1.6, 3.0); // 2.5 is max

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Rates and efficiencies for ditau 
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (DEBUG) std::cout << "\n=== Rates and efficiencies for ditau" << endl;

    ////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////
    ApplyDiTauZMatching(L1TkEGTauCandidates);
    ApplyDiTauZMatching(L1TkEGTaus_VtxIso); 
    ApplyDiTauZMatching(L1TkEGTaus_RelIso);
    ApplyDiTauZMatching(L1TkEGTaus_VtxIsoLoose);
    ApplyDiTauZMatching(L1TkEGTaus_VtxIsoTight);
    ApplyDiTauZMatching(L1TkEGTaus_RelIsoLoose);
    ApplyDiTauZMatching(L1TkEGTaus_RelIsoTight);

    FillDiTau_(L1TkEGTauCandidates, hDiTau_Rate_TkEG  , hDiTau_Eff_TkEG);
    FillDiTau_(L1TkEGTauCandidates, hDiTau_Rate_TkEG_C, hDiTau_Eff_TkEG_C, 0.0, 1.0);
    FillDiTau_(L1TkEGTauCandidates, hDiTau_Rate_TkEG_I, hDiTau_Eff_TkEG_I, 1.0, 1.6);
    FillDiTau_(L1TkEGTauCandidates, hDiTau_Rate_TkEG_F, hDiTau_Eff_TkEG_F, 1.6, 3.0);
    
    FillDiTau_(L1TkEGTaus_VtxIso, hDiTau_Rate_VtxIso  , hDiTau_Eff_VtxIso);
    FillDiTau_(L1TkEGTaus_VtxIso, hDiTau_Rate_VtxIso_C, hDiTau_Eff_VtxIso_C, 0.0, 1.0);
    FillDiTau_(L1TkEGTaus_VtxIso, hDiTau_Rate_VtxIso_I, hDiTau_Eff_VtxIso_I, 1.0, 1.6);
    FillDiTau_(L1TkEGTaus_VtxIso, hDiTau_Rate_VtxIso_F, hDiTau_Eff_VtxIso_F, 1.6, 3.0);

    FillDiTau_(L1TkEGTaus_RelIso, hDiTau_Rate_RelIso  , hDiTau_Eff_RelIso);
    FillDiTau_(L1TkEGTaus_RelIso, hDiTau_Rate_RelIso_C, hDiTau_Eff_RelIso_C, 0.0, 1.0);
    FillDiTau_(L1TkEGTaus_RelIso, hDiTau_Rate_RelIso_I, hDiTau_Eff_RelIso_I, 1.0, 1.6);
    FillDiTau_(L1TkEGTaus_RelIso, hDiTau_Rate_RelIso_F, hDiTau_Eff_RelIso_F, 1.6, 3.0);

    FillDiTau_(L1TkEGTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose  , hDiTau_Eff_VtxIsoLoose);
    FillDiTau_(L1TkEGTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose_C, hDiTau_Eff_VtxIsoLoose_C, 0.0, 1.0);
    FillDiTau_(L1TkEGTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose_I, hDiTau_Eff_VtxIsoLoose_I, 1.0, 1.6);
    FillDiTau_(L1TkEGTaus_VtxIsoLoose, hDiTau_Rate_VtxIsoLoose_F, hDiTau_Eff_VtxIsoLoose_F, 1.6, 3.0);

    FillDiTau_(L1TkEGTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight  , hDiTau_Eff_VtxIsoTight);
    FillDiTau_(L1TkEGTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight_C, hDiTau_Eff_VtxIsoTight_C, 0.0, 1.0);
    FillDiTau_(L1TkEGTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight_I, hDiTau_Eff_VtxIsoTight_I, 1.0, 1.6);
    FillDiTau_(L1TkEGTaus_VtxIsoTight, hDiTau_Rate_VtxIsoTight_F, hDiTau_Eff_VtxIsoTight_F, 1.6, 3.0);

    FillDiTau_(L1TkEGTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose  , hDiTau_Eff_RelIsoLoose);
    FillDiTau_(L1TkEGTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose_C, hDiTau_Eff_RelIsoLoose_C, 0.0, 1.0);
    FillDiTau_(L1TkEGTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose_I, hDiTau_Eff_RelIsoLoose_I, 1.0, 1.6);
    FillDiTau_(L1TkEGTaus_RelIsoLoose, hDiTau_Rate_RelIsoLoose_F, hDiTau_Eff_RelIsoLoose_F, 1.6, 3.0);

    FillDiTau_(L1TkEGTaus_RelIsoTight, hDiTau_Rate_RelIsoTight  , hDiTau_Eff_RelIsoTight);
    FillDiTau_(L1TkEGTaus_RelIsoTight, hDiTau_Rate_RelIsoTight_C, hDiTau_Eff_RelIsoTight_C, 0.0, 1.0);
    FillDiTau_(L1TkEGTaus_RelIsoTight, hDiTau_Rate_RelIsoTight_I, hDiTau_Eff_RelIsoTight_I, 1.0, 1.6);
    FillDiTau_(L1TkEGTaus_RelIsoTight, hDiTau_Rate_RelIsoTight_F, hDiTau_Eff_RelIsoTight_F, 1.6, 3.0);
    


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Rates and efficiencies for the best performing WP
    ////////////////////////////////////////////////////////////////////////////////////////////////    

    // Single-Tau 
    FillSingleTau_(L1TkEGTaus_VtxIsoTight, hL1TkEGTaus_SingleTau_Rate  , hL1TkEGTaus_SingleTau_Eff);
    
    // Di-Tau
    FillDiTau_(L1TkEGTaus_VtxIsoTight, hL1TkEGTaus_DiTau_Rate  , hL1TkEGTaus_DiTau_Eff);

    
    // Progress bar
    if (!DEBUG) auxTools_.ProgressBar(jentry, nEntries, 100, 100);
   
  } // For-loop: Entries


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Fill counters
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (DEBUG) std::cout << "\n=== Fill counters" << endl;


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
  h_Counters->SetBinContent( 10, nEvtsVtxIsoLoose);
  h_Counters->SetBinContent( 11, nEvtsVtxIsoTight);
  h_Counters->SetBinContent( 12, nEvtsRelIsoLoose);
  h_Counters->SetBinContent( 13, nEvtsRelIsoTight);

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
    

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Convert/Finalise histograms
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (DEBUG) std::cout << "\n=== Convert/Finalise histograms" << endl;

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

  // Turn-ons for the best performing WP
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn25, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn25_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn25_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn25_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);

  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn50, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn50_1pr, hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn50_3pr, hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hL1TkEGTaus_TurnOn50_noNeutrals, hMcHadronicTau_VisEt_noNeutrals);
  

    // SingleTau
  double N = nEntries;
  histoTools_.ConvertToRateHisto_1D(hTkEG_Rate  , N);
  histoTools_.ConvertToRateHisto_1D(hTkEG_Rate_C, N);
  histoTools_.ConvertToRateHisto_1D(hTkEG_Rate_I, N);
  histoTools_.ConvertToRateHisto_1D(hTkEG_Rate_F, N);
      
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

  FinaliseEffHisto_( hTkEG_Eff  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hTkEG_Eff_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hTkEG_Eff_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hTkEG_Eff_F, nEvtsWithMaxHTaus);
  
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
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_TkEG  , N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_TkEG_C, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_TkEG_I, N);
  histoTools_.ConvertToRateHisto_1D(hDiTau_Rate_TkEG_F, N);
  
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
  
  FinaliseEffHisto_( hDiTau_Eff_TkEG  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_TkEG_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_TkEG_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_( hDiTau_Eff_TkEG_F, nEvtsWithMaxHTaus);
  
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


  // Best performing WP

  // Single-Tau
  histoTools_.ConvertToRateHisto_1D( hL1TkEGTaus_SingleTau_Rate, N );
  FinaliseEffHisto_( hL1TkEGTaus_SingleTau_Eff, nEvtsWithMaxHTaus);
  
  // Di-Tau
  histoTools_.ConvertToRateHisto_1D( hL1TkEGTaus_DiTau_Rate, N );
  FinaliseEffHisto_( hL1TkEGTaus_DiTau_Eff,nEvtsWithMaxHTaus);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (DEBUG) std::cout << "\n=== Write the histograms to the file" << endl;
  
  WriteHistos_();
  auxTools_.StopwatchStop(5, "minutes", "Total Time");
  
}

//============================================================================
//void SortL1TkEGs()
//============================================================================
//{
//  sort(L1TkEGTauCandidates.begin(), L1TkEGTauCandidates.end(), [](L1TkEGParticle& a, L1TkEGParticle& b) {return a.GetTotalEt() > b.GetTotalEt();});
//  return;
//}


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
double TkEG::GetDonutRatio(L1TkEGParticle &L1TkEG,
			     vector<TTTrack> isoTTTracks, 
			     bool bUseCone)
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

  TTTrack seedTk = L1TkEG.GetLeadingTrack();  
  vector<TTTrack> isoConeTks; 
  isoConeTks = L1TkEG.GetIsoConeTracks();
  double sumPt_smallRAnnulus = 0.0;
  double sumPt_largeRAnnulus = 0.0;
  double smallR = L1TkEG.GetIsoConeMax();
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
void TkEG::BookHistos_(void)
//============================================================================
{

  // Binning
  const unsigned int nEt = 300;
  const float minEt = 0.0;
  const float maxEt = +300.0;

  const unsigned int nEta = 60;
  const float minEta = -3.0;
  const float maxEta = +3.0;

  const unsigned int nG = 1000;
  const float minG      =    0.0;
  const float maxG      =   10.0;


  const char* tEt   = ";E_{T} (GeV); Entries / %.0f GeV";
  const char* tEta  = ";#eta;Entries / %.2f";
  const char* tG    = ";#gamma;Entries / %.2f";


  // Event-Type Histograms
  histoTools_.BookHisto_1D(h_Counters, "Counters",  "", 13, 0.0, +13.0);
  h_Counters->GetXaxis()->SetBinLabel( 1, "All Evts");
  h_Counters->GetXaxis()->SetBinLabel( 2, "Evts");
  h_Counters->GetXaxis()->SetBinLabel( 3, "Seed Pt");
  h_Counters->GetXaxis()->SetBinLabel( 4, "Seed Eta");
  h_Counters->GetXaxis()->SetBinLabel( 5, "Seed ChiSq");
  h_Counters->GetXaxis()->SetBinLabel( 6, "Seed Stubs");
  h_Counters->GetXaxis()->SetBinLabel( 7, "Seed ldg in SigCone");
  h_Counters->GetXaxis()->SetBinLabel( 8, "RelIso");
  h_Counters->GetXaxis()->SetBinLabel( 9, "VtxIso");
  h_Counters->GetXaxis()->SetBinLabel(10, "VtxIso (L)");
  h_Counters->GetXaxis()->SetBinLabel(11, "VtxIso (T)");
  h_Counters->GetXaxis()->SetBinLabel(12, "RelIso (L)");
  h_Counters->GetXaxis()->SetBinLabel(13, "RelIso (T)");

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
  histoTools_.BookHisto_2D(h_genTau_chargedDaugh_PtLead_dRmax, "genTau_chargedDaugh_PtLead_dRmax", "; #DeltaR_{max} (#pi^{#pm}_{ldg}, #pi^{#pm}) ;p_{T} (#pi^{#pm}_{ldg}) (GeV)", 100, 0.0, 1.0, 100, 0.0, 100.0);

  // Pt of leading charged decay product  (leading pt>10GeV) Vs dRmax with  other neutral products
  histoTools_.BookHisto_2D(h_genTau_neutralDaugh_PtLead_dRmax, "genTau_neutralDaugh_PtLead_dRmax", "; #DeltaR_{max} (#pi^{#pm}_{ldg}, #pi^{0});p_{T} (#pi^{#pm}_{ldg}) (GeV)",  100, 0.0, 1.0, 100, 0.0, 100.0);

  // Et of pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Pion0_Et, "Pion0_Et", ";E_{T} (GeV); Particles / bin", 100, +0.0, +100.0);

  // Et of the photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Photons_Et, "Photons_Et", ";E_{T} (GeV); Particles / bin", 100, +0.0, +100.0);

  // dR of the two photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Photons_dR  , "Photons_dR"  , ";#Delta R; Entries / bin",   200,  0.0,   +2.0);

  // dPhi of the two photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Photons_dPhi  , "Photons_dPhi"  , ";#Delta#phi; Entries / bin",  200,  -1.0,  1.0);
  
  // dEta of the two photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_1D(h_Photons_dEta  , "Photons_dEta"  , ";#Delta#eta; Entries / bin",  400,  -2.0,  2.0);
  
  histoTools_.BookHisto_2D(h_Photons_dEtaVsdPhi, "Photons_dEtaVsdPhi", ";#Delta#eta;#Delta#phi; Entries / bin", 400,  -2.0,  2.0, 200,  -1.0,  1.0);
  // Et of pi0 Vs dR of the two photons from pi0 (1 prong decay, 1 pi0 cases)
  histoTools_.BookHisto_2D(h_Pion0Et_Vs_PhotonsDR, "Pion0Et_Vs_PhotonsDR", ";E_{T} (#pi^{0})  (GeV); #DeltaR (photon_{1}, photon_{2})", 100, 0.0, 100.0, 200, 0.0, 2.0);

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
  histoTools_.BookHisto_1D(h_leadTrks_Multiplicity, "leadTrks_Multiplicity", ";Number of lead tracks in event;Events / bin", 16, -0.5, +15.5);

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

  // Opening of signal cone 
  histoTools_.BookHisto_1D(h_SigCone_DeltaR, "SigCone_DeltaR", ";#DeltaR; Entries/ bin", 150, 0, 0.3 );

  // Pt of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_Pt, "trkClusters_Pt", ";p_{T} (GeV);Clusters / bin", 150, +0.0, +300.0);

  // Invariant mass of track clusters (before cut)
  histoTools_.BookHisto_1D(h_trkClusters_M_beforeCut, "trkClusters_M_beforeCut", ";Invariant mass;Clusters / bin", 40, 0.0, +4.0);

  // Invariant mass of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_M, "trkClusters_M", ";Invariant mass;Clusters / bin", 40, 0.0, +4.0);

  // Number of EGs                                                                                                                                                       
  histoTools_.BookHisto_1D(h_EGs_N, "EGs_Multiplicity", ";Number of EGs in the event ; entries / bin", 201, -0.5, +200.5);

  // MCmatched EGs Et 
  histoTools_.BookHisto_1D(h_EGs_MCmatched_Et, "EGs_MCmatched_Et", ";E_{T} (GeV) ;EGs / bin", 50, 0.0, +100.0);

  // EGs Properties
  histoTools_.BookHisto_1D(h_EGs_Et, "EGs_Et", ";E_{T} (GeV) ;EGs / bin", 200, 0.0, +100.0);
  histoTools_.BookHisto_1D(h_EGs_Eta, "EGs_Eta", tEta  , nEta  , minEta  , maxEta  );
  histoTools_.BookHisto_1D(h_EGs_Phi, "EGs_Phi", ";#phi (rads);EGs / bin", 36,  -3.15,  +3.15);
  // histoTools_.BookHisto_1D(h_EGs_IEta, "EGs_IEta", ";i#eta;EGs / bin", 70, -35, +35); 
  // histoTools_.BookHisto_1D(h_EGs_IPhi, "EGs_IPhi", ";i#phi;EGs / bin", 36,  0,  145);
  histoTools_.BookHisto_1D(h_EGs_HwQual, "EGs_HwQual", ";Quality bit; EGs / bin", 32, 0.5, 32.5); 
  histoTools_.BookHisto_2D(h_EGs_EtaVsEt, "EGs_EtaVsEt", ";#eta ;E_{T} (GeV); Entries/bin"  , nEta  , minEta  , maxEta, nEt  , minEt  , maxEt);

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

  // Invariant mass of track+eg clusters (before cut)
  histoTools_.BookHisto_1D(h_TkEGClusters_M_beforeCut, "TkEGClusters_M_beforeCut", ";Invariant mass;Clusters / bin", 40, 0.0, +4.0);


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
  
  
  histoTools_.BookHisto_1D(h_ldgTkEG_ET, "TkEG_ldgET", ";E_{T} (GeV); TkEGs / bin", 300, 0.0,  +300.0); 

  // DeltaR in MC matching of the lead track of tau candidates
  histoTools_.BookHisto_1D(hTkEG_DeltaRmatch  , "TkEG_DeltaRmatch"  , ";#Delta R;Particles / bin",   100,  0.0,   +1.0);
  
  // Visible Et of the matched generator taus
  histoTools_.BookHisto_1D(hTkEG_genVisEt, "TkEG_genVisEt", ";E_{T}^{vis} (gen) (GeV);Particles / bin", 40, 0.0,  +200.0);

  // Visible Et of the matched generator taus (when there is at least one EG clustered)
  histoTools_.BookHisto_1D(hTkEG_genVisEt_clustEG, "TkEG_genVisEt_clustEG", ";E_{T}^{vis} (gen) (GeV);Particles / bin", 40, 0.0,  +200.0);
  
  // Visible Pt of the matched generator taus (when there is at least one EG clustered)
  histoTools_.BookHisto_1D(hTkEG_genVisPt_clustEG, "TkEG_genVisPt_clustEG", ";p_{T}^{vis} (gen) (GeV);Particles / bin", 40, 0.0,  +200.0);

  // Multiplicity of tau candidates 
  histoTools_.BookHisto_1D(h_TkEG_N, "TkEG_N", ";#tau-candidate multiplicity;Entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_TkEG_RelIso_N, "TkEG_RelIso_N", ";#tau-candidate multiplicity;Entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_TkEG_VtxIso_N, "TkEG_VtxIso_N", ";#tau-candidate multiplicity;Entries / bin", 16, -0.5, +15.5); 
  histoTools_.BookHisto_1D(h_TkEG_RelIsoLoose_N, "TkEG_RelIsoLoose_N", ";#tau-candidate multiplicity;Entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_TkEG_VtxIsoLoose_N, "TkEG_VtxIsoLoose_N", ";#tau-candidate multiplicity;Entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_TkEG_RelIsoTight_N, "TkEG_RelIsoTight_N", ";#tau-candidate multiplicity;Entries / bin", 16, -0.5, +15.5);
  histoTools_.BookHisto_1D(h_TkEG_VtxIsoTight_N, "TkEG_VtxIsoTight_N", ";#tau-candidate multiplicity;Entries / bin", 16, -0.5, +15.5);

  // Properties of tau candidate 
  histoTools_.BookHisto_1D(h_TkEG_Pt, "TkEG_Pt", ";p_{T} (GeV);Entries / bin", 300, +0.0, +300.0);
  histoTools_.BookHisto_1D(h_TkEG_ET, "TkEG_ET", ";E_{T} (GeV);Entries / bin", 300, +0.0, +300.0);
  histoTools_.BookHisto_1D(h_TkEG_Eta, "TkEG_Eta", ";#eta;Entries / bin", 60, -3.0, +3.0);
  histoTools_.BookHisto_1D(h_TkEG_Phi, "TkEG_Phi", ";#phi (rads); Entries / bin", 21,  -3.15,  +3.15);
  histoTools_.BookHisto_1D(h_TkEG_InvMass, "TkEG_InvMass", ";Mass (GeV); Entries / bin", 40, 0.0, +4.0);
  histoTools_.BookHisto_1D(h_TkEG_NEGs, "TkEG_NEGs", ";N_{EGs};Entries / bin", 6,-0.5, 5.5);
  histoTools_.BookHisto_1D(h_TkEG_NEGs_withEGs, "TkEG_NEGs_withEGs", ";N_{EGs};Entries / bin", 6,-0.5, 5.5);
  histoTools_.BookHisto_1D(h_TkEG_NTks, "TkEG_NTks", ";N_{Tks};Entries / bin", 6,-0.5, 5.5);
  histoTools_.BookHisto_1D(h_TkEG_NEGs_C, "TkEG_NEGs_C", ";N_{EGs};Entries / bin", 6,-0.5, 5.5);
  histoTools_.BookHisto_1D(h_TkEG_NEGs_I, "TkEG_NEGs_I", ";N_{EGs};Entries / bin", 6,-0.5, 5.5);
  histoTools_.BookHisto_1D(h_TkEG_NEGs_F, "TkEG_NEGs_F", ";N_{EGs};Entries / bin", 6,-0.5, 5.5);
  histoTools_.BookHisto_1D(h_TkEG_CHF, "TkEG_CHF", ";CHF; Entries / bin", 200,  0.0,   +2.0);
  histoTools_.BookHisto_1D(h_TkEG_NHF, "TkEG_NHF", ";NHF; Entries / bin", 200,  0.0,   +2.0);
  histoTools_.BookHisto_1D(h_TkEG_CHF_withNeutrals, "TkEG_CHF_withNeutrals", ";CHF; Entries / bin", 100,  0.0,   +1.0);
  histoTools_.BookHisto_1D(h_TkEG_NHF_withNeutrals, "TkEG_NHF_withNeutrals", ";NHF; Entries / bin", 100,  0.0,   +1.0);

  histoTools_.BookHisto_1D(h_TkEG_clustEGs_MCMatch, "TkEG_clustEGs_MCMatch", ";;Entries / bin", 2, 0, 2);
  h_TkEG_clustEGs_MCMatch->GetXaxis()->SetBinLabel(1,"not MC Matched EG");
  h_TkEG_clustEGs_MCMatch->GetXaxis()->SetBinLabel(2,"MC Matched EG");

  histoTools_.BookHisto_1D(h_TkEG_clustEGs_Matched_HwQual, "TkEG_clustEGs_Matched_HwQual", ";Quality bit; Entries / bin", 32, 0.5, 32.5);
  histoTools_.BookHisto_1D(h_TkEG_clustEGs_nonMatched_HwQual, "TkEG_clustEGs_nonMatched_HwQual", ";Quality bit; Entries / bin", 32, 0.5, 32.5);
  histoTools_.BookHisto_1D(h_TkEG_clustEGs_dET_matchPion0, "TkEG_clustEGs_dET_matchPion0", ";#deltaE_{T}; Entries / bin", 100, -50, 50); 
  histoTools_.BookHisto_1D(h_TkEG_clustEGs_ETResolution, "TkEG_clustEGs_ETResolution", ";E_{T} resolution (GeV);Clusters / bin", 20000, -10.0, +10.0); 
  histoTools_.BookHisto_1D(h_TkEG_isoTracks_InvMass, "TkEG_isoTracks_InvMass", ";Invariant mass (iso-cone); Entries / bin", 5000, +0.0, +5.0);
  histoTools_.BookHisto_1D(h_TkEG_isoTracks_Multiplicity, "TkEG_isoTracks_Multiplicity", ";N_{iso-tks}; Entries / bin",11, -0.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_isoTracks_Et, "TkEG_isoTracks_Et", tEt  , nEt  , minEt  , maxEt  );
  histoTools_.BookHisto_1D(h_TkEG_isoTracks_Eta, "TkEG_isoTracks_Eta", tEta  , nEta  , minEta  , maxEta  );
  histoTools_.BookHisto_1D(h_TkEG_DonutRatio, "TkEG_DonutRatio", tG   , nG   , minG   , maxG   );
  histoTools_.BookHisto_1D(h_TkEG_signalEGs_Multiplicity, "TkEG_signalEGs_Multiplicity", ";N_{signal-EGs}; Entries / bin",11, -0.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_isoEGs_Multiplicity, "TkEG_isoEGs_Multiplicity", ";N_{iso-EGs}; Entries / bin",11, -0.5, 10.5);

  // Poor Et Resolution Candidates Properties
  histoTools_.BookHisto_1D(h_PoorEtResolCand_InvMass, "PoorEtResolCand_InvMass", ";Mass (GeV); Entries / bin", 40, 0.0, +4.0);
  histoTools_.BookHisto_1D(h_PoorEtResolCand_RelIso, "PoorEtResolCand_RelIso", ";Relative isolation;Clusters / bin", 100, 0.0, +5.0);
  histoTools_.BookHisto_1D(h_PoorEtResolCand_VtxIso, "PoorEtResolCand_VtxIso", ";Vertex isolation;Clusters / bin", 100, 0.0, +5.0);
  histoTools_.BookHisto_1D(h_PoorEtResolCand_CHF, "PoorEtResolCand_CHF", ";CHF; Entries / bin", 100,  0.0,   +1.0);
  histoTools_.BookHisto_1D(h_PoorEtResolCand_IsoTracks_N, "PoorEtResolCand_IsoTracks_N", ";N_{iso-tks}; Entries / bin",11, -0.5, 10.5);
  histoTools_.BookHisto_1D(h_PoorEtResolCand_dR_EG_Seed, "PoorEtResolCand_dR_EG_Seed", ";#DeltaR(#tau-seed, EG); Entries / bin",30,  0.0,   +0.3);

  histoTools_.BookHisto_1D(h_GoodEtResolCand_InvMass, "GoodEtResolCand_InvMass", ";Mass (GeV); Entries / bin", 40, 0.0, +4.0);
  histoTools_.BookHisto_1D(h_GoodEtResolCand_RelIso, "GoodEtResolCand_RelIso", ";Relative isolation;Clusters / bin", 100, 0.0, +5.0);
  histoTools_.BookHisto_1D(h_GoodEtResolCand_VtxIso, "GoodEtResolCand_VtxIso", ";Vertex isolation;Clusters / bin", 100, 0.0, +5.0);
  histoTools_.BookHisto_1D(h_GoodEtResolCand_CHF, "GoodEtResolCand_CHF", ";CHF; Entries / bin", 100,  0.0,   +1.0);
  histoTools_.BookHisto_1D(h_GoodEtResolCand_IsoTracks_N, "GoodEtResolCand_IsoTracks_N", ";N_{iso-tks}; Entries / bin",11, -0.5, 10.5);
  histoTools_.BookHisto_1D(h_GoodEtResolCand_dR_EG_Seed, "GoodEtResolCand_dR_EG_Seed", ";#DeltaR(#tau-seed, EG); Entries / bin",30,  0.0,   +0.3);

  // Poor Neutral Resolution candidates 
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_NeuMultiplicity, "TkEG_PoorNeuResol_NeuMultiplicity", ";N_{#pi^{0}}; Entries / bin", 11, -0.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_dR_Pi0_visTau, "TkEG_PoorNeuResol_dR_Pi0_visTau", ";#DeltaR (#pi^{0}, #tau_{vis})", 50, 0.0, 0.5);
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_dEta_Pi0_visTau, "TkEG_PoorNeuResol_dEta_Pi0_visTau", ";#Delta#eta (#pi^{0}, #tau_{vis})", 100,  -10.0,  10.0);
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_dPhi_Pi0_visTau, "TkEG_PoorNeuResol_dPhi_Pi0_visTau", ";#Delta#phi (#pi^{0}, #tau_{vis})",  21,  -3.15,  +3.15);
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_Pi0_ET, "TkEG_PoorNeuResol_Pi0_ET", ";E_{T}^{#pi^{0}};Entries/bin", nEt, minEt, maxEt);
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_Pi0_closestEG_ET, "TkEG_PoorNeuResol_Pi0_closestEG_ET", ";E_{T}^{EG}; Entries/bin", nEt, minEt, maxEt);
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_dRmin_Pi0_EG, "TkEG_PoorNeuResol_dRmin_Pi0_EG", ";#DeltaR_{min} (#pi^{0}, EG)", 120, 0.0, 6.0);
  histoTools_.BookHisto_1D(h_TkEG_PoorNeuResol_dRmin_Seed_closestEG, "TkEG_PoorNeuResol_dRmin_Seed_closestEG", ";#DeltaR_{min} (#tau-seed, EG)", 120, 0.0, 6.0);
  histoTools_.BookHisto_2D(h_TkEG_PoorNeuResol_Pi0_closestEG_ET_Vs_dRmin_Pi0_EG, "TkEG_PoorNeuResol_Pi0_closestEG_ET_Vs_dRmin_Pi0_EG",";E_{T}^{EG};#DeltaR_{min} (#pi^{0}, EG); Entries/bin", nEt, minEt, maxEt, 120, 0.0, 6.0);
  histoTools_.BookHisto_2D(h_TkEG_PoorNeuResol_Pi0_ET_Vs_closestEG_ET, "TkEG_PoorNeuResol_Pi0_ET_Vs_closestEG_ET", ";E_{T}^{#pi^{0}};E_{T}^{EG};Entries/bin", nEt, minEt, maxEt, nEt, minEt, maxEt);

  // Pt resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_PtResolution, "TkEG_PtResolution", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_C, "TkEG_PtResolution_C", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_I, "TkEG_PtResolution_I", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F, "TkEG_PtResolution_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals, "TkEG_PtResolution_noNeutrals", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals, "TkEG_PtResolution_withNeutrals", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr, "TkEG_PtResolution_1pr", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr, "TkEG_PtResolution_3pr", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F_withEGs, "TkEG_PtResolution_F_withEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F_withEGs_posEta, "TkEG_PtResolution_F_withEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F_withEGs_negEta, "TkEG_PtResolution_F_withEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F_noEGs, "TkEG_PtResolution_F_noEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F_noEGs_posEta, "TkEG_PtResolution_F_noEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_F_noEGs_negEta, "TkEG_PtResolution_F_noEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F, "TkEG_PtResolution_noNeutrals_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F, "TkEG_PtResolution_withNeutrals_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F, "TkEG_PtResolution_1pr_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F, "TkEG_PtResolution_3pr_F", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F_withEGs, "TkEG_PtResolution_noNeutrals_F_withEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F_withEGs, "TkEG_PtResolution_withNeutrals_F_withEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F_withEGs, "TkEG_PtResolution_1pr_F_withEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F_withEGs, "TkEG_PtResolution_3pr_F_withEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F_noEGs, "TkEG_PtResolution_noNeutrals_F_noEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F_noEGs, "TkEG_PtResolution_withNeutrals_F_noEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F_noEGs, "TkEG_PtResolution_1pr_F_noEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F_noEGs, "TkEG_PtResolution_3pr_F_noEGs", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // eta > 0
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F_withEGs_posEta, "TkEG_PtResolution_noNeutrals_F_withEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F_withEGs_posEta, "TkEG_PtResolution_withNeutrals_F_withEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F_withEGs_posEta, "TkEG_PtResolution_1pr_F_withEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F_withEGs_posEta, "TkEG_PtResolution_3pr_F_withEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F_noEGs_posEta, "TkEG_PtResolution_noNeutrals_F_noEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F_noEGs_posEta, "TkEG_PtResolution_withNeutrals_F_noEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F_noEGs_posEta, "TkEG_PtResolution_1pr_F_noEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F_noEGs_posEta, "TkEG_PtResolution_3pr_F_noEGs_posEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // eta < 0
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F_withEGs_negEta, "TkEG_PtResolution_noNeutrals_F_withEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F_withEGs_negEta, "TkEG_PtResolution_withNeutrals_F_withEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F_withEGs_negEta, "TkEG_PtResolution_1pr_F_withEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F_withEGs_negEta, "TkEG_PtResolution_3pr_F_withEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PtResolution_noNeutrals_F_noEGs_negEta, "TkEG_PtResolution_noNeutrals_F_noEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_withNeutrals_F_noEGs_negEta, "TkEG_PtResolution_withNeutrals_F_noEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_1pr_F_noEGs_negEta, "TkEG_PtResolution_1pr_F_noEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PtResolution_3pr_F_noEGs_negEta, "TkEG_PtResolution_3pr_F_noEGs_negEta", ";p_{T} resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // Et resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_EtResolution, "TkEG_EtResolution", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_C, "TkEG_EtResolution_C", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_I, "TkEG_EtResolution_I", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F, "TkEG_EtResolution_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals, "TkEG_EtResolution_noNeutrals", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals, "TkEG_EtResolution_withNeutrals", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr, "TkEG_EtResolution_1pr", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr, "TkEG_EtResolution_3pr", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withEGs, "TkEG_EtResolution_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noEGs, "TkEG_EtResolution_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_1pr, "TkEG_EtResolution_withNeutrals_1pr", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_3pr, "TkEG_EtResolution_withNeutrals_3pr", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_1pion0, "TkEG_EtResolution_withNeutrals_1pion0", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_2pion0, "TkEG_EtResolution_withNeutrals_2pion0", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_3pion0, "TkEG_EtResolution_withNeutrals_3pion0", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_4pion0, "TkEG_EtResolution_withNeutrals_4pion0", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
    
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_withEGs, "TkEG_EtResolution_noNeutrals_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs, "TkEG_EtResolution_withNeutrals_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_withEGs, "TkEG_EtResolution_1pr_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_withEGs, "TkEG_EtResolution_3pr_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs_0to5GeV,  "TkEG_EtResolution_withNeutrals_withEGs_0to5GeV", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs_5to10GeV,  "TkEG_EtResolution_withNeutrals_withEGs_5to10GeV", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs_10to15GeV,  "TkEG_EtResolution_withNeutrals_withEGs_10to15GeV", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs_15to20GeV,  "TkEG_EtResolution_withNeutrals_withEGs_15to20GeV", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs_20to30GeV,  "TkEG_EtResolution_withNeutrals_withEGs_20to30GeV", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs_30to40GeV,  "TkEG_EtResolution_withNeutrals_withEGs_30to40GeV", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_withEGs_40to50GeV,  "TkEG_EtResolution_withNeutrals_withEGs_40to50GeV", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);


  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_noEGs, "TkEG_EtResolution_noNeutrals_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_noEGs, "TkEG_EtResolution_withNeutrals_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_noEGs, "TkEG_EtResolution_1pr_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_noEGs, "TkEG_EtResolution_3pr_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F_withEGs, "TkEG_EtResolution_F_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F_withEGs_posEta, "TkEG_EtResolution_F_withEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F_withEGs_negEta, "TkEG_EtResolution_F_withEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F_noEGs, "TkEG_EtResolution_F_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F_noEGs_posEta, "TkEG_EtResolution_F_noEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_F_noEGs_negEta, "TkEG_EtResolution_F_noEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F, "TkEG_EtResolution_noNeutrals_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F, "TkEG_EtResolution_withNeutrals_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F, "TkEG_EtResolution_1pr_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F, "TkEG_EtResolution_3pr_F", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F_withEGs, "TkEG_EtResolution_noNeutrals_F_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F_withEGs, "TkEG_EtResolution_withNeutrals_F_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F_withEGs, "TkEG_EtResolution_1pr_F_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F_withEGs, "TkEG_EtResolution_3pr_F_withEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F_noEGs, "TkEG_EtResolution_noNeutrals_F_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F_noEGs, "TkEG_EtResolution_withNeutrals_F_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F_noEGs, "TkEG_EtResolution_1pr_F_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F_noEGs, "TkEG_EtResolution_3pr_F_noEGs", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // eta > 0
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F_withEGs_posEta, "TkEG_EtResolution_noNeutrals_F_withEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F_withEGs_posEta, "TkEG_EtResolution_withNeutrals_F_withEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F_withEGs_posEta, "TkEG_EtResolution_1pr_F_withEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F_withEGs_posEta, "TkEG_EtResolution_3pr_F_withEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F_noEGs_posEta, "TkEG_EtResolution_noNeutrals_F_noEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F_noEGs_posEta, "TkEG_EtResolution_withNeutrals_F_noEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F_noEGs_posEta, "TkEG_EtResolution_1pr_F_noEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F_noEGs_posEta, "TkEG_EtResolution_3pr_F_noEGs_posEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  
    // eta < 0
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F_withEGs_negEta, "TkEG_EtResolution_noNeutrals_F_withEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F_withEGs_negEta, "TkEG_EtResolution_withNeutrals_F_withEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F_withEGs_negEta, "TkEG_EtResolution_1pr_F_withEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F_withEGs_negEta, "TkEG_EtResolution_3pr_F_withEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_noNeutrals_F_noEGs_negEta, "TkEG_EtResolution_noNeutrals_F_noEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_withNeutrals_F_noEGs_negEta, "TkEG_EtResolution_withNeutrals_F_noEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_1pr_F_noEGs_negEta, "TkEG_EtResolution_1pr_F_noEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtResolution_3pr_F_noEGs_negEta, "TkEG_EtResolution_3pr_F_noEGs_negEta", ";Et resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);


  // Eta resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution, "TkEG_EtaResolution", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_C, "TkEG_EtaResolution_C", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_I, "TkEG_EtaResolution_I", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F, "TkEG_EtaResolution_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals, "TkEG_EtaResolution_noNeutrals", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals, "TkEG_EtaResolution_withNeutrals", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr, "TkEG_EtaResolution_1pr", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr, "TkEG_EtaResolution_3pr", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F_withEGs, "TkEG_EtaResolution_F_withEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F_withEGs_posEta, "TkEG_EtaResolution_F_withEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F_withEGs_negEta, "TkEG_EtaResolution_F_withEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F_noEGs, "TkEG_EtaResolution_F_noEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F_noEGs_posEta, "TkEG_EtaResolution_F_noEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_F_noEGs_negEta, "TkEG_EtaResolution_F_noEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F, "TkEG_EtaResolution_noNeutrals_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F, "TkEG_EtaResolution_withNeutrals_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F, "TkEG_EtaResolution_1pr_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F, "TkEG_EtaResolution_3pr_F", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F_withEGs, "TkEG_EtaResolution_noNeutrals_F_withEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F_withEGs, "TkEG_EtaResolution_withNeutrals_F_withEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F_withEGs, "TkEG_EtaResolution_1pr_F_withEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F_withEGs, "TkEG_EtaResolution_3pr_F_withEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F_noEGs, "TkEG_EtaResolution_noNeutrals_F_noEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F_noEGs, "TkEG_EtaResolution_withNeutrals_F_noEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F_noEGs, "TkEG_EtaResolution_1pr_F_noEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F_noEGs, "TkEG_EtaResolution_3pr_F_noEGs", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // eta > 0
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F_withEGs_posEta, "TkEG_EtaResolution_noNeutrals_F_withEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F_withEGs_posEta, "TkEG_EtaResolution_withNeutrals_F_withEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F_withEGs_posEta, "TkEG_EtaResolution_1pr_F_withEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F_withEGs_posEta, "TkEG_EtaResolution_3pr_F_withEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F_noEGs_posEta, "TkEG_EtaResolution_noNeutrals_F_noEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F_noEGs_posEta, "TkEG_EtaResolution_withNeutrals_F_noEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F_noEGs_posEta, "TkEG_EtaResolution_1pr_F_noEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F_noEGs_posEta, "TkEG_EtaResolution_3pr_F_noEGs_posEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // eta < 0
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F_withEGs_negEta, "TkEG_EtaResolution_noNeutrals_F_withEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F_withEGs_negEta, "TkEG_EtaResolution_withNeutrals_F_withEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F_withEGs_negEta, "TkEG_EtaResolution_1pr_F_withEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F_withEGs_negEta, "TkEG_EtaResolution_3pr_F_withEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_noNeutrals_F_noEGs_negEta, "TkEG_EtaResolution_noNeutrals_F_noEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_withNeutrals_F_noEGs_negEta, "TkEG_EtaResolution_withNeutrals_F_noEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_1pr_F_noEGs_negEta, "TkEG_EtaResolution_1pr_F_noEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_EtaResolution_3pr_F_noEGs_negEta, "TkEG_EtaResolution_3pr_F_noEGs_negEta", ";#eta resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  
  // Phi resolution of TkEG candidate
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution, "TkEG_PhiResolution", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_C, "TkEG_PhiResolution_C", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_I, "TkEG_PhiResolution_I", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F, "TkEG_PhiResolution_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals, "TkEG_PhiResolution_noNeutrals", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals, "TkEG_PhiResolution_withNeutrals", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr, "TkEG_PhiResolution_1pr", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr, "TkEG_PhiResolution_3pr", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F_withEGs, "TkEG_PhiResolution_F_withEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F_withEGs_posEta, "TkEG_PhiResolution_F_withEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F_withEGs_negEta, "TkEG_PhiResolution_F_withEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F_noEGs, "TkEG_PhiResolution_F_noEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F_noEGs_posEta, "TkEG_PhiResolution_F_noEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_F_noEGs_negEta, "TkEG_PhiResolution_F_noEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F, "TkEG_PhiResolution_noNeutrals_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F, "TkEG_PhiResolution_withNeutrals_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F, "TkEG_PhiResolution_1pr_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F, "TkEG_PhiResolution_3pr_F", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F_withEGs, "TkEG_PhiResolution_noNeutrals_F_withEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F_withEGs, "TkEG_PhiResolution_withNeutrals_F_withEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F_withEGs, "TkEG_PhiResolution_1pr_F_withEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F_withEGs, "TkEG_PhiResolution_3pr_F_withEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F_noEGs, "TkEG_PhiResolution_noNeutrals_F_noEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F_noEGs, "TkEG_PhiResolution_withNeutrals_F_noEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F_noEGs, "TkEG_PhiResolution_1pr_F_noEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F_noEGs, "TkEG_PhiResolution_3pr_F_noEGs", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // eta > 0
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F_withEGs_posEta, "TkEG_PhiResolution_noNeutrals_F_withEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F_withEGs_posEta, "TkEG_PhiResolution_withNeutrals_F_withEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F_withEGs_posEta, "TkEG_PhiResolution_1pr_F_withEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F_withEGs_posEta, "TkEG_PhiResolution_3pr_F_withEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F_noEGs_posEta, "TkEG_PhiResolution_noNeutrals_F_noEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F_noEGs_posEta, "TkEG_PhiResolution_withNeutrals_F_noEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F_noEGs_posEta, "TkEG_PhiResolution_1pr_F_noEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F_noEGs_posEta, "TkEG_PhiResolution_3pr_F_noEGs_posEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // eta < 0
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F_withEGs_negEta, "TkEG_PhiResolution_noNeutrals_F_withEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F_withEGs_negEta, "TkEG_PhiResolution_withNeutrals_F_withEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F_withEGs_negEta, "TkEG_PhiResolution_1pr_F_withEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F_withEGs_negEta, "TkEG_PhiResolution_3pr_F_withEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_noNeutrals_F_noEGs_negEta, "TkEG_PhiResolution_noNeutrals_F_noEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_withNeutrals_F_noEGs_negEta, "TkEG_PhiResolution_withNeutrals_F_noEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_1pr_F_noEGs_negEta, "TkEG_PhiResolution_1pr_F_noEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(h_TkEG_PhiResolution_3pr_F_noEGs_negEta, "TkEG_PhiResolution_3pr_F_noEGs_negEta", ";#phi resolution (GeV);Clusters / bin", 2000, -1.0, +1.0);

  // Charged Resolution
  histoTools_.BookHisto_1D(h_TkEG_ChargedResolution, "TkEG_ChargedResolution", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_ChargedResolution_noNeutrals, "TkEG_ChargedResolution_noNeutrals", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_ChargedResolution_withNeutrals, "TkEG_ChargedResolution_withNeutrals", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_ChargedResolution_1pr, "TkEG_ChargedResolution_1pr", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_ChargedResolution_3pr, "TkEG_ChargedResolution_3pr", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);

  // Neutrals Resolution
  histoTools_.BookHisto_1D(h_TkEG_NeutralsResolution, "TkEG_NeutralsResolution", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_NeutralsResolution_noNeutrals, "TkEG_NeutralsResolution_noNeutrals", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_NeutralsResolution_withNeutrals, "TkEG_NeutralsResolution_withNeutrals", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_NeutralsResolution_1pr, "TkEG_NeutralsResolution_1pr", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  histoTools_.BookHisto_1D(h_TkEG_NeutralsResolution_3pr, "TkEG_NeutralsResolution_3pr", ";#pi^{0} resolution (GeV);Clusters / bin", 21, -10.5, 10.5);
  
  // Decay Mode of non MCmatched Candidates 
  histoTools_.BookHisto_1D(h_nonMCmatchedCandidates_decayMode, "nonMCmatchedCandidates_decayMode", ";;Events", 3, -0.5 , 2.5);
  h_nonMCmatchedCandidates_decayMode->GetXaxis()->SetBinLabel(1, "e");
  h_nonMCmatchedCandidates_decayMode->GetXaxis()->SetBinLabel(2, "#mu");
  h_nonMCmatchedCandidates_decayMode->GetXaxis()->SetBinLabel(3, "other");

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

  // Turn-ons for the best performing WP
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn25, "L1TkEGTaus_TurnOn25", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn25_1pr, "L1TkEGTaus_TurnOn25_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn25_3pr, "L1TkEGTaus_TurnOn25_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn25_withNeutrals, "L1TkEGTaus_TurnOn25_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn25_noNeutrals, "L1TkEGTaus_TurnOn25_noNeutrals" , "", 60 , minEt , maxEt );

  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn50, "L1TkEGTaus_TurnOn50", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn50_1pr, "L1TkEGTaus_TurnOn50_1pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn50_3pr, "L1TkEGTaus_TurnOn50_3pr" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn50_withNeutrals, "L1TkEGTaus_TurnOn50_withNeutrals" , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hL1TkEGTaus_TurnOn50_noNeutrals, "L1TkEGTaus_TurnOn50_noNeutrals" , "", 60 , minEt , maxEt );


  // Single tau rates
  histoTools_.BookHisto_1D(hTkEG_Rate      , "TkEG_Rate"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_Rate_C    , "TkEG_Rate_C"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_Rate_I    , "TkEG_Rate_I"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_Rate_F    , "TkEG_Rate_F"    , "", nEt , minEt , maxEt );
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

  histoTools_.BookHisto_1D(hL1TkEGTaus_SingleTau_Rate, "L1TkEGTaus_SingleTau_Rate", "", nEt , minEt , maxEt );

  // Di-tau rates
    histoTools_.BookHisto_1D(hDiTau_Rate_TkEG      , "DiTau_Rate_TkEG"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_TkEG_C    , "DiTau_Rate_TkEG_C"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_TkEG_I    , "DiTau_Rate_TkEG_I"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Rate_TkEG_F    , "DiTau_Rate_TkEG_F"    , "", nEt , minEt , maxEt );
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

  histoTools_.BookHisto_1D(hL1TkEGTaus_DiTau_Rate, "L1TkEGTaus_DiTau_Rate", "", nEt , minEt , maxEt );

  // Single-tau efficiencies
  histoTools_.BookHisto_1D(hTkEG_Eff       , "TkEG_Eff"       , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_Eff_C     , "TkEG_Eff_C"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_Eff_I     , "TkEG_Eff_I"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hTkEG_Eff_F     , "TkEG_Eff_F"     , "", nEt , minEt , maxEt );
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

  histoTools_.BookHisto_1D(hL1TkEGTaus_SingleTau_Eff, "L1TkEGTaus_SingleTau_Eff", "", nEt , minEt , maxEt );

  // Di-tau efficiencies
  histoTools_.BookHisto_1D(hDiTau_Eff_TkEG       , "DiTau_Eff_TkEG"      , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_TkEG_C     , "DiTau_Eff_TkEG_C"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_TkEG_I     , "DiTau_Eff_TkEG_I"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hDiTau_Eff_TkEG_F     , "DiTau_Eff_TkEG_F"    , "", nEt , minEt , maxEt );
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

  histoTools_.BookHisto_1D(hL1TkEGTaus_DiTau_Eff, "L1TkEGTaus_DiTau_Eff", "", nEt , minEt , maxEt );

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
  h_Photons_dEta->Write();
  h_Photons_dPhi->Write();
  h_Photons_dEtaVsdPhi->Write();
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

  h_SigCone_DeltaR->Write();

  h_trkClusters_MultiplicityPerCluster->Write();
  h_MCmatch_chargedDaugh_N->Write();
  h_MCmatch_neutralDaugh_N->Write();
  h_trkClusters_Pt->Write();
  h_trkClusters_M->Write();
  h_trkClusters_M_beforeCut->Write();

  h_clustEGs_Et->Write();
  h_clustEGs_Eta->Write();
  h_clustEGs_Phi->Write();
  h_clustEGs_M->Write();

  h_EGs_N -> Write();
  h_EGs_MCmatched_Et->Write();
  h_EGs_Et->Write();
  h_EGs_Eta->Write();
  h_EGs_Phi->Write();
//  h_EGs_IEta->Write();
//  h_EGs_IPhi->Write();
  h_EGs_HwQual->Write();
  h_EGs_EtaVsEt->Write();

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

  h_TkEGClusters_M_beforeCut->Write();
  h_EGClusters_MultiplicityPerCluster->Write();
  h_EGClusters_Et->Write();
  h_EGClusters_M->Write();

  h_TkEG_relIso->Write();
  h_TkEG_vtxIso->Write();

  h_TkEG_N->Write();
  h_TkEG_RelIso_N->Write();
  h_TkEG_VtxIso_N->Write();
  h_TkEG_RelIsoLoose_N->Write();
  h_TkEG_VtxIsoLoose_N->Write();
  h_TkEG_RelIsoTight_N->Write();
  h_TkEG_VtxIsoTight_N->Write();

  h_TkEG_Pt->Write();
  h_TkEG_ET->Write();
  h_TkEG_Eta->Write();
  h_TkEG_Phi->Write();
  h_TkEG_InvMass->Write();
  h_TkEG_NEGs->Write();
  h_TkEG_NEGs_withEGs->Write();
  h_TkEG_NTks->Write();
  h_TkEG_NEGs_C->Write();
  h_TkEG_NEGs_I->Write();
  h_TkEG_NEGs_F->Write();
  //h_TkEG_CHF->Write();
  //h_TkEG_NHF->Write();
  h_TkEG_CHF_withNeutrals->Write();
  h_TkEG_NHF_withNeutrals->Write();
  h_TkEG_clustEGs_MCMatch->Write();
  h_TkEG_clustEGs_Matched_HwQual->Write();
  h_TkEG_clustEGs_nonMatched_HwQual->Write();
  h_TkEG_clustEGs_dET_matchPion0->Write();
  h_TkEG_clustEGs_ETResolution->Write();
  h_TkEG_isoTracks_InvMass->Write();
  h_TkEG_isoTracks_Multiplicity->Write();
  h_TkEG_isoTracks_Et->Write();
  h_TkEG_isoTracks_Eta->Write();
  h_TkEG_DonutRatio->Write();
  h_TkEG_signalEGs_Multiplicity->Write();
  h_TkEG_isoEGs_Multiplicity->Write();

  // Poor Et Resolution Candidates Properties
  h_PoorEtResolCand_InvMass->Write();
  h_PoorEtResolCand_RelIso->Write();
  h_PoorEtResolCand_VtxIso->Write();
  h_PoorEtResolCand_CHF->Write();
  h_PoorEtResolCand_IsoTracks_N->Write();
  h_PoorEtResolCand_dR_EG_Seed->Write();
  // Good Et resolution candidates
  h_GoodEtResolCand_InvMass->Write();
  h_GoodEtResolCand_RelIso->Write();
  h_GoodEtResolCand_VtxIso->Write();
  h_GoodEtResolCand_CHF->Write();
  h_GoodEtResolCand_IsoTracks_N->Write();
  h_GoodEtResolCand_dR_EG_Seed->Write();

  // Poor Neutral Resolution candidates
  h_TkEG_PoorNeuResol_NeuMultiplicity->Write();
  h_TkEG_PoorNeuResol_dR_Pi0_visTau->Write();
  h_TkEG_PoorNeuResol_dEta_Pi0_visTau->Write();
  h_TkEG_PoorNeuResol_dPhi_Pi0_visTau->Write();
  h_TkEG_PoorNeuResol_Pi0_ET->Write();
  h_TkEG_PoorNeuResol_dRmin_Pi0_EG->Write();
  h_TkEG_PoorNeuResol_dRmin_Seed_closestEG->Write();
  h_TkEG_PoorNeuResol_Pi0_closestEG_ET->Write();
  h_TkEG_PoorNeuResol_Pi0_closestEG_ET_Vs_dRmin_Pi0_EG->Write();
  h_TkEG_PoorNeuResol_Pi0_ET_Vs_closestEG_ET->Write();

  h_TkEG_PtResolution->Write();
  h_TkEG_PtResolution_C->Write();
  h_TkEG_PtResolution_I->Write();
  h_TkEG_PtResolution_F->Write();
  h_TkEG_PtResolution_noNeutrals->Write();
  h_TkEG_PtResolution_withNeutrals->Write();
  h_TkEG_PtResolution_1pr->Write();
  h_TkEG_PtResolution_3pr->Write();

  h_TkEG_PtResolution_F_withEGs->Write();
  h_TkEG_PtResolution_F_withEGs_posEta->Write();
  h_TkEG_PtResolution_F_withEGs_negEta->Write();
  h_TkEG_PtResolution_F_noEGs->Write();
  h_TkEG_PtResolution_F_noEGs_posEta->Write();
  h_TkEG_PtResolution_F_noEGs_negEta->Write();
  
  h_TkEG_PtResolution_withNeutrals_F->Write();
  h_TkEG_PtResolution_1pr_F->Write();
  h_TkEG_PtResolution_3pr_F->Write();

  h_TkEG_PtResolution_noNeutrals_F_withEGs->Write();
  h_TkEG_PtResolution_withNeutrals_F_withEGs->Write();
  h_TkEG_PtResolution_1pr_F_withEGs->Write();
  h_TkEG_PtResolution_3pr_F_withEGs->Write();
  h_TkEG_PtResolution_noNeutrals_F_noEGs->Write();
  h_TkEG_PtResolution_withNeutrals_F_noEGs->Write();
  h_TkEG_PtResolution_1pr_F_noEGs->Write();
  h_TkEG_PtResolution_3pr_F_noEGs->Write();

  // eta > 0
  h_TkEG_PtResolution_noNeutrals_F_withEGs_posEta->Write();
  h_TkEG_PtResolution_withNeutrals_F_withEGs_posEta->Write();
  h_TkEG_PtResolution_1pr_F_withEGs_posEta->Write();
  h_TkEG_PtResolution_3pr_F_withEGs_posEta->Write();
  h_TkEG_PtResolution_noNeutrals_F_noEGs_posEta->Write();
  h_TkEG_PtResolution_withNeutrals_F_noEGs_posEta->Write();
  h_TkEG_PtResolution_1pr_F_noEGs_posEta->Write();
  h_TkEG_PtResolution_3pr_F_noEGs_posEta->Write();

  // eta < 0
  h_TkEG_PtResolution_noNeutrals_F_withEGs_negEta->Write();
  h_TkEG_PtResolution_withNeutrals_F_withEGs_negEta->Write();
  h_TkEG_PtResolution_1pr_F_withEGs_negEta->Write();
  h_TkEG_PtResolution_3pr_F_withEGs_negEta->Write();
  h_TkEG_PtResolution_noNeutrals_F_noEGs_negEta->Write();
  h_TkEG_PtResolution_withNeutrals_F_noEGs_negEta->Write();
  h_TkEG_PtResolution_1pr_F_noEGs_negEta->Write();
  h_TkEG_PtResolution_3pr_F_noEGs_negEta->Write();

  h_TkEG_EtResolution->Write();
  h_TkEG_EtResolution_C->Write();
  h_TkEG_EtResolution_I->Write();
  h_TkEG_EtResolution_F->Write();
  h_TkEG_EtResolution_noNeutrals->Write();
  h_TkEG_EtResolution_withNeutrals->Write();
  h_TkEG_EtResolution_1pr->Write();
  h_TkEG_EtResolution_3pr->Write();
  h_TkEG_EtResolution_noEGs->Write();
  h_TkEG_EtResolution_withEGs->Write();

  h_TkEG_EtResolution_withNeutrals_1pr->Write();
  h_TkEG_EtResolution_withNeutrals_3pr->Write();
  h_TkEG_EtResolution_withNeutrals_1pion0->Write();
  h_TkEG_EtResolution_withNeutrals_2pion0->Write();
  h_TkEG_EtResolution_withNeutrals_3pion0->Write();
  h_TkEG_EtResolution_withNeutrals_4pion0->Write();
    
  h_TkEG_EtResolution_noNeutrals_withEGs->Write();
  h_TkEG_EtResolution_withNeutrals_withEGs->Write();
  h_TkEG_EtResolution_1pr_withEGs->Write();
  h_TkEG_EtResolution_3pr_withEGs->Write();
  h_TkEG_EtResolution_noNeutrals_noEGs->Write();
  h_TkEG_EtResolution_withNeutrals_noEGs->Write();
  h_TkEG_EtResolution_1pr_noEGs->Write();
  h_TkEG_EtResolution_3pr_noEGs->Write();

  h_TkEG_EtResolution_withNeutrals_withEGs_0to5GeV->Write();
  h_TkEG_EtResolution_withNeutrals_withEGs_5to10GeV->Write();
  h_TkEG_EtResolution_withNeutrals_withEGs_10to15GeV->Write();
  h_TkEG_EtResolution_withNeutrals_withEGs_15to20GeV->Write();
  h_TkEG_EtResolution_withNeutrals_withEGs_20to30GeV->Write();
  h_TkEG_EtResolution_withNeutrals_withEGs_30to40GeV->Write();
  h_TkEG_EtResolution_withNeutrals_withEGs_40to50GeV->Write();

  h_TkEG_EtResolution_F_withEGs->Write();
  h_TkEG_EtResolution_F_withEGs_posEta->Write();
  h_TkEG_EtResolution_F_withEGs_negEta->Write();
  h_TkEG_EtResolution_F_noEGs->Write();
  h_TkEG_EtResolution_F_noEGs_posEta->Write();
  h_TkEG_EtResolution_F_noEGs_negEta->Write();

  h_TkEG_EtResolution_noNeutrals_F->Write();
  h_TkEG_EtResolution_withNeutrals_F->Write();
  h_TkEG_EtResolution_1pr_F->Write();
  h_TkEG_EtResolution_3pr_F->Write();

  h_TkEG_EtResolution_noNeutrals_F_withEGs->Write();
  h_TkEG_EtResolution_withNeutrals_F_withEGs->Write();
  h_TkEG_EtResolution_1pr_F_withEGs->Write();
  h_TkEG_EtResolution_3pr_F_withEGs->Write();
  h_TkEG_EtResolution_noNeutrals_F_noEGs->Write();
  h_TkEG_EtResolution_withNeutrals_F_noEGs->Write();
  h_TkEG_EtResolution_1pr_F_noEGs->Write();
  h_TkEG_EtResolution_3pr_F_noEGs->Write();

  // eta > 0
  h_TkEG_EtResolution_noNeutrals_F_withEGs_posEta->Write();
  h_TkEG_EtResolution_withNeutrals_F_withEGs_posEta->Write();
  h_TkEG_EtResolution_1pr_F_withEGs_posEta->Write();
  h_TkEG_EtResolution_3pr_F_withEGs_posEta->Write();
  h_TkEG_EtResolution_noNeutrals_F_noEGs_posEta->Write();
  h_TkEG_EtResolution_withNeutrals_F_noEGs_posEta->Write();
  h_TkEG_EtResolution_1pr_F_noEGs_posEta->Write();
  h_TkEG_EtResolution_3pr_F_noEGs_posEta->Write();

  // eta < 0
  h_TkEG_EtResolution_noNeutrals_F_withEGs_negEta->Write();
  h_TkEG_EtResolution_withNeutrals_F_withEGs_negEta->Write();
  h_TkEG_EtResolution_1pr_F_withEGs_negEta->Write();
  h_TkEG_EtResolution_3pr_F_withEGs_negEta->Write();
  h_TkEG_EtResolution_noNeutrals_F_noEGs_negEta->Write();
  h_TkEG_EtResolution_withNeutrals_F_noEGs_negEta->Write();
  h_TkEG_EtResolution_1pr_F_noEGs_negEta->Write();
  h_TkEG_EtResolution_3pr_F_noEGs_negEta->Write();

  h_TkEG_EtaResolution->Write();
  h_TkEG_EtaResolution_C->Write();
  h_TkEG_EtaResolution_I->Write();
  h_TkEG_EtaResolution_F->Write();
  h_TkEG_EtaResolution_noNeutrals->Write();
  h_TkEG_EtaResolution_withNeutrals->Write();
  h_TkEG_EtaResolution_1pr->Write();
  h_TkEG_EtaResolution_3pr->Write();

  h_TkEG_EtaResolution_F_withEGs->Write();
  h_TkEG_EtaResolution_F_withEGs_posEta->Write();
  h_TkEG_EtaResolution_F_withEGs_negEta->Write();
  h_TkEG_EtaResolution_F_noEGs->Write();
  h_TkEG_EtaResolution_F_noEGs_posEta->Write();
  h_TkEG_EtaResolution_F_noEGs_negEta->Write();

  h_TkEG_EtaResolution_noNeutrals_F->Write();
  h_TkEG_EtaResolution_withNeutrals_F->Write();
  h_TkEG_EtaResolution_1pr_F->Write();
  h_TkEG_EtaResolution_3pr_F->Write();

  h_TkEG_EtaResolution_noNeutrals_F_withEGs->Write();
  h_TkEG_EtaResolution_withNeutrals_F_withEGs->Write();
  h_TkEG_EtaResolution_1pr_F_withEGs->Write();
  h_TkEG_EtaResolution_3pr_F_withEGs->Write();
  h_TkEG_EtaResolution_noNeutrals_F_noEGs->Write();
  h_TkEG_EtaResolution_withNeutrals_F_noEGs->Write();
  h_TkEG_EtaResolution_1pr_F_noEGs->Write();
  h_TkEG_EtaResolution_3pr_F_noEGs->Write();

  // eta > 0
  h_TkEG_EtaResolution_noNeutrals_F_withEGs_posEta->Write();
  h_TkEG_EtaResolution_withNeutrals_F_withEGs_posEta->Write();
  h_TkEG_EtaResolution_1pr_F_withEGs_posEta->Write();
  h_TkEG_EtaResolution_3pr_F_withEGs_posEta->Write();
  h_TkEG_EtaResolution_noNeutrals_F_noEGs_posEta->Write();
  h_TkEG_EtaResolution_withNeutrals_F_noEGs_posEta->Write();
  h_TkEG_EtaResolution_1pr_F_noEGs_posEta->Write();
  h_TkEG_EtaResolution_3pr_F_noEGs_posEta->Write();

  // eta < 0
  h_TkEG_EtaResolution_noNeutrals_F_withEGs_negEta->Write();
  h_TkEG_EtaResolution_withNeutrals_F_withEGs_negEta->Write();
  h_TkEG_EtaResolution_1pr_F_withEGs_negEta->Write();
  h_TkEG_EtaResolution_3pr_F_withEGs_negEta->Write();
  h_TkEG_EtaResolution_noNeutrals_F_noEGs_negEta->Write();
  h_TkEG_EtaResolution_withNeutrals_F_noEGs_negEta->Write();
  h_TkEG_EtaResolution_1pr_F_noEGs_negEta->Write();
  h_TkEG_EtaResolution_3pr_F_noEGs_negEta->Write();

  h_TkEG_PhiResolution->Write();
  h_TkEG_PhiResolution_C->Write();
  h_TkEG_PhiResolution_I->Write();
  h_TkEG_PhiResolution_F->Write();
  h_TkEG_PhiResolution_noNeutrals->Write();
  h_TkEG_PhiResolution_withNeutrals->Write();
  h_TkEG_PhiResolution_1pr->Write();
  h_TkEG_PhiResolution_3pr->Write();

  h_TkEG_PhiResolution_F_withEGs->Write();
  h_TkEG_PhiResolution_F_withEGs_posEta->Write();
  h_TkEG_PhiResolution_F_withEGs_negEta->Write();
  h_TkEG_PhiResolution_F_noEGs->Write();
  h_TkEG_PhiResolution_F_noEGs_posEta->Write();
  h_TkEG_PhiResolution_F_noEGs_negEta->Write();

  h_TkEG_PhiResolution_noNeutrals_F->Write();
  h_TkEG_PhiResolution_withNeutrals_F->Write();
  h_TkEG_PhiResolution_1pr_F->Write();
  h_TkEG_PhiResolution_3pr_F->Write();

  h_TkEG_PhiResolution_noNeutrals_F_withEGs->Write();
  h_TkEG_PhiResolution_withNeutrals_F_withEGs->Write();
  h_TkEG_PhiResolution_1pr_F_withEGs->Write();
  h_TkEG_PhiResolution_3pr_F_withEGs->Write();
  h_TkEG_PhiResolution_noNeutrals_F_noEGs->Write();
  h_TkEG_PhiResolution_withNeutrals_F_noEGs->Write();
  h_TkEG_PhiResolution_1pr_F_noEGs->Write();
  h_TkEG_PhiResolution_3pr_F_noEGs->Write();

  // eta > 0
  h_TkEG_PhiResolution_noNeutrals_F_withEGs_posEta->Write();
  h_TkEG_PhiResolution_withNeutrals_F_withEGs_posEta->Write();
  h_TkEG_PhiResolution_1pr_F_withEGs_posEta->Write();
  h_TkEG_PhiResolution_3pr_F_withEGs_posEta->Write();
  h_TkEG_PhiResolution_noNeutrals_F_noEGs_posEta->Write();
  h_TkEG_PhiResolution_withNeutrals_F_noEGs_posEta->Write();
  h_TkEG_PhiResolution_1pr_F_noEGs_posEta->Write();
  h_TkEG_PhiResolution_3pr_F_noEGs_posEta->Write();

  // eta < 0
  h_TkEG_PhiResolution_noNeutrals_F_withEGs_negEta->Write();
  h_TkEG_PhiResolution_withNeutrals_F_withEGs_negEta->Write();
  h_TkEG_PhiResolution_1pr_F_withEGs_negEta->Write();
  h_TkEG_PhiResolution_3pr_F_withEGs_negEta->Write();
  h_TkEG_PhiResolution_noNeutrals_F_noEGs_negEta->Write();
  h_TkEG_PhiResolution_withNeutrals_F_noEGs_negEta->Write();
  h_TkEG_PhiResolution_1pr_F_noEGs_negEta->Write();
  h_TkEG_PhiResolution_3pr_F_noEGs_negEta->Write();

  // Charged Resolution
  h_TkEG_ChargedResolution->Write();
  h_TkEG_ChargedResolution_noNeutrals->Write();
  h_TkEG_ChargedResolution_withNeutrals->Write();
  h_TkEG_ChargedResolution_1pr->Write();
  h_TkEG_ChargedResolution_3pr->Write();

  // Neutrals Resolution
  h_TkEG_NeutralsResolution->Write();
  h_TkEG_NeutralsResolution_noNeutrals->Write();
  h_TkEG_NeutralsResolution_withNeutrals->Write();
  h_TkEG_NeutralsResolution_1pr->Write();
  h_TkEG_NeutralsResolution_3pr->Write();

  h_ldgTkEG_ET->Write();
  hTkEG_DeltaRmatch->Write();
  hTkEG_genVisEt->Write();
  hTkEG_genVisEt_clustEG->Write();
  hTkEG_genVisPt_clustEG->Write();

  
  h_MCmatch_dR->Write(); 
  h_leadTrk_MCmatch->Write(); 
  h_leadTrk4stubs_MCmatch->Write();
  
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

  // Turn-ons for the best performing WP
  hL1TkEGTaus_TurnOn25->Write();
  hL1TkEGTaus_TurnOn25_1pr->Write();
  hL1TkEGTaus_TurnOn25_3pr->Write();
  hL1TkEGTaus_TurnOn25_withNeutrals->Write();
  hL1TkEGTaus_TurnOn25_noNeutrals->Write();

  hL1TkEGTaus_TurnOn50->Write();
  hL1TkEGTaus_TurnOn50_1pr->Write();
  hL1TkEGTaus_TurnOn50_3pr->Write();
  hL1TkEGTaus_TurnOn50_withNeutrals->Write();
  hL1TkEGTaus_TurnOn50_noNeutrals->Write();

  // SingleTau: Efficiencies
  hTkEG_Rate->Write();
  hTkEG_Rate_C->Write();
  hTkEG_Rate_I->Write();
  hTkEG_Rate_F->Write();
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
  
  hL1TkEGTaus_SingleTau_Rate->Write();

  // DiTau: Rates
  hDiTau_Rate_TkEG->Write();
  hDiTau_Rate_TkEG_C->Write();
  hDiTau_Rate_TkEG_I->Write();
  hDiTau_Rate_TkEG_F->Write();
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

  hL1TkEGTaus_DiTau_Rate->Write();

  // SingleTau: Efficiencies
  hTkEG_Eff->Write();
  hTkEG_Eff_C->Write();
  hTkEG_Eff_I->Write();
  hTkEG_Eff_F->Write();
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

  hL1TkEGTaus_SingleTau_Eff->Write();

  // DiTau: Efficiencies
  hDiTau_Eff_TkEG->Write();
  hDiTau_Eff_TkEG_C->Write();
  hDiTau_Eff_TkEG_I->Write();
  hDiTau_Eff_TkEG_F->Write();
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

  hL1TkEGTaus_DiTau_Eff->Write();

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
void TkEG::ApplyDiTauZMatching(vector<L1TkEGParticle> &L1TkEGs)
//============================================================================
{
  
  //
  // NOTE:
  // Is this function needed? It was introduced to remove from the collection
  // L1 Tau candidates that have a Z0 value which is more than X cm away from 
  // the leading-in=ET L1 Tau candidate. The objects are permanently removedf
  // rom the collection input, which is dangerous!
  //

  // Sanity check
  if (L1TkEGs.size() < 2) return;

  // Initialise variables
  double deltaPOCAz = 9999.9;

  TTTrack match_tk0 = L1TkEGs.at(0).GetLeadingTrack(); 

  // For-loop: L1TkEGs
  for (size_t i = 1; i < L1TkEGs.size(); i++)
    {

      TTTrack match_tk = L1TkEGs.at(i).GetLeadingTrack();
      deltaPOCAz = abs( match_tk0.getZ0() - match_tk.getZ0() );
    
    // If the Trigger objects is not within x-cm reject it
    if (deltaPOCAz > diTau_deltaPOCAz) L1TkEGs.erase ( L1TkEGs.begin()+i );
   
    }  // For-loop: L1TkEGs
  
  return;
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
  vector<L1TkEGParticle> L1TkEGs_mcMatched = GetMcMatchedL1TkEGs(L1TkEGs);
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
  /*
  TH1D *h = new TH1D("TotalEvents", "nEvtsTotal", 300, 0.0, 300);
  for (int i = 1; i<= h->GetNbinsX(); i++) {
    h->SetBinContent(i, nEvtsTotal);
  }

  // for (int i = 1; i<= h->GetNbinsX(); i++) { 
  //   cout <<  h->GetBinContent(i)<<endl;
  // }

  TEfficiency* hEff = 0;
  hEff = new TEfficiency(*histo, *h);
  `*/
  const int nBins = histo->GetNbinsX()+1;
  double eff, err;

  // For-loop: Histogram bins
  for (int i = 1; i<= nBins; i++){
    
    const int nPass = histo->GetBinContent(i);
    auxTools_.Efficiency(nPass, nEvtsTotal, "binomial", eff, err ); //fixme: use TEfficiency?

    // Update current histo bin to true eff value and error
    histo->SetBinContent(i, eff);
    histo->SetBinError  (i, err);

    //cout<<"previous = "<< histo->GetBinContent(i)<<" +- "<<err<<endl;
    //cout<<"now      = "<< hEff->GetEfficiency(i)<<" +- "<<hEff->GetEfficiencyErrorUp(i)<<endl;
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

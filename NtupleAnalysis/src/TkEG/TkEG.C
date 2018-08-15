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
  cfg_tk_maxEta      = 1.44;        // Default: 1e6
  cfg_tk_maxChiSqRed = 1e6;         // Default: 1e6
  cfg_tk_minStubs    =   4;         // Default: 0
  cfg_tk_minStubsPS  =   0;         // Default: 0
  cfg_tk_maxStubsPS  = 1e6;         // Default: 1e6

  // TkEG algorithm parameters
  ETmin = 5.0; // min ET in GeV of L1EG objects
  ZMAX = 10.0; // |z_track| < ZMAX in cm
  CHI2MAX = 200;		
  DRmin=0.1;
  DRmax=0.3;
  minPt_leading_track = 10.0;
  PrimaryVtxConstrain = false; // use the primary vertex
  DeltaZMax = 1.0; // | z_track - z_primaryvtx | < DeltaZMax in cm. 
  IsoCut=0.15; 
  RelativeIsolation = true;

  // New parameters
  minStubs_trk      = 4; 
  maxChi2_trk       = 40.0; // GeV
  maxChi2_trk_alt   = 20.0; // GeV
  minPt_leadtrk     = 10.0; // GeV
  maxEta_leadtrk    = 1.44;
  minDeltaR_leadtrk = 0.0;
  maxDeltaR_leadtrk = 0.3;
  maxDeltaZ_trk     = 0.8;  // cm
  maxInvMass_trk    = 1.40; // GeV 
  minEt_EG          = 2.0;  // GeV
  minDeltaR_EG      = 0.0;
  maxDeltaR_EG      = 0.25;
  maxInvMass_EG     = 1.77; // GeV
  maxDeltaR_MCmatch = 0.2;    
  minDeltaR_iso     = 0.0;
  maxDeltaR_iso     = 5.0;
  maxDeltaZ_iso     = 0.8;  // cm
  useRelIso         = true;
  maxRelIso         = 0.20; //0.15
  
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
  
  cout << "=== TkEG:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;
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
  unsigned int counter_allTrkTauCand = 0;
  unsigned int counter_hasGenTau = 0;
  unsigned int counter_hasEG = 0;
  
  for (int jentry = 0; jentry < nEntries; jentry++, nEvts++){
    
    if(cfg_DEBUG) cout << "\t------------ Entry = " << jentry << endl;
    
    // Init variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    
    ////////////////////////////////////////////////
    // Build collections
    ////////////////////////////////////////////////       

    // Tracks
    if (cfg_AddL1Tks) {
      TTTracks = GetTTTracks(cfg_tk_minPt, cfg_tk_minEta, cfg_tk_maxEta, cfg_tk_maxChiSqRed, cfg_tk_minStubs, cfg_tk_nFitParams, false);
      sort( TTTracks.begin(), TTTracks.end() ); // Sort from highest Pt to lowest (not done by default)
      if (cfg_DEBUG*0) cout << "\n=== TTracks (" << L1Tks_Pt->size() << ")" << endl;
      if (cfg_DEBUG*0) PrintTTTrackCollection(TTTracks);
    }    
    
    // EGs
    if (cfg_AddEGs) {
      L1EGs = GetL1EGs(false);
      sort( L1EGs.begin(), L1EGs.end() ); // Sort from highest Et to lowest Et (should be already done by default)
      if (cfg_DEBUG*0) cout << "\n=== L1EGs (" << L1EGs.size() << ")" << endl;
      if (cfg_DEBUG*0) PrintL1EGCollection(L1EGs);
    }
    
    // GenParticles (skip for MinBias samples as no real taus exist)
    vector<GenParticle> GenTaus;
    if (!isMinBias) GenTaus = GetGenParticles(15, true);
    if (cfg_DEBUG*0) PrintGenParticleCollection(GenTaus);
    
    // Hadronic GenTaus (skip for MinBias samples)
    GenTausHadronic.clear();
    if (cfg_AddGenP) {
      if (cfg_DEBUG*0) cout << "\n=== GenParticles (" << GenP_Pt.size() << ")" << endl;
      if (!isMinBias) GenTausHadronic = GetHadronicGenTaus(GenTaus, 00.0, 1.4);
      if (cfg_DEBUG*0) PrintGenParticleCollection(GenTausHadronic);
    }
    
    // Triggred GenTaus (skip for MinBias samples)
    vector<GenParticle> GenTausTrigger;  
    if (!isMinBias) GenTausTrigger = GetHadronicGenTaus(GenTaus, 20.0, 1.3);  
    
    // Ensure that all taus are found, needed by the current efficiency definition 
    // E.g. for ttbar, only events with two taus within trigger acceptance are considered for efficiency calculation)
    bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;

    // Fill histos with genTau info
    if (!isMinBias){
      h_genTausAll_N->Fill(GenTaus.size());
      for (vector<GenParticle>::iterator genTau = GenTaus.begin(); genTau != GenTaus.end(); genTau++) {
	h_genTausAll_Pt  -> Fill(genTau->pt());
	h_genTausAll_Eta -> Fill(genTau->eta());
	h_genTausAll_Phi -> Fill(genTau->phi());
      }
      h_genTausAll_Eta1VsEta2 -> Fill (GenTaus.at(0).eta(), GenTaus.at(1).eta());
      h_genTausAll_Phi1VsPhi2 -> Fill (GenTaus.at(0).phi(), GenTaus.at(1).phi());
    }

    ////////////////////////////////////////////////
    // TkEG algorithm
    ////////////////////////////////////////////////    
    
    // Consider only events with at least one genuine hadronic tau (except for MinBias sample)
    h_genTausHad_N->Fill( GenTausHadronic.size() );
    if (GenTausHadronic.size() < 1 && !isMinBias) continue;


    //================ Gen-Level ===================

    // For-loop: All hadronic gen-taus

    for (vector<GenParticle>::iterator genTau = GenTausHadronic.begin(); genTau != GenTausHadronic.end(); genTau++) {
 
      // Fill histo with the number of the decay products of hadronic taus
      h_genTausHad_Daughters_N-> Fill(genTau->finalDaughters().size());
      h_genTausHad_chargedDaugh_N-> Fill(genTau->finalDaughtersCharged().size());
      h_genTausHad_neutralDaugh_N-> Fill(genTau->finalDaughtersNeutral().size());

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
	  if ( d_2.pt() > 2.0 ) {
	    // Find the maximum deltaR of the leading charged product and the other charged products
	    deltaRcharged = auxTools_.DeltaR(d.eta(), d.phi(), d_2.eta(), d_2.phi());
	    if (deltaRcharged > deltaRcharged_max) deltaRcharged_max = deltaRcharged;
	  }

	}
	if (notLeading) continue;
	
	// Sum of ET of neutral products should be greater than 2 GeV
	sumET =0;
	for (unsigned int i = 0; i < genTau->finalDaughtersNeutral().size(); i++) sumET += genTau->finalDaughtersNeutral().at(i).et();
	if (sumET > 2.0) {
	  // For-loop: All neutral daughters
	  for (unsigned int i = 0; i < genTau->finalDaughtersNeutral().size(); i++) {
	    GenParticle d_n      = genTau->finalDaughtersNeutral().at(i);
	    
	    // Find the maximum deltaR of the leading charged product and the other neutral products
	    deltaRneutral = auxTools_.DeltaR(d.eta(), d.phi(), d_n.eta(), d_n.phi());
	    if (deltaRneutral > deltaRneutral_max) deltaRneutral_max = deltaRneutral;
	  }
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
	
      }

    }

    //============= Leading-Tracks =================
    
    // Select lead tracks  
    vector<unsigned short> leadTrackIndices;
    trackTauCandidates.clear();
    float deltaR;
    bool highPtNeighbourFound;
    vector< TTTrack> newCandidateTracks; // temporary container
    for (vector<TTTrack>::iterator leadTrkIter = TTTracks.begin(); leadTrkIter != TTTracks.end(); leadTrkIter++) {
      if (cfg_DEBUG*0) std:: cout << "Considering leading track candidate " << leadTrkIter->index() << std::endl;

      // General cut on chi^2 
      if (leadTrkIter->getChi2() > 100.0) continue;

      // Pt and eta cuts
      if (leadTrkIter->getPt() < minPt_leadtrk) continue;
      if (leadTrkIter->getEta() > maxEta_leadtrk) continue;
      
      // Plot quality quantities for the TTTracks 
      if (leadTrkIter->getNumOfStubs() == minStubs_trk) h_trk_Chi2_all_4stubs -> Fill (leadTrkIter->getChi2());
      if (leadTrkIter->getNumOfStubs() > minStubs_trk) h_trk_Chi2_all_5stubs -> Fill (leadTrkIter->getChi2());

      h_trk_NStubs_all -> Fill (leadTrkIter->getNumOfStubs());
      h_trk_Chi2_all   -> Fill (leadTrkIter->getChi2());
      h_trk_Chi2Red_all-> Fill (leadTrkIter->getChi2Red());

      // Only use high quality tracks
      if (leadTrkIter->getNumOfStubs() < minStubs_trk) continue;
      if (leadTrkIter->getChi2() > maxChi2_trk) continue;
      if ((leadTrkIter->getNumOfStubs() > minStubs_trk) && (leadTrkIter->getChi2() > maxChi2_trk_alt)) continue;
      
      if (cfg_DEBUG*0) std:: cout << "Leading track candidate " << leadTrkIter->index() << "passed cuts!" << std::endl;
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
      if (cfg_DEBUG) std::cout << "leading track index = " << trackTauCandidates[i][0].index() << ", Pt = " << trackTauCandidates[i][0].getPt() << std::endl;
    }  
    // TODO: fill counters in h_leadTrkSelection
   

    
    // Cluster surrounding tracks with lead tracks
    TTTrack* leadTrackPtr = NULL;
    float invMass;
    float pT;
    bool stopClustering;
    //TLorentzVector p4sum; //marina
    for (size_t i=0; i<trackTauCandidates.size(); i++) {
        leadTrackPtr = &(trackTauCandidates[i][0]);
        int Ntrks = 1;
        // Track clustering counters
        unsigned int trkcounter_allTracks = 0;
        unsigned int trkcounter_allNonLeading = 0;
        unsigned int trkcounter_passZ = 0;
        unsigned int trkcounter_passDRmax = 0;
        unsigned int trkcounter_passDRmin = 0;
        unsigned int trkcounter_passInvMass = 0;    


        if (cfg_DEBUG*0) cout << "Starting to cluster lead track " << leadTrackPtr->index();
        // Loop over other tracks
        stopClustering = false;
        for (vector<TTTrack>::iterator trackIter = TTTracks.begin(); ( !stopClustering && (trackIter != TTTracks.end()) ); trackIter++) {
	  // General cut on chi^2 
	  if (trackIter->getChi2() > 100.0) continue;

          trkcounter_allTracks += 1;
          // Do not double-counts the lead track
          if (trackIter->index() == leadTrackPtr->index()) continue;
          trkcounter_allNonLeading += 1;
          // Skip tracks that are not close in terms of deltaR and deltaZ
          if (abs (trackIter->getZ0() - leadTrackPtr->getZ0() ) > maxDeltaZ_trk) continue;
          trkcounter_passZ += 1;
          deltaR = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), trackIter->getEta(), trackIter->getPhi());
          if (deltaR > maxDeltaR_leadtrk) continue;
          trkcounter_passDRmax += 1;
          if (deltaR < minDeltaR_leadtrk) continue;
          trkcounter_passDRmin += 1;
          // Calculate invariant mass of already clustered tracks + new track
          invMass = 0.0;
          pT = 0.0;
	  TLorentzVector p4sum; // initialized to (0,0,0,0)  //marina
          for(size_t j=0; j < trackTauCandidates[i].size(); j++){
            p4sum += trackTauCandidates[i][j].p4();
          }
          invMass = p4sum.M();
          pT = p4sum.Pt();
          p4sum += trackIter->p4();
          if (p4sum.M() < maxInvMass_trk) {
            trackTauCandidates[i].push_back(*trackIter);
            h_clustTrks_Pt->Fill( trackIter->getPt() );
            h_clustTrks_Eta->Fill( trackIter->getEta() );
            h_clustTrks_Phi->Fill( trackIter->getPhi() );
            h_clustTrks_Phi_Eta->Fill( trackIter->getPhi(),trackIter->getEta() );
            trkcounter_passInvMass += 1;
          }  
          else {
            // Stop clustering and fill track cluster histograms
            stopClustering = true;
            Ntrks = trackTauCandidates[i].size();
            h_trkClusters_MultiplicityPerCluster->Fill( Ntrks );
            h_trkClusters_Pt->Fill( pT );
            h_trkClusters_M->Fill( invMass );
          }
        }
        if (cfg_DEBUG) cout << "After clustering trackTauCandidates[" << i << "] with leadTrack " << leadTrackPtr->index() << ", it has " 
                            << trackTauCandidates[i].size() << " tracks and a mass of " << invMass << endl;

        // FIX ME: Fill track clustering counter histogram 
        //h_clustTrks_counter->Fill(1,trkcounter_allTracks,1);
        //h_clustTrks_counter->Fill(2,trkcounter_allNonLeading,1);
        //h_clustTrks_counter->Fill(3,trkcounter_passZ,1);
        //h_clustTrks_counter->Fill(4,trkcounter_passDRmax,1);
        //h_clustTrks_counter->Fill(5,trkcounter_passDRmin,1);
        //h_clustTrks_counter->Fill(6,trkcounter_passInvMass,1);
    }
    
    // EGs properties 
    // Fill histo with number of EGs in each event
    h_EGs_N -> Fill( L1EGs.size() );
    
    // For-loop: All the EGs in the event
    for (auto eg = L1EGs.begin(); eg != L1EGs.end() ; eg++) {
      
      // Fill histos with EGs properties 
      h_EGs_Et  -> Fill( eg->getEt() );
      h_EGs_Eta -> Fill( eg->getEta() );
      h_EGs_Phi -> Fill( eg->getPhi() );
      h_EGs_IEta -> Fill( eg->getIEta() );
      h_EGs_IPhi -> Fill( eg->getIPhi() );

    }// For-loop: All the EGs in the event
    
    // Build EG clusters and create tau candidates
    vector<L1EG> EGcluster;
    TauCandidates.clear();
    double ET;
    float deltaRmin;
    float deltaR_beforecorrection;
    //deltaR = 999.0; //marina
    for (size_t i=0; i<trackTauCandidates.size(); i++) {
      //cout << "Clustering EGs around trackTauCandidates[" << i << "]" << endl;
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

      //Increase counter for all trackTauCandidates 
      counter_allTrkTauCand++;

      // Calculate the p4sum of the clustered tracks for each trackTau Candidate 
      for(size_t j=0; j < trackTauCandidates[i].size(); j++){
	p4sum += trackTauCandidates[i][j].p4();
      }
      
      // Initialize the ET (scalar sum of the ET of EGs) and the InvMass as InvMass of the clustered tracks
      ET      = 0.0; //= p4sum.Et();
      invMass = p4sum.M();
      deltaRmin = 1000000.0;
      // For-loop: All the EGs in the event
      for (auto eg = L1EGs.begin(); ( !stopClustering && (eg != L1EGs.end()) ); eg++) {
	counter_allEG++;
	// Skip small-Et EGs and those not matching to lead track (in terms of DeltaR)
	// TODO: Correct eta of EG based on vertex position of the leading track?
	if (eg->getEt() < minEt_EG) continue;
	counter_passEt++;
	deltaR = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), CorrectedEta(eg->getEta(), leadTrackPtr->getZ0()), eg->getPhi());
	deltaR_beforecorrection = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), eg->getEta(), eg->getPhi());
	//cout << "eta = "<<eg->getEta()<<"  corrected eta = " <<CorrectedEta(eg->getEta(), leadTrackPtr->getZ0())<<endl;
	// Fill  dR histogram
	h_leadTrk_EG_dR -> Fill (deltaR); 
	h_leadTrk_EG_dR_beforecorrection -> Fill(deltaR_beforecorrection);
	if (deltaR < deltaRmin) deltaRmin=deltaR;
	if (cfg_DEBUG) std::cout << "deltaR = " << deltaR << endl;
	// Skip if EG is not within the signal cone 
	if (deltaR > maxDeltaR_EG) continue;
	counter_passDRmax++;
	if (deltaR < minDeltaR_EG) continue;
	counter_passDRmin++;

	p4sum += eg->p4();
	
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
	  //cout << "Invariant mass of cluster is " << invMass << " and ET is " << ET << endl; //marina
	}
	else {
	  stopClustering = true;
	  // Fill EG cluster variables to histograms (if >=1 EG is clustered)
	  if (EGcluster.size() >=1) {
	    h_EGClusters_Et->Fill( ET );            
	    h_EGClusters_M->Fill( invMass );
	  }
	}
      }// For-loop: All the EGs in the event
      
      h_leadTrk_EG_dRmin -> Fill (deltaRmin);

      // Fill EG clustering counter histogram
      h_clustEGs_allEGs      -> Fill(counter_allEG);
      h_clustEGs_passEt      -> Fill(counter_passEt);
      h_clustEGs_passDRmax   -> Fill(counter_passDRmax);
      h_clustEGs_passDRmin   -> Fill(counter_passDRmin);
      h_clustEGs_passInvMass -> Fill(counter_passInvMass);
      
      // Fill the number of EGs in cluster
      h_EGClusters_MultiplicityPerCluster->Fill( EGcluster.size() );
      
      // Find the genTau matching to the lead track
      GenParticle genTau;
      double deltaR_match = GetMatchingGenParticle(trackTauCandidates[i][0], &genTau);
      h_MCmatch_dR->Fill(deltaR_match);

      bool hasGenTau = false;
      if (deltaR_match < maxDeltaR_MCmatch) hasGenTau = true;

      // marina
      if (hasGenTau) {
	h_leadTrk_MCmatch -> Fill (hasGenTau);
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

      /*
      // Fill histograms with the matched genTau
      if (hasGenTau){
	//Increase counter for trackTauCandidates which have a genTau
	counter_hasGenTau++;
	if (EGcluster.size() >=1) counter_hasEG++;
	
	// ET of TkEG when it has a matched genTau
	hTkEG_matched_Et  -> Fill( ET );

	// Properties of the genTau
	hTkEG_genVisEt    -> Fill( genTau.p4vis().Et() );
	hTkEG_DeltaRmatch -> Fill( deltaR_match );
	// Pt resolution
	TLorentzVector p4tracks; // initialized to (0,0,0,0) 
	for (auto tk = trackTauCandidates[i].begin(); tk != trackTauCandidates[i].end(); tk++) {
	  p4tracks += tk->p4();
	}
	double pTresolution = ( p4tracks.Pt()-genTau.p4vis().Pt() ) / genTau.p4vis().Pt();
	h_trkClusters_PtResolution->Fill(pTresolution);
	// ET resolution
	double ETresolution = ( ET -genTau.p4vis().Et() ) / genTau.p4vis().Et(); //marina
	h_EGClusters_EtResolution->Fill(ETresolution);
      }
      
      // Build a tau candidate from tracks and EGs
      L1TkEGParticle newTauCandidate(trackTauCandidates[i], EGcluster, genTau, hasGenTau);
      if (cfg_DEBUG) { //marina
	cout << "Constructed a new tau candidate with a track invariant mass of " << newTauCandidate.GetTrackInvMass() 
	     << ", an EG invariant mass of " << newTauCandidate.GetEGInvMass();
	if (hasGenTau) cout  << " and a generator tau with visible pT " << newTauCandidate.GetGenTauPt() << endl;
	else cout << " and does not have a matching generator tau" << endl;
      }
      
      
      TauCandidates.push_back(newTauCandidate);
      */
    }

    /*
    // Check relative isolation of tau candidates
    float relIso = -1.0;
    float ptSum;
    TTTrack leadingTrack;
    TauCandidatesIsolated.clear();
    
    for (auto tkeg = TauCandidates.begin(); tkeg != TauCandidates.end(); tkeg++) {
      leadingTrack = tkeg->GetLeadingTrack();
      
      // Sum Pt's inside the isolation cone
      ptSum = 0.0;
      for (auto tk = TTTracks.begin(); tk != TTTracks.end(); tk++) {
	bool clustered = false;
	//for (auto clusttk = tkeg->GetTracks().begin(); clusttk != tkeg->GetTracks().end(); clusttk++){
	//if (clusttk->index() == tk->index()) clustered = true;
	//}
	//if (clustered) continue;
	if (abs (tk->getZ0() - leadingTrack.getZ0() ) > maxDeltaZ_iso) continue;
	deltaR = auxTools_.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), tk->getEta(), tk->getPhi());
        if (deltaR > minDeltaR_iso && deltaR < maxDeltaR_iso)
          ptSum += tk -> getPt();
      }
      // Calculate relative isolation
      relIso = ptSum / tkeg->GetTrackBasedPt();
      // Fill relative isolation histogram
      h_trkClusters_relIso->Fill( relIso );
      //cout << "relIso = " << relIso << endl;
      if(relIso < maxRelIso) {
        TauCandidatesIsolated.push_back(*tkeg);
      }
    }
    */
    /*
    // Fill MCmatch counters 
    h_MCmatch_counters -> SetBinContent(1, counter_allTrkTauCand);
    h_MCmatch_counters -> SetBinContent(2, counter_hasGenTau    );
    h_MCmatch_counters -> SetBinContent(3, counter_hasEG        );
    */
    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////
    
    // Fill visible Et of hadronic taus to a histogram
    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++) hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() );
    //for (vector<GenParticle>::iterator tau = GenTausTrigger.begin(); tau != GenTausTrigger.end(); tau++) hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() );


    // Fill turn-on numerators
    FillTurnOn_Numerator_(TauCandidates, 25.0, hTurnOn25_all);
    FillTurnOn_Numerator_(TauCandidatesIsolated, 25.0, hTurnOn25_relIso);
    FillTurnOn_Numerator_(TauCandidates, 50.0, hTurnOn50_all);
    FillTurnOn_Numerator_(TauCandidatesIsolated, 50.0, hTurnOn50_relIso);

    ////////////////////////////////////////////////
    // Rates and efficiencies for single tau
    ////////////////////////////////////////////////
    FillSingleTau_(TauCandidates, hRateSingleTau_all, hEffSingleTau_all);
    FillSingleTau_(TauCandidates, hRateSingleTau_C, hEffSingleTau_C, 0.0, 1.0);
    FillSingleTau_(TauCandidates, hRateSingleTau_I, hEffSingleTau_I, 1.0, 1.6);
    FillSingleTau_(TauCandidates, hRateSingleTau_F, hEffSingleTau_F, 1.6, 3.0); // 2.5 is max

    FillSingleTau_(TauCandidatesIsolated, hRateSingleTau_relIso, hEffSingleTau_relIso);
    

    ////////////////////////////////////////////////
    // Rates and efficiencies for ditau
    ////////////////////////////////////////////////
    FillDiTau_(TauCandidates, hRateDiTau_all, hEffDiTau_all);
    FillDiTau_(TauCandidates, hRateDiTau_C, hEffDiTau_C, 0.0, 1.0);
    FillDiTau_(TauCandidates, hRateDiTau_I, hEffDiTau_I, 1.0, 1.6);
    FillDiTau_(TauCandidates, hRateDiTau_F, hEffDiTau_F, 1.6, 3.0); // 2.5 is max

    FillDiTau_(TauCandidatesIsolated, hRateDiTau_relIso, hEffDiTau_relIso);

    ////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////
    //ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_Tk);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk  , hDiTau_Eff_Tk);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_C, hDiTau_Eff_Tk_C, 0.0, 1.0);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_I, hDiTau_Eff_Tk_I, 1.0, 1.6);
    //FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_F, hDiTau_Eff_Tk_F, 1.6, 3.0); 
    
    // Progress bar
    if (!cfg_DEBUG) auxTools_.ProgressBar(jentry, nEntries, 100, 100);
   
  } // For-loop: Entries

  ////////////////////////////////////////////////
  // Fill counters
  ////////////////////////////////////////////////

  hCounters->SetBinContent(1, nAllEvts);
  hCounters->SetBinContent(2, nEvts);

  
  
  ////////////////////////////////////////////////
  // Convert/Finalise histograms
  ////////////////////////////////////////////////

  // Divide turn-on numerators by the denumerator
  histoTools_.DivideHistos_1D(hTurnOn25_all, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTurnOn25_relIso, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTurnOn50_all, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTurnOn50_relIso, hMcHadronicTau_VisEt); 
  
  // Convert rate histograms
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_all, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_F, nEntries);

  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_relIso, nEntries);

  histoTools_.ConvertToRateHisto_1D(hRateDiTau_all, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateDiTau_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateDiTau_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateDiTau_F, nEntries);

  histoTools_.ConvertToRateHisto_1D(hRateDiTau_relIso, nEntries);

  // Finalise efficiency histograms
  FinaliseEffHisto_(hEffSingleTau_all  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_(hEffSingleTau_relIso  , nEvtsWithMaxHTaus);

  FinaliseEffHisto_(hRateDiTau_all  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hRateDiTau_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hRateDiTau_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hRateDiTau_F, nEvtsWithMaxHTaus);

  FinaliseEffHisto_(hRateDiTau_relIso  , nEvtsWithMaxHTaus);

  ////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////
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
void TkEG::BookHistos_(void)
//============================================================================
{
  // Event-Type Histograms
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 2, 0.0, +2.0);

  // Number of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_N, "genTausAll_N", ";Number of all genuine taus in event;Events / bin", 5, -0.5, +4.5);

  // Pt of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_Pt, "genTausAll_Pt", ";p_{T} (GeV);Events / bin",100, +0.0, +100.0);

  // Eta of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_Eta, "genTausAll_Eta", ";#eta ;Events / bin", 60, -3.0, +3.0);

  // Phi of all genuine  taus (per event)
  histoTools_.BookHisto_1D(h_genTausAll_Phi, "genTausAll_Phi", ";#phi (rads) ;Events / bin", 21,  -3.15,  +3.15);

  // Eta of first tau Vs Eta of the second tau
  histoTools_.BookHisto_2D(h_genTausAll_Eta1VsEta2, "genTausAll_Eta1VsEta2", ";#eta_{1}; #eta_{2};Events / bin", 60, -3.0, +3.0, 60, -3.0, +3.0);

  // Phi of first tau Vs Phi of the second tau
  histoTools_.BookHisto_2D(h_genTausAll_Phi1VsPhi2, "genTausAll_Phi1VsPhi2", ";#phi_{1}; #phi_{2};Events / bin", 21,  -3.15,  +3.15, 21,  -3.15,  +3.15);

  // Number of daughters of genuine taus
  histoTools_.BookHisto_1D(h_genTausHad_Daughters_N, "genTausHad_Daughters_N", ";N_{daughters} (Hadronic Gen taus); Events / bin", 11, -0.5, +10.5);

  // Number of charged daughters of genuine taus
  histoTools_.BookHisto_1D(h_genTausHad_chargedDaugh_N, "genTausHad_chargedDaugh_N", ";N_{daughters^{#pm}} (Hadronic Gen taus); Events / bin", 11, -0.5, +10.5);

  // Number of neutral daughters of genuine taus
  histoTools_.BookHisto_1D(h_genTausHad_neutralDaugh_N, "genTausHad_neutralDaugh_N", ";N_{daughters^{0}} (Hadronic Gen taus); Events / bin", 11, -0.5, +10.5);

  // Number of genuine hadronic taus (per event)
  histoTools_.BookHisto_1D(h_genTausHad_N, "genTausHadronic_N", ";Number of genuine hadronic taus in event;Events / bin", 5, -0.5, +4.5);

  // Pt of the charged products of the tau decay
  histoTools_.BookHisto_1D(h_genTau_chargedDaugh_Pt, "genTau_chargedDaugh_Pt", ";p_{T} (charged decay products) (GeV);Particles / bin", 300, +0.0, +300.0);

  // TotalMass of the charged products of the tau decay 
  histoTools_.BookHisto_1D(h_genTau_chargedDaugh_totalMass, "genTau_chargedDaugh_totalMass", ";M_{total} (charged decay products); GenTaus/ bin", 40, 0.0, +4.0);

  // TotalMass of the neutral products of the tau decay 
  histoTools_.BookHisto_1D(h_genTau_neutralDaugh_totalMass, "genTau_neutralDaugh_totalMass", ";M_{total} (neutral decay products); GenTaus/ bin", 40, 0.0, +4.0);

  // Et of the neutral products of the tau decay
  histoTools_.BookHisto_1D(h_genTau_neutralDaugh_Et, "genTau_neutralDaugh_Et", ";E_{T} (neutral decay products) (GeV);Particles / bin", 300, +0.0, +300.0);

  // Pt of charged products Vs Et of neutral products of hadronic tau
  histoTools_.BookHisto_2D(h_genTauHad_chargedPtVsneutralET, "genTauHad_chargedPtVsneutralET", ";p_{T}^{charged} (GeV); E_{T}^{neutral} (GeV);Events / bin", 120, +0.0, +120.0, 20, +0.0, +120.0);
  
  // Visible pt of charged decay products (leading pt>10GeV) Vs dRmax
  histoTools_.BookHisto_2D(h_genTau_chargedDaugh_visPt_dRmax, "genTau_chargedDaugh_visPt_dRmax", ";p_{T}^{vis} (#tau)  (GeV); #DeltaR_{max} (charged^{ldg}, charged)", 100, 0.0, 100.0, 50, 0.0, 0.5);

  // Pt of leading charged product (leading pt>10GeV) Vs dRmax with other charged products
  histoTools_.BookHisto_2D(h_genTau_chargedDaugh_PtLead_dRmax, "genTau_chargedDaugh_PtLead_dRmax", ";p_{T} (charged^{ldg}) (GeV); #DeltaR_{max} (charged^{ldg}, charged)", 100, 0.0, 100.0, 50, 0.0, 0.5);

  // Pt of leading charged decay product  (leading pt>10GeV) Vs dRmax with  other neutral products
  histoTools_.BookHisto_2D(h_genTau_neutralDaugh_PtLead_dRmax, "genTau_neutralDaugh_PtLead_dRmax", ";p_{T} (charged^{ldg}) (GeV); #DeltaR_{max} (charged^{ldg}, neutral)", 100, 0.0, 100.0, 50, 0.0, 0.5);

  // Number of stubs of all tracks
  histoTools_.BookHisto_1D(h_trk_NStubs_all, "trk_NStubs_all", ";N_{stubs} (all tracks); Events / bin", 16, -0.5, 15.5);

  // Chi Squared of the tracks (4 stubs)
  histoTools_.BookHisto_1D(h_trk_Chi2_all_4stubs, "trk_Chi2_all_4stubs", ";#chi^{2} (4-stubs tracks); Events / bin", 100, 0.0, 100.0);

  // Chi Squared of the tracks (>=5 stubs)
  histoTools_.BookHisto_1D(h_trk_Chi2_all_5stubs, "trk_Chi2_all_5stubs", ";#chi^{2} (5-stubs tracks); Events / bin", 100, 0.0, 100.0);

  // Chi Squared of all tracks
  histoTools_.BookHisto_1D(h_trk_Chi2_all, "trk_Chi2_all", ";#chi^{2} (all tracks); Events / bin", 100, 0.0, 100.0);

  // Chi Squared/ dof of all tracks
  histoTools_.BookHisto_1D(h_trk_Chi2Red_all, "trk_Chi2Red_all", ";#chi^{2}/dof (all tracks); Events / bin", 50, 0.0, 50.0);

  // Number of lead tracks (per event)
  histoTools_.BookHisto_1D(h_leadTrks_Multiplicity, "leadTrks_Multiplicity", ";Number of lead tracks in event;Events / bin", 15, -0.5, +14.5);

  // Lead track Pt
  histoTools_.BookHisto_1D(h_leadTrks_Pt, "leadTrks_Pt", ";p_{T} (GeV);Tracks / bin", 300, +0.0, +300.0);

  // Lead track Eta
  histoTools_.BookHisto_1D(h_leadTrks_Eta, "leadTrks_Eta", ";#eta;Tracks / bin", 60, -3.0, +3.0);

  // Lead track Phi
  histoTools_.BookHisto_1D(h_leadTrks_Phi, "leadTrks_Phi", ";#phi (rads);Tracks / bin", 21,  -3.15,  +3.15);
  
  // Lead tracks in Eta-Phi plane
  // (Syntax: BookHisto_2D(histogram, hName, hTitle, binsX, xMin, xMax, binsY, yMin, yMax)
  histoTools_.BookHisto_2D(h_leadTrks_Phi_Eta, "leadTrks_Phi_Eta",  ";#phi (rads);#eta", 210,  -3.15,  +3.15, 600,  -3.0,  +3.0);  

  // Counters for the selection of lead tracks
  histoTools_.BookHisto_1D(h_leadTrkSelection, "leadTrkSelection", "", 5,  0.0,  5.0);  

  // Clustered tracks Pt
  histoTools_.BookHisto_1D(h_clustTrks_Pt, "clustTrks_Pt", ";p_{T} (GeV);Tracks / bin", 300, +0.0, +300.0);

  // Clustered tracks Eta
  histoTools_.BookHisto_1D(h_clustTrks_Eta, "clustTrks_Eta", ";#eta;Tracks / bin", 60, -3.0, +3.0);

  // Clustered tracks Phi
  histoTools_.BookHisto_1D(h_clustTrks_Phi, "clustTrks_Phi", ";#phi (rads);Tracks / bin", 21,  -3.15,  +3.15);
  
  // Clustered tracks in Eta-Phi plane
  histoTools_.BookHisto_2D(h_clustTrks_Phi_Eta, "clustTrks_Phi_Eta",  ";#phi (rads);#eta", 210,  -3.15,  +3.15, 300,  -3.0,  +3.0);  			      

  // Track clustering counter
  histoTools_.BookHisto_2D(h_clustTrks_counter, "clustTrks_counter", ";Selection steps;clustered tracks / selection step", 6, 0., 6., 150, 0., 300.);
  h_clustTrks_counter->GetXaxis()->SetBinLabel(1,"allTracks");
  h_clustTrks_counter->GetXaxis()->SetBinLabel(2,"allNonLeading");
  h_clustTrks_counter->GetXaxis()->SetBinLabel(3,"passZ");
  h_clustTrks_counter->GetXaxis()->SetBinLabel(4,"passDRmax");
  h_clustTrks_counter->GetXaxis()->SetBinLabel(5,"passDRmin");
  h_clustTrks_counter->GetXaxis()->SetBinLabel(6,"passInvMass");
  h_clustTrks_counter->SetMarkerStyle(24);

  // Number of clustered tracks (per cluster)
  histoTools_.BookHisto_1D(h_trkClusters_MultiplicityPerCluster, "trkClusters_MultiplicityPerCluster", ";Number of tracks in cluster;Clusters / bin", 16, -0.5, +15.5);

  histoTools_.BookHisto_1D(h_MCmatch_chargedDaugh_N, "MCmatch_chargedDaugh_N", ";Number of daughters^{+} of matched gen-#tau; entries / bin", 16, -0.5, +15.5);

  histoTools_.BookHisto_1D(h_MCmatch_neutralDaugh_N, "MCmatch_neutralDaugh_N", ";Number of daughters^{0} of matched gen-#tau; entries / bin", 16, -0.5, +15.5);

  // Pt of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_Pt, "trkClusters_Pt", ";p_{T} (GeV);Clusters / bin", 300, +0.0, +300.0);

  // Pt resolution of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_PtResolution, "trkClusters_PtResolution", ";p_{T} resolution (GeV);Clusters / bin", 50, -5.0, +5.0);
  
  // Invariant mass of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_M, "trkClusters_M", ";Invariant mass;Clusters / bin", 40, 0.0, +4.0);

  // Number of EGs                                                                                                                                                       
  histoTools_.BookHisto_1D(h_EGs_N, "EGs_Multiplicity", ";Number of EGs in the event ; entries / bin", 21, -0.5, +20.5);

  // EGs Et                                                                                                                                                       
  histoTools_.BookHisto_1D(h_EGs_Et, "EGs_Et", ";E_{T} (GeV) ;EGs / bin", 300, 0.0, +300.0);

  // EGs Eta 
  histoTools_.BookHisto_1D(h_EGs_Eta, "EGs_Eta", ";#eta ;EGs / bin", 64, -1.6, +1.6);

  // EGs Phi
  histoTools_.BookHisto_1D(h_EGs_Phi, "EGs_Phi", ";#phi (rads);EGs / bin", 36,  -3.15,  +3.15);

  // EGs IEta
  histoTools_.BookHisto_1D(h_EGs_IEta, "EGs_IEta", ";i#eta;EGs / bin", 70, -35, +35); 

  // EGs IPhi
  histoTools_.BookHisto_1D(h_EGs_IPhi, "EGs_IPhi", ";i#phi;EGs / bin", 36,  0,  145);

  // DR of lead trk and EGs
  histoTools_.BookHisto_1D(h_leadTrk_EG_dR, "leadTrk_EG_dR", ";#DeltaR(trk_{ldg}, EG); entries / bin",   60,  0.0,   +6.0);

  // DR of lead trk and EGs (before the correction of EG's eta
  histoTools_.BookHisto_1D(h_leadTrk_EG_dR_beforecorrection, "leadTrk_EG_dR_beforecorrection", ";#DeltaR(trk_{ldg}, EG); entries / bin",   60,  0.0,   +6.0);

  // minimum DR of lead trk and EGs
  histoTools_.BookHisto_1D(h_leadTrk_EG_dRmin, "leadTrk_EG_dRmin", ";#DeltaR_{min}(trk_{ldg}, EG); entries / bin",   60,  0.0,   +6.0);

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
  histoTools_.BookHisto_1D(h_EGClusters_Et, "EGClusters_Et", ";Et (GeV);Clusters / bin", 60, +0.0, +300.0);

  // Et resolution of EG clusters
  histoTools_.BookHisto_1D(h_EGClusters_EtResolution, "EGClusters_EtResolution", ";Et resolution (GeV);Clusters / bin", 50, -5.0, +5.0);
  
  // Invariant mass of EG clusters
  histoTools_.BookHisto_1D(h_EGClusters_M, "EGClusters_M", ";Invariant mass (GeV);Clusters / bin", 40, +0.0, +4.0);

  // Track-based relative isolation of tau candidates
  histoTools_.BookHisto_1D(h_trkClusters_relIso, "trkClusters_relIso", ";Relative isolation;Clusters / bin", 100, 0.0, +5.0);
  
  // Et of the tkeg when it has a matched genTau
  histoTools_.BookHisto_1D(hTkEG_matched_Et, "TkEG_matched_Et", ";E_{T} (GeV); TkEGs / bin", 40, 0.0,  +200.0);

  // DeltaR in MC matching of the lead track of tau candidates
  histoTools_.BookHisto_1D(hTkEG_DeltaRmatch  , "TkEG_DeltaRmatch"  , ";#Delta R;Particles / bin",   100,  0.0,   +1.0);
  
  // Visible Et of the matched generator taus
  histoTools_.BookHisto_1D(hTkEG_genVisEt, "TkEG_genVisEt", ";E_{T}^{vis} (gen) (GeV);Particles / bin", 40, 0.0,  +200.0);
  
  // Visible Et of all generator taus
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt, "hMcHadronicTau_VisEt", ";VisEt (GeV);Particles / bin", 40, 0.0,  +200.0);

  histoTools_.BookHisto_1D(h_MCmatch_counters, "MCmatch_counters", ";;Events", 3, 0, 3);
  h_MCmatch_counters->GetXaxis()->SetBinLabel(1, "All TrkClusters");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(2, "MC Matched");
  h_MCmatch_counters->GetXaxis()->SetBinLabel(3, "EG associated");
  
  histoTools_.BookHisto_1D(h_leadTrk_MCmatch, "leadTrk_MCmatch", ";;Events", 2, 0, 2);
  h_leadTrk_MCmatch->GetXaxis()->SetBinLabel(1, "NOT MC Matched");
  h_leadTrk_MCmatch->GetXaxis()->SetBinLabel(2, "MC Matched");

  // minimum Delta R of leading track and gen particle 
  histoTools_.BookHisto_1D(h_MCmatch_dR, "MCmatch_dR", ";#DeltaR_{min}(trk_{ldg}, gen-#tau); entries / bin",   60,  0.0,   +6.0);

  // Turn-on histograms
  histoTools_.BookHisto_1D(hTurnOn25_all, "TurnOn25_all", ";track cluster p_{T} (GeV); Efficiency / bin", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTurnOn25_relIso, "TurnOn25_relIso", ";track cluster p_{T} (GeV); Efficiency / bin", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTurnOn50_all, "TurnOn50_all", ";track cluster p_{T} (GeV); Efficiency / bin", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTurnOn50_relIso, "TurnOn50_relIso", ";track cluster p_{T} (GeV); Efficiency / bin", 40, 0.0,  +200.0);
  
  // Single-tau rates
  histoTools_.BookHisto_1D(hRateSingleTau_all, "Rate_SingleTau_all"    , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateSingleTau_C  , "Rate_SingleTau_C"  , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateSingleTau_I  , "Rate_SingleTau_I"  , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateSingleTau_F  , "Rate_SingleTau_F"  , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);  
  
  histoTools_.BookHisto_1D(hRateSingleTau_relIso  , "Rate_SingleTau_RelIso"    , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);

  // Di-tau rates
  histoTools_.BookHisto_1D(hRateDiTau_all, "Rate_DiTau_all"    , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateDiTau_C  , "Rate_DiTau_C"  , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateDiTau_I  , "Rate_DiTau_I"  , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateDiTau_F  , "Rate_DiTau_F"  , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);  

  histoTools_.BookHisto_1D(hRateDiTau_relIso    , "Rate_DiTau_RelIso"    , ";track cluster p_{T} threshold (GeV); Rate (kHz) / bin", 200, 0.0,  +200.0);
  
  // Single-tau efficiencies
  histoTools_.BookHisto_1D(hEffSingleTau_all , "Eff_SingleTau_all"     , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffSingleTau_C   , "Eff_SingleTau_C"   , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffSingleTau_I   , "Eff_SingleTau_I"   , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffSingleTau_F   , "Eff_SingleTau_F"   , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);

  histoTools_.BookHisto_1D(hEffSingleTau_relIso     , "Eff_SingleTau_RelIso"     , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);

  // Di-tau efficiencies
  histoTools_.BookHisto_1D(hEffDiTau_all , "Eff_DiTau_all"     , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffDiTau_C   , "Eff_DiTau_C"   , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffDiTau_I   , "Eff_DiTau_I"   , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffDiTau_F   , "Eff_DiTau_F"   , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);

  histoTools_.BookHisto_1D(hEffDiTau_relIso     , "Eff_DiTau_RelIso"     , ";track cluster p_{T} threshold (GeV); Efficiency / bin", 200, 0.0,  +200.0);
  
  return;
}


//============================================================================
void TkEG::WriteHistos_(void)
//============================================================================
{
  // Location -> outFile
  outFile->cd();
  
  // Write Counters
  hCounters->Write();

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
  h_genTau_chargedDaugh_visPt_dRmax->Write();
  h_genTau_chargedDaugh_PtLead_dRmax->Write();
  h_genTau_neutralDaugh_PtLead_dRmax->Write();

  h_trk_NStubs_all->Write();
  h_trk_Chi2_all_4stubs->Write();
  h_trk_Chi2_all_5stubs->Write();
  h_trk_Chi2_all->Write();
  h_trk_Chi2Red_all->Write();
  
  h_leadTrks_Multiplicity->Write();

  h_leadTrks_Pt->Write();
  h_leadTrks_Eta->Write();
  h_leadTrks_Phi->Write();

  h_leadTrkSelection->Write();

  h_clustTrks_Pt->Write();
  h_clustTrks_Eta->Write();
  h_clustTrks_Phi->Write();
  //h_clustTrks_counter->Write();

  h_trkClusters_MultiplicityPerCluster->Write();
  h_MCmatch_chargedDaugh_N->Write();
  h_MCmatch_neutralDaugh_N->Write();
  h_trkClusters_Pt->Write();
  h_trkClusters_PtResolution->Write();
  h_trkClusters_M->Write();

  h_clustEGs_Et->Write();
  h_clustEGs_Eta->Write();
  h_clustEGs_Phi->Write();
  h_clustEGs_M->Write();

  h_EGs_N -> Write();
  h_EGs_Et->Write();
  h_EGs_Eta->Write();
  h_EGs_Phi->Write();
  h_EGs_IEta->Write();
  h_EGs_IPhi->Write();

  h_clustEGs_allEGs      -> Write();
  h_clustEGs_passEt      -> Write();
  h_clustEGs_passDRmax   -> Write();
  h_clustEGs_passDRmin   -> Write();
  h_clustEGs_passInvMass -> Write();

  h_leadTrk_EG_dR->Write();
  h_leadTrk_EG_dR_beforecorrection->Write();
  h_leadTrk_EG_dRmin->Write();

  h_EGClusters_MultiplicityPerCluster->Write();
  h_EGClusters_Et->Write();
  h_EGClusters_EtResolution->Write();
  h_EGClusters_M->Write();

  h_trkClusters_relIso->Write();

  hTkEG_matched_Et->Write();
  hTkEG_DeltaRmatch->Write();
  hTkEG_genVisEt->Write();
  hMcHadronicTau_VisEt->Write();
  
  h_MCmatch_counters->Write(); 
  h_MCmatch_dR->Write(); 
  h_leadTrk_MCmatch->Write(); 

  hTurnOn25_all->Write();
  hTurnOn25_relIso->Write();
  hTurnOn50_all->Write();
  hTurnOn50_relIso->Write();

  hRateSingleTau_all->Write(); // Inclusive = C+I+F                                                                                                                   
  hRateSingleTau_C->Write();
  hRateSingleTau_I->Write();
  hRateSingleTau_F->Write();

  hRateSingleTau_relIso->Write(); // Inclusive = C+I+F                                                                                                                   

  hRateDiTau_all->Write(); // Inclusive = C+I+F                                                                                                                       
  hRateDiTau_C->Write();
  hRateDiTau_I->Write();
  hRateDiTau_F->Write();

  hRateDiTau_relIso->Write(); // Inclusive = C+I+F                                                                                                                       

  hEffSingleTau_all->Write();  // Inclusive = C+I+F                                                                                                                   
  hEffSingleTau_C->Write();
  hEffSingleTau_I->Write();
  hEffSingleTau_F->Write();

  hEffSingleTau_relIso->Write();  // Inclusive = C+I+F                                                                                                                   

  hEffDiTau_all->Write();  // Inclusive = C+I+F                                                                                                                       
  hEffDiTau_C->Write();
  hEffDiTau_I->Write();
  hEffDiTau_F->Write();

  hEffDiTau_relIso->Write();  // Inclusive = C+I+F                                                                                                                       
  
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
  // if (GenTausHadronic.size() < 1 ) return; // FIXME: move outside this function

  // Initialise the GenParticle (to be returned)
  GenParticle match_GenParticle;
  double deltaR;
  double match_dR = 999;
  
  // For-loop: All hadronic GenTaus
  for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++) {
  //for (vector<GenParticle>::iterator tau = GenTausTrigger.begin(); tau != GenTausTrigger.end(); tau++) {
      // If no hadronic decay products found (pi+/-, pi0, K+/-, K0, K0L), skip this tau
      if (tau->finalDaughtersCharged().size() < 1) continue;
      deltaR = auxTools_.DeltaR( track.getEta(), track.getPhi(), tau->eta(), tau->phi() );
      //if (deltaR > maxDeltaR_MCmatch) continue;
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
  double ldgEt = L1TkEGs.at(0).GetTrackBasedPt();

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
  double ldgEt_mcMatched = L1TkEGs_mcMatched.at(0).GetTrackBasedPt();
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
  double subLdgEt = L1TkEG.GetTrackBasedPt();
  
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
  double subLdgEt_mcMatched = L1TkEGs_mcMatched.at(1).GetTrackBasedPt();
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
  double ldgEt1 = L1TkEGs1.at(0).GetTrackBasedPt();
  double ldgEt2 = L1TkEGs2.at(0).GetTrackBasedPt();

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
  double ldgEt1_mcMatched = L1TkEGs1_mcMatched.at(0).GetTrackBasedPt();
  double ldgEt2_mcMatched = L1TkEGs2_mcMatched.at(0).GetTrackBasedPt();

  
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
				   const double minPt,
				   TH1D *hTurnOn)
//============================================================================
{
  
  // For-loop: L1TkEGs
  int iterated_particles = 0;
  int matched_particles = 0;
  int passed_cut = 0;
  for (vector<L1TkEGParticle>::iterator tau = L1TkEGs.begin(); tau != L1TkEGs.end(); tau++)
    {
      
      iterated_particles++;
      // Skip if trigger object is not MC matched
      if (!tau->HasMatchingGenParticle()) continue;	 
      
      matched_particles++;
      // Skip if trigger object has Pt < minPt
      if (tau->GetTrackBasedPt() < minPt) continue;
      passed_cut++;      
            
      // Fill the turn-on
      hTurnOn->Fill( tau->GetGenTauEt() );
      // DEBUG

    } // For-loop: L1TkEGs
//      cout << "Filling turn-on histogram: iterated particles=" << iterated_particles << ", matched=" << matched_particles << ", passed_cut=" << passed_cut << endl;

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

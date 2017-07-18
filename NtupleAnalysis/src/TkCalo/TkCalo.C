#ifndef TkCalo_cxx
#define TkCalo_cxx

// User
//#include "../Auxiliary/interface/constants.h"
#include "TkCalo.h"

// ROOT
//#include "TFitResult.h"
//#include "TF1.h"

// C++
#include <stdexcept>

//============================================================================
void TkCalo::InitObjects(void)
//============================================================================
{
  return; 
}


//============================================================================
void TkCalo::InitVars_()
//============================================================================
{

  datasets_  = datasets_.GetDataset(mcSample);
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;
  
  cfg_DEBUG = true;
  if (cfg_DEBUG) std::cout << "=== TkCalo::InitVars_()" << std::endl;
  
  cfg_AddL1Tks   = true;
  cfg_AddEGs     = true;
  cfg_AddGenP    = true;
  
  // Track parameters
  cfg_tk_Collection  =  "TTTracks"; // Default: "TTTracks" (not "TTPixelTracks")
  cfg_tk_nFitParams  = 4;           // Default: 4
  cfg_tk_minPt       = 2.00;        // Default: 2.0
  cfg_tk_minEta      = 0.0;         // Default: 0.0
  cfg_tk_maxEta      = 1e6;         // Default: 1e6
  cfg_tk_maxChiSqRed = 1e6;         // Default: 1e6
  cfg_tk_minStubs    =   0;         // Default: 0
  cfg_tk_minStubsPS  =   0;         // Default: 0
  cfg_tk_maxStubsPS  = 1e6;         // Default: 1e6

  // TkCalo algorithm parameters
  ETmin = 5.0; // min ET in GeV of L1EG objects
  ZMAX = 10.0; // |z_track| < ZMAX in cm
  CHI2MAX = 200;		
  DRmin=0.1;
  DRmax=0.3;
  minPt_leading_track = 10.0;
  PrimaryVtxConstrain = false; // use the primary vertex
  DeltaZMax = 1.0; // | z_track - z_primaryvtx | < DeltaZMax in cm. 
  IsoCut=0.15; 
  RelativeIsolation=true;

  // New parameters
  minStubs_trk      = 5; 
  maxChi2_trk       = 40.0; // GeV
  minPt_leadtrk     = 10.0; // GeV
  maxEta_leadtrk    = 2.3;
  minDeltaR_leadtrk = 0.0;
  maxDeltaR_leadtrk = 0.3;
  maxDeltaZ_trk     = 1.0;  // cm
  maxInvMass_trk    = 1.77; // GeV 
  minEt_EG          = 5.0;  // GeV
  minDeltaR_EG      = 0.0;
  maxDeltaR_EG      = 0.3;
  maxInvMass_EG     = 1.77; // GeV
  maxDeltaR_MCmatch = 0.1;    
  minDeltaR_iso     = 0.0;
  maxDeltaR_iso     = 5.0;
  maxDeltaZ_iso     = 0.8;  // cm
  useRelIso         = true;
  maxRelIso         = 0.15;
  
  return;

}


//============================================================================
void TkCalo::PrintSettings(void)
//============================================================================
{

  // TODO
  
  return;
}


//============================================================================
void TkCalo::Loop()
//============================================================================
{

  // Sanity check
  if (fChain == 0) return;
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  cout << "=== TkCalo:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;
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
  
    if(cfg_DEBUG) cout << "\tEntry = " << jentry << endl;
    
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
	  if (cfg_DEBUG) cout << "\n=== TTracks (" << L1Tks_Pt->size() << ")" << endl;
	  TTTracks = GetTTTracks(cfg_tk_minPt, cfg_tk_minEta, cfg_tk_maxEta, cfg_tk_maxChiSqRed, cfg_tk_minStubs, cfg_tk_minStubsPS, cfg_tk_maxStubsPS, cfg_tk_nFitParams, false);
	  sort( TTTracks.begin(), TTTracks.end() ); // Sort from highest Pt to lowest (not done by default)
	}    
	
	// EGs
	if (cfg_AddEGs) {
	  L1EGs = GetL1EGs(false);
	  sort( L1EGs.begin(), L1EGs.end() ); // Sort from highest Et to lowest Et (should be already done by default)
	  cout << "Sorted EGs!" << endl;
	  if (cfg_DEBUG) cout << "\n=== L1EGs (" << L1EGs.size() << ")" << endl;
	}

    // GenParticles
    vector<GenParticle> GenTaus;
	GenTaus = GetGenParticles(15, false);
	
	// Hadronic GenTaus
    GenTausHadronic.clear();
    if (cfg_AddGenP) {
      if (cfg_DEBUG) cout << "\n=== GenParticles (" << GenP_Pt->size() << ")" << endl;
	  GenTausHadronic = GetHadronicGenTaus(GenTaus, 00.0, 999.9);
	}

    // Triggred GenTaus
    vector<GenParticle> GenTausTrigger;  
    GenTausTrigger      = GetHadronicGenTaus(GenTaus, 20.0, 2.4);  
    
    // Ensure that all taus are found //TODO: what is the point?
    bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;

    ////////////////////////////////////////////////
    // TkCalo algorithm
    ////////////////////////////////////////////////    

    // Consider only events with at least one genuine hadronic tau
    h_genTausHad_N->Fill( GenTausHadronic.size() );
    if (GenTausHadronic.size() < 1) continue;

    // Select lead tracks  
    vector<unsigned short> leadTrackIndices;
    trackTauCandidates.clear();
    float deltaR;
    bool highPtNeighbourFound;
    vector< TTTrack> newCandidateTracks; // temporary container
    for (vector<TTTrack>::iterator leadTrkIter = TTTracks.begin(); leadTrkIter != TTTracks.end(); leadTrkIter++) {
      if (cfg_DEBUG*0) std:: cout << "Considering leading track candidate " << leadTrkIter->index() << std::endl;
      // Only use high quality tracks
      if (leadTrkIter->getNumOfStubs() < minStubs_trk) continue;
      if (leadTrkIter->getChi2() > maxChi2_trk) continue;
      // Pt and eta cuts
      if (leadTrkIter->getPt() < minPt_leadtrk) continue;
      if (leadTrkIter->getEta() > maxEta_leadtrk) continue;
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
    for (size_t i=0; i<trackTauCandidates.size(); i++) {
        leadTrackPtr = &(trackTauCandidates[i][0]);
        if (cfg_DEBUG*0) cout << "Starting to cluster lead track " << leadTrackPtr->index();
        // Loop over other tracks
        stopClustering = false;
        for (vector<TTTrack>::iterator trackIter = TTTracks.begin(); !stopClustering && trackIter != TTTracks.end(); trackIter++) {
          // Skip tracks that are not close in terms of deltaR and deltaZ
          if (abs (trackIter->getZ0() - leadTrackPtr->getZ0() ) > maxDeltaZ_trk) continue;
          deltaR = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), trackIter->getEta(), trackIter->getPhi());
          if (deltaR > maxDeltaR_leadtrk) continue;
          if (deltaR < minDeltaR_leadtrk) continue;
          // Calculate invariant mass of already clustered tracks + new track
          invMass = 0.0;
          pT = 0.0;
          TLorentzVector p4sum; // initialized to (0,0,0,0)
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
          }  
          else {
            // Stop clustering and fill track cluster histograms
            stopClustering = true;
            h_trkClusters_MultiplicityPerCluster->Fill( trackTauCandidates[i].size() );
            h_trkClusters_Pt->Fill( pT );
            h_trkClusters_M->Fill( invMass );
          }
        }
        if (cfg_DEBUG) cout << "After clustering trackTauCandidates[" << i << "] with leadTrack " << leadTrackPtr->index() << ", it has " 
                            << trackTauCandidates[i].size() << " tracks and a mass of " << invMass << endl;

    }

    // Build EG clusters and create tau candidates
    vector<L1EG> EGcluster;
    TauCandidates.clear();
    for (size_t i=0; i<trackTauCandidates.size(); i++) {
        std::cout << "Clustering EGs around trackTauCandidates[" << i << "]" << endl;
        EGcluster.clear();
        pT = 0.0;
        invMass = 0.0;
        leadTrackPtr = &(trackTauCandidates[i][0]);
        stopClustering = false;
        TLorentzVector p4sum; // initialized to (0,0,0,0)        
	    for (auto eg = L1EGs.begin(); !stopClustering && eg != L1EGs.end(); eg++) {
          // Skip small-Et EGs and those not matching to lead track (in terms of DeltaR)
          // TODO: Correct eta of EG based on vertex position of the leading track?
          if (eg->getEt() < minEt_EG) continue;
          deltaR = auxTools_.DeltaR(leadTrackPtr->getEta(), leadTrackPtr->getPhi(), eg->getEta(), eg->getPhi());
          if (cfg_DEBUG*0) std::cout << "deltaR = " << deltaR << endl;
          if (deltaR > maxDeltaR_EG) continue;
          if (deltaR < minDeltaR_EG) continue;
          p4sum += eg->p4();
          std::cout << "invMassTmp = " << p4sum.M() << endl;
          if (p4sum.M() < maxInvMass_EG){
            EGcluster.push_back(*eg);
            pT = p4sum.Pt();
            invMass = p4sum.M();
            // Fill EG histograms
            h_clustEGs_Et->Fill( eg->getEt() );
            h_clustEGs_Eta->Fill( eg->getEta() );
            h_clustEGs_Phi->Fill( eg->getPhi() );
            h_clustEGs_Phi_Eta->Fill( eg->getPhi(),eg->getEta() );
            // Fill EG cluster histograms
            h_EGClusters_Pt->Fill( pT );            
            h_EGClusters_M->Fill( invMass );            
          }
          else {
            stopClustering = true;
          }
        }
        // Fill the number of EGs in cluster
        h_clustEGs_MultiplicityPerCluster->Fill( EGcluster.size() );

        // Find the genTau matching to the lead track
        GenParticle genTau;
        double deltaR_match = GetMatchingGenParticle(trackTauCandidates[i][0], &genTau);

        bool hasGenTau = false;
        if (deltaR_match < maxDeltaR_MCmatch) hasGenTau = true;

        // Fill histograms with the matched genTau
        if (hasGenTau){
          hTkEG_VisEt->Fill( genTau.p4vis().Et() );
          hTkEG_DeltaRmatch->Fill( deltaR_match );
        }

        // Build a tau candidate from tracks and EGs
        L1TkEGParticle newTauCandidate(trackTauCandidates[i], EGcluster, genTau, hasGenTau);
        cout << "Constructed a new tau candidate with a track invariant mass of " << newTauCandidate.GetTrackInvMass() 
             << ", an EG invariant mass of " << newTauCandidate.GetEGInvMass();
        if (hasGenTau) cout  << " and a generator tau with visible pT " << newTauCandidate.GetGenTauPt() << endl;
        else cout << " and does not have a matching generator tau" << endl;
        TauCandidates.push_back(newTauCandidate);
    }  
       
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
        if (abs (tk->getZ0() - leadingTrack.getZ0() ) > maxDeltaZ_iso) continue;
        deltaR = auxTools_.DeltaR(leadingTrack.getEta(), leadingTrack.getPhi(), tk->getEta(), tk->getPhi());
        if (deltaR > minDeltaR_iso && deltaR < maxDeltaR_iso)
          ptSum += tk -> getPt();
      }
      // Calculate relative isolation
      relIso = tkeg->GetTrackPtSum() / ptSum;
      // Fill relative isolation histogram
      h_trkClusters_relIso->Fill( relIso );
      //cout << "relIso = " << relIso << endl;
      if(relIso < maxRelIso) {
        TauCandidatesIsolated.push_back(*tkeg);
      }
    }

    ////////////////////////////////////////////////
    // Print the properties of L1TkEMParticles
    ////////////////////////////////////////////////

    // TODO
    
    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////

    // Fill visible Et of hadronic taus to a histogram
    for (vector<GenParticle>::iterator tau = GenTausHadronic.begin(); tau != GenTausHadronic.end(); tau++) hMcHadronicTau_VisEt->Fill( tau->p4vis().Et() );

    // Fill turn-on numerators
    FillTurnOn_Numerator_(TauCandidates, 25.0, hTurnOn25_all);
    FillTurnOn_Numerator_(TauCandidatesIsolated, 25.0, hTurnOn25_relIso);
    FillTurnOn_Numerator_(TauCandidates, 50.0, hTurnOn50_all);
    FillTurnOn_Numerator_(TauCandidatesIsolated, 50.0, hTurnOn50_relIso);

    ////////////////////////////////////////////////
    // Rates and efficiencies for single tau
    ////////////////////////////////////////////////
    FillSingleTau_(TauCandidates, hRateSingleTau, hEffSingleTau);
    FillSingleTau_(TauCandidates, hRateSingleTau_C, hEffSingleTau_C, 0.0, 1.0);
    FillSingleTau_(TauCandidates, hRateSingleTau_I, hEffSingleTau_I, 1.0, 1.6);
    FillSingleTau_(TauCandidates, hRateSingleTau_F, hEffSingleTau_F, 1.6, 3.0); // 2.5 is max

    ////////////////////////////////////////////////
    // Rates and efficiencies for ditau
    ////////////////////////////////////////////////
/*    FillDiTau_(L1TkTaus_Calo, L1TkTaus_Tk    , hDiTau_Rate_Calo_Tk    , hDiTau_Eff_Calo_Tk    );
    FillDiTau_(L1TkTaus_Calo, L1TkTaus_VtxIso, hDiTau_Rate_Calo_VtxIso, hDiTau_Eff_Calo_VtxIso);
    FillDiTau_(L1TkTaus_Tk  , L1TkTaus_VtxIso, hDiTau_Rate_Tk_VtxIso  , hDiTau_Eff_Tk_VtxIso  );

    ////////////////////////////////////////////////
    // WARNING: Erases L1TkTaus from vector!
    ////////////////////////////////////////////////
    ApplyDiTauZMatching(matchTk_Collection, L1TkTaus_Tk);
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk  , hDiTau_Eff_Tk);
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_C, hDiTau_Eff_Tk_C, 0.0, 1.0);
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_I, hDiTau_Eff_Tk_I, 1.0, 1.6);
    FillDiTau_(L1TkTaus_Tk, hDiTau_Rate_Tk_F, hDiTau_Eff_Tk_F, 1.6, 3.0); */

    // Progress bar
    if (!cfg_DEBUG) auxTools_.ProgressBar(jentry, nEntries, 100, 100);
    
  } // For-loop: Entries

  ////////////////////////////////////////////////
  // Fill counters
  ////////////////////////////////////////////////

  // TODO
  
  
  ////////////////////////////////////////////////
  // Convert/Finalise histogramss
  ////////////////////////////////////////////////

  // Divide turn-on numerators by the denumerator
  histoTools_.DivideHistos_1D(hTurnOn25_all, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTurnOn25_relIso, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTurnOn50_all, hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hTurnOn50_relIso, hMcHadronicTau_VisEt); 
  
  // Convert rate histograms
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_C, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_I, nEntries);
  histoTools_.ConvertToRateHisto_1D(hRateSingleTau_F, nEntries);

  // Finalise efficiency histograms
  FinaliseEffHisto_(hEffSingleTau  , nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_C, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_I, nEvtsWithMaxHTaus);
  FinaliseEffHisto_(hEffSingleTau_F, nEvtsWithMaxHTaus);


  ////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////

  ////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////
  outFile->cd();
  outFile->Write();
  auxTools_.StopwatchStop(5, "minutes", "Total Time");

}


//============================================================================
float TkCalo::DeltaPhi(float phi1, float phi2) 
//============================================================================
{  float dphi = phi1 - phi2;
   if (dphi < 0) dphi = dphi + 2.*TMath::Pi();
   if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
   return dphi;
}


//============================================================================
float TkCalo::deltaR(float eta1, float eta2, float phi1, float phi2) 
//============================================================================
{   float deta = eta1 - eta2;
    float dphi = DeltaPhi(phi1, phi2);
    float DR= sqrt( deta*deta + dphi*dphi );
    return DR;
}


//============================================================================
void TkCalo::BookHistos_(void)
//============================================================================
{

  // Number of genuine hadronic taus (per event)
  histoTools_.BookHisto_1D(h_genTausHad_N, "genTausHadronic_N", ";Number of genuine hadronic taus in event;Events / bin", 4, -0.5, +4.5);

  // Number of lead tracks (per event)
  histoTools_.BookHisto_1D(h_leadTrks_Multiplicity, "leadTrks_Multiplicity", ";Number of lead tracks in event;Events / bin", 30, -0.5, +14.5);

  // Lead track Pt
  histoTools_.BookHisto_1D(h_leadTrks_Pt, "leadTrks_Pt", ";Pt (GeV);Tracks / bin", 300, +0.0, +300.0);

  // Lead track Eta
  histoTools_.BookHisto_1D(h_leadTrks_Eta, "leadTrks_Eta", ";#eta;Tracks / bin", 600, -3.0, +3.0);

  // Lead track Phi
  histoTools_.BookHisto_1D(h_leadTrks_Phi, "leadTrks_Phi", ";#phi (rads);Tracks / bin", 23,  -3.15,  +3.15);
  
  // Lead tracks in Eta-Phi plane
  // (Syntax: BookHisto_2D(histogram, hName, hTitle, binsX, xMin, xMax, binsY, yMin, yMax)
  histoTools_.BookHisto_2D(h_leadTrks_Phi_Eta, "leadTrks_Phi_Eta",  ";#phi (rads);#eta", 230,  -3.15,  +3.15, 600,  -3.0,  +3.0);  

  // Counters for the selection of lead tracks
  histoTools_.BookHisto_1D(h_leadTrkSelection, "leadTrkSelection", "", 5,  0.0,  5.0);  

  // Clustered tracks Pt
  histoTools_.BookHisto_1D(h_clustTrks_Pt, "clustTrks_Pt", ";Pt (GeV);Tracks / bin", 300, +0.0, +300.0);

  // Clustered tracks Eta
  histoTools_.BookHisto_1D(h_clustTrks_Eta, "clustTrks_Eta", ";#eta;Tracks / bin", 600, -3.0, +3.0);

  // Clustered tracks Phi
  histoTools_.BookHisto_1D(h_clustTrks_Phi, "clustTrks_Phi", ";#phi (rads);Tracks / bin", 23,  -3.15,  +3.15);
  
  // Clustered tracks in Eta-Phi plane
  histoTools_.BookHisto_2D(h_clustTrks_Phi_Eta, "clustTrks_Phi_Eta",  ";#phi (rads);#eta", 230,  -3.15,  +3.15, 600,  -3.0,  +3.0);  			      

  // Number of clustered tracks (per cluster)
  histoTools_.BookHisto_1D(h_trkClusters_MultiplicityPerCluster, "trkClusters_MultiplicityPerCluster", ";Number of tracks in cluster;Clusters / bin", 30, -0.5, +14.5);

  // Pt of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_Pt, "trkClusters_Pt", ";Pt (GeV);Clusters / bin", 300, +0.0, +300.0);
  
  // Invariant mass of track clusters
  histoTools_.BookHisto_1D(h_trkClusters_M, "trkClusters_M", ";Invariant mass;Clusters / bin", 100, +0.0, +5.0);

  // Number of EGs (per cluster)
  histoTools_.BookHisto_1D(h_clustEGs_MultiplicityPerCluster, "clustEGs_MultiplicityPerCluster", ";Number of EGs in cluster;Clusters / bin", 30, -0.5, +14.5);

  // Clustered EGs Pt
  histoTools_.BookHisto_1D(h_clustEGs_Et, "clustEGs_Pt", ";Pt (GeV);Events / bin", 300, +0.0, +300.0);

  // Clustered EGs Eta
  histoTools_.BookHisto_1D(h_clustEGs_Eta, "clustEGs_Eta", ";#eta;Events / bin", 600, -3.0, +3.0);

  // Clustered EGs Phi
  histoTools_.BookHisto_1D(h_clustEGs_Phi, "clustEGs_Phi", ";#phi (rads);Events / bin", 23,  -3.15,  +3.15);

  // Clustered EGs in Eta-Phi plane
  histoTools_.BookHisto_2D(h_clustEGs_Phi_Eta, "clustEGs_Phi_Eta",  ";#phi (rads);#eta", 230,  -3.15,  +3.15, 600,  -3.0,  +3.0);		      

  // Pt of EG clusters
  histoTools_.BookHisto_1D(h_EGClusters_Pt, "EGClusters_Pt", ";Pt (GeV);Clusters / bin", 300, +0.0, +300.0);
  
  // Invariant mass of EG clusters
  histoTools_.BookHisto_1D(h_EGClusters_M, "EGClusters_M", ";Invariant mass;Clusters / bin", 100, +0.0, +5.0);

  // Track-based relative isolation of tau candidates
  histoTools_.BookHisto_1D(h_trkClusters_relIso, "trkClusters_relIso", ";Relative isolation;Clusters / bin", 100, 0.0, +5.0);

  // DeltaR in MC matching of the lead track of tau candidates
  histoTools_.BookHisto_1D(hTkEG_DeltaRmatch  , "TkEG_DeltaRmatch"  , ";#Delta R;Particles / bin",   100,  0.0,   +1.0);

  // Visible Et of the matched generator taus
  histoTools_.BookHisto_1D(hTkEG_VisEt, "TkEG_VisEt", ";VisEt (GeV);Particles / bin", 40, 0.0,  +200.0);

  // Visible Et of all generator taus
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt, "hMcHadronicTau_VisEt", ";VisEt (GeV);Particles / bin", 40, 0.0,  +200.0);

  // Turn-on histograms
  histoTools_.BookHisto_1D(hTurnOn25_all, "TurnOn25_all", "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTurnOn25_relIso, "TurnOn25_relIso", "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTurnOn50_all, "TurnOn50_all", "", 40, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hTurnOn50_relIso, "TurnOn50_relIso", "", 40, 0.0,  +200.0);
  
  // Single-tau rates
  histoTools_.BookHisto_1D(hRateSingleTau    , "Rate_SingleTau"    , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateSingleTau_C  , "Rate_SingleTau_C"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateSingleTau_I  , "Rate_SingleTau_I"  , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hRateSingleTau_F  , "Rate_SingleTau_F"  , "", 200, 0.0,  +200.0);  
  
  // Single-tau efficiencies
  histoTools_.BookHisto_1D(hEffSingleTau     , "Eff_SingleTau"     , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffSingleTau_C   , "Eff_SingleTau_C"   , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffSingleTau_I   , "Eff_SingleTau_I"   , "", 200, 0.0,  +200.0);
  histoTools_.BookHisto_1D(hEffSingleTau_F   , "Eff_SingleTau_F"   , "", 200, 0.0,  +200.0);
  
  return;
}



//============================================================================
vector<L1TkEGParticle> TkCalo::GetMcMatchedL1TkEGs(vector<L1TkEGParticle> L1TkEGs)
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
double TkCalo::GetMatchingGenParticle(TTTrack track, GenParticle *genParticlePtr)					    
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
void TkCalo::FillSingleTau_(vector<L1TkEGParticle> L1TkEGs, 
			    TH1D *hRate,
			    TH1D *hEfficiency,
			    double minEta,
			    double maxEta)
//============================================================================
{

  // Sanity check
  if( L1TkEGs.size() < 1 ) return;
  
  // Fill rate
  double ldgEt = L1TkEGs.at(0).GetTrackBasedEt();

  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) < minEta) return;
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) > maxEta) return;
    
  FillRate_(hRate, ldgEt);
  
  // Get MC-matched trigger objects
  vector<L1TkEGParticle> L1TkEGs_mcMatched = GetMcMatchedL1TkEGs(TauCandidates);
  if (L1TkEGs_mcMatched.size() < 1) return;
  
  // Check that all taus were found
  //if(!bFoundAllTaus_) return; //TODO what is the point?

  // Fill efficiency
  double ldgEt_mcMatched = L1TkEGs_mcMatched.at(0).GetTrackBasedEt();
  FillEfficiency_(hEfficiency, ldgEt_mcMatched);

  return;
}


//============================================================================
void TkCalo::FillDiTau_(vector<L1TkEGParticle> L1TkEGs, 
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
  double subLdgEt = L1TkEG.GetTrackBasedEt();
  
  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) < minEta) return;
  if ( abs(L1TkEGs.at(0).GetLeadingTrack().getEta()) > maxEta) return;
  if ( abs(L1TkEGs.at(1).GetLeadingTrack().getEta()) < minEta) return;
  if ( abs(L1TkEGs.at(1).GetLeadingTrack().getEta()) > maxEta) return;

  FillRate_(hRate, subLdgEt);

  // Get MC-Matched trigger objects
  vector<L1TkEGParticle> L1TkEGs_mcMatched = GetMcMatchedL1TkEGs(L1TkEGs);
  if (L1TkEGs_mcMatched.size() < 2) return;
    
  // Check that all taus were found
 // if(!bFoundAllTaus_) return; //TODO

  // Fill efficiency
  double subLdgEt_mcMatched = L1TkEGs_mcMatched.at(1).GetTrackBasedEt();
  FillEfficiency_(hEfficiency, subLdgEt_mcMatched);

  return;
}


//============================================================================
void TkCalo::FillDiTau_(vector<L1TkEGParticle> L1TkEGs1,
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
  double ldgEt1 = L1TkEGs1.at(0).GetTrackBasedEt();
  double ldgEt2 = L1TkEGs2.at(0).GetTrackBasedEt();

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
  double ldgEt1_mcMatched = L1TkEGs1_mcMatched.at(0).GetTrackBasedEt();
  double ldgEt2_mcMatched = L1TkEGs2_mcMatched.at(0).GetTrackBasedEt();

  
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
void TkCalo::FillRate_(TH1D *hRate,
			   const double ldgEt)
//============================================================================
{
  
  if (ldgEt < 0) return;
  hRate ->Fill( ldgEt );
  return;
}


//============================================================================
void TkCalo::FillRate_(TH2D *hRate,
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
void TkCalo::FillEfficiency_(TH1D *hEfficiency,
			     const double ldgEt)
//============================================================================
{
  
  histoTools_.FillAllBinsUpToValue_1D(hEfficiency, ldgEt);

  return;
}


//============================================================================
void TkCalo::FillTurnOn_Numerator_(vector<L1TkEGParticle> L1TkEGs, 
				   const double minEt,
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
      // Skip if trigger object has eT < minEt
      if (tau->GetTrackBasedEt() < minEt) continue;
      passed_cut++;      
            
      // Fill the turn-on
      hTurnOn->Fill( tau->GetGenTauEt() );
      // DEBUG

    } // For-loop: L1TkEGs
//      cout << "Filling turn-on histogram: iterated particles=" << iterated_particles << ", matched=" << matched_particles << ", passed_cut=" << passed_cut << endl;

  return;
   
}

//============================================================================
void TkCalo::FinaliseEffHisto_(TH1D *histo, 
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
void TkCalo::FinaliseEffHisto_(TH2D *histo, 
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

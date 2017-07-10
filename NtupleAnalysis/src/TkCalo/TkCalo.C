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
  
  cfg_DEBUG = true;
  if (cfg_DEBUG) std::cout << "=== TkCalo::InitVars_()" << std::endl;
  
  cfg_AddL1Tks   = true;
  cfg_AddEGs     = true;
  
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
    if (true) { //TODO if (cfg_AddL1Tks) {
	  if (cfg_DEBUG) cout << "\n=== TTracks (" << L1Tks_Pt->size() << ")" << endl;
	  TTTracks = GetTTTracks(cfg_tk_minPt, cfg_tk_minEta, cfg_tk_maxEta, cfg_tk_maxChiSqRed, cfg_tk_minStubs, cfg_tk_minStubsPS, cfg_tk_maxStubsPS, cfg_tk_nFitParams, false);
	  sort( TTTracks.rbegin(), TTTracks.rend() ); // Sort from highest Pt to lowest (not done by default)
	}    
	
	// EGs
	if (cfg_AddEGs) {
	  L1EGs = GetL1EGs(cfg_DEBUG);
	  sort( L1EGs.rbegin(), L1EGs.rend() ); // Sort from highest Et to lowest Et (should be already done by default)
	  if (cfg_DEBUG) cout << "\n=== L1EGs (" << L1EGs.size() << ")" << endl;
	}

    ////////////////////////////////////////////////
    // TkCalo algorithm
    ////////////////////////////////////////////////    

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

    // Build EG clusters
    vector<L1EG> EGcluster;
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

        // Build a tau candidate from tracks and EGs
        L1TkEGParticle newTauCandidate(trackTauCandidates[i], EGcluster);
        cout << "Constructed a new tau candidate with track invariant mass of " << newTauCandidate.GetTrackInvMass() 
             << " and EG invariant mass of " << newTauCandidate.GetEGInvMass() << endl;
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

   
    ////////////////////////////////////////////////
    // Fill Turn-On histograms
    ////////////////////////////////////////////////
    
    // TODO
    
    // Progress bar
    if (!cfg_DEBUG) auxTools_.ProgressBar(jentry, nEntries, 100, 100);
    
  }// For-loop: Entries

  ////////////////////////////////////////////////
  // Fill counters
  ////////////////////////////////////////////////

  // TODO
  
  
  ////////////////////////////////////////////////
  // Convert/Finalise Histos
  ////////////////////////////////////////////////

  // TODO


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

  // Number of lead tracks (per event)
  histoTools_.BookHisto_1D(h_leadTrks_Multiplicity, "leadTrks_Multiplicity", ";Number of lead tracks in event;Events / bin", 30, -0.5, +14.5);

  // Lead track Pt
  histoTools_.BookHisto_1D(h_leadTrks_Pt, "leadTrks_Pt", ";Pt (GeV);Events / bin", 300, +0.0, +300.0);

  // Lead track Eta
  histoTools_.BookHisto_1D(h_leadTrks_Eta, "leadTrks_Eta", ";#eta;Events / bin", 600, -3.0, +3.0);

  // Lead track Phi
  histoTools_.BookHisto_1D(h_leadTrks_Phi, "leadTrks_Phi", ";#phi (rads);Events / bin", 23,  -3.15,  +3.15);
  
  // Lead tracks in Eta-Phi plane
  // (Syntax: BookHisto_2D(histogram, hName, hTitle, binsX, xMin, xMax, binsY, yMin, yMax)
  histoTools_.BookHisto_2D(h_leadTrks_Phi_Eta, "leadTrks_Phi_Eta",  ";#phi (rads);#eta", 230,  -3.15,  +3.15, 600,  -3.0,  +3.0);  

  // Counters for the selection of lead tracks
  histoTools_.BookHisto_1D(h_leadTrkSelection, "leadTrkSelection", "", 5,  0.0,  5.0);  

  // Clustered tracks Pt
  histoTools_.BookHisto_1D(h_clustTrks_Pt, "clustTrks_Pt", ";Pt (GeV);Events / bin", 300, +0.0, +300.0);

  // Clustered tracks Eta
  histoTools_.BookHisto_1D(h_clustTrks_Eta, "clustTrks_Eta", ";#eta;Events / bin", 600, -3.0, +3.0);

  // Clustered tracks Phi
  histoTools_.BookHisto_1D(h_clustTrks_Phi, "clustTrks_Phi", ";#phi (rads);Events / bin", 23,  -3.15,  +3.15);
  
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

  return;
}


#endif

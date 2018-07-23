#ifndef TauTrigger_cxx
#define TauTrigger_cxx

#include "TauTrigger.h"
#include "../utilities/constants.h"

//****************************************************************************
void TauTrigger::InitObjects(void)
//****************************************************************************
{

  pvProducer = new L1TkPrimaryVertex(this->s);

  return; 
}


//****************************************************************************
void TauTrigger::InitVars_()
//****************************************************************************
{
  
  // MC matching
  DEBUG                    = true;
  mcMatching_maxDeltaR     = +0.10;
  mcMatching_unique        = true;
  tauDecayMode             = "all";

  // L1TkTau 
  tk_CollectionType        = "TTPixelTracks"; // "TTTracks"; // "TTPixelTracks";
  tauAlgorithmMode         = "ShrinkingCone";
  matchTk_minPt_           = +5.00;

  // Dataset-related
  datasets_  = datasets_.GetDataset(mcSample);
  realTauMom = datasets_.McTauMomPdgId_;
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;

  // Signal cone 
  sigCone_tkCollectionWP   = "vLooseWP";
  matchTk_caloDeltaR_      = +0.12;
  sigCone_ShrinkConstant   = +1.50; // (+0.0 = signal cone, not signal annulus)
  sigCone_maxDeltaR        = +0.15;
  sigCone_maxTkInvMass     = +1.77;  // 3-pr
  sigCone_maxTkDeltaPOCAz  = +0.20;  // 3-pr
  cout << "*** WARNING! Why sigCone_ShrinkConstant = +1.50 ??????" << endl;

  // Isolation cone
  isoCone_tkCollectionWP   = "vLooseWP"; 
  isoCone_ShrinkConstant   = +3.50; // TP: 3.50 GeV
  isoCone_VtxIsoWP         = +0.40; // TP: 1.0cm
  isoCone_minDeltaR        = sigCone_maxDeltaR;
  isoCone_maxDeltaR        = +0.40;

  // Di-Tau
  diTau_deltaPOCAz         = +0.50;

  PrintSettings();

  return;
}


//****************************************************************************
void TauTrigger::PrintSettings(void)
//****************************************************************************
{

  // Inform user of settings
  Table settings("Variable | Value | Units", "Text");  // Table settingsTable("Variable & Value & Units", "LaTeX", "l l l");
  settings.AddRowColumn(0, "MC Sample");
  settings.AddRowColumn(0, mcSample );
  settings.AddRowColumn(1, "Tau Algorithm");
  settings.AddRowColumn(1, tauAlgorithmMode);
  settings.AddRowColumn(2, "Track Collection");
  settings.AddRowColumn(2, tk_CollectionType);
  settings.AddRowColumn(3, "SigTks WP");
  settings.AddRowColumn(3, sigCone_tkCollectionWP);
  settings.AddRowColumn(4, "Calo-Matching Tk Pt (min)");
  settings.AddRowColumn(4, auxTools_.ToString(matchTk_minPt_) );
  settings.AddRowColumn(4, "GeVc^{-1}");
  settings.AddRowColumn(5, "Calo-Matching Tk DeltaR (max)");
  settings.AddRowColumn(5, auxTools_.ToString(matchTk_caloDeltaR_) );
  settings.AddRowColumn(6, "SigCone Shrink Constant");
  settings.AddRowColumn(6, auxTools_.ToString(sigCone_ShrinkConstant) );
  settings.AddRowColumn(6, "GeV");
  settings.AddRowColumn(7, "SigCone DeltaR (max)");
  settings.AddRowColumn(7, auxTools_.ToString(sigCone_maxDeltaR) );
  settings.AddRowColumn(8, "SigCone InvMass (max)");
  settings.AddRowColumn(8, auxTools_.ToString(sigCone_maxTkInvMass) );
  settings.AddRowColumn(8, "GeVc^{-2}");
  settings.AddRowColumn(9, "SigCone TkDeltaPOCAz (max)");
  settings.AddRowColumn(9, auxTools_.ToString(sigCone_maxTkDeltaPOCAz) );
  settings.AddRowColumn(9, "cm");
  settings.AddRowColumn(10, "IsoTks WP");
  settings.AddRowColumn(10, isoCone_tkCollectionWP);
  settings.AddRowColumn(11, "IsoCone Shrink Constant");
  settings.AddRowColumn(11, auxTools_.ToString(isoCone_ShrinkConstant) );
  settings.AddRowColumn(11, "GeV");
  settings.AddRowColumn(12, "IsoCone DeltaR (min)");
  settings.AddRowColumn(12, auxTools_.ToString(isoCone_minDeltaR) );
  settings.AddRowColumn(13, "IsoCone DeltaR (max)");
  settings.AddRowColumn(13, auxTools_.ToString(isoCone_maxDeltaR) );
  settings.AddRowColumn(14, "VtxIso WP" );
  settings.AddRowColumn(14, auxTools_.ToString(isoCone_VtxIsoWP) );
  settings.AddRowColumn(14, "cm");
  settings.AddRowColumn(15, "MC-Matching DeltaR (max)");
  settings.AddRowColumn(15, auxTools_.ToString(mcMatching_maxDeltaR) );
  settings.AddRowColumn(16, "Unique MC-Matching");
  settings.AddRowColumn(16, auxTools_.ToString(mcMatching_unique) );
  settings.AddRowColumn(17, "MC-Tau Decay Mode (all, 1p, 3p) ");
  settings.AddRowColumn(17, tauDecayMode);
  settings.AddRowColumn(18, "MC-Tau Mom (PdgId)");
  settings.AddRowColumn(18, auxTools_.ToString(realTauMom));
  settings.AddRowColumn(19, "Num of MC-Taus");
  settings.AddRowColumn(19, auxTools_.ToString(nMaxNumOfHTausPossible));
  settings.AddRowColumn(20, "Di-TAu Delta POCA-z");
  settings.AddRowColumn(20, auxTools_.ToString(diTau_deltaPOCAz));
  settings.AddRowColumn(20, "cm");  

  settings.Print();
  
  return;
}


//****************************************************************************
void TauTrigger::Loop()
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


  // For-loop: Entries
  for (int jentry=0; jentry < nEntries; jentry++){
  
    if(DEBUG) cout << "Entry = " << jentry << endl;
    
    // Init variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // MC Tau containers
    vector<int> mcHadronicTaus         = GetTriggerHadronicTaus_(realTauMom, 0.0);
    int nMcHadronicTaus_VisEt20Eta2p3  = (int) GetTriggerHadronicTaus_(realTauMom, 20.0).size();
    bFoundAllTaus_ = ( nMcHadronicTaus_VisEt20Eta2p3 >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;

    // 1-prong or 3-prong or all?
    bool bForbiddenTauDecayMode = EventWithForbiddenTauDecayMode(tauDecayMode, mcHadronicTaus);
    if(bForbiddenTauDecayMode) continue;
    
    // Sanity checks
    auxTools_.EnsureVectorIsSorted(*L1CaloTau_Et, true);
    auxTools_.EnsureVectorIsSorted(*L1Tks_Pt, true);
    auxTools_.EnsureVectorIsSorted(*L1PixTks_Pt, true);

    // Get PV z-position
    vector<int> pvTks_index;
    pv_z = 0.0;
    pv_z = pvProducer->GetPrimaryVertexZ("default", pvTks_index);

    // Event-Type Histograms
    hHepMCEvt_VtxX_VtxY     ->Fill(HepMCEvt_VtxX, HepMCEvt_VtxY);
    hHepMCEvt_VtxZ          ->Fill(HepMCEvt_VtxZ);
    hL1TkPV_VtxZ            ->Fill(pv_z);
    hPrimaryVertex_DeltaZ   ->Fill(pv_z - HepMCEvt_VtxZ);
    hPrimaryVertex_AbsDeltaZ->Fill( abs(pv_z - HepMCEvt_VtxZ) );


    // Initialisations
    vector<L1TkTauParticle> L1TkTaus_Calo;
    vector<L1TkTauParticle> L1TkTaus_Tk;
    vector<L1TkTauParticle> L1TkTaus_VtxIso;

    // For-loop: L1CaloTaus
    for (Size_t iCalo = 0; iCalo < L1CaloTau_Et->size(); iCalo++) {            

      // Variable declaration
      int matchTk_index     = -1;
      int matchGenp_Index   = -1;
      double sigCone_minDeltaR_tmp =  0.0;
      double sigCone_maxDeltaR_tmp =  0.0;
      double isoCone_minDeltaR_tmp =  0.0;
      double matchTk_deltaR_tmp    = -1.0;
      double matchGenp_deltaR;
      vector<int> sigTks_index;
      vector<int> isoTks_index;
      double relIso;
      double vtxIso;

      // Get L1TkTauParticles
      if ( tk_CollectionType.compare("TTPixelTracks") == 0 )
	{
	  GetL1TkTauMatchPixTrack(matchTk_index, iCalo, sigCone_tkCollectionWP, matchTk_minPt_, matchTk_caloDeltaR_, matchTk_deltaR_tmp);
	  GetL1TkTauAlgoConeSizes(iCalo, tauAlgorithmMode, sigCone_maxDeltaR, sigCone_minDeltaR_tmp, sigCone_maxDeltaR_tmp, isoCone_minDeltaR_tmp);

	  GetL1TkTauSigConePixTracks(matchTk_index, sigCone_minDeltaR_tmp, sigCone_maxDeltaR_tmp, sigCone_tkCollectionWP, sigTks_index);
	  GetL1TkTauIsoConePixTracks(matchTk_index, isoCone_minDeltaR_tmp, isoCone_maxDeltaR    , isoCone_tkCollectionWP, isoTks_index);

	  GetL1TkTauPixIsolation(matchTk_index, isoTks_index, relIso, vtxIso);

	}      
      else if ( tk_CollectionType.compare("TTTracks") == 0 )
	{
	  GetL1TkTauMatchTrack(matchTk_index, iCalo, sigCone_tkCollectionWP, matchTk_minPt_, matchTk_caloDeltaR_, matchTk_deltaR_tmp);
	  GetL1TkTauAlgoConeSizes(iCalo, tauAlgorithmMode, sigCone_maxDeltaR, sigCone_minDeltaR_tmp, sigCone_maxDeltaR_tmp, isoCone_minDeltaR_tmp);
	  
	  GetL1TkTauSigConeTracks(matchTk_index, sigCone_minDeltaR_tmp, sigCone_maxDeltaR_tmp, sigCone_tkCollectionWP, sigTks_index);
	  GetL1TkTauIsoConeTracks(matchTk_index, isoCone_minDeltaR_tmp, isoCone_maxDeltaR    , isoCone_tkCollectionWP, isoTks_index);
	  
	  GetL1TkTauIsolation(matchTk_index, isoTks_index, relIso, vtxIso);
	}
      else
	{
	  cout << "E R R O R ! TauTrigger::Loop(...) - Invalid Track Collection Type \"" << tk_CollectionType << "\". EXIT" << endl;
	  exit(1);
	}

    
      // Construct L1TkTaus and store them
      L1TkTauParticle L1TkTau(iCalo, matchTk_index, matchTk_deltaR_tmp, sigTks_index, isoTks_index, vtxIso, relIso, 
			      sigCone_minDeltaR_tmp, sigCone_maxDeltaR_tmp, isoCone_minDeltaR_tmp, isoCone_maxDeltaR);
      GetL1CaloTauMatchGenp(iCalo, mcHadronicTaus, matchGenp_Index, matchGenp_deltaR);
      L1TkTau.SetMatchGenp( matchGenp_Index, matchGenp_deltaR);
      L1TkTaus_Calo.push_back(L1TkTau);
      PrintL1TkTauProperties(DEBUG, L1TkTau);

    } // For-loop: L1CaloTaus


    // Associate MC Taus ONLY to the closest L1CaloTau match (unique matching). Makes difference for large values of mcMatching_maxDeltaR (e.g. 0.8)
    if (mcMatching_unique) GetL1CaloTauUniqueMatchGenp(L1TkTaus_Calo);


    ////////////////////////////////////////////////
    /// Construct additional collections of L1TkTaus 
    ////////////////////////////////////////////////
    for( int i=0; i < (int) L1TkTaus_Calo.size(); i++){

      L1TkTauParticle L1TkTau = L1TkTaus_Calo.at(i);

      // Apply track-matching
      if ( L1TkTau.matchTk_Index_ < 0) continue;
      L1TkTaus_Tk.push_back(L1TkTau);

      // Apply isolation
      if ( abs(L1TkTau.vtxIso_) < isoCone_VtxIsoWP) continue;
      L1TkTaus_VtxIso.push_back(L1TkTau);

    }
      

    ////////////////////////////////////////////////
    /// L1Tktau Properties 
    ////////////////////////////////////////////////
    vector<L1TkTauParticle> myL1TkTaus = L1TkTaus_VtxIso; // L1TkTaus_VtxIso; // L1TkTaus_Tk; 
    hL1TkTau_Multiplicity ->Fill( myL1TkTaus.size() );

    /// For-loop: L1TkTaus
    for( int i=0; i < (int) myL1TkTaus.size(); i++){
           
      // Variables
      L1TkTauParticle L1TkTau  = myL1TkTaus.at(i);
      TLorentzVector sigTks_p4 = GetL1TkTauSigTksP4(tk_CollectionType, L1TkTau);
      TLorentzVector isoTks_p4 = GetL1TkTauIsoTksP4(tk_CollectionType, L1TkTau);
    
      // Do not skip if using MinBias sample as no real taus exist!
      if ( (L1TkTau.matchGenp_Index_ < 0) && (mcSample.compare("MinBias") !=0) ) continue;
  
      // Fill Histos
      hL1TkTau_Rtau        ->Fill( GetL1TkTauLdgTkPt(tk_CollectionType, L1TkTau)  / L1CaloTau_Et->at(L1TkTau.caloTau_Index_) );
      hL1TkTau_CHF         ->Fill( L1CaloTau_Et->at(L1TkTau.caloTau_Index_) / sigTks_p4.Et() );
      hL1TkTau_NHF         ->Fill( (L1CaloTau_Et->at(L1TkTau.caloTau_Index_) - sigTks_p4.Et())/L1CaloTau_Et->at(L1TkTau.caloTau_Index_) );
      hL1TkTau_NHFAbs      ->Fill( abs(L1CaloTau_Et->at(L1TkTau.caloTau_Index_) - sigTks_p4.Et())/L1CaloTau_Et->at(L1TkTau.caloTau_Index_) );
      hL1TkTau_NSigTks     ->Fill( L1TkTau.sigTks_Index_.size() );
      hL1TkTau_NIsoTks     ->Fill( L1TkTau.isoTks_Index_.size() );
      hL1TkTau_InvMass     ->Fill( sigTks_p4.M()                ); 
      hL1TkTau_InvMassIncl ->Fill( (sigTks_p4 + isoTks_p4).M()  );
      hL1TkTau_SigConeRMin ->Fill( L1TkTau.sigCone_minDeltaR_   );
      hL1TkTau_SigConeRMax ->Fill( L1TkTau.sigCone_maxDeltaR_   );
      hL1TkTau_IsoConeRMin ->Fill( L1TkTau.isoCone_minDeltaR_   );
      hL1TkTau_IsoConeRMax ->Fill( L1TkTau.isoCone_maxDeltaR_   );
      hL1TkTau_DeltaRGenP  ->Fill( L1TkTau.matchGenp_deltaR_    );           
      // if (L1TkTau.isoTks_Index_.size() > 0) {
      hL1TkTau_RelIso      ->Fill( L1TkTau.relIso_            );
      hL1TkTau_VtxIso      ->Fill( L1TkTau.vtxIso_            );
      hL1TkTau_VtxIsoAbs   ->Fill( abs(L1TkTau.vtxIso_)       );
	//}
    

      // Matching Tk
      const int mTk     =  L1TkTau.matchTk_Index_;
      const int mCalo   =  L1TkTau.caloTau_Index_;
      TLorentzVector mCalo_p4 = auxTools_.GetTLorentzVector( L1CaloTau_Et->at(mCalo), L1CaloTau_Eta->at(mCalo), L1CaloTau_Phi->at(mCalo), L1CaloTau_E->at(mCalo) );
      if ( tk_CollectionType.compare("TTPixelTracks") == 0 ){

	// Get the transverse component of this track with respect to the matching track
	TVector3 mPixTk_p3      = auxTools_.GetTVector3( L1PixTks_Px->at(mTk), L1PixTks_Py->at(mTk), L1PixTks_Pz->at(mTk) );
	double matchPixTk_PtRel = mPixTk_p3.Perp(mCalo_p4.Vect());

	hL1TkTau_Charge                   ->Fill( s->GetChargeFromPixTracks(L1TkTau.sigTks_Index_) );
	hL1TkTau_MatchPixTk_PtRel         ->Fill( matchPixTk_PtRel              );
	hL1TkTau_MatchPixTk_DeltaR        ->Fill( L1TkTau.matchTk_deltaR_       );
	hL1TkTau_MatchPixTk_StubPtCons    ->Fill( L1Tks_StubPtConsistency->at( L1PixTks_TTTrackIndex->at(mTk) ) );
	hL1TkTau_MatchPixTk_NStubs        ->Fill( s->GetNumOfStubs  ( L1PixTks_TTTrackIndex->at(mTk) ) );
	hL1TkTau_MatchPixTk_NPsStubs      ->Fill( s->GetNumOfPSStubs( L1PixTks_TTTrackIndex->at(mTk) ) );
	hL1TkTau_MatchPixTk_NBarrelStubs  ->Fill( s->GetNumOfBarrelStubs( L1PixTks_TTTrackIndex->at(mTk) ) );
	hL1TkTau_MatchPixTk_NEndcapStubs  ->Fill( s->GetNumOfEndcapStubs( L1PixTks_TTTrackIndex->at(mTk) ) );
	hL1TkTau_MatchPixTk_Pt            ->Fill( L1PixTks_Pt->at(mTk)          );
	hL1TkTau_MatchPixTk_NPixHits      ->Fill( L1PixTks_NPixHits->at(mTk)    );
	hL1TkTau_MatchPixTk_Eta           ->Fill( L1PixTks_Eta->at(mTk)         );
	hL1TkTau_MatchPixTk_POCAz         ->Fill( L1PixTks_POCAz->at(mTk)       );
	hL1TkTau_MatchPixTk_POCAzSig      ->Fill( s->GetPixTrackPOCAzSig( mTk ) );
	hL1TkTau_MatchPixTk_d0            ->Fill( s->GetPixTrackD0( mTk )       );
	hL1TkTau_MatchPixTk_d0Abs         ->Fill( abs(s->GetPixTrackD0( mTk ))  );
	hL1TkTau_MatchPixTk_d0Sig         ->Fill( s->GetPixTrackD0Sig( mTk )    );
	hL1TkTau_MatchPixTk_d0SigAbs      ->Fill( abs(s->GetPixTrackD0Sig(mTk)) );
	hL1TkTau_MatchPixTk_ChiSquared    ->Fill( L1PixTks_ChiSquared->at(mTk)  );
	hL1TkTau_MatchPixTk_RedChiSquared ->Fill( s->GetPixTrackRedChiSq(mTk)   ); // L1PixTks_RedChiSquared->at(mTk)
	hL1TkTau_MatchPixTk_SigmaRInv     ->Fill( L1PixTks_SigmaRInv->at(mTk)   );
	hL1TkTau_MatchPixTk_SigmaPhi0     ->Fill( L1PixTks_SigmaPhi0->at(mTk)   );
	hL1TkTau_MatchPixTk_SigmaD0       ->Fill( L1PixTks_SigmaD0->at(mTk)     );
	hL1TkTau_MatchPixTk_SigmaT        ->Fill( L1PixTks_SigmaT->at(mTk)      );
	hL1TkTau_MatchPixTk_SigmaZ0       ->Fill( L1PixTks_SigmaZ0->at(mTk)     );

	// For-loop: SigCone Tks
	for (int j = 0; j < (int) L1TkTau.sigTks_Index_.size(); j++ ) {

	  // Get indices and skip if matching track
	  int sigTk = L1TkTau.sigTks_Index_.at(j);
	  if (sigTk == mTk) continue;

	  // Get the transverse component of this track with respect to the matching track
	  TVector3 sigTk_p3  = auxTools_.GetTVector3( L1PixTks_Px->at(sigTk), L1PixTks_Py->at(sigTk), L1PixTks_Pz->at(sigTk) );
	  double sigTk_PtRel = sigTk_p3.Perp(mPixTk_p3); // mPixTk_p3.Perp(sigTk_p3);

	  // Fill Histograms
	  hL1TkTau_SigTks_Pt        ->Fill( L1PixTks_Pt->at(sigTk)          );
	  hL1TkTau_SigTks_Eta       ->Fill( L1PixTks_Eta->at(sigTk)         );
	  hL1TkTau_SigTks_POCAz     ->Fill( L1PixTks_POCAz->at(sigTk)       );
	  hL1TkTau_SigTks_DeltaPOCAz->Fill( abs( L1PixTks_POCAz->at(mTk) - L1PixTks_POCAz->at(sigTk) ) );
	  hL1TkTau_SigTks_d0        ->Fill( s->GetPixTrackD0(sigTk)         );
	  hL1TkTau_SigTks_d0Abs     ->Fill( abs(s->GetPixTrackD0(sigTk))    );
	  hL1TkTau_SigTks_d0Sig     ->Fill( s->GetPixTrackD0Sig(sigTk)      );
	  hL1TkTau_SigTks_d0SigAbs  ->Fill( abs(s->GetPixTrackD0Sig(sigTk)) );
	  hL1TkTau_SigTks_PtRel     ->Fill( sigTk_PtRel                     );
	  hL1TkTau_SigTks_StubPtCons->Fill( L1Tks_StubPtConsistency->at( L1PixTks_TTTrackIndex->at(sigTk) ) );

	}	  

	// For-loop: IsoCone Tracks
	for (int k = 0; k < (int) L1TkTau.isoTks_Index_.size(); k++ ) {
	  
	  // Get indices and skip if matching track
	  int isoTk = L1TkTau.isoTks_Index_.at(k);
	  
	  // Get the transverse component of this track with respect to the matching track
	  TVector3 isoTk_p3  = auxTools_.GetTVector3( L1PixTks_Px->at(isoTk), L1PixTks_Py->at(isoTk), L1PixTks_Pz->at(isoTk) );
	  double isoTk_PtRel = isoTk_p3.Perp(mPixTk_p3); // mPixTk_p3.Perp(isoTk_p3);
	  
	  // Fill Histograms
	  hL1TkTau_IsoTks_Pt        ->Fill( L1PixTks_Pt->at(isoTk)          );
	  hL1TkTau_IsoTks_Eta       ->Fill( L1PixTks_Eta->at(isoTk)         );
	  hL1TkTau_IsoTks_POCAz     ->Fill( L1PixTks_POCAz->at(isoTk)       );
	  hL1TkTau_IsoTks_DeltaPOCAz->Fill( abs( L1PixTks_POCAz->at(mTk) - L1PixTks_POCAz->at(isoTk) ) );
	  hL1TkTau_IsoTks_d0        ->Fill( s->GetPixTrackD0(isoTk)         );
	  hL1TkTau_IsoTks_d0Abs     ->Fill( abs(s->GetPixTrackD0(isoTk))    );
	  hL1TkTau_IsoTks_d0Sig     ->Fill( s->GetPixTrackD0Sig(isoTk)      );
	  hL1TkTau_IsoTks_d0SigAbs  ->Fill( abs(s->GetPixTrackD0Sig(isoTk)) );
	  hL1TkTau_IsoTks_PtRel     ->Fill( isoTk_PtRel                     );
	  hL1TkTau_IsoTks_StubPtCons->Fill( L1Tks_StubPtConsistency->at( L1PixTks_TTTrackIndex->at(isoTk) ) );
	  
	}
      }      
      else if ( tk_CollectionType.compare("TTTracks") == 0 )
	{

	  // Get the transverse component of this track with respect to the matching track
	  TVector3 mTk_p3      = auxTools_.GetTVector3( L1Tks_Px->at(mTk), L1Tks_Py->at(mTk), L1Tks_Pz->at(mTk) );
	  double matchTk_PtRel = mCalo_p4.Vect().Perp(mTk_p3); // mTk_p3.Perp(mCalo_p4.Vect());
	  
	  hL1TkTau_Charge                ->Fill( s->GetChargeFromTracks(L1TkTau.sigTks_Index_) );
	  hL1TkTau_MatchTk_PtRel         ->Fill( matchTk_PtRel                    );
	  hL1TkTau_MatchTk_DeltaR        ->Fill( L1TkTau.matchTk_deltaR_          );
	  hL1TkTau_MatchTk_Pt            ->Fill( L1Tks_Pt->at(mTk)                );
	  hL1TkTau_MatchTk_Eta           ->Fill( L1Tks_Eta->at(mTk)               );
	  hL1TkTau_MatchTk_POCAz         ->Fill( L1Tks_POCAz->at(mTk)             );
	  hL1TkTau_MatchTk_d0            ->Fill( s->GetTrackD0(mTk)               );
	  hL1TkTau_MatchTk_d0Abs         ->Fill( abs(s->GetTrackD0(mTk))          );
	  hL1TkTau_MatchTk_NStubs        ->Fill( s->GetNumOfStubs(mTk)            );
	  hL1TkTau_MatchTk_NPsStubs      ->Fill( s->GetNumOfPSStubs(mTk)          );
	  hL1TkTau_MatchTk_NBarrelStubs  ->Fill( s->GetNumOfBarrelStubs(mTk)      );
	  hL1TkTau_MatchTk_NEndcapStubs  ->Fill( s->GetNumOfEndcapStubs(mTk)      );
	  hL1TkTau_MatchTk_StubPtCons    ->Fill( L1Tks_StubPtConsistency->at(mTk) );
	  hL1TkTau_MatchTk_ChiSquared    ->Fill( L1Tks_ChiSquared->at(mTk)        );
	  // hL1TkTau_MatchTk_RedChiSquared ->Fill( L1Tks_RedChiSquared->at(mTk)     );
	  hL1TkTau_MatchTk_IsGenuine     ->Fill( L1Tks_IsGenuine->at(mTk)         );
	  hL1TkTau_MatchTk_IsUnknown     ->Fill( L1Tks_IsUnknown->at(mTk)         );
	  hL1TkTau_MatchTk_IsCombinatoric->Fill( L1Tks_IsCombinatoric->at(mTk)    );

	  // For-loop: SigCone Tracks 
	  for (int j = 0; j < (int) L1TkTau.sigTks_Index_.size(); j++ ) {

	    // Get indices and skip if matching track
	    int sigTk = L1TkTau.sigTks_Index_.at(j);
	    if (sigTk == mTk) continue;

	    // Get the transverse component of this track with respect to the matching track
	    TVector3 sigTk_p3  = auxTools_.GetTVector3( L1Tks_Px->at(sigTk), L1Tks_Py->at(sigTk), L1Tks_Pz->at(sigTk) );
	    double sigTk_PtRel = sigTk_p3.Perp(mTk_p3); // mTk_p3.Perp(sigTk_p3);

	    // Fill Histograms
	    hL1TkTau_SigTks_Pt        ->Fill( L1Tks_Pt->at(sigTk)       );
	    hL1TkTau_SigTks_Eta       ->Fill( L1Tks_Eta->at(sigTk)      );
	    hL1TkTau_SigTks_POCAz     ->Fill( L1Tks_POCAz->at(sigTk)    );
	    hL1TkTau_SigTks_DeltaPOCAz->Fill( abs( L1Tks_POCAz->at(mTk) - L1Tks_POCAz->at(sigTk) ) );
	    hL1TkTau_SigTks_d0        ->Fill( s->GetTrackD0(sigTk)      );
	    hL1TkTau_SigTks_d0Abs     ->Fill( abs(s->GetTrackD0(sigTk)) );
	    hL1TkTau_SigTks_PtRel     ->Fill( sigTk_PtRel   );
	    hL1TkTau_SigTks_StubPtCons->Fill( L1Tks_StubPtConsistency->at(sigTk) );
	  }
	  
	  // For-loop: IsoCone Tracks 
	  for (int k = 0; k < (int) L1TkTau.isoTks_Index_.size(); k++ ) {
	    
	    // Get indices and skip if matching track
	    int isoTk = L1TkTau.isoTks_Index_.at(k);
	    
	    // Get the transverse component of this track with respect to the matching track
	    TVector3 isoTk_p3  = auxTools_.GetTVector3( L1Tks_Px->at(isoTk), L1Tks_Py->at(isoTk), L1Tks_Pz->at(isoTk) );
	    double isoTk_PtRel = isoTk_p3.Perp(mTk_p3);  // mTk_p3.Perp(isoTk_p3);
		  
	    // Fill Histograms
	    hL1TkTau_IsoTks_Pt        ->Fill( L1Tks_Pt->at(isoTk)       );
	    hL1TkTau_IsoTks_Eta       ->Fill( L1Tks_Eta->at(isoTk)      );
	    hL1TkTau_IsoTks_POCAz     ->Fill( L1Tks_POCAz->at(isoTk)    );
	    hL1TkTau_IsoTks_DeltaPOCAz->Fill( abs( L1Tks_POCAz->at(mTk) - L1Tks_POCAz->at(isoTk) ) );
	    hL1TkTau_IsoTks_d0        ->Fill( s->GetTrackD0(isoTk)      );
	    hL1TkTau_IsoTks_d0Abs     ->Fill( abs(s->GetTrackD0(isoTk)) );
	    hL1TkTau_IsoTks_PtRel     ->Fill( isoTk_PtRel               );
	    hL1TkTau_IsoTks_StubPtCons->Fill( L1Tks_StubPtConsistency->at(isoTk) );
	  }	  

	}
      
    } /// For-loop: Tk-Matched L1TkTaus



    // Fill Turn-On histograms
    FillTurnOn_Denominator_(mcHadronicTaus, hMcHadronicTau_VisEt);
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 25.0, hCalo_TurnOn25  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 25.0, hTk_TurnOn25    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn25);
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 50.0, hCalo_TurnOn50  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 50.0, hTk_TurnOn50    ); 
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn50);
    //
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 66.0, hCalo_TurnOn_SingleTau50KHz  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 65.0, hTk_TurnOn_SingleTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 50.0, hVtxIso_TurnOn_SingleTau50KHz);
    //
    FillTurnOn_Numerator_(L1TkTaus_Calo   , 42.0, hCalo_TurnOn_DiTau50KHz  );
    FillTurnOn_Numerator_(L1TkTaus_Tk     , 40.0, hTk_TurnOn_DiTau50KHz    );
    FillTurnOn_Numerator_(L1TkTaus_VtxIso , 25.0, hVtxIso_TurnOn_DiTau50KHz);


    // SingleTau
    FillSingleTau_(L1TkTaus_Calo  , hCalo_Rate  , hCalo_Eff  );
    FillSingleTau_(L1TkTaus_Tk    , hTk_Rate    , hTk_Eff    );
    FillSingleTau_(L1TkTaus_VtxIso, hVtxIso_Rate, hVtxIso_Eff);
    
    // DiTau
    FillDiTau_(L1TkTaus_Calo  , L1TkTaus_Tk       , hDiTau_Rate_Calo_Tk    , hDiTau_Eff_Calo_Tk    );
    FillDiTau_(L1TkTaus_Calo  , L1TkTaus_VtxIso   , hDiTau_Rate_Calo_VtxIso, hDiTau_Eff_Calo_VtxIso);
    FillDiTau_(L1TkTaus_Tk    , L1TkTaus_VtxIso   , hDiTau_Rate_Tk_VtxIso  , hDiTau_Eff_Tk_VtxIso  );
    
    // WARNING: Removal of non Z-matching should be just before I need it!
    ApplyDiTauZMatching(tk_CollectionType, L1TkTaus_Tk);
    ApplyDiTauZMatching(tk_CollectionType, L1TkTaus_VtxIso);
    FillDiTau_(L1TkTaus_Calo  , hDiTau_Rate_Calo  , hDiTau_Eff_Calo  );
    FillDiTau_(L1TkTaus_Tk    , hDiTau_Rate_Tk    , hDiTau_Eff_Tk    );
    FillDiTau_(L1TkTaus_VtxIso, hDiTau_Rate_VtxIso, hDiTau_Eff_VtxIso);

    // Progress bar
    auxTools_.ProgressBar(jentry, nEntries, 100, 100);
  
  }// For-loop: Entries


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

  // Write the histograms to the file
  outFile->cd();
  outFile->Write();
  auxTools_.StopwatchStop(5, "minutes");

}


//****************************************************************************
void TauTrigger::BookHistos_(void)
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

  // VtxIsolated L1TkTaus, Signal TTTracks
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_DeltaR       , "L1TkTau_MatchPixTk_DeltaR"       , 2000,    +0.0,    +2.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_PtRel        , "L1TkTau_MatchPixTk_PtRel"        , 2000,    +0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_StubPtCons   , "L1TkTau_MatchPixTk_StubPtCons"   ,  200,    +0.0,  +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_NStubs       , "L1TkTau_MatchPixTk_NStubs"       ,   30,    -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_NPsStubs     , "L1TkTau_MatchPixTk_NPsStubs"     ,   30,    -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_NBarrelStubs , "L1TkTau_MatchPixTk_NBarrelStubs" ,   30,    -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_NEndcapStubs , "L1TkTau_MatchPixTk_NEndcapStubs" ,   30,    -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_NPixHits     , "L1TkTau_MatchPixTk_NPixHits"     ,   30,    -0.5,   +14.5  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_Pt           , "L1TkTau_MatchPixTk_Pt"           ,  300,    +0.0,  +300.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_Eta          , "L1TkTau_MatchPixTk_Eta"          ,  600,    -3.0,    +3.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_POCAz        , "L1TkTau_MatchPixTk_POCAz"        ,  600,   -30.0,   +30.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_POCAzSig     , "L1TkTau_MatchPixTk_POCAzSig"     , 2000, -1000.0, +1000.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_d0           , "L1TkTau_MatchPixTk_d0"           , 2000,   -10.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_d0Abs        , "L1TkTau_MatchPixTk_d0Abs"        , 1000,     0.0,   +10.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_d0Sig        , "L1TkTau_MatchPixTk_d0Sig"        , 4000,   -20.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_d0SigAbs     , "L1TkTau_MatchPixTk_d0SigAbs"     , 2000,     0.0,   +20.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_ChiSquared   , "L1TkTau_MatchPixTk_ChiSquared"   ,  200,    +0.0,  +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_RedChiSquared, "L1TkTau_MatchPixTk_RedChiSquared", 2000,    +0.0,  +200.0  );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_SigmaRInv    , "L1TkTau_MatchPixTk_SigmaRInv"    , 2000,    +0.0,    +0.01 );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_SigmaPhi0    , "L1TkTau_MatchPixTk_SigmaPhi0"    , 5000,    +0.0,    +0.05 ); 
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_SigmaD0      , "L1TkTau_MatchPixTk_SigmaD0"      , 5000,    +0.0,    +0.05 );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_SigmaT       , "L1TkTau_MatchPixTk_SigmaT"       , 5000,    +0.0,    +0.05 );
  histoTools_.BookHisto_1D(hL1TkTau_MatchPixTk_SigmaZ0      , "L1TkTau_MatchPixTk_SigmaZ0"      , 5000,    +0.0,    +0.05 );
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

  // VtxIsolated L1TkTaus, Signal TTPixelTracks
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
void TauTrigger::FinaliseEffHisto_(TH1D *histo, 
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
void TauTrigger::FinaliseEffHisto_(TH2D *histo, 
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
void TauTrigger::ApplyDiTauZMatching(const string tkCollectionType, 
				     vector<L1TkTauParticle> &L1TkTaus)
//****************************************************************************
{
  
  // cout << "FIXME! more complicated function needed" << endl;  //xenios 

  if (L1TkTaus.size() < 2) return;
  const int matchTk_index0 = L1TkTaus.at(0).matchTk_Index_;
  double deltaPOCAz = 9999.9;

  // For-loop: L1TkTaus
  for( int i=1; i < (int) L1TkTaus.size(); i++){
    
    // Variables
    const int matchTk_index = L1TkTaus.at(i).matchTk_Index_;

    if ( tkCollectionType.compare("TTPixelTracks") == 0 ) {
      deltaPOCAz = abs( L1PixTks_POCAz->at(matchTk_index) - L1PixTks_POCAz->at(matchTk_index0) );
    }
    else if ( tkCollectionType.compare("TTTracks") == 0 ) {
      deltaPOCAz = abs( L1Tks_POCAz->at(matchTk_index) - L1Tks_POCAz->at(matchTk_index0) );
    }
    else{
      cout << "E R R O R ! TauTrigger::ApplyDiTauZMatching(...) - Unknown sample \"" << mcSample << "\". EXIT" << endl;
      exit(1);
    }
    
    // If the Trigger objects is not within 1.0 cm reject it
    if (deltaPOCAz > diTau_deltaPOCAz) L1TkTaus.erase ( L1TkTaus.begin()+i );
  }  // For-loop: L1TkTaus
  
  return;
}


//****************************************************************************
void TauTrigger::FillSingleTau_(const vector<L1TkTauParticle> L1TkTaus, 
				TH1D *hRate,
				TH1D *hEfficiency)
//****************************************************************************
{

  if( L1TkTaus.size() == 0 ) return;
  
  // Fill rate 
  const Double_t ldgEt = L1CaloTau_Et->at(L1TkTaus.at(0).caloTau_Index_);
  FillRate_(hRate, ldgEt);
  
  // Get MC-matched trigger objects above given Et Threshold
  vector<L1TkTauParticle> L1TkTaus_mcMatched;
  GetMcMatchedL1TkTaus(L1TkTaus, L1TkTaus_mcMatched);
  if (L1TkTaus_mcMatched.size() < 1) return;

  // Fill efficiency
  if(!bFoundAllTaus_) return;
  Double_t ldgEt_mcMatched = L1CaloTau_Et->at(L1TkTaus_mcMatched.at(0).caloTau_Index_);
  FillEfficiency_(hEfficiency, ldgEt_mcMatched);

  return;
}


//****************************************************************************
void TauTrigger::FillDiTau_(const vector<L1TkTauParticle> L1TkTaus, 
			    TH1D *hRate,
			    TH1D *hEfficiency)
//****************************************************************************
{

  if( L1TkTaus.size() < 2 ) return;  

  // Fill rate 
  const Double_t subLdgEt = L1CaloTau_Et->at(L1TkTaus.at(1).caloTau_Index_);
  FillRate_(hRate, subLdgEt);

  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1TkTaus_mcMatched;
  GetMcMatchedL1TkTaus(L1TkTaus, L1TkTaus_mcMatched);
  if (L1TkTaus_mcMatched.size() < 2) return;
    
  // Fill efficiency
  if(!bFoundAllTaus_) return;
  Double_t subLdgEt_mcMatched = L1CaloTau_Et->at(L1TkTaus_mcMatched.at(1).caloTau_Index_);
  FillEfficiency_(hEfficiency, subLdgEt_mcMatched);

  return;
}


//****************************************************************************
void TauTrigger::FillDiTau_(const vector<L1TkTauParticle> L1TkTaus1,
			    const vector<L1TkTauParticle> L1TkTaus2, 
			    TH2D *hRate,
			    TH2D *hEfficiency)
//****************************************************************************
{

  if( L1TkTaus1.size() < 1 ) return;
  if( L1TkTaus2.size() < 1 ) return;
  
  vector<L1TkTauParticle> L1TkTaus1_ = L1TkTaus1;
  
  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1TkTaus_mcMatched1;
  vector<L1TkTauParticle> L1TkTaus_mcMatched2;
  GetMcMatchedL1TkTaus(L1TkTaus1_, L1TkTaus_mcMatched1);
  GetMcMatchedL1TkTaus(L1TkTaus2, L1TkTaus_mcMatched2);

  // Fill rate 
  RemoveDuplicates(L1TkTaus2, L1TkTaus1_);
  if( L1TkTaus1_.size() < 1 ) return;
  
  Double_t ldgEt1 = L1CaloTau_Et->at(L1TkTaus1_.at(0).caloTau_Index_);
  Double_t ldgEt2 = L1CaloTau_Et->at(L1TkTaus2.at(0).caloTau_Index_);

  // Make x-axis the ldgEt axis
  if (ldgEt1 > ldgEt2)      FillRate_(hRate, ldgEt1, ldgEt2); 
  else if (ldgEt1 < ldgEt2) FillRate_(hRate, ldgEt2, ldgEt1);
  else{
    if (L1TkTaus1_.at(0).caloTau_Index_ == L1TkTaus2.at(0).caloTau_Index_ ) {
      cout << "E R R O R ! TauTrigger::FillDiTau_(...) - Leading and Sub-Leading Trigger objects have the same eT (" << ldgEt1 << " , " << ldgEt2 << ")! "
	   << "This cannot happen. Exit." << endl;
      exit(1);
    }
  }

  // Get MC-matched trigger objects above given Et Threshold
  if (L1TkTaus_mcMatched1.size() < 1) return;
  if (L1TkTaus_mcMatched2.size() < 1) return;
  RemoveDuplicates(L1TkTaus_mcMatched2, L1TkTaus_mcMatched1);
  if (L1TkTaus_mcMatched1.size() < 1) return;

  Double_t ldgEt_mcMatched1 = L1CaloTau_Et->at(L1TkTaus_mcMatched1.at(0).caloTau_Index_);
  Double_t ldgEt_mcMatched2 = L1CaloTau_Et->at(L1TkTaus_mcMatched2.at(0).caloTau_Index_);

  // Only fill if all signal taus are within acceptance
  if(!bFoundAllTaus_) return;
  if (ldgEt_mcMatched1 > ldgEt_mcMatched2)      histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt_mcMatched1, ldgEt_mcMatched2);
  else if (ldgEt_mcMatched1 < ldgEt_mcMatched2) histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt_mcMatched2, ldgEt_mcMatched1);
  else{
    if (L1TkTaus_mcMatched1.at(0).caloTau_Index_ == L1TkTaus_mcMatched2.at(0).caloTau_Index_ ) {
      cout << "E R R O R ! TauTrigger::FillDiTau_(...) - Leading and Sub-Leading MC-Matched Trigger objects have the same eT (" << ldgEt_mcMatched1 << " , " << ldgEt_mcMatched2 << ")! "
	   << "This cannot happen. Exit." << endl;
      exit(1);
    }
  }

  return;
}


//****************************************************************************
void TauTrigger::RemoveDuplicates(const vector<L1TkTauParticle> L1TkTaus2, 
				  vector<L1TkTauParticle> &L1TkTaus1)
//****************************************************************************
{
  
  if (L1TkTaus1.size() < 1) return;
  if (L1TkTaus2.size() < 1) return;

  if ( L1TkTaus1.size() < L1TkTaus2.size() )
    {
      cout << "E R R O R ! TauTrigger::RemoveDuplicates(...) - Size of L1TkTaus1 (" << L1TkTaus1.size() << ") is smaller than the size of L1TkTaus2 (" << L1TkTaus2.size() << "). "
	   <<  "Perhaps you should reverse their order when passing them as arguments to this function? EXIT" << endl;
      exit(1);
    }

  int calo_index1 = L1TkTaus1.at(0).caloTau_Index_;
  int calo_index2 = L1TkTaus2.at(0).caloTau_Index_;
  if( calo_index1 == calo_index2) L1TkTaus1.erase ( L1TkTaus1.begin()+0 );

  return;
}


//****************************************************************************
void TauTrigger::GetMcMatchedL1TkTaus(const vector<L1TkTauParticle> L1TkTaus, 
				      vector<L1TkTauParticle> &L1TkTaus_McMatched)
//****************************************************************************
{
  
  if (L1TkTaus.size() < 1) return;
  for (int iTau = 0; iTau < (int) L1TkTaus.size(); iTau++)
    {      
      L1TkTauParticle L1TkTau = L1TkTaus.at(iTau);
      if (L1TkTau.matchGenp_Index_ < 0) continue;
      L1TkTaus_McMatched.push_back(L1TkTau);
    }
  
  return;
}


//****************************************************************************
void TauTrigger::FillRate_(TH1D *hRate,
			   const Double_t ldgEt)
//****************************************************************************
{
  
  if (ldgEt < 0) return;
  hRate ->Fill( ldgEt );
  
  return;
}


//****************************************************************************
void TauTrigger::FillRate_(TH2D *hRate,
			   const Double_t ldgEt1,
			   const Double_t ldgEt2)
//****************************************************************************
{
  
  if (ldgEt1 < 0) return;
  if (ldgEt2 < 0) return;

  hRate ->Fill( ldgEt1, ldgEt2 );
  
  return;
}


//****************************************************************************
void TauTrigger::FillEfficiency_(TH1D *hEfficiency,
				 const Double_t ldgEt)
//****************************************************************************
{
  
  histoTools_.FillAllBinsUpToValue_1D(hEfficiency, ldgEt);

  return;
}


//****************************************************************************
vector<int> TauTrigger::GetTriggerHadronicTaus_(const int tauMom, 					  
						const double hTauVisEtCut)
//****************************************************************************
{

  // Variable declaration
  vector<int> hTaus_index;
  const bool bApplyEtaCut      = true;
  const bool bUseAbsoluteMomId = true;

  // For-loop: GenParticles
  for (Size_t iGenP = 0; iGenP < GenP_PdgId->size(); iGenP++) {
    
    // Get hadronically decaying taus within acceptance
    bool bIsHTau = IsTriggerHadronicTau(iGenP, bApplyEtaCut, tauMom, bUseAbsoluteMomId, hTauVisEtCut);
    if(!bIsHTau) continue;
    
    // Push back GenP index hadronic tau 
    hTaus_index.push_back(iGenP);
    
  } // For-loop: GenParticles
  
  return hTaus_index;
}


//****************************************************************************
TLorentzVector TauTrigger::GetHadronicTauVisP4_(const int genP_index)
//****************************************************************************
{

  // Get visible 4-momentum of hadronic tau
  vector<unsigned short> hTau_Prods;
  GetHadronicTauFinalDaughters(genP_index, hTau_Prods);
  TLorentzVector visP4 = GetVisibleP4(hTau_Prods);

  return visP4;
}


//****************************************************************************
void TauTrigger::GetL1CaloTauMatchGenp(const int calo_index,
				       vector<int> hTaus_index,
				       int &matchGenp_index, 
				       double &matchGenp_deltaR)
//****************************************************************************
{
  
  // Initialise return values
  matchGenp_deltaR = 9999.9;

  if( calo_index < 0 ) return;
  if( hTaus_index.size() < 1) return;

  // For-loop: VisP4 of hTaus
  for (size_t i = 0; i < hTaus_index.size(); i++) {
    
    const int iTau = hTaus_index.at(i);

    // Is the CaloTau deltaR matched to the hTau?
    TLorentzVector hTau_VisP4 = GetHadronicTauVisP4_(iTau);
    double deltaR = auxTools_.DeltaR( L1CaloTau_Eta->at(calo_index), L1CaloTau_Phi->at(calo_index), hTau_VisP4.Eta(), hTau_VisP4.Phi() );
    if (deltaR > mcMatching_maxDeltaR) continue;
    matchGenp_index  = iTau;
    matchGenp_deltaR = deltaR;

    if(DEBUG) {
      PrintL1CaloTau(calo_index);
      PrintGenpMinimalInfo(iTau);
      cout << "    matchGenp_index = " << matchGenp_index << "\tmatchGenp_deltaR = " << matchGenp_deltaR << endl;
    }

  
  } // For-loop: VisP4 of hTaus

  return;
}      


//****************************************************************************
void TauTrigger::GetL1CaloTauUniqueMatchGenp(vector<L1TkTauParticle> &L1TkTaus)
//****************************************************************************
{

  if(DEBUG) {  
    Table table1("CaloTau Index | GenP-Match Index | Delta-R", "Text");
    // For-loop: L1TkTaus
    for(int i = 0; i < (int) L1TkTaus.size(); i++){

      int index     = L1TkTaus.at(i).matchGenp_Index_;
      double deltaR = L1TkTaus.at(i).matchGenp_deltaR_;
      table1.AddRowColumn(i, auxTools_.ToString(i) );
      table1.AddRowColumn(i, auxTools_.ToString(index) );
      table1.AddRowColumn(i, auxTools_.ToString(deltaR) );
    }
    table1.Print();
  }


  // For-loop: L1TkTaus
  for(int i = 0; i < (int) L1TkTaus.size(); i++){

    int mcMatch_index     = L1TkTaus.at(i).matchGenp_Index_;
    double mcMatch_deltaR = L1TkTaus.at(i).matchGenp_deltaR_;
    
    // Skip L1TkTau if it has no MC Tau associated
    if (mcMatch_index < 0) continue;

    // For-loop: L1TkTaus
    for(int j = 0; j < (int) L1TkTaus.size(); j++){

      int mcMatch_index_tmp     = L1TkTaus.at(j).matchGenp_Index_;
      double mcMatch_deltaR_tmp = L1TkTaus.at(j).matchGenp_deltaR_;

      // Skip if self or if it has no MC Tau associated
      if (i==j) continue;
      if (mcMatch_index_tmp < 0) continue;
      
      // De-associated the MC-tau with this L1TkTau
      if (mcMatch_index_tmp != mcMatch_index) continue;
      if (mcMatch_deltaR_tmp > mcMatch_deltaR) continue;
      L1TkTaus.at(i).SetMatchGenp(-1.0, 9999.9);

    } // For-loop: L1TkTaus

  } // For-loop: L1TkTaus

  if(DEBUG) {
    Table table2("CaloTau Index | GenP-Match Index (Unique) | Delta-R (Unique)", "Text");
    // For-loop: L1TkTaus
    for(int i = 0; i < (int) L1TkTaus.size(); i++){

      int index     = L1TkTaus.at(i).matchGenp_Index_;
      double deltaR = L1TkTaus.at(i).matchGenp_deltaR_;
      table2.AddRowColumn(i, auxTools_.ToString(i) );
      table2.AddRowColumn(i, auxTools_.ToString(index) );
      table2.AddRowColumn(i, auxTools_.ToString(deltaR) );
    }
    table2.Print();
  }
    
  return;
}      


//****************************************************************************
void TauTrigger::FillTurnOn_Numerator_(const vector<L1TkTauParticle> L1TkTaus, 
				       const double minEt,
				       TH1D *hTurnOn)
//****************************************************************************
{

  if (L1TkTaus.size() < 1) return;

  // For-loop: L1TkTaus
  for(int iTau = 0; iTau < (int) L1TkTaus.size(); iTau++){
      
    // Skip if trigger object is not MC matched
    const L1TkTauParticle L1TkTau = L1TkTaus.at(iTau);
    if ( L1TkTau.matchGenp_Index_ < 0) continue;
    
    // Skip if trigger object has eT < minEt
    double calo_et = L1CaloTau_Et->at(L1TkTau.caloTau_Index_);
    if ( calo_et < minEt) continue;
    
    // Fill histo with visible eT of mathing MC Tau
    const TLorentzVector hTau_visP4 = GetHadronicTauVisP4_( L1TkTau.matchGenp_Index_ );
    // const TLorentzVector hTau_P4 = GetP4( L1TkTau.matchGenp_Index_ );
    hTurnOn->Fill( hTau_visP4.Et() );
    
    if(DEBUG) {
      PrintGenpMinimalInfo(L1TkTau.matchGenp_Index_);
      PrintAllDaughtersMinimalInfo(L1TkTau.matchGenp_Index_);
      PrintL1CaloTau(L1TkTau.caloTau_Index_);
    }

  } // For-loop: L1TkTaus
  
  return;
   
}


//****************************************************************************
void TauTrigger::FillTurnOn_Denominator_(vector<int> mcHadronicTaus,
					 TH1D *hVisEt)
//****************************************************************************
{

  if (mcHadronicTaus.size() < 1) return;

  // For-loop: MC Haronic Taus
  for (Size_t i = 0; i < mcHadronicTaus.size(); i++){
    
    int iTau  = mcHadronicTaus.at(i);
    TLorentzVector hTau_visP4 = GetHadronicTauVisP4_(iTau);
    hVisEt->Fill( hTau_visP4.Et() );    

  }// For-loop: MC Haronic Taus

  return;

}

//**************************************************
void TauTrigger::PrintL1CaloTau(int Indx)
//**************************************************
{

  double E   = L1CaloTau_E  ->at(Indx);
  double Et  = L1CaloTau_Et ->at(Indx);
  double Eta = L1CaloTau_Eta->at(Indx);
  double Phi = L1CaloTau_Phi->at(Indx);
  // Table caloTauInfo("Event | CaloTau Index | Energy | Et | Eta | Phi", "Text");
  // caloTauInfo.AddRowColumn(0, auxTools_.ToString( EvtNumber ) );
  // caloTauInfo.AddRowColumn(0, auxTools_.ToString( Indx ) );
  // caloTauInfo.AddRowColumn(0, auxTools_.ToString( E ) );
  // caloTauInfo.AddRowColumn(0, auxTools_.ToString( Et ) );
  // caloTauInfo.AddRowColumn(0, auxTools_.ToString( Eta ) );
  // caloTauInfo.AddRowColumn(0, auxTools_.ToString( Phi ) );
  // caloTauInfo.Print();
  Table caloTauInfo("iCalo | Energy | Et | Eta | Phi", "Text");
  // caloTauInfo.AddRowColumn(0, auxTools_.ToString( EvtNumber ) );
  caloTauInfo.AddRowColumn(0, auxTools_.ToString( Indx ) );
  caloTauInfo.AddRowColumn(0, auxTools_.ToString( E ) );
  caloTauInfo.AddRowColumn(0, auxTools_.ToString( Et ) );
  caloTauInfo.AddRowColumn(0, auxTools_.ToString( Eta ) );
  caloTauInfo.AddRowColumn(0, auxTools_.ToString( Phi ) );
  caloTauInfo.Print();
  

  return;
}

//**************************************************
bool TauTrigger::EventWithForbiddenTauDecayMode(string allowedTauDecayMode, 
						vector<int> hTaus_index)
//**************************************************
{

  // 1-prong or 3-prong?
  if (mcSample.compare("nugun") == 0) return false;  
  if (allowedTauDecayMode.compare("all") == 0) return false;

  bool bFoundForbiddenMode = false;
  // For-loop: All hadronic taus
  for (Size_t i = 0; i < hTaus_index.size(); i++) {
    
    // Get charged pions
    const int tauIndex = hTaus_index.at(i);
    vector<unsigned short> chargedPions;
    GetHadronicTauChargedPions(tauIndex, chargedPions);
	
    if (allowedTauDecayMode.compare("1p") == 0){ if (chargedPions.size() != 1) bFoundForbiddenMode = true;}
    else if (allowedTauDecayMode.compare("3p") == 0){ if (chargedPions.size() != 3) bFoundForbiddenMode =  true;}
    else{
      cout << "E R R O R ! TauTrigger::EventWithForbiddenTauDecayMode(...) - Unknown tau decay mode \"" << allowedTauDecayMode << "\". EXIT" << endl;
      exit(1);
    }

  } // For-loop: All hadronic taus

  return bFoundForbiddenMode;
}


//****************************************************************************
void TauTrigger::GetL1TkTauMatchTrack(int &matchTk_index,
				      const int calo_index,
				      const string matchTk_CollectionType,
				      const double matchTk_minPt,
				      const double matchTk_maxDeltaR,
				      double &matchTk_deltaR)
//****************************************************************************
{

  // Init variables
  matchTk_deltaR       = 999.99;
  double matchTk_pt    = 0.0;
  double matchTk_PtRel = 0.0;

  // Get calo properties
  const double calo_eta = L1CaloTau_Eta->at(calo_index);
  const double calo_phi = L1CaloTau_Phi->at(calo_index);
  const TLorentzVector calo_p4 = auxTools_.GetTLorentzVector( L1CaloTau_Et->at(calo_index), L1CaloTau_Eta->at(calo_index), L1CaloTau_Phi->at(calo_index), L1CaloTau_E->at(calo_index) );

  // For-Loop: L1Tks  
  for (Size_t iTk = 0; iTk < L1Tks_Pt->size(); iTk++) {
    
    // Get tk qualities
    const double tk_pt     = L1Tks_Pt->at(iTk);
    const double tk_eta    = L1Tks_Eta->at(iTk);
    const double tk_phi    = L1Tks_Phi->at(iTk);
    const TVector3 tk_p3   = auxTools_.GetTVector3( L1Tks_Px->at(iTk), L1Tks_Py->at(iTk), L1Tks_Pz->at(iTk) );
    double tk_PtRel        = tk_p3.Perp( calo_p4.Vect() );

    // Apply pT requirement
    bool bIsBelowMinPt   = (tk_pt < matchTk_minPt);
    if (bIsBelowMinPt) continue;

    // Apply track quality criteria
    bool bIsSignalTk = s->SelectTracks(iTk, matchTk_CollectionType);
    if (!bIsSignalTk) continue;

    // Attempt to deltaR-match CaloTau with a track
    const double deltaR = auxTools_.DeltaR(calo_eta, calo_phi, tk_eta, tk_phi);

    // Only look at tracks within a matching cone of:  0 <= DeltaR < matchTk_maxDeltaR
    if (deltaR > matchTk_maxDeltaR) continue;

    // Find track which is closest to the calo. Save deltaR and index
    if (deltaR < matchTk_deltaR) 
      {
	matchTk_pt     = tk_pt;
	matchTk_deltaR = deltaR;
	matchTk_index  = iTk;
	matchTk_PtRel  = tk_PtRel;
      }
    
  }// For-Loop: L1Tks

  return;
}


//****************************************************************************
void TauTrigger::GetL1TkTauMatchPixTrack(int &matchTk_index,
					 const int calo_index,
					 const string matchTk_CollectionType,
					 const double matchTk_minPt,
					 const double matchTk_maxDeltaR,
					 double &matchTk_deltaR)
//****************************************************************************
{

  // Init variables
  matchTk_deltaR       = 999.99;
  double matchTk_pt    = 0.0;
  double matchTk_PtRel = 0.0;

  // Get calo properties
  const double calo_eta = L1CaloTau_Eta->at(calo_index);
  const double calo_phi = L1CaloTau_Phi->at(calo_index);
  const TLorentzVector calo_p4 = auxTools_.GetTLorentzVector( L1CaloTau_Et->at(calo_index), L1CaloTau_Eta->at(calo_index), L1CaloTau_Phi->at(calo_index), L1CaloTau_E->at(calo_index) );

  // For-Loop: L1PixTks  
  for (Size_t iTk = 0; iTk < L1PixTks_Pt->size(); iTk++) {
    
    // Get tk qualities
    const double tk_pt     = L1PixTks_Pt->at(iTk);
    const double tk_eta    = L1PixTks_Eta->at(iTk);
    const double tk_phi    = L1PixTks_Phi->at(iTk);
    const TVector3 tk_p3   = auxTools_.GetTVector3( L1PixTks_Px->at(iTk), L1PixTks_Py->at(iTk), L1PixTks_Pz->at(iTk) );
    double tk_PtRel        = tk_p3.Perp( calo_p4.Vect() );

    // Apply pT requirement
    bool bIsBelowMinPt   = (tk_pt < matchTk_minPt);
    if (bIsBelowMinPt) continue;

    // Apply track quality criteria
    bool bIsSignalTk = s->SelectPixTracks(iTk, matchTk_CollectionType);
    if (!bIsSignalTk) continue;

    // Attempt to deltaR-match CaloTau with a track
    const double deltaR = auxTools_.DeltaR(calo_eta, calo_phi, tk_eta, tk_phi);

    // Only look at tracks within a matching cone of:  0 <= DeltaR < matchTk_maxDeltaR
    if (deltaR > matchTk_maxDeltaR) continue;

    // Find track which is closest to the calo. Save deltaR and index
    if (deltaR < matchTk_deltaR) 
      {
	matchTk_pt     = tk_pt;
	matchTk_deltaR = deltaR;
	matchTk_index  = iTk;
	matchTk_PtRel  = tk_PtRel;
      }
    
  }// For-Loop: L1PixTks

  return;
}


//****************************************************************************
void TauTrigger::GetL1TkTauIsolation(const int matchTk_index,
				     const vector<int> isoTks_index,
				     double &relIso,
				     double &vtxIso)
//****************************************************************************
{

  relIso = 0.0;
  vtxIso = 999.9;
  if (isoTks_index.size() < 1) return;

  // Initialise variables
  vector< double > v_deltaPOCAz;  
  v_deltaPOCAz.push_back(999.9); // for cases where no track is found
  double isoTk_scalarSumPt = 0.0; 

  // Get ldgTk properties
  const double matchTk_Pt    = L1Tks_Pt   ->at(matchTk_index);
  const double matchTk_POCAz = L1Tks_POCAz->at(matchTk_index);

  // For-loop: IsoTks
  for (Size_t iTk = 0; iTk < isoTks_index.size(); iTk++) { 
	
    const int isoTk_index    = isoTks_index.at(iTk);
    const double isoTk_Pt    = L1Tks_Pt->at(isoTk_index);
    const double isoTk_POCAz = L1Tks_POCAz->at(isoTk_index);

    // Add up the pT of all (selected) tracks in isolation cone/annulus
    isoTk_scalarSumPt += isoTk_Pt;

    // Store delta-POCAz of isoTk from matchTk
    v_deltaPOCAz.push_back(matchTk_POCAz - isoTk_POCAz);
    
  }// For-loop: IsoTks


  // Calculated vertex isolation
  std::sort(v_deltaPOCAz.begin(), v_deltaPOCAz.end(), SortDescendingAbs() );
  vtxIso = v_deltaPOCAz[0];

  // Calculated relative isolation
  relIso = isoTk_scalarSumPt/matchTk_Pt;

  return;
}


//****************************************************************************
void TauTrigger::GetL1TkTauPixIsolation(const int matchTk_index,
					const vector<int> isoTks_index,
					double &relIso,
					double &vtxIso)
//****************************************************************************
{
  
  relIso = 0.0;
  vtxIso = 999.9;
  if (isoTks_index.size() < 1) return;

  // Initialise variables
  vector<double> v_deltaPOCAz;
  v_deltaPOCAz.push_back(999.9); // for cases where no track is found
  double isoTk_scalarSumPt = 0.0; 

  // Get matching track properties
  const double matchTk_Pt    = L1PixTks_Pt   ->at(matchTk_index);
  const double matchTk_POCAz = L1PixTks_POCAz->at(matchTk_index);

  // For-loop: IsoTks
  for (Size_t iTk = 0; iTk < isoTks_index.size(); iTk++) { 
	
    const int isoTk_index    = isoTks_index.at(iTk);
    const double isoTk_Pt    = L1PixTks_Pt->at(isoTk_index);
    const double isoTk_POCAz = L1PixTks_POCAz->at(isoTk_index);
    
    // Add up the pT of all (selected) tracks in isolation cone/annulus
    isoTk_scalarSumPt += isoTk_Pt;

    // Store delta-POCAz of isoTk from matchTk
    v_deltaPOCAz.push_back(matchTk_POCAz - isoTk_POCAz);
    
  }// For-loop: IsoTks


  // Calculated vertex isolation
  std::sort(v_deltaPOCAz.begin(), v_deltaPOCAz.end(), SortDescendingAbs() );
  vtxIso = v_deltaPOCAz[0];

  // Calculated relative isolation
  relIso = isoTk_scalarSumPt/matchTk_Pt;

  return;
}


//****************************************************************************
void TauTrigger::GetL1TkTauAlgoConeSizes(const int calo_index,
					 const string algoType,
					 const double sigCone_dRMax,
					 double &sig_minDeltaR_tmp,
					 double &sig_maxDeltaR_tmp,
					 double &iso_minDeltaR_tmp)
//****************************************************************************
{

  double deltaR_sigMin = -1.0;
  double deltaR_sigMax = -1.0;

  if (algoType.compare("FixedCone") == 0){ 
    deltaR_sigMin = 0.0;
    deltaR_sigMax = sigCone_maxDeltaR;
  }  
  else if (algoType.compare("ShrinkingCone") == 0){ 
    deltaR_sigMin = (sigCone_ShrinkConstant)/(L1CaloTau_Et->at(calo_index));
    deltaR_sigMax = (isoCone_ShrinkConstant)/(L1CaloTau_Et->at(calo_index));
    if (deltaR_sigMax > sigCone_dRMax) deltaR_sigMax = sigCone_dRMax;
  }
  else{
    cout << "E R R O R ! TauTrigger::GetL1TkTauAlgoConeSizes(...) - Unknown algorithm type \"" << algoType << "\". EXIT" << endl;
    exit(1);
  }

  // Assign signal and isolation cone sizes
  sig_minDeltaR_tmp = deltaR_sigMin;
  sig_maxDeltaR_tmp = deltaR_sigMax;
  iso_minDeltaR_tmp = deltaR_sigMax;
#ifdef DEBUG
  cout << "I N F O ! TauTrigger::GetL1TkTauAlgoConeSizes(...) - New max signal cone (min isolation cone ) size is " << sig_maxDeltaR << " (" << iso_minDeltaR << ")" << endl;
#endif

  return;
}


//****************************************************************************
void TauTrigger::GetL1TkTauSigConeTracks(const int matchTk_index, 
					 const double deltaR_min,
					 const double deltaR_max,
					 const string selectionType, 
					 vector<int> &sigTks_index)
//****************************************************************************
{
  // Add the matchTk to the sigTks
  if (matchTk_index == -1) return;
  else sigTks_index.push_back(matchTk_index);

  // Get matchTk properties
  const double matchTk_pt    = L1Tks_Pt ->at(matchTk_index);
  const double matchTk_eta   = L1Tks_Eta->at(matchTk_index);
  const double matchTk_phi   = L1Tks_Phi->at(matchTk_index);
  const double matchTk_mass  = constants::pionMass;
  const double matchTk_POCAz = L1Tks_POCAz->at( matchTk_index );
  TLorentzVector tau_p4;
  tau_p4.SetPtEtaPhiM(matchTk_pt, matchTk_eta, matchTk_phi, matchTk_mass);

  // For-Loop: L1Tks  
  for (Size_t iTk = 0; iTk < L1Tks_Pt->size(); iTk++) {
    
    bool bIsQualityTk = s->SelectTracks(iTk, selectionType);
    if (!bIsQualityTk) continue;

    // Check if track has already being used
    bool bIsOwnSigTk  = find(sigTks_index.begin(), sigTks_index.end(), iTk) != sigTks_index.end();
    if (bIsOwnSigTk) continue;

    // Get Track properties
    double tk_pt      = L1Tks_Pt->at(iTk);
    double tk_eta     = L1Tks_Eta->at(iTk);
    double tk_phi     = L1Tks_Phi->at(iTk);
    double tk_POCAz   = L1Tks_POCAz->at(iTk);
    double deltaR     = auxTools_.DeltaR(matchTk_eta, matchTk_phi, tk_eta, tk_phi);
    double deltaPOCAz = abs(matchTk_POCAz - tk_POCAz);
    TLorentzVector tmp_p4;
    tmp_p4.SetPtEtaPhiM(tk_pt, tk_eta, tk_phi, constants::pionMass);
    
    // Is track within the signal cone?
    bool bIsInsideCone = ( (deltaR <= deltaR_max) && (deltaR >= deltaR_min) );
    if (!bIsInsideCone) continue;

    // Is track close enough in POCA-z?
    if (deltaPOCAz > sigCone_maxTkDeltaPOCAz) continue;

    // Invariant mass cut
    if ( (tau_p4 + tmp_p4).M() > sigCone_maxTkInvMass) continue;

    // Save the track as a signal-cone track
    tau_p4 = tau_p4 + tmp_p4;
    sigTks_index.push_back(iTk);
    
  } // For-Loop: L1Tks 


#ifdef DEBUG
  auxTools.PrintVector( sigTks_index, "\tSigTks: ");
#endif
  
  return;
}



//****************************************************************************
void TauTrigger::GetL1TkTauSigConePixTracks(const int matchTk_index, 
					    const double deltaR_min,
					    const double deltaR_max,
					    const string selectionType, 
					    vector<int> &sigTks_index)
//****************************************************************************
{

  // Add the matchTk to the sigTks
  if (matchTk_index == -1) return;
  else sigTks_index.push_back(matchTk_index);

  // Get matchTk properties
  const double matchTk_pt    = L1PixTks_Pt ->at(matchTk_index);
  const double matchTk_eta   = L1PixTks_Eta->at(matchTk_index);
  const double matchTk_phi   = L1PixTks_Phi->at(matchTk_index);
  const double matchTk_mass  = constants::pionMass;
  const double matchTk_POCAz = L1PixTks_POCAz->at( matchTk_index );
  TVector3 matchTk_p3        = auxTools_.GetTVector3( L1PixTks_Px->at(matchTk_index), L1PixTks_Py->at(matchTk_index), L1PixTks_Pz->at(matchTk_index) );
  TLorentzVector tau_p4;
  tau_p4.SetPtEtaPhiM(matchTk_pt, matchTk_eta, matchTk_phi, matchTk_mass);

  // For-Loop: L1PixTks  
  for (Size_t iTk = 0; iTk < L1PixTks_Pt->size(); iTk++) {
    
    bool bIsQualityTk = s->SelectPixTracks(iTk, selectionType);
    if (!bIsQualityTk) continue;

    // Check if track has already being used
    bool bIsOwnSigTk  = find(sigTks_index.begin(), sigTks_index.end(), iTk) != sigTks_index.end();
    if (bIsOwnSigTk) continue;

    // Get Track properties
    double tk_pt      = L1PixTks_Pt->at(iTk);
    double tk_eta     = L1PixTks_Eta->at(iTk);
    double tk_phi     = L1PixTks_Phi->at(iTk);
    double tk_POCAz   = L1PixTks_POCAz->at(iTk);
    double deltaR     = auxTools_.DeltaR(matchTk_eta, matchTk_phi, tk_eta, tk_phi);
    double deltaPOCAz = abs(matchTk_POCAz - tk_POCAz);
    TVector3 tk_p3    = auxTools_.GetTVector3( L1PixTks_Px->at(iTk), L1PixTks_Py->at(iTk), L1PixTks_Pz->at(iTk) );
    TLorentzVector tmp_p4;
    tmp_p4.SetPtEtaPhiM(tk_pt, tk_eta, tk_phi, constants::pionMass);
    
    // Is track within the signal cone?
    bool bIsInsideCone = ( (deltaR <= deltaR_max) && (deltaR >= deltaR_min) );
    if (!bIsInsideCone) continue;

    // Is track close enough in POCA-z?
    if (deltaPOCAz > sigCone_maxTkDeltaPOCAz) continue;
    
    // Invariant mass cut
    if ( (tau_p4 + tmp_p4).M() > sigCone_maxTkInvMass) continue;

    // Save the track as a signal-cone track
    tau_p4 = tau_p4 + tmp_p4;
    sigTks_index.push_back(iTk);
    
  } // For-Loop: L1PixTks 


#ifdef DEBUG
  auxTools.PrintVector( sigTks_index, "\tSigTks: ");
#endif
  
  return;
}


//****************************************************************************
void TauTrigger::GetL1TkTauIsoConeTracks(const int matchTk_index, 
					 const double deltaR_min,
					 const double deltaR_max,
					 const string selectionType,
					 vector<int> &isoTks_index)
//****************************************************************************
{
  
#ifdef DEBUG
  for (Size_t l = 0; l < isoTks_index.size(); l++)  PrintTrackProperties(isoTks_index.at(l));
#endif

  if (matchTk_index == -1) return;
  const double matchTk_eta   = L1Tks_Eta->at(matchTk_index);
  const double matchTk_phi   = L1Tks_Phi->at(matchTk_index);
  const double matchTk_POCAz = L1Tks_POCAz->at(matchTk_index);

  // For-Loop: L1Tks
  for (Size_t iTk = 0; iTk < L1Tks_Pt->size(); iTk++) {
    bool bIsQualityTk = s->SelectTracks(iTk, selectionType);
    if (!bIsQualityTk) continue;
    bool bIsUsedIsoTk = find(isoTks_index.begin(), isoTks_index.end(), iTk) != isoTks_index.end();
    if (bIsUsedIsoTk) continue;

    // Is tk within or isolation annulus
    double tk_eta    = L1Tks_Eta->at(iTk);
    double tk_phi    = L1Tks_Phi->at(iTk);
    double tk_POCAz  = L1Tks_POCAz->at(iTk);
    double deltaR    = auxTools_.DeltaR(matchTk_eta, matchTk_phi, tk_eta, tk_phi);
    bool bIsInsideCone = (deltaR <= deltaR_max) && (deltaR >= deltaR_min);
    if (!bIsInsideCone) continue;

    // Consider only tracks coming from the matching track vertex?
    double deltaPOCAz = abs(matchTk_POCAz - tk_POCAz);
    if (deltaPOCAz < 0.0) continue; // disabled

    // Save the track as an isolation-cone track
    isoTks_index.push_back( int(iTk) );

  } // For-Loop: L1Tks
  
  return;
}


//****************************************************************************
double TauTrigger::GetL1TkTauLdgTkPt(const string tkCollectionType, 
				     const L1TkTauParticle L1TkTau)
//****************************************************************************
{
  
  double ldgTk_Pt    = -1.0;
  vector<int> sigTks = L1TkTau.sigTks_Index_;
  std::sort(sigTks.begin(), sigTks.end(), SortDescendingAbs() );

  if ( tkCollectionType.compare("TTPixelTracks") == 0 ) ldgTk_Pt = L1PixTks_Pt->at( sigTks.at(0) );
  else if ( tkCollectionType.compare("TTTracks") == 0 ) ldgTk_Pt = L1Tks_Pt   ->at( sigTks.at(0) );
  else{
    cout << "E R R O R ! TauTrigger::GetL1TkTauLdgTkPt(...) - Unknown sample \"" << mcSample << "\". EXIT" << endl;
    exit(1);
  }

  return ldgTk_Pt;
}



//****************************************************************************
TLorentzVector TauTrigger::GetL1TkTauSigTksP4(const string tkCollectionType, 
					      const L1TkTauParticle L1TkTau)
//****************************************************************************
{
   
  TLorentzVector sigTks_p4;
  if ( tkCollectionType.compare("TTPixelTracks") == 0 ) sigTks_p4 = s->GetP4FromPixTracks(L1TkTau.sigTks_Index_);
  else if ( tkCollectionType.compare("TTTracks") == 0 ) sigTks_p4 = s->GetP4FromTracks(L1TkTau.sigTks_Index_);
  else{
    cout << "E R R O R ! TauTrigger::GetL1TkTauSigTksP4(...) - Unknown sample \"" << mcSample << "\". EXIT" << endl;
    exit(1);
  }

  return sigTks_p4;
}


//****************************************************************************
TLorentzVector TauTrigger::GetL1TkTauIsoTksP4(const string tkCollectionType, 
					      const L1TkTauParticle L1TkTau)
//****************************************************************************
{
   
  TLorentzVector isoTks_p4;
  if ( tkCollectionType.compare("TTPixelTracks") == 0 ) isoTks_p4 = s->GetP4FromPixTracks(L1TkTau.isoTks_Index_);
  else if ( tkCollectionType.compare("TTTracks") == 0 ) isoTks_p4 = s->GetP4FromTracks(L1TkTau.isoTks_Index_);
  else{
    cout << "E R R O R ! TauTrigger::GetL1TkTauSigTksP4(...) - Unknown sample \"" << mcSample << "\". EXIT" << endl;
    exit(1);
  }

  return isoTks_p4;
}


//****************************************************************************
void TauTrigger::PrintTrackProperties(int tk_Index)
//****************************************************************************
{

  if ( tk_CollectionType.compare("TTPixelTracks") == 0 ) s->PrintPixTrackProperties(tk_Index);
  else if ( tk_CollectionType.compare("TTTracks") == 0 ) s->PrintTrackProperties(tk_Index);
  else
    {
      cout << "E R R O R ! TauTrigger::Loop(...) - Invalid Track Collection Type \"" << tk_CollectionType << "\". EXIT" << endl;
      exit(1);
    }
  
  return;
}


//****************************************************************************
void TauTrigger::PrintL1TkTauProperties(bool bDebug, 
					L1TkTauParticle L1TkTau)
//****************************************************************************
{

  if(!bDebug) return;
  L1TkTau.PrintProperties();
  PrintL1CaloTau(L1TkTau.caloTau_Index_);
  PrintTrackProperties(L1TkTau.matchTk_Index_);

  return;
}


//****************************************************************************
void TauTrigger::GetL1TkTauIsoConePixTracks(const int matchTk_index, 
					    const double deltaR_min,
					    const double deltaR_max,
					    const string selectionType,
					    vector<int> &isoTks_index)
//****************************************************************************
{
  
#ifdef DEBUG
  for (Size_t l = 0; l < isoTks_index.size(); l++)  PrintPixTrackProperties(isoTks_index.at(l));
#endif

  if (matchTk_index == -1) return;
  const double matchTk_eta   = L1PixTks_Eta->at(matchTk_index);
  const double matchTk_phi   = L1PixTks_Phi->at(matchTk_index);
  const double matchTk_POCAz = L1PixTks_POCAz->at(matchTk_index);

  // For-Loop: L1PixTks
  for (Size_t iTk = 0; iTk < L1PixTks_Pt->size(); iTk++) {
    bool bIsQualityTk = s->SelectPixTracks(iTk, selectionType);
    if (!bIsQualityTk) continue;
    bool bIsUsedIsoTk = find(isoTks_index.begin(), isoTks_index.end(), iTk) != isoTks_index.end();
    if (bIsUsedIsoTk) continue;

    // Is tk within or isolation annulus
    double tk_eta      = L1PixTks_Eta->at(iTk);
    double tk_phi      = L1PixTks_Phi->at(iTk);
    double tk_POCAz    = L1PixTks_POCAz->at(iTk);
    double deltaR      = auxTools_.DeltaR(matchTk_eta, matchTk_phi, tk_eta, tk_phi);
    bool bIsInsideCone = (deltaR <= deltaR_max) && (deltaR >= deltaR_min);
    if (!bIsInsideCone) continue;

    // Consider only tracks coming from the matching track vertex?
    double deltaPOCAz = abs(matchTk_POCAz - tk_POCAz);
    if (deltaPOCAz < 0.0) continue; // disabled

    // Save the track as an isolation-cone track
    isoTks_index.push_back( int(iTk) );

  } // For-Loop: L1PixTks
  
  return;
}

#endif

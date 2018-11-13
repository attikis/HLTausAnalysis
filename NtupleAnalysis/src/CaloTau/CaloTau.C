#ifndef CaloTau_cxx
#define CaloTau_cxx

// User
#include "../Auxiliary/interface/constants.h"
#include "CaloTau.h"

// ROOT
#include "TFitResult.h"
#include "TF1.h"
#include "Math/VectorUtil.h"

// C++
#include <stdexcept>

//============================================================================
void CaloTau::InitObjects(void)
//============================================================================
{
  return; 
}


//============================================================================
void CaloTau::InitVars_()
//============================================================================
{

  // Public & Private Variables
  DEBUG = false;
  datasets_  = datasets_.GetDataset(mcSample);
  realTauMom = datasets_.McTauMomPdgId_;
  nMaxNumOfHTausPossible = datasets_.nMcTaus_;
  mcMatching_dRMax  = +0.05;
  dR_match_min      = 0.00;
  dR_match_max      = 0.10;
  dR_sigCone_min    = 0.00;
  dR_sigCone_max    = 0.25; // arbitrarily selected!
  dR_isoCone_min    = 0.00;
  dR_isoCone_max    = 0.435; // 1 x 1 TT = 5 x 5 crystals = 5 x 0.087 = 0.435 (for central region)
  mcMatching_unique = true;
  _eta_C = 0.8;
  _eta_F = 1.6;

  return;
}


//============================================================================
void CaloTau::PrintSettings(void)
//============================================================================
{

  // if (!DEBUG) return;

  // Inform user of settings
  Table settings("Variable | Cut | Value | Default | Units", "Text");
  settings.AddRowColumn(0, "MC Sample");
  settings.AddRowColumn(0, "==" );
  settings.AddRowColumn(0, mcSample );
 
  settings.AddRowColumn(1, "MC-Matching DeltaR");
  settings.AddRowColumn(1, "<=");
  settings.AddRowColumn(1, auxTools_.ToString(mcMatching_dRMax) );
  settings.AddRowColumn(1, "0.05" );
  settings.AddRowColumn(1, "" );

  settings.AddRowColumn(2, "MC-Matching IsUnique");
  settings.AddRowColumn(2, "==");
  settings.AddRowColumn(2, auxTools_.ToString(mcMatching_unique) );
  settings.AddRowColumn(2, "1" );
  settings.AddRowColumn(2, "" );
  
  settings.AddRowColumn(3, "MC-Taus: Mom PdgId");
  settings.AddRowColumn(3, "==");
  settings.AddRowColumn(3, auxTools_.ToString(realTauMom));
  settings.AddRowColumn(3, "N/A" );
  settings.AddRowColumn(3, "" );

  settings.AddRowColumn(3, "MC-Taus: Number Expected");
  settings.AddRowColumn(3, ">=");
  settings.AddRowColumn(3, auxTools_.ToString(nMaxNumOfHTausPossible));
  settings.AddRowColumn(3, "N/A" );
  settings.AddRowColumn(3, "" );

  settings.Print();
  
  return;
}


//============================================================================
void CaloTau::Loop()
//============================================================================
{
  // Sanity check
  if (fChain == 0) return;
  
  const Long64_t nEntries   = (MaxEvents == -1) ? fChain->GetEntries() : min((int)fChain->GetEntries(), MaxEvents);
  
  if (DEBUG) cout << "=== CaloTau:\n\tAnalyzing: " << nEntries << "/" << fChain->GetEntries() << " events" << endl;

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
  
  if (isMinBias) PrintSettings();
  
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
	// GenTausHadronic     = GetHadronicGenTaus(GenTaus, 00.0, 1.479); // CaloTaus currently eta-restricted
	// GenTausTrigger      = GetHadronicGenTaus(GenTaus, 20.0, 1.479); // CaloTaus currently eta-restricted
	GenTausHadronic     = GetHadronicGenTaus(GenTaus, 00.0, 1.3); // CaloTaus currently eta-restricted
	GenTausTrigger      = GetHadronicGenTaus(GenTaus, 20.0, 1.3); // CaloTaus currently eta-restricted
      }

    if (DEBUG)
    {
	cout << "\tPrinting all GenParticle Collections" << endl;
	if (0) PrintGenParticleCollection(GenParticles);
	// PrintGenParticleCollection(GenTaus);
	// PrintGenParticleCollection(GenTausHadronic);    
	PrintGenParticleCollection(GenTausTrigger);
      }

        // Tau Collections
    vector<L1Tau> L1Taus = GetL1Taus(false);
    sort(L1Taus.begin(), L1Taus.end(), [](L1Tau& a, L1Tau& b) {return a.et()  > b.et();});
    vector<L1TkTauParticle> CaloTaus;
    vector<L1TkTauParticle> CaloTausIso;
    vector<L1TkTauParticle> CaloTaus_MC;
    vector<L1TkTauParticle> CaloTausIso_MC;

    // Ensure that all taus are found
    bFoundAllTaus_ = ( (int) GenTausTrigger.size() >= nMaxNumOfHTausPossible);
    if (bFoundAllTaus_) nEvtsWithMaxHTaus++;
    else if (!isMinBias) continue;

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
    // CaloTaus
    // ======================================================================================
    int iCalo = -1;
    bool bFoundMC = false;

    // For-loop: All calo taus
    for (vector<L1Tau>::iterator calo = L1Taus.begin(); calo != L1Taus.end(); calo++)
      {
	iCalo += 1;
	if (0) calo->PrintProperties(iCalo==0);

	// Fill histograms
	hCalo_Et        ->Fill( calo->getEt()    );
	hCalo_Eta       ->Fill( calo->getEta()   );
	hCalo_Phi       ->Fill( calo->getPhi()   );
	hCalo_IEt       ->Fill( calo->getIEt()   );
	hCalo_IEta      ->Fill( calo->getIEta()  );
	hCalo_IPhi      ->Fill( calo->getIPhi()  );
	hCalo_Iso       ->Fill( calo->getIso()   );
	hCalo_TowerIPhi ->Fill( calo->getIPhi()  );
	hCalo_TowerIEta ->Fill( calo->getIEta()  );
	hCalo_RawEt     ->Fill( calo->getRawEt() );
	hCalo_IsoEt     ->Fill( calo->getIsoEt() );
	hCalo_NTT       ->Fill( calo->getNTT()   );
	hCalo_HasEM     ->Fill( calo->getHasEM() );
	hCalo_IsMerged  ->Fill( calo->getIsMerged() );

	// Construct the CaloTk candidate       
        L1TkTauParticle tau(dR_match_min  , dR_match_max,
			    dR_match_min  , dR_match_max,
			    dR_isoCone_min, dR_isoCone_max);
	
	// Assign values to the L1TkTau
	// std::cout << "calo->index() = " << calo->index() << std::endl;
	//tau.SetCaloTau(calo->index());
	tau.SetCaloTau(*calo);
	// tau.SetMatchingTk(matchTk);
	// tau.SetMatchTkDeltaRNew(matchTk_dR);

	// Get the smatching gen-particle
	GetMatchingGenParticle(tau, GenTausTrigger);

	// Save Candidates
	CaloTaus.push_back(tau);
	if ( tau.HasMatchingGenParticle() )
	  {
	  CaloTaus_MC.push_back(tau);
	  bFoundMC = true;
	  }

	// Save Isolated Candidates
	if (!calo->getIso()) continue;
	// if (!calo->getHasEM()) continue;
	// if (calo->getIsMerged() ) continue;
	CaloTausIso.push_back(tau);

	// Fill histograms
	hCaloIso_Et        ->Fill( calo->getEt()    );
	hCaloIso_Eta       ->Fill( calo->getEta()   );
	hCaloIso_Phi       ->Fill( calo->getPhi()   );
	hCaloIso_IEt       ->Fill( calo->getIEt()   );
	hCaloIso_IEta      ->Fill( calo->getIEta()  );
	hCaloIso_IPhi      ->Fill( calo->getIPhi()  );
	hCaloIso_Iso       ->Fill( calo->getIso()   );
	hCaloIso_TowerIPhi ->Fill( calo->getIPhi()  );
	hCaloIso_TowerIEta ->Fill( calo->getIEta()  );
	hCaloIso_RawEt     ->Fill( calo->getRawEt() );
	hCaloIso_IsoEt     ->Fill( calo->getIsoEt() );
	hCaloIso_NTT       ->Fill( calo->getNTT()   );
	hCaloIso_HasEM     ->Fill( calo->getHasEM() );
	hCaloIso_IsMerged  ->Fill( calo->getIsMerged() );

	if ( tau.HasMatchingGenParticle() ) CaloTausIso_MC.push_back(tau);

	// Print?
	if (DEBUG) tau.PrintProperties(false, false, true, true);
      }
    if (bFoundMC) nEvtsMcMatch++;

    // Sort CaloTaus by Et
    sort(CaloTaus.begin(), CaloTaus.end(), [](L1TkTauParticle&a, L1TkTauParticle& b) {return a.GetCaloTau().et()  > b.GetCaloTau().et();});
    sort(CaloTaus_MC.begin(), CaloTaus_MC.end(), [](L1TkTauParticle&a, L1TkTauParticle& b) {return a.GetCaloTau().et()  > b.GetCaloTau().et();});
    sort(CaloTausIso.begin(), CaloTausIso.end(), [](L1TkTauParticle&a, L1TkTauParticle& b) {return a.GetCaloTau().et()  > b.GetCaloTau().et();});
    sort(CaloTausIso_MC.begin(), CaloTausIso_MC.end(), [](L1TkTauParticle&a, L1TkTauParticle& b) {return a.GetCaloTau().et()  > b.GetCaloTau().et();});

    // Fill histos
    hCalo_Multiplicity       ->Fill( CaloTaus.size() );
    hCalo_Multiplicity_MC    ->Fill( CaloTaus_MC.size() );
    hCaloIso_Multiplicity    ->Fill( CaloTausIso.size() );
    hCaloIso_Multiplicity_MC ->Fill( CaloTausIso_MC.size() );
  
    ////////////////////////////////////////////////
    /// CaloTau Properties 
    ////////////////////////////////////////////////
    // for (vector<L1TkTauParticle>::iterator tau = CaloTausIso.begin(); tau != CaloTausIso.end(); tau++)
    for (vector<L1TkTauParticle>::iterator tau = CaloTaus.begin(); tau != CaloTaus.end(); tau++)
      {
	
	if (DEBUG) tau->PrintProperties(true, true, true, true);

	// Do not skip if using MinBias sample as no real taus exist!
	if (!tau->HasMatchingGenParticle() && (isMinBias == false) ) continue;
	
	// Get matching gen particle
	GenParticle p = tau->GetMatchingGenParticle();
	double tauEt  = tau->GetCaloTau().et();
	double tauEta = tau->GetCaloTau().eta();
	double tauPhi = tau->GetCaloTau().phi();
	double visEt  = p.p4vis().Et();
	double visEta = p.p4vis().Eta();
	double visPhi = p.p4vis().Phi();
	double etRes  = (tauEt - visEt)/(visEt);
	double etaRes = (tauEta - visEta)/(visEta);
	double phiRes = (tauPhi - visPhi)/(visPhi);

	// Resolution
	hCalo_ResolutionEt ->Fill( etRes  );
	hCalo_ResolutionEta->Fill( etaRes );
	hCalo_ResolutionPhi->Fill( phiRes );
	if (!tau->GetCaloTau().getIso()) 
	  {
	    hCaloIso_ResolutionEt ->Fill( etRes  );
	    hCaloIso_ResolutionEta->Fill( etaRes );
	    hCaloIso_ResolutionPhi->Fill( phiRes );
	  }

	if (p.finalDaughtersNeutral().size() > 0)
	  {
	    hCalo_ResolutionEt_withNeutrals ->Fill( etRes  );
	    hCalo_ResolutionEta_withNeutrals->Fill( etaRes );
	    hCalo_ResolutionPhi_withNeutrals->Fill( phiRes );
	    if (!tau->GetCaloTau().getIso()) 
	      {
		hCaloIso_ResolutionEt_withNeutrals ->Fill( etRes  );
		hCaloIso_ResolutionEta_withNeutrals->Fill( etaRes );
		hCaloIso_ResolutionPhi_withNeutrals->Fill( phiRes );
	      }
	  }
	else{
	  hCalo_ResolutionEt_noNeutrals ->Fill( etRes  );
	  hCalo_ResolutionEta_noNeutrals->Fill( etaRes );
	  hCalo_ResolutionPhi_noNeutrals->Fill( phiRes );
	    if (!tau->GetCaloTau().getIso()) 
	      {
		hCaloIso_ResolutionEt_noNeutrals ->Fill( etRes  );
		hCaloIso_ResolutionEta_noNeutrals->Fill( etaRes );
		hCaloIso_ResolutionPhi_noNeutrals->Fill( phiRes );
	      }
	}

	if (p.finalDaughtersCharged().size() == 1) 
	  {
	    hCalo_ResolutionEt_1pr ->Fill( etRes  );
	    hCalo_ResolutionEta_1pr->Fill( etaRes );
	    hCalo_ResolutionPhi_1pr->Fill( phiRes );
	    if (!tau->GetCaloTau().getIso()) 
	      {
		hCaloIso_ResolutionEt_1pr ->Fill( etRes  );
		hCaloIso_ResolutionEta_1pr->Fill( etaRes );
		hCaloIso_ResolutionPhi_1pr->Fill( phiRes );
	      }
	  }
	else if (p.finalDaughtersCharged().size() == 3) 
	  {
	    hCalo_ResolutionEt_3pr ->Fill( etRes  );
	    hCalo_ResolutionEta_3pr->Fill( etaRes );
	    hCalo_ResolutionPhi_3pr->Fill( phiRes );
	    if (!tau->GetCaloTau().getIso()) 
	      {
		hCaloIso_ResolutionEt_3pr ->Fill( etRes  );
		hCaloIso_ResolutionEta_3pr->Fill( etaRes );
		hCaloIso_ResolutionPhi_3pr->Fill( phiRes );
	      }
	  }


	if ( IsWithinEtaRegion("Central", tauEta) )
	  {
	    hCalo_ResolutionEt_C ->Fill( etRes  );
	    hCalo_ResolutionEta_C->Fill( etaRes );
	    hCalo_ResolutionPhi_C->Fill( phiRes ); 
	    if (!tau->GetCaloTau().getIso()) 
	      {
		hCaloIso_ResolutionEt_C ->Fill( etRes  );
		hCaloIso_ResolutionEta_C->Fill( etaRes );
		hCaloIso_ResolutionPhi_C->Fill( phiRes ); 
	      }
	  }
	else if ( IsWithinEtaRegion("Intermediate", tauEta) )
	  {
	    hCalo_ResolutionEt_I ->Fill( etRes  );
	    hCalo_ResolutionEta_I->Fill( etaRes );
	    hCalo_ResolutionPhi_I->Fill( phiRes ); 
	    if (!tau->GetCaloTau().getIso()) 
	      {
		hCaloIso_ResolutionEt_I ->Fill( etRes  );
		hCaloIso_ResolutionEta_I->Fill( etaRes );
		hCaloIso_ResolutionPhi_I->Fill( phiRes ); 
	      }
	  }
	// currently no L1Taus in forward eta region
	else if ( IsWithinEtaRegion("Forward", tauEta) )
	  {
	    hCalo_ResolutionEt_F ->Fill( etRes  );
	    hCalo_ResolutionEta_F->Fill( etaRes );
	    hCalo_ResolutionPhi_F->Fill( phiRes ); 
	    if (!tau->GetCaloTau().getIso()) 
	      {
		hCaloIso_ResolutionEt_F ->Fill( etRes  );
		hCaloIso_ResolutionEta_F->Fill( etaRes );
		hCaloIso_ResolutionPhi_F->Fill( phiRes ); 
	      }
	  }
	else{
	  cout << "=== CaloTau::Loop() - Unexpected Eta value of \"" << tauEta << "\". EXIT" << endl;
	  exit(1);
	}
	
	// Calibrations
	hCaloEta_Vs_CaloEtOverVisEt->Fill(visEta, tauEt / visEt );
	if ( tauEt >= 20.0 && tauEt <  40.0)
	  {
	    hCaloEta_Vs_CaloEtOverVisEt_Pt20to40->Fill(visEta, tauEt / visEt );
	    //pCaloEta_Vs_CaloEtOverVisEt_Pt20to40->Fill(visEta, tauEt / visEt, 1);
	  }
	
	if ( tauEt >= 40.0 && tauEt <  60.0) 
	  {
	    hCaloEta_Vs_CaloEtOverVisEt_Pt40to60->Fill(visEta, tauEt / visEt );
	    // pCaloEta_Vs_CaloEtOverVisEt_Pt40to60->Fill(visEta, tauEt / visEt, 1);
	  }
	if ( tauEt >= 60.0 && tauEt <  80.0)
	  {
	    hCaloEta_Vs_CaloEtOverVisEt_Pt60to80  ->Fill(visEta, tauEt / visEt );
	    //pCaloEta_Vs_CaloEtOverVisEt_Pt60to80  ->Fill(visEta, tauEt / visEt, 1);
	  }
	if ( tauEt >= 80.0 && tauEt < 100.0)
	  {
	    hCaloEta_Vs_CaloEtOverVisEt_Pt80to100 ->Fill(visEta, tauEt / visEt );
	    //pCaloEta_Vs_CaloEtOverVisEt_Pt80to100 ->Fill(visEta, tauEt / visEt, 1);	   
	  }
	if ( tauEt >= 100.0 && tauEt < 150.0) 
	  {
	    hCaloEta_Vs_CaloEtOverVisEt_Pt100to150->Fill(visEta, tauEt / visEt );
	    //pCaloEta_Vs_CaloEtOverVisEt_Pt100to150->Fill(visEta, tauEt / visEt, 1);
	  }
	if ( tauEt >= 150.0 && tauEt < 200.0)
	  {
	    hCaloEta_Vs_CaloEtOverVisEt_Pt150to200->Fill(visEta, tauEt / visEt );
	    //pCaloEta_Vs_CaloEtOverVisEt_Pt150to200->Fill(visEta, tauEt / visEt, 1);
	  }

	if (tauEt >= 200.0)
	  {
	    hCaloEta_Vs_CaloEtOverVisEt_PtGE200   ->Fill(visEta, tauEt / visEt );
	    //pCaloEta_Vs_CaloEtOverVisEt_PtGE200   ->Fill(visEta, tauEt / visEt, 1);
	  }

      }// CaloTaus


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

    ////////////////////////////////////////////////
    // Turn-Ons
    ////////////////////////////////////////////////
    FillTurnOn_Numerator_(CaloTaus   , 25.0, hCalo_TurnOn25, hCalo_TurnOn25_1pr, hCalo_TurnOn25_3pr, hCalo_TurnOn25_withNeutrals, hCalo_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(CaloTaus   , 50.0, hCalo_TurnOn50, hCalo_TurnOn50_1pr, hCalo_TurnOn50_3pr, hCalo_TurnOn50_withNeutrals, hCalo_TurnOn50_noNeutrals);
    //
    FillTurnOn_Numerator_(CaloTausIso, 25.0, hCaloIso_TurnOn25, hCaloIso_TurnOn25_1pr, hCaloIso_TurnOn25_3pr, hCaloIso_TurnOn25_withNeutrals, hCaloIso_TurnOn25_noNeutrals);
    FillTurnOn_Numerator_(CaloTausIso, 50.0, hCaloIso_TurnOn25, hCaloIso_TurnOn25_1pr, hCaloIso_TurnOn25_3pr, hCaloIso_TurnOn25_withNeutrals, hCaloIso_TurnOn25_noNeutrals);

    ////////////////////////////////////////////////
    // SingleTau
    ////////////////////////////////////////////////
    FillSingleTau_(CaloTaus, hCalo_Rate  , hCalo_Eff  );
    FillSingleTau_(CaloTaus, hCalo_Rate_C, hCalo_Eff_C, 0.0, 1.0);
    FillSingleTau_(CaloTaus, hCalo_Rate_I, hCalo_Eff_I, 1.0, 1.6);
    FillSingleTau_(CaloTaus, hCalo_Rate_F, hCalo_Eff_F, 1.6, 3.0); // 2.5 is max
    //
    FillSingleTau_(CaloTausIso, hCaloIso_Rate  , hCaloIso_Eff  );
    FillSingleTau_(CaloTausIso, hCaloIso_Rate_C, hCaloIso_Eff_C, 0.0, 1.0);
    FillSingleTau_(CaloTausIso, hCaloIso_Rate_I, hCaloIso_Eff_I, 1.0, 1.6);
    FillSingleTau_(CaloTausIso, hCaloIso_Rate_F, hCaloIso_Eff_F, 1.6, 3.0); // 2.5 is max

    ////////////////////////////////////////////////
    // DiTau
    ////////////////////////////////////////////////
    FillDiTau_(CaloTaus, hCalo_Rate_DiTau  , hCalo_Eff_DiTau);
    FillDiTau_(CaloTaus, hCalo_Rate_DiTau_C, hCalo_Eff_DiTau_C, 0.0, 1.0);
    FillDiTau_(CaloTaus, hCalo_Rate_DiTau_I, hCalo_Eff_DiTau_I, 1.0, 1.6);
    FillDiTau_(CaloTaus, hCalo_Rate_DiTau_F, hCalo_Eff_DiTau_F, 1.6, 3.0);
    //
    FillDiTau_(CaloTausIso, hCaloIso_Rate_DiTau  , hCaloIso_Eff_DiTau);
    FillDiTau_(CaloTausIso, hCaloIso_Rate_DiTau_C, hCaloIso_Eff_DiTau_C, 0.0, 1.0);
    FillDiTau_(CaloTausIso, hCaloIso_Rate_DiTau_I, hCaloIso_Eff_DiTau_I, 1.0, 1.6);
    FillDiTau_(CaloTausIso, hCaloIso_Rate_DiTau_F, hCaloIso_Eff_DiTau_F, 1.6, 3.0);

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
  double N = nEntries;
  if (isMinBias) // new: speed things up a bit (Rate plots only make sense for Neutrino gun!)
    {
      // SingleTau
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate  , N);
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate_C, N);
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate_I, N);
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate_F, N);
      //
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate  , N);
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate_C, N);
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate_I, N);
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate_F, N);

      // DiTau
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate_DiTau  , N);
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate_DiTau_C, N);
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate_DiTau_I, N);
      histoTools_.ConvertToRateHisto_1D(hCalo_Rate_DiTau_F, N);
      //
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate_DiTau  , N);
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate_DiTau_C, N);
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate_DiTau_I, N);
      histoTools_.ConvertToRateHisto_1D(hCaloIso_Rate_DiTau_F, N);

    }
  else // new: speed things up a bit
    {
      // Single Tau
      FinaliseEffHisto_( hCalo_Eff  , nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCalo_Eff_C, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCalo_Eff_I, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCalo_Eff_F, nEvtsWithMaxHTaus);
      //
      FinaliseEffHisto_( hCaloIso_Eff  , nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCaloIso_Eff_C, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCaloIso_Eff_I, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCaloIso_Eff_F, nEvtsWithMaxHTaus);

      // DiTau
      FinaliseEffHisto_( hCalo_Eff_DiTau  , nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCalo_Eff_DiTau_C, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCalo_Eff_DiTau_I, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCalo_Eff_DiTau_F, nEvtsWithMaxHTaus);
      //
      FinaliseEffHisto_( hCaloIso_Eff_DiTau  , nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCaloIso_Eff_DiTau_C, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCaloIso_Eff_DiTau_I, nEvtsWithMaxHTaus);
      FinaliseEffHisto_( hCaloIso_Eff_DiTau_F, nEvtsWithMaxHTaus);

    }
      
  ////////////////////////////////////////////////
  // Turn-Ons
  ////////////////////////////////////////////////
  histoTools_.DivideHistos_1D(hCalo_TurnOn25             , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hCalo_TurnOn25_1pr         , hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hCalo_TurnOn25_3pr         , hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hCalo_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hCalo_TurnOn25_noNeutrals  , hMcHadronicTau_VisEt_noNeutrals);
  histoTools_.DivideHistos_1D(hCalo_TurnOn50             , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hCalo_TurnOn50_1pr         , hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hCalo_TurnOn50_3pr         , hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hCalo_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hCalo_TurnOn50_noNeutrals  , hMcHadronicTau_VisEt_noNeutrals);
  //
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn25             , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn25_1pr         , hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn25_3pr         , hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn25_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn25_noNeutrals  , hMcHadronicTau_VisEt_noNeutrals);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn50             , hMcHadronicTau_VisEt);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn50_1pr         , hMcHadronicTau_VisEt_1pr);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn50_3pr         , hMcHadronicTau_VisEt_3pr);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn50_withNeutrals, hMcHadronicTau_VisEt_withNeutrals);
  histoTools_.DivideHistos_1D(hCaloIso_TurnOn50_noNeutrals  , hMcHadronicTau_VisEt_noNeutrals);

  ////////////////////////////////////////////////
  // Profile Histograms
  ////////////////////////////////////////////////
  pCaloEta_Vs_CaloEtOverVisEt            = hCaloEta_Vs_CaloEtOverVisEt           ->ProfileX("pCaloEta_Vs_CaloEtOverVisEt");
  pCaloEta_Vs_CaloEtOverVisEt_Pt20to40   = hCaloEta_Vs_CaloEtOverVisEt_Pt20to40  ->ProfileX("pCaloEta_Vs_CaloEtOverVisEt_Pt20to40");
  pCaloEta_Vs_CaloEtOverVisEt_Pt40to60   = hCaloEta_Vs_CaloEtOverVisEt_Pt40to60  ->ProfileX("pCaloEta_Vs_CaloEtOverVisEt_Pt40to60");
  pCaloEta_Vs_CaloEtOverVisEt_Pt60to80   = hCaloEta_Vs_CaloEtOverVisEt_Pt60to80  ->ProfileX("pCaloEta_Vs_CaloEtOverVisEt_Pt60to80");
  pCaloEta_Vs_CaloEtOverVisEt_Pt80to100  = hCaloEta_Vs_CaloEtOverVisEt_Pt80to100 ->ProfileX("pCaloEta_Vs_CaloEtOverVisEt_Pt80to100");
  pCaloEta_Vs_CaloEtOverVisEt_Pt100to150 = hCaloEta_Vs_CaloEtOverVisEt_Pt100to150->ProfileX("pCaloEta_Vs_CaloEtOverVisEt_Pt100to150");
  pCaloEta_Vs_CaloEtOverVisEt_Pt150to200 = hCaloEta_Vs_CaloEtOverVisEt_Pt150to200->ProfileX("pCaloEta_Vs_CaloEtOverVisEt_Pt150to200");
  pCaloEta_Vs_CaloEtOverVisEt_PtGE200    = hCaloEta_Vs_CaloEtOverVisEt_PtGE200   ->ProfileX("pCaloEta_Vs_CaloEtOverVisEt_PtGE200");

  ////////////////////////////////////////////////
  // Write the histograms to the file
  ////////////////////////////////////////////////
  WriteHistos_();
  if (DEBUG) auxTools_.StopwatchStop(5, "minutes", "Total Time");
  
}


//============================================================================
void CaloTau::BookHistos_(void)
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

  // Calibration (Binning according to EB granularity at L1: i.e. 17 iEta for range 0.0 < eta < 1.479
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt           , "CaloEta_Vs_CaloEtOverVisEt"           , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt20to40  , "CaloEta_Vs_CaloEtOverVisEt_Pt20to40"  , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt40to60  , "CaloEta_Vs_CaloEtOverVisEt_Pt40to60"  , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt60to80  , "CaloEta_Vs_CaloEtOverVisEt_Pt60to80"  , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt80to100 , "CaloEta_Vs_CaloEtOverVisEt_Pt80to100" , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt100to150, "CaloEta_Vs_CaloEtOverVisEt_Pt100to150", ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt150to200, "CaloEta_Vs_CaloEtOverVisEt_Pt150to200", ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_PtGE200   , "CaloEta_Vs_CaloEtOverVisEt_PtGE200"   , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", 17*2,  -1.479, +1.479, 100, 0.0, +5.0);
  //
  // (const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, Option_t* option = "")
  pCaloEta_Vs_CaloEtOverVisEt            = new TProfile("pCaloEta_Vs_CaloEtOverVisEt"           , ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  pCaloEta_Vs_CaloEtOverVisEt_Pt20to40   = new TProfile("pCaloEta_Vs_CaloEtOverVisEt_Pt20to40"  , ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  pCaloEta_Vs_CaloEtOverVisEt_Pt40to60   = new TProfile("pCaloEta_Vs_CaloEtOverVisEt_Pt40to60"  , ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  pCaloEta_Vs_CaloEtOverVisEt_Pt60to80   = new TProfile("pCaloEta_Vs_CaloEtOverVisEt_Pt60to80"  , ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  pCaloEta_Vs_CaloEtOverVisEt_Pt80to100  = new TProfile("pCaloEta_Vs_CaloEtOverVisEt_Pt80to100" , ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  pCaloEta_Vs_CaloEtOverVisEt_Pt100to150 = new TProfile("pCaloEta_Vs_CaloEtOverVisEt_Pt100to150", ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  pCaloEta_Vs_CaloEtOverVisEt_Pt150to200 = new TProfile("pCaloEta_Vs_CaloEtOverVisEt_Pt150to200", ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  pCaloEta_Vs_CaloEtOverVisEt_PtGE200    = new TProfile("pCaloEta_Vs_CaloEtOverVisEt_PtGE200"   , ";<#eta^{L1}>; <E_{T}^{L1} / E_{T}^{vis}>", 17*2,  -1.479, +1.479, 0.0, +5.0, "");
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt             , "hCaloEta_Vs_CaloEtOverVisEt"           , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt20to40    , "hCaloEta_Vs_CaloEtOverVisEt_Pt20to40"  , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt40to60    , "hCaloEta_Vs_CaloEtOverVisEt_Pt40to60"  , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt60to80    , "hCaloEta_Vs_CaloEtOverVisEt_Pt60to80"  , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_Pt80to100   , "hCaloEta_Vs_CaloEtOverVisEt_Pt80to100" , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_PtGE100to150, "hCaloEta_Vs_CaloEtOverVisEt_Pt100to150", ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_PtGE150to200, "hCaloEta_Vs_CaloEtOverVisEt_Pt150to200", ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);
  // histoTools_.BookHisto_2D(hCaloEta_Vs_CaloEtOverVisEt_PtGE200     , "hCaloEta_Vs_CaloEtOverVisEt_PtGE200"   , ";#eta^{L1}; E_{T}^{L1} / E_{T}^{vis}", nEta,  minEta, maxEta, 100, 0.0, +5.0);

  // Counters
  histoTools_.BookHisto_1D(hCounters, "Counters",  "", 15, 0.0, +15.0);

  // CaloTaus
  histoTools_.BookHisto_1D(hCalo_Multiplicity   , "Calo_Multiplicity"   , tN  ,  nN  , minN   , maxN    );
  histoTools_.BookHisto_1D(hCalo_Multiplicity_MC, "Calo_Multiplicity_MC", tN  ,  nN  , minN   , maxN    );
  histoTools_.BookHisto_1D(hCalo_Et             , "Calo_Et"             , tEt ,  nEt , minEt  , maxEt   );
  histoTools_.BookHisto_1D(hCalo_Eta            , "Calo_Eta"            , tEta,  nEta, minEta , maxEta  );
  histoTools_.BookHisto_1D(hCalo_Phi            , "Calo_Phi"            , tPhi,  nPhi, minPhi , maxPhi  );
  histoTools_.BookHisto_1D(hCalo_IEt            , "Calo_IEt"            , tEt ,  nEt , minEt  , maxEt   );
  histoTools_.BookHisto_1D(hCalo_IEta           , "Calo_IEta"           , tEta, nIEta, minIEta, maxIEta );
  histoTools_.BookHisto_1D(hCalo_IPhi           , "Calo_IPhi"           , tPhi, nIPhi, minIPhi, maxIPhi );
  histoTools_.BookHisto_1D(hCalo_Iso            , "Calo_Iso"            ,   "", nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hCalo_TowerIEta      , "Calo_TowerIEta"      , tEta, nIEta, minIEta, maxIEta );
  histoTools_.BookHisto_1D(hCalo_TowerIPhi      , "Calo_TowerIPhi"      , tPhi, nIPhi, minIPhi, maxIPhi );
  histoTools_.BookHisto_1D(hCalo_RawEt          , "Calo_RawEt"          , tEt , nEt*3, minEt  , maxEt*3 );
  histoTools_.BookHisto_1D(hCalo_IsoEt          , "Calo_IsoEt"          , tEt , nEt  , minEt  , maxEt   );
  histoTools_.BookHisto_1D(hCalo_NTT            , "Calo_NTT"            , tN  ,nN*100, minN   , maxN*100);
  histoTools_.BookHisto_1D(hCalo_HasEM          , "Calo_HasEM"          , ""  , nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hCalo_IsMerged       , "Calo_IsMerged"       , ""  , nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hCalo_DeltaRGenP     , "Calo_DeltaRGenP"     , tDR , nR   , minR   , maxR    );
  //
  histoTools_.BookHisto_1D(hCaloIso_Multiplicity   , "CaloIso_Multiplicity"   , tN  ,  nN  , minN   , maxN    );
  histoTools_.BookHisto_1D(hCaloIso_Multiplicity_MC, "CaloIso_Multiplicity_MC", tN  ,  nN  , minN   , maxN    );
  histoTools_.BookHisto_1D(hCaloIso_Et             , "CaloIso_Et"             , tEt ,  nEt , minEt  , maxEt   );
  histoTools_.BookHisto_1D(hCaloIso_Eta            , "CaloIso_Eta"            , tEta,  nEta, minEta , maxEta  );
  histoTools_.BookHisto_1D(hCaloIso_Phi            , "CaloIso_Phi"            , tPhi,  nPhi, minPhi , maxPhi  );
  histoTools_.BookHisto_1D(hCaloIso_IEt            , "CaloIso_IEt"            , tEt ,  nEt , minEt  , maxEt   );
  histoTools_.BookHisto_1D(hCaloIso_IEta           , "CaloIso_IEta"           , tEta, nIEta, minIEta, maxIEta );
  histoTools_.BookHisto_1D(hCaloIso_IPhi           , "CaloIso_IPhi"           , tPhi, nIPhi, minIPhi, maxIPhi );
  histoTools_.BookHisto_1D(hCaloIso_Iso            , "CaloIso_Iso"            ,   "", nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hCaloIso_TowerIEta      , "CaloIso_TowerIEta"      , tEta, nIEta, minIEta, maxIEta );
  histoTools_.BookHisto_1D(hCaloIso_TowerIPhi      , "CaloIso_TowerIPhi"      , tPhi, nIPhi, minIPhi, maxIPhi );
  histoTools_.BookHisto_1D(hCaloIso_RawEt          , "CaloIso_RawEt"          , tEt , nEt*3, minEt  , maxEt*3 );
  histoTools_.BookHisto_1D(hCaloIso_IsoEt          , "CaloIso_IsoEt"          , tEt , nEt  , minEt  , maxEt   );
  histoTools_.BookHisto_1D(hCaloIso_NTT            , "CaloIso_NTT"            , tN  ,nN*100, minN   , maxN*100);
  histoTools_.BookHisto_1D(hCaloIso_HasEM          , "CaloIso_HasEM"          , ""  , nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hCaloIso_IsMerged       , "CaloIso_IsMerged"       , ""  , nBool, minBool, maxBool );
  histoTools_.BookHisto_1D(hCaloIso_DeltaRGenP     , "CaloIso_DeltaRGenP"     , tDR , nR   , minR   , maxR    );

  // Resolutions 
  histoTools_.BookHisto_1D(hCalo_ResolutionEt              , "Calo_ResolutionEt"              , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEt_1pr          , "Calo_ResolutionEt_1pr"          , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEt_3pr          , "Calo_ResolutionEt_3pr"          , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEt_withNeutrals , "Calo_ResolutionEt_withNeutrals" , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEt_noNeutrals   , "Calo_ResolutionEt_noNeutrals"   , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta             , "Calo_ResolutionEta"             , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta_1pr         , "Calo_ResolutionEta_1pr"         , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta_3pr         , "Calo_ResolutionEta_3pr"         , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta_withNeutrals, "Calo_ResolutionEta_withNeutrals", ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta_noNeutrals  , "Calo_ResolutionEta_noNeutrals"  , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi             , "Calo_ResolutionPhi"             , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi_1pr         , "Calo_ResolutionPhi_1pr"         , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi_3pr         , "Calo_ResolutionPhi_3pr"         , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi_withNeutrals, "Calo_ResolutionPhi_withNeutrals", ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi_noNeutrals  , "Calo_ResolutionPhi_noNeutrals"  , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEt_C            , "Calo_ResolutionEt_C"            , ";E_{T} (GeV);Events / %.0f GeV"               , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEt_I            , "Calo_ResolutionEt_I"            , ";E_{T} (GeV);Events / %.0f GeV"               , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEt_F            , "Calo_ResolutionEt_F"            , ";E_{T} (GeV);Events / %.0f GeV"               , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta_C           , "Calo_ResolutionEta_C"           , ";#eta;Events / %.2f"                          , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta_I           , "Calo_ResolutionEta_I"           , ";#eta;Events / %.2f"                          , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionEta_F           , "Calo_ResolutionEta_F"           , ";#eta;Events / %.2f"                          , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi_C           , "Calo_ResolutionPhi_C"           , ";#phi (rads);Events / %.2f rads"              , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi_I           , "Calo_ResolutionPhi_I"           , ";#phi (rads);Events / %.2f rads"              , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCalo_ResolutionPhi_F           , "Calo_ResolutionPhi_F"           , ";#phi (rads);Events / %.2f rads"              , 2000, -10.0, +10.0);
  //
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt              , "CaloIso_ResolutionEt"              , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_1pr          , "CaloIso_ResolutionEt_1pr"          , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_3pr          , "CaloIso_ResolutionEt_3pr"          , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_withNeutrals , "CaloIso_ResolutionEt_withNeutrals" , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_noNeutrals   , "CaloIso_ResolutionEt_noNeutrals"   , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta             , "CaloIso_ResolutionEta"             , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_1pr         , "CaloIso_ResolutionEta_1pr"         , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_3pr         , "CaloIso_ResolutionEta_3pr"         , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_withNeutrals, "CaloIso_ResolutionEta_withNeutrals", ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_noNeutrals  , "CaloIso_ResolutionEta_noNeutrals"  , ";#delta #eta / #eta^{vis};Events / %.2f"      , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi             , "CaloIso_ResolutionPhi"             , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_1pr         , "CaloIso_ResolutionPhi_1pr"         , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_3pr         , "CaloIso_ResolutionPhi_3pr"         , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_withNeutrals, "CaloIso_ResolutionPhi_withNeutrals", ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_noNeutrals  , "CaloIso_ResolutionPhi_noNeutrals"  , ";#delta #phi / #phi^{vis};Events / %.2f rads" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_C            , "CaloIso_ResolutionEt_C"            , ";E_{T} (GeV);Events / %.0f GeV"               , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_I            , "CaloIso_ResolutionEt_I"            , ";E_{T} (GeV);Events / %.0f GeV"               , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_F            , "CaloIso_ResolutionEt_F"            , ";E_{T} (GeV);Events / %.0f GeV"               , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_C           , "CaloIso_ResolutionEta_C"           , ";#eta;Events / %.2f"                          , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_I           , "CaloIso_ResolutionEta_I"           , ";#eta;Events / %.2f"                          , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_F           , "CaloIso_ResolutionEta_F"           , ";#eta;Events / %.2f"                          , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_C           , "CaloIso_ResolutionPhi_C"           , ";#phi (rads);Events / %.2f rads"              , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_I           , "CaloIso_ResolutionPhi_I"           , ";#phi (rads);Events / %.2f rads"              , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_F           , "CaloIso_ResolutionPhi_F"           , ";#phi (rads);Events / %.2f rads"              , 2000, -10.0, +10.0);

  // Resolutions (Iso)
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt    , "CaloIso_ResolutionEt"     , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_1pr, "CaloIso_ResolutionEt_1pr" , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_3pr, "CaloIso_ResolutionEt_3pr" , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_withNeutrals, "CaloIso_ResolutionEt_withNeutrals", ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_noNeutrals  , "CaloIso_ResolutionEt_noNeutrals"  , ";#delta E_{T} / E_{T}^{vis};Events / %.0f GeV", 2000, -10.0, +10.0);

  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta    , "CaloIso_ResolutionEta"    , ";#delta #eta / #eta^{vis};Events / %.2f", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_1pr, "CaloIso_ResolutionEta_1pr", ";#delta #eta / #eta^{vis};Events / %.2f", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_3pr, "CaloIso_ResolutionEta_3pr", ";#delta #eta / #eta^{vis};Events / %.2f", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_withNeutrals, "CaloIso_ResolutionEta_withNeutrals", ";#delta #eta / #eta^{vis};Events / %.2f", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_noNeutrals  , "CaloIso_ResolutionEta_noNeutrals"  , ";#delta #eta / #eta^{vis};Events / %.2f", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi    , "CaloIso_ResolutionPhi"    , ";#delta #phi / #phi^{vis};Events / %.2f rads", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_1pr, "CaloIso_ResolutionPhi_1pr", ";#delta #phi / #phi^{vis};Events / %.2f rads", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_3pr, "CaloIso_ResolutionPhi_3pr", ";#delta #phi / #phi^{vis};Events / %.2f rads", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_withNeutrals, "CaloIso_ResolutionPhi_withNeutrals", ";#phi (rads);Events / %.2f rads", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_noNeutrals  , "CaloIso_ResolutionPhi_noNeutrals"  , ";#phi (rads);Events / %.2f rads", 2000, -1.0, +1.0);

  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_C , "CaloIso_ResolutionEt_C" , ";E_{T} (GeV);Events / %.0f GeV" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_C, "CaloIso_ResolutionEta_C", ";#eta;Events / %.2f"            , 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_C, "CaloIso_ResolutionPhi_C", ";#phi (rads);Events / %.2f rads", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_I , "CaloIso_ResolutionEt_I" , ";E_{T} (GeV);Events / %.0f GeV" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_I, "CaloIso_ResolutionEta_I", ";#eta;Events / %.2f"            , 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_I, "CaloIso_ResolutionPhi_I", ";#phi (rads);Events / %.2f rads", 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEt_F , "CaloIso_ResolutionEt_F" , ";E_{T} (GeV);Events / %.0f GeV" , 2000, -10.0, +10.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionEta_F, "CaloIso_ResolutionEta_F", ";#eta;Events / %.2f"            , 2000, -1.0, +1.0);
  histoTools_.BookHisto_1D(hCaloIso_ResolutionPhi_F, "CaloIso_ResolutionPhi_F", ";#phi (rads);Events / %.2f rads", 2000, -1.0, +1.0);

  // SingleTau
  histoTools_.BookHisto_1D(hCalo_Rate    , "Calo_Rate"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_C  , "Calo_Rate_C"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_I  , "Calo_Rate_I"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_F  , "Calo_Rate_F"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff     , "Calo_Eff"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_C   , "Calo_Eff_C"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_I   , "Calo_Eff_I"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_F   , "Calo_Eff_F"   , "", nEt , minEt , maxEt );
  //
  histoTools_.BookHisto_1D(hCaloIso_Rate    , "CaloIso_Rate"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Rate_C  , "CaloIso_Rate_C"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Rate_I  , "CaloIso_Rate_I"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Rate_F  , "CaloIso_Rate_F"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff     , "CaloIso_Eff"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff_C   , "CaloIso_Eff_C"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff_I   , "CaloIso_Eff_I"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff_F   , "CaloIso_Eff_F"   , "", nEt , minEt , maxEt );

  // DiTau
  histoTools_.BookHisto_1D(hCalo_Rate_DiTau    , "Calo_Rate_DiTau"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_DiTau_C  , "Calo_Rate_DiTau_C" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_DiTau_I  , "Calo_Rate_DiTau_I" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Rate_DiTau_F  , "Calo_Rate_DiTau_F" , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_DiTau     , "Calo_Eff_DiTau"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_DiTau_C   , "Calo_Eff_DiTau_C"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_DiTau_I   , "Calo_Eff_DiTau_I"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_Eff_DiTau_F   , "Calo_Eff_DiTau_F"  , "", nEt , minEt , maxEt );
  //
  histoTools_.BookHisto_1D(hCaloIso_Rate_DiTau  , "CaloIso_Rate_DiTau"    , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Rate_DiTau_C, "CaloIso_Rate_DiTau_C"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Rate_DiTau_I, "CaloIso_Rate_DiTau_I"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Rate_DiTau_F, "CaloIso_Rate_DiTau_F"  , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff_DiTau   , "CaloIso_Eff_DiTau"     , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff_DiTau_C , "CaloIso_Eff_DiTau_C"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff_DiTau_I , "CaloIso_Eff_DiTau_I"   , "", nEt , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_Eff_DiTau_F , "CaloIso_Eff_DiTau_F"   , "", nEt , minEt , maxEt );

  // Turn-Ons
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt             , "McHadronicTau_VisEt"             , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_1pr         , "McHadronicTau_VisEt_1pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_3pr         , "McHadronicTau_VisEt_3pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_withNeutrals, "McHadronicTau_VisEt_withNeutrals", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hMcHadronicTau_VisEt_noNeutrals  , "McHadronicTau_VisEt_noNeutrals"  , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn25             , "Calo_TurnOn25"             , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn25_1pr         , "Calo_TurnOn25_1pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn25_3pr         , "Calo_TurnOn25_3pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn25_withNeutrals, "Calo_TurnOn25_withNeutrals", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn25_noNeutrals  , "Calo_TurnOn25_noNeutrals"  , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn50             , "Calo_TurnOn50"             , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn50_1pr         , "Calo_TurnOn50_1pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn50_3pr         , "Calo_TurnOn50_3pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn50_withNeutrals, "Calo_TurnOn50_withNeutrals", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCalo_TurnOn50_noNeutrals  , "Calo_TurnOn50_noNeutrals"  , "", 60 , minEt , maxEt );
  //
  histoTools_.BookHisto_1D(hCaloIso_TurnOn25             , "CaloIso_TurnOn25"             , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn25_1pr         , "CaloIso_TurnOn25_1pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn25_3pr         , "CaloIso_TurnOn25_3pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn25_withNeutrals, "CaloIso_TurnOn25_withNeutrals", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn25_noNeutrals  , "CaloIso_TurnOn25_noNeutrals"  , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn50             , "CaloIso_TurnOn50"             , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn50_1pr         , "CaloIso_TurnOn50_1pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn50_3pr         , "CaloIso_TurnOn50_3pr"         , "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn50_withNeutrals, "CaloIso_TurnOn50_withNeutrals", "", 60 , minEt , maxEt );
  histoTools_.BookHisto_1D(hCaloIso_TurnOn50_noNeutrals, "CaloIso_TurnOn50_noNeutrals"    , "", 60 , minEt , maxEt );

  return;
}

//============================================================================
void CaloTau::WriteHistos_(void)
//============================================================================
{
  // Location -> outFile 
  outFile->cd();
  
  // GenParticles Histograms
  hGenP_VisEt_Vs_dRMaxLdgPion->Write();
  hGenP_PtLdg_Vs_dRMaxLdgPion->Write();

  // Calibration
  hCaloEta_Vs_CaloEtOverVisEt->Write();
  hCaloEta_Vs_CaloEtOverVisEt_Pt20to40  ->Write();
  hCaloEta_Vs_CaloEtOverVisEt_Pt40to60  ->Write();
  hCaloEta_Vs_CaloEtOverVisEt_Pt60to80  ->Write();
  hCaloEta_Vs_CaloEtOverVisEt_Pt80to100 ->Write();
  hCaloEta_Vs_CaloEtOverVisEt_Pt100to150->Write();
  hCaloEta_Vs_CaloEtOverVisEt_Pt150to200->Write();
  hCaloEta_Vs_CaloEtOverVisEt_PtGE200  ->Write();
  //
  pCaloEta_Vs_CaloEtOverVisEt->Write();
  pCaloEta_Vs_CaloEtOverVisEt_Pt20to40  ->Write();
  pCaloEta_Vs_CaloEtOverVisEt_Pt40to60  ->Write();
  pCaloEta_Vs_CaloEtOverVisEt_Pt60to80  ->Write();
  pCaloEta_Vs_CaloEtOverVisEt_Pt80to100 ->Write();
  pCaloEta_Vs_CaloEtOverVisEt_Pt100to150->Write();
  pCaloEta_Vs_CaloEtOverVisEt_Pt150to200->Write();
  pCaloEta_Vs_CaloEtOverVisEt_PtGE200  ->Write();

  // Counters
  hCounters->Write();

  // Calos
  hCalo_Multiplicity->Write();
  hCalo_Multiplicity_MC->Write();
  hCalo_Et->Write();
  hCalo_Eta->Write();
  hCalo_Phi->Write();
  hCalo_IEt->Write();
  hCalo_IEta->Write();
  hCalo_IPhi->Write();
  hCalo_Iso->Write();
  hCalo_TowerIEta->Write();
  hCalo_TowerIPhi->Write();
  hCalo_RawEt->Write();
  hCalo_IsoEt->Write();
  hCalo_NTT->Write();
  hCalo_HasEM->Write();
  hCalo_IsMerged->Write();
  //
  hCaloIso_Multiplicity->Write();
  hCaloIso_Multiplicity_MC->Write();
  hCaloIso_Et->Write();
  hCaloIso_Eta->Write();
  hCaloIso_Phi->Write();
  hCaloIso_IEt->Write();
  hCaloIso_IEta->Write();
  hCaloIso_IPhi->Write();
  hCaloIso_Iso->Write();
  hCaloIso_TowerIEta->Write();
  hCaloIso_TowerIPhi->Write();
  hCaloIso_RawEt->Write();
  hCaloIso_IsoEt->Write();
  hCaloIso_NTT->Write();
  hCaloIso_HasEM->Write();
  hCaloIso_IsMerged->Write();

  // Resolutions
  hCalo_ResolutionEt->Write();
  hCalo_ResolutionEt_1pr->Write();
  hCalo_ResolutionEt_3pr->Write();
  hCalo_ResolutionEt_withNeutrals->Write();
  hCalo_ResolutionEt_noNeutrals->Write();
  hCalo_ResolutionEta->Write();
  hCalo_ResolutionEta_1pr->Write();
  hCalo_ResolutionEta_3pr->Write();
  hCalo_ResolutionEta_withNeutrals->Write();
  hCalo_ResolutionEta_noNeutrals->Write();
  hCalo_ResolutionPhi->Write();
  hCalo_ResolutionPhi_1pr->Write();
  hCalo_ResolutionPhi_3pr->Write();
  hCalo_ResolutionPhi_withNeutrals->Write();
  hCalo_ResolutionPhi_noNeutrals->Write();
  hCalo_ResolutionEt_C->Write();
  hCalo_ResolutionEta_C->Write();
  hCalo_ResolutionPhi_C->Write();
  hCalo_ResolutionEt_I->Write();
  hCalo_ResolutionEta_I->Write();
  hCalo_ResolutionPhi_I->Write();
  hCalo_ResolutionEt_F->Write();
  hCalo_ResolutionEta_F->Write();
  hCalo_ResolutionPhi_F->Write();
  //
  hCaloIso_ResolutionEt->Write();
  hCaloIso_ResolutionEt_1pr->Write();
  hCaloIso_ResolutionEt_3pr->Write();
  hCaloIso_ResolutionEt_withNeutrals->Write();
  hCaloIso_ResolutionEt_noNeutrals->Write();
  hCaloIso_ResolutionEta->Write();
  hCaloIso_ResolutionEta_1pr->Write();
  hCaloIso_ResolutionEta_3pr->Write();
  hCaloIso_ResolutionEta_withNeutrals->Write();
  hCaloIso_ResolutionEta_noNeutrals->Write();
  hCaloIso_ResolutionPhi->Write();
  hCaloIso_ResolutionPhi_1pr->Write();
  hCaloIso_ResolutionPhi_3pr->Write();
  hCaloIso_ResolutionPhi_withNeutrals->Write();
  hCaloIso_ResolutionPhi_noNeutrals->Write();
  hCaloIso_ResolutionEt_C->Write();
  hCaloIso_ResolutionEta_C->Write();
  hCaloIso_ResolutionPhi_C->Write();
  hCaloIso_ResolutionEt_I->Write();
  hCaloIso_ResolutionEta_I->Write();
  hCaloIso_ResolutionPhi_I->Write();
  hCaloIso_ResolutionEt_F->Write();
  hCaloIso_ResolutionEta_F->Write();
  hCaloIso_ResolutionPhi_F->Write();

  // SingleTau
  hCalo_Rate->Write();
  hCalo_Rate_C->Write();
  hCalo_Rate_I->Write();
  hCalo_Rate_F->Write();
  hCalo_Eff->Write();
  hCalo_Eff_C->Write();
  hCalo_Eff_I->Write();
  hCalo_Eff_F->Write();
  //
  hCaloIso_Rate->Write();
  hCaloIso_Rate_C->Write();
  hCaloIso_Rate_I->Write();
  hCaloIso_Rate_F->Write();
  hCaloIso_Eff->Write();
  hCaloIso_Eff_C->Write();
  hCaloIso_Eff_I->Write();
  hCaloIso_Eff_F->Write();

  // DiTau
  hCalo_Rate_DiTau->Write();
  hCalo_Rate_DiTau_C->Write();
  hCalo_Rate_DiTau_I->Write();
  hCalo_Rate_DiTau_F->Write();
  hCalo_Eff_DiTau->Write();
  hCalo_Eff_DiTau_C->Write();
  hCalo_Eff_DiTau_I->Write();
  hCalo_Eff_DiTau_F->Write();
  //
  hCaloIso_Rate_DiTau->Write();
  hCaloIso_Rate_DiTau_C->Write();
  hCaloIso_Rate_DiTau_I->Write();
  hCaloIso_Rate_DiTau_F->Write();
  hCaloIso_Eff_DiTau->Write();
  hCaloIso_Eff_DiTau_C->Write();
  hCaloIso_Eff_DiTau_I->Write();
  hCaloIso_Eff_DiTau_F->Write();

  // Turn-Ons
  hMcHadronicTau_VisEt->Write();
  hMcHadronicTau_VisEt_1pr->Write();
  hMcHadronicTau_VisEt_3pr->Write();
  hMcHadronicTau_VisEt_withNeutrals->Write();
  hMcHadronicTau_VisEt_noNeutrals->Write();
  hCalo_TurnOn25->Write();
  hCalo_TurnOn25_1pr->Write();
  hCalo_TurnOn25_3pr->Write();
  hCalo_TurnOn25_withNeutrals->Write();
  hCalo_TurnOn25_noNeutrals->Write();
  hCalo_TurnOn50->Write();
  hCalo_TurnOn50_1pr->Write();
  hCalo_TurnOn50_3pr->Write();
  hCalo_TurnOn50_withNeutrals->Write();
  hCalo_TurnOn50_noNeutrals->Write();
  //
  hCaloIso_TurnOn25->Write();
  hCaloIso_TurnOn25_1pr->Write();
  hCaloIso_TurnOn25_3pr->Write();
  hCaloIso_TurnOn25_withNeutrals->Write();
  hCaloIso_TurnOn25_noNeutrals->Write();
  hCaloIso_TurnOn50->Write();
  hCaloIso_TurnOn50_1pr->Write();
  hCaloIso_TurnOn50_3pr->Write();
  hCaloIso_TurnOn50_withNeutrals->Write();
  hCaloIso_TurnOn50_noNeutrals->Write();

  outFile->Write();

  return;
}


//****************************************************************************
bool CaloTau::IsWithinEtaRegion(string etaRegion,
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
void CaloTau::FinaliseEffHisto_(TH1D *histo, 
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
void CaloTau::FinaliseEffHisto_(TH2D *histo, 
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
void CaloTau::FillSingleTau_(vector<L1TkTauParticle> L1Taus, 
			    TH1D *hRate,
			    TH1D *hEfficiency,
			    double minEta,
			    double maxEta)
//============================================================================
{

  int iLdg = -1;
  int iSubldg = -1;
  GetLdgAndSubldgIndices(L1Taus, iLdg, iSubldg);

  // Sanity check
  if( L1Taus.size() < 1 ) return;
  
  // Fill rate
  if (0) std::cout << "1) iLdg = " << iLdg << ", L1Taus.size() = " << L1Taus.size() << std::endl;
  // TLorentzVector sigTks_p4 = GetSigConeTTTracksP4();
  double ldgEt  = L1Taus.at(iLdg).GetCaloTau().et();
  double ldgEta = L1Taus.at(iLdg).GetCaloTau().eta();
  
  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(ldgEta) < minEta) return;
  if ( abs(ldgEta) > maxEta) return;
    
  FillRate_(hRate, ldgEt);
  
  // Get MC-matched trigger objects
  vector<L1TkTauParticle> L1Taus_mcMatched = GetMcMatchedL1Taus(L1Taus);
  if (L1Taus_mcMatched.size() < 1) return;

  // Check that all taus were found
  if(!bFoundAllTaus_) return;

  int iLdgMC = -1;
  int iSubldgMC = -1;
  GetLdgAndSubldgIndices(L1Taus_mcMatched, iLdgMC, iSubldgMC);  

  // Fill efficiency
  if (0) std::cout << "3) iLdgMC = " << iLdgMC << ", L1Taus_mcMatched.size() = " << L1Taus_mcMatched.size() << std::endl;
  if (0) L1Taus_mcMatched.at(iLdgMC).PrintProperties(false, true, false, false);
  double ldgEt_mcMatched = L1Taus_mcMatched.at(iLdg).GetCaloTau().et();
  FillEfficiency_(hEfficiency, ldgEt_mcMatched);

  return;
}


//============================================================================
void CaloTau::FillDiTau_(vector<L1TkTauParticle> L1Taus, 
			TH1D *hRate,
			TH1D *hEfficiency,
			double minEta,
			double maxEta)
//============================================================================
{

  // Sanity check
  if( L1Taus.size() < 2 ) return;  

  // Get Ldg and Subldg object indices
  int iLdg = -1;
  int iSubldg = -1;
  GetLdgAndSubldgIndices(L1Taus, iLdg, iSubldg);

  // Fill rate
  if (0) std::cout << "4) iSubldg = " << iSubldg << ", L1Taus.size() = " << L1Taus.size() << std::endl;
  // TLorentzVector sigTks_p4 = L1Taus.at(iSubldg).GetSigConeTTTracksP4();
  double subLdgEt  = L1Taus.at(iSubldg).GetCaloTau().et();
  double subLdgEta = L1Taus.at(iSubldg).GetCaloTau().eta();

  // Inclusive or Eta slice in Central/Intermedieate/Forward Tracker region?
  if ( abs(subLdgEta) < minEta) return;
  if ( abs(subLdgEta) > maxEta) return;
  if ( abs(subLdgEta) < minEta) return;
  if ( abs(subLdgEta) > maxEta) return;

  FillRate_(hRate, subLdgEt);

  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1Taus_mcMatched = GetMcMatchedL1Taus(L1Taus);
  if (L1Taus_mcMatched.size() < 2) return;
    
  // Check that all taus were found
  if(!bFoundAllTaus_) return;

  // Get Ldg and Subldg object indices
  int iLdgMC = -1;
  int iSubldgMC = -1;
  GetLdgAndSubldgIndices(L1Taus_mcMatched, iLdgMC, iSubldgMC);

  // Fill efficiency
  if (0) std::cout << "5) iSubldgMC = " << iSubldgMC << ", L1Taus_mcMatched.size() = " << L1Taus_mcMatched.size() << std::endl;
  // TLorentzVector sigTks_p4_mc = L1Taus_mcMatched.at(iSubldgMC).GetSigConeTTTracksP4(); // fixme
  double subLdgEt_mcMatched  = L1Taus_mcMatched.at(iSubldg).GetCaloTau().et();
  FillEfficiency_(hEfficiency, subLdgEt_mcMatched);

  return;
}


//============================================================================
void CaloTau::FillDiTau_(vector<L1TkTauParticle> L1Taus1,
			vector<L1TkTauParticle> L1Taus2, 
			TH2D *hRate,
			TH2D *hEfficiency)
//============================================================================
{

  // Sanity check
  if( L1Taus1.size() < 1 ) return;
  if( L1Taus2.size() < 1 ) return;
  
  // Get MC-Matched trigger objects
  vector<L1TkTauParticle> L1Taus1_mcMatched = GetMcMatchedL1Taus(L1Taus1);
  vector<L1TkTauParticle> L1Taus2_mcMatched = GetMcMatchedL1Taus(L1Taus2);

  // Get Ldg and Subldg object indices
  int iLdg1, iLdg2, iLdg1MC, iLdg2MC = -1;
  int iSubldg1, iSubldg2, iSubldg1MC, iSubldg2MC = -1;
  GetLdgAndSubldgIndices(L1Taus1, iLdg1, iSubldg1);
  GetLdgAndSubldgIndices(L1Taus2, iLdg2, iSubldg2);
  GetLdgAndSubldgIndices(L1Taus1_mcMatched, iLdg1MC, iSubldg1MC);
  GetLdgAndSubldgIndices(L1Taus2_mcMatched, iLdg2MC, iSubldg2MC);

  // Fill rate 
  if (0) std::cout << "6) iLdg1 = " << iLdg1 << ", L1Taus1.size() = " << L1Taus1.size() << std::endl;
  if (0) std::cout << "7) iLdg2 = " << iLdg2 << ", L1Taus2.size() = " << L1Taus2.size() << std::endl;
  // TLorentzVector sigTks1_p4 = L1Taus1.at(iLdg1).GetSigConeTTTracksP4();
  // TLorentzVector sigTks2_p4 = L1Taus2.at(iLdg2).GetSigConeTTTracksP4();
  double ldgEt1 = L1Taus1.at(iLdg1).GetCaloTau().et();
  double ldgEt2 = L1Taus2.at(iLdg2).GetCaloTau().et();

  // Ensure that different calo objects are used
  double eta1 = L1Taus1.at(iLdg1).GetMatchingTk().getEta();
  double phi1 = L1Taus1.at(iLdg1).GetMatchingTk().getPhi();
  double eta2 = L1Taus2.at(iLdg2).GetMatchingTk().getEta();
  double phi2 = L1Taus2.at(iLdg2).GetMatchingTk().getPhi();
  double dR = auxTools_.DeltaR(eta1, phi1, eta2, phi2);
  if (dR < 0.4)
    {
      if (L1Taus2.size() < 2) return;
      if (0) std::cout << "8) iSubldg1 = " << iSubldg1 << ", L1Taus2.size() = " << L1Taus2.size() << std::endl;
      ldgEt2 = L1Taus2.at(iSubldg2).GetCaloTau().et(); //fixme! better clean-up needed
    }

  // Make x-axis the ldgEt axis
  if (ldgEt1 > ldgEt2) FillRate_(hRate, ldgEt1, ldgEt2); 
  else FillRate_(hRate, ldgEt2, ldgEt1);
  
  // Get MC-matched trigger objects
  if (L1Taus1_mcMatched.size() < 1) return;
  if (L1Taus2_mcMatched.size() < 1) return;

  // Get MC-matched Et
  if (0) std::cout << "9) iSubldg1 = " << iLdg1MC << ", L1Taus1_mcMatched.size() = " << L1Taus1_mcMatched.size() << std::endl;
  if (0) std::cout << "10) iSubldg2 = " << iLdg2MC << ", L1Taus2_mcMatched.size() = " << L1Taus2_mcMatched.size() << std::endl;
  // TLorentzVector sigTks1_p4_mc = L1Taus1_mcMatched.at(iLdg1MC).GetSigConeTTTracksP4();
  // TLorentzVector sigTks2_p4_mc = L1Taus2_mcMatched.at(iLdg2MC).GetSigConeTTTracksP4();
  double ldgEt1_mcMatched = L1Taus1.at(iLdg1MC).GetCaloTau().et();
  double ldgEt2_mcMatched = L1Taus2.at(iLdg2MC).GetCaloTau().et();

  // Ensure that different calo objects are used
  eta1 = L1Taus1_mcMatched.at(iLdg1MC).GetMatchingTk().getEta();
  phi1 = L1Taus1_mcMatched.at(iLdg1MC).GetMatchingTk().getPhi();
  eta2 = L1Taus2_mcMatched.at(iLdg2MC).GetMatchingTk().getEta();
  phi2 = L1Taus2_mcMatched.at(iLdg2MC).GetMatchingTk().getPhi();
  dR    = auxTools_.DeltaR(eta1, phi1, eta2, phi2);
  if (dR < 0.4)
    {
      if (L1Taus2_mcMatched.size() < 2) return;
      ldgEt2_mcMatched = L1Taus2_mcMatched.at(iSubldg2MC).GetCaloTau().et(); //fixme! double-check this
    }
  
  // Check that all taus were found
  if(!bFoundAllTaus_) return;

  // Make x-axis the ldgEt axis
  if (ldgEt1_mcMatched > ldgEt2_mcMatched) histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt1_mcMatched, ldgEt2_mcMatched);
  else histoTools_.FillAllBinsUpToValue_2D(hEfficiency, ldgEt2_mcMatched, ldgEt1_mcMatched);

  return;
}



//============================================================================
void CaloTau::FillRate_(TH1D *hRate,
		       const double ldgEt)
//============================================================================
{
  
  if (ldgEt < 0) return;
  hRate ->Fill( ldgEt );
  
  return;
}


//============================================================================
void CaloTau::FillRate_(TH2D *hRate,
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
void CaloTau::FillEfficiency_(TH1D *hEfficiency,
			     const double ldgEt)
//============================================================================
{
  
  histoTools_.FillAllBinsUpToValue_1D(hEfficiency, ldgEt);

  return;
}


//============================================================================
void CaloTau::FillTurnOn_Numerator_(vector<L1TkTauParticle> L1Taus, 
				   const double minEt,
				   TH1D *hTurnOn,
				   TH1D *hTurnOn_1pr,
				   TH1D *hTurnOn_3pr,
				   TH1D *hTurnOn_withNeutrals,
				   TH1D *hTurnOn_noNeutrals)
//============================================================================
{
  
  // For-loop: L1Taus
  for (vector<L1TkTauParticle>::iterator L1TkTau = L1Taus.begin(); L1TkTau != L1Taus.end(); L1TkTau++)
    {
      
      // Skip if trigger object is not MC matched
      if (!L1TkTau->HasMatchingGenParticle()) continue;	 
      
      // Skip if trigger object has eT < minEt
      // TLorentzVector sigTks2_p4 = L1TkTau->GetSigConeTTTracksP4();
      double tkTau_et = L1TkTau->GetCaloTau().et();
      if ( tkTau_et < minEt) continue;
      
      // Get MC-match
      GenParticle p = L1TkTau->GetMatchingGenParticle();
      
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
      
    } // For-loop: L1Taus

  return;
   
}


//============================================================================
void CaloTau::GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
				     vector<GenParticle> hadGenTaus)					    
//============================================================================
{

  // Sanity check
  if (hadGenTaus.size() < 1 ) return;

  // Initialise
  double caloEta   = L1TkTau.GetCaloTau().eta();
  double caloPhi   = L1TkTau.GetCaloTau().phi();
  double match_dR  = 9999.9;
  GenParticle match_GenParticle;
  
  // For-loop: All hadronic GenTaus
  for (vector<GenParticle>::iterator tau = hadGenTaus.begin(); tau != hadGenTaus.end(); tau++)
    {
      // If no hadronic decay products found (pi+/-, pi0, K+/-, K0, K0L), skip this tau
      if (tau->finalDaughtersCharged().size() < 1) continue;

      if (0) tau->PrintFinalDaughtersCharged();

      TLorentzVector p4charged = tau->p4charged(false);

      double deltaR = auxTools_.DeltaR( p4charged.Eta(), p4charged.Phi(), caloEta, caloPhi ); 
      // std::cout << "dR = " << deltaR << ", caloEta = " << caloEta << ", caloPhi = " << caloPhi << ", p4charged.Eta() = "<< p4charged.Eta() << ", p4charged.Phi() = " << p4charged.Phi() << std::endl;
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
vector<L1TkTauParticle> CaloTau::GetMcMatchedL1Taus(vector<L1TkTauParticle> L1Taus)
//============================================================================
{

  // Get all MC-matched trigger objects
  vector<L1TkTauParticle> matchedL1Taus;
  for (vector<L1TkTauParticle>::iterator tau = L1Taus.begin(); tau != L1Taus.end(); tau++)
    {
      if (!tau->HasMatchingGenParticle()) continue;
      matchedL1Taus.push_back(*tau);
    }
  
  return matchedL1Taus;
}

//============================================================================
void CaloTau::GetLdgAndSubldgIndices(vector<L1TkTauParticle> myTaus,
				    int &iLdg,
				    int &iSubldg)
{
  //============================================================================
  
  // Declarations
  iLdg = 0.0;
  iSubldg = 0.0;
  int index = -1;
  double EtLdg = 0.0;
  double EtSubldg = 0.0;
  
  // Skip if empty
  if (myTaus.size() == 0) return;

  // For-loop: All taus candidates
  for (vector<L1TkTauParticle>::iterator tau = myTaus.begin(); tau != myTaus.end(); tau++)
    {
      index++;
      
      double Et = tau->GetCaloTau().et(); // tau->GetSigConeTTTracksP4().Et();
      if (0) std::cout << index << ") Et = " << Et << std::endl;
      
      if (Et > EtLdg)
	{
	  EtLdg = Et;
	  iLdg  = index;
	}
      else if (Et > EtSubldg)
	{
	  EtSubldg = Et;
	  iSubldg  = index;
	}
      else{}
	
      // std::cout << "EtLdg = " << EtLdg << ", EtSubldg = " << EtSubldg << std::endl;
    }

  if (0)
    {
      std::cout << "================================================" << std::endl;
      std::cout << "EtLdg = " << EtLdg << ", EtSubldg = " << EtSubldg << std::endl;
      std::cout << "iLDg = " << iLdg << ", iSubldg = " << iSubldg << std::endl;
      std::cout << "================================================" << std::endl;
    }

  // Sanity check 
  if (iLdg > int(myTaus.size()-1))
    {
      std::cout << "=== CaloTk::GetLdgAndSubldgIndices() ERROR!\n\tiLdg = " << iLdg << ", myTaus.size() = " << myTaus.size() << std::endl;
    }

  if (iSubldg > int(myTaus.size()-1))
    {
      std::cout << "=== CaloTk::GetLdgAndSubldgIndices() WARNING!\n\tiSubldg = " << iSubldg << ", myTaus.size() = " << myTaus.size() << std::endl;
    }

  return;
}

#endif


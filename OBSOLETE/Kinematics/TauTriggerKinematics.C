#ifndef TauTriggerKinematics_cxx
#define TauTriggerKinematics_cxx

#include "TauTriggerKinematics.h"
#include "../utilities/constants.h"

//#define DEBUG

//****************************************************************************
void TauTriggerKinematics::InitVars()
//****************************************************************************
{

  // Variable initialisation
  tauVisEtCut  = 20.0;


  // Sample-depenend settings: 
  if (mcSample.find("MinBias") != string::npos){
    realTauMom = 0;
    nMaxNumOfHTausPossible = 0;
  }
  else if (mcSample.find("VBF") != string::npos){
    realTauMom = 25;
    nMaxNumOfHTausPossible = 2;
  }
  else if( (mcSample.find("HPlus160") != string::npos) || (mcSample.find("HPlus200") != string::npos) ) {
    realTauMom = 37;
    nMaxNumOfHTausPossible = 1;
  }
  else if (mcSample.find("TTBar") != string::npos){
    realTauMom = 24;
    nMaxNumOfHTausPossible = 2;
  }
  else if (mcSample.find("SingleTauGun1p") != string::npos){
    realTauMom = 0;
    nMaxNumOfHTausPossible = 1;
  }
  else if (mcSample.find("DiTauGun3p") != string::npos){
    realTauMom = 0;
    nMaxNumOfHTausPossible = 2;
  }
  else{
    cout << "E R R O R ! TauTrigger::InitVars_(...) - Unknown sample \"" << mcSample << "\". EXIT" << endl;
    exit(1);
  }

  return;
}


//****************************************************************************
void TauTriggerKinematics::Loop()
//****************************************************************************
{
  
  TauTriggerKinematics::InitVars();
  
  if (fChain == 0) return;
  
  const Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : std::min((int)fChain->GetEntries(), MaxEvents);

  std::cout << "Analyzing " << nEntries << " events.\n";

  ///////////////////
  // Initialisations
  ///////////////////
  BookHistos();
  Long64_t nbytes       = 0;
  Long64_t nb           = 0;

  // For-loop: Entries
  for (Int_t jentry=0; jentry < nEntries; jentry++){

    // Init loop variables + Sanity checks
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Declarations
    int nTaus_1pr = 0;
    int nTaus_3pr = 0;
    std::vector<int> hTaus_index;
    TLorentzVector diTau_P4;
    TLorentzVector diTau_VisP4;
    
    // Tau Properties
    hTaus_index = GetTriggerHadronicTaus(realTauMom, tauVisEtCut);
   
    // For-loop: All hadronic taus
    for (Size_t i = 0; i < hTaus_index.size(); i++) {

      const int tauIndex = hTaus_index.at(i);
      // Declarations
      TLorentzVector hTau_P4;
      TLorentzVector hTau_VisP4;
      TLorentzVector hTau_ChargedPionsP4;
      TLorentzVector hTau_NeutralPionsP4;
      TLorentzVector hTau_LdgChPion_P4;
      std::vector<unsigned short> hTau_Daughters;
      std::vector<unsigned short> chargedPions;
      std::vector<unsigned short> neutralPions;

      // Get charged and neutral pions. Check tauDecayMode
      GetHadronicTauChargedPions(tauIndex, chargedPions);
      GetHadronicTauNeutralPions(tauIndex, neutralPions);      

      if (tauDecayMode.compare("1pr") == 0){ if(chargedPions.size() != 1)  continue;}
      else if (tauDecayMode.compare("3pr") == 0){ if (chargedPions.size() != 3) continue;}
      else if (tauDecayMode.compare("all") == 0) {}
      else{
	std::cout << "E R R O R ! TauTriggerKinematics::InitVars(...) - Unknown tau decay mode \"" << tauDecayMode << "\". EXIT" << std::endl;
	exit(1);
      }

      // Tau properties
      hTau_P4.SetPtEtaPhiM(GenP_Pt->at(tauIndex), GenP_Eta->at(tauIndex), GenP_Phi->at(tauIndex), GenP_Mass->at(tauIndex));
      diTau_P4 = diTau_P4 + hTau_P4;
      
      hHadronicTau_Pt     ->Fill( GenP_Pt     ->at(tauIndex) );
      hHadronicTau_Eta    ->Fill( GenP_Eta    ->at(tauIndex) );
      hHadronicTau_Phi    ->Fill( GenP_Phi    ->at(tauIndex) );
      hHadronicTau_Mass   ->Fill( GenP_Mass   ->at(tauIndex) );
      hHadronicTau_Charge ->Fill( GenP_Charge ->at(tauIndex) );
      hHadronicTau_PdgId  ->Fill( GenP_PdgId  ->at(tauIndex) );
      hHadronicTau_Status ->Fill( GenP_Status ->at(tauIndex) );
      hHadronicTau_VertexX->Fill( GenP_VertexX->at(tauIndex) );
      hHadronicTau_VertexY->Fill( GenP_VertexY->at(tauIndex) );
      hHadronicTau_VertexZ->Fill( GenP_VertexZ->at(tauIndex) );


      // Visible-Tau properties
      GetHadronicTauFinalDaughters(tauIndex, hTau_Daughters);
      hTau_VisP4 = GetVisibleP4(hTau_Daughters);
      diTau_VisP4 = diTau_VisP4 + hTau_VisP4;
      
      hHadronicTau_VisEt    ->Fill( hTau_VisP4.Et() );
      hHadronicTau_VisEta   ->Fill( hTau_VisP4.Eta() );
      hHadronicTau_VisPhi   ->Fill( hTau_VisP4.Phi() );
      hHadronicTau_VisMass  ->Fill( hTau_VisP4.M() );
      hHadronicTau_DecayMode->Fill( GetGenTauDecayMode(hTau_Daughters) );
      
      // For-loop: Tau charged pions
      hHadronicTau_ChargedPion_N->Fill( chargedPions.size() );
      if (chargedPions.size() == 1 ) nTaus_1pr++;
      if (chargedPions.size() == 3 ) nTaus_3pr++;
    
      for (Size_t j = 0; j < chargedPions.size(); j++) {
      
	const int pionIndex = chargedPions.at(j);
	TLorentzVector tmp_P4;
	tmp_P4.SetPtEtaPhiM(GenP_Pt->at(pionIndex), GenP_Eta->at(pionIndex), GenP_Phi->at(pionIndex), GenP_Mass->at(pionIndex));
	hTau_ChargedPionsP4 += tmp_P4;

	hHadronicTau_ChargedPion_Pt     ->Fill( GenP_Pt     ->at(pionIndex) );
	hHadronicTau_ChargedPion_Eta    ->Fill( GenP_Eta    ->at(pionIndex) );
	hHadronicTau_ChargedPion_Phi    ->Fill( GenP_Phi    ->at(pionIndex) );
	hHadronicTau_ChargedPion_Mass   ->Fill( GenP_Mass   ->at(pionIndex) );
	hHadronicTau_ChargedPion_Charge ->Fill( GenP_Charge ->at(pionIndex) );
	hHadronicTau_ChargedPion_PdgId  ->Fill( GenP_PdgId  ->at(pionIndex) );
	hHadronicTau_ChargedPion_Status ->Fill( GenP_Status ->at(pionIndex) );
	hHadronicTau_ChargedPion_VertexX->Fill( GenP_VertexX->at(pionIndex) );
	hHadronicTau_ChargedPion_VertexY->Fill( GenP_VertexY->at(pionIndex) );
	hHadronicTau_ChargedPion_VertexZ->Fill( GenP_VertexZ->at(pionIndex) );

      } // For-loop: Tau charged pions
    

      // For-loop: Tau neutral pions
      hHadronicTau_NeutralPion_N->Fill( neutralPions.size() );
      for (Size_t k = 0; k < neutralPions.size(); k++) {
      
	const int piZeroIndex = neutralPions.at(k);
	TLorentzVector tmp2_P4;
	tmp2_P4.SetPtEtaPhiM(GenP_Pt->at(piZeroIndex), GenP_Eta->at(piZeroIndex), GenP_Phi->at(piZeroIndex), GenP_Mass->at(piZeroIndex));
	hTau_NeutralPionsP4 += tmp2_P4;

	hHadronicTau_NeutralPion_Pt     ->Fill( GenP_Pt     ->at(piZeroIndex) );
	hHadronicTau_NeutralPion_Eta    ->Fill( GenP_Eta    ->at(piZeroIndex) );
	hHadronicTau_NeutralPion_Phi    ->Fill( GenP_Phi    ->at(piZeroIndex) );
	hHadronicTau_NeutralPion_Mass   ->Fill( GenP_Mass   ->at(piZeroIndex) );
	hHadronicTau_NeutralPion_Charge ->Fill( GenP_Charge ->at(piZeroIndex) );
	hHadronicTau_NeutralPion_PdgId  ->Fill( GenP_PdgId  ->at(piZeroIndex) );
	hHadronicTau_NeutralPion_Status ->Fill( GenP_Status ->at(piZeroIndex) );
	hHadronicTau_NeutralPion_VertexX->Fill( GenP_VertexX->at(piZeroIndex) );
	hHadronicTau_NeutralPion_VertexY->Fill( GenP_VertexY->at(piZeroIndex) );
	hHadronicTau_NeutralPion_VertexZ->Fill( GenP_VertexZ->at(piZeroIndex) );

      } // For-loop: Tau neutral pions


      // Leading charged pion
      Int_t LdgChPionIndex = GetLdgDaughterIndex(hTau_Daughters, true);
      hTau_LdgChPion_P4.SetPtEtaPhiM(GenP_Pt->at(LdgChPionIndex), GenP_Eta->at(LdgChPionIndex), GenP_Phi->at(LdgChPionIndex), GenP_Mass->at(LdgChPionIndex));

      hHadronicTau_LdgChPion_Pt      ->Fill( hTau_LdgChPion_P4.Pt()  );
      hHadronicTau_LdgChPion_Eta     ->Fill( hTau_LdgChPion_P4.Eta() );
      hHadronicTau_LdgChPion_Phi     ->Fill( hTau_LdgChPion_P4.Phi() );
      hHadronicTau_LdgChPion_Mass    ->Fill( GenP_Mass   ->at(LdgChPionIndex) );
      hHadronicTau_LdgChPion_Charge  ->Fill( GenP_Charge ->at(LdgChPionIndex) );
      hHadronicTau_LdgChPion_PdgId   ->Fill( GenP_PdgId  ->at(LdgChPionIndex) );
      hHadronicTau_LdgChPion_Status  ->Fill( GenP_Status ->at(LdgChPionIndex) );
      hHadronicTau_LdgChPion_VertexX ->Fill( GenP_VertexX->at(LdgChPionIndex) );
      hHadronicTau_LdgChPion_VertexY ->Fill( GenP_VertexY->at(LdgChPionIndex) );
      hHadronicTau_LdgChPion_VertexZ ->Fill( GenP_VertexZ->at(LdgChPionIndex) );
      hHadronicTau_LdgChPion_Rtau    ->Fill( hTau_LdgChPion_P4.P() / hTau_VisP4.P() );
      hHadronicTau_LdgChPion_Pt_VisEt->Fill( hTau_LdgChPion_P4.Pt() , hTau_VisP4.Et() );     

      // 3-prong
      if (chargedPions.size() >= 3)
	{
	  double maxDeltaR = GetHadronicTauMaxSignalCone(hTau_Daughters);
	  hHadronicTau_LdgChPion_DeltaRMax         ->Fill(maxDeltaR);
	  hHadronicTau_LdgChPion_DeltaRMax_Pt      ->Fill(maxDeltaR, hTau_LdgChPion_P4.Pt() );
	  hHadronicTau_LdgChPion_DeltaRMax_TauEt   ->Fill(maxDeltaR, hTau_P4.Et() );
	  hHadronicTau_LdgChPion_DeltaRMax_VisTauEt->Fill(maxDeltaR, hTau_VisP4.Et() );

	  double maxDeltaR_MinPt = GetHadronicTauMaxSignalCone(hTau_Daughters, true, 2.0);
	  hHadronicTau_LdgChPion_DeltaRMax_MinPt         ->Fill(maxDeltaR_MinPt);
	  hHadronicTau_LdgChPion_DeltaRMax_Pt_MinPt      ->Fill(maxDeltaR_MinPt, hTau_LdgChPion_P4.Pt() );
	  hHadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt   ->Fill(maxDeltaR_MinPt, hTau_P4.Et() );
	  hHadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt->Fill(maxDeltaR_MinPt, hTau_VisP4.Et() );
	}

      // Visible-Tau properties (again)
      hHadronicTau_EcalFraction->Fill( hTau_NeutralPionsP4.Et() / hTau_VisP4.Et() );
      hHadronicTau_HcalFraction->Fill( hTau_ChargedPionsP4.Et() / hTau_VisP4.Et() );

    } // For-loop: All hadronic taus
    
    // Fill Histograms outside loop
    hHadronicTau_N->Fill( hTaus_index.size() );    

    bool bFillDiTaus = false;
    if ( (nMaxNumOfHTausPossible > 1) && (hTaus_index.size() > 1) ) {

      if ( (tauDecayMode.compare("1pr") == 0) && (nTaus_1pr == 2) ) bFillDiTaus = true;
      else if ( (tauDecayMode.compare("3pr") == 0) && (nTaus_3pr == 2) ) bFillDiTaus = true;
      else if ( (tauDecayMode.compare("all") == 0) )bFillDiTaus = true;
      else bFillDiTaus = false;

      if (bFillDiTaus==true)
	{
	  if (diTau_P4.Pt() < 1 ) std::cout << "tauDecayMode = " << tauDecayMode << ", nTaus_1pr = " << nTaus_1pr << ", nTaus_3pr = " << nTaus_3pr << std::endl;
	  hDiTau_Pt     ->Fill( diTau_P4.Pt()   );
	  hDiTau_InvMass->Fill( diTau_P4.M()    );
	  hDiTau_VisMass->Fill( diTau_VisP4.M() );
	}

    }
    
    // Update progress bar
#ifndef DEBUG
    tools.ProgressBar(jentry, nEntries, 50, 50);
#endif

  }// For-loop: Entries



  ///////////////////////////////////////////////////////////
  // Normalisation/Finalisation of Histos outside loop
  ///////////////////////////////////////////////////////////

  // Write the histograms to the file
  outFile->cd();
  outFile->Write();
  std::cout << "\nDone.\n" << std::endl;

}


//****************************************************************************
void TauTriggerKinematics::BookHistos(void)
//****************************************************************************
{
  
  // Tau
  histos.BookHisto_1D( hHadronicTau_N      , "HadronicTau_N"      ,  20,  +0.00,  + 20.0);
  histos.BookHisto_1D( hHadronicTau_Pt     , "HadronicTau_Pt"     , 400,  +0.00,  +400.0);
  histos.BookHisto_1D( hHadronicTau_Eta    , "HadronicTau_Eta"    , 300,  -3.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_Phi    , "HadronicTau_Phi"    , 160,  -4.00,  +  4.00);
  histos.BookHisto_1D( hHadronicTau_Mass   , "HadronicTau_Mass"   , 300,  +0.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_Charge , "HadronicTau_Charge" ,   5,  -2.50,  +  2.50);
  histos.BookHisto_1D( hHadronicTau_PdgId  , "HadronicTau_PdgId"  ,  20,  +0.00,  + 20.00);
  histos.BookHisto_1D( hHadronicTau_Status , "HadronicTau_Status" ,  10,  +0.00,  + 10.00);
  histos.BookHisto_1D( hHadronicTau_VertexX, "HadronicTau_VertexX",  60, -30.00,  + 30.00);
  histos.BookHisto_1D( hHadronicTau_VertexY, "HadronicTau_VertexY",  60, -30.00,  + 30.00);
  histos.BookHisto_1D( hHadronicTau_VertexZ, "HadronicTau_VertexZ",  60, -30.00,  + 30.00);

  // Visible-Tau
  histos.BookHisto_1D( hHadronicTau_VisEt       , "HadronicTau_VisEt"       , 400,  +0.00,  +400.00);
  histos.BookHisto_1D( hHadronicTau_VisEta      , "HadronicTau_VisEta"      , 300,  -3.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_VisPhi      , "HadronicTau_VisPhi"      , 160,  -4.00,  +  4.00);
  histos.BookHisto_1D( hHadronicTau_VisMass     , "HadronicTau_VisMass"     , 300,  +0.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_DecayMode   , "HadronicTau_DecayMode"   , 100,  +0.00,  +100.00);
  histos.BookHisto_1D( hHadronicTau_EcalFraction, "HadronicTau_EcalFraction",  110,  +0.00,  + 1.10);
  histos.BookHisto_1D( hHadronicTau_HcalFraction, "HadronicTau_HcalFraction",  110,  +0.00,  + 1.10);

  // Tau Charged Pions
  histos.BookHisto_1D( hHadronicTau_ChargedPion_N       , "HadronicTau_ChargedPion_N"      ,   10, +   0.00, +  10.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_Pt      , "HadronicTau_ChargedPion_Pt"     ,  400, +   0.00, + 400.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_Eta     , "HadronicTau_ChargedPion_Eta"    ,  300, -   3.00, +   3.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_Phi     , "HadronicTau_ChargedPion_Phi"    ,  160, -   4.00, +   4.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_Mass    , "HadronicTau_ChargedPion_Mass"   ,  300, +   0.00, +   3.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_Charge  , "HadronicTau_ChargedPion_Charge" ,    5, -   2.50, +   2.50);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_PdgId   , "HadronicTau_ChargedPion_PdgId"  , 4000, -2000.00, +2000.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_Status  , "HadronicTau_ChargedPion_Status" ,   10, +   0.00, +  10.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_VertexX , "HadronicTau_ChargedPion_VertexX",   60, -  30.00, +  30.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_VertexY , "HadronicTau_ChargedPion_VertexY",   60, -  30.00, +  30.00);
  histos.BookHisto_1D( hHadronicTau_ChargedPion_VertexZ , "HadronicTau_ChargedPion_VertexZ",   60, -  30.00, +  30.00);

  // Tau Neutral pions
  histos.BookHisto_1D( hHadronicTau_NeutralPion_N       , "HadronicTau_NeutralPion_N"      ,   10, -   0.00,  + 10.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_Pt      , "HadronicTau_NeutralPion_Pt"     ,  400, -   0.00,  +400.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_Eta     , "HadronicTau_NeutralPion_Eta"    ,  300, -   3.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_Phi     , "HadronicTau_NeutralPion_Phi"    ,  160, -   4.00,  +  4.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_Mass    , "HadronicTau_NeutralPion_Mass"   ,  300, -   0.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_Charge  , "HadronicTau_NeutralPion_Charge" ,    5, -   2.50,  +  2.50);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_PdgId   , "HadronicTau_NeutralPion_PdgId"  , 4000, -2000.00, +2000.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_Status  , "HadronicTau_NeutralPion_Status" ,   10, +   0.00,  + 10.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_VertexX , "HadronicTau_NeutralPion_VertexX",   60, -  30.00,  + 30.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_VertexY , "HadronicTau_NeutralPion_VertexY",   60, -  30.00,  + 30.00);
  histos.BookHisto_1D( hHadronicTau_NeutralPion_VertexZ , "HadronicTau_NeutralPion_VertexZ",   60, -  30.00,  + 30.00);

  // Leading Charged Pion
  histos.BookHisto_1D( hHadronicTau_LdgChPion_Pt        , "HadronicTau_LdgChPion_Pt"       ,  400, -   0.00,  +400.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_Eta       , "HadronicTau_LdgChPion_Eta"      ,  600, -   3.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_Phi       , "HadronicTau_LdgChPion_Phi"      ,  800, -   4.00,  +  4.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_Mass      , "HadronicTau_LdgChPion_Mass"     ,  300, +   0.00,  +  3.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_Charge    , "HadronicTau_LdgChPion_Charge"   ,    5, -   2.50,  +  2.50);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_PdgId     , "HadronicTau_LdgChPion_PdgId"    , 4000, -2000.00, +2000.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_Status    , "HadronicTau_LdgChPion_Status"   ,   10, +   0.00,  + 10.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_VertexX   , "HadronicTau_LdgChPion_VertexX"  ,   60, -  30.00,  + 30.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_VertexY   , "HadronicTau_LdgChPion_VertexY"  ,   60, -  30.00,  + 30.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_VertexZ   , "HadronicTau_LdgChPion_VertexZ"  ,   60, -  30.00,  + 30.00);
  histos.BookHisto_1D( hHadronicTau_LdgChPion_Rtau      , "HadronicTau_LdgChPion_Rtau"     ,  24, +   0.00,  +  1.20);
  histos.BookHisto_2D( hHadronicTau_LdgChPion_Pt_VisEt  , "HadronicTau_LdgChPion_Pt_VisEt" ,  400, +   0.00,  + 400.00, 400, 0.00,  +400.00);

  // 3-prong
  histos.BookHisto_1D( hHadronicTau_LdgChPion_DeltaRMax         , "HadronicTau_LdgChPion_DeltaRMax"         , 200, + 0.00,  +1.00);
  histos.BookHisto_2D( hHadronicTau_LdgChPion_DeltaRMax_Pt      , "HadronicTau_LdgChPion_DeltaRMax_Pt"      , 200, + 0.00,  +1.00, 400, +0.00,  +400.00);
  histos.BookHisto_2D( hHadronicTau_LdgChPion_DeltaRMax_TauEt   , "HadronicTau_LdgChPion_DeltaRMax_TauEt"   , 200, + 0.00,  +1.00, 400, +0.00,  +400.00);
  histos.BookHisto_2D( hHadronicTau_LdgChPion_DeltaRMax_VisTauEt, "HadronicTau_LdgChPion_DeltaRMax_VisTauEt", 200, + 0.00,  +1.00, 400, +0.00,  +400.00);

  histos.BookHisto_1D( hHadronicTau_LdgChPion_DeltaRMax_MinPt         , "HadronicTau_LdgChPion_DeltaRMax_MinPt"         , 200, +0.00,  +1.00);
  histos.BookHisto_2D( hHadronicTau_LdgChPion_DeltaRMax_Pt_MinPt      , "HadronicTau_LdgChPion_DeltaRMax_Pt_MinPt"      , 200, +0.00,  +1.00, 400, +0.0, +400.00);
  histos.BookHisto_2D( hHadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt   , "HadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt"   , 200, +0.00,  +1.00, 400, +0.0, +400.00);
  histos.BookHisto_2D( hHadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt, "HadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt", 200, +0.00,  +1.00, 400, +0.0, +400.00);

  // DiTau System
  histos.BookHisto_1D( hDiTau_Pt     , "DiTau_Pt"     , 400, +0.00 , +400.0);
  histos.BookHisto_1D( hDiTau_InvMass, "DiTau_InvMass", 200, +0.00 , +200.0);
  histos.BookHisto_1D( hDiTau_VisMass, "DiTau_VisMass", 200, +0.00 , +200.0);

  return;
}


//****************************************************************************
std::vector<int> TauTriggerKinematics::GetTriggerHadronicTaus(const int tauMom, 					  
						    const double hTauVisEtCut)
//****************************************************************************
{

  std::vector<int> hTaus_index;

  // For-loop: GenParticles
  for (Size_t iGenP = 0; iGenP < GenP_PdgId->size(); iGenP++) {
    
    // Get hadronically decaying taus within acceptance
    bool bIsHTau = IsTriggerHadronicTau(iGenP, true, tauMom, true, hTauVisEtCut);
    if(!bIsHTau) continue;
    
    // Push back GenP index hadronic tau 
    hTaus_index.push_back(iGenP);
    
  } // For-loop: GenParticles
  
  return hTaus_index;
}


//****************************************************************************
TLorentzVector TauTriggerKinematics::GetHadronicTauVisP4(const int genP_index)
//****************************************************************************
{

  // Get visible 4-momentum of hadronic tau
  std::vector<unsigned short> hTau_Prods;
  GetHadronicTauFinalDaughters(genP_index, hTau_Prods);
  TLorentzVector hTau_VisP4 = GetVisibleP4(hTau_Prods);
  return hTau_VisP4;
}

#endif // TauTriggerKinematics_cxx

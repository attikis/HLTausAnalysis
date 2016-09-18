#ifndef Tracking_cxx
#define Tracking_cxx

#include "Tracking.h"
#include "../utilities/constants.h"
#include "TFitResult.h"
#include "TF1.h"

//#define bDebug

//****************************************************************************
void Tracking::InitVars_(void)
//****************************************************************************
{
  
  // Options
  bVerbose           = false;
  bDoPixelTracks     = true;
  bPrintEfficiencies = false;
  bPrintResolutions  = true;
  bPrintFitInfo      = false;
  bSaveFitInfo       = false;
  tp_dxy_Max         = 1.0;

  // Pixel Cuts
  pixTk_NPixHits_Min = 3;
  pixTk_NPixHits_Max = 1e6; 
  // pixTk_MustHit_Type.push_back(+1);
  // pixTk_MustHit_Type.push_back(-1);
  
  // Initialisations 
  PT_MAX       = 50.0;
  PT_BINS      = 10;
  PT_BINWIDTH  = PT_MAX/(double)PT_BINS;
  ETA_MAX      = 2.4; // default: 2.5 
  ETA_BINS     = 12;  // default: 25
  ETA_BINWIDTH = ETA_MAX/(double)ETA_BINS;

  n_all_eta_C        = 0;
  n_all_eta_I        = 0;
  n_all_eta_F        = 0;
  n_match_eta_C      = 0;
  n_match_eta_I      = 0;
  n_match_eta_F      = 0;
  //
  n_all_eta_C_pt_L   = 0;
  n_all_eta_I_pt_L   = 0;
  n_all_eta_F_pt_L   = 0;
  n_match_eta_C_pt_L = 0;
  n_match_eta_I_pt_L = 0;
  n_match_eta_F_pt_L = 0;
  //
  n_all_eta_C_pt_M   = 0;
  n_all_eta_I_pt_M   = 0;
  n_all_eta_F_pt_M   = 0;
  n_match_eta_C_pt_M = 0;
  n_match_eta_I_pt_M = 0;
  n_match_eta_F_pt_M = 0;
  //
  n_all_eta_C_pt_H   = 0;
  n_all_eta_I_pt_H   = 0;
  n_all_eta_F_pt_H   = 0;
  n_match_eta_C_pt_H = 0;
  n_match_eta_I_pt_H = 0;
  n_match_eta_F_pt_H = 0;

  // Event Pixel Hits
  pixHits_XYZ.clear();
  pixHits_TTPixelTrackIndex.clear();  

  _pt_acceptance     =   0.2;
  _eta_acceptance    =   2.5;
  _z0_acceptance     =  30.0;
  _pt_L              =   6.0; //  5.0
  _pt_H              =  16.0; // 15.0
  _eta_C             =   0.8; 
  _eta_F             =   1.6;
  _pt_tkMin          =   2.0;
  _nStubs_tkMin      =   4;
  _chiSq_overflow    = 100.0;
  _redChiSq_overflow =  15.0;
  _fromCentimetresToMicrometers = 10000.0;  // 1) cm->m: 100 2) m->micrometers*1000000;
  //
  s_eta_all          = "\\abs{\\eta} \\leq 2.5";
  s_eta_C            = "\\abs{\\eta} < " + auxTools_.ToString(_eta_C);
  s_eta_I            = auxTools_.ToString(_eta_C) + " < \\abs{\\eta} < " + auxTools_.ToString(_eta_F);
  s_eta_F            = "\\abs{\\eta} \\geq " + auxTools_.ToString(_eta_F);
  s_pt_tkMin         = "\\pT > " + auxTools_.ToString(_pt_tkMin) + " \\GeVc{-1}";
  s_pt_L             = "\\pT \\leq " + auxTools_.ToString(_pt_L) + " \\GeVc{-1}";
  s_pt_M             = auxTools_.ToString(_pt_L) + " \\GeVc{-1} < \\pT \\leq " + auxTools_.ToString(_pt_H) + " \\GeVc{-1}";
  s_pt_H             = "\\pT \\geq " + auxTools_.ToString(_pt_H) + " \\GeVc{-1}";
  
  return;
}



//****************************************************************************
void Tracking::Loop(){
//****************************************************************************

  if (fChain == 0) return;
  
  Long64_t nEntries      = (MaxEvents == -1) ? fChain->GetEntries() : std::min((int)fChain->GetEntries(), MaxEvents);
  Int_t  NPBarDivisions  = 100;
  Int_t PBarWidth        = 150;
  
  std::cout << "Analyzing " << nEntries << " events.\n";

  InitVars_();
  BookHistos_();  
  auxTools_.PrintPSets(fChain);
  PrintSettings_();
  PrintPixelTrackCuts_();

  Long64_t nbytes = 0, nb = 0;  
  ///////////////////////////////////////////////////////////
  // For-loop: Entries
  ///////////////////////////////////////////////////////////
  for (Int_t jentry=0; jentry < nEntries; jentry++){
      
    // Init loop variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    // Get all pixel hits for given event
    pixHits_XYZ.clear();
    pixHits_TTPixelTrackIndex.clear();  
    s->GetPixTrackAllHits(pixHits_XYZ, pixHits_TTPixelTrackIndex);    

  
    ///////////////////////////////////////////////////////////
    // For-loop: TPs
    ///////////////////////////////////////////////////////////
    for (Int_t t = 0; t < (int) TP_Pt->size(); t++) {

      
      // Get TP properties
      int tp_NMatch = TP_NMatch->at(t);
      double tp_Pt  = TP_Pt->at(t);
      double tp_Eta = TP_Eta->at(t);
      double tp_Phi = TP_Phi->at(t);
      double tp_x0  = TP_POCAx->at(t);
      double tp_y0  = TP_POCAy->at(t);
      double tp_z0  = TP_POCAz->at(t);
      double tp_d0  = tp->GetD0(t);
#ifdef bDebug
      tp->PrintTPProperties(t);
#endif
      
      // Ignore TPs outside our acceptance
      if ( !IsWithinTrackerAcceptance(tp_Pt, tp_Eta, tp_z0) ) continue;

      // Only consider TPs coming from near the IP (Louise has this enabled)
      // N.B.: Significantly affect d0 Resolution if disabled!!!
      double tp_dxy = sqrt(tp_x0*tp_x0 + tp_y0*tp_y0);
      if (tp_dxy > tp_dxy_Max) continue;

      
      // Fill Histograms
      h_tp_pt ->Fill(tp_Pt);
      if (tp_Pt < _pt_L) h_tp_pt_L->Fill(tp_Pt);

      // Counters
      if (tp_Pt > _pt_tkMin) {
	if ( IsWithinEtaRegion("Central", tp_Eta) ){
	  n_all_eta_C++;
	  h_tp_pt_C->Fill(tp_Pt);
	}
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ){
	  n_all_eta_I++;
	  h_tp_pt_I->Fill(tp_Pt);
	}
	else if ( IsWithinEtaRegion("Forward", tp_Eta) ){
	  n_all_eta_F++;
	  h_tp_pt_F->Fill(tp_Pt);
	}
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}
	
	h_tp_eta->Fill(tp_Eta);
	h_tp_phi->Fill(tp_Phi);
	h_tp_z0 ->Fill(tp_z0);
	h_tp_d0 ->Fill(tp_d0 * _fromCentimetresToMicrometers);
      }          

      // Low Pt Range
      if ( IsWithinPtRange("Low", tp_Pt) ){
	if ( IsWithinEtaRegion("Central"          , tp_Eta) ) n_all_eta_C_pt_L++;
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ) n_all_eta_I_pt_L++;
	else if ( IsWithinEtaRegion("Forward"     , tp_Eta) ) n_all_eta_F_pt_L++;
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}

	h_tp_eta_L->Fill(tp_Eta);
      }

      // Middle Pt Range
      if ( IsWithinPtRange("Middle", tp_Pt) ){
	if ( IsWithinEtaRegion("Central"          , tp_Eta) ) n_all_eta_C_pt_M++;
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ) n_all_eta_I_pt_M++;
	else if ( IsWithinEtaRegion("Forward"     , tp_Eta) ) n_all_eta_F_pt_M++;
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}

	h_tp_eta_M->Fill(tp_Eta);
      }

      // High Pt Range
      if ( IsWithinPtRange("High", tp_Pt) ){
	if ( IsWithinEtaRegion("Central"          , tp_Eta) ) n_all_eta_C_pt_H++;
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ) n_all_eta_I_pt_H++;
	else if ( IsWithinEtaRegion("Forward"     , tp_Eta) ) n_all_eta_F_pt_H++;
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}

	h_tp_eta_H->Fill(tp_Eta);
      }
    
      // Was the TP matched to a TTTrack? If not continue
      if (tp_NMatch < 1) continue;


      // NEW! Was the TP matched uniquely to a TTTrack? Inserted by myself on 31 July 2015  (xenios)
      if (tp_NMatch > 1) continue;
     
      ///////////////////////////////////////////////////////////
      // Matching TTTracks (overwritten later on by pixel values if bDoPixelTracks is enabled)
      ///////////////////////////////////////////////////////////
      int tk_Index         = TP_TTTrackIndex->at(t); //t is the TP index
      double tk_Pt         = L1Tks_Pt->at(tk_Index);
      double tk_Eta        = L1Tks_Eta->at(tk_Index);
      double tk_Phi        = L1Tks_Phi->at(tk_Index);
      double tk_z0         = L1Tks_POCAz->at(tk_Index);
      double tk_d0         = s->GetTrackD0(tk_Index); // L1Tks_d0->at(tk_Index);
      double tk_ChiSq      = L1Tks_ChiSquared->at(tk_Index);
      double tk_RedChiSq   = s->GetRedChiSquared(tk_Index, 5);  // L1Tks_RedChiSquared->at(tk_Index);
      int tk_nStubs        = s->GetNumOfStubs(tk_Index);
      if (bVerbose && !bDoPixelTracks) s->PrintTrackProperties(tk_Index);

      // Use only TTTracks with a minimum number of stubs
      if (tk_nStubs < _nStubs_tkMin) continue; 
      
      ///////////////////////////////////////////////////////////
      // Matching TTPixelTracks: OVERWRITE track properties by pixel re-fitted L1 Pixel Tracks
      ///////////////////////////////////////////////////////////
      if (bDoPixelTracks){

	int pixTk_Index = s->GetPixelIndexOfTrack(tk_Index);
	if (pixTk_Index < 0) continue;
	double pixTk_Pt         = L1PixTks_Pt ->at(pixTk_Index);
	double pixTk_Eta        = L1PixTks_Eta->at(pixTk_Index);
	double pixTk_Phi        = L1PixTks_Phi->at(pixTk_Index);
	double pixTk_z0         = L1PixTks_POCAz->at(pixTk_Index);
	double pixTk_d0         = s->GetPixTrackD0(pixTk_Index);
	double pixTk_ChiSq      = L1PixTks_ChiSquared->at(pixTk_Index);
	double pixTk_RedChiSq   = s->GetPixTrackRedChiSq(pixTk_Index);
	int pixTk_PixHits       = L1PixTks_NPixHits->at(pixTk_Index);
	//
	vector<double> pixTk_PixHits_X   = L1PixTks_PixHits_X->at(pixTk_Index);
	vector<double> pixTk_PixHits_Y   = L1PixTks_PixHits_Y->at(pixTk_Index);
	vector<double> pixTk_PixHits_Z   = L1PixTks_PixHits_Z->at(pixTk_Index);
	vector<double> pixTk_PixHits_R   = L1PixTks_PixHits_R->at(pixTk_Index);
	vector<double> pixTk_PixHits_Phi = L1PixTks_PixHits_Phi->at(pixTk_Index);
	vector< int > pixTk_PixHits_Type = L1PixTks_PixHits_Type->at(pixTk_Index);
	// Temporary assignement - handled later on
	int pixTk_SharedPixHits = 0;
	vector<int> pixTk_SharedPixHits_PixTkRefIndex;
	vector<double> pixTk_SharedPixHits_X = pixTk_PixHits_X;
	vector<double> pixTk_SharedPixHits_Y   = pixTk_PixHits_Y;
	vector<double> pixTk_SharedPixHits_Z   = pixTk_PixHits_Z;
	vector<double> pixTk_SharedPixHits_R   = pixTk_PixHits_R;
	vector<double> pixTk_SharedPixHits_Phi = pixTk_PixHits_Phi;
	vector< int > pixTk_SharedPixHits_Type = pixTk_PixHits_Type;
	//
	int pixTk_CandPixHits_N              = L1PixTks_CandPixHits_X->at(pixTk_Index).size();
	vector<double> pixTk_CandPixHits_X   = L1PixTks_CandPixHits_X->at(pixTk_Index);
	vector<double> pixTk_CandPixHits_Y   = L1PixTks_CandPixHits_Y->at(pixTk_Index);
	vector<double> pixTk_CandPixHits_Z   = L1PixTks_CandPixHits_Z->at(pixTk_Index);
	vector<double> pixTk_CandPixHits_R   = L1PixTks_CandPixHits_R->at(pixTk_Index);
	vector<double> pixTk_CandPixHits_Phi = L1PixTks_CandPixHits_Phi->at(pixTk_Index);
	vector< int > pixTk_CandPixHits_Type = L1PixTks_CandPixHits_Type->at(pixTk_Index);
	
	// Apply TTPixelTracks cuts: Number of Pixel Hits
	if (pixTk_PixHits < pixTk_NPixHits_Min) continue;
	if (pixTk_PixHits >= pixTk_NPixHits_Max) continue;
	// if (pixTk_PixHits != pixTk_NPixHits_Min) continue; // xenios
	
	// Apply TTPixelTracks cuts: Required Pixel Hits Type
	bool bHasMustHitLayerOrDisk = s->HasMustPixelHitInLayeOrDisk(pixTk_Index, pixTk_MustHit_Type);
	// bool bHasMustHitLayerOrDisk = s->HasMustPixelHitInLayerAndDisk(pixTk_Index, pixTk_Barrel_Type, pixTk_Endcap_Type, pixTk_Type_EndcapEtaBoundary);
	if (!bHasMustHitLayerOrDisk) continue;

	// Overwrite TTTrack values with TTPixelTrack values!
	tk_Pt         = pixTk_Pt;
	tk_Eta        = pixTk_Eta;
	tk_Phi        = pixTk_Phi;
	tk_z0         = pixTk_z0;
	tk_d0         = pixTk_d0;
	tk_ChiSq      = pixTk_ChiSq;
	tk_RedChiSq   = pixTk_RedChiSq;

	if (bVerbose) s->PrintPixTrackProperties(pixTk_Index, false);

	// Pixel Hits of TTPixelTrack with index pixTk_Index
	h_pixTk_pixHits_N->Fill(pixTk_PixHits);
	for(int i = 0; i < pixTk_PixHits; i++){
	  
	  h_pixTk_pixHits_Rho   ->Fill( pixTk_PixHits_R.at(i) );
	  h_pixTk_pixHits_Z     ->Fill( pixTk_PixHits_Z.at(i) );
	  h_pixTk_pixHits_Type  ->Fill( pixTk_PixHits_Type.at(i) );
	  h_pixTk_pixHits_ZVsRho->Fill( pixTk_PixHits_Z.at(i), pixTk_PixHits_R.at(i) );
	  if (fabs(tk_Eta) < _eta_C) h_pixTk_pixHits_ZVsRho_C->Fill( pixTk_PixHits_Z.at(i), pixTk_PixHits_R.at(i) );
	  else if (fabs(tk_Eta) < _eta_F && fabs(tk_Eta) >= _eta_C) h_pixTk_pixHits_ZVsRho_I->Fill( pixTk_PixHits_Z.at(i), pixTk_PixHits_R.at(i) );
	  else if (fabs(tk_Eta) >= _eta_F) h_pixTk_pixHits_ZVsRho_F->Fill( pixTk_PixHits_Z.at(i), pixTk_PixHits_R.at(i) );
	  else{}
	  if (fabs(tk_Eta) >= 2.0) h_pixTk_pixHits_ZVsRho_EtaGE2->Fill( pixTk_PixHits_Z.at(i), pixTk_PixHits_R.at(i) );	  
	}// for(int i = 0; i < pixTk_PixHits; i++){

	
	// Shared Pixel Hits
	s->GetPixTrackSharedHits(pixTk_Index, pixTk_SharedPixHits_PixTkRefIndex, pixTk_SharedPixHits_X, pixTk_SharedPixHits_Y, pixTk_PixHits_Z,
				 pixTk_PixHits_R, pixTk_PixHits_Phi, pixTk_PixHits_Type, pixHits_XYZ, pixHits_TTPixelTrackIndex);
	
	pixTk_SharedPixHits = pixTk_SharedPixHits_X.size();
	h_pixTk_sharedPixHits_N->Fill(pixTk_SharedPixHits);
	
	if(pixTk_SharedPixHits > 0) {
	  cout << "Event " << jentry << "/" << nEntries << endl;
	  cout << "TP_Pt->size() = " << TP_Pt->size() << ", tp_NMatch = " << tp_NMatch << ", L1PixTks_Pt->size() = " << L1PixTks_Pt->size() << ", L1Tks_Pt->size() = " << L1Tks_Pt->size() << ", pixHits_XYZ = " << pixHits_XYZ.size() << ", pixTk_PixHits = " << pixTk_PixHits << ", pixTk_SharedPixHits = " << pixTk_SharedPixHits << endl;
	  // cout << "pixTk_PixHits = " << pixTk_SharedPixHits << endl;
	  // auxTools_.PrintVector(pixTk_PixHits_X);
	  // auxTools_.PrintVector(pixTk_PixHits_Y);
	  // auxTools_.PrintVector(pixTk_PixHits_Z);
	  // cout << "pixTk_SharedPixHits = " << pixTk_SharedPixHits << endl;
	  // auxTools_.PrintVector(pixTk_SharedPixHits_X);
	  // auxTools_.PrintVector(pixTk_SharedPixHits_Y);
	  // auxTools_.PrintVector(pixTk_SharedPixHits_Z);

	  // tp->PrintTPProperties(t);
	  // for (int i = 0; i < (int) L1PixTks_Pt->size(); i++) s->PrintPixTrackProperties(i, true);

	  // std::cout << "\n" << std::endl;
	  // auxTools_.PrintVector(pixTk_SharedPixHits_PixTkRefIndex);
	  // for (int i = 0; i < (int) pixTk_PixHits_PixTkRefIndex.size(); i++){
	  //   s->PrintPixTrackProperties(pixTk_PixHits_PixTkRefIndex.at(i), false);

	  // }

	  
	  
	  // s->PrintPixTrackProperties(pixTk_Index);
	  
	}
	   for(int i = 0; i < pixTk_SharedPixHits; i++){
	  h_pixTk_sharedPixHits_Rho   ->Fill( pixTk_SharedPixHits_R.at(i) );
	  h_pixTk_sharedPixHits_Z     ->Fill( pixTk_SharedPixHits_Z.at(i) );
	  h_pixTk_sharedPixHits_Type  ->Fill( pixTk_SharedPixHits_Type.at(i) );
	  h_pixTk_sharedPixHits_ZVsRho->Fill( pixTk_SharedPixHits_Z.at(i), pixTk_SharedPixHits_R.at(i) );
	  if (fabs(tk_Eta) < _eta_C) h_pixTk_sharedPixHits_ZVsRho_C->Fill( pixTk_SharedPixHits_Z.at(i), pixTk_SharedPixHits_R.at(i) );
	  else if (fabs(tk_Eta) < _eta_F && fabs(tk_Eta) >= _eta_C) h_pixTk_sharedPixHits_ZVsRho_I->Fill( pixTk_SharedPixHits_Z.at(i), pixTk_SharedPixHits_R.at(i) );
	  else if (fabs(tk_Eta) >= _eta_F) h_pixTk_sharedPixHits_ZVsRho_F->Fill( pixTk_SharedPixHits_Z.at(i), pixTk_SharedPixHits_R.at(i) );
	  else{}
	  if (fabs(tk_Eta) >= 2.0) h_pixTk_sharedPixHits_ZVsRho_EtaGE2->Fill( pixTk_SharedPixHits_Z.at(i), pixTk_SharedPixHits_R.at(i) );
	}// for(int i = 0; i < pixTk_SharedPixHits; i++){
	
	    
	// Candidate (Projected) Pixel Hits
	h_pixTk_candPixHits_N->Fill(pixTk_CandPixHits_N);
	for(int j = 0; j < pixTk_CandPixHits_N; j++){
	  h_pixTk_candPixHits_Rho   ->Fill( pixTk_CandPixHits_R.at(j) );
	  h_pixTk_candPixHits_Z     ->Fill( pixTk_CandPixHits_Z.at(j) );
	  h_pixTk_candPixHits_Type  ->Fill( pixTk_CandPixHits_Type.at(j) );
	  h_pixTk_candPixHits_ZVsRho->Fill( pixTk_CandPixHits_Z.at(j), pixTk_CandPixHits_R.at(j) );
	  if (fabs(tk_Eta) < _eta_C) h_pixTk_candPixHits_ZVsRho_C->Fill( pixTk_CandPixHits_Z.at(j), pixTk_CandPixHits_R.at(j) );
	  else if (fabs(tk_Eta) < _eta_F && fabs(tk_Eta) >= _eta_C) h_pixTk_candPixHits_ZVsRho_I->Fill( pixTk_CandPixHits_Z.at(j), pixTk_CandPixHits_R.at(j) );
	  else if (fabs(tk_Eta) >= _eta_F) h_pixTk_candPixHits_ZVsRho_F->Fill( pixTk_CandPixHits_Z.at(j), pixTk_CandPixHits_R.at(j) );
	  else{}
	  if (fabs(tk_Eta) >= 2.0) h_pixTk_candPixHits_ZVsRho_EtaGE2->Fill( pixTk_CandPixHits_Z.at(j), pixTk_CandPixHits_R.at(j) );
	}

      } // if (bDoPixelTracks){


      
      // Overflow bins for chi-square-related distributions
      if (tk_ChiSq > _chiSq_overflow)  tk_ChiSq = _chiSq_overflow - 0.01;
      if (tk_RedChiSq > _redChiSq_overflow) tk_RedChiSq = _redChiSq_overflow - 0.01;

      // Fill Histos	  
      h_2d_logchi2_eta    ->Fill( tp_Eta, log(tk_ChiSq)    );
      h_2d_logchi2_dof_eta->Fill( tp_Eta, log(tk_RedChiSq) );
      h_match_trk_chi2    ->Fill(tk_ChiSq   );
      h_match_trk_chi2_dof->Fill(tk_RedChiSq);

      // Central eta
      if (fabs(tk_Eta) < _eta_C) {
	// Low pT
	if (tk_Pt < _pt_L) {
	  h_match_trk_chi2_C_L    ->Fill(tk_ChiSq   );
          h_match_trk_chi2_dof_C_L->Fill(tk_RedChiSq);
	  h_match_trk_chi2_L    ->Fill(tk_ChiSq   );
          h_match_trk_chi2_dof_L->Fill(tk_RedChiSq);
	}

	// Medium pT
	else if (tk_Pt < _pt_H && tk_Pt >= _pt_L) {
	  h_match_trk_chi2_C_M    ->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_C_M->Fill(tk_RedChiSq);
	  h_match_trk_chi2_M    ->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_M->Fill(tk_RedChiSq);
	}

	// High pT
	else {
	  h_match_trk_chi2_C_H->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_C_H->Fill(tk_RedChiSq);
	  h_match_trk_chi2_H    ->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_H->Fill(tk_RedChiSq);
	}
      }
      // Intermediate eta
      else if (fabs(tk_Eta) < _eta_F && fabs(tk_Eta) >= _eta_C) {
	// Low pT
	if (tk_Pt < _pt_L) {
          h_match_trk_chi2_I_L    ->Fill(tk_ChiSq);
          h_match_trk_chi2_dof_I_L->Fill(tk_RedChiSq);
	  h_match_trk_chi2_L    ->Fill(tk_ChiSq   );
          h_match_trk_chi2_dof_L->Fill(tk_RedChiSq);	 
	}

	// Medium pT
	else if (tk_Pt < _pt_H && tk_Pt >= _pt_L) {
          h_match_trk_chi2_I_M->Fill(tk_ChiSq);
          h_match_trk_chi2_dof_I_M->Fill(tk_RedChiSq);
	  h_match_trk_chi2_M    ->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_M->Fill(tk_RedChiSq);
	}

	// High pT
	else {
          h_match_trk_chi2_I_H->Fill(tk_ChiSq);
          h_match_trk_chi2_dof_I_H->Fill(tk_RedChiSq);
	  h_match_trk_chi2_H    ->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_H->Fill(tk_RedChiSq);
	}
      }

      // Forward eta
      else if (fabs(tk_Eta) >= _eta_F) {
	// Low pT
	if (tk_Pt < _pt_L) {
          h_match_trk_chi2_F_L    ->Fill(tk_ChiSq   );
          h_match_trk_chi2_dof_F_L->Fill(tk_RedChiSq);
	  h_match_trk_chi2_L    ->Fill(tk_ChiSq   );
          h_match_trk_chi2_dof_L->Fill(tk_RedChiSq);
	}
	// Medium pT
	else if (tk_Pt < _pt_H && tk_Pt >= _pt_L) {
          h_match_trk_chi2_F_M    ->Fill(tk_ChiSq);
          h_match_trk_chi2_dof_F_M->Fill(tk_RedChiSq);
	  h_match_trk_chi2_M    ->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_M->Fill(tk_RedChiSq);
       	}
	// High pT
	else {
          h_match_trk_chi2_F_H    ->Fill(tk_ChiSq);
          h_match_trk_chi2_dof_F_H->Fill(tk_RedChiSq);
	  h_match_trk_chi2_H    ->Fill(tk_ChiSq);
	  h_match_trk_chi2_dof_H->Fill(tk_RedChiSq);
	}
      }

      // Cut on chi2
      if (tk_ChiSq > _chiSq_overflow) exit(1);
      
      // Fill matched track histograms
      h_match_tp_pt->Fill(tp_Pt);
      if (tp_Pt < _pt_L) h_match_tp_pt_L->Fill(tp_Pt);
      
      if(tp_Pt > _pt_tkMin) {
	h_match_tp_eta->Fill(tp_Eta);
	h_match_tp_phi->Fill(tp_Phi);
	h_match_tp_z0 ->Fill(tp_z0);
	h_match_tp_d0 ->Fill(tp_d0 * _fromCentimetresToMicrometers);
	
	if ( IsWithinEtaRegion("Central"          , tp_Eta) ){
	  n_match_eta_C++;
	  h_match_tp_pt_C->Fill(tp_Pt);
	}
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ){
	  n_match_eta_I++;
	  h_match_tp_pt_I->Fill(tp_Pt);
	}
	else if ( IsWithinEtaRegion("Forward"     , tp_Eta) ){
	  n_match_eta_F++;
	  h_match_tp_pt_F->Fill(tp_Pt);
	}
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}
      }

      // Low pT
      if ( IsWithinPtRange("Low", tp_Pt) ){
	if ( IsWithinEtaRegion("Central"          , tp_Eta) ) n_match_eta_C_pt_L++;
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ) n_match_eta_I_pt_L++;
	else if ( IsWithinEtaRegion("Forward"     , tp_Eta) ) n_match_eta_F_pt_L++;
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}

	h_match_tp_eta_L->Fill(tp_Eta);
      }

      // Medium pT
      if ( IsWithinPtRange("Middle", tp_Pt) ){
	if ( IsWithinEtaRegion("Central"          , tp_Eta) ) n_match_eta_C_pt_M++;
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ) n_match_eta_I_pt_M++;
	else if ( IsWithinEtaRegion("Forward"     , tp_Eta) ) n_match_eta_F_pt_M++;
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}

	h_match_tp_eta_M->Fill(tp_Eta);
      }
      // High pT
      if ( IsWithinPtRange("High", tp_Pt) ){
	if ( IsWithinEtaRegion("Central"          , tp_Eta) ) n_match_eta_C_pt_H++;
	else if ( IsWithinEtaRegion("Intermediate", tp_Eta) ) n_match_eta_I_pt_H++;
	else if ( IsWithinEtaRegion("Forward"     , tp_Eta) ) n_match_eta_F_pt_H++;
	else{
	  cout << "E R R O R ! Tracking::Loop(...) - Unexpected Eta value of \"" << tp_Eta << "\". EXIT" << endl;
	  exit(1);
	}

	h_match_tp_eta_H->Fill(tp_Eta);
      }

      // Fill nstub histograms
      h_match_trk_nstub->Fill(tk_nStubs);
      if (fabs(tk_Eta) < _eta_C) h_match_trk_nstub_C->Fill(tk_nStubs);
      else if (fabs(tk_Eta) < _eta_F && fabs(tk_Eta) >= _eta_C) h_match_trk_nstub_I->Fill(tk_nStubs);
      else if (fabs(tk_Eta) >= _eta_F) h_match_trk_nstub_F->Fill(tk_nStubs);

      // Fill residual histograms
      h_res_pt   ->Fill( tk_Pt  - tp_Pt);
      h_res_ptRel->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_res_eta  ->Fill( tk_Eta - tp_Eta);
      h_res_phi  ->Fill( tk_Phi - tp_Phi);
      h_res_z0   ->Fill( tk_z0  - tp_z0);
      h_res_d0   ->Fill( tk_d0  - tp_d0);
      // Central Region
      if ( IsWithinEtaRegion("Central", tk_Eta) ){
	h_res_pt_C   ->Fill( tk_Pt  - tp_Pt);
	h_res_ptRel_C->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	h_res_eta_C  ->Fill( tk_Eta - tp_Eta);
	h_res_phi_C  ->Fill( tk_Phi - tp_Phi);
	h_res_z0_C   ->Fill( tk_z0 - tp_z0);
	h_res_d0_C   ->Fill( tk_d0 - tp_d0);
      }
      // Intermediate Region
      else if ( IsWithinEtaRegion("Intermediate", tk_Eta) ){
	h_res_pt_I   ->Fill( tk_Pt  - tp_Pt);
	h_res_ptRel_I->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	h_res_eta_I  ->Fill( tk_Eta - tp_Eta);
	h_res_phi_I  ->Fill( tk_Phi - tp_Phi);
	h_res_z0_I   ->Fill( tk_z0 - tp_z0);
	h_res_d0_I   ->Fill( tk_d0 - tp_d0);
      }
      // Forward Region
      else if ( IsWithinEtaRegion("Forward", tk_Eta) ){
	h_res_pt_F   ->Fill( tk_Pt  - tp_Pt);
	h_res_ptRel_F->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	h_res_eta_F  ->Fill( tk_Eta - tp_Eta);
	h_res_phi_F  ->Fill( tk_Phi - tp_Phi);
	h_res_z0_F   ->Fill( tk_z0 - tp_z0);
	h_res_d0_F   ->Fill( tk_d0 - tp_d0);
      }
      else{
	cout << "E R R O R ! Tracking::Loop(...) - Unexpected Track-Eta of \"" << tk_Eta << "\". EXIT" << endl;
	exit(1);
      }


      // For-loop: 0-nRange 
      for (int j=0; j < PT_BINS; j++) {

	// Fill residual vs. pt histograms in pT-steps of 5 GeV/c
	double ptLow  = (double)j * PT_BINWIDTH;
	double ptHigh = ptLow + 5.0;
	
	if ( (tp_Pt >= ptLow) && (tp_Pt < ptHigh) ) {
	  h_resVsPt_pt[j]   ->Fill( tk_Pt  - tp_Pt);
	  h_resVsPt_ptRel[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	  h_resVsPt_eta[j]  ->Fill( tk_Eta - tp_Eta);
	  h_resVsPt_phi[j]  ->Fill( tk_Phi - tp_Phi);
	  h_resVsPt_z0[j]   ->Fill( tk_z0  - tp_z0);
	  h_resVsPt_d0[j]   ->Fill( tk_d0  - tp_d0);

	  // Central Region
	  if (fabs(tk_Eta) < _eta_C) {
	    h_resVsPt_pt_C[j]   ->Fill( tk_Pt  - tp_Pt);
	    h_resVsPt_ptRel_C[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsPt_eta_C[j]  ->Fill( tk_Eta - tp_Eta);
	    h_resVsPt_phi_C[j]  ->Fill( tk_Phi - tp_Phi);
	    h_resVsPt_z0_C[j]   ->Fill( tk_z0  - tp_z0);
	    h_resVsPt_d0_C[j]   ->Fill( tk_d0  - tp_d0);
	  }
	  // Intermedieate Region
	  else if (fabs(tk_Eta) < _eta_F && fabs(tk_Eta) >= _eta_C) {
	    h_resVsPt_pt_I[j]   ->Fill( tk_Pt  - tp_Pt);
	    h_resVsPt_ptRel_I[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsPt_eta_I[j]  ->Fill( tk_Eta - tp_Eta);
	    h_resVsPt_phi_I[j]  ->Fill( tk_Phi - tp_Phi);
	    h_resVsPt_z0_I[j]   ->Fill( tk_z0  - tp_z0);
	    h_resVsPt_d0_I[j]   ->Fill( tk_d0  - tp_d0);
	  }
	  // Forward Region
	  else if (fabs(tk_Eta) >= _eta_F) {
	    h_resVsPt_pt_F[j]   ->Fill( tk_Pt  - tp_Pt);
	    h_resVsPt_ptRel_F[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsPt_eta_F[j]  ->Fill( tk_Eta - tp_Eta);
	    h_resVsPt_phi_F[j]  ->Fill( tk_Phi - tp_Phi);
	    h_resVsPt_z0_F[j]   ->Fill( tk_z0  - tp_z0);
	    h_resVsPt_d0_F[j]   ->Fill( tk_d0  - tp_d0);
	  }
	  
	} // if ( (tp_Pt >= ptLow) && (tp_Pt < ptHigh) ) {
      } // for (int j=0; j < PT_BINS; j++) {



      // For-loop: 0-nEtaRange 
      for (int k=0; k < ETA_BINS; k++) {

	// Fill residual vs. eta histogram in eta-steps of DeltaEta = 0.1
	double etaLow  = (double)k * ETA_BINWIDTH;
	double etaHigh = etaLow + ETA_BINWIDTH;

	//  eta_low <= |eta | < eta_high
	if ( (fabs(tp_Eta) >= etaLow) && (fabs(tp_Eta) < etaHigh) ) {
	  h_resVsEta_pt[k]   ->Fill( tk_Pt  - tp_Pt);
	  h_resVsEta_ptRel[k]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	  h_resVsEta_eta[k]  ->Fill( tk_Eta - tp_Eta);
	  h_resVsEta_phi[k]  ->Fill( tk_Phi - tp_Phi);
	  h_resVsEta_z0[k]   ->Fill( tk_z0  - tp_z0);
	  h_resVsEta_d0[k]   ->Fill( tk_d0  - tp_d0);

	  if ( IsWithinPtRange("Low", tk_Pt) ){
	    h_resVsEta_pt_L[k]   ->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsEta_ptRel_L[k]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsEta_eta_L[k]  ->Fill( tk_Eta - tp_Eta);
	    h_resVsEta_phi_L[k]  ->Fill( tk_Phi - tp_Phi);
	    h_resVsEta_z0_L[k]   ->Fill( tk_z0  - tp_z0);
	    h_resVsEta_d0_L[k]   ->Fill( tk_d0  - tp_d0);
	  }
	  else if ( IsWithinPtRange("Middle", tk_Pt) ){
	    h_resVsEta_pt_M[k]   ->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsEta_ptRel_M[k]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsEta_eta_M[k]  ->Fill( tk_Eta - tp_Eta);
	    h_resVsEta_phi_M[k]  ->Fill( tk_Phi - tp_Phi);
	    h_resVsEta_z0_M[k]   ->Fill( tk_z0  - tp_z0);
	    h_resVsEta_d0_M[k]   ->Fill( tk_d0  - tp_d0);
	  }
	  else if ( IsWithinPtRange("High", tk_Pt) ){
	    h_resVsEta_pt_H[k]   ->Fill((tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsEta_ptRel_H[k]->Fill((tk_Pt - tp_Pt)/tp_Pt);
	    h_resVsEta_eta_H[k]  ->Fill(tk_Eta - tp_Eta);
	    h_resVsEta_phi_H[k]  ->Fill(tk_Phi - tp_Phi);
	    h_resVsEta_z0_H[k]   ->Fill(tk_z0 - tp_z0);
	    h_resVsEta_d0_H[k]   ->Fill(tk_d0 - tp_d0);
	  }
	  else{
	    cout << "E R R O R ! Tracking::Loop(...) - Unexpected Track-Pt of \"" << tk_Pt << "\". EXIT" << endl;
	    exit(1);
	  }
	  
	} // if ( (fabs(tp_Eta) >= etaLow) && (fabs(tp_Eta) < etaHigh) ) {
      } // for (int k=0; k < ETA_BINS; k++) {
      
      // Fill 2D histograms
      h_2d_dpt_eta   ->Fill(tp_Eta, fabs((tk_Pt - tp_Pt)));
      h_2d_dptRel_eta->Fill(tp_Eta, fabs((tk_Pt - tp_Pt)/tp_Pt));
      h_2d_deta_eta  ->Fill(tp_Eta, fabs(tk_Eta - tp_Eta));
      h_2d_dphi_eta  ->Fill(tp_Eta, fabs(tk_Phi - tp_Phi));
      h_2d_dz0_eta   ->Fill(tp_Eta, fabs(tk_z0  - tp_z0));
      h_2d_dd0_eta   ->Fill(tp_Eta, fabs(tk_d0  - tp_d0));    

#ifdef bDebug	       
      tp->PrintTPProperties(t);
#endif
      
      } // For-loop: TPs

#ifndef bDebug
  auxTools_.ProgressBar(jentry, nEntries, NPBarDivisions, PBarWidth);
#endif
  }// For-loop: Entries


  // Fill efficiency histograms
  MakeEfficiencyHisto(h_eff_pt   , h_match_tp_pt   , h_tp_pt   );
  MakeEfficiencyHisto(h_eff_pt_L , h_match_tp_pt_L , h_tp_pt_L );
  MakeEfficiencyHisto(h_eff_pt_C , h_match_tp_pt_C , h_tp_pt_C );
  MakeEfficiencyHisto(h_eff_pt_I , h_match_tp_pt_I , h_tp_pt_I );
  MakeEfficiencyHisto(h_eff_pt_F , h_match_tp_pt_F , h_tp_pt_F );
  MakeEfficiencyHisto(h_eff_eta  , h_match_tp_eta  , h_tp_eta  );
  MakeEfficiencyHisto(h_eff_eta_L, h_match_tp_eta_L, h_tp_eta_L);
  MakeEfficiencyHisto(h_eff_eta_M, h_match_tp_eta_M, h_tp_eta_M);
  MakeEfficiencyHisto(h_eff_eta_H, h_match_tp_eta_H, h_tp_eta_H);
  MakeEfficiencyHisto(h_eff_phi  , h_match_tp_phi  , h_tp_phi  );
  MakeEfficiencyHisto(h_eff_z0   , h_match_tp_z0   , h_tp_z0   );
  MakeEfficiencyHisto(h_eff_d0   , h_match_tp_d0   , h_tp_d0   );

  FinaliseHistos_();
  PrintEfficiencies_();
  PrintResolutions_();

  // Write the histograms to the file
  outFile->cd();
  outFile->Write();

}


//****************************************************************************
void Tracking::BookHistos_(void)
//****************************************************************************
{

  // const int nEtaBins = 52;
  const int nEtaBins = 26;
  histoTools_.BookHisto_1D( h_tp_pt   , "tp_pt"   ,   100,   0.0,  +100.0 );
  histoTools_.BookHisto_1D( h_tp_pt_L , "tp_pt_L" ,    50,   0.0,    +5.0 );
  histoTools_.BookHisto_1D( h_tp_pt_C , "tp_pt_C" ,   100,   0.0,  +100.0 );
  histoTools_.BookHisto_1D( h_tp_pt_I , "tp_pt_I" ,   100,   0.0,  +100.0 );
  histoTools_.BookHisto_1D( h_tp_pt_F , "tp_pt_F" ,   100,   0.0,  +100.0 );
  histoTools_.BookHisto_1D( h_tp_eta  , "tp_eta"  ,    26,  -2.6,    +2.6 );
  histoTools_.BookHisto_1D( h_tp_eta_L, "tp_eta_L",    26,  -2.6,    +2.6 );
  histoTools_.BookHisto_1D( h_tp_eta_M, "tp_eta_M",    26,  -2.6,    +2.6 );
  histoTools_.BookHisto_1D( h_tp_eta_H, "tp_eta_H",    26,  -2.6,    +2.6 );
  histoTools_.BookHisto_1D( h_tp_phi  , "tp_phi"  ,    64,  -3.2,    +3.2 );
  histoTools_.BookHisto_1D( h_tp_z0   , "tp_z0"   ,    80, -20.0,   +20.0 );
  histoTools_.BookHisto_1D( h_tp_d0   , "tp_d0"   ,   100, -50.0,   +50.0 ); //micrometers

  histoTools_.BookHisto_1D( h_match_tp_pt   , "match_tp_pt"   , 100,   0.0, +100.0 );
  histoTools_.BookHisto_1D( h_match_tp_pt_L , "match_tp_pt_L" ,  50,   0.0,   +5.0 );
  histoTools_.BookHisto_1D( h_match_tp_pt_C , "match_tp_pt_C" , 100,   0.0, +100.0 );
  histoTools_.BookHisto_1D( h_match_tp_pt_I , "match_tp_pt_I" , 100,   0.0, +100.0 );
  histoTools_.BookHisto_1D( h_match_tp_pt_F , "match_tp_pt_F" , 100,   0.0, +100.0 );
  histoTools_.BookHisto_1D( h_match_tp_eta  , "match_tp_eta"  ,  26,  -2.6,   +2.6 );
  histoTools_.BookHisto_1D( h_match_tp_eta_L, "match_tp_eta_L",  26,  -2.6,   +2.6 );
  histoTools_.BookHisto_1D( h_match_tp_eta_M, "match_tp_eta_M",  26,  -2.6,   +2.6 );
  histoTools_.BookHisto_1D( h_match_tp_eta_H, "match_tp_eta_H",  26,  -2.6,   +2.6 );
  histoTools_.BookHisto_1D( h_match_tp_phi  , "match_tp_phi"  ,  64,  -3.2,   +3.2 );
  histoTools_.BookHisto_1D( h_match_tp_z0   , "match_tp_z0"   ,  80, -20.0,  +20.0 );
  histoTools_.BookHisto_1D( h_match_tp_d0   , "match_tp_d0"   , 100, -50.0,  +50.0 ); //micrometers

  // efficiency histograms
  histoTools_.BookHisto_1D( h_eff_pt   , "eff_pt"   , 100,   0.0,  +100.0  );
  histoTools_.BookHisto_1D( h_eff_pt_L , "eff_pt_L" ,  50,   0.0,    +5.0  );
  histoTools_.BookHisto_1D( h_eff_pt_C , "eff_pt_C" , 100,   0.0,  +100.0  );
  histoTools_.BookHisto_1D( h_eff_pt_I , "eff_pt_I" , 100,   0.0,  +100.0  ); 
  histoTools_.BookHisto_1D( h_eff_pt_F , "eff_pt_F" , 100,   0.0,  +100.0  );
  histoTools_.BookHisto_1D( h_eff_eta  , "eff_eta"  ,  26,  -2.6,    +2.6  );
  histoTools_.BookHisto_1D( h_eff_eta_L, "eff_eta_L",  26,  -2.6,    +2.6  );
  histoTools_.BookHisto_1D( h_eff_eta_M, "eff_eta_M",  26,  -2.6,    +2.6  );
  histoTools_.BookHisto_1D( h_eff_eta_H, "eff_eta_H",  26,  -2.6,    +2.6  );
  histoTools_.BookHisto_1D( h_eff_phi  , "eff_phi"  ,  64,  -3.2,    +3.2  );
  histoTools_.BookHisto_1D( h_eff_z0   , "eff_z0"   ,  80, -20.0,   +20.0  );
  histoTools_.BookHisto_1D( h_eff_d0   , "eff_d0"   , 100,-50.0,   +50.0  );

  histoTools_.BookHisto_1D( h_match_trk_nstub  , "match_trk_nstub"  , 15, -0.5, 14.5);
  histoTools_.BookHisto_1D( h_match_trk_nstub_C, "match_trk_nstub_C", 15, -0.5, 14.5);
  histoTools_.BookHisto_1D( h_match_trk_nstub_I, "match_trk_nstub_I", 15, -0.5, 14.5);
  histoTools_.BookHisto_1D( h_match_trk_nstub_F, "match_trk_nstub_F", 15, -0.5, 14.5);

  // chi2 histograms (last bin is an overflow bin)
  histoTools_.BookHisto_1D( h_match_trk_chi2     , "match_trk_chi2"    , 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_L   , "match_trk_chi2_L"  , 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_M   , "match_trk_chi2_M"  , 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_H   , "match_trk_chi2_H"  , 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_C_L , "match_trk_chi2_C_L", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_I_L , "match_trk_chi2_I_L", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_F_L , "match_trk_chi2_F_L", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_C_M , "match_trk_chi2_C_M", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_I_M , "match_trk_chi2_I_M", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_F_M , "match_trk_chi2_F_M", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_C_H , "match_trk_chi2_C_H", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_I_H , "match_trk_chi2_I_H", 100, 0, 100);
  histoTools_.BookHisto_1D( h_match_trk_chi2_F_H , "match_trk_chi2_F_H", 100, 0, 100);

  // chi2/dof histograms (lastbin is an overflow bin)
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof    , "match_trk_chi2_dof"    , 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_L  , "match_trk_chi2_dof_L"  , 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_M  , "match_trk_chi2_dof_M"  , 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_H  , "match_trk_chi2_dof_H"  , 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_C_L, "match_trk_chi2_dof_C_L", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_I_L, "match_trk_chi2_dof_I_L", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_F_L, "match_trk_chi2_dof_F_L", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_C_M, "match_trk_chi2_dof_C_M", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_I_M, "match_trk_chi2_dof_I_M", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_F_M, "match_trk_chi2_dof_F_M", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_C_H, "match_trk_chi2_dof_C_H", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_I_H, "match_trk_chi2_dof_I_H", 150, 0, 15);
  histoTools_.BookHisto_1D( h_match_trk_chi2_dof_F_H, "match_trk_chi2_dof_F_H", 150, 0, 15);

  // resolution histograms
  histoTools_.BookHisto_1D( h_res_pt     , "res_pt"     , 100, -5.0  , 5.0  );
  histoTools_.BookHisto_1D( h_res_pt_C   , "res_pt_C"   , 100, -5.0  , 5.0  );
  histoTools_.BookHisto_1D( h_res_pt_I   , "res_pt_I"   , 100, -5.0  , 5.0  );
  histoTools_.BookHisto_1D( h_res_pt_F   , "res_pt_F"   , 100, -5.0  , 5.0  );

  histoTools_.BookHisto_1D( h_res_ptRel  , "res_ptRel"  , 1000, -1.0  , 1.0  );
  histoTools_.BookHisto_1D( h_res_ptRel_C, "res_ptRel_C", 1000, -1.0  , 1.0  );
  histoTools_.BookHisto_1D( h_res_ptRel_I, "res_ptRel_I", 1000, -1.0  , 1.0  );
  histoTools_.BookHisto_1D( h_res_ptRel_F, "res_ptRel_F", 1000, -1.0  , 1.0  );

  histoTools_.BookHisto_1D( h_res_eta    , "res_eta"    , 400, -0.01 , 0.01 );
  histoTools_.BookHisto_1D( h_res_eta_C  , "res_eta_C"  , 200, -0.01 , 0.01 );
  histoTools_.BookHisto_1D( h_res_eta_I  , "res_eta_I"  , 200, -0.01 , 0.01 );
  histoTools_.BookHisto_1D( h_res_eta_F  , "res_eta_F"  , 200, -0.01 , 0.01 );

  histoTools_.BookHisto_1D( h_res_phi    , "res_phi"    , 100, -0.005, 0.005);
  histoTools_.BookHisto_1D( h_res_phi_C  , "res_phi_C"  , 100, -0.005, 0.005);
  histoTools_.BookHisto_1D( h_res_phi_I  , "res_phi_I"  , 100, -0.005, 0.005);
  histoTools_.BookHisto_1D( h_res_phi_F  , "res_phi_F"  , 100, -0.005, 0.005);

  histoTools_.BookHisto_1D( h_res_z0     , "res_z0"     , 1000, -1.0  , 1.0  );
  histoTools_.BookHisto_1D( h_res_z0_C   , "res_z0_C"   , 1000, -1.0  , 1.0  );
  histoTools_.BookHisto_1D( h_res_z0_I   , "res_z0_I"   , 1000, -1.0  , 1.0  );
  histoTools_.BookHisto_1D( h_res_z0_F   , "res_z0_F"   , 1000, -1.0  , 1.0  );

  histoTools_.BookHisto_1D( h_res_d0     , "res_d0"     , 500, -0.5  , 0.5  );
  histoTools_.BookHisto_1D( h_res_d0_C   , "res_d0_C"   , 500, -0.5  , 0.5  );
  histoTools_.BookHisto_1D( h_res_d0_I   , "res_d0_I"   , 500, -0.5  , 0.5  );
  histoTools_.BookHisto_1D( h_res_d0_F   , "res_d0_F"   , 500, -0.5  , 0.5  );

  // resolution vs. pt histograms  
  const int nPtRange = 20;
  TString ptrange[nPtRange] = {"0-5","5-10", "10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50",
			       "50-55", "55-60","60-65","65-70","70-75","75-80","80-85","85-90","90-95","95-100"};

  // For-loop: All pt bins
  for (int i=0; i < PT_BINS; i++) {

    double pt = PT_BINWIDTH*(i);
    int nBinsXFactor   = 1;
    int nBinsXFactor_C = 2;
    int nBinsXFactor_I = 2;
    int nBinsXFactor_F = 2;

    histoTools_.BookHisto_1D( h_resVsPt_pt[i]  , "resVsPt_pt_"   + ptrange[i], 100, -5.0, 5.0);
    histoTools_.BookHisto_1D( h_resVsPt_pt_C[i], "resVsPt_pt_C_" + ptrange[i], 100, -5.0, 5.0);
    histoTools_.BookHisto_1D( h_resVsPt_pt_I[i], "resVsPt_pt_I_" + ptrange[i], 100, -5.0, 5.0);
    histoTools_.BookHisto_1D( h_resVsPt_pt_F[i], "resVsPt_pt_F_" + ptrange[i], 100, -5.0, 5.0);

    // restictive range: -0.15 to 0.15
    histoTools_.BookHisto_1D( h_resVsPt_ptRel[i]  , "resVsPt_ptRel_"  + ptrange[i], 400/nBinsXFactor, -0.50, +0.50);
    histoTools_.BookHisto_1D( h_resVsPt_ptRel_C[i], "resVsPt_ptRel_C_"+ ptrange[i], 200/nBinsXFactor, -0.50, +0.50);
    histoTools_.BookHisto_1D( h_resVsPt_ptRel_I[i], "resVsPt_ptRel_I_"+ ptrange[i], 200/nBinsXFactor, -0.50, +0.50);
    histoTools_.BookHisto_1D( h_resVsPt_ptRel_F[i], "resVsPt_ptRel_F_"+ ptrange[i], 200/nBinsXFactor, -0.50, +0.50);

    histoTools_.BookHisto_1D( h_resVsPt_eta[i]  , "resVsPt_eta_"    + ptrange[i], 400/nBinsXFactor, -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsPt_eta_C[i], "resVsPt_eta_C_"  + ptrange[i], 200/nBinsXFactor, -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsPt_eta_I[i], "resVsPt_eta_I_"  + ptrange[i], 200/nBinsXFactor, -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsPt_eta_F[i], "resVsPt_eta_F_"  + ptrange[i], 200/nBinsXFactor, -0.01, +0.01);

    histoTools_.BookHisto_1D( h_resVsPt_phi[i]  , "resVsPt_phi_"   + ptrange[i], 100, -0.005, +0.005);
    histoTools_.BookHisto_1D( h_resVsPt_phi_C[i], "resVsPt_phi_C_" + ptrange[i], 100, -0.005, +0.005);
    histoTools_.BookHisto_1D( h_resVsPt_phi_I[i], "resVsPt_phi_I_" + ptrange[i], 100, -0.005, +0.005);
    histoTools_.BookHisto_1D( h_resVsPt_phi_F[i], "resVsPt_phi_F_" + ptrange[i], 100, -0.005, +0.005);

    histoTools_.BookHisto_1D( h_resVsPt_z0[i]  , "resVsPt_z0_"   + ptrange[i], 500/nBinsXFactor, -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsPt_z0_C[i], "resVsPt_z0_C_" + ptrange[i], 500/nBinsXFactor, -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsPt_z0_I[i], "resVsPt_z0_I_" + ptrange[i], 500/nBinsXFactor, -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsPt_z0_F[i], "resVsPt_z0_F_" + ptrange[i], 250/nBinsXFactor, -0.5, +0.5);

    histoTools_.BookHisto_1D( h_resVsPt_d0[i]  , "resVsPt_d0_"   + ptrange[i], 200/nBinsXFactor, -0.1, +0.1);
    histoTools_.BookHisto_1D( h_resVsPt_d0_C[i], "resVsPt_d0_C_" + ptrange[i], 200/nBinsXFactor, -0.1, +0.1);
    histoTools_.BookHisto_1D( h_resVsPt_d0_I[i], "resVsPt_d0_I_" + ptrange[i], 200/nBinsXFactor, -0.1, +0.1);
    histoTools_.BookHisto_1D( h_resVsPt_d0_F[i], "resVsPt_d0_F_" + ptrange[i], 100/nBinsXFactor, -0.1, +0.1);
  }

  // resolution vs. eta histograms
  //  const int nEtaRange = 25;
  //  TString etarange[nEtaRange] = {"0.1", "0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0",
  //  				 "1.1", "1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2.0",
  //  				 "2.1", "2.2","2.3","2.4","2.5"};
  const int nEtaRange = 12;
  TString etarange[nEtaRange] = {"0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4"};
 
  for (int i=0; i < ETA_BINS; i++) {  

    double eta = ETA_BINWIDTH*(i+1);
    int nBinsXFactor   = 2;
    int nBinsXFactor_L = 4;
    int nBinsXFactor_M = 2;
    int nBinsXFactor_H = 2;
    // Increase bin-width at large eta (last bin)
    if (eta >= 2.2){
      nBinsXFactor   = nBinsXFactor   * 2;
      nBinsXFactor_L = nBinsXFactor_L * 2;
      nBinsXFactor_M = nBinsXFactor_M * 2;
      nBinsXFactor_H = nBinsXFactor_H * 2;
    }
    if (eta >= 2.4){
      nBinsXFactor   = nBinsXFactor   * 2;
      // nBinsXFactor_L = nBinsXFactor_L * 2;
      nBinsXFactor_M = nBinsXFactor_M * 2;
      nBinsXFactor_H = nBinsXFactor_H * 2;
    }
	
    histoTools_.BookHisto_1D( h_resVsEta_pt[i]  , "h_resVsEta_pt_"   + etarange[i], 100, -5.00, +5.00);
    histoTools_.BookHisto_1D( h_resVsEta_pt_L[i], "h_resVsEta_pt_L_" + etarange[i], 100, -0.25, +0.25);
    histoTools_.BookHisto_1D( h_resVsEta_pt_M[i], "h_resVsEta_pt_M_" + etarange[i], 100, -0.25, +0.25);
    histoTools_.BookHisto_1D( h_resVsEta_pt_H[i], "h_resVsEta_pt_H_" + etarange[i], 100, -0.25, +0.25);

    histoTools_.BookHisto_1D( h_resVsEta_ptRel[i]  , "h_resVsEta_ptRel_"  + etarange[i], 200, -0.50, +0.50);
    histoTools_.BookHisto_1D( h_resVsEta_ptRel_L[i], "h_resVsEta_ptRel_L_"+ etarange[i], 200, -0.50, +0.50);
    histoTools_.BookHisto_1D( h_resVsEta_ptRel_M[i], "h_resVsEta_ptRel_M_"+ etarange[i], 200, -0.50, +0.50);
    histoTools_.BookHisto_1D( h_resVsEta_ptRel_H[i], "h_resVsEta_ptRel_H_"+ etarange[i], 200, -0.50, +0.50);

    histoTools_.BookHisto_1D( h_resVsEta_eta[i]  , "h_resVsEta_eta_"   + etarange[i], 400/nBinsXFactor  , -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsEta_eta_L[i], "h_resVsEta_eta_L_" + etarange[i], 200/nBinsXFactor_L, -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsEta_eta_M[i], "h_resVsEta_eta_M_" + etarange[i], 400/nBinsXFactor_M, -0.01, +0.01); 
    histoTools_.BookHisto_1D( h_resVsEta_eta_H[i], "h_resVsEta_eta_H_" + etarange[i], 400/nBinsXFactor_H, -0.01, +0.01);

    histoTools_.BookHisto_1D( h_resVsEta_phi[i]  , "h_resVsEta_phi_"   + etarange[i], 100, -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsEta_phi_L[i], "h_resVsEta_phi_L_" + etarange[i], 100, -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsEta_phi_M[i], "h_resVsEta_phi_M_" + etarange[i], 100, -0.01, +0.01);
    histoTools_.BookHisto_1D( h_resVsEta_phi_H[i], "h_resVsEta_phi_H_" + etarange[i], 100, -0.01, +0.01);

    histoTools_.BookHisto_1D( h_resVsEta_z0[i]  , "h_resVsEta_z0_"   + etarange[i], 1000/nBinsXFactor  , -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsEta_z0_L[i], "h_resVsEta_z0_L_" + etarange[i], 1000/nBinsXFactor_L, -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsEta_z0_M[i], "h_resVsEta_z0_M_" + etarange[i], 1000/nBinsXFactor_M, -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsEta_z0_H[i], "h_resVsEta_z0_H_" + etarange[i], 1000/nBinsXFactor_H, -0.5, +0.5);

    histoTools_.BookHisto_1D( h_resVsEta_d0[i]  , "h_resVsEta_d0_"   + etarange[i], 2000/nBinsXFactor  , -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsEta_d0_L[i], "h_resVsEta_d0_L_" + etarange[i], 1000/nBinsXFactor_L, -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsEta_d0_M[i], "h_resVsEta_d0_M_" + etarange[i], 1000/nBinsXFactor_M, -0.5, +0.5);
    histoTools_.BookHisto_1D( h_resVsEta_d0_H[i], "h_resVsEta_d0_H_" + etarange[i], 2000/nBinsXFactor_H, -0.5, +0.5);

  }

  // 2D histograms
  histoTools_.BookHisto_2D( h_2d_logchi2_eta    , "2d_logchi2_eta"    , 50, -2.5, +2.5, 120, -4.0, +8.0  );
  histoTools_.BookHisto_2D( h_2d_logchi2_dof_eta, "2d_logchi2_dof_eta", 50, -2.5, +2.5, 120, -4.0, +8.0  );
  histoTools_.BookHisto_2D( h_2d_dz0_eta        , "2d_dz0_eta"        , 50, -2.5, +2.5, 120,  0.0, +1.2  );
  histoTools_.BookHisto_2D( h_2d_dd0_eta        , "2d_dd0_eta"        , 50, -2.5, +2.5, 120,  0.0, +1.2  );
  histoTools_.BookHisto_2D( h_2d_deta_eta       , "2d_deta_eta"       , 50, -2.5, +2.5, 120,  0.0, +0.012);
  histoTools_.BookHisto_2D( h_2d_dphi_eta       , "2d_dphi_eta"       , 50, -2.5, +2.5, 100,  0.0, +0.010);
  histoTools_.BookHisto_2D( h_2d_dpt_eta        , "2d_dpt_eta"        , 50, -2.5, +2.5, 500,  0.0, +5.0  );
  histoTools_.BookHisto_2D( h_2d_dptRel_eta     , "2d_dptRel_eta"     , 50, -2.5, +2.5, 500,  0.0, +5.0  );

  // resolution vs. pT histograms
  histoTools_.BookHisto_1D( h2_resVsPt_pt  , "resVsPt_pt"  , PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_pt_C, "resVsPt_pt_C", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_pt_I, "resVsPt_pt_I", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_pt_F, "resVsPt_pt_F", PT_BINS, 0.0, PT_MAX);

  histoTools_.BookHisto_1D( h2_resVsPt_ptRel  , "resVsPt_ptRel"  , PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_ptRel_C, "resVsPt_ptRel_C", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_ptRel_I, "resVsPt_ptRel_I", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_ptRel_F, "resVsPt_ptRel_F", PT_BINS, 0.0, PT_MAX);

  histoTools_.BookHisto_1D( h2_resVsPt_eta  , "resVsPt_eta"  , PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_eta_C, "resVsPt_eta_C", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_eta_I, "resVsPt_eta_I", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_eta_F, "resVsPt_eta_F", PT_BINS, 0.0, PT_MAX);

  histoTools_.BookHisto_1D( h2_resVsPt_phi  , "resVsPt_phi"  , PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_phi_C, "resVsPt_phi_C", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_phi_I, "resVsPt_phi_I", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_phi_F, "resVsPt_phi_F", PT_BINS, 0.0, PT_MAX);

  histoTools_.BookHisto_1D( h2_resVsPt_z0  , "resVsPt_z0"  , PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_z0_C, "resVsPt_z0_C", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_z0_I, "resVsPt_z0_I", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_z0_F, "resVsPt_z0_F", PT_BINS, 0.0, PT_MAX);

  histoTools_.BookHisto_1D( h2_resVsPt_d0  , "resVsPt_d0"  , PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_d0_C, "resVsPt_d0_C", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_d0_I, "resVsPt_d0_I", PT_BINS, 0.0, PT_MAX);
  histoTools_.BookHisto_1D( h2_resVsPt_d0_F, "resVsPt_d0_F", PT_BINS, 0.0, PT_MAX);

  // resolution vs. eta histograms
  histoTools_.BookHisto_1D( h2_resVsEta_pt   , "resVsEta_pt"   , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_pt_L , "resVsEta_pt_L" , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_pt_M , "resVsEta_pt_M" , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_pt_H , "resVsEta_pt_H" , ETA_BINS, 0.0, ETA_MAX);

  histoTools_.BookHisto_1D( h2_resVsEta_ptRel  , "resVsEta_ptRel"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_ptRel_L, "resVsEta_ptRel_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_ptRel_M, "resVsEta_ptRel_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_ptRel_H, "resVsEta_ptRel_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools_.BookHisto_1D( h2_resVsEta_eta  , "resVsEta_eta"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_eta_L, "resVsEta_eta_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_eta_M, "resVsEta_eta_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_eta_H, "resVsEta_eta_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools_.BookHisto_1D( h2_resVsEta_phi  , "resVsEta_phi"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_phi_L, "resVsEta_phi_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_phi_M, "resVsEta_phi_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_phi_H, "resVsEta_phi_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools_.BookHisto_1D( h2_resVsEta_z0  , "resVsEta_z0"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_z0_L, "resVsEta_z0_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_z0_M, "resVsEta_z0_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_z0_H, "resVsEta_z0_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools_.BookHisto_1D( h2_resVsEta_d0  , "resVsEta_d0"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_d0_L, "resVsEta_d0_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_d0_M, "resVsEta_d0_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools_.BookHisto_1D( h2_resVsEta_d0_H, "resVsEta_d0_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools_.BookHisto_1D( h_pixTk_pixHits_N            , "pixTk_pixHits_N"            ,   30,   -0.5,  +29.5 );
  histoTools_.BookHisto_1D( h_pixTk_pixHits_Rho          , "pixTk_pixHits_Rho"          , 1000,   +0.0,  +20.0 );
  histoTools_.BookHisto_1D( h_pixTk_pixHits_Z            , "pixTk_pixHits_Z"            ,  240,  -60.0,  +60.0 );
  histoTools_.BookHisto_1D( h_pixTk_pixHits_Type         , "pixTk_pixHits_Type"         ,   10,   -4.5,   +5.5 );
  histoTools_.BookHisto_2D( h_pixTk_pixHits_ZVsRho       , "pixTk_pixHits_ZVsRho"       ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_pixHits_ZVsRho_C     , "pixTk_pixHits_ZVsRho_C"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_pixHits_ZVsRho_I     , "pixTk_pixHits_ZVsRho_I"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_pixHits_ZVsRho_F     , "pixTk_pixHits_ZVsRho_F"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_pixHits_ZVsRho_EtaGE2, "pixTk_pixHits_ZVsRho_EtaGE2",  600,  -60.0,  +60.0, 100, +0.0, +20.0 );

  histoTools_.BookHisto_1D( h_pixTk_sharedPixHits_N            , "pixTk_sharedPixHits_N"            ,   30,   -0.5,  +29.5 );
  histoTools_.BookHisto_1D( h_pixTk_sharedPixHits_Rho          , "pixTk_sharedPixHits_Rho"          , 1000,   +0.0,  +20.0 );
  histoTools_.BookHisto_1D( h_pixTk_sharedPixHits_Z            , "pixTk_sharedPixHits_Z"            ,  240,  -60.0,  +60.0 );
  histoTools_.BookHisto_1D( h_pixTk_sharedPixHits_Type         , "pixTk_sharedPixHits_Type"         ,   10,   -4.5,   +5.5 );
  histoTools_.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho       , "pixTk_sharedPixHits_ZVsRho"       ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_C     , "pixTk_sharedPixHits_ZVsRho_C"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_I     , "pixTk_sharedPixHits_ZVsRho_I"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_F     , "pixTk_sharedPixHits_ZVsRho_F"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_EtaGE2, "pixTk_sharedPixHits_ZVsRho_EtaGE2",  600,  -60.0,  +60.0, 100, +0.0, +20.0 );

  histoTools_.BookHisto_1D( h_pixTk_candPixHits_N            , "pixTk_candPixHits_N"            ,   30,   -0.5,  +29.5 );
  histoTools_.BookHisto_1D( h_pixTk_candPixHits_Rho          , "pixTk_candPixHits_Rho"          , 1000,   +0.0,  +20.0 );
  histoTools_.BookHisto_1D( h_pixTk_candPixHits_Z            , "pixTk_candPixHits_Z"            ,  240,  -60.0,  +60.0 );
  histoTools_.BookHisto_1D( h_pixTk_candPixHits_Type         , "pixTk_candPixHits_Type"         ,   10,   -4.5,   +5.5 );
  histoTools_.BookHisto_2D( h_pixTk_candPixHits_ZVsRho       , "pixTk_candPixHits_ZVsRho"       ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_C     , "pixTk_candPixHits_ZVsRho_C"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_I     , "pixTk_candPixHits_ZVsRho_I"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_F     , "pixTk_candPixHits_ZVsRho_F"     ,  600,  -60.0,  +60.0, 100, +0.0, +20.0 );
  histoTools_.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_EtaGE2, "pixTk_candPixHits_ZVsRho_EtaGE2",  600,  -60.0,  +60.0, 100, +0.0, +20.0 );

  return;
}


//****************************************************************************
void Tracking::FillHistoBinWithRMS_(TH1D *hToFill,
				    TH1D *hToExtractRMS,
				    int binNumber)
//****************************************************************************
{
  
  hToFill->SetBinContent(binNumber, hToExtractRMS->GetRMS() );
  hToFill->SetBinError  (binNumber, hToExtractRMS->GetRMSError() );

  return;
}


  
//****************************************************************************
void Tracking::FitAndFillHistoBin_(TH1D *hToFill,
				   TH1D *hToFit,
				   int binNumber,
				   double fitRangeLow,
				   double fitRangeHigh,
				   double redChiSqMax,
				   std::vector<double> v_fitArea)
//****************************************************************************
{
  
  // Treat red-chi-sq values less than 0 as disabling the gaussian fit. Extract resolutions directly from histogram with GetRMS()
  if(redChiSqMax < 0){
    FillHistoBinWithRMS_(hToFill, hToFit, binNumber);
    return;    
  }

  // Histogramming
  if (hToFit->GetEntries() < 2) return;

  // Perform the function fit with the desired options in the requested fit range. Options:
  TF1 *fitFunction = new TF1(hToFit->GetName(), "gaus", fitRangeLow, fitRangeHigh);
  fitFunction->SetLineColor(kBlue);
  fitFunction->SetLineWidth(3);
  fitFunction->SetLineStyle(kSolid);
  fitFunction->SetParNames ("const.: ","#mu: ","#sigma: ");
  TFitResultPtr r = histoTools_.FitFunction( hToFit, fitFunction, "Q S ", "", fitRangeLow, fitRangeHigh, v_fitArea, redChiSqMax);
  
  // Get the standard deviation (sigma) of the gaussian fit and its related error
  Double_t chi2       = r->Chi2();
  Double_t dof        = r->Ndf();
  Double_t redChi2    = chi2/dof;
  Double_t sigmaValue = 0;
  Double_t sigmaError = 0;

  // Use fit-extracted value only if the fit was reasonably good. Otherwise get rms and rms-error from entire histo  
  if ( (chi2 > -1) && (redChi2 <= redChiSqMax) ){
    sigmaValue = r->Parameter(2);
    sigmaError = r->ParError(2);
    hToFill->SetBinContent(binNumber,  sigmaValue);
    hToFill->SetBinError  (binNumber,  sigmaError);
  }
  else{
    // If fit fails then get RMS of 80% of entire histogram range (ignore the 10% of the tails from each side)
    histoTools_.FindSymmetricFitRange(hToFit, 0.8, fitRangeLow, fitRangeHigh);
    hToFit->GetXaxis()->SetRange(fitRangeLow, fitRangeHigh);
    
    hToFill->SetBinContent(binNumber, hToFit->GetRMS() );
    hToFill->SetBinError  (binNumber, hToFit->GetRMSError() );
  }

   
  return;
}


//****************************************************************************
void Tracking::FitAndFillHistoBinSL_(TH1D *hToFill,
				     TH1D *hToFit,
				     int binNumber,
				     double significanceLevel)
//****************************************************************************
{
  // NOTE: The function name derives from Significance Level (SL)
  if (hToFit->GetEntries() < 2) return;

  // Perform the function fit with the desired options in the requested fit range. Options:  
  Double_t fitRangeLow  = hToFit->GetXaxis()->GetBinLowEdge(0);
  Double_t fitRangeHigh = hToFit->GetXaxis()->GetBinLowEdge( hToFit->GetNbinsX()+1 );
  TF1 *fitFunction      = new TF1(hToFit->GetName(), "gaus", fitRangeLow, fitRangeHigh);
  fitFunction->SetLineColor(kBlue);
  fitFunction->SetLineWidth(3);
  fitFunction->SetLineStyle(kSolid);
  fitFunction->SetParNames ("const.: ","#mu: ","#sigma: ");
  TFitResultPtr r = histoTools_.FitFunctionSL( hToFit, fitFunction, "Q S ", "", GetFitAreaVector(), significanceLevel, bPrintFitInfo, bSaveFitInfo);

  // Get the standard deviation (sigma) of the gaussian fit and its related error
  Double_t chi2       = r->Chi2();
  Double_t dof        = r->Ndf();
  Double_t redChi2    = chi2/dof;
  Double_t sigmaValue = 0;
  Double_t sigmaError = 0;

  // Use fit-extracted value only if the fit was reasonably good. Otherwise get rms and rms-error from entire histo  
  if (chi2 > -1){
    sigmaValue = r->Parameter(2);
    sigmaError = r->ParError(2);
    hToFill->SetBinContent(binNumber,  sigmaValue);
    hToFill->SetBinError  (binNumber,  sigmaError);
  }
  else{
    // If fit fails then get RMS of 80% of entire histogram range (ignore the 10% of the tails from each side)
    histoTools_.FindSymmetricFitRange(hToFit, 0.8, fitRangeLow, fitRangeHigh);
    hToFit->GetXaxis()->SetRange(fitRangeLow, fitRangeHigh);
    
    hToFill->SetBinContent(binNumber, hToFit->GetRMS() );
    hToFill->SetBinError  (binNumber, hToFit->GetRMSError() );
  }

   
  return;
}


//****************************************************************************
std::vector<double> Tracking::GetFitAreaVector(const double minTotalArea)
//****************************************************************************
{

  // Fit area from 0.5 -> minTotalArea/2
  std::vector<double> v_fitArea;  
  for (int i = 0; i <= 20; i++){
    double area = 0.5 - i*0.02;
    if (2*area < minTotalArea) break;
    v_fitArea.push_back(area);
  }
  // auxTools_.PrintVector(v_fitArea);
  
  return v_fitArea;
}


//****************************************************************************
void Tracking::FinaliseHistos_(void)
//****************************************************************************
{
  
  // Loop over all pt bins
  for (int i=0; i < PT_BINS; i++) {

    // Cannot fit "gaus" to phi due to asymmetric distribution (lorentz-drift)
    FitAndFillHistoBinSL_( h2_resVsPt_pt  , h_resVsPt_pt[i]  , i+1, +0.005);
    FitAndFillHistoBinSL_( h2_resVsPt_pt_C, h_resVsPt_pt_C[i], i+1, +0.005);
    FitAndFillHistoBinSL_( h2_resVsPt_pt_I, h_resVsPt_pt_I[i], i+1, +0.005);
    FitAndFillHistoBinSL_( h2_resVsPt_pt_F, h_resVsPt_pt_F[i], i+1, +0.005);

    FitAndFillHistoBinSL_( h2_resVsPt_ptRel  , h_resVsPt_ptRel[i]  , i+1, +0.005);
    FitAndFillHistoBinSL_( h2_resVsPt_ptRel_C, h_resVsPt_ptRel_C[i], i+1, +0.005);
    FitAndFillHistoBinSL_( h2_resVsPt_ptRel_I, h_resVsPt_ptRel_I[i], i+1, +0.005);
    FitAndFillHistoBinSL_( h2_resVsPt_ptRel_F, h_resVsPt_ptRel_F[i], i+1, +0.005);
    
    // Fit "gaus" to eta to get resolution
    FitAndFillHistoBinSL_( h2_resVsPt_eta  , h_resVsPt_eta[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_eta_C, h_resVsPt_eta_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_eta_I, h_resVsPt_eta_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_eta_F, h_resVsPt_eta_F[i], i+1, +0.005 );

    // Cannot fit "gaus" to phi due to asymmetric distribution (lorentz-drift)
    FitAndFillHistoBinSL_( h2_resVsPt_phi  , h_resVsPt_phi[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_phi_C, h_resVsPt_phi_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_phi_I, h_resVsPt_phi_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_phi_F, h_resVsPt_phi_F[i], i+1, +0.005 );

    // Fit "gaus" to z0 to get resolution
    FitAndFillHistoBinSL_( h2_resVsPt_z0  , h_resVsPt_z0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_z0_C, h_resVsPt_z0_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_z0_I, h_resVsPt_z0_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_z0_F, h_resVsPt_z0_F[i], i+1, +0.005 );

    // Fit "gaus" to d0 to get resolution
    FitAndFillHistoBinSL_( h2_resVsPt_d0  , h_resVsPt_d0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_d0_C, h_resVsPt_d0_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_d0_I, h_resVsPt_d0_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL_( h2_resVsPt_d0_F, h_resVsPt_d0_F[i], i+1, +0.005 );

  }

  // Loop over all eta bins
  for (int i=0; i < ETA_BINS; i++) {

    // double eta = i*ETA_BINWIDTH;
    
    FitAndFillHistoBinSL_(h2_resVsEta_pt  , h_resVsEta_pt[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_pt_L, h_resVsEta_pt_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_pt_M, h_resVsEta_pt_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_pt_H, h_resVsEta_pt_H[i], i+1, +0.005 );
  
    FitAndFillHistoBinSL_(h2_resVsEta_ptRel  , h_resVsEta_ptRel[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_ptRel_L, h_resVsEta_ptRel_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_ptRel_M, h_resVsEta_ptRel_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_ptRel_H, h_resVsEta_ptRel_H[i], i+1, +0.005 );

    FitAndFillHistoBinSL_(h2_resVsEta_eta  , h_resVsEta_eta[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_eta_L, h_resVsEta_eta_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_eta_M, h_resVsEta_eta_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_eta_H, h_resVsEta_eta_H[i], i+1, +0.005 );
    
    FitAndFillHistoBinSL_(h2_resVsEta_phi  , h_resVsEta_phi[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_phi_L, h_resVsEta_phi_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_phi_M, h_resVsEta_phi_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_phi_H, h_resVsEta_phi_H[i], i+1, +0.005 );

    FitAndFillHistoBinSL_(h2_resVsEta_z0  , h_resVsEta_z0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_z0_L, h_resVsEta_z0_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_z0_M, h_resVsEta_z0_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_z0_H, h_resVsEta_z0_H[i], i+1, +0.005 );

    FitAndFillHistoBinSL_(h2_resVsEta_d0  , h_resVsEta_d0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_d0_L, h_resVsEta_d0_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_d0_M, h_resVsEta_d0_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL_(h2_resVsEta_d0_H, h_resVsEta_d0_H[i], i+1, +0.005 );
    
  }
  
  return;
}

    

//****************************************************************************
void Tracking::MakeEfficiencyHisto(TH1D* h_eff, 
				   TH1D* h_match_tp, 
				   TH1D* h_tp)
//****************************************************************************
{

  h_eff->Divide(h_match_tp, h_tp, 1.0, 1.0, "B");
  
  return;
}
 
 

//****************************************************************************
void Tracking::PrintSettings_(void)
//****************************************************************************
{

  // Inform user of settings  
  Table settings("Symbol | Definition | Units | Description", "Text");
  settings.AddRowColumn(0, "_eta_C");
  settings.AddRowColumn(0, "abs(eta) < " + auxTools_.ToString(_eta_C) );
  settings.AddRowColumn(0, "-");
  settings.AddRowColumn(0, "Central Eta Region");
  settings.AddRowColumn(1, "_eta_I");
  settings.AddRowColumn(1, auxTools_.ToString(_eta_C) + " < abs(eta) < " + auxTools_.ToString(_eta_F));
  settings.AddRowColumn(1, "-");
  settings.AddRowColumn(1, "Intermediate Eta Region");
  settings.AddRowColumn(2, "_eta_F");
  settings.AddRowColumn(2, "abs(eta) > " + auxTools_.ToString(_eta_F) );
  settings.AddRowColumn(2, "-");
  settings.AddRowColumn(2, "Forward Eta Region");
  settings.AddRowColumn(3, "_pt_L");
  settings.AddRowColumn(3, "pT < " + auxTools_.ToString(_pt_L));
  settings.AddRowColumn(3, "GeV/c");
  settings.AddRowColumn(3, "Low pT Range");
  settings.AddRowColumn(4, "_pt_M");
  settings.AddRowColumn(4, auxTools_.ToString(_pt_L) + " GeV/c < pT < " + auxTools_.ToString(_pt_H) );
  settings.AddRowColumn(4, "GeV/c");
  settings.AddRowColumn(4, "Middle pT Range");
  settings.AddRowColumn(5, "_pt_H");
  settings.AddRowColumn(5, "pT > " + auxTools_.ToString(_pt_H) );
  settings.AddRowColumn(5, "GeV/c");
  settings.AddRowColumn(5, "High pT Range");
  settings.AddRowColumn(6, "_pt_acceptance");
  settings.AddRowColumn(6, auxTools_.ToString(_pt_acceptance) );
  settings.AddRowColumn(6, "GeV/c");
  settings.AddRowColumn(6, "TPs with smaller pt are outside acceptance" );
  settings.AddRowColumn(7, "_eta_acceptance");
  settings.AddRowColumn(7, auxTools_.ToString(_eta_acceptance) );
  settings.AddRowColumn(7, "" );
  settings.AddRowColumn(7, "TPs with larger abs(eta) are outside acceptance" );
  settings.AddRowColumn(8, "_z0_acceptance");
  settings.AddRowColumn(8, auxTools_.ToString(_z0_acceptance) );
  settings.AddRowColumn(8, "cm" );
  settings.AddRowColumn(8, "TPs with larger abs(z0) are outside acceptance" );
  settings.AddRowColumn(9, "_pt_tkMin");
  settings.AddRowColumn(9, auxTools_.ToString(_pt_tkMin) );
  settings.AddRowColumn(9, "GeV/c" );
  settings.AddRowColumn(9, "Minimum allowable L1 track pt" );
  settings.AddRowColumn(10, "_nStubs_tkMin");
  settings.AddRowColumn(10, auxTools_.ToString(_nStubs_tkMin) );
  settings.AddRowColumn(10, "-" );
  settings.AddRowColumn(10, "Minimum number of stubs for a L1 track" );
  settings.AddRowColumn(11, "_chiSq_overflow");
  settings.AddRowColumn(11, auxTools_.ToString(_chiSq_overflow) );
  settings.AddRowColumn(11, "-" );
  settings.AddRowColumn(11, "Value of overflow bin for chi-squared of L1 tracks" );
  settings.AddRowColumn(12, "_redChiSq_overflow");
  settings.AddRowColumn(12, auxTools_.ToString(_redChiSq_overflow) );
  settings.AddRowColumn(12, "-" );
  settings.AddRowColumn(12, "Value of overflow bin for reduced-chi-squared of L1 tracks" );
  settings.AddRowColumn(13, "bDoPixelTracks");
  settings.AddRowColumn(13, auxTools_.ToString(bDoPixelTracks) );
  settings.AddRowColumn(13, "-" );
  settings.AddRowColumn(13, "Enable (set to true) to make plots for TTPixelTracks instead of TTTracks");
  settings.AddRowColumn(14, "bVerbose");
  settings.AddRowColumn(14, auxTools_.ToString(bVerbose) );
  settings.AddRowColumn(14, "-" );
  settings.AddRowColumn(14, "Enable (set to true) to print info");
  settings.AddRowColumn(15, "tp_dxy_Max");
  settings.AddRowColumn(15, auxTools_.ToString(tp_dxy_Max) );
  settings.AddRowColumn(15, "cm");
  settings.AddRowColumn(15, "Max dxy of TPs considered (Look at those from near the IP)");

  settings.Print();
  
  return;
}


//****************************************************************************
void Tracking::PrintPixelTrackCuts_(void)
//****************************************************************************
{

  // Inform user of settings  
  Table pixCuts("Variable | Units | Cut Direction | Cut Value | Explanation", "Text");
  pixCuts.AddRowColumn(0, "pixTk_NPixHits_Min");
  pixCuts.AddRowColumn(0, "" );
  pixCuts.AddRowColumn(0, " >= " );
  pixCuts.AddRowColumn(0, auxTools_.ToString(pixTk_NPixHits_Min) );
  pixCuts.AddRowColumn(0, "Minimum number of Pixel Fit-Hits for TTPixelTracks" );

  pixCuts.AddRowColumn(1, "pixTk_NPixHits_Max");
  pixCuts.AddRowColumn(1, "" );
  pixCuts.AddRowColumn(1, " < " );
  pixCuts.AddRowColumn(1, auxTools_.ToString(pixTk_NPixHits_Max) );
  pixCuts.AddRowColumn(1, "Maximum number of Pixel Fit-Hits per TTPixelTracks" );
  
  pixCuts.AddRowColumn(2, "pixTk_MustHit_Type");
  pixCuts.AddRowColumn(2, "" );
  pixCuts.AddRowColumn(2, "==" );
  pixCuts.AddRowColumn(2, auxTools_.ConvertIntVectorToString(pixTk_MustHit_Type) );
  pixCuts.AddRowColumn(2, "Pixel Hit Types that TTPixelTracks must at least satisfy ONE of them" );

  // pixCuts.AddRowColumn(2, "pixTk-Barrel-Type");
  // pixCuts.AddRowColumn(2, "" );
  // pixCuts.AddRowColumn(2, "==" );
  // pixCuts.AddRowColumn(2, auxTools_.ConvertIntVectorToString(pixTk_Barrel_Type) );
  // pixCuts.AddRowColumn(2, "Pixel Hit Types that TTPixelTracks must have (eta <= pixTk_Type_EndcapEtaBoundary)" );
  
  // pixCuts.AddRowColumn(3, "pixTk_Endcap_Type");
  // pixCuts.AddRowColumn(3, "" );
  // pixCuts.AddRowColumn(3, "==" );
  // pixCuts.AddRowColumn(3, auxTools_.ConvertIntVectorToString(pixTk_Endcap_Type) );
  // pixCuts.AddRowColumn(3, "Pixel Hit Types that TTPixelTracks must have (eta > pixTk_Type_EndcapEtaBoundary" );

  // pixCuts.AddRowColumn(4, "pixTk_Type_EndcapEtaBoundary");
  // pixCuts.AddRowColumn(4, "" );
  // pixCuts.AddRowColumn(4, "==" );
  // pixCuts.AddRowColumn(4, auxTools_.ToString(pixTk_Type_EndcapEtaBoundary) );
  // pixCuts.AddRowColumn(4, "Where the barrel ends and endcap  starts (when demanding pixel fit-hit types)" );

  //
  pixCuts.Print();
  
  return;
}


//****************************************************************************
void Tracking::PrintEfficiencies_(void)
//****************************************************************************
{

  // Calculate efficiencies
  double eff_eta_C   = 0.0;
  double err_eta_C   = 0.0;
  double eff_eta_I   = 0.0;
  double err_eta_I   = 0.0;
  double eff_eta_F   = 0.0;
  double err_eta_F   = 0.0;
  double eff_overall = 0.0;
  double err_overall = 0.0;
  auxTools_.Efficiency(n_match_eta_C, n_all_eta_C, "binomial", eff_eta_C, err_eta_C);
  auxTools_.Efficiency(n_match_eta_I, n_all_eta_I, "binomial", eff_eta_I, err_eta_I);
  auxTools_.Efficiency(n_match_eta_F, n_all_eta_F, "binomial", eff_eta_F, err_eta_F);
  auxTools_.Efficiency(n_match_eta_C + n_match_eta_I + n_match_eta_F, n_all_eta_C +  n_all_eta_I +  n_all_eta_F, "binomial", eff_overall, err_overall);

  double eff_eta_C_pt_L   = 0.0;
  double err_eta_C_pt_L   = 0.0;
  double eff_eta_I_pt_L   = 0.0;
  double err_eta_I_pt_L   = 0.0;
  double eff_eta_F_pt_L   = 0.0;
  double err_eta_F_pt_L   = 0.0;
  double eff_overall_pt_L = 0.0;
  double err_overall_pt_L = 0.0;
  auxTools_.Efficiency(n_match_eta_C_pt_L, n_all_eta_C_pt_L, "binomial", eff_eta_C_pt_L, err_eta_C_pt_L);
  auxTools_.Efficiency(n_match_eta_I_pt_L, n_all_eta_I_pt_L, "binomial", eff_eta_I_pt_L, err_eta_I_pt_L);
  auxTools_.Efficiency(n_match_eta_F_pt_L, n_all_eta_F_pt_L, "binomial", eff_eta_F_pt_L, err_eta_F_pt_L);
  auxTools_.Efficiency(n_match_eta_C_pt_L + n_match_eta_I_pt_L + n_match_eta_F_pt_L, n_all_eta_C_pt_L +  n_all_eta_I_pt_L +  n_all_eta_F_pt_L, "binomial", eff_overall_pt_L, err_overall_pt_L);

  double eff_eta_C_pt_M   = 0.0;
  double err_eta_C_pt_M   = 0.0;
  double eff_eta_I_pt_M   = 0.0;
  double err_eta_I_pt_M   = 0.0;
  double eff_eta_F_pt_M   = 0.0;
  double err_eta_F_pt_M   = 0.0;
  double eff_overall_pt_M = 0.0;
  double err_overall_pt_M = 0.0;
  auxTools_.Efficiency(n_match_eta_C_pt_M, n_all_eta_C_pt_M, "binomial", eff_eta_C_pt_M, err_eta_C_pt_M);
  auxTools_.Efficiency(n_match_eta_I_pt_M, n_all_eta_I_pt_M, "binomial", eff_eta_I_pt_M, err_eta_I_pt_M);
  auxTools_.Efficiency(n_match_eta_F_pt_M, n_all_eta_F_pt_M, "binomial", eff_eta_F_pt_M, err_eta_F_pt_M);
  auxTools_.Efficiency(n_match_eta_C_pt_M + n_match_eta_I_pt_M + n_match_eta_F_pt_M, n_all_eta_C_pt_M +  n_all_eta_I_pt_M +  n_all_eta_F_pt_M, "binomial", eff_overall_pt_M, err_overall_pt_M);

  double eff_eta_C_pt_H   = 0.0;
  double err_eta_C_pt_H   = 0.0;
  double eff_eta_I_pt_H   = 0.0;
  double err_eta_I_pt_H   = 0.0;
  double eff_eta_F_pt_H   = 0.0;
  double err_eta_F_pt_H   = 0.0;
  double eff_overall_pt_H = 0.0;
  double err_overall_pt_H = 0.0;
  auxTools_.Efficiency(n_match_eta_C_pt_H, n_all_eta_C_pt_H, "binomial", eff_eta_C_pt_H, err_eta_C_pt_H);
  auxTools_.Efficiency(n_match_eta_I_pt_H, n_all_eta_I_pt_H, "binomial", eff_eta_I_pt_H, err_eta_I_pt_H);
  auxTools_.Efficiency(n_match_eta_F_pt_H, n_all_eta_F_pt_H, "binomial", eff_eta_F_pt_H, err_eta_F_pt_H);
  auxTools_.Efficiency(n_match_eta_C_pt_H + n_match_eta_I_pt_H + n_match_eta_F_pt_H, n_all_eta_C_pt_H +  n_all_eta_I_pt_H +  n_all_eta_F_pt_H, "binomial", eff_overall_pt_H, err_overall_pt_H);

  // Create efficiencies table
  Table efficiencies("\\pT Range | " + s_eta_all + " | " + s_eta_C + " | " + s_eta_I + " | " + s_eta_F, "LaTeX", "r R R R R");
  efficiencies.AddRowColumn(0, "$" + s_pt_tkMin + "$");
  efficiencies.AddRowColumn(0, auxTools_.ToString(eff_overall, 2) + " \\pm " + auxTools_.ToString(err_overall, 2) );
  efficiencies.AddRowColumn(0, auxTools_.ToString(eff_eta_C  , 2) + " \\pm " + auxTools_.ToString(err_eta_C  , 2) );
  efficiencies.AddRowColumn(0, auxTools_.ToString(eff_eta_I  , 2) + " \\pm " + auxTools_.ToString(err_eta_I  , 2) );
  efficiencies.AddRowColumn(0, auxTools_.ToString(eff_eta_F  , 2) + " \\pm " + auxTools_.ToString(err_eta_F  , 2) );

  efficiencies.AddRowColumn(1, "$" + s_pt_L + "$");
  efficiencies.AddRowColumn(1, auxTools_.ToString(eff_overall_pt_L, 2) + " \\pm " + auxTools_.ToString(err_overall_pt_L, 2) );
  efficiencies.AddRowColumn(1, auxTools_.ToString(eff_eta_C_pt_L  , 2) + " \\pm " + auxTools_.ToString(err_eta_C_pt_L  , 2) );
  efficiencies.AddRowColumn(1, auxTools_.ToString(eff_eta_I_pt_L  , 2) + " \\pm " + auxTools_.ToString(err_eta_I_pt_L  , 2) );
  efficiencies.AddRowColumn(1, auxTools_.ToString(eff_eta_F_pt_L  , 2) + " \\pm " + auxTools_.ToString(err_eta_F_pt_L  , 2) );

  efficiencies.AddRowColumn(2, "$" + s_pt_M + "$");
  efficiencies.AddRowColumn(2, auxTools_.ToString(eff_overall_pt_M, 2) + " \\pm " + auxTools_.ToString(err_overall_pt_M, 2) );
  efficiencies.AddRowColumn(2, auxTools_.ToString(eff_eta_C_pt_M  , 2) + " \\pm " + auxTools_.ToString(err_eta_C_pt_M  , 2) );
  efficiencies.AddRowColumn(2, auxTools_.ToString(eff_eta_I_pt_M  , 2) + " \\pm " + auxTools_.ToString(err_eta_I_pt_M  , 2) );
  efficiencies.AddRowColumn(2, auxTools_.ToString(eff_eta_F_pt_M  , 2) + " \\pm " + auxTools_.ToString(err_eta_F_pt_M  , 2) );

  efficiencies.AddRowColumn(3, "$" + s_pt_H + "$");
  efficiencies.AddRowColumn(3, auxTools_.ToString(eff_overall_pt_H, 2) + " \\pm " + auxTools_.ToString(err_overall_pt_H, 2) );
  efficiencies.AddRowColumn(3, auxTools_.ToString(eff_eta_C_pt_H  , 2) + " \\pm " + auxTools_.ToString(err_eta_C_pt_H  , 2) );
  efficiencies.AddRowColumn(3, auxTools_.ToString(eff_eta_I_pt_H  , 2) + " \\pm " + auxTools_.ToString(err_eta_I_pt_H  , 2) );
  efficiencies.AddRowColumn(3, auxTools_.ToString(eff_eta_F_pt_H  , 2) + " \\pm " + auxTools_.ToString(err_eta_F_pt_H  , 2) );
  
  if(bPrintEfficiencies) efficiencies.Print();

return;
}



//****************************************************************************
void Tracking::PrintResolutions_(void)
 //****************************************************************************
 {
   
   // Calculate resolutions for table
   double ptResolution_L    = GetHistoMeanInRange(h2_resVsEta_pt_L , 0.0, 2.5);
   double ptResolution_L_C  = GetHistoMeanInRange(h2_resVsEta_pt_L , 0.0, 0.8);
   double ptResolution_L_I  = GetHistoMeanInRange(h2_resVsEta_pt_L , 0.8, 1.6);
   double ptResolution_L_F  = GetHistoMeanInRange(h2_resVsEta_pt_L , 1.6, 2.5);

   double etaResolution_L   = GetHistoMeanInRange(h2_resVsEta_eta_L, 0.0, 2.5);
   double etaResolution_L_C = GetHistoMeanInRange(h2_resVsEta_eta_L, 0.0, 0.8);
   double etaResolution_L_I = GetHistoMeanInRange(h2_resVsEta_eta_L, 0.8, 1.6);
   double etaResolution_L_F = GetHistoMeanInRange(h2_resVsEta_eta_L, 1.6, 2.5);

   double phiResolution_L   = GetHistoMeanInRange(h2_resVsEta_phi_L, 0.0, 2.5);
   double phiResolution_L_C = GetHistoMeanInRange(h2_resVsEta_phi_L, 0.0, 0.8);
   double phiResolution_L_I = GetHistoMeanInRange(h2_resVsEta_phi_L, 0.8, 1.6);
   double phiResolution_L_F = GetHistoMeanInRange(h2_resVsEta_phi_L, 1.6, 2.5);

   double z0Resolution_L    = GetHistoMeanInRange(h2_resVsEta_z0_L , 0.0, 2.5);
   double z0Resolution_L_C  = GetHistoMeanInRange(h2_resVsEta_z0_L , 0.0, 0.8);
   double z0Resolution_L_I  = GetHistoMeanInRange(h2_resVsEta_z0_L , 0.8, 1.6);
   double z0Resolution_L_F  = GetHistoMeanInRange(h2_resVsEta_z0_L , 1.6, 2.5);

   double d0Resolution_L    = GetHistoMeanInRange(h2_resVsEta_d0_L , 0.0, 2.5);
   double d0Resolution_L_C  = GetHistoMeanInRange(h2_resVsEta_d0_L , 0.0, 0.8);
   double d0Resolution_L_I  = GetHistoMeanInRange(h2_resVsEta_d0_L , 0.8, 1.6);
   double d0Resolution_L_F  = GetHistoMeanInRange(h2_resVsEta_d0_L , 1.6, 2.5);

   double ptResolution_M    = GetHistoMeanInRange(h2_resVsEta_pt_M, 0.0, 2.5);
   double ptResolution_M_C  = GetHistoMeanInRange(h2_resVsEta_pt_M, 0.0, 0.8);
   double ptResolution_M_I  = GetHistoMeanInRange(h2_resVsEta_pt_M, 0.8, 1.6);
   double ptResolution_M_F  = GetHistoMeanInRange(h2_resVsEta_pt_M, 1.6, 2.5);

   double etaResolution_M   = GetHistoMeanInRange(h2_resVsEta_eta_M, 0.0, 2.5);
   double etaResolution_M_C = GetHistoMeanInRange(h2_resVsEta_eta_M, 0.0, 0.8);
   double etaResolution_M_I = GetHistoMeanInRange(h2_resVsEta_eta_M, 0.8, 1.6);
   double etaResolution_M_F = GetHistoMeanInRange(h2_resVsEta_eta_M, 1.6, 2.5);

   double phiResolution_M   = GetHistoMeanInRange(h2_resVsEta_phi_M, 0.0, 2.5);
   double phiResolution_M_C = GetHistoMeanInRange(h2_resVsEta_phi_M, 0.0, 0.8);
   double phiResolution_M_I = GetHistoMeanInRange(h2_resVsEta_phi_M, 0.8, 1.6);
   double phiResolution_M_F = GetHistoMeanInRange(h2_resVsEta_phi_M, 1.6, 2.5);

   double z0Resolution_M    = GetHistoMeanInRange(h2_resVsEta_z0_M, 0.0, 2.5);
   double z0Resolution_M_C  = GetHistoMeanInRange(h2_resVsEta_z0_M, 0.0, 0.8);
   double z0Resolution_M_I  = GetHistoMeanInRange(h2_resVsEta_z0_M, 0.8, 1.6);
   double z0Resolution_M_F  = GetHistoMeanInRange(h2_resVsEta_z0_M, 1.6, 2.5);

   double d0Resolution_M    = GetHistoMeanInRange(h2_resVsEta_d0_M, 0.0, 2.5);
   double d0Resolution_M_C  = GetHistoMeanInRange(h2_resVsEta_d0_M, 0.0, 0.8);
   double d0Resolution_M_I  = GetHistoMeanInRange(h2_resVsEta_d0_M, 0.8, 1.6);
   double d0Resolution_M_F  = GetHistoMeanInRange(h2_resVsEta_d0_M, 1.6, 2.5);

   double ptResolution_H    = GetHistoMeanInRange(h2_resVsEta_pt_H, 0.0, 2.5);
   double ptResolution_H_C  = GetHistoMeanInRange(h2_resVsEta_pt_H, 0.0, 0.8);
   double ptResolution_H_I  = GetHistoMeanInRange(h2_resVsEta_pt_H, 0.8, 1.6);
   double ptResolution_H_F  = GetHistoMeanInRange(h2_resVsEta_pt_H, 1.6, 2.5);

   double etaResolution_H   = GetHistoMeanInRange(h2_resVsEta_eta_H, 0.0, 2.5);
   double etaResolution_H_C = GetHistoMeanInRange(h2_resVsEta_eta_H, 0.0, 0.8);
   double etaResolution_H_I = GetHistoMeanInRange(h2_resVsEta_eta_H, 0.8, 1.6);
   double etaResolution_H_F = GetHistoMeanInRange(h2_resVsEta_eta_H, 1.6, 2.5);

   double phiResolution_H   = GetHistoMeanInRange(h2_resVsEta_phi_H, 0.0, 2.5);
   double phiResolution_H_C = GetHistoMeanInRange(h2_resVsEta_phi_H, 0.0, 0.8);
   double phiResolution_H_I = GetHistoMeanInRange(h2_resVsEta_phi_H, 0.8, 1.6);
   double phiResolution_H_F = GetHistoMeanInRange(h2_resVsEta_phi_H, 1.6, 2.5);

   double z0Resolution_H    = GetHistoMeanInRange(h2_resVsEta_z0_H, 0.0, 2.5);
   double z0Resolution_H_C  = GetHistoMeanInRange(h2_resVsEta_z0_H, 0.0, 0.8);
   double z0Resolution_H_I  = GetHistoMeanInRange(h2_resVsEta_z0_H, 0.8, 1.6);
   double z0Resolution_H_F  = GetHistoMeanInRange(h2_resVsEta_z0_H, 1.6, 2.5);

   double d0Resolution_H    = GetHistoMeanInRange(h2_resVsEta_d0_H, 0.0, 0.8);
   double d0Resolution_H_C  = GetHistoMeanInRange(h2_resVsEta_d0_H, 0.0, 0.8);
   double d0Resolution_H_I  = GetHistoMeanInRange(h2_resVsEta_d0_H, 0.8, 1.6);
   double d0Resolution_H_F  = GetHistoMeanInRange(h2_resVsEta_d0_H, 1.6, 2.5);

   
   // Create resolutions table
   Table resolutions("Resolution | Units | \\abs{\\eta} < 0.8 | 0.8 < \\abs{\\eta} < 1.6 | \\abs{\\eta} > 1.6 | \\abs{\\eta} \\leq 2.5", "LaTeX", "c c R R R R");
   resolutions.AddRowColumn(0, "\\multicolumn{6}{c}{$\\pT < " + auxTools_.ToString(_pt_L) + " \\GeVc{-1}$}");
   resolutions.AddRowColumn(1, "$\\pT$");
   resolutions.AddRowColumn(1, "$\\MeVc{-1}$");
   resolutions.AddRowColumn(1, auxTools_.ToString(1000 * ptResolution_L_C, 2) );
   resolutions.AddRowColumn(1, auxTools_.ToString(1000 * ptResolution_L_I, 2) );
   resolutions.AddRowColumn(1, auxTools_.ToString(1000 * ptResolution_L_F, 2) );
   resolutions.AddRowColumn(1, auxTools_.ToString(1000 * ptResolution_L  , 2) );
   resolutions.AddRowColumn(2, "$\\eta$");
   resolutions.AddRowColumn(2, "$10^{-3}$");
   resolutions.AddRowColumn(2, auxTools_.ToString(1000 * etaResolution_L_C, 2) );
   resolutions.AddRowColumn(2, auxTools_.ToString(1000 * etaResolution_L_I, 2) );
   resolutions.AddRowColumn(2, auxTools_.ToString(1000 * etaResolution_L_F, 2) );
   resolutions.AddRowColumn(2, auxTools_.ToString(1000 * etaResolution_L  , 2) );
   resolutions.AddRowColumn(3, "$\\phi$");
   resolutions.AddRowColumn(3, "$\\sMilliRads$");
   resolutions.AddRowColumn(3, auxTools_.ToString(1000 * phiResolution_L_C, 2) );
   resolutions.AddRowColumn(3, auxTools_.ToString(1000 * phiResolution_L_I, 2) );
   resolutions.AddRowColumn(3, auxTools_.ToString(1000 * phiResolution_L_F, 2) );
   resolutions.AddRowColumn(3, auxTools_.ToString(1000 * phiResolution_L  , 2) );
   resolutions.AddRowColumn(4, "$\\zPOCA{}$");
   resolutions.AddRowColumn(4, "$\\sMilliMeter$");
   resolutions.AddRowColumn(4, auxTools_.ToString(10 * z0Resolution_L_C, 2) );
   resolutions.AddRowColumn(4, auxTools_.ToString(10 * z0Resolution_L_I, 2) );
   resolutions.AddRowColumn(4, auxTools_.ToString(10 * z0Resolution_L_F, 2) );
   resolutions.AddRowColumn(4, auxTools_.ToString(10 * z0Resolution_L  , 2) );
   resolutions.AddRowColumn(5, "$\\dZero{}$");
   resolutions.AddRowColumn(5, "$\\sMicroMeter$");
   resolutions.AddRowColumn(5, auxTools_.ToString(10000 * d0Resolution_L_C, 3) );
   resolutions.AddRowColumn(5, auxTools_.ToString(10000 * d0Resolution_L_I, 3) );
   resolutions.AddRowColumn(5, auxTools_.ToString(10000 * d0Resolution_L_F, 3) );
   resolutions.AddRowColumn(5, auxTools_.ToString(10000 * d0Resolution_L  , 3) );

   resolutions.AddRowColumn(6, "\\multicolumn{6}{c}{$" + auxTools_.ToString(_pt_L) + " \\GeVc{-1}$} < \\pT < " + auxTools_.ToString(_pt_H) + " \\GeVc{-1}$}");
   resolutions.AddRowColumn(7, "$\\pT$");
   resolutions.AddRowColumn(7, "$\\MeVc{-1}$");
   resolutions.AddRowColumn(7, auxTools_.ToString(1000 * ptResolution_M_C, 2) );
   resolutions.AddRowColumn(7, auxTools_.ToString(1000 * ptResolution_M_I, 2) );
   resolutions.AddRowColumn(7, auxTools_.ToString(1000 * ptResolution_M_F, 2) );
   resolutions.AddRowColumn(7, auxTools_.ToString(1000 * ptResolution_M  , 2) );
   resolutions.AddRowColumn(8, "$\\eta$");
   resolutions.AddRowColumn(8, "$10^{-3}$");   
   resolutions.AddRowColumn(8, auxTools_.ToString(1000 * etaResolution_M_C, 2) );
   resolutions.AddRowColumn(8, auxTools_.ToString(1000 * etaResolution_M_I, 2) );
   resolutions.AddRowColumn(8, auxTools_.ToString(1000 * etaResolution_M_F, 2) );
   resolutions.AddRowColumn(8, auxTools_.ToString(1000 * etaResolution_M  , 2) );
   resolutions.AddRowColumn(9, "$\\phi$");
   resolutions.AddRowColumn(9, "$\\sMilliRads$");  
   resolutions.AddRowColumn(9, auxTools_.ToString(1000 * phiResolution_M_C, 2) );
   resolutions.AddRowColumn(9, auxTools_.ToString(1000 * phiResolution_M_I, 2) );
   resolutions.AddRowColumn(9, auxTools_.ToString(1000 * phiResolution_M_F, 2) );
   resolutions.AddRowColumn(9, auxTools_.ToString(1000 * phiResolution_M  , 2) );
   resolutions.AddRowColumn(10, "$\\zPOCA{}$");
   resolutions.AddRowColumn(10, "$\\sMilliMeter$");
   resolutions.AddRowColumn(10, auxTools_.ToString(10 * z0Resolution_M_C, 2) );
   resolutions.AddRowColumn(10, auxTools_.ToString(10 * z0Resolution_M_I, 2) );
   resolutions.AddRowColumn(10, auxTools_.ToString(10 * z0Resolution_M_F, 2) );
   resolutions.AddRowColumn(10, auxTools_.ToString(10 * z0Resolution_M  , 2) );
   resolutions.AddRowColumn(11, "$\\dZero{}$");
   resolutions.AddRowColumn(11, "$\\sMicroMeter$");   
   resolutions.AddRowColumn(11, auxTools_.ToString(10000 * d0Resolution_M_C, 3) );
   resolutions.AddRowColumn(11, auxTools_.ToString(10000 * d0Resolution_M_I, 3) );
   resolutions.AddRowColumn(11, auxTools_.ToString(10000 * d0Resolution_M_F, 3) );
   resolutions.AddRowColumn(11, auxTools_.ToString(10000 * d0Resolution_M  , 3) );

   resolutions.AddRowColumn(12, "\\multicolumn{6}{c}{$\\pT \\geq " + auxTools_.ToString(_pt_H) + " \\GeVc{-1}$}");
   resolutions.AddRowColumn(13, "$\\pT$");
   resolutions.AddRowColumn(13, "$\\MeVc{-1}$");  
   resolutions.AddRowColumn(13, auxTools_.ToString(1000 * ptResolution_H_C, 2) );
   resolutions.AddRowColumn(13, auxTools_.ToString(1000 * ptResolution_H_I, 2) );
   resolutions.AddRowColumn(13, auxTools_.ToString(1000 * ptResolution_H_F, 2) );
   resolutions.AddRowColumn(13, auxTools_.ToString(1000 * ptResolution_H  , 2) );
   resolutions.AddRowColumn(14, "$\\eta$");
   resolutions.AddRowColumn(14, "$10^{-3}$");  
   resolutions.AddRowColumn(14, auxTools_.ToString(1000 * etaResolution_H_C, 2) );
   resolutions.AddRowColumn(14, auxTools_.ToString(1000 * etaResolution_H_I, 2) );
   resolutions.AddRowColumn(14, auxTools_.ToString(1000 * etaResolution_H_F, 2) );
   resolutions.AddRowColumn(14, auxTools_.ToString(1000 * etaResolution_H  , 2) );
   resolutions.AddRowColumn(15, "$\\phi$");
   resolutions.AddRowColumn(15, "$\\sMilliRads$");
   resolutions.AddRowColumn(15, auxTools_.ToString(1000 * phiResolution_H_C, 2) );
   resolutions.AddRowColumn(15, auxTools_.ToString(1000 * phiResolution_H_I, 2) );
   resolutions.AddRowColumn(15, auxTools_.ToString(1000 * phiResolution_H_F, 2) );
   resolutions.AddRowColumn(15, auxTools_.ToString(1000 * phiResolution_H  , 2) );
   resolutions.AddRowColumn(16, "$\\zPOCA{}$");
   resolutions.AddRowColumn(16, "$\\sMilliMeter$");
   resolutions.AddRowColumn(16, auxTools_.ToString(10 * z0Resolution_H_C, 2) );
   resolutions.AddRowColumn(16, auxTools_.ToString(10 * z0Resolution_H_I, 2) );
   resolutions.AddRowColumn(16, auxTools_.ToString(10 * z0Resolution_H_F, 2) );
   resolutions.AddRowColumn(16, auxTools_.ToString(10 * z0Resolution_H  , 2) );
   resolutions.AddRowColumn(17, "$\\dZero{}$");
   resolutions.AddRowColumn(17, "$\\sMicroMeter$");
   resolutions.AddRowColumn(17, auxTools_.ToString(10000 * d0Resolution_H_C, 3) );
   resolutions.AddRowColumn(17, auxTools_.ToString(10000 * d0Resolution_H_I, 3) );
   resolutions.AddRowColumn(17, auxTools_.ToString(10000 * d0Resolution_H_F, 3) );
   resolutions.AddRowColumn(17, auxTools_.ToString(10000 * d0Resolution_H  , 3) );

   if(bPrintResolutions) resolutions.Print();

   return;
 }



//****************************************************************************
 double Tracking::GetHistoMeanInRange(TH1D* histo, 
				      const double xMin, 
				      const double xMax)
 //****************************************************************************
 {
   
  int xBinMin = histo->FindBin(xMin);
  int xBinMax = histo->FindBin(xMax);
  double arithmeticMean = 0.0;
  int nValues = 0;

  for(int i = xBinMin; i < xBinMax; i++, nValues++){
    
    double tmp = histo->GetBinContent(i);
    arithmeticMean = arithmeticMean + tmp;
  }
  arithmeticMean = arithmeticMean/nValues;

  return arithmeticMean;
 }



 //****************************************************************************
bool Tracking::IsWithinTrackerAcceptance(double pt, 
					 double eta, 
					 double z0)
 //****************************************************************************
 {
   
   // cout << "pt | eta | z0 = " << pt << " | " << eta << " | " << z0 << endl;
   bool bWithinPt    = ( pt >= _pt_acceptance );    //fixme: needed?
   bool bWithinEta   = ( fabs( eta ) <= _eta_acceptance );
   bool bWithinZ0    = ( fabs( z0 ) <= _z0_acceptance );

   return bWithinPt * bWithinEta * bWithinZ0;
 }



//****************************************************************************
bool Tracking::IsWithinEtaRegion(string etaRegion, 
				 double eta)
//****************************************************************************
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


//****************************************************************************
bool Tracking::IsWithinPtRange(string ptRange, 
			       double pt)
//****************************************************************************
{

  bool bWithinPtRange = false;
  if ( ptRange.compare("Low") == 0 )         bWithinPtRange = (pt <= _pt_L);
  else if ( ptRange.compare("Middle") == 0 ) bWithinPtRange = (pt <= _pt_H && pt > _pt_L);
  else if ( ptRange.compare("High") == 0 )   bWithinPtRange = (pt > _pt_H);
  else{
    cout << "E R R O R ! Tracking::IsWithinPtRange(...) - Invalid pt range type \"" << ptRange << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinPtRange;  
}

#endif // Tracking_cxx

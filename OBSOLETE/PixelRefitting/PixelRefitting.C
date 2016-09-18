#ifndef PixelRefitting_cxx
#define PixelRefitting_cxx

#include "PixelRefitting.h"
#include "../utilities/constants.h"
#include "TFitResult.h"
#include "TF1.h"

// Re-fitted tracks. Description:
//  1) Loop over all Tracking Particles (TPs)
//  2) For current TP get the index of the corresponding TTrack (if any)
//  3) Demand that the corresponding TTrack match is unique
//  4) Get the TTPixelTrack corresponding the the TTrack which was uniquely matched to the current TP
//  5) Get the Candidate Pixel Hits of the TTPixelTrack
//  6) Re-fit the Candidate Pixel Hits of the TTPixelTrack to get a Re-fitted TTPixelTrack
//  7) Apply cuts to TTPixelTracks (kinematics, pixel hits)
//  8) Fill Candidate Pixel Hit Histos, Pixel Hit Histos, Shared Pixel Hit Histos
//  9) Fill matched refitted TTPixelTrack property Histograms
// 10) Fill matched refitted PixTrackMatchedTPHistos(pixFitTk, tp_Index);
// 11) Fill residual histograms (TP value - TTPixelTrack value), 

//****************************************************************************
void PixelRefitting::InitVars(void)
//****************************************************************************
{
  
  // Options
  ///////////////////
  bVerbose              = false;
  bPrintResolutions     = false;
  bPrintResolutions_Incl= false;
  bPrintResidualFitInfo = false;
  bSaveResidualFitInfo  = false;
  

  // TP Cuts
  ///////////////////
  TP_PT_MIN  =  0.2;  // default =  +0.2
  TP_ETA_MAX =  2.5;  // default =  +2.5
  TP_Z0_MAX  = 30.0;  // default = +30.0
  TP_DXY_MAX =  1.0;  // default =  +1.0
  TP_NMATCH  =  1.0;  // default =  +1.0

  
  // Pixel Hit Cuts
  ///////////////////
  TK_INDEX_MIN   = 0;
  TK_PT_MIN      = 1.0; // default = +1.5
  TK_PT_MAX      = 1e6; // default = +1e6
  TK_ABSETA_MIN  = 0.0; // default = +0.0
  TK_ABSETA_MAX  = 1e6; // default = +1e6
  TK_PIXHITS_MIN = -1;  // default = -1
  TK_PIXHITS_MAX = 1e6; // default = +1e6
  // Hit Requirement: L1, L1/D1 etc..
  // TK_PIXHITS_TYPE.push_back(+1);
  // TK_PIXHITS_TYPE.push_back(-1);

  
  // Pixel Hits/Patterns
  /////////////////// 
  TK_CANDPIXHITS_MIN = -1;   // default = -1
  TK_CANDPIXHITS_MAX = 1e6;  // default = +1e6
  // TK_PIXHITS_PATTERNS.push_back(7);  // L1+L2+L3
  // TK_PIXHITS_PATTERNS.push_back(11); // L1+L2+L4
  // TK_PIXHITS_PATTERNS.push_back(13); // L1+L3+L4
  // TK_PIXHITS_PATTERNS.push_back(14); // L2+L3+L4
  // TK_PIXHITS_PATTERNS.push_back(15); // L1+L2+L3+L4
  // TK_PIXHITS_PATTERNS.push_back(23);
  // TK_PIXHITS_PATTERNS.push_back(51);
  // TK_PIXHITS_PATTERNS.push_back(112);
  // TK_PIXHITS_PATTERNS.push_back(113);
  
  
  // Settings
  ///////////////////
  _pt_L              =   6.0; //  5.0
  _pt_H              =  16.0; // 15.0
  _eta_C             =   0.8; 
  _eta_F             =   1.6;
  _chiSq_overflow    = 100.0;
  _redChiSq_overflow =  15.0;
  _fromCentimetresToMicrometers = 10000.0;  // 1) cm->m: 100 2) m->micrometers*1000000;
  _fitSignificanceLevel = 0.01;
  

  // Histogram Binning
  ///////////////////
  PT_MAX       = 50.0;
  PT_BINS      = 10;
  PT_BINWIDTH  = PT_MAX/(double)PT_BINS;
  ETA_MAX      = 2.4; // default: 2.5 
  ETA_BINS     = 12;  // default: 25
  ETA_BINWIDTH = ETA_MAX/(double)ETA_BINS;
  
  return;
}



//****************************************************************************
void PixelRefitting::Loop()
//****************************************************************************
{
  
  if (fChain == 0) return;
  Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : std::min((int)fChain->GetEntries(), MaxEvents);
  
  ///////////////////////////////////////////////////////////
  /// Initialisations
  ///////////////////////////////////////////////////////////
  Long64_t nbytes = 0, nb = 0;
  L1PixelTrackFit f(3.8112); // Bz in Tesla
  InitVars();
  BookHistos();  


  ///////////////////////////////////////////////////////////
  /// Inform User of Settings & Cuts
  ///////////////////////////////////////////////////////////
  auxTools.PrintPSets(fChain);
  PrintSettings();
  PrintTPCuts();
  PrintPixTrackCuts();
  

  ///////////////////////////////////////////////////////////
  // For-loop: Entries
  ///////////////////////////////////////////////////////////
  for (Int_t jentry=0; jentry < nEntries; jentry++){
      
    // Init loop variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Get all pixel hits for each event
    vector<TVector3> pixHits_XYZ;
    vector<int> pixHits_TTPixelTrackIndex;
    s->GetPixTrackAllHits(pixHits_XYZ, pixHits_TTPixelTrackIndex);

    if (bVerbose) cout << "Event " << jentry << "/" << nEntries << ":" << endl;   

    ///////////////////////////////////////////////////////////
    // For-loop: TPs
    ///////////////////////////////////////////////////////////
    for (Int_t tp_Index = 0; tp_Index < (int) TP_Pt->size(); tp_Index++) {
     
      // Get TP properties
      double tp_Pt  = TP_Pt->at(tp_Index);
      double tp_Eta = TP_Eta->at(tp_Index);
      double tp_Phi = TP_Phi->at(tp_Index);
      double tp_x0  = TP_POCAx->at(tp_Index);
      double tp_y0  = TP_POCAy->at(tp_Index);
      double tp_z0  = TP_POCAz->at(tp_Index);
      double tp_d0  = tp->GetD0(tp_Index);
      int tp_NMatch = TP_NMatch->at(tp_Index);
      double tp_dxy = sqrt(tp_x0*tp_x0 + tp_y0*tp_y0);      

      // Apply TP cuts
      if ( !IsWithinTrackerAcceptance(tp_Pt, tp_Eta, tp_z0, tp_dxy) ) continue;

      // Fill Histograms
      h_tp_pt ->Fill(tp_Pt);
      h_tp_eta->Fill(tp_Eta);
      if (tp_Pt < _pt_L) h_tp_pt_L->Fill(tp_Pt);
      h_tp_phi->Fill(tp_Phi);
      h_tp_z0 ->Fill(tp_z0);
      h_tp_d0 ->Fill(tp_d0 * _fromCentimetresToMicrometers);
      h_tp_dxy->Fill(tp_dxy * _fromCentimetresToMicrometers);
      FillPtRegionsHistos(tp_Pt, tp_Eta, h_tp_eta_L, h_tp_eta_M, h_tp_eta_H);
      FillEtaRegionsHistos(tp_Eta, tp_Pt , h_tp_pt_C, h_tp_pt_I, h_tp_pt_F);

      
      ///////////////////////////////////////////////////////////
      // TPs for TTTrack Matching
      ///////////////////////////////////////////////////////////
      // Was the match unique? If you don't want to use this set a negative value
      if (TP_NMATCH > -1){
	if(tp_NMatch != TP_NMATCH) continue;
      }
      
          
      ///////////////////////////////////////////////////////////
          // TTPixelTracks (with unique matching with Tracking Particles)
      ///////////////////////////////////////////////////////////
      double sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0;
      // Get the TP-Corresponding TTTrack
      int tk_Index        = TP_TTTrackIndex->at(tp_Index);
      double tk_Pt        = L1Tks_Pt ->at(tk_Index);
      double tk_Eta       = L1Tks_Eta->at(tk_Index);
      double tk_Phi       = L1Tks_Phi->at(tk_Index);
      double tk_RInv      = L1Tks_RInv->at(tk_Index);
      double tk_t         = sinh(tk_Eta);
      double tk_d0        = s->GetTrackD0 (tk_Index);
      double tk_z0        = L1Tks_POCAz->at(tk_Index);
      double tk_ChiSq     = L1Tks_ChiSquared->at(tk_Index);
      double tk_RedChiSq  = s->GetTrackRedChiSq(tk_Index);

      // Get the TP-Corresponding TTPixelTrack (from the TP-Corresponding TTTrack)
      int pixTk_Index     = s->GetPixelIndexOfTrack(tk_Index);
      if (pixTk_Index < TK_INDEX_MIN) continue;
            
      // Get TTPixelTrack variables
      if (0) s->PrintPixTrackProperties(pixTk_Index, false);

      // Get the Candidate Pixel Hits of TP-Corresponding TTPixelTracks
      vector<double> candPixHits_X    = L1PixTks_CandPixHits_X->at(pixTk_Index);
      vector<double> candPixHits_Y    = L1PixTks_CandPixHits_Y->at(pixTk_Index);
      vector<double> candPixHits_Z    = L1PixTks_CandPixHits_Z->at(pixTk_Index);
      vector<double> candPixHits_R    = L1PixTks_CandPixHits_R->at(pixTk_Index);
      vector<double> candPixHits_Phi  = L1PixTks_CandPixHits_Phi->at(pixTk_Index);
      vector< int >  candPixHits_Type = L1PixTks_CandPixHits_Type->at(pixTk_Index);
      const int pixTk_nCandPixHits    = (int) candPixHits_X.size();

      // Refit the TTPixelTrack using the Candidate Pixel Hits of the TP-Corresponding TTPixelTracks
      TTPixelTrack pixFitTk =  f.FitPixelTrack(tk_RInv, tk_Phi, tk_d0, tk_t, tk_z0, candPixHits_X, candPixHits_Y, candPixHits_Z, candPixHits_Type);
      
      ///////////////////////////////////////////////////////////      
      // Re-fitted TTPixelTrack
      ///////////////////////////////////////////////////////////
      double pixTk_Pt       = pixFitTk.getMomentum().Perp();
      double pixTk_Eta      = pixFitTk.getMomentum().Eta();
      int tk_PixHits        = pixFitTk.getNhit();
      int pixHits_Pattern   = pixFitTk.getPixelHitsPattern();      

      // Acceptance
      if(pixTk_Pt < TK_PT_MIN) continue;
      if(pixTk_Pt > TK_PT_MAX) continue;
      if(fabs(pixTk_Eta) < TK_ABSETA_MIN) continue;
      if(fabs(pixTk_Eta) > TK_ABSETA_MAX) continue;
      
      // Candidate pixel hits?     
      if( pixTk_nCandPixHits < TK_CANDPIXHITS_MIN ) continue;
      if( pixTk_nCandPixHits > TK_CANDPIXHITS_MAX ) continue;

      // Fitted pixel hits
      if(tk_PixHits < TK_PIXHITS_MIN) continue;
      if(tk_PixHits > TK_PIXHITS_MAX) continue;
      if(TK_PIXHITS_TYPE.size() > 0){
	bool bHasMustHitLayerOrDisk = s->HasMustPixelHitInLayeOrDisk(pixTk_Index, TK_PIXHITS_TYPE);
	if (!bHasMustHitLayerOrDisk) continue;
      }

      // Hit pattern Cuts
      bool bFoundAllowedHitPattern =  s->HasMustPixelHitPattern(pixHits_Pattern, TK_PIXHITS_PATTERNS);
      if(!bFoundAllowedHitPattern) continue;
      
      // Fill Pixel Hit Histograms
      FillCandPixelHitHistos(pixFitTk, false);
      FillPixelHitHistos(pixFitTk, false);
      FillSharedPixelHitHistos(pixTk_Index, tp_Index, pixHits_XYZ, pixHits_TTPixelTrackIndex, false); // convert to accomodate TTPixelTrack

      // Fill TP<-->Track Matched Histograms
      FillTPMatchedPixTrackHistos(pixFitTk);
      FillPixTrackMatchedTPHistos(pixFitTk, tp_Index);

      // Print TP & TTPixelTrack properties
      if (bVerbose) tp->PrintTPProperties(tp_Index);
      if (0) s->PrintPixTrackProperties(pixTk_Index, false);
      if (0) pixFitTk.PrintProperties();
      
    } // For-loop: TPs

    if (!bVerbose) auxTools.ProgressBar(jentry, nEntries, 100, 150);
    
  }// For-loop: Entries
  cout << "\n" << endl;


  ///////////////////  ///////////////////  ///////////////////
  // After ALL cuts (TP: acceptance & unique matching. TTPixelTrack: sanity, pT, eta, candidate pixel hits, used pixel hits, 
  ///////////////////  ///////////////////  ///////////////////
  
  // Efficiency Histograms
  MakeEfficiencyHisto(h_efficiency_pt   , h_match_tp_pt   , h_tp_pt   );
  MakeEfficiencyHisto(h_efficiency_pt_L , h_match_tp_pt_L , h_tp_pt_L );
  MakeEfficiencyHisto(h_efficiency_pt_C , h_match_tp_pt_C , h_tp_pt_C );
  MakeEfficiencyHisto(h_efficiency_pt_I , h_match_tp_pt_I , h_tp_pt_I );
  MakeEfficiencyHisto(h_efficiency_pt_F , h_match_tp_pt_F , h_tp_pt_F );
  MakeEfficiencyHisto(h_efficiency_eta  , h_match_tp_eta  , h_tp_eta  );
  MakeEfficiencyHisto(h_efficiency_eta_L, h_match_tp_eta_L, h_tp_eta_L);
  MakeEfficiencyHisto(h_efficiency_eta_M, h_match_tp_eta_M, h_tp_eta_M);
  MakeEfficiencyHisto(h_efficiency_eta_H, h_match_tp_eta_H, h_tp_eta_H);
  MakeEfficiencyHisto(h_efficiency_phi  , h_match_tp_phi  , h_tp_phi  );
  MakeEfficiencyHisto(h_efficiency_z0   , h_match_tp_z0   , h_tp_z0   );
  MakeEfficiencyHisto(h_efficiency_d0   , h_match_tp_d0   , h_tp_d0   );

  // Finalise Histos, Print Resolutions & Write to ROOT file
  FinaliseHistos();
  PrintResolutions();
  outFile->cd();
  outFile->Write();

}


//****************************************************************************
void PixelRefitting::BookHistos(void)
//****************************************************************************
{

  // const int nEtaBins = 52;
  const int nEtaBins = 26;
  histoTools.BookHisto_1D( h_tp_pt   , "tp_pt"   ,   100,   0.0,  +100.0 );
  histoTools.BookHisto_1D( h_tp_pt_L , "tp_pt_L" ,    50,   0.0,    +5.0 );
  histoTools.BookHisto_1D( h_tp_pt_C , "tp_pt_C" ,   100,   0.0,  +100.0 );
  histoTools.BookHisto_1D( h_tp_pt_I , "tp_pt_I" ,   100,   0.0,  +100.0 );
  histoTools.BookHisto_1D( h_tp_pt_F , "tp_pt_F" ,   100,   0.0,  +100.0 );
  histoTools.BookHisto_1D( h_tp_eta  , "tp_eta"  ,    26,  -2.6,    +2.6 );
  histoTools.BookHisto_1D( h_tp_eta_L, "tp_eta_L",    26,  -2.6,    +2.6 );
  histoTools.BookHisto_1D( h_tp_eta_M, "tp_eta_M",    26,  -2.6,    +2.6 );
  histoTools.BookHisto_1D( h_tp_eta_H, "tp_eta_H",    26,  -2.6,    +2.6 );
  histoTools.BookHisto_1D( h_tp_phi  , "tp_phi"  ,    64,  -3.2,    +3.2 );
  histoTools.BookHisto_1D( h_tp_z0   , "tp_z0"   ,    80, -20.0,   +20.0 );
  histoTools.BookHisto_1D( h_tp_d0   , "tp_d0"   ,   100, -50.0,   +50.0 ); //micrometers
  histoTools.BookHisto_1D( h_tp_dxy  , "tp_dxy"  ,   100, -50.0,   +50.0 ); //micrometers

  histoTools.BookHisto_1D( h_match_tp_pt   , "match_tp_pt"   , 100,   0.0, +100.0 );
  histoTools.BookHisto_1D( h_match_tp_pt_L , "match_tp_pt_L" ,  50,   0.0,   +5.0 );
  histoTools.BookHisto_1D( h_match_tp_pt_C , "match_tp_pt_C" , 100,   0.0, +100.0 );
  histoTools.BookHisto_1D( h_match_tp_pt_I , "match_tp_pt_I" , 100,   0.0, +100.0 );
  histoTools.BookHisto_1D( h_match_tp_pt_F , "match_tp_pt_F" , 100,   0.0, +100.0 );
  histoTools.BookHisto_1D( h_match_tp_eta  , "match_tp_eta"  ,  26,  -2.6,   +2.6 );
  histoTools.BookHisto_1D( h_match_tp_eta_L, "match_tp_eta_L",  26,  -2.6,   +2.6 );
  histoTools.BookHisto_1D( h_match_tp_eta_M, "match_tp_eta_M",  26,  -2.6,   +2.6 );
  histoTools.BookHisto_1D( h_match_tp_eta_H, "match_tp_eta_H",  26,  -2.6,   +2.6 );
  histoTools.BookHisto_1D( h_match_tp_phi  , "match_tp_phi"  ,  64,  -3.2,   +3.2 );
  histoTools.BookHisto_1D( h_match_tp_z0   , "match_tp_z0"   ,  80, -20.0,  +20.0 );
  histoTools.BookHisto_1D( h_match_tp_d0   , "match_tp_d0"   , 100, -50.0,  +50.0 ); //micrometers

  // efficiency histograms
  histoTools.BookHisto_1D( h_efficiency_pt   , "efficiency_pt"   , 100,   0.0,  +100.0  );
  histoTools.BookHisto_1D( h_efficiency_pt_L , "efficiency_pt_L" ,  50,   0.0,    +5.0  );
  histoTools.BookHisto_1D( h_efficiency_pt_C , "efficiency_pt_C" , 100,   0.0,  +100.0  );
  histoTools.BookHisto_1D( h_efficiency_pt_I , "efficiency_pt_I" , 100,   0.0,  +100.0  ); 
  histoTools.BookHisto_1D( h_efficiency_pt_F , "efficiency_pt_F" , 100,   0.0,  +100.0  );
  histoTools.BookHisto_1D( h_efficiency_eta  , "efficiency_eta"  ,  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_eta_L, "efficiency_eta_L",  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_eta_M, "efficiency_eta_M",  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_eta_H, "efficiency_eta_H",  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_phi  , "efficiency_phi"  ,  64,  -3.2,    +3.2  );
  histoTools.BookHisto_1D( h_efficiency_z0   , "efficiency_z0"   ,  80, -20.0,   +20.0  );
  histoTools.BookHisto_1D( h_efficiency_d0   , "efficiency_d0"   , 100,-50.0,   +50.0  );

  histoTools.BookHisto_1D( h_match_trk_nstub  , "match_trk_nstub"  , 15, -0.5, 14.5);
  histoTools.BookHisto_1D( h_match_trk_nstub_C, "match_trk_nstub_C", 15, -0.5, 14.5);
  histoTools.BookHisto_1D( h_match_trk_nstub_I, "match_trk_nstub_I", 15, -0.5, 14.5);
  histoTools.BookHisto_1D( h_match_trk_nstub_F, "match_trk_nstub_F", 15, -0.5, 14.5);

  // chi2 histograms (last bin is an overflow bin)
  histoTools.BookHisto_1D( h_match_trk_chi2     , "match_trk_chi2"    , 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_L   , "match_trk_chi2_L"  , 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_M   , "match_trk_chi2_M"  , 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_H   , "match_trk_chi2_H"  , 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_C_L , "match_trk_chi2_C_L", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_I_L , "match_trk_chi2_I_L", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_F_L , "match_trk_chi2_F_L", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_C_M , "match_trk_chi2_C_M", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_I_M , "match_trk_chi2_I_M", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_F_M , "match_trk_chi2_F_M", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_C_H , "match_trk_chi2_C_H", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_I_H , "match_trk_chi2_I_H", 100, 0, 100);
  histoTools.BookHisto_1D( h_match_trk_chi2_F_H , "match_trk_chi2_F_H", 100, 0, 100);

  // chi2/dof histograms (lastbin is an overflow bin)
  histoTools.BookHisto_1D( h_match_trk_chi2_dof    , "match_trk_chi2_dof"    , 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_L  , "match_trk_chi2_dof_L"  , 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_M  , "match_trk_chi2_dof_M"  , 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_H  , "match_trk_chi2_dof_H"  , 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_C_L, "match_trk_chi2_dof_C_L", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_I_L, "match_trk_chi2_dof_I_L", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_F_L, "match_trk_chi2_dof_F_L", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_C_M, "match_trk_chi2_dof_C_M", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_I_M, "match_trk_chi2_dof_I_M", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_F_M, "match_trk_chi2_dof_F_M", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_C_H, "match_trk_chi2_dof_C_H", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_I_H, "match_trk_chi2_dof_I_H", 150, 0, 15);
  histoTools.BookHisto_1D( h_match_trk_chi2_dof_F_H, "match_trk_chi2_dof_F_H", 150, 0, 15);

  // resolution histograms
  histoTools.BookHisto_1D( h_residual_pt     , "residual_pt"     , 100, -5.0  , 5.0  );
  histoTools.BookHisto_1D( h_residual_pt_C   , "residual_pt_C"   , 100, -5.0  , 5.0  );
  histoTools.BookHisto_1D( h_residual_pt_I   , "residual_pt_I"   , 100, -5.0  , 5.0  );
  histoTools.BookHisto_1D( h_residual_pt_F   , "residual_pt_F"   , 100, -5.0  , 5.0  );

  histoTools.BookHisto_1D( h_residual_ptRel  , "residual_ptRel"  , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_ptRel_C, "residual_ptRel_C", 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_ptRel_I, "residual_ptRel_I", 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_ptRel_F, "residual_ptRel_F", 1000, -1.0  , 1.0  );

  histoTools.BookHisto_1D( h_residual_eta    , "residual_eta"    , 400, -0.01 , 0.01 );
  histoTools.BookHisto_1D( h_residual_eta_C  , "residual_eta_C"  , 200, -0.01 , 0.01 );
  histoTools.BookHisto_1D( h_residual_eta_I  , "residual_eta_I"  , 200, -0.01 , 0.01 );
  histoTools.BookHisto_1D( h_residual_eta_F  , "residual_eta_F"  , 200, -0.01 , 0.01 );

  histoTools.BookHisto_1D( h_residual_phi    , "residual_phi"    , 100, -0.005, 0.005);
  histoTools.BookHisto_1D( h_residual_phi_C  , "residual_phi_C"  , 100, -0.005, 0.005);
  histoTools.BookHisto_1D( h_residual_phi_I  , "residual_phi_I"  , 100, -0.005, 0.005);
  histoTools.BookHisto_1D( h_residual_phi_F  , "residual_phi_F"  , 100, -0.005, 0.005);

  histoTools.BookHisto_1D( h_residual_z0     , "residual_z0"     , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_z0_C   , "residual_z0_C"   , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_z0_I   , "residual_z0_I"   , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_z0_F   , "residual_z0_F"   , 1000, -1.0  , 1.0  );

  histoTools.BookHisto_1D( h_residual_d0     , "residual_d0"     , 500, -0.5  , 0.5  );
  histoTools.BookHisto_1D( h_residual_d0_C   , "residual_d0_C"   , 500, -0.5  , 0.5  );
  histoTools.BookHisto_1D( h_residual_d0_I   , "residual_d0_I"   , 500, -0.5  , 0.5  );
  histoTools.BookHisto_1D( h_residual_d0_F   , "residual_d0_F"   , 500, -0.5  , 0.5  );

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
    
    histoTools.BookHisto_1D( h_residualVsPt_pt[i]  , "residualVsPt_pt_"   + ptrange[i], 100, -5.0, 5.0);
    histoTools.BookHisto_1D( h_residualVsPt_pt_C[i], "residualVsPt_pt_C_" + ptrange[i], 100, -5.0, 5.0);
    histoTools.BookHisto_1D( h_residualVsPt_pt_I[i], "residualVsPt_pt_I_" + ptrange[i], 100, -5.0, 5.0);
    histoTools.BookHisto_1D( h_residualVsPt_pt_F[i], "residualVsPt_pt_F_" + ptrange[i], 100, -5.0, 5.0);

    // restictive range: -0.15 to 0.15
    histoTools.BookHisto_1D( h_residualVsPt_ptRel[i]  , "residualVsPt_ptRel_"  + ptrange[i], 400/nBinsXFactor, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsPt_ptRel_C[i], "residualVsPt_ptRel_C_"+ ptrange[i], 200/nBinsXFactor, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsPt_ptRel_I[i], "residualVsPt_ptRel_I_"+ ptrange[i], 200/nBinsXFactor, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsPt_ptRel_F[i], "residualVsPt_ptRel_F_"+ ptrange[i], 200/nBinsXFactor, -0.50, +0.50);

    histoTools.BookHisto_1D( h_residualVsPt_eta[i]  , "residualVsPt_eta_"    + ptrange[i], 400/nBinsXFactor, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsPt_eta_C[i], "residualVsPt_eta_C_"  + ptrange[i], 200/nBinsXFactor, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsPt_eta_I[i], "residualVsPt_eta_I_"  + ptrange[i], 200/nBinsXFactor, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsPt_eta_F[i], "residualVsPt_eta_F_"  + ptrange[i], 200/nBinsXFactor, -0.01, +0.01);

    histoTools.BookHisto_1D( h_residualVsPt_phi[i]  , "residualVsPt_phi_"   + ptrange[i], 100, -0.005, +0.005);
    histoTools.BookHisto_1D( h_residualVsPt_phi_C[i], "residualVsPt_phi_C_" + ptrange[i], 100, -0.005, +0.005);
    histoTools.BookHisto_1D( h_residualVsPt_phi_I[i], "residualVsPt_phi_I_" + ptrange[i], 100, -0.005, +0.005);
    histoTools.BookHisto_1D( h_residualVsPt_phi_F[i], "residualVsPt_phi_F_" + ptrange[i], 100, -0.005, +0.005);

    histoTools.BookHisto_1D( h_residualVsPt_z0[i]  , "residualVsPt_z0_"   + ptrange[i], 500/nBinsXFactor, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsPt_z0_C[i], "residualVsPt_z0_C_" + ptrange[i], 500/nBinsXFactor, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsPt_z0_I[i], "residualVsPt_z0_I_" + ptrange[i], 500/nBinsXFactor, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsPt_z0_F[i], "residualVsPt_z0_F_" + ptrange[i], 250/nBinsXFactor, -0.5, +0.5);

    histoTools.BookHisto_1D( h_residualVsPt_d0[i]  , "residualVsPt_d0_"   + ptrange[i], 200/nBinsXFactor, -0.1, +0.1);
    histoTools.BookHisto_1D( h_residualVsPt_d0_C[i], "residualVsPt_d0_C_" + ptrange[i], 200/nBinsXFactor, -0.1, +0.1);
    histoTools.BookHisto_1D( h_residualVsPt_d0_I[i], "residualVsPt_d0_I_" + ptrange[i], 200/nBinsXFactor, -0.1, +0.1);
    histoTools.BookHisto_1D( h_residualVsPt_d0_F[i], "residualVsPt_d0_F_" + ptrange[i], 100/nBinsXFactor, -0.1, +0.1);
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
	
    histoTools.BookHisto_1D( h_residualVsEta_pt[i]  , "residualVsEta_pt_"   + etarange[i], 100, -5.00, +5.00);
    histoTools.BookHisto_1D( h_residualVsEta_pt_L[i], "residualVsEta_pt_L_" + etarange[i], 100, -0.25, +0.25);
    histoTools.BookHisto_1D( h_residualVsEta_pt_M[i], "residualVsEta_pt_M_" + etarange[i], 100, -0.25, +0.25);
    histoTools.BookHisto_1D( h_residualVsEta_pt_H[i], "residualVsEta_pt_H_" + etarange[i], 100, -0.25, +0.25);

    histoTools.BookHisto_1D( h_residualVsEta_ptRel[i]  , "residualVsEta_ptRel_"  + etarange[i], 200, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsEta_ptRel_L[i], "residualVsEta_ptRel_L_"+ etarange[i], 200, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsEta_ptRel_M[i], "residualVsEta_ptRel_M_"+ etarange[i], 200, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsEta_ptRel_H[i], "residualVsEta_ptRel_H_"+ etarange[i], 200, -0.50, +0.50);

    histoTools.BookHisto_1D( h_residualVsEta_eta[i]  , "residualVsEta_eta_"   + etarange[i], 400/nBinsXFactor  , -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_eta_L[i], "residualVsEta_eta_L_" + etarange[i], 200/nBinsXFactor_L, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_eta_M[i], "residualVsEta_eta_M_" + etarange[i], 400/nBinsXFactor_M, -0.01, +0.01); 
    histoTools.BookHisto_1D( h_residualVsEta_eta_H[i], "residualVsEta_eta_H_" + etarange[i], 400/nBinsXFactor_H, -0.01, +0.01);

    histoTools.BookHisto_1D( h_residualVsEta_phi[i]  , "residualVsEta_phi_"   + etarange[i], 100, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_phi_L[i], "residualVsEta_phi_L_" + etarange[i], 100, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_phi_M[i], "residualVsEta_phi_M_" + etarange[i], 100, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_phi_H[i], "residualVsEta_phi_H_" + etarange[i], 100, -0.01, +0.01);

    histoTools.BookHisto_1D( h_residualVsEta_z0[i]  , "residualVsEta_z0_"   + etarange[i], 1000/nBinsXFactor  , -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_z0_L[i], "residualVsEta_z0_L_" + etarange[i], 1000/nBinsXFactor_L, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_z0_M[i], "residualVsEta_z0_M_" + etarange[i], 1000/nBinsXFactor_M, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_z0_H[i], "residualVsEta_z0_H_" + etarange[i], 1000/nBinsXFactor_H, -0.5, +0.5);

    histoTools.BookHisto_1D( h_residualVsEta_d0[i]  , "residualVsEta_d0_"   + etarange[i], 2000/nBinsXFactor  , -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_d0_L[i], "residualVsEta_d0_L_" + etarange[i], 1000/nBinsXFactor_L, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_d0_M[i], "residualVsEta_d0_M_" + etarange[i], 1000/nBinsXFactor_M, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_d0_H[i], "residualVsEta_d0_H_" + etarange[i], 2000/nBinsXFactor_H, -0.5, +0.5);

  }

  // 2D histograms
  histoTools.BookHisto_2D( h_2d_logchi2_eta    , "2d_logchi2_eta"    , 50, -2.5, +2.5, 120, -4.0, +8.0  );
  histoTools.BookHisto_2D( h_2d_logchi2_dof_eta, "2d_logchi2_dof_eta", 50, -2.5, +2.5, 120, -4.0, +8.0  );
  histoTools.BookHisto_2D( h_2d_dz0_eta        , "2d_dz0_eta"        , 50, -2.5, +2.5, 120,  0.0, +1.2  );
  histoTools.BookHisto_2D( h_2d_dd0_eta        , "2d_dd0_eta"        , 50, -2.5, +2.5, 120,  0.0, +1.2  );
  histoTools.BookHisto_2D( h_2d_deta_eta       , "2d_deta_eta"       , 50, -2.5, +2.5, 120,  0.0, +0.012);
  histoTools.BookHisto_2D( h_2d_dphi_eta       , "2d_dphi_eta"       , 50, -2.5, +2.5, 100,  0.0, +0.010);
  histoTools.BookHisto_2D( h_2d_dpt_eta        , "2d_dpt_eta"        , 50, -2.5, +2.5, 500,  0.0, +5.0  );
  histoTools.BookHisto_2D( h_2d_dptRel_eta     , "2d_dptRel_eta"     , 50, -2.5, +2.5, 500,  0.0, +5.0  );

  // Resolution Vs pT histograms
  histoTools.BookHisto_1D( h_resolutionVsPt_pt  , "resolutionVsPt_pt"  , PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_pt_C, "resolutionVsPt_pt_C", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_pt_I, "resolutionVsPt_pt_I", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_pt_F, "resolutionVsPt_pt_F", PT_BINS, 0.0, PT_MAX);

  histoTools.BookHisto_1D( h_resolutionVsPt_ptRel  , "resolutionVsPt_ptRel"  , PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_ptRel_C, "resolutionVsPt_ptRel_C", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_ptRel_I, "resolutionVsPt_ptRel_I", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_ptRel_F, "resolutionVsPt_ptRel_F", PT_BINS, 0.0, PT_MAX);

  histoTools.BookHisto_1D( h_resolutionVsPt_eta  , "resolutionVsPt_eta"  , PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_eta_C, "resolutionVsPt_eta_C", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_eta_I, "resolutionVsPt_eta_I", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_eta_F, "resolutionVsPt_eta_F", PT_BINS, 0.0, PT_MAX);

  histoTools.BookHisto_1D( h_resolutionVsPt_phi  , "resolutionVsPt_phi"  , PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_phi_C, "resolutionVsPt_phi_C", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_phi_I, "resolutionVsPt_phi_I", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_phi_F, "resolutionVsPt_phi_F", PT_BINS, 0.0, PT_MAX);

  histoTools.BookHisto_1D( h_resolutionVsPt_z0  , "resolutionVsPt_z0"  , PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_z0_C, "resolutionVsPt_z0_C", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_z0_I, "resolutionVsPt_z0_I", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_z0_F, "resolutionVsPt_z0_F", PT_BINS, 0.0, PT_MAX);

  histoTools.BookHisto_1D( h_resolutionVsPt_d0  , "resolutionVsPt_d0"  , PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_d0_C, "resolutionVsPt_d0_C", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_d0_I, "resolutionVsPt_d0_I", PT_BINS, 0.0, PT_MAX);
  histoTools.BookHisto_1D( h_resolutionVsPt_d0_F, "resolutionVsPt_d0_F", PT_BINS, 0.0, PT_MAX);

  // Resolution Vs eta histograms
  histoTools.BookHisto_1D( h_resolutionVsEta_pt   , "resolutionVsEta_pt"   , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_pt_L , "resolutionVsEta_pt_L" , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_pt_M , "resolutionVsEta_pt_M" , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_pt_H , "resolutionVsEta_pt_H" , ETA_BINS, 0.0, ETA_MAX);

  histoTools.BookHisto_1D( h_resolutionVsEta_ptRel  , "resolutionVsEta_ptRel"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_ptRel_L, "resolutionVsEta_ptRel_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_ptRel_M, "resolutionVsEta_ptRel_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_ptRel_H, "resolutionVsEta_ptRel_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools.BookHisto_1D( h_resolutionVsEta_eta  , "resolutionVsEta_eta"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_eta_L, "resolutionVsEta_eta_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_eta_M, "resolutionVsEta_eta_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_eta_H, "resolutionVsEta_eta_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools.BookHisto_1D( h_resolutionVsEta_phi  , "resolutionVsEta_phi"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_phi_L, "resolutionVsEta_phi_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_phi_M, "resolutionVsEta_phi_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_phi_H, "resolutionVsEta_phi_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools.BookHisto_1D( h_resolutionVsEta_z0  , "resolutionVsEta_z0"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_z0_L, "resolutionVsEta_z0_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_z0_M, "resolutionVsEta_z0_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_z0_H, "resolutionVsEta_z0_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools.BookHisto_1D( h_resolutionVsEta_d0  , "resolutionVsEta_d0"  , ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_d0_L, "resolutionVsEta_d0_L", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_d0_M, "resolutionVsEta_d0_M", ETA_BINS, 0.0, ETA_MAX);
  histoTools.BookHisto_1D( h_resolutionVsEta_d0_H, "resolutionVsEta_d0_H", ETA_BINS, 0.0, ETA_MAX);

  histoTools.BookHisto_1D( h_pixTk_pixHits_N            , "pixTk_pixHits_N"            ,   30,   -0.5,  +29.5 );
  histoTools.BookHisto_1D( h_pixTk_pixHits_Rho          , "pixTk_pixHits_Rho"          , 1000,   +0.0,  +20.0 );
  histoTools.BookHisto_1D( h_pixTk_pixHits_Z            , "pixTk_pixHits_Z"            ,  240,  -60.0,  +60.0 );
  histoTools.BookHisto_1D( h_pixTk_pixHits_Type         , "pixTk_pixHits_Type"         ,   10,   -4.5,   +5.5 );
  histoTools.BookHisto_1D( h_pixTk_pixHits_Pattern      , "pixTk_pixHits_Pattern"      ,  150,   -0.5, +149.5 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_EtaVsN       , "pixTk_pixHits_EtaVsN"       ,  500,   -5.0,   +5.0,   30, -0.5, +29.5 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_ZVsRho       , "pixTk_pixHits_ZVsRho"       ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_ZVsRho_C     , "pixTk_pixHits_ZVsRho_C"     ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_ZVsRho_I     , "pixTk_pixHits_ZVsRho_I"     ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_ZVsRho_F     , "pixTk_pixHits_ZVsRho_F"     ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_ZVsRho_EtaGE2, "pixTk_pixHits_ZVsRho_EtaGE2",  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_XVsY         , "pixTk_pixHits_XVsY"         ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_XVsY_C       , "pixTk_pixHits_XVsY_C"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_XVsY_I       , "pixTk_pixHits_XVsY_I"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_XVsY_F       , "pixTk_pixHits_XVsY_F"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_XVsY_EtaGE2  , "pixTk_pixHits_XVsY_EtaGE2"  ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_pixHits_PatternVsEta , "pixTk_pixHits_PatternVsEta" ,  150,   -0.5, +149.5, 2*ETA_BINS, 0.0, ETA_MAX); 
  
  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_N            , "pixTk_sharedPixHits_N"            ,   30,   -0.5,  +29.5 );
  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_Rho          , "pixTk_sharedPixHits_Rho"          , 1000,   +0.0,  +20.0 );
  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_Z            , "pixTk_sharedPixHits_Z"            ,  240,  -60.0,  +60.0 );
  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_Type         , "pixTk_sharedPixHits_Type"         ,   10,   -4.5,   +5.5 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_EtaVsN       , "pixTk_sharedPixHits_EtaVsN"       ,  500,   -5.0,   +5.0,   30, -0.5, +29.5 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho       , "pixTk_sharedPixHits_ZVsRho"       ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_C     , "pixTk_sharedPixHits_ZVsRho_C"     ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_I     , "pixTk_sharedPixHits_ZVsRho_I"     ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_F     , "pixTk_sharedPixHits_ZVsRho_F"     ,  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_ZVsRho_EtaGE2, "pixTk_sharedPixHits_ZVsRho_EtaGE2",  600,  -60.0,  +60.0,  100, +0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_XVsY         , "pixTk_sharedPixHits_XVsY"         ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_XVsY_C       , "pixTk_sharedPixHits_XVsY_C"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_XVsY_I       , "pixTk_sharedPixHits_XVsY_I"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_XVsY_F       , "pixTk_sharedPixHits_XVsY_F"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_sharedPixHits_XVsY_EtaGE2  , "pixTk_sharedPixHits_XVsY_EtaGE2"  ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );

  histoTools.BookHisto_1D( h_pixTk_candPixHits_N            , "pixTk_candPixHits_N"            ,   30,   -0.5,  +29.5 );
  histoTools.BookHisto_1D( h_pixTk_candPixHits_Rho          , "pixTk_candPixHits_Rho"          , 1000,   +0.0,  +20.0 );
  histoTools.BookHisto_1D( h_pixTk_candPixHits_Z            , "pixTk_candPixHits_Z"            ,  240,  -60.0,  +60.0 );
  histoTools.BookHisto_1D( h_pixTk_candPixHits_Type         , "pixTk_candPixHits_Type"         ,   10,   -4.5,   +5.5 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_EtaVsN       , "pixTk_candPixHits_EtaVsN"       ,  500,   -5.0,   +5.0,   30, -0.5, +29.5  );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_ZVsRho       , "pixTk_candPixHits_ZVsRho"       ,  600,  -60.0,  +60.0,  100, + 0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_C     , "pixTk_candPixHits_ZVsRho_C"     ,  600,  -60.0,  +60.0,  100, + 0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_I     , "pixTk_candPixHits_ZVsRho_I"     ,  600,  -60.0,  +60.0,  100, + 0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_F     , "pixTk_candPixHits_ZVsRho_F"     ,  600,  -60.0,  +60.0,  100, + 0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_ZVsRho_EtaGE2, "pixTk_candPixHits_ZVsRho_EtaGE2",  600,  -60.0,  +60.0,  100, + 0.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_XVsY         , "pixTk_candPixHits_XVsY"         ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_XVsY_C       , "pixTk_candPixHits_XVsY_C"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_XVsY_I       , "pixTk_candPixHits_XVsY_I"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_XVsY_F       , "pixTk_candPixHits_XVsY_F"       ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );
  histoTools.BookHisto_2D( h_pixTk_candPixHits_XVsY_EtaGE2  , "pixTk_candPixHits_XVsY_EtaGE2"  ,  400,  -20.0,  +20.0,  400, -20.0, +20.0 );

  return;
}


//****************************************************************************
void PixelRefitting::FillHistoBinWithRMS(TH1D *hToFill,
					TH1D *hToExtractRMS,
					int binNumber)
//****************************************************************************
{
  
  hToFill->SetBinContent(binNumber, hToExtractRMS->GetRMS() );
  hToFill->SetBinError  (binNumber, hToExtractRMS->GetRMSError() );

  return;
}
  

//****************************************************************************
void PixelRefitting::FitAndFillHistoBinSL(TH1D *hToFill,
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

  // "Q"= Quiet mode (minimum printing), "S" = The result of the fit is returned in the TFitResultPtr (see below Access to the Fit Result)
  TFitResultPtr r = histoTools.FitFunctionSL( hToFit, fitFunction, "Q S ", "", GetFitAreaVector(), significanceLevel, bPrintResidualFitInfo, bSaveResidualFitInfo);

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
    histoTools.FindSymmetricFitRange(hToFit, 0.8, fitRangeLow, fitRangeHigh);
    hToFit->GetXaxis()->SetRange(fitRangeLow, fitRangeHigh);
    
    hToFill->SetBinContent(binNumber, hToFit->GetRMS() );
    hToFill->SetBinError  (binNumber, hToFit->GetRMSError() );
  }

   
  return;
}


//****************************************************************************
std::vector<double> PixelRefitting::GetFitAreaVector(const double minTotalArea)
//****************************************************************************
{

  // Maximum possible area on ONE side of gaussian is 50% (symmetric)
  double area_max = 0.5;
  
  // Fit area from 0.5 -> minTotalArea/2
  std::vector<double> v_fitArea;  
  for (int i = 0; i <= 20; i++){
 
    double area = area_max - i*0.02;

    // Break if total area (both sides) is bigger than the minimum area allowed (minTotalArea=0.50 by default)
    if (2*area < minTotalArea) break;

    v_fitArea.push_back(area);
  }

  // auxTools.PrintVector(v_fitArea);
  
  return v_fitArea;
}


//****************************************************************************
void PixelRefitting::FinaliseHistos(void)
//****************************************************************************
{
  
  // Loop over all pt bins
  for (int i=0; i < PT_BINS; i++) {

    // Cannot fit "gaus" to phi due to asymmetric distribution (lorentz-drift)
    FitAndFillHistoBinSL( h_resolutionVsPt_pt  , h_residualVsPt_pt[i]  , i+1, _fitSignificanceLevel);
    FitAndFillHistoBinSL( h_resolutionVsPt_pt_C, h_residualVsPt_pt_C[i], i+1, _fitSignificanceLevel);
    FitAndFillHistoBinSL( h_resolutionVsPt_pt_I, h_residualVsPt_pt_I[i], i+1, _fitSignificanceLevel);
    FitAndFillHistoBinSL( h_resolutionVsPt_pt_F, h_residualVsPt_pt_F[i], i+1, _fitSignificanceLevel);

    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel  , h_residualVsPt_ptRel[i]  , i+1, _fitSignificanceLevel);
    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel_C, h_residualVsPt_ptRel_C[i], i+1, _fitSignificanceLevel);
    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel_I, h_residualVsPt_ptRel_I[i], i+1, _fitSignificanceLevel);
    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel_F, h_residualVsPt_ptRel_F[i], i+1, _fitSignificanceLevel);
    
    // Fit "gaus" to eta to get resolution
    FitAndFillHistoBinSL( h_resolutionVsPt_eta  , h_residualVsPt_eta[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_eta_C, h_residualVsPt_eta_C[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_eta_I, h_residualVsPt_eta_I[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_eta_F, h_residualVsPt_eta_F[i], i+1, _fitSignificanceLevel );

    // Cannot fit "gaus" to phi due to asymmetric distribution (lorentz-drift)
    FitAndFillHistoBinSL( h_resolutionVsPt_phi  , h_residualVsPt_phi[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_phi_C, h_residualVsPt_phi_C[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_phi_I, h_residualVsPt_phi_I[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_phi_F, h_residualVsPt_phi_F[i], i+1, _fitSignificanceLevel );

    // Fit "gaus" to z0 to get resolution
    FitAndFillHistoBinSL( h_resolutionVsPt_z0  , h_residualVsPt_z0[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_z0_C, h_residualVsPt_z0_C[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_z0_I, h_residualVsPt_z0_I[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_z0_F, h_residualVsPt_z0_F[i], i+1, _fitSignificanceLevel );

    // Fit "gaus" to d0 to get resolution
    FitAndFillHistoBinSL( h_resolutionVsPt_d0  , h_residualVsPt_d0[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_d0_C, h_residualVsPt_d0_C[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_d0_I, h_residualVsPt_d0_I[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL( h_resolutionVsPt_d0_F, h_residualVsPt_d0_F[i], i+1, _fitSignificanceLevel );

  }

  // Loop over all eta bins
  for (int i=0; i < ETA_BINS; i++) {

    FitAndFillHistoBinSL(h_resolutionVsEta_pt  , h_residualVsEta_pt[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_pt_L, h_residualVsEta_pt_L[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_pt_M, h_residualVsEta_pt_M[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_pt_H, h_residualVsEta_pt_H[i], i+1, _fitSignificanceLevel );
  
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel  , h_residualVsEta_ptRel[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel_L, h_residualVsEta_ptRel_L[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel_M, h_residualVsEta_ptRel_M[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel_H, h_residualVsEta_ptRel_H[i], i+1, _fitSignificanceLevel );

    FitAndFillHistoBinSL(h_resolutionVsEta_eta  , h_residualVsEta_eta[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_eta_L, h_residualVsEta_eta_L[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_eta_M, h_residualVsEta_eta_M[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_eta_H, h_residualVsEta_eta_H[i], i+1, _fitSignificanceLevel );
    
    FitAndFillHistoBinSL(h_resolutionVsEta_phi  , h_residualVsEta_phi[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_phi_L, h_residualVsEta_phi_L[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_phi_M, h_residualVsEta_phi_M[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_phi_H, h_residualVsEta_phi_H[i], i+1, _fitSignificanceLevel );

    FitAndFillHistoBinSL(h_resolutionVsEta_z0  , h_residualVsEta_z0[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_z0_L, h_residualVsEta_z0_L[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_z0_M, h_residualVsEta_z0_M[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_z0_H, h_residualVsEta_z0_H[i], i+1, _fitSignificanceLevel );

    FitAndFillHistoBinSL(h_resolutionVsEta_d0  , h_residualVsEta_d0[i]  , i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_d0_L, h_residualVsEta_d0_L[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_d0_M, h_residualVsEta_d0_M[i], i+1, _fitSignificanceLevel );
    FitAndFillHistoBinSL(h_resolutionVsEta_d0_H, h_residualVsEta_d0_H[i], i+1, _fitSignificanceLevel );
    
  }
  
  return;
}

    

//****************************************************************************
void PixelRefitting::MakeEfficiencyHisto(TH1D* h_eff, 
					TH1D* h_match_tp, 
					TH1D* h_tp)
//****************************************************************************
{

  h_eff->Divide(h_match_tp, h_tp, 1.0, 1.0, "B");
  
  return;
}
 
 

//****************************************************************************
void PixelRefitting::PrintSettings(void)
//****************************************************************************
{

  // Inform user of settings  
  Table settings("Symbol | Definition | Units | Description | Comments", "Text");
  settings.AddRowColumn(0, "bVerbose");
  settings.AddRowColumn(0, auxTools.ToString(bVerbose) );
  settings.AddRowColumn(0, "-" );
  settings.AddRowColumn(0, "Enable to print event-by-event info");
  settings.AddRowColumn(0, "Default is 0");

  settings.AddRowColumn(1, "bPrintResolutions");
  settings.AddRowColumn(1, auxTools.ToString(bPrintResolutions) );
  settings.AddRowColumn(1, "-" );
  settings.AddRowColumn(1, "Enable to print resolutions for pT regions");
  settings.AddRowColumn(1, "Default is 0");

  settings.AddRowColumn(2, "bPrintResolutions_Incl");
  settings.AddRowColumn(2, auxTools.ToString(bPrintResolutions_Incl) );
  settings.AddRowColumn(2, "-" );
  settings.AddRowColumn(2, "Enable to print inclusive resolutions");
  settings.AddRowColumn(2, "Default is 0");

  settings.AddRowColumn(3, "bPrintResidualFitInfo");
  settings.AddRowColumn(3, auxTools.ToString(bPrintResidualFitInfo) );
  settings.AddRowColumn(3, "-" );
  settings.AddRowColumn(3, "Enable to print all info on residual fits attempted");
  settings.AddRowColumn(3, "Default is 0");

  settings.AddRowColumn(4, "bSaveResidualFitInfo");
  settings.AddRowColumn(4, auxTools.ToString(bSaveResidualFitInfo) );
  settings.AddRowColumn(4, "-" );
  settings.AddRowColumn(4, "Enable to save to txt file all info for the residual fits attempted");
  settings.AddRowColumn(4, "Default is 0");

  settings.AddRowColumn(5, "_eta_C");
  settings.AddRowColumn(5, "abs(eta) < " + auxTools.ToString(_eta_C) );
  settings.AddRowColumn(5, "-");
  settings.AddRowColumn(5, "Central Eta Region");
  settings.AddRowColumn(5, "-");

  settings.AddRowColumn(6, "_eta_I");
  settings.AddRowColumn(6, auxTools.ToString(_eta_C) + " < abs(eta) < " + auxTools.ToString(_eta_F));
  settings.AddRowColumn(6, "-");
  settings.AddRowColumn(6, "Intermediate Eta Region");
  settings.AddRowColumn(6, "-");
    
  settings.AddRowColumn(7, "_eta_F");
  settings.AddRowColumn(7, "abs(eta) > " + auxTools.ToString(_eta_F) );
  settings.AddRowColumn(7, "-");
  settings.AddRowColumn(7, "Forward Eta Region");
  settings.AddRowColumn(7, "-");
  
  settings.AddRowColumn(8, "_pt_L");
  settings.AddRowColumn(8, "pT < " + auxTools.ToString(_pt_L));
  settings.AddRowColumn(8, "GeV/c");
  settings.AddRowColumn(8, "Low pT Range");
  settings.AddRowColumn(8, "Default is 5");
    
  settings.AddRowColumn(9, "_pt_M");
  settings.AddRowColumn(9, auxTools.ToString(_pt_L) + " GeV/c < pT < " + auxTools.ToString(_pt_H) );
  settings.AddRowColumn(9, "GeV/c");
  settings.AddRowColumn(9, "Middle pT Range");
  settings.AddRowColumn(9, "Default is 5-15");
    
  settings.AddRowColumn(10, "_pt_H");
  settings.AddRowColumn(10, "pT > " + auxTools.ToString(_pt_H) );
  settings.AddRowColumn(10, "GeV/c");
  settings.AddRowColumn(10, "High pT Range");
  settings.AddRowColumn(10, "Default is 15");
  
  settings.AddRowColumn(11, "_chiSq_overflow");
  settings.AddRowColumn(11, auxTools.ToString(_chiSq_overflow) );
  settings.AddRowColumn(11, "-" );
  settings.AddRowColumn(11, "Value of overflow bin for chi-squared of L1 tracks" );
  settings.AddRowColumn(11, "Default is 100");
  
  settings.AddRowColumn(12, "_redChiSq_overflow");
  settings.AddRowColumn(12, auxTools.ToString(_redChiSq_overflow) );
  settings.AddRowColumn(12, "-" );
  settings.AddRowColumn(12, "Value of overflow bin for reduced-chi-squared of L1 tracks" );
  settings.AddRowColumn(12, "Default is 15");

  settings.AddRowColumn(13, "_fitSignificanceLevel");
  settings.AddRowColumn(13, auxTools.ToString(_fitSignificanceLevel) );
  settings.AddRowColumn(13, "-" );
  settings.AddRowColumn(13, "Value of the significance level in Gaussian Residual Fits = P(Reject Hypothesis)" );
  settings.AddRowColumn(13, "Default is 0.01");
  
  settings.Print();
  
  return;
}



//****************************************************************************
void PixelRefitting::PrintTPCuts(void)
//****************************************************************************
{

  Table tpCuts("# | Variable | Cut | Cut Value | Units | Description | Comments", "Text");
  tpCuts.AddRowColumn(0, "1");
  tpCuts.AddRowColumn(0, "TP_PT_MIN");
  tpCuts.AddRowColumn(0, " >= " );
  tpCuts.AddRowColumn(0, auxTools.ToString(TP_PT_MIN) );
  tpCuts.AddRowColumn(0, "GeV/c" );
  tpCuts.AddRowColumn(0, "Minimum pT for TP to be within acceptance" );
  tpCuts.AddRowColumn(0, "Default is 0.2 GeV/c" );

  tpCuts.AddRowColumn(1, "2");
  tpCuts.AddRowColumn(1, "TP_ETA_MAX");
  tpCuts.AddRowColumn(1, " <= " );
  tpCuts.AddRowColumn(1, auxTools.ToString(TP_ETA_MAX) );
  tpCuts.AddRowColumn(1, "-" );
  tpCuts.AddRowColumn(1, "Maximum abs(eta) value TP to be within acceptance" );
  tpCuts.AddRowColumn(1, "Default is 2.5" );

  tpCuts.AddRowColumn(2, "3");
  tpCuts.AddRowColumn(2, "TP_Z0_MAX");
  tpCuts.AddRowColumn(2, " <= " );
  tpCuts.AddRowColumn(2, auxTools.ToString(TP_Z0_MAX) );
  tpCuts.AddRowColumn(2, "cm" );
  tpCuts.AddRowColumn(2, "Maximum abs(z0) value TP to be within acceptance" );
  tpCuts.AddRowColumn(2, "Default is 30cm" );

  tpCuts.AddRowColumn(3, "4");
  tpCuts.AddRowColumn(3, "TP_DXY_MAX");
  tpCuts.AddRowColumn(3, " <= " );
  tpCuts.AddRowColumn(3, auxTools.ToString(TP_DXY_MAX) );
  tpCuts.AddRowColumn(3, "cm" );
  tpCuts.AddRowColumn(3, "Maximum dxy of TPs considered (Choose near IP");
  tpCuts.AddRowColumn(3, "Default is 1cm. Affects d0 Resolution!" );

  tpCuts.AddRowColumn(4, "5");  
  tpCuts.AddRowColumn(4, "TP_NMATCH");
  tpCuts.AddRowColumn(4, " == " );
  tpCuts.AddRowColumn(4, auxTools.ToString(TP_DXY_MAX) );
  tpCuts.AddRowColumn(4, "-" );
  tpCuts.AddRowColumn(4, "Number of matches with a Track");
  tpCuts.AddRowColumn(4, "Default is +1 (Unique). Disable is -1");
  tpCuts.Print();
  
  return;
}


//****************************************************************************
void PixelRefitting::PrintPixTrackCuts(void)
//****************************************************************************
{

  // Inform user of settings  
  Table pixTkCuts("# | Variable | Cut | Cut Value | Units | Description | Comments", "Text");
  pixTkCuts.AddRowColumn(0, "1");
  pixTkCuts.AddRowColumn(0, "TK_INDEX_MIN");
  pixTkCuts.AddRowColumn(0, " >= " );
  pixTkCuts.AddRowColumn(0, auxTools.ToString(TK_INDEX_MIN) );
  pixTkCuts.AddRowColumn(0, "-" );
  pixTkCuts.AddRowColumn(0, "Sanity check for TTPixelTracks considered" );
  pixTkCuts.AddRowColumn(0, "Default is 0" );

  pixTkCuts.AddRowColumn(1, "2");
  pixTkCuts.AddRowColumn(1, "TK_PT_MIN");
  pixTkCuts.AddRowColumn(1, " >= " );
  pixTkCuts.AddRowColumn(1, auxTools.ToString(TK_PT_MIN) );
  pixTkCuts.AddRowColumn(1, "-" );
  pixTkCuts.AddRowColumn(1, "Minimum pT for TTPixelTracks considered" );
  pixTkCuts.AddRowColumn(1, "Default is 0" );

  pixTkCuts.AddRowColumn(2, "3");
  pixTkCuts.AddRowColumn(2, "TK_PT_MAX");
  pixTkCuts.AddRowColumn(2, " <= " );
  pixTkCuts.AddRowColumn(2, auxTools.ToString(TK_PT_MAX) );
  pixTkCuts.AddRowColumn(2, "-" );
  pixTkCuts.AddRowColumn(2, "Maximum pT for TTPixelTracks considered" );
  pixTkCuts.AddRowColumn(2, "Default is +1e6" );

  pixTkCuts.AddRowColumn(3, "4");
  pixTkCuts.AddRowColumn(3, "TK_ABSETA_MIN");
  pixTkCuts.AddRowColumn(3, " >= " );
  pixTkCuts.AddRowColumn(3, auxTools.ToString(TK_ABSETA_MIN) );
  pixTkCuts.AddRowColumn(3, "-" );
  pixTkCuts.AddRowColumn(3, "Minimum abs{eta} for TTPixelTracks considered" );
  pixTkCuts.AddRowColumn(3, "Default is 0.0" );

  pixTkCuts.AddRowColumn(4, "5");
  pixTkCuts.AddRowColumn(4, "TK_ABSETA_MAX");
  pixTkCuts.AddRowColumn(4, " <= " );
  pixTkCuts.AddRowColumn(4, auxTools.ToString(TK_ABSETA_MAX) );
  pixTkCuts.AddRowColumn(4, "-" );
  pixTkCuts.AddRowColumn(4, "Maximum abs{eta} for TTPixelTracks considered" );
  pixTkCuts.AddRowColumn(4, "Default is +1e6" );
  
  pixTkCuts.AddRowColumn(5, "6");
  pixTkCuts.AddRowColumn(5, "TK_PIXHITS_MIN");
  pixTkCuts.AddRowColumn(5, " >= " );
  pixTkCuts.AddRowColumn(5, auxTools.ToString(TK_PIXHITS_MIN) );
  pixTkCuts.AddRowColumn(5, "-" );
  pixTkCuts.AddRowColumn(5, "Minimum number of Pixel Hits for TTPixelTracks" );
  pixTkCuts.AddRowColumn(5, "Default is -1" );

  pixTkCuts.AddRowColumn(6, "7");
  pixTkCuts.AddRowColumn(6, "TK_PIXHITS_MAX"); 
  pixTkCuts.AddRowColumn(6, " <= " );
  pixTkCuts.AddRowColumn(6, auxTools.ToString(TK_PIXHITS_MAX) );
  pixTkCuts.AddRowColumn(6, "-" );
  pixTkCuts.AddRowColumn(6, "Maximum number of Pixel Hits for TTPixelTracks" );
  pixTkCuts.AddRowColumn(6, "Default is +1e6" );

  pixTkCuts.AddRowColumn(7, "8");
  pixTkCuts.AddRowColumn(7, "TK_PIXHITS_TYPE");
  pixTkCuts.AddRowColumn(7, " == " );
  pixTkCuts.AddRowColumn(7, auxTools.ConvertIntVectorToString(TK_PIXHITS_TYPE) );
  pixTkCuts.AddRowColumn(7, "-" );
  pixTkCuts.AddRowColumn(7, "Pixel Hit Types for TTPixelTracks (OR)");
  pixTkCuts.AddRowColumn(7, "Default is empty" );

  pixTkCuts.AddRowColumn(8, "9");
  pixTkCuts.AddRowColumn(8, "TK_PIXHITS_PATTERNS");
  pixTkCuts.AddRowColumn(8, " == " );
  pixTkCuts.AddRowColumn(8, auxTools.ConvertIntVectorToString(TK_PIXHITS_PATTERNS) );
  pixTkCuts.AddRowColumn(8, "-" );
  pixTkCuts.AddRowColumn(8, "Bit-based Hit pattern variable" );
  pixTkCuts.AddRowColumn(8, "L1=1,L2=2,L3=4,L4=8,D1=16,D2=32,D3=64" );

  pixTkCuts.AddRowColumn(9, "10");
  pixTkCuts.AddRowColumn(9, "TK_CANDPIXHITS_MIN");
  pixTkCuts.AddRowColumn(9, " >= " );
  pixTkCuts.AddRowColumn(9, auxTools.ToString(TK_CANDPIXHITS_MIN) );
  pixTkCuts.AddRowColumn(9, "-" );
  pixTkCuts.AddRowColumn(9, "Minimum number of candidate Pixel Hits" );
  pixTkCuts.AddRowColumn(9, "Default is -1" );


  pixTkCuts.AddRowColumn(10, "11");
  pixTkCuts.AddRowColumn(10, "TK_CANDPIXHITS_MAX");
  pixTkCuts.AddRowColumn(10, " <= " );
  pixTkCuts.AddRowColumn(10, auxTools.ToString(TK_CANDPIXHITS_MAX) );
  pixTkCuts.AddRowColumn(10, "-" );
  pixTkCuts.AddRowColumn(10, "Maximum number of candidate Pixel Hits" );
  pixTkCuts.AddRowColumn(10, "Default is +1e6" );
  
  pixTkCuts.Print();
  
  return;
}


//****************************************************************************
void PixelRefitting::PrintResolutions(void)
//****************************************************************************
{
   
  // Low pT
  double ptResolution_L_C  = GetHistoMeanInRange(h_resolutionVsEta_pt_L , 0.0, 0.8);
  double ptResolution_L_I  = GetHistoMeanInRange(h_resolutionVsEta_pt_L , 0.8, 1.6);
  double ptResolution_L_F  = GetHistoMeanInRange(h_resolutionVsEta_pt_L , 1.6, 2.5);

  double etaResolution_L_C = GetHistoMeanInRange(h_resolutionVsEta_eta_L, 0.0, 0.8);
  double etaResolution_L_I = GetHistoMeanInRange(h_resolutionVsEta_eta_L, 0.8, 1.6);
  double etaResolution_L_F = GetHistoMeanInRange(h_resolutionVsEta_eta_L, 1.6, 2.5);

  double phiResolution_L_C = GetHistoMeanInRange(h_resolutionVsEta_phi_L, 0.0, 0.8);
  double phiResolution_L_I = GetHistoMeanInRange(h_resolutionVsEta_phi_L, 0.8, 1.6);
  double phiResolution_L_F = GetHistoMeanInRange(h_resolutionVsEta_phi_L, 1.6, 2.5);

  double z0Resolution_L_C  = GetHistoMeanInRange(h_resolutionVsEta_z0_L , 0.0, 0.8);
  double z0Resolution_L_I  = GetHistoMeanInRange(h_resolutionVsEta_z0_L , 0.8, 1.6);
  double z0Resolution_L_F  = GetHistoMeanInRange(h_resolutionVsEta_z0_L , 1.6, 2.5);

  double d0Resolution_L_C  = GetHistoMeanInRange(h_resolutionVsEta_d0_L , 0.0, 0.8);
  double d0Resolution_L_I  = GetHistoMeanInRange(h_resolutionVsEta_d0_L , 0.8, 1.6);
  double d0Resolution_L_F  = GetHistoMeanInRange(h_resolutionVsEta_d0_L , 1.6, 2.5);

  // Middle pT
  double ptResolution_M_C  = GetHistoMeanInRange(h_resolutionVsEta_pt_M, 0.0, 0.8);
  double ptResolution_M_I  = GetHistoMeanInRange(h_resolutionVsEta_pt_M, 0.8, 1.6);
  double ptResolution_M_F  = GetHistoMeanInRange(h_resolutionVsEta_pt_M, 1.6, 2.5);

  double etaResolution_M_C = GetHistoMeanInRange(h_resolutionVsEta_eta_M, 0.0, 0.8);
  double etaResolution_M_I = GetHistoMeanInRange(h_resolutionVsEta_eta_M, 0.8, 1.6);
  double etaResolution_M_F = GetHistoMeanInRange(h_resolutionVsEta_eta_M, 1.6, 2.5);

  double phiResolution_M_C = GetHistoMeanInRange(h_resolutionVsEta_phi_M, 0.0, 0.8);
  double phiResolution_M_I = GetHistoMeanInRange(h_resolutionVsEta_phi_M, 0.8, 1.6);
  double phiResolution_M_F = GetHistoMeanInRange(h_resolutionVsEta_phi_M, 1.6, 2.5);

  double z0Resolution_M_C  = GetHistoMeanInRange(h_resolutionVsEta_z0_M, 0.0, 0.8);
  double z0Resolution_M_I  = GetHistoMeanInRange(h_resolutionVsEta_z0_M, 0.8, 1.6);
  double z0Resolution_M_F  = GetHistoMeanInRange(h_resolutionVsEta_z0_M, 1.6, 2.5);

  double d0Resolution_M_C  = GetHistoMeanInRange(h_resolutionVsEta_d0_M, 0.0, 0.8);
  double d0Resolution_M_I  = GetHistoMeanInRange(h_resolutionVsEta_d0_M, 0.8, 1.6);
  double d0Resolution_M_F  = GetHistoMeanInRange(h_resolutionVsEta_d0_M, 1.6, 2.5);

  // High pT
  double ptResolution_H_C  = GetHistoMeanInRange(h_resolutionVsEta_pt_H, 0.0, 0.8);
  double ptResolution_H_I  = GetHistoMeanInRange(h_resolutionVsEta_pt_H, 0.8, 1.6);
  double ptResolution_H_F  = GetHistoMeanInRange(h_resolutionVsEta_pt_H, 1.6, 2.5);

  double etaResolution_H_C = GetHistoMeanInRange(h_resolutionVsEta_eta_H, 0.0, 0.8);
  double etaResolution_H_I = GetHistoMeanInRange(h_resolutionVsEta_eta_H, 0.8, 1.6);
  double etaResolution_H_F = GetHistoMeanInRange(h_resolutionVsEta_eta_H, 1.6, 2.5);

  double phiResolution_H_C = GetHistoMeanInRange(h_resolutionVsEta_phi_H, 0.0, 0.8);
  double phiResolution_H_I = GetHistoMeanInRange(h_resolutionVsEta_phi_H, 0.8, 1.6);
  double phiResolution_H_F = GetHistoMeanInRange(h_resolutionVsEta_phi_H, 1.6, 2.5);

  double z0Resolution_H_C  = GetHistoMeanInRange(h_resolutionVsEta_z0_H, 0.0, 0.8);
  double z0Resolution_H_I  = GetHistoMeanInRange(h_resolutionVsEta_z0_H, 0.8, 1.6);
  double z0Resolution_H_F  = GetHistoMeanInRange(h_resolutionVsEta_z0_H, 1.6, 2.5);

  double d0Resolution_H_C  = GetHistoMeanInRange(h_resolutionVsEta_d0_H, 0.0, 0.8);
  double d0Resolution_H_I  = GetHistoMeanInRange(h_resolutionVsEta_d0_H, 0.8, 1.6);
  double d0Resolution_H_F  = GetHistoMeanInRange(h_resolutionVsEta_d0_H, 1.6, 2.5);

  // Inclusive
  double ptResolution_C  = GetHistoMeanInRange(h_resolutionVsEta_pt , 0.0, 0.8);
  double ptResolution_I  = GetHistoMeanInRange(h_resolutionVsEta_pt , 0.8, 1.6);
  double ptResolution_F  = GetHistoMeanInRange(h_resolutionVsEta_pt , 1.6, 2.5);

  double etaResolution_C = GetHistoMeanInRange(h_resolutionVsEta_eta, 0.0, 0.8);
  double etaResolution_I = GetHistoMeanInRange(h_resolutionVsEta_eta, 0.8, 1.6);
  double etaResolution_F = GetHistoMeanInRange(h_resolutionVsEta_eta, 1.6, 2.5);

  double phiResolution_C = GetHistoMeanInRange(h_resolutionVsEta_phi, 0.0, 0.8);
  double phiResolution_I = GetHistoMeanInRange(h_resolutionVsEta_phi, 0.8, 1.6);
  double phiResolution_F = GetHistoMeanInRange(h_resolutionVsEta_phi, 1.6, 2.5);

  double z0Resolution_C  = GetHistoMeanInRange(h_resolutionVsEta_z0 , 0.0, 0.8);
  double z0Resolution_I  = GetHistoMeanInRange(h_resolutionVsEta_z0 , 0.8, 1.6);
  double z0Resolution_F  = GetHistoMeanInRange(h_resolutionVsEta_z0 , 1.6, 2.5);

  double d0Resolution_C  = GetHistoMeanInRange(h_resolutionVsEta_d0 , 0.0, 0.8);
  double d0Resolution_I  = GetHistoMeanInRange(h_resolutionVsEta_d0 , 0.8, 1.6);
  double d0Resolution_F  = GetHistoMeanInRange(h_resolutionVsEta_d0 , 1.6, 2.5);

    
  // Create resolutions table
  
  Table resolutions_incl("Quantity | \\abs{\\eta} \\leq 0.8 | 0.8 < \\abs{\\eta} \\leq 1.6 | \\abs{\\eta} > 1.6 | Units", "LaTeX", "c R R R c");
  resolutions_incl.AddRowColumn(0, "$\\pT$");
  resolutions_incl.AddRowColumn(0, auxTools.ToString(1000 * ptResolution_H_C, 3) );
  resolutions_incl.AddRowColumn(0, auxTools.ToString(1000 * ptResolution_H_I, 3) );
  resolutions_incl.AddRowColumn(0, auxTools.ToString(1000 * ptResolution_H_F, 3) );
  resolutions_incl.AddRowColumn(0, "$\\MeVc{-1}$");
  
  resolutions_incl.AddRowColumn(1, "$\\eta$");  
  resolutions_incl.AddRowColumn(1, auxTools.ToString(1000 * etaResolution_H_C, 3) );
  resolutions_incl.AddRowColumn(1, auxTools.ToString(1000 * etaResolution_H_I, 3) );
  resolutions_incl.AddRowColumn(1, auxTools.ToString(1000 * etaResolution_H_F, 3) );
  resolutions_incl.AddRowColumn(1, "$10^{-3}$");
  
  resolutions_incl.AddRowColumn(2, "$\\phi$");
  resolutions_incl.AddRowColumn(2, auxTools.ToString(1000 * phiResolution_H_C, 3) );
  resolutions_incl.AddRowColumn(2, auxTools.ToString(1000 * phiResolution_H_I, 3) );
  resolutions_incl.AddRowColumn(2, auxTools.ToString(1000 * phiResolution_H_F, 3) );
  resolutions_incl.AddRowColumn(2, "$\\sMilliRads$");
  
  resolutions_incl.AddRowColumn(3, "$\\zPOCA{}$");  
  resolutions_incl.AddRowColumn(3, auxTools.ToString(10 * z0Resolution_H_C, 2) );
  resolutions_incl.AddRowColumn(3, auxTools.ToString(10 * z0Resolution_H_I, 2) );
  resolutions_incl.AddRowColumn(3, auxTools.ToString(10 * z0Resolution_H_F, 2) );
  resolutions_incl.AddRowColumn(3, "$\\sMilliMeter$");
  
  resolutions_incl.AddRowColumn(4, "$\\dZero{}$");
  resolutions_incl.AddRowColumn(4, auxTools.ToString(10000 * d0Resolution_H_C, 3) );
  resolutions_incl.AddRowColumn(4, auxTools.ToString(10000 * d0Resolution_H_I, 3) );
  resolutions_incl.AddRowColumn(4, auxTools.ToString(10000 * d0Resolution_H_F, 3) );
  resolutions_incl.AddRowColumn(4, "$\\sMicroMeter$");


  
  Table resolutions("Quantity | \\abs{\\eta} \\leq 0.8 | 0.8 < \\abs{\\eta} \\leq 1.6 | \\abs{\\eta} > 1.6 | Units", "LaTeX", "c R R R c");
  resolutions.AddRowColumn(0, "\\multicolumn{5}{c}{$\\pT \\leq " + auxTools.ToString(_pt_L) + "\\, \\GeVc{-1}$}");
  resolutions.AddRowColumn(1, "$\\pT$");
  resolutions.AddRowColumn(1, auxTools.ToString(1000 * ptResolution_L_C, 3) );
  resolutions.AddRowColumn(1, auxTools.ToString(1000 * ptResolution_L_I, 3) );
  resolutions.AddRowColumn(1, auxTools.ToString(1000 * ptResolution_L_F, 3) );
  resolutions.AddRowColumn(1, "$\\MeVc{-1}$");
  
  resolutions.AddRowColumn(2, "$\\eta$");
  resolutions.AddRowColumn(2, auxTools.ToString(1000 * etaResolution_L_C, 3) );
  resolutions.AddRowColumn(2, auxTools.ToString(1000 * etaResolution_L_I, 3) );
  resolutions.AddRowColumn(2, auxTools.ToString(1000 * etaResolution_L_F, 3) );
  resolutions.AddRowColumn(2, "$10^{-3}$");
 
  resolutions.AddRowColumn(3, "$\\phi$");
  resolutions.AddRowColumn(3, auxTools.ToString(1000 * phiResolution_L_C, 3) );
  resolutions.AddRowColumn(3, auxTools.ToString(1000 * phiResolution_L_I, 3) );
  resolutions.AddRowColumn(3, auxTools.ToString(1000 * phiResolution_L_F, 3) );
  resolutions.AddRowColumn(3, "$\\sMilliRads$");
  
  resolutions.AddRowColumn(4, "$\\zPOCA{}$");
  resolutions.AddRowColumn(4, auxTools.ToString(10 * z0Resolution_L_C, 2) );
  resolutions.AddRowColumn(4, auxTools.ToString(10 * z0Resolution_L_I, 2) );
  resolutions.AddRowColumn(4, auxTools.ToString(10 * z0Resolution_L_F, 2) );
  resolutions.AddRowColumn(4, "$\\sMilliMeter$");
 
  resolutions.AddRowColumn(5, "$\\dZero{}$");  
  resolutions.AddRowColumn(5, auxTools.ToString(10000 * d0Resolution_L_C, 3) );
  resolutions.AddRowColumn(5, auxTools.ToString(10000 * d0Resolution_L_I, 3) );
  resolutions.AddRowColumn(5, auxTools.ToString(10000 * d0Resolution_L_F, 3) );
  resolutions.AddRowColumn(5, "$\\sMicroMeter$");
  
  resolutions.AddRowColumn(6, "\\multicolumn{5}{c}{$" + auxTools.ToString(_pt_L) + "\\, \\GeVc{-1} < \\pT \\leq " + auxTools.ToString(_pt_H) + "\\, \\GeVc{-1}$}");
  resolutions.AddRowColumn(7, "$\\pT$");
  resolutions.AddRowColumn(7, auxTools.ToString(1000 * ptResolution_M_C, 3) );
  resolutions.AddRowColumn(7, auxTools.ToString(1000 * ptResolution_M_I, 3) );
  resolutions.AddRowColumn(7, auxTools.ToString(1000 * ptResolution_M_F, 3) );
  resolutions.AddRowColumn(7, "$\\MeVc{-1}$");
  
  resolutions.AddRowColumn(8, "$\\eta$");  
  resolutions.AddRowColumn(8, auxTools.ToString(1000 * etaResolution_M_C, 3) );
  resolutions.AddRowColumn(8, auxTools.ToString(1000 * etaResolution_M_I, 3) );
  resolutions.AddRowColumn(8, auxTools.ToString(1000 * etaResolution_M_F, 3) );
  resolutions.AddRowColumn(8, "$10^{-3}$");
 
  resolutions.AddRowColumn(9, "$\\phi$");  
  resolutions.AddRowColumn(9, auxTools.ToString(1000 * phiResolution_M_C, 2) );
  resolutions.AddRowColumn(9, auxTools.ToString(1000 * phiResolution_M_I, 2) );
  resolutions.AddRowColumn(9, auxTools.ToString(1000 * phiResolution_M_F, 2) );
  resolutions.AddRowColumn(9, "$\\sMilliRads$");
  
  resolutions.AddRowColumn(10, "$\\zPOCA{}$");  
  resolutions.AddRowColumn(10, auxTools.ToString(10 * z0Resolution_M_C, 2) );
  resolutions.AddRowColumn(10, auxTools.ToString(10 * z0Resolution_M_I, 2) );
  resolutions.AddRowColumn(10, auxTools.ToString(10 * z0Resolution_M_F, 2) );
  resolutions.AddRowColumn(10, "$\\sMilliMeter$");
  
  resolutions.AddRowColumn(11, "$\\dZero{}$");  
  resolutions.AddRowColumn(11, auxTools.ToString(10000 * d0Resolution_M_C, 3) );
  resolutions.AddRowColumn(11, auxTools.ToString(10000 * d0Resolution_M_I, 3) );
  resolutions.AddRowColumn(11, auxTools.ToString(10000 * d0Resolution_M_F, 3) );
  resolutions.AddRowColumn(11, "$\\sMicroMeter$");
  
  resolutions.AddRowColumn(12, "\\multicolumn{5}{c}{$\\pT > " + auxTools.ToString(_pt_H) + "\\, \\GeVc{-1}$}");
  resolutions.AddRowColumn(13, "$\\pT$");
  resolutions.AddRowColumn(13, auxTools.ToString(1000 * ptResolution_H_C, 3) );
  resolutions.AddRowColumn(13, auxTools.ToString(1000 * ptResolution_H_I, 3) );
  resolutions.AddRowColumn(13, auxTools.ToString(1000 * ptResolution_H_F, 3) );
  resolutions.AddRowColumn(13, "$\\MeVc{-1}$");
  
  resolutions.AddRowColumn(14, "$\\eta$");  
  resolutions.AddRowColumn(14, auxTools.ToString(1000 * etaResolution_H_C, 3) );
  resolutions.AddRowColumn(14, auxTools.ToString(1000 * etaResolution_H_I, 3) );
  resolutions.AddRowColumn(14, auxTools.ToString(1000 * etaResolution_H_F, 3) );
  resolutions.AddRowColumn(14, "$10^{-3}$");
  
  resolutions.AddRowColumn(15, "$\\phi$");
  resolutions.AddRowColumn(15, auxTools.ToString(1000 * phiResolution_H_C, 3) );
  resolutions.AddRowColumn(15, auxTools.ToString(1000 * phiResolution_H_I, 3) );
  resolutions.AddRowColumn(15, auxTools.ToString(1000 * phiResolution_H_F, 3) );
  resolutions.AddRowColumn(15, "$\\sMilliRads$");
  
  resolutions.AddRowColumn(16, "$\\zPOCA{}$");  
  resolutions.AddRowColumn(16, auxTools.ToString(10 * z0Resolution_H_C, 2) );
  resolutions.AddRowColumn(16, auxTools.ToString(10 * z0Resolution_H_I, 2) );
  resolutions.AddRowColumn(16, auxTools.ToString(10 * z0Resolution_H_F, 2) );
  resolutions.AddRowColumn(16, "$\\sMilliMeter$");
  
  resolutions.AddRowColumn(17, "$\\dZero{}$");
  resolutions.AddRowColumn(17, auxTools.ToString(10000 * d0Resolution_H_C, 3) );
  resolutions.AddRowColumn(17, auxTools.ToString(10000 * d0Resolution_H_I, 3) );
  resolutions.AddRowColumn(17, auxTools.ToString(10000 * d0Resolution_H_F, 3) );
  resolutions.AddRowColumn(17, "$\\sMicroMeter$");

  if(bPrintResolutions_Incl) resolutions_incl.Print();
  if(bPrintResolutions) resolutions.Print();

  
  return;
}



//****************************************************************************
double PixelRefitting::GetHistoMeanInRange(TH1D* histo, 
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
bool PixelRefitting::IsWithinTrackerAcceptance(double pt, 
					      double eta, 
					      double z0, 
					      double dxy)
//****************************************************************************
{

  bool bWithin_pt  = ( fabs( pt  ) >= TP_PT_MIN );
  bool bWithin_eta = ( fabs( eta ) <= TP_ETA_MAX );
  bool bWithin_z0  = ( fabs( z0  ) <= TP_Z0_MAX  );
  bool bWithin_dxy = ( fabs( dxy ) <= TP_DXY_MAX );  // Only TPs coming from near the IP (Louise) - Significantly affects d0 Resolution.

  return bWithin_pt * bWithin_eta * bWithin_z0 * bWithin_dxy;
}



//****************************************************************************
bool PixelRefitting::IsWithinEtaRegion(string etaRegion, 
				      double eta)
//****************************************************************************
{

  bool bWithinEtaRegion = false;
  if ( etaRegion.compare("Central") == 0 )           bWithinEtaRegion = (fabs(eta) <= _eta_C);
  else if ( etaRegion.compare("Intermediate") == 0 ) bWithinEtaRegion = (fabs(eta) <= _eta_F && fabs(eta) > _eta_C);
  else if ( etaRegion.compare("Forward") == 0 )      bWithinEtaRegion = (fabs(eta) > _eta_F);
  else{
    cout << "E R R O R ! PixelRefitting::IsWithinEtaRegion(...) - Invalid eta region type \"" << etaRegion << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinEtaRegion;  
}


//****************************************************************************
bool PixelRefitting::IsWithinPtRange(string ptRange, 
				    double pt)
//****************************************************************************
{

  bool bWithinPtRange = false;
  if ( ptRange.compare("Low") == 0 )         bWithinPtRange = (pt <= _pt_L);
  else if ( ptRange.compare("Middle") == 0 ) bWithinPtRange = (pt <= _pt_H && pt > _pt_L);
  else if ( ptRange.compare("High") == 0 )   bWithinPtRange = (pt > _pt_H);
  else{
    cout << "E R R O R ! PixelRefitting::IsWithinPtRange(...) - Invalid pt range type \"" << ptRange << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinPtRange;  
}


//****************************************************************************
void PixelRefitting::FillPixelHitHistos(TTPixelTrack pixFitTk,
				       bool bPrintPixInfo)
//****************************************************************************
{
  
  if (pixFitTk.getNhit() < 3) return;

  // Create Pixel-Hit Table
  Table pixInfo("Hits | X | Y | Z | R | Phi | Type | Pattern", "Text");
  
  // Get variables
  double tk_Eta = pixFitTk.getMomentum().Eta();
  vector<TVector3> pixHits = pixFitTk.getPixelHits();
  vector< int >  pixHits_Type;
  int pixHits_Pattern = pixFitTk.getPixelHitsPattern();
  int tk_PixHits      = pixFitTk.getNhit();    
  for (int i = 0; i < (int) pixHits.size(); i++) pixHits_Type.push_back( pixFitTk.getPixelHitType( pixHits.at(i) ) );

 
  // Fill Histograms
  h_pixTk_pixHits_N->Fill(tk_PixHits);
 
  // For-loop: Pixel Hits (best-fit)  
  for(int i = 0; i < tk_PixHits; i++){

    pixInfo.AddRowColumn(i, auxTools.ToString( tk_PixHits) );
    pixInfo.AddRowColumn(i, auxTools.ToString( pixHits.at(i).X() ) );
    pixInfo.AddRowColumn(i, auxTools.ToString( pixHits.at(i).Y() ) );
    pixInfo.AddRowColumn(i, auxTools.ToString( pixHits.at(i).Z() ) );
    pixInfo.AddRowColumn(i, auxTools.ToString( pixHits.at(i).Perp() ) );
    pixInfo.AddRowColumn(i, auxTools.ToString( pixHits.at(i).Phi() ) );
    pixInfo.AddRowColumn(i, auxTools.ToString( pixHits_Type.at(i)) );
    pixInfo.AddRowColumn(i, auxTools.ToString( pixHits_Pattern ) );

    // Fill Histos
    h_pixTk_pixHits_Rho   ->Fill( pixHits.at(i).Perp() );
    h_pixTk_pixHits_Z     ->Fill( pixHits.at(i).Z()    );
    h_pixTk_pixHits_Type  ->Fill( pixHits_Type.at(i)   );
    h_pixTk_pixHits_Pattern->Fill( pixHits_Pattern );
    h_pixTk_pixHits_PatternVsEta->Fill( pixHits_Pattern, fabs(tk_Eta) );
    
    h_pixTk_pixHits_EtaVsN->Fill( tk_Eta, tk_PixHits);
    h_pixTk_pixHits_ZVsRho->Fill( pixHits.at(i).Z(), pixHits.at(i).Perp() );   
    FillEtaRegionsHistos2D(tk_Eta, pixHits.at(i).Z(), pixHits.at(i).Perp(), h_pixTk_pixHits_ZVsRho_C, h_pixTk_pixHits_ZVsRho_I, h_pixTk_pixHits_ZVsRho_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_pixHits_ZVsRho_EtaGE2->Fill( pixHits.at(i).Z(), pixHits.at(i).Perp() );

    h_pixTk_pixHits_XVsY->Fill( pixHits.at(i).X(), pixHits.at(i).Y() );
    FillEtaRegionsHistos2D(tk_Eta, pixHits.at(i).X(), pixHits.at(i).Y(), h_pixTk_pixHits_XVsY_C, h_pixTk_pixHits_XVsY_I, h_pixTk_pixHits_XVsY_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_pixHits_XVsY_EtaGE2->Fill( pixHits.at(i).X(), pixHits.at(i).Y() );
    
  }// for(int i = 0; i < tk_PixHit; i++){
  
  if (bPrintPixInfo) pixInfo.Print();
  return;
}


//****************************************************************************
void PixelRefitting::FillCandPixelHitHistos(TTPixelTrack pixFitTk,
					   bool bPrintPixInfo)
//****************************************************************************
{
  if (pixFitTk.getNhit() < 3) return;

  // Create Pixel-Hit Table
  Table pixInfo("Candidate Hits | X | Y | Z | R | Phi | Type | Pattern", "Text");
  
  // Get variables
  double tk_Eta = pixFitTk.getMomentum().Eta();
  vector<TVector3> candPixHits = pixFitTk.getCandidatePixelHits();
  vector< int >  candPixHits_Type;
  for (int i = 0; i < (int) candPixHits.size(); i++) candPixHits_Type.push_back( pixFitTk.getPixelHitType( candPixHits.at(i) ) );
  int tk_PixHits = pixFitTk.getNcandidatehit();
  
  // Fill Histograms 
  h_pixTk_candPixHits_N->Fill(tk_PixHits);

  // For-loop: Pixel Hits (candidates)
  for(int j = 0; j < tk_PixHits; j++){

    pixInfo.AddRowColumn(j, auxTools.ToString( tk_PixHits) );
    pixInfo.AddRowColumn(j, auxTools.ToString( candPixHits.at(j).X() )    );
    pixInfo.AddRowColumn(j, auxTools.ToString( candPixHits.at(j).Y() )    );
    pixInfo.AddRowColumn(j, auxTools.ToString( candPixHits.at(j).Z() )    );
    pixInfo.AddRowColumn(j, auxTools.ToString( candPixHits.at(j).Perp() ) );
    pixInfo.AddRowColumn(j, auxTools.ToString( candPixHits.at(j).Phi() )  );
    pixInfo.AddRowColumn(j, auxTools.ToString( candPixHits_Type.at(j) )   );

    // Fill Histos
    h_pixTk_candPixHits_Rho   ->Fill( candPixHits.at(j).Perp() );
    h_pixTk_candPixHits_Z     ->Fill( candPixHits.at(j).Z() );
    h_pixTk_candPixHits_Type  ->Fill( candPixHits_Type.at(j) );
    h_pixTk_candPixHits_EtaVsN->Fill( tk_Eta, tk_PixHits);
    h_pixTk_candPixHits_ZVsRho->Fill( candPixHits.at(j).Z(), candPixHits.at(j).Perp() );
    FillEtaRegionsHistos2D(tk_Eta, candPixHits.at(j).Z(), candPixHits.at(j).Perp(), h_pixTk_candPixHits_ZVsRho_C, h_pixTk_candPixHits_ZVsRho_I, h_pixTk_candPixHits_ZVsRho_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_candPixHits_ZVsRho_EtaGE2->Fill( candPixHits.at(j).Z(), candPixHits.at(j).Perp() );
    
    h_pixTk_candPixHits_XVsY->Fill( candPixHits.at(j).X(), candPixHits.at(j).Y() );
    FillEtaRegionsHistos2D(tk_Eta, candPixHits.at(j).X(), candPixHits.at(j).Y(), h_pixTk_candPixHits_XVsY_C, h_pixTk_candPixHits_XVsY_I, h_pixTk_candPixHits_XVsY_F);    
    if (fabs(tk_Eta) >= 2.0) h_pixTk_candPixHits_XVsY_EtaGE2->Fill( candPixHits.at(j).X(), candPixHits.at(j).Y() );
      
  }// for(int j = 0; j < pixTk_CandPixHits; j++){

  if (bPrintPixInfo) pixInfo.Print();
  return;
}


//****************************************************************************
void PixelRefitting::FillSharedPixelHitHistos(int pixTk_Index,
					     int tp_Index,
					     vector<TVector3> pixHits_XYZ,
					     vector<int> pixHits_TTPixelTrackIndex,
					     bool bPrintPixInfo)
//****************************************************************************
{
  if (pixTk_Index < 0) return;

  // Create Pixel-Hit Table
  Table pixInfo("Index | Hit Index | Type | Hit R | Hit Z | Hit Phi | TTPixelTrack | TTTrack", "Text");

  // Get variables
  double tk_Eta = L1PixTks_Eta->at(pixTk_Index);
  const int TTTrackIndex = L1PixTks_TTTrackIndex->at(pixTk_Index);
  vector<int>    pixHits_PixTkRefIndex;
  vector<double> pixHits_X;
  vector<double> pixHits_Y;
  vector<double> pixHits_Z;
  vector<double> pixHits_R;
  vector<double> pixHits_Phi;
  vector< int >  pixHits_Type;
  
  // Get shared pixel hits!
  s->GetPixTrackSharedHits(pixTk_Index, pixHits_PixTkRefIndex, pixHits_X, pixHits_Y, pixHits_Z, pixHits_R, pixHits_Phi, pixHits_Type, pixHits_XYZ, pixHits_TTPixelTrackIndex);

  // Fill Histos
  int tk_PixHits = (int) pixHits_X.size();
  h_pixTk_sharedPixHits_N->Fill(tk_PixHits);

  // For-loop: Pixel Hits (candidates)
  for(int j = 0; j < tk_PixHits; j++){

    pixInfo.AddRowColumn(j, auxTools.ToString(pixTk_Index) );
    pixInfo.AddRowColumn(j, auxTools.ToString(j) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_Type.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_R.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_Z.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_Phi.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_PixTkRefIndex.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(TTTrackIndex) );
    
    // Fill Histos
    h_pixTk_sharedPixHits_Rho   ->Fill( pixHits_R.at(j) );
    h_pixTk_sharedPixHits_Z     ->Fill( pixHits_Z.at(j) );
    h_pixTk_sharedPixHits_Type  ->Fill( pixHits_Type.at(j) );
    h_pixTk_sharedPixHits_EtaVsN->Fill( tk_Eta, tk_PixHits);
    h_pixTk_sharedPixHits_ZVsRho->Fill( pixHits_Z.at(j), pixHits_R.at(j) );
    FillEtaRegionsHistos2D(tk_Eta, pixHits_Z.at(j), pixHits_R.at(j), h_pixTk_sharedPixHits_ZVsRho_C, h_pixTk_sharedPixHits_ZVsRho_I, h_pixTk_sharedPixHits_ZVsRho_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_sharedPixHits_ZVsRho_EtaGE2->Fill( pixHits_Z.at(j), pixHits_R.at(j) );

    h_pixTk_sharedPixHits_XVsY->Fill( pixHits_X.at(j), pixHits_Y.at(j) );
    FillEtaRegionsHistos2D(tk_Eta, pixHits_X.at(j), pixHits_Y.at(j), h_pixTk_sharedPixHits_XVsY_C, h_pixTk_sharedPixHits_XVsY_I, h_pixTk_sharedPixHits_XVsY_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_sharedPixHits_XVsY_EtaGE2->Fill( pixHits_X.at(j), pixHits_Y.at(j) );
          
  }// for(int j = 0; j < pixTk_CandPixHits; j++){

  if (bPrintPixInfo){
    if (pixHits_PixTkRefIndex.size()>0){

      std::cout << "Found Shared Hits!" << std::endl;
      // pixInfo.Print();
      tp->PrintTPProperties(tp_Index);
      s->PrintPixTrackProperties(pixTk_Index, true);
      s->PrintPixTrackProperties(pixHits_PixTkRefIndex.at(0), true);

    }// if (pixHits_PixTkRefIndex.size()>0){

  }// if (bPrintPixInfo){

  return;
}


//****************************************************************************
void PixelRefitting::FillResidualHistos(int tp_Index,
				       TTPixelTrack pixFitTk)
//****************************************************************************
{

  // Get TP variables
  double tp_Pt  = TP_Pt->at(tp_Index);
  double tp_Eta = TP_Eta->at(tp_Index);
  double tp_Phi = TP_Phi->at(tp_Index);
  double tp_z0  = TP_POCAz->at(tp_Index);
  double tp_d0  = tp->GetD0(tp_Index);

  // Get TTPixelTrack variables
  double tk_Pt         = pixFitTk.getMomentum().Perp();
  double tk_Eta        = pixFitTk.getMomentum().Eta();
  double tk_Phi        = pixFitTk.getMomentum().Phi();
  double tk_z0         = pixFitTk.getZ0();
  double tk_d0         = pixFitTk.getD0();
  
  // Fill Resolutions (Total)
  h_residual_pt   ->Fill( tk_Pt  - tp_Pt);
  h_residual_ptRel->Fill( (tk_Pt - tp_Pt)/tp_Pt);
  h_residual_eta  ->Fill( tk_Eta - tp_Eta);
  h_residual_phi  ->Fill( tk_Phi - tp_Phi);
  h_residual_z0   ->Fill( tk_z0  - tp_z0);
  h_residual_d0   ->Fill( tk_d0  - tp_d0);
  
  // Fill Resolutions (Eta Regions)
  if ( IsWithinEtaRegion("Central", tk_Eta) ){
    h_residual_pt_C   ->Fill( tk_Pt  - tp_Pt);
    h_residual_ptRel_C->Fill( (tk_Pt - tp_Pt)/tp_Pt);
    h_residual_eta_C  ->Fill( tk_Eta - tp_Eta);
    h_residual_phi_C  ->Fill( tk_Phi - tp_Phi);
    h_residual_z0_C   ->Fill( tk_z0 - tp_z0);
    h_residual_d0_C   ->Fill( tk_d0 - tp_d0);
  }
  else if ( IsWithinEtaRegion("Intermediate", tk_Eta) ){
    h_residual_pt_I   ->Fill( tk_Pt  - tp_Pt);
    h_residual_ptRel_I->Fill( (tk_Pt - tp_Pt)/tp_Pt);
    h_residual_eta_I  ->Fill( tk_Eta - tp_Eta);
    h_residual_phi_I  ->Fill( tk_Phi - tp_Phi);
    h_residual_z0_I   ->Fill( tk_z0 - tp_z0);
    h_residual_d0_I   ->Fill( tk_d0 - tp_d0);
  }
  else if ( IsWithinEtaRegion("Forward", tk_Eta) ){
    h_residual_pt_F   ->Fill( tk_Pt  - tp_Pt);
    h_residual_ptRel_F->Fill( (tk_Pt - tp_Pt)/tp_Pt);
    h_residual_eta_F  ->Fill( tk_Eta - tp_Eta);
    h_residual_phi_F  ->Fill( tk_Phi - tp_Phi);
    h_residual_z0_F   ->Fill( tk_z0 - tp_z0);
    h_residual_d0_F   ->Fill( tk_d0 - tp_d0);
  }
  else{
    cout << "E R R O R ! PixelRefitting::FillResidualHistos(...) - Unexpected eta of \"" << tk_Eta << "\". EXIT" << endl;
    exit(1);
  }

  // For-loop: Pt Bins
  for (int j=0; j < PT_BINS; j++) {

    // Fill residual vs. pt histograms in pT-steps of 5 GeV/c
    double ptLow  = (double)j * PT_BINWIDTH;
    double ptHigh = ptLow + 5.0;
    bool bPtIsWithinPtBin = ( (tp_Pt >= ptLow) && (tp_Pt < ptHigh) );
    if (!bPtIsWithinPtBin) continue;

    // Fill Residuals
    h_residualVsPt_pt[j]   ->Fill( tk_Pt  - tp_Pt);
    h_residualVsPt_ptRel[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
    h_residualVsPt_eta[j]  ->Fill( tk_Eta - tp_Eta);
    h_residualVsPt_phi[j]  ->Fill( tk_Phi - tp_Phi);
    h_residualVsPt_z0[j]   ->Fill( tk_z0  - tp_z0);
    h_residualVsPt_d0[j]   ->Fill( tk_d0  - tp_d0);

    // Fill Residuals (Eta Regions)    
    if ( IsWithinEtaRegion("Central", tk_Eta) ){
      h_residualVsPt_pt_C[j]   ->Fill( tk_Pt  - tp_Pt);
      h_residualVsPt_ptRel_C[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsPt_eta_C[j]  ->Fill( tk_Eta - tp_Eta);
      h_residualVsPt_phi_C[j]  ->Fill( tk_Phi - tp_Phi);
      h_residualVsPt_z0_C[j]   ->Fill( tk_z0  - tp_z0);
      h_residualVsPt_d0_C[j]   ->Fill( tk_d0  - tp_d0);
    }
    else if ( IsWithinEtaRegion("Intermediate", tk_Eta) ){
      h_residualVsPt_pt_I[j]   ->Fill( tk_Pt  - tp_Pt);
      h_residualVsPt_ptRel_I[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsPt_eta_I[j]  ->Fill( tk_Eta - tp_Eta);
      h_residualVsPt_phi_I[j]  ->Fill( tk_Phi - tp_Phi);
      h_residualVsPt_z0_I[j]   ->Fill( tk_z0  - tp_z0);
      h_residualVsPt_d0_I[j]   ->Fill( tk_d0  - tp_d0);
    }	
    else if ( IsWithinEtaRegion("Forward", tk_Eta) ){
      h_residualVsPt_pt_F[j]   ->Fill( tk_Pt  - tp_Pt);
      h_residualVsPt_ptRel_F[j]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsPt_eta_F[j]  ->Fill( tk_Eta - tp_Eta);
      h_residualVsPt_phi_F[j]  ->Fill( tk_Phi - tp_Phi);
      h_residualVsPt_z0_F[j]   ->Fill( tk_z0  - tp_z0);
      h_residualVsPt_d0_F[j]   ->Fill( tk_d0  - tp_d0);
    }
    else{
      cout << "E R R O R ! PixelRefitting::FillResidualHistos(...) - Unexpected Eta value of \"" << tk_Eta << "\". EXIT" << endl;
      exit(1);
    }
  } // for (int j=0; j < PT_BINS; j++) {
  


  // For-loop: Eta Bin
  for (int k=0; k < ETA_BINS; k++) {

    // Fill residual vs. eta histogram in eta-steps of DeltaEta = 0.1
    double etaLow  = (double)k * ETA_BINWIDTH;
    double etaHigh = etaLow + ETA_BINWIDTH;
    bool bEtaIsWithinEtaBin = ( (fabs(tp_Eta) >= etaLow) && (fabs(tp_Eta) < etaHigh) );
    if (!bEtaIsWithinEtaBin) continue;

    // Fill Residuals
    h_residualVsEta_pt[k]   ->Fill( tk_Pt  - tp_Pt);
    h_residualVsEta_ptRel[k]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
    h_residualVsEta_eta[k]  ->Fill( tk_Eta - tp_Eta);
    h_residualVsEta_phi[k]  ->Fill( tk_Phi - tp_Phi);
    h_residualVsEta_z0[k]   ->Fill( tk_z0  - tp_z0);
    h_residualVsEta_d0[k]   ->Fill( tk_d0  - tp_d0);
    
    // Fill Residuals (Pt Regions)    
    if ( IsWithinPtRange("Low", tk_Pt) ){
      h_residualVsEta_pt_L[k]   ->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsEta_ptRel_L[k]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsEta_eta_L[k]  ->Fill( tk_Eta - tp_Eta);
      h_residualVsEta_phi_L[k]  ->Fill( tk_Phi - tp_Phi);
      h_residualVsEta_z0_L[k]   ->Fill( tk_z0  - tp_z0);
      h_residualVsEta_d0_L[k]   ->Fill( tk_d0  - tp_d0);
    }
    else if ( IsWithinPtRange("Middle", tk_Pt) ){
      h_residualVsEta_pt_M[k]   ->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsEta_ptRel_M[k]->Fill( (tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsEta_eta_M[k]  ->Fill( tk_Eta - tp_Eta);
      h_residualVsEta_phi_M[k]  ->Fill( tk_Phi - tp_Phi);
      h_residualVsEta_z0_M[k]   ->Fill( tk_z0  - tp_z0);
      h_residualVsEta_d0_M[k]   ->Fill( tk_d0  - tp_d0);
    }
    else if ( IsWithinPtRange("High", tk_Pt) ){
      h_residualVsEta_pt_H[k]   ->Fill((tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsEta_ptRel_H[k]->Fill((tk_Pt - tp_Pt)/tp_Pt);
      h_residualVsEta_eta_H[k]  ->Fill(tk_Eta - tp_Eta);
      h_residualVsEta_phi_H[k]  ->Fill(tk_Phi - tp_Phi);
      h_residualVsEta_z0_H[k]   ->Fill(tk_z0 - tp_z0);
      h_residualVsEta_d0_H[k]   ->Fill(tk_d0 - tp_d0);
    }
    else{
      cout << "E R R O R ! PixelRefitting::Loop(...) - Unexpected Track-Pt of \"" << tk_Pt << "\". EXIT" << endl;
      exit(1);
    }	  

  } // for (int k=0; k < ETA_BINS; k++) {

  return;
}


//****************************************************************************
void PixelRefitting::FillEtaRegionsHistos(double eta,
					 double variable,
					 TH1D *hCentral,
					 TH1D *hIntermediate,
					 TH1D *hForward)
//****************************************************************************
{

  if (IsWithinEtaRegion("Central", eta) )            hCentral->Fill(variable);
  else if ( IsWithinEtaRegion("Intermediate", eta) ) hIntermediate->Fill(variable);
  else if ( IsWithinEtaRegion("Forward", eta) )      hForward->Fill(variable);
  else{
    cout << "E R R O R ! PixelRefitting::FillEtaRegionsHistos(...) - Unexpected Eta value of \"" << eta << "\". EXIT" << endl;
    exit(1);
  }
  
  return;
}


//****************************************************************************
void PixelRefitting::FillEtaRegionsHistos2D(double eta,
					   double variableX,
					   double variableY,
					   TH2D *hCentral,
					   TH2D *hIntermediate,
					   TH2D *hForward)
//****************************************************************************
{

  if (IsWithinEtaRegion("Central", eta) )            hCentral->Fill(variableX, variableY);
  else if ( IsWithinEtaRegion("Intermediate", eta) ) hIntermediate->Fill(variableX, variableY);
  else if ( IsWithinEtaRegion("Forward", eta) )      hForward->Fill(variableX, variableY);
  else{
    cout << "E R R O R ! PixelRefitting::FillEtaRegionsHistos(...) - Unexpected Eta value of \"" << eta << "\". EXIT" << endl;
    exit(1);
  }
  
  return;
}


//****************************************************************************
void PixelRefitting::FillPtRegionsHistos(double pt,
					double variable,
					TH1D *hLow,
					TH1D *hMiddle,
					TH1D *hHigh)
//****************************************************************************
{
  
  if ( IsWithinPtRange("Low", pt) )         hLow->Fill(variable);
  else if ( IsWithinPtRange("Middle", pt) ) hMiddle->Fill(variable);
  else if ( IsWithinPtRange("High", pt) )   hHigh->Fill(variable);
  else{
    cout << "E R R O R ! PixelRefitting::FillPtRegionsHistos(...) - Unexpected pT value of \"" << pt << "\". EXIT" << endl;
    exit(1);
  }
  
  return;  
}


//****************************************************************************
void PixelRefitting::FillTPMatchedPixTrackHistos(TTPixelTrack pixFitTk)
//****************************************************************************
{

  // Get TTPixelTrack variables
  double tk_Pt         = pixFitTk.getMomentum().Perp();
  double tk_Eta        = pixFitTk.getMomentum().Eta();
  double tk_ChiSq      = pixFitTk.getChi2();
  double tk_RedChiSq   = pixFitTk.getChi2Red();

  // Overflow bins for chi-square-related distributions
  if (tk_ChiSq > _chiSq_overflow)  tk_ChiSq = _chiSq_overflow - 0.01;
  if (tk_RedChiSq > _redChiSq_overflow) tk_RedChiSq = _redChiSq_overflow - 0.01;
  
  // Generic Histograms
  h_match_trk_chi2    ->Fill(tk_ChiSq   );
  h_match_trk_chi2_dof->Fill(tk_RedChiSq);    

  if ( IsWithinEtaRegion("Central", tk_Eta) ){
    
    if ( IsWithinPtRange("Low", tk_Pt) ){
      h_match_trk_chi2_C_L    ->Fill(tk_ChiSq   );
      h_match_trk_chi2_dof_C_L->Fill(tk_RedChiSq);
      h_match_trk_chi2_L    ->Fill(tk_ChiSq   );
      h_match_trk_chi2_dof_L->Fill(tk_RedChiSq);
    }
    else if ( IsWithinPtRange("Middle", tk_Pt) ) {
      h_match_trk_chi2_C_M    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_C_M->Fill(tk_RedChiSq);
      h_match_trk_chi2_M    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_M->Fill(tk_RedChiSq);
    }
    else if ( IsWithinPtRange("High", tk_Pt) ){
      h_match_trk_chi2_C_H->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_C_H->Fill(tk_RedChiSq);
      h_match_trk_chi2_H    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_H->Fill(tk_RedChiSq);
    }
    else{
      cout << "E R R O R ! PixelRefitting::Loop(...) - Unexpected tk pt value of \"" << tk_Pt << "\". EXIT" << endl;
      exit(1);
    }
  }// Central eta
  else if ( IsWithinEtaRegion("Intermediate", tk_Eta) ){
    
    if ( IsWithinPtRange("Low", tk_Pt) ){
      h_match_trk_chi2_I_L    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_I_L->Fill(tk_RedChiSq);
      h_match_trk_chi2_L    ->Fill(tk_ChiSq   );
      h_match_trk_chi2_dof_L->Fill(tk_RedChiSq);	 
    }
    else if ( IsWithinPtRange("Middle", tk_Pt) ) {
      h_match_trk_chi2_I_M->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_I_M->Fill(tk_RedChiSq);
      h_match_trk_chi2_M    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_M->Fill(tk_RedChiSq);
    }
    else if ( IsWithinPtRange("High", tk_Pt) ){
      h_match_trk_chi2_I_H->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_I_H->Fill(tk_RedChiSq);
      h_match_trk_chi2_H    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_H->Fill(tk_RedChiSq);
	}
    else{
      cout << "E R R O R ! PixelRefitting::Loop(...) - Unexpected tk pt value of \"" << tk_Pt << "\". EXIT" << endl;
      exit(1);
    }
  }// Intermediate Eta Region
  else if ( IsWithinEtaRegion("Forward", tk_Eta) ) {
    if ( IsWithinPtRange("Low", tk_Pt) ){
      h_match_trk_chi2_F_L    ->Fill(tk_ChiSq   );
      h_match_trk_chi2_dof_F_L->Fill(tk_RedChiSq);
      h_match_trk_chi2_L    ->Fill(tk_ChiSq   );
      h_match_trk_chi2_dof_L->Fill(tk_RedChiSq);
    }
    else if ( IsWithinPtRange("Middle", tk_Pt) ) {
      h_match_trk_chi2_F_M    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_F_M->Fill(tk_RedChiSq);
      h_match_trk_chi2_M    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_M->Fill(tk_RedChiSq);
    }
    else if ( IsWithinPtRange("High", tk_Pt) ){
      h_match_trk_chi2_F_H    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_F_H->Fill(tk_RedChiSq);
      h_match_trk_chi2_H    ->Fill(tk_ChiSq);
      h_match_trk_chi2_dof_H->Fill(tk_RedChiSq);
    }
  }// Forward Eta Region
  else{
    cout << "E R R O R ! PixelRefitting::Loop(...) - Unexpected eta value of \"" << tk_Eta << "\". EXIT" << endl;
    exit(1);
  }

  return;
}


//****************************************************************************
void PixelRefitting::FillPixTrackMatchedTPHistos(TTPixelTrack pixFitTk,
						 int tp_Index)
//****************************************************************************
{

  // Get TP variables
  double tp_Pt  = TP_Pt->at(tp_Index);
  double tp_Eta = TP_Eta->at(tp_Index);
  double tp_Phi = TP_Phi->at(tp_Index);
  double tp_z0  = TP_POCAz->at(tp_Index);
  double tp_d0  = tp->GetD0(tp_Index);

  // Get TTPixelTrack variables
  double tk_Pt         = pixFitTk.getMomentum().Perp();
  double tk_Eta        = pixFitTk.getMomentum().Eta();
  double tk_Phi        = pixFitTk.getMomentum().Phi();
  double tk_z0         = pixFitTk.getZ0();
  double tk_d0         = pixFitTk.getD0();
  double tk_ChiSq      = pixFitTk.getChi2();
  double tk_RedChiSq   = pixFitTk.getChi2Red();
  
  // Generic Histograms
  h_match_tp_pt->Fill(tp_Pt);
  if (tp_Pt < _pt_L) h_match_tp_pt_L->Fill(tp_Pt);      
  h_match_tp_eta->Fill(tp_Eta);
  h_match_tp_phi->Fill(tp_Phi);
  h_match_tp_z0 ->Fill(tp_z0);
  h_match_tp_d0 ->Fill(tp_d0 * _fromCentimetresToMicrometers);

  FillEtaRegionsHistos(tp_Eta, tp_Pt, h_match_tp_pt_C, h_match_tp_pt_I, h_match_tp_pt_F);
  FillPtRegionsHistos(tp_Pt, tp_Eta, h_match_tp_eta_L, h_match_tp_eta_M, h_match_tp_eta_H);
  FillResidualHistos(tp_Index, pixFitTk);

  // 2D Histograms
  h_2d_dpt_eta   ->Fill(tp_Eta, fabs((tk_Pt - tp_Pt)));
  h_2d_dptRel_eta->Fill(tp_Eta, fabs((tk_Pt - tp_Pt)/tp_Pt));
  h_2d_deta_eta  ->Fill(tp_Eta, fabs(tk_Eta - tp_Eta));
  h_2d_dphi_eta  ->Fill(tp_Eta, fabs(tk_Phi - tp_Phi));
  h_2d_dz0_eta   ->Fill(tp_Eta, fabs(tk_z0  - tp_z0));
  h_2d_dd0_eta   ->Fill(tp_Eta, fabs(tk_d0  - tp_d0));    
  h_2d_logchi2_eta    ->Fill( tp_Eta, log(tk_ChiSq)    );
  h_2d_logchi2_dof_eta->Fill( tp_Eta, log(tk_RedChiSq) );

  return;
}

#endif // PixelRefitting_cxx

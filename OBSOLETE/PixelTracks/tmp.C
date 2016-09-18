#ifndef PixelTracks_cxx
#define PixelTracks_cxx

#include "PixelTracks.h"
#include "../utilities/constants.h"
#include "TFitResult.h"
#include "TF1.h"

//****************************************************************************
void PixelTracks::InitVars(void)
//****************************************************************************
{
  
  // Options
  bVerbose              = false;
  bPrintResolutions     = false;
  bPrintResolutions_Incl= false;
  bPrintResidualFitInfo = false;
  bSaveResidualFitInfo  = false;
  
  // TP Cuts
  TP_PT_MIN  =  0.2;
  TP_ETA_MAX =  2.5;
  TP_Z0_MAX  = 30.0;
  TP_DXY_MAX =  1.0;

  // TTPixelTrack Cuts
  TK_PIXHITS_MIN = 3;
  TK_PIXHITS_MAX = 1e6; 
  // TK_PIXHITS_TYPE.push_back(+1);
  // TK_PIXHITS_TYPE.push_back(-1);
  
  // Settings
  _pt_L              =   6.0; //  5.0
  _pt_H              =  16.0; // 15.0
  _eta_C             =   0.8; 
  _eta_F             =   1.6;
  _chiSq_overflow    = 100.0;
  _redChiSq_overflow =  15.0;
  _fromCentimetresToMicrometers = 10000.0;  // 1) cm->m: 100 2) m->micrometers*1000000;
  
  // Histogram Binning
  PT_MAX       = 50.0;
  PT_BINS      = 10;
  PT_BINWIDTH  = PT_MAX/(double)PT_BINS;
  ETA_MAX      = 2.4; // default: 2.5 
  ETA_BINS     = 12;  // default: 25
  ETA_BINWIDTH = ETA_MAX/(double)ETA_BINS;
  
  return;
}



//****************************************************************************
void PixelTracks::Loop()
//****************************************************************************
{
  
  if (fChain == 0) return;
  Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : std::min((int)fChain->GetEntries(), MaxEvents);
  
  ///////////////////////////////////////////////////////////
  /// Initialisations
  ///////////////////////////////////////////////////////////
  Long64_t nbytes = 0, nb = 0;
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
      bool bIsTkMatched     = (tp_NMatch > 0); // Was the TP matched to a TTTrack?       
      bool bIsTkMatchUnique = (tp_NMatch < 2); // Was the match unique? (31 July 2015) FIXME
      if ( !(bIsTkMatched * bIsTkMatchUnique) ) continue; 

      ///////////////////////////////////////////////////////////
      // TTPixelTracks
      ///////////////////////////////////////////////////////////
      int tk_Index    = TP_TTTrackIndex->at(tp_Index);
      int pixTk_Index = s->GetPixelIndexOfTrack(tk_Index);
      if (pixTk_Index < 0) continue;
      
      // Get TTPixelTrack variables
      double tk_Eta        = L1PixTks_Eta->at(pixTk_Index);
      int tk_PixHits       = L1PixTks_NPixHits->at(pixTk_Index);

      
      ///////////////////////////////////////////////////////////
      // TTPixelTracks Cuts
      ///////////////////////////////////////////////////////////
      bool bPassNumOfHits_Min     = (tk_PixHits >=  TK_PIXHITS_MIN);
      bool bPassNumOfHits_Max     = (tk_PixHits <= TK_PIXHITS_MAX);
      bool bPassNumOfHits_Exact   = (tk_PixHits == TK_PIXHITS_MIN);
      bool bHasMustHitLayerOrDisk = s->HasMustPixelHitInLayeOrDisk(pixTk_Index, TK_PIXHITS_TYPE);
      // bool bPassAllCuts = bPassNumOfHits_Min * bPassNumOfHits_Max * bPassNumOfHits_Exact * bHasMustHitLayerOrDisk;
      bool bPassAllCuts = bPassNumOfHits_Min * bPassNumOfHits_Max * bHasMustHitLayerOrDisk;
      if( !bPassAllCuts ) continue;
                 
      // Fill Pixel Hit Histograms
      FillCandPixelHitHistos(pixTk_Index, false);
      FillPixelHitHistos(pixTk_Index, false);
      FillSharedPixelHitHistos(pixTk_Index, tp_Index, pixHits_XYZ, pixHits_TTPixelTrackIndex, false);

      // Fill TP<-->Track Matched Histograms
      FillTPMatchedPixTrackHistos(pixTk_Index);
      FillPixTrackMatchedTPHistos(pixTk_Index, tp_Index);

      // Print TP & TTPixelTrack properties
      if (bVerbose) tp->PrintTPProperties(tp_Index);
      if (bVerbose) s->PrintPixTrackProperties(pixTk_Index, true);
      
      
    } // For-loop: TPs

    if (!bVerbose) auxTools.ProgressBar(jentry, nEntries, 100, 150);
    
  }// For-loop: Entries


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
void PixelTracks::BookHistos(void)
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
  histoTools.BookHisto_1D( h_efficiency_pt   , "eff_pt"   , 100,   0.0,  +100.0  );
  histoTools.BookHisto_1D( h_efficiency_pt_L , "eff_pt_L" ,  50,   0.0,    +5.0  );
  histoTools.BookHisto_1D( h_efficiency_pt_C , "eff_pt_C" , 100,   0.0,  +100.0  );
  histoTools.BookHisto_1D( h_efficiency_pt_I , "eff_pt_I" , 100,   0.0,  +100.0  ); 
  histoTools.BookHisto_1D( h_efficiency_pt_F , "eff_pt_F" , 100,   0.0,  +100.0  );
  histoTools.BookHisto_1D( h_efficiency_eta  , "eff_eta"  ,  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_eta_L, "eff_eta_L",  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_eta_M, "eff_eta_M",  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_eta_H, "eff_eta_H",  26,  -2.6,    +2.6  );
  histoTools.BookHisto_1D( h_efficiency_phi  , "eff_phi"  ,  64,  -3.2,    +3.2  );
  histoTools.BookHisto_1D( h_efficiency_z0   , "eff_z0"   ,  80, -20.0,   +20.0  );
  histoTools.BookHisto_1D( h_efficiency_d0   , "eff_d0"   , 100,-50.0,   +50.0  );

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
  histoTools.BookHisto_1D( h_residual_pt     , "res_pt"     , 100, -5.0  , 5.0  );
  histoTools.BookHisto_1D( h_residual_pt_C   , "res_pt_C"   , 100, -5.0  , 5.0  );
  histoTools.BookHisto_1D( h_residual_pt_I   , "res_pt_I"   , 100, -5.0  , 5.0  );
  histoTools.BookHisto_1D( h_residual_pt_F   , "res_pt_F"   , 100, -5.0  , 5.0  );

  histoTools.BookHisto_1D( h_residual_ptRel  , "res_ptRel"  , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_ptRel_C, "res_ptRel_C", 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_ptRel_I, "res_ptRel_I", 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_ptRel_F, "res_ptRel_F", 1000, -1.0  , 1.0  );

  histoTools.BookHisto_1D( h_residual_eta    , "res_eta"    , 400, -0.01 , 0.01 );
  histoTools.BookHisto_1D( h_residual_eta_C  , "res_eta_C"  , 200, -0.01 , 0.01 );
  histoTools.BookHisto_1D( h_residual_eta_I  , "res_eta_I"  , 200, -0.01 , 0.01 );
  histoTools.BookHisto_1D( h_residual_eta_F  , "res_eta_F"  , 200, -0.01 , 0.01 );

  histoTools.BookHisto_1D( h_residual_phi    , "res_phi"    , 100, -0.005, 0.005);
  histoTools.BookHisto_1D( h_residual_phi_C  , "res_phi_C"  , 100, -0.005, 0.005);
  histoTools.BookHisto_1D( h_residual_phi_I  , "res_phi_I"  , 100, -0.005, 0.005);
  histoTools.BookHisto_1D( h_residual_phi_F  , "res_phi_F"  , 100, -0.005, 0.005);

  histoTools.BookHisto_1D( h_residual_z0     , "res_z0"     , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_z0_C   , "res_z0_C"   , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_z0_I   , "res_z0_I"   , 1000, -1.0  , 1.0  );
  histoTools.BookHisto_1D( h_residual_z0_F   , "res_z0_F"   , 1000, -1.0  , 1.0  );

  histoTools.BookHisto_1D( h_residual_d0     , "res_d0"     , 500, -0.5  , 0.5  );
  histoTools.BookHisto_1D( h_residual_d0_C   , "res_d0_C"   , 500, -0.5  , 0.5  );
  histoTools.BookHisto_1D( h_residual_d0_I   , "res_d0_I"   , 500, -0.5  , 0.5  );
  histoTools.BookHisto_1D( h_residual_d0_F   , "res_d0_F"   , 500, -0.5  , 0.5  );

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
	
    histoTools.BookHisto_1D( h_residualVsEta_pt[i]  , "h_residualVsEta_pt_"   + etarange[i], 100, -5.00, +5.00);
    histoTools.BookHisto_1D( h_residualVsEta_pt_L[i], "h_residualVsEta_pt_L_" + etarange[i], 100, -0.25, +0.25);
    histoTools.BookHisto_1D( h_residualVsEta_pt_M[i], "h_residualVsEta_pt_M_" + etarange[i], 100, -0.25, +0.25);
    histoTools.BookHisto_1D( h_residualVsEta_pt_H[i], "h_residualVsEta_pt_H_" + etarange[i], 100, -0.25, +0.25);

    histoTools.BookHisto_1D( h_residualVsEta_ptRel[i]  , "h_residualVsEta_ptRel_"  + etarange[i], 200, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsEta_ptRel_L[i], "h_residualVsEta_ptRel_L_"+ etarange[i], 200, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsEta_ptRel_M[i], "h_residualVsEta_ptRel_M_"+ etarange[i], 200, -0.50, +0.50);
    histoTools.BookHisto_1D( h_residualVsEta_ptRel_H[i], "h_residualVsEta_ptRel_H_"+ etarange[i], 200, -0.50, +0.50);

    histoTools.BookHisto_1D( h_residualVsEta_eta[i]  , "h_residualVsEta_eta_"   + etarange[i], 400/nBinsXFactor  , -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_eta_L[i], "h_residualVsEta_eta_L_" + etarange[i], 200/nBinsXFactor_L, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_eta_M[i], "h_residualVsEta_eta_M_" + etarange[i], 400/nBinsXFactor_M, -0.01, +0.01); 
    histoTools.BookHisto_1D( h_residualVsEta_eta_H[i], "h_residualVsEta_eta_H_" + etarange[i], 400/nBinsXFactor_H, -0.01, +0.01);

    histoTools.BookHisto_1D( h_residualVsEta_phi[i]  , "h_residualVsEta_phi_"   + etarange[i], 100, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_phi_L[i], "h_residualVsEta_phi_L_" + etarange[i], 100, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_phi_M[i], "h_residualVsEta_phi_M_" + etarange[i], 100, -0.01, +0.01);
    histoTools.BookHisto_1D( h_residualVsEta_phi_H[i], "h_residualVsEta_phi_H_" + etarange[i], 100, -0.01, +0.01);

    histoTools.BookHisto_1D( h_residualVsEta_z0[i]  , "h_residualVsEta_z0_"   + etarange[i], 1000/nBinsXFactor  , -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_z0_L[i], "h_residualVsEta_z0_L_" + etarange[i], 1000/nBinsXFactor_L, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_z0_M[i], "h_residualVsEta_z0_M_" + etarange[i], 1000/nBinsXFactor_M, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_z0_H[i], "h_residualVsEta_z0_H_" + etarange[i], 1000/nBinsXFactor_H, -0.5, +0.5);

    histoTools.BookHisto_1D( h_residualVsEta_d0[i]  , "h_residualVsEta_d0_"   + etarange[i], 2000/nBinsXFactor  , -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_d0_L[i], "h_residualVsEta_d0_L_" + etarange[i], 1000/nBinsXFactor_L, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_d0_M[i], "h_residualVsEta_d0_M_" + etarange[i], 1000/nBinsXFactor_M, -0.5, +0.5);
    histoTools.BookHisto_1D( h_residualVsEta_d0_H[i], "h_residualVsEta_d0_H_" + etarange[i], 2000/nBinsXFactor_H, -0.5, +0.5);

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

  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_N            , "pixTk_sharedPixHits_N"            ,   30,   -0.5,  +29.5 );
  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_Rho          , "pixTk_sharedPixHits_Rho"          , 1000,   +0.0,  +20.0 );
  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_Z            , "pixTk_sharedPixHits_Z"            ,  240,  -60.0,  +60.0 );
  histoTools.BookHisto_1D( h_pixTk_sharedPixHits_Type         , "pixTk_sharedPixHits_Type"         ,   10,   -4.5,   +5.5 );
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
void PixelTracks::FillHistoBinWithRMS(TH1D *hToFill,
					TH1D *hToExtractRMS,
					int binNumber)
//****************************************************************************
{
  
  hToFill->SetBinContent(binNumber, hToExtractRMS->GetRMS() );
  hToFill->SetBinError  (binNumber, hToExtractRMS->GetRMSError() );

  return;
}


  
//****************************************************************************
void PixelTracks::FitAndFillHistoBin(TH1D *hToFill,
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
    FillHistoBinWithRMS(hToFill, hToFit, binNumber);
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
  TFitResultPtr r = histoTools.FitFunction( hToFit, fitFunction, "Q S ", "", fitRangeLow, fitRangeHigh, v_fitArea, redChiSqMax);
  
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
    histoTools.FindSymmetricFitRange(hToFit, 0.8, fitRangeLow, fitRangeHigh);
    hToFit->GetXaxis()->SetRange(fitRangeLow, fitRangeHigh);
    
    hToFill->SetBinContent(binNumber, hToFit->GetRMS() );
    hToFill->SetBinError  (binNumber, hToFit->GetRMSError() );
  }

   
  return;
}


//****************************************************************************
void PixelTracks::FitAndFillHistoBinSL(TH1D *hToFill,
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
std::vector<double> PixelTracks::GetFitAreaVector(const double minTotalArea)
//****************************************************************************
{

  // Fit area from 0.5 -> minTotalArea/2
  std::vector<double> v_fitArea;  
  for (int i = 0; i <= 20; i++){
    double area = 0.5 - i*0.02;
    if (2*area < minTotalArea) break;
    v_fitArea.push_back(area);
  }
  // auxTools.PrintVector(v_fitArea);
  
  return v_fitArea;
}


//****************************************************************************
void PixelTracks::FinaliseHistos(void)
//****************************************************************************
{
  
  // Loop over all pt bins
  for (int i=0; i < PT_BINS; i++) {

    // Cannot fit "gaus" to phi due to asymmetric distribution (lorentz-drift)
    FitAndFillHistoBinSL( h_resolutionVsPt_pt  , h_residualVsPt_pt[i]  , i+1, +0.005);
    FitAndFillHistoBinSL( h_resolutionVsPt_pt_C, h_residualVsPt_pt_C[i], i+1, +0.005);
    FitAndFillHistoBinSL( h_resolutionVsPt_pt_I, h_residualVsPt_pt_I[i], i+1, +0.005);
    FitAndFillHistoBinSL( h_resolutionVsPt_pt_F, h_residualVsPt_pt_F[i], i+1, +0.005);

    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel  , h_residualVsPt_ptRel[i]  , i+1, +0.005);
    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel_C, h_residualVsPt_ptRel_C[i], i+1, +0.005);
    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel_I, h_residualVsPt_ptRel_I[i], i+1, +0.005);
    FitAndFillHistoBinSL( h_resolutionVsPt_ptRel_F, h_residualVsPt_ptRel_F[i], i+1, +0.005);
    
    // Fit "gaus" to eta to get resolution
    FitAndFillHistoBinSL( h_resolutionVsPt_eta  , h_residualVsPt_eta[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_eta_C, h_residualVsPt_eta_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_eta_I, h_residualVsPt_eta_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_eta_F, h_residualVsPt_eta_F[i], i+1, +0.005 );

    // Cannot fit "gaus" to phi due to asymmetric distribution (lorentz-drift)
    FitAndFillHistoBinSL( h_resolutionVsPt_phi  , h_residualVsPt_phi[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_phi_C, h_residualVsPt_phi_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_phi_I, h_residualVsPt_phi_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_phi_F, h_residualVsPt_phi_F[i], i+1, +0.005 );

    // Fit "gaus" to z0 to get resolution
    FitAndFillHistoBinSL( h_resolutionVsPt_z0  , h_residualVsPt_z0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_z0_C, h_residualVsPt_z0_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_z0_I, h_residualVsPt_z0_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_z0_F, h_residualVsPt_z0_F[i], i+1, +0.005 );

    // Fit "gaus" to d0 to get resolution
    FitAndFillHistoBinSL( h_resolutionVsPt_d0  , h_residualVsPt_d0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_d0_C, h_residualVsPt_d0_C[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_d0_I, h_residualVsPt_d0_I[i], i+1, +0.005 );
    FitAndFillHistoBinSL( h_resolutionVsPt_d0_F, h_residualVsPt_d0_F[i], i+1, +0.005 );

  }

  // Loop over all eta bins
  for (int i=0; i < ETA_BINS; i++) {

    FitAndFillHistoBinSL(h_resolutionVsEta_pt  , h_residualVsEta_pt[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_pt_L, h_residualVsEta_pt_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_pt_M, h_residualVsEta_pt_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_pt_H, h_residualVsEta_pt_H[i], i+1, +0.005 );
  
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel  , h_residualVsEta_ptRel[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel_L, h_residualVsEta_ptRel_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel_M, h_residualVsEta_ptRel_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_ptRel_H, h_residualVsEta_ptRel_H[i], i+1, +0.005 );

    FitAndFillHistoBinSL(h_resolutionVsEta_eta  , h_residualVsEta_eta[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_eta_L, h_residualVsEta_eta_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_eta_M, h_residualVsEta_eta_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_eta_H, h_residualVsEta_eta_H[i], i+1, +0.005 );
    
    FitAndFillHistoBinSL(h_resolutionVsEta_phi  , h_residualVsEta_phi[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_phi_L, h_residualVsEta_phi_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_phi_M, h_residualVsEta_phi_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_phi_H, h_residualVsEta_phi_H[i], i+1, +0.005 );

    FitAndFillHistoBinSL(h_resolutionVsEta_z0  , h_residualVsEta_z0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_z0_L, h_residualVsEta_z0_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_z0_M, h_residualVsEta_z0_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_z0_H, h_residualVsEta_z0_H[i], i+1, +0.005 );

    FitAndFillHistoBinSL(h_resolutionVsEta_d0  , h_residualVsEta_d0[i]  , i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_d0_L, h_residualVsEta_d0_L[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_d0_M, h_residualVsEta_d0_M[i], i+1, +0.005 );
    FitAndFillHistoBinSL(h_resolutionVsEta_d0_H, h_residualVsEta_d0_H[i], i+1, +0.005 );
    
  }
  
  return;
}

    

//****************************************************************************
void PixelTracks::MakeEfficiencyHisto(TH1D* h_eff, 
					TH1D* h_match_tp, 
					TH1D* h_tp)
//****************************************************************************
{

  h_eff->Divide(h_match_tp, h_tp, 1.0, 1.0, "B");
  
  return;
}
 
 

//****************************************************************************
void PixelTracks::PrintSettings(void)
//****************************************************************************
{

  // Inform user of settings  
  Table settings("Symbol | Definition | Units | Description | Comments |", "Text");
  settings.AddRowColumn(0, "bVerbose");
  settings.AddRowColumn(0, auxTools.ToString(bVerbose) );
  settings.AddRowColumn(0, "-" );
  settings.AddRowColumn(0, "Enable to print event-by-event info");
  settings.AddRowColumn(0, "-");

  settings.AddRowColumn(1, "bPrintResolutions");
  settings.AddRowColumn(1, auxTools.ToString(bPrintResolutions) );
  settings.AddRowColumn(1, "-" );
  settings.AddRowColumn(1, "Enable to print resolutions for pT regions");
  settings.AddRowColumn(1, "-");

  settings.AddRowColumn(2, "bPrintResolutions_Incl");
  settings.AddRowColumn(2, auxTools.ToString(bPrintResolutions_Incl) );
  settings.AddRowColumn(2, "-" );
  settings.AddRowColumn(2, "Enable to print inclusive resolutions");
  settings.AddRowColumn(2, "-");

  settings.AddRowColumn(3, "bPrintResidualFitInfo");
  settings.AddRowColumn(3, auxTools.ToString(bPrintResidualFitInfo) );
  settings.AddRowColumn(3, "-" );
  settings.AddRowColumn(3, "Enable to print all info on residual fits attempted");
  settings.AddRowColumn(3, "-");

  settings.AddRowColumn(4, "bSaveResidualFitInfo");
  settings.AddRowColumn(4, auxTools.ToString(bSaveResidualFitInfo) );
  settings.AddRowColumn(4, "-" );
  settings.AddRowColumn(4, "Enable to save to txt file all info for the residual fits attempted");
  settings.AddRowColumn(4, "-");

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
  settings.AddRowColumn(8, "-");
    
  settings.AddRowColumn(9, "_pt_M");
  settings.AddRowColumn(9, auxTools.ToString(_pt_L) + " GeV/c < pT < " + auxTools.ToString(_pt_H) );
  settings.AddRowColumn(9, "GeV/c");
  settings.AddRowColumn(9, "Middle pT Range");
  settings.AddRowColumn(9, "-");
    
  settings.AddRowColumn(10, "_pt_H");
  settings.AddRowColumn(10, "pT > " + auxTools.ToString(_pt_H) );
  settings.AddRowColumn(10, "GeV/c");
  settings.AddRowColumn(10, "High pT Range");
  settings.AddRowColumn(10, "-");
  
  settings.AddRowColumn(11, "_chiSq_overflow");
  settings.AddRowColumn(11, auxTools.ToString(_chiSq_overflow) );
  settings.AddRowColumn(11, "-" );
  settings.AddRowColumn(11, "Value of overflow bin for chi-squared of L1 tracks" );
  settings.AddRowColumn(11, "-");
  
  settings.AddRowColumn(12, "_redChiSq_overflow");
  settings.AddRowColumn(12, auxTools.ToString(_redChiSq_overflow) );
  settings.AddRowColumn(12, "-" );
  settings.AddRowColumn(12, "Value of overflow bin for reduced-chi-squared of L1 tracks" );
  settings.AddRowColumn(12, "-");
  
  settings.Print();
  
  return;
}



//****************************************************************************
void PixelTracks::PrintTPCuts(void)
//****************************************************************************
{

  Table tpCuts("Variable | Cut Direction | Cut Value | Units | Description | Comments |", "Text");
  tpCuts.AddRowColumn(0, "TP_PT_MIN");
  tpCuts.AddRowColumn(0, " >= " );
  tpCuts.AddRowColumn(0, auxTools.ToString(TP_PT_MIN) );
  tpCuts.AddRowColumn(0, "GeV/c" );
  tpCuts.AddRowColumn(0, "Minimum pT for TP to be considered within acceptance" );
  tpCuts.AddRowColumn(0, "-" );

  tpCuts.AddRowColumn(1, "TP_ETA_MAX");
  tpCuts.AddRowColumn(1, " <= " );
  tpCuts.AddRowColumn(1, auxTools.ToString(TP_ETA_MAX) );
  tpCuts.AddRowColumn(1, "-" );
  tpCuts.AddRowColumn(1, "Maximum absolute eta value TP to be considered within acceptance" );
  tpCuts.AddRowColumn(1, "-" );
  
  tpCuts.AddRowColumn(2, "TP_Z0_MAX");
  tpCuts.AddRowColumn(2, " <= " );
  tpCuts.AddRowColumn(2, auxTools.ToString(TP_Z0_MAX) );
  tpCuts.AddRowColumn(2, "cm" );
  tpCuts.AddRowColumn(2, "Maximum absolute z0 value TP to be considered within acceptance" );
  tpCuts.AddRowColumn(2, "-" );
  
  tpCuts.AddRowColumn(3, "TP_DXY_MAX");
  tpCuts.AddRowColumn(3, " <= " );
  tpCuts.AddRowColumn(3, auxTools.ToString(TP_DXY_MAX) );
  tpCuts.AddRowColumn(3, "cm" );
  tpCuts.AddRowColumn(3, "Maximum dxy of TPs considered - Select only those from near the IP");
  tpCuts.AddRowColumn(3, "Affects d0 Resolution!" );

  tpCuts.Print();
  
  return;
}


//****************************************************************************
void PixelTracks::PrintPixTrackCuts(void)
//****************************************************************************
{

  // Inform user of settings  
  Table pixTkCuts("Variable | Cut Direction | Cut Value | Units | Description | Comments |", "Text");
  pixTkCuts.AddRowColumn(0, "TK_PIXHITS_MIN");
  pixTkCuts.AddRowColumn(0, " >= " );
  pixTkCuts.AddRowColumn(0, auxTools.ToString(TK_PIXHITS_MIN) );
  pixTkCuts.AddRowColumn(0, "-" );
  pixTkCuts.AddRowColumn(0, "Minimum number of Pixel Fit-Hits for TTPixelTracks" );
  pixTkCuts.AddRowColumn(0, "-" );

  pixTkCuts.AddRowColumn(1, "TK_PIXHITS_MAX"); 
  pixTkCuts.AddRowColumn(1, " <= " );
  pixTkCuts.AddRowColumn(1, auxTools.ToString(TK_PIXHITS_MAX) );
  pixTkCuts.AddRowColumn(1, "-" );
  pixTkCuts.AddRowColumn(1, "Maximum number of Pixel Fit-Hits per TTPixelTracks" );
  pixTkCuts.AddRowColumn(1, "-" );
  
  pixTkCuts.AddRowColumn(2, "TK_PIXHITS_TYPE");
  pixTkCuts.AddRowColumn(2, " == " );
  pixTkCuts.AddRowColumn(2, auxTools.ConvertIntVectorToString(TK_PIXHITS_TYPE) );
  pixTkCuts.AddRowColumn(2, "-" );
  pixTkCuts.AddRowColumn(2, "Pixel Hit Types that TTPixelTracks must at least satisfy ONE of them" );
  pixTkCuts.AddRowColumn(2, "-" );
  
  pixTkCuts.Print();
  
  return;
}


//****************************************************************************
void PixelTracks::PrintResolutions(void)
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
double PixelTracks::GetHistoMeanInRange(TH1D* histo, 
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
bool PixelTracks::IsWithinTrackerAcceptance(double pt, 
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
bool PixelTracks::IsWithinEtaRegion(string etaRegion, 
				      double eta)
//****************************************************************************
{

  bool bWithinEtaRegion = false;
  if ( etaRegion.compare("Central") == 0 )           bWithinEtaRegion = (fabs(eta) <= _eta_C);
  else if ( etaRegion.compare("Intermediate") == 0 ) bWithinEtaRegion = (fabs(eta) <= _eta_F && fabs(eta) > _eta_C);
  else if ( etaRegion.compare("Forward") == 0 )      bWithinEtaRegion = (fabs(eta) > _eta_F);
  else{
    cout << "E R R O R ! PixelTracks::IsWithinEtaRegion(...) - Invalid eta region type \"" << etaRegion << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinEtaRegion;  
}


//****************************************************************************
bool PixelTracks::IsWithinPtRange(string ptRange, 
				    double pt)
//****************************************************************************
{

  bool bWithinPtRange = false;
  if ( ptRange.compare("Low") == 0 )         bWithinPtRange = (pt <= _pt_L);
  else if ( ptRange.compare("Middle") == 0 ) bWithinPtRange = (pt <= _pt_H && pt > _pt_L);
  else if ( ptRange.compare("High") == 0 )   bWithinPtRange = (pt > _pt_H);
  else{
    cout << "E R R O R ! PixelTracks::IsWithinPtRange(...) - Invalid pt range type \"" << ptRange << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinPtRange;  
}

//****************************************************************************
void PixelTracks::FillPixelHitHistos(int pixTk_Index,
				       bool bPrintPixInfo)
//****************************************************************************
{
  
  if (pixTk_Index < 0) return;

  // Create Pixel-Hit Table
  Table pixInfo("Hits | X | Y | Z | R | Phi | Type", "Text");
  
  // Get variables
  double tk_Eta  = L1PixTks_Eta->at(pixTk_Index);
  vector<double> pixHits_X    = L1PixTks_PixHits_X->at(pixTk_Index);
  vector<double> pixHits_Y    = L1PixTks_PixHits_Y->at(pixTk_Index);
  vector<double> pixHits_Z    = L1PixTks_PixHits_Z->at(pixTk_Index);
  vector<double> pixHits_R    = L1PixTks_PixHits_R->at(pixTk_Index);
  vector<double> pixHits_Phi  = L1PixTks_PixHits_Phi->at(pixTk_Index);
  vector< int >  pixHits_Type = L1PixTks_PixHits_Type->at(pixTk_Index);
  int tk_PixHits     = (int) pixHits_X.size();
  
  // Fill Histograms
  h_pixTk_pixHits_N->Fill(tk_PixHits);
 
  // For-loop: Pixel Hits (best-fit)  
  for(int i = 0; i < tk_PixHits; i++){

    pixInfo.AddRowColumn(i, auxTools.ToString(tk_PixHits) );
    pixInfo.AddRowColumn(i, auxTools.ToString(pixHits_X.at(i)) );
    pixInfo.AddRowColumn(i, auxTools.ToString(pixHits_Y.at(i)) );
    pixInfo.AddRowColumn(i, auxTools.ToString(pixHits_Z.at(i)) );
    pixInfo.AddRowColumn(i, auxTools.ToString(pixHits_R.at(i)) );
    pixInfo.AddRowColumn(i, auxTools.ToString(pixHits_Phi.at(i)) );
    pixInfo.AddRowColumn(i, auxTools.ToString(pixHits_Type.at(i)) );

    
    // Fill Histos
    h_pixTk_pixHits_Rho   ->Fill( pixHits_R.at(i) );
    h_pixTk_pixHits_Z     ->Fill( pixHits_Z.at(i) );
    h_pixTk_pixHits_Type  ->Fill( pixHits_Type.at(i) );

    h_pixTk_pixHits_ZVsRho->Fill( pixHits_Z.at(i), pixHits_R.at(i) );   
    FillEtaRegionsHistos2D(tk_Eta, pixHits_Z.at(i), pixHits_R.at(i), h_pixTk_pixHits_ZVsRho_C, h_pixTk_pixHits_ZVsRho_I, h_pixTk_pixHits_ZVsRho_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_pixHits_ZVsRho_EtaGE2->Fill( pixHits_Z.at(i), pixHits_R.at(i) );

    h_pixTk_pixHits_XVsY->Fill( pixHits_X.at(i), pixHits_Y.at(i) );
    FillEtaRegionsHistos2D(tk_Eta, pixHits_X.at(i), pixHits_Y.at(i), h_pixTk_pixHits_XVsY_C, h_pixTk_pixHits_XVsY_I, h_pixTk_pixHits_XVsY_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_pixHits_XVsY_EtaGE2->Fill( pixHits_X.at(i), pixHits_Y.at(i) );
    
  }// for(int i = 0; i < tk_PixHits; i++){
  
  if (bPrintPixInfo) pixInfo.Print();
  return;
}


//****************************************************************************
void PixelTracks::FillCandPixelHitHistos(int pixTk_Index,
					   bool bPrintPixInfo)
//****************************************************************************
{
  if (pixTk_Index < 0) return;

  // Create Pixel-Hit Table
  Table pixInfo("Candidate Hits | X | Y | Z | R | Phi | Type", "Text");
  
  // Get variables
  double tk_Eta = L1PixTks_Eta->at(pixTk_Index);
  vector<double> pixHits_X    = L1PixTks_CandPixHits_X->at(pixTk_Index);
  vector<double> pixHits_Y    = L1PixTks_CandPixHits_Y->at(pixTk_Index);
  vector<double> pixHits_Z    = L1PixTks_CandPixHits_Z->at(pixTk_Index);
  vector<double> pixHits_R    = L1PixTks_CandPixHits_R->at(pixTk_Index);
  vector<double> pixHits_Phi  = L1PixTks_CandPixHits_Phi->at(pixTk_Index);
  vector< int >  pixHits_Type = L1PixTks_CandPixHits_Type->at(pixTk_Index);
  int tk_PixHits = (int) pixHits_X.size();  
  
  // Fill Histograms 
  h_pixTk_candPixHits_N->Fill(tk_PixHits);

  // For-loop: Pixel Hits (candidates)
  for(int j = 0; j < tk_PixHits; j++){

    pixInfo.AddRowColumn(j, auxTools.ToString(tk_PixHits) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_X.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_Y.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_Z.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_R.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_Phi.at(j)) );
    pixInfo.AddRowColumn(j, auxTools.ToString(pixHits_Type.at(j)) );

    // Fill Histos
    h_pixTk_candPixHits_Rho   ->Fill( pixHits_R.at(j) );
    h_pixTk_candPixHits_Z     ->Fill( pixHits_Z.at(j) );
    h_pixTk_candPixHits_Type  ->Fill( pixHits_Type.at(j) );

    h_pixTk_candPixHits_ZVsRho->Fill( pixHits_Z.at(j), pixHits_R.at(j) );
    FillEtaRegionsHistos2D(tk_Eta, pixHits_Z.at(j), pixHits_R.at(j), h_pixTk_candPixHits_ZVsRho_C, h_pixTk_candPixHits_ZVsRho_I, h_pixTk_candPixHits_ZVsRho_F);
    if (fabs(tk_Eta) >= 2.0) h_pixTk_candPixHits_ZVsRho_EtaGE2->Fill( pixHits_Z.at(j), pixHits_R.at(j) );
    
    h_pixTk_candPixHits_XVsY->Fill( pixHits_X.at(j), pixHits_Y.at(j) );
    FillEtaRegionsHistos2D(tk_Eta, pixHits_X.at(j), pixHits_Y.at(j), h_pixTk_candPixHits_XVsY_C, h_pixTk_candPixHits_XVsY_I, h_pixTk_candPixHits_XVsY_F);    
    if (fabs(tk_Eta) >= 2.0) h_pixTk_candPixHits_XVsY_EtaGE2->Fill( pixHits_X.at(j), pixHits_Y.at(j) );
      
  }// for(int j = 0; j < pixTk_CandPixHits; j++){

  if (bPrintPixInfo) pixInfo.Print();
  return;
}


//****************************************************************************
void PixelTracks::FillSharedPixelHitHistos(int pixTk_Index,
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
void PixelTracks::FillResidualHistos(int tp_Index,
				       int pixTk_Index)
//****************************************************************************
{

  // Get TP variables
  double tp_Pt  = TP_Pt->at(tp_Index);
  double tp_Eta = TP_Eta->at(tp_Index);
  double tp_Phi = TP_Phi->at(tp_Index);
  double tp_z0  = TP_POCAz->at(tp_Index);
  double tp_d0  = tp->GetD0(tp_Index);

  // Get TTPixelTrack variables
  double tk_Pt  = L1PixTks_Pt ->at(pixTk_Index);
  double tk_Eta = L1PixTks_Eta->at(pixTk_Index);
  double tk_Phi = L1PixTks_Phi->at(pixTk_Index);
  double tk_z0  = L1PixTks_POCAz->at(pixTk_Index);
  double tk_d0  = s->GetPixTrackD0(pixTk_Index);
  
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
    cout << "E R R O R ! PixelTracks::FillResidualHistos(...) - Unexpected eta of \"" << tk_Eta << "\". EXIT" << endl;
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
      cout << "E R R O R ! PixelTracks::FillResidualHistos(...) - Unexpected Eta value of \"" << tk_Eta << "\". EXIT" << endl;
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
      cout << "E R R O R ! PixelTracks::Loop(...) - Unexpected Track-Pt of \"" << tk_Pt << "\". EXIT" << endl;
      exit(1);
    }	  

  } // for (int k=0; k < ETA_BINS; k++) {

  return;
}


//****************************************************************************
void PixelTracks::FillEtaRegionsHistos(double eta,
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
    cout << "E R R O R ! PixelTracks::FillEtaRegionsHistos(...) - Unexpected Eta value of \"" << eta << "\". EXIT" << endl;
    exit(1);
  }
  
  return;
}


//****************************************************************************
void PixelTracks::FillEtaRegionsHistos2D(double eta,
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
    cout << "E R R O R ! PixelTracks::FillEtaRegionsHistos(...) - Unexpected Eta value of \"" << eta << "\". EXIT" << endl;
    exit(1);
  }
  
  return;
}


//****************************************************************************
void PixelTracks::FillPtRegionsHistos(double pt,
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
    cout << "E R R O R ! PixelTracks::FillPtRegionsHistos(...) - Unexpected pT value of \"" << pt << "\". EXIT" << endl;
    exit(1);
  }
  
  return;  
}


//****************************************************************************
void PixelTracks::FillTPMatchedPixTrackHistos(int pixTk_Index)
//****************************************************************************
{

  // Get TTPixelTrack variables
  double tk_Pt         = L1PixTks_Pt ->at(pixTk_Index);
  double tk_Eta        = L1PixTks_Eta->at(pixTk_Index);
  double tk_ChiSq      = L1PixTks_ChiSquared->at(pixTk_Index);
  double tk_RedChiSq   = s->GetPixTrackRedChiSq(pixTk_Index);

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
      cout << "E R R O R ! PixelTracks::Loop(...) - Unexpected tk pt value of \"" << tk_Pt << "\". EXIT" << endl;
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
      cout << "E R R O R ! PixelTracks::Loop(...) - Unexpected tk pt value of \"" << tk_Pt << "\". EXIT" << endl;
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
    cout << "E R R O R ! PixelTracks::Loop(...) - Unexpected eta value of \"" << tk_Eta << "\". EXIT" << endl;
    exit(1);
  }

  return;
}


//****************************************************************************
void PixelTracks::FillPixTrackMatchedTPHistos(int pixTk_Index,
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
  double tk_Pt         = L1PixTks_Pt ->at(pixTk_Index);
  double tk_Eta        = L1PixTks_Eta->at(pixTk_Index);
  double tk_Phi        = L1PixTks_Phi->at(pixTk_Index);
  double tk_z0         = L1PixTks_POCAz->at(pixTk_Index);
  double tk_d0         = s->GetPixTrackD0(pixTk_Index);
  double tk_ChiSq      = L1PixTks_ChiSquared->at(pixTk_Index);
  double tk_RedChiSq   = s->GetPixTrackRedChiSq(pixTk_Index);
  
  // Generic Histograms
  h_match_tp_pt->Fill(tp_Pt);
  if (tp_Pt < _pt_L) h_match_tp_pt_L->Fill(tp_Pt);      
  h_match_tp_eta->Fill(tp_Eta);
  h_match_tp_phi->Fill(tp_Phi);
  h_match_tp_z0 ->Fill(tp_z0);
  h_match_tp_d0 ->Fill(tp_d0 * _fromCentimetresToMicrometers);

  FillEtaRegionsHistos(tp_Eta, tp_Pt, h_match_tp_pt_C, h_match_tp_pt_I, h_match_tp_pt_F);
  FillPtRegionsHistos(tp_Pt, tp_Eta, h_match_tp_eta_L, h_match_tp_eta_M, h_match_tp_eta_H);
  FillResidualHistos(tp_Index, pixTk_Index);

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

#endif // PixelTracks_cxx

#ifndef PixelRefitting_h
#define PixelRefitting_h

#include <iostream>

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/MCTools.C"
#include "../utilities/AuxTools.C"
#include "../utilities/HistoTools.C"
#include "../utilities/L1Tracks.C"
#include "../utilities/TrackingParticles.C"
#include "../utilities/L1PixelTrackFit.C"
#include "../utilities/TTPixelTrack.C"

class PixelRefitting : public TreeAnalyserMC{
 public:
 PixelRefitting(const string SamplePath,
	       const std::string SampleName,
	       const std::string text_, 
	       const int maxEvents_ = -1, 
	       TTree* tree=0) : 
  TreeAnalyserMC("PixelRefitting", SamplePath, SampleName, text_, maxEvents_, tree) {};
  virtual void Loop();
  
 private:
  void InitVars(void);

  void BookHistos(void);

  void FinaliseHistos(void);

  void FillTPMatchedPixTrackHistos(TTPixelTrack pixFitTk);

  void FillPixTrackMatchedTPHistos(TTPixelTrack pixFitTk,
				   int tp_Index);
			   
  void FillPtRegionsHistos(double pt,
			   double variable,
			   TH1D *hLow,
			   TH1D *hMiddle,
			   TH1D *hHigh);

  void FillEtaRegionsHistos(double eta,
			    double variable,
			    TH1D *hCentral,
			    TH1D *hIntermediate,
			    TH1D *hForward);

  void FillEtaRegionsHistos2D(double eta,
			      double variableX,
			      double variableY,
			      TH2D *hCentral,
			      TH2D *hIntermediate,
			      TH2D *hForward);

  void FillPixelHitHistos(TTPixelTrack pixFitTk,
			  bool bPrintPixInfo = false);
  
  void FillCandPixelHitHistos(TTPixelTrack pixFitTk,
			      bool bPrintPixInfo = false);
  
  void FillSharedPixelHitHistos(int pixTk_Index,
				int tp_Index,
				vector<TVector3> pixHits_XYZ,
				vector<int> pixHits_TTPixelTrackIndex,  
				bool bPrintPixInfo = false);

  void FillResidualHistos(int tp_Index,
			  TTPixelTrack pixFitTk);
    
  void FillHistoBinWithRMS(TH1D *hToFill,
			   TH1D *hToExtractRMS,
			   int binNumber);

  void FitAndFillHistoBin(TH1D *hToFill,
			  TH1D *hToFit,
			  int binNumber,
			  double fitRangeLow,
			  double fitRangeHigh,
			  double redChiSqMax,
			  std::vector<double> v_fitArea);

  void FitAndFillHistoBinSL(TH1D *hToFill,
			    TH1D *hToFit,
			    int binNumber,
			    double significanceLevel=0.05);

 
  std::vector<double> GetFitAreaVector(const double minTotalArea=0.50);
    
  void PrintSettings(void);

  void PrintTPCuts(void);
  
  void PrintPixTrackCuts(void);

  void PrintEfficiencies(void);

  void PrintResolutions(void);

  void MakeEfficiencyHisto(TH1D* h_eff, 
			   TH1D* h_match_tp, 
			   TH1D* h_tp);

  int GetTkMatchingGenP(const int tk_index, 
			const int pdgId);

  int GetPixTkMatchingGenP(const int tk_index, 
			   const int pdgId);
  
  double GetHistoMeanInRange(TH1D* histo, 
			     const double xMin, 
			     const double xMax);
  
  bool IsWithinTrackerAcceptance(double pt, 
				 double eta, 
				 double z0,
				 double dxy);

  bool IsWithinEtaRegion(string etaRegion, 
			 double eta);

  bool IsWithinPtRange(string ptRange, 
		       double pt);


  // Objects/Variables
  AuxTools auxTools;
  HistoTools histoTools;
  TrackingParticles *TPs;
  
  // Options
  bool bVerbose;
  bool bPrintResolutions;
  bool bPrintResolutions_Incl;
  bool bPrintResidualFitInfo;
  bool bSaveResidualFitInfo;
  
  // TP Cuts
  double TP_PT_MIN;
  double TP_ETA_MAX;
  double TP_Z0_MAX;
  double TP_DXY_MAX;
  int TP_NMATCH;
 
  // TTPixelTrack Cuts
  int TK_INDEX_MIN;
  int TK_PIXHITS_MIN;
  int TK_PIXHITS_MAX;
  double TK_PT_MIN;
  double TK_PT_MAX;
  double TK_ABSETA_MIN;
  double TK_ABSETA_MAX;
  int TK_CANDPIXHITS_MIN;
  int TK_CANDPIXHITS_MAX;
  std::vector<int> TK_PIXHITS_TYPE;
  std::vector<int> TK_PIXHITS_PATTERNS;
  
  // Settings
  double _pt_L;
  double _pt_H;
  double _eta_C;
  double _eta_F;
  double _chiSq_overflow;
  double _redChiSq_overflow;
  double _fromCentimetresToMicrometers;
  double _fitSignificanceLevel;  // Gaussian Fits

  // Histogram Binning
  double PT_MAX;
  double PT_BINS;
  double PT_BINWIDTH;
  double ETA_MAX;
  int ETA_BINS;
  double ETA_BINWIDTH;

  // Histograms
  TH1D* h_tp_pt;
  TH1D* h_tp_pt_L;
  TH1D* h_tp_pt_C;
  TH1D* h_tp_pt_I;
  TH1D* h_tp_pt_F;
  TH1D* h_tp_eta;
  TH1D* h_tp_eta_L;
  TH1D* h_tp_eta_M;
  TH1D* h_tp_eta_H;
  TH1D* h_tp_phi;
  TH1D* h_tp_z0;
  TH1D* h_tp_d0;
  TH1D* h_tp_dxy;

  TH1D* h_match_tp_pt;
  TH1D* h_match_tp_pt_L;
  TH1D* h_match_tp_pt_C;
  TH1D* h_match_tp_pt_I;
  TH1D* h_match_tp_pt_F;
  TH1D* h_match_tp_eta;
  TH1D* h_match_tp_eta_L;
  TH1D* h_match_tp_eta_M;
  TH1D* h_match_tp_eta_H;
  TH1D* h_match_tp_phi;
  TH1D* h_match_tp_z0;
  TH1D* h_match_tp_d0;

  TH1D* h_match_trk_nstub;
  TH1D* h_match_trk_nstub_C;
  TH1D* h_match_trk_nstub_I;
  TH1D* h_match_trk_nstub_F;

  // chi2 histograms (last bin is an overflow bin)
  TH1D* h_match_trk_chi2;
  TH1D* h_match_trk_chi2_L;
  TH1D* h_match_trk_chi2_M;
  TH1D* h_match_trk_chi2_H;
  TH1D* h_match_trk_chi2_C_L;
  TH1D* h_match_trk_chi2_I_L;
  TH1D* h_match_trk_chi2_F_L;
  TH1D* h_match_trk_chi2_C_M;
  TH1D* h_match_trk_chi2_I_M;
  TH1D* h_match_trk_chi2_F_M;
  TH1D* h_match_trk_chi2_C_H;
  TH1D* h_match_trk_chi2_I_H;
  TH1D* h_match_trk_chi2_F_H;

  // chi2/dof histograms (lastbin is an overflow bin)
  TH1D* h_match_trk_chi2_dof;
  TH1D* h_match_trk_chi2_dof_L;
  TH1D* h_match_trk_chi2_dof_M;
  TH1D* h_match_trk_chi2_dof_H;
  TH1D* h_match_trk_chi2_dof_C_L;
  TH1D* h_match_trk_chi2_dof_I_L;
  TH1D* h_match_trk_chi2_dof_F_L;
  TH1D* h_match_trk_chi2_dof_C_M;
  TH1D* h_match_trk_chi2_dof_I_M;
  TH1D* h_match_trk_chi2_dof_F_M;
  TH1D* h_match_trk_chi2_dof_C_H;
  TH1D* h_match_trk_chi2_dof_I_H;
  TH1D* h_match_trk_chi2_dof_F_H;

  // Efficiency histograms
  TH1D* h_efficiency_pt;
  TH1D* h_efficiency_pt_L;
  TH1D* h_efficiency_pt_C;
  TH1D* h_efficiency_pt_I;
  TH1D* h_efficiency_pt_F;
  TH1D* h_efficiency_eta;
  TH1D* h_efficiency_eta_L;
  TH1D* h_efficiency_eta_M;
  TH1D* h_efficiency_eta_H;
  TH1D* h_efficiency_phi;
  TH1D* h_efficiency_z0;
  TH1D* h_efficiency_d0;

  // Resolution histograms
  TH1D* h_residual_pt;
  TH1D* h_residual_pt_C;
  TH1D* h_residual_pt_I;
  TH1D* h_residual_pt_F;
  TH1D* h_residual_ptRel;
  TH1D* h_residual_ptRel_C;
  TH1D* h_residual_ptRel_I;
  TH1D* h_residual_ptRel_F;
  TH1D* h_residual_eta;
  TH1D* h_residual_eta_C;
  TH1D* h_residual_eta_I;
  TH1D* h_residual_eta_F;
  TH1D* h_residual_phi ;
  TH1D* h_residual_phi_C;
  TH1D* h_residual_phi_I;
  TH1D* h_residual_phi_F;
  TH1D* h_residual_z0  ;
  TH1D* h_residual_z0_C;
  TH1D* h_residual_z0_I;
  TH1D* h_residual_z0_F;
  TH1D* h_residual_d0  ;
  TH1D* h_residual_d0_C;
  TH1D* h_residual_d0_I;
  TH1D* h_residual_d0_F;

  // Residual Vs pT histograms
  TH1D* h_residualVsPt_pt[20];
  TH1D* h_residualVsPt_pt_C[20];
  TH1D* h_residualVsPt_pt_I[20];
  TH1D* h_residualVsPt_pt_F[20];

  TH1D* h_residualVsPt_ptRel[20];
  TH1D* h_residualVsPt_ptRel_C[20];
  TH1D* h_residualVsPt_ptRel_I[20];
  TH1D* h_residualVsPt_ptRel_F[20];

  TH1D* h_residualVsPt_eta[20];
  TH1D* h_residualVsPt_eta_C[20];
  TH1D* h_residualVsPt_eta_I[20];
  TH1D* h_residualVsPt_eta_F[20];

  TH1D* h_residualVsPt_phi[20];
  TH1D* h_residualVsPt_phi_C[20];
  TH1D* h_residualVsPt_phi_I[20];
  TH1D* h_residualVsPt_phi_F[20];

  TH1D* h_residualVsPt_z0[20];
  TH1D* h_residualVsPt_z0_C[20];
  TH1D* h_residualVsPt_z0_I[20];
  TH1D* h_residualVsPt_z0_F[20];

  TH1D* h_residualVsPt_d0[20];
  TH1D* h_residualVsPt_d0_C[20];
  TH1D* h_residualVsPt_d0_I[20];
  TH1D* h_residualVsPt_d0_F[20];

  // Residual Vs eta histograms
  TH1D* h_residualVsEta_pt[25];
  TH1D* h_residualVsEta_pt_L[25];
  TH1D* h_residualVsEta_pt_M[25];
  TH1D* h_residualVsEta_pt_H[25];

  TH1D* h_residualVsEta_ptRel[25];
  TH1D* h_residualVsEta_ptRel_L[25];
  TH1D* h_residualVsEta_ptRel_M[25];
  TH1D* h_residualVsEta_ptRel_H[25];

  TH1D* h_residualVsEta_eta[25];
  TH1D* h_residualVsEta_eta_L[25];
  TH1D* h_residualVsEta_eta_M[25];
  TH1D* h_residualVsEta_eta_H[25];

  TH1D* h_residualVsEta_phi[25];
  TH1D* h_residualVsEta_phi_L[25];
  TH1D* h_residualVsEta_phi_M[25];
  TH1D* h_residualVsEta_phi_H[25];

  TH1D* h_residualVsEta_z0[25];
  TH1D* h_residualVsEta_z0_L[25];
  TH1D* h_residualVsEta_z0_M[25];
  TH1D* h_residualVsEta_z0_H[25];

  TH1D* h_residualVsEta_d0[25];
  TH1D* h_residualVsEta_d0_L[25];
  TH1D* h_residualVsEta_d0_M[25];
  TH1D* h_residualVsEta_d0_H[25];

  // 2D histograms
  TH2D* h_2d_logchi2_eta;
  TH2D* h_2d_logchi2_dof_eta;
  TH2D* h_2d_dz0_eta;
  TH2D* h_2d_dd0_eta;
  TH2D* h_2d_deta_eta;
  TH2D* h_2d_dphi_eta;
  TH2D* h_2d_dpt_eta; 
  TH2D* h_2d_dptRel_eta; 

  // Resolution Vs pT histograms
  TH1D* h_resolutionVsPt_pt;
  TH1D* h_resolutionVsPt_pt_C;
  TH1D* h_resolutionVsPt_pt_I;
  TH1D* h_resolutionVsPt_pt_F;

  TH1D* h_resolutionVsPt_ptRel;
  TH1D* h_resolutionVsPt_ptRel_C;
  TH1D* h_resolutionVsPt_ptRel_I;
  TH1D* h_resolutionVsPt_ptRel_F;

  TH1D* h_resolutionVsPt_eta;
  TH1D* h_resolutionVsPt_eta_C;
  TH1D* h_resolutionVsPt_eta_I;
  TH1D* h_resolutionVsPt_eta_F;

  TH1D* h_resolutionVsPt_phi;
  TH1D* h_resolutionVsPt_phi_C;
  TH1D* h_resolutionVsPt_phi_I;
  TH1D* h_resolutionVsPt_phi_F;

  TH1D* h_resolutionVsPt_z0;
  TH1D* h_resolutionVsPt_z0_C;
  TH1D* h_resolutionVsPt_z0_I; 
  TH1D* h_resolutionVsPt_z0_F;

  TH1D* h_resolutionVsPt_d0;
  TH1D* h_resolutionVsPt_d0_C;
  TH1D* h_resolutionVsPt_d0_I; 
  TH1D* h_resolutionVsPt_d0_F;

  // Resolution Vs eta histograms
  TH1D* h_resolutionVsEta_pt;
  TH1D* h_resolutionVsEta_pt_L;
  TH1D* h_resolutionVsEta_pt_M;
  TH1D* h_resolutionVsEta_pt_H;

  TH1D* h_resolutionVsEta_ptRel;
  TH1D* h_resolutionVsEta_ptRel_L;
  TH1D* h_resolutionVsEta_ptRel_M;
  TH1D* h_resolutionVsEta_ptRel_H;

  TH1D* h_resolutionVsEta_eta;
  TH1D* h_resolutionVsEta_eta_L;
  TH1D* h_resolutionVsEta_eta_M;
  TH1D* h_resolutionVsEta_eta_H;

  TH1D* h_resolutionVsEta_phi;
  TH1D* h_resolutionVsEta_phi_L;
  TH1D* h_resolutionVsEta_phi_M;
  TH1D* h_resolutionVsEta_phi_H;

  TH1D* h_resolutionVsEta_z0;
  TH1D* h_resolutionVsEta_z0_L;
  TH1D* h_resolutionVsEta_z0_M;
  TH1D* h_resolutionVsEta_z0_H;

  TH1D* h_resolutionVsEta_d0;
  TH1D* h_resolutionVsEta_d0_L;
  TH1D* h_resolutionVsEta_d0_M;
  TH1D* h_resolutionVsEta_d0_H;

  // Pixel Hits
  TH1D* h_pixTk_pixHits_N;
  TH1D* h_pixTk_pixHits_Rho;
  TH1D* h_pixTk_pixHits_Z;
  TH1D* h_pixTk_pixHits_Type;
  TH1D* h_pixTk_pixHits_Pattern;
  TH2D* h_pixTk_pixHits_EtaVsN;
  TH2D* h_pixTk_pixHits_ZVsRho;
  TH2D* h_pixTk_pixHits_ZVsRho_C;
  TH2D* h_pixTk_pixHits_ZVsRho_I;
  TH2D* h_pixTk_pixHits_ZVsRho_F;
  TH2D* h_pixTk_pixHits_ZVsRho_EtaGE2;
  TH2D* h_pixTk_pixHits_XVsY;
  TH2D* h_pixTk_pixHits_XVsY_C;
  TH2D* h_pixTk_pixHits_XVsY_I;
  TH2D* h_pixTk_pixHits_XVsY_F;
  TH2D* h_pixTk_pixHits_XVsY_EtaGE2;
  TH2D* h_pixTk_pixHits_PatternVsEta;
  
  // Shared Pixel Hits
  TH1D* h_pixTk_sharedPixHits_N;
  TH1D* h_pixTk_sharedPixHits_Rho;
  TH1D* h_pixTk_sharedPixHits_Z;
  TH1D* h_pixTk_sharedPixHits_Type;
  TH2D* h_pixTk_sharedPixHits_EtaVsN;
  TH2D* h_pixTk_sharedPixHits_ZVsRho;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_C;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_I;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_F;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_EtaGE2;
  TH2D* h_pixTk_sharedPixHits_XVsY;
  TH2D* h_pixTk_sharedPixHits_XVsY_C;
  TH2D* h_pixTk_sharedPixHits_XVsY_I;
  TH2D* h_pixTk_sharedPixHits_XVsY_F;
  TH2D* h_pixTk_sharedPixHits_XVsY_EtaGE2;
  
  // Candidate Pixel Hits
  TH1D* h_pixTk_candPixHits_N;
  TH1D* h_pixTk_candPixHits_Rho;
  TH1D* h_pixTk_candPixHits_Z;
  TH1D* h_pixTk_candPixHits_Type;
  TH2D* h_pixTk_candPixHits_EtaVsN;
  TH2D* h_pixTk_candPixHits_ZVsRho;
  TH2D* h_pixTk_candPixHits_ZVsRho_C;
  TH2D* h_pixTk_candPixHits_ZVsRho_I;
  TH2D* h_pixTk_candPixHits_ZVsRho_F;
  TH2D* h_pixTk_candPixHits_ZVsRho_EtaGE2;
  TH2D* h_pixTk_candPixHits_XVsY;
  TH2D* h_pixTk_candPixHits_XVsY_C;
  TH2D* h_pixTk_candPixHits_XVsY_I;
  TH2D* h_pixTk_candPixHits_XVsY_F;
  TH2D* h_pixTk_candPixHits_XVsY_EtaGE2;


};

#endif // PixelRefitting_h

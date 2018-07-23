#ifndef Tracking_h
#define Tracking_h

#include <iostream>

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/MCTools.C"
#include "../utilities/AuxTools.C"
#include "../utilities/HistoTools.C"
#include "../utilities/L1Tracks.C"
#include "../utilities/TrackingParticles.C"

class Tracking : public TreeAnalyserMC{
 public:
 Tracking(const string SamplePath,
	    const std::string SampleName,
	    const std::string text_, 
	    const int maxEvents_ = -1, 
	    TTree* tree=0) : 
  TreeAnalyserMC("Tracking", SamplePath, SampleName, text_, maxEvents_, tree) {};
  virtual void Loop();
  
 private:
  void InitVars_(void);

  void BookHistos_(void);

  void FinaliseHistos_(void);

  void FillHistoBinWithRMS_(TH1D *hToFill,
			    TH1D *hToExtractRMS,
			    int binNumber);

  void FitAndFillHistoBin_(TH1D *hToFill,
			   TH1D *hToFit,
			   int binNumber,
			   double fitRangeLow,
			   double fitRangeHigh,
			   double redChiSqMax,
			   std::vector<double> v_fitArea);

  void FitAndFillHistoBinSL_(TH1D *hToFill,
			     TH1D *hToFit,
			     int binNumber,
			     double significanceLevel=0.05);

  std::vector<double> GetFitAreaVector(const double minTotalArea=0.50);
    
  void PrintSettings_(void);

  void PrintPixelTrackCuts_(void);

  void PrintEfficiencies_(void);

  void PrintResolutions_(void);

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
				 double z0);

  bool IsWithinEtaRegion(string etaRegion, 
			 double eta);

  bool IsWithinPtRange(string ptRange, 
			 double pt);


  // Objects/Variables
  AuxTools auxTools_;
  HistoTools histoTools_;
  TrackingParticles *TPs;
 
  // Counters for integrated efficiencies
  bool bDoPixelTracks;
  bool bVerbose;
  bool bPrintEfficiencies;
  bool bPrintResolutions;
  bool bPrintFitInfo;
  bool bSaveFitInfo;
  int pixTk_NPixHits_Min;
  int pixTk_NPixHits_Max;
  std::vector<int> pixTk_MustHit_Type;
  double tp_dxy_Max;
  double _pt_acceptance;
  double _pt_acceptanceMax;
  double _eta_acceptance;
  double _z0_acceptance;
  double _pt_L;
  double _pt_H;
  double _eta_C;
  double _eta_F;
  double _pt_tkMin;  
  int _nStubs_tkMin;
  double _chiSq_overflow;
  double _redChiSq_overflow;
  double _fromCentimetresToMicrometers;
  string s_eta_all;
  string s_eta_C;
  string s_eta_I;
  string s_eta_F;
  string s_pt_tkMin;
  string s_pt_L;
  string s_pt_M;
  string s_pt_H;

  double PT_MAX;
  double PT_BINS;
  double PT_BINWIDTH;
  double ETA_MAX;
  int ETA_BINS;
  double ETA_BINWIDTH;

  int n_all_eta_F;
  int n_all_eta_I;
  int n_all_eta_C;
  int n_match_eta_F;
  int n_match_eta_I;
  int n_match_eta_C;
  
  int n_all_eta_F_pt_L;
  int n_all_eta_I_pt_L;
  int n_all_eta_C_pt_L;
  int n_match_eta_F_pt_L;
  int n_match_eta_I_pt_L;
  int n_match_eta_C_pt_L;

  int n_all_eta_F_pt_M;
  int n_all_eta_I_pt_M;
  int n_all_eta_C_pt_M;
  int n_match_eta_F_pt_M;
  int n_match_eta_I_pt_M;
  int n_match_eta_C_pt_M;

  int n_all_eta_F_pt_H;
  int n_all_eta_I_pt_H;
  int n_all_eta_C_pt_H;
  int n_match_eta_F_pt_H;
  int n_match_eta_I_pt_H;
  int n_match_eta_C_pt_H;
 
  // Event Pixel Hits
  vector<TVector3> pixHits_XYZ;
  vector<int> pixHits_TTPixelTrackIndex;

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

  // efficiency histograms
  TH1D* h_eff_pt;
  TH1D* h_eff_pt_L;
  TH1D* h_eff_pt_C;
  TH1D* h_eff_pt_I;
  TH1D* h_eff_pt_F;
  TH1D* h_eff_eta;
  TH1D* h_eff_eta_L;
  TH1D* h_eff_eta_M;
  TH1D* h_eff_eta_H;
  TH1D* h_eff_phi;
  TH1D* h_eff_z0;
  TH1D* h_eff_d0;

  // resolution histograms
  TH1D* h_res_pt;
  TH1D* h_res_pt_C;
  TH1D* h_res_pt_I;
  TH1D* h_res_pt_F;
  TH1D* h_res_ptRel;
  TH1D* h_res_ptRel_C;
  TH1D* h_res_ptRel_I;
  TH1D* h_res_ptRel_F;
  TH1D* h_res_eta;
  TH1D* h_res_eta_C;
  TH1D* h_res_eta_I;
  TH1D* h_res_eta_F;
  TH1D* h_res_phi ;
  TH1D* h_res_phi_C;
  TH1D* h_res_phi_I;
  TH1D* h_res_phi_F;
  TH1D* h_res_z0  ;
  TH1D* h_res_z0_C;
  TH1D* h_res_z0_I;
  TH1D* h_res_z0_F;
  TH1D* h_res_d0  ;
  TH1D* h_res_d0_C;
  TH1D* h_res_d0_I;
  TH1D* h_res_d0_F;

  // resolution vs. pt histograms
  TH1D* h_resVsPt_pt[20];
  TH1D* h_resVsPt_pt_C[20];
  TH1D* h_resVsPt_pt_I[20];
  TH1D* h_resVsPt_pt_F[20];

  TH1D* h_resVsPt_ptRel[20];
  TH1D* h_resVsPt_ptRel_C[20];
  TH1D* h_resVsPt_ptRel_I[20];
  TH1D* h_resVsPt_ptRel_F[20];

  TH1D* h_resVsPt_eta[20];
  TH1D* h_resVsPt_eta_C[20];
  TH1D* h_resVsPt_eta_I[20];
  TH1D* h_resVsPt_eta_F[20];

  TH1D* h_resVsPt_phi[20];
  TH1D* h_resVsPt_phi_C[20];
  TH1D* h_resVsPt_phi_I[20];
  TH1D* h_resVsPt_phi_F[20];

  TH1D* h_resVsPt_z0[20];
  TH1D* h_resVsPt_z0_C[20];
  TH1D* h_resVsPt_z0_I[20];
  TH1D* h_resVsPt_z0_F[20];

  TH1D* h_resVsPt_d0[20];
  TH1D* h_resVsPt_d0_C[20];
  TH1D* h_resVsPt_d0_I[20];
  TH1D* h_resVsPt_d0_F[20];

  // resolution vs. eta histograms
  TH1D* h_resVsEta_pt[25];
  TH1D* h_resVsEta_pt_L[25];
  TH1D* h_resVsEta_pt_M[25];
  TH1D* h_resVsEta_pt_H[25];

  TH1D* h_resVsEta_ptRel[25];
  TH1D* h_resVsEta_ptRel_L[25];
  TH1D* h_resVsEta_ptRel_M[25];
  TH1D* h_resVsEta_ptRel_H[25];

  TH1D* h_resVsEta_eta[25];
  TH1D* h_resVsEta_eta_L[25];
  TH1D* h_resVsEta_eta_M[25];
  TH1D* h_resVsEta_eta_H[25];

  TH1D* h_resVsEta_phi[25];
  TH1D* h_resVsEta_phi_L[25];
  TH1D* h_resVsEta_phi_M[25];
  TH1D* h_resVsEta_phi_H[25];

  TH1D* h_resVsEta_z0[25];
  TH1D* h_resVsEta_z0_L[25];
  TH1D* h_resVsEta_z0_M[25];
  TH1D* h_resVsEta_z0_H[25];

  TH1D* h_resVsEta_d0[25];
  TH1D* h_resVsEta_d0_L[25];
  TH1D* h_resVsEta_d0_M[25];
  TH1D* h_resVsEta_d0_H[25];

  // 2D histograms
  TH2D* h_2d_logchi2_eta;
  TH2D* h_2d_logchi2_dof_eta;
  TH2D* h_2d_dz0_eta;
  TH2D* h_2d_dd0_eta;
  TH2D* h_2d_deta_eta;
  TH2D* h_2d_dphi_eta;
  TH2D* h_2d_dpt_eta; 
  TH2D* h_2d_dptRel_eta; 

  // resolution vs. pT histograms
  TH1D* h2_resVsPt_pt;
  TH1D* h2_resVsPt_pt_C;
  TH1D* h2_resVsPt_pt_I;
  TH1D* h2_resVsPt_pt_F;

  TH1D* h2_resVsPt_ptRel;
  TH1D* h2_resVsPt_ptRel_C;
  TH1D* h2_resVsPt_ptRel_I;
  TH1D* h2_resVsPt_ptRel_F;

  TH1D* h2_resVsPt_eta;
  TH1D* h2_resVsPt_eta_C;
  TH1D* h2_resVsPt_eta_I;
  TH1D* h2_resVsPt_eta_F;

  TH1D* h2_resVsPt_phi;
  TH1D* h2_resVsPt_phi_C;
  TH1D* h2_resVsPt_phi_I;
  TH1D* h2_resVsPt_phi_F;

  TH1D* h2_resVsPt_z0;
  TH1D* h2_resVsPt_z0_C;
  TH1D* h2_resVsPt_z0_I; 
  TH1D* h2_resVsPt_z0_F;

  TH1D* h2_resVsPt_d0;
  TH1D* h2_resVsPt_d0_C;
  TH1D* h2_resVsPt_d0_I; 
  TH1D* h2_resVsPt_d0_F;

  // resolution vs. eta histograms
  TH1D* h2_resVsEta_pt;
  TH1D* h2_resVsEta_pt_L;
  TH1D* h2_resVsEta_pt_M;
  TH1D* h2_resVsEta_pt_H;

  TH1D* h2_resVsEta_ptRel;
  TH1D* h2_resVsEta_ptRel_L;
  TH1D* h2_resVsEta_ptRel_M;
  TH1D* h2_resVsEta_ptRel_H;

  TH1D* h2_resVsEta_eta;
  TH1D* h2_resVsEta_eta_L;
  TH1D* h2_resVsEta_eta_M;
  TH1D* h2_resVsEta_eta_H;

  TH1D* h2_resVsEta_phi;
  TH1D* h2_resVsEta_phi_L;
  TH1D* h2_resVsEta_phi_M;
  TH1D* h2_resVsEta_phi_H;

  TH1D* h2_resVsEta_z0;
  TH1D* h2_resVsEta_z0_L;
  TH1D* h2_resVsEta_z0_M;
  TH1D* h2_resVsEta_z0_H;

  TH1D* h2_resVsEta_d0;
  TH1D* h2_resVsEta_d0_L;
  TH1D* h2_resVsEta_d0_M;
  TH1D* h2_resVsEta_d0_H;

  // Pixel Hits
  TH1D* h_pixTk_pixHits_N;
  TH1D* h_pixTk_pixHits_Rho;
  TH1D* h_pixTk_pixHits_Z;
  TH1D* h_pixTk_pixHits_Type;
  TH2D* h_pixTk_pixHits_ZVsRho;
  TH2D* h_pixTk_pixHits_ZVsRho_C;
  TH2D* h_pixTk_pixHits_ZVsRho_I;
  TH2D* h_pixTk_pixHits_ZVsRho_F;
  TH2D* h_pixTk_pixHits_ZVsRho_EtaGE2;

  // Shared Pixel Hits
  TH1D* h_pixTk_sharedPixHits_N;
  TH1D* h_pixTk_sharedPixHits_Rho;
  TH1D* h_pixTk_sharedPixHits_Z;
  TH1D* h_pixTk_sharedPixHits_Type;
  TH2D* h_pixTk_sharedPixHits_ZVsRho;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_C;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_I;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_F;
  TH2D* h_pixTk_sharedPixHits_ZVsRho_EtaGE2;
  
  // Candidate Pixel Hits
  TH1D* h_pixTk_candPixHits_N;
  TH1D* h_pixTk_candPixHits_Rho;
  TH1D* h_pixTk_candPixHits_Z;
  TH1D* h_pixTk_candPixHits_Type;
  TH2D* h_pixTk_candPixHits_ZVsRho;
  TH2D* h_pixTk_candPixHits_ZVsRho_C;
  TH2D* h_pixTk_candPixHits_ZVsRho_I;
  TH2D* h_pixTk_candPixHits_ZVsRho_F;
  TH2D* h_pixTk_candPixHits_ZVsRho_EtaGE2;


};

#endif // Tracking_h

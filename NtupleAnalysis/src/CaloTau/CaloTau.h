#ifndef CaloTau_h
#define CaloTau_h

// System 
#include <iostream>
#include <stdlib.h> 

// User
#include "../Framework/src/TreeAnalyserMC.C"
#include "../Framework/src/TreeReaderMC.C"


#include "../Auxiliary/src/AuxTools.C"
#include "../Auxiliary/src/Table.C"
#include "../Auxiliary/src/MCTools.C"
#include "../Auxiliary/src/HistoTools.C"
#include "../Auxiliary/src/Datasets.C" 

#include "../DataFormat/src/L1TkTauParticle.C"
#include "../DataFormat/src/GenParticle.C"
#include "../DataFormat/src/L1Jet.C"
#include "../DataFormat/src/L1Tau.C"
#include "../DataFormat/src/L1Sum.C"
#include "../DataFormat/src/L1CaloTP.C"

// ROOT
//#include "TEfficiency.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

//marina
#include <vector>
#include <iostream>

using namespace std;

class CaloTau : public TreeAnalyserMC{
 public:
  // Constructors/Destructors
  ~CaloTau(){};
 CaloTau(const string SamplePath,
	const string SampleName,
	const string text_, 
	const int maxEvents_ = -1, 
	TChain* chain=0) :

  TreeAnalyserMC("", SamplePath, SampleName, text_, maxEvents_, chain)
    { 
      auxTools_.StopwatchStart();
      mcSample = SampleName;
      InitObjects();
    };
  
  // Public Variables
  virtual void Loop();

  void PrintSettings(void);

  void GetMatchingGenParticle(L1TkTauParticle &L1TkTau,
			      vector<GenParticle> hadGenTaus);			    

  void GetLdgAndSubldgIndices(vector<L1TkTauParticle> myTaus,
			      int &iLdg,
			      int &iSubldg);

  // Public Variables
  bool DEBUG;
  int realTauMom;
  int nMaxNumOfHTausPossible;
  bool mcMatching_unique;
  double mcMatching_dRMax;
  double dR_match_min;
  double dR_match_max;
  double dR_sigCone_min;
  double dR_sigCone_max;
  double dR_isoCone_min;
  double dR_isoCone_max; // 1 x 1 TT = 5 x 5 crystals = 5 x 0.087 = 0.435 (for central region)
  
  double _eta_C;
  double _eta_F;
  string caloCorrectionType;
  string mcSample;

  struct SortAscendingAbs{ bool operator() (double a, double b) const { return abs(a) > abs(b); } }; 
  struct SortDescendingAbs{ bool operator() (double a, double b) const { return abs(a) < abs(b); } }; 
  
 private:
  // Function declaration
  void BookHistos_(void);
  void WriteHistos_(void);
  void InitObjects(void);
  void InitVars_(void);

  void FillEfficiency_(TH1D *hSignalEfficiency,
		       const double ldgEt);

  void FillEfficiency_(TH2D *hSignalEfficiency,
		       const double ldgEt1,
		       const double ldgEt2);  

  void FillRate_(TH1D *hRate,
		 const double ldgEt);

  void FillRate_(TH2D *hRate,
		 const double ldgEt1,
		 const double ldgEt2);  

  void FinaliseEffHisto_(TH1D *histo, 
			 const int nEvtsTotal);

  void FinaliseEffHisto_(TH2D *histo, 
			 const int nEvtsTotal);  

  void FillTurnOn_Numerator_(vector<L1TkTauParticle> L1Taus,
			     const double minEt,
			     TH1D *hTurnOn, 
			     TH1D *hTurnOn_1pr, 
			     TH1D *hTurnOn_3pr, 
			     TH1D *hTurnOn_withNeutrals, 
			     TH1D *hTurnOn_noNeutrals);

  void FillSingleTau_(vector<L1TkTauParticle> L1Taus,
		      TH1D *hRate,
		      TH1D *hEfficiency,
		      double minEta=0.0,
		      double maxEta=999.9);

  void FillDiTau_(vector<L1TkTauParticle> L1Taus, 
		  TH1D *hRate,
		  TH1D *hEfficiency,
		  double minEta=0.0,
		  double maxEta=999.9);

  void FillDiTau_(vector<L1TkTauParticle> L1Taus1,
		  vector<L1TkTauParticle> L1Taus2,
		  TH2D *hRate,
		  TH2D *hEfficiency);

  vector<L1TkTauParticle> GetMcMatchedL1Taus(vector<L1TkTauParticle> L1Taus);

  bool IsWithinEtaRegion(string etaRegion,
			 double eta);

  
  // Variable declaration
  AuxTools auxTools_;
  MCTools mcTools_;
  Datasets datasets_;
  HistoTools histoTools_;
  bool bFoundAllTaus_;

  // GenParticles Histograms
  TH2D* hGenP_VisEt_Vs_dRMaxLdgPion;
  TH2D* hGenP_PtLdg_Vs_dRMaxLdgPion;

  // Counters
  TH1D* hCounters;

  // CaloTaus (http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_3_7/doc/html/d1/dc8/classl1t_1_1CaloTools.html)
  TH1D* hCalo_Multiplicity;
  TH1D* hCalo_Multiplicity_MC;
  TH1D* hCalo_Et; 
  TH1D* hCalo_Eta;
  TH1D* hCalo_Phi;
  TH1D* hCalo_IEt;
  TH1D* hCalo_IEta; // ieta of seed tower
  TH1D* hCalo_IPhi; // iphi of seed tower
  TH1D* hCalo_Iso;
  TH1D* hCalo_TowerIEta;
  TH1D* hCalo_TowerIPhi;
  TH1D* hCalo_RawEt; // raw (uncalibrated) cluster sum
  TH1D* hCalo_IsoEt; // raw isolation sum - cluster sum
  TH1D* hCalo_NTT;   // n towers above threshold
  TH1D* hCalo_HasEM;
  TH1D* hCalo_IsMerged;
  TH1D* hCalo_DeltaRGenP;
  //
  TH1D* hCaloIso_Multiplicity;
  TH1D* hCaloIso_Multiplicity_MC;
  TH1D* hCaloIso_Et; 
  TH1D* hCaloIso_Eta;
  TH1D* hCaloIso_Phi;
  TH1D* hCaloIso_IEt;
  TH1D* hCaloIso_IEta; // ieta of seed tower
  TH1D* hCaloIso_IPhi; // iphi of seed tower
  TH1D* hCaloIso_Iso;
  TH1D* hCaloIso_TowerIEta;
  TH1D* hCaloIso_TowerIPhi;
  TH1D* hCaloIso_RawEt; // raw (uncalibrated) cluster sum
  TH1D* hCaloIso_IsoEt; // raw isolation sum - cluster sum
  TH1D* hCaloIso_NTT;   // n towers above threshold
  TH1D* hCaloIso_HasEM;
  TH1D* hCaloIso_IsMerged;
  TH1D* hCaloIso_DeltaRGenP;

  // Resolutions
  TH1D* hCalo_ResolutionEt;
  TH1D* hCalo_ResolutionEt_1pr;
  TH1D* hCalo_ResolutionEt_3pr;
  TH1D* hCalo_ResolutionEt_withNeutrals;
  TH1D* hCalo_ResolutionEt_noNeutrals;
  TH1D* hCalo_ResolutionEta;
  TH1D* hCalo_ResolutionEta_1pr;
  TH1D* hCalo_ResolutionEta_3pr;
  TH1D* hCalo_ResolutionEta_withNeutrals;
  TH1D* hCalo_ResolutionEta_noNeutrals;
  TH1D* hCalo_ResolutionPhi;
  TH1D* hCalo_ResolutionPhi_1pr;
  TH1D* hCalo_ResolutionPhi_3pr;
  TH1D* hCalo_ResolutionPhi_withNeutrals;
  TH1D* hCalo_ResolutionPhi_noNeutrals;
  TH1D* hCalo_ResolutionEt_C;
  TH1D* hCalo_ResolutionEta_C;
  TH1D* hCalo_ResolutionPhi_C;
  TH1D* hCalo_ResolutionEt_I;
  TH1D* hCalo_ResolutionEta_I;
  TH1D* hCalo_ResolutionPhi_I;
  TH1D* hCalo_ResolutionEt_F;
  TH1D* hCalo_ResolutionEta_F;
  TH1D* hCalo_ResolutionPhi_F;
  //
  TH1D* hCaloIso_ResolutionEt;
  TH1D* hCaloIso_ResolutionEt_1pr;
  TH1D* hCaloIso_ResolutionEt_3pr;
  TH1D* hCaloIso_ResolutionEt_withNeutrals;
  TH1D* hCaloIso_ResolutionEt_noNeutrals;
  TH1D* hCaloIso_ResolutionEta;
  TH1D* hCaloIso_ResolutionEta_1pr;
  TH1D* hCaloIso_ResolutionEta_3pr;
  TH1D* hCaloIso_ResolutionEta_withNeutrals;
  TH1D* hCaloIso_ResolutionEta_noNeutrals;
  TH1D* hCaloIso_ResolutionPhi;
  TH1D* hCaloIso_ResolutionPhi_1pr;
  TH1D* hCaloIso_ResolutionPhi_3pr;
  TH1D* hCaloIso_ResolutionPhi_withNeutrals;
  TH1D* hCaloIso_ResolutionPhi_noNeutrals;
  TH1D* hCaloIso_ResolutionEt_C;
  TH1D* hCaloIso_ResolutionEta_C;
  TH1D* hCaloIso_ResolutionPhi_C;
  TH1D* hCaloIso_ResolutionEt_I;
  TH1D* hCaloIso_ResolutionEta_I;
  TH1D* hCaloIso_ResolutionPhi_I;
  TH1D* hCaloIso_ResolutionEt_F;
  TH1D* hCaloIso_ResolutionEta_F;
  TH1D* hCaloIso_ResolutionPhi_F;

  // SingleTau
  TH1D* hCalo_Rate;
  TH1D* hCalo_Rate_C;
  TH1D* hCalo_Rate_I;
  TH1D* hCalo_Rate_F;
  TH1D* hCalo_Eff;
  TH1D* hCalo_Eff_C;
  TH1D* hCalo_Eff_I;
  TH1D* hCalo_Eff_F;
  //
  TH1D* hCaloIso_Rate;
  TH1D* hCaloIso_Rate_C;
  TH1D* hCaloIso_Rate_I;
  TH1D* hCaloIso_Rate_F;
  TH1D* hCaloIso_Eff;
  TH1D* hCaloIso_Eff_C;
  TH1D* hCaloIso_Eff_I;
  TH1D* hCaloIso_Eff_F;

  // DiTau
  TH1D* hCalo_Rate_DiTau;
  TH1D* hCalo_Rate_DiTau_C;
  TH1D* hCalo_Rate_DiTau_I;
  TH1D* hCalo_Rate_DiTau_F;
  TH1D* hCalo_Eff_DiTau;
  TH1D* hCalo_Eff_DiTau_C;
  TH1D* hCalo_Eff_DiTau_I;
  TH1D* hCalo_Eff_DiTau_F;
  //
  TH1D* hCaloIso_Rate_DiTau;
  TH1D* hCaloIso_Rate_DiTau_C;
  TH1D* hCaloIso_Rate_DiTau_I;
  TH1D* hCaloIso_Rate_DiTau_F;
  TH1D* hCaloIso_Eff_DiTau;
  TH1D* hCaloIso_Eff_DiTau_C;
  TH1D* hCaloIso_Eff_DiTau_I;
  TH1D* hCaloIso_Eff_DiTau_F;

  // Turn-Ons
  TH1D* hMcHadronicTau_VisEt;
  TH1D* hMcHadronicTau_VisEt_1pr;
  TH1D* hMcHadronicTau_VisEt_3pr;
  TH1D* hMcHadronicTau_VisEt_withNeutrals;
  TH1D* hMcHadronicTau_VisEt_noNeutrals;

  TH1D* hCalo_TurnOn25;
  TH1D* hCalo_TurnOn25_1pr;
  TH1D* hCalo_TurnOn25_3pr;
  TH1D* hCalo_TurnOn25_withNeutrals;
  TH1D* hCalo_TurnOn25_noNeutrals;
  TH1D* hCalo_TurnOn50;
  TH1D* hCalo_TurnOn50_1pr;
  TH1D* hCalo_TurnOn50_3pr;
  TH1D* hCalo_TurnOn50_withNeutrals;
  TH1D* hCalo_TurnOn50_noNeutrals;
  //
  TH1D* hCaloIso_TurnOn25;
  TH1D* hCaloIso_TurnOn25_1pr;
  TH1D* hCaloIso_TurnOn25_3pr;
  TH1D* hCaloIso_TurnOn25_withNeutrals;
  TH1D* hCaloIso_TurnOn25_noNeutrals;
  TH1D* hCaloIso_TurnOn50;
  TH1D* hCaloIso_TurnOn50_1pr;
  TH1D* hCaloIso_TurnOn50_3pr;
  TH1D* hCaloIso_TurnOn50_withNeutrals;
  TH1D* hCaloIso_TurnOn50_noNeutrals;

};

#endif

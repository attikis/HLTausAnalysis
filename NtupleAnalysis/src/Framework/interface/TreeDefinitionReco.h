#ifndef TreeDefinitionReco_h
#define TreeDefinitionReco_h

// User
#include "TreeDefinitionBase.h"
#include <TFriendElement.h>
#include <TList.h>
#include <TString.h>

using namespace std;

class TreeDefinitionReco : public virtual TreeDefinitionBase
{
 public:


  // Event
  UInt_t          run;
  ULong64_t       event;
  UInt_t          lumi;
  UInt_t          bx;
  ULong64_t       orbit;
  ULong64_t       time;
  Int_t           nPV;
  Int_t           nPV_True;
  vector<TString> hlt;
  Double_t        puWeight;


  // L1PhaseII
  UInt_t          EG_N;
  vector<double>  EG_Et;
  vector<double>  EG_Eta;
  vector<double>  EG_Phi;
  vector<int>     EG_Bx;
  vector<double>  EG_Iso;
  vector<double>  EG_zVtx;
  vector<int>     EG_HwQual;
  
  UInt_t           tkEM_N;
  vector<double>   tkEM_Et;
  vector<double>   tkEM_Eta;
  vector<double>   tkEM_Phi;
  vector<int>      tkEM_Bx;
  vector<double>   tkEM_TrkIso;
  vector<double>   tkEM_zVtx;
  vector<double>   tkEM_HwQual;
  vector<double>   tkEM_EGRefPt;
  vector<double>   tkEM_EGRefEta;
  vector<double>   tkEM_EGRefPhi;


  // Calo Towers 
  //==Calo TPs
  Short_t         hcalTP_N;
  vector<short>   hcalTP_ieta;
  vector<short>   hcalTP_iphi;
  vector<short>   hcalTP_Caliphi;
  vector<float>   hcalTP_et;
  vector<short>   hcalTP_compEt;
  vector<short>   hcalTP_fineGrain;

  Short_t         ecalTP_N;
  vector<short>   ecalTP_ieta;
  vector<short>   ecalTP_iphi;
  vector<short>   ecalTP_Caliphi;
  vector<float>   ecalTP_et;
  vector<short>   ecalTP_compEt;
  vector<short>   ecalTP_fineGrain;

  Short_t         ecalEBTP_N;
  vector<short>   ecalEBTP_ieta;
  vector<short>   ecalEBTP_iphi;
  vector<short>   ecalEBTP_Caliphi;
  vector<float>   ecalEBTP_et;
  vector<short>   ecalEBTP_encodEt;
  vector<short>   ecalEBTP_spike;
  vector<short>   ecalEBTP_time;

  //==Calo Towers 
  Short_t         caloTower_N;
  vector<short>   caloTower_ieta;
  vector<short>   caloTower_iphi;
  vector<short>   caloTower_iet;
  vector<short>   caloTower_iem;
  vector<short>   caloTower_ihad;
  vector<short>   caloTower_iratio;
  vector<short>   caloTower_iqual;
  vector<float>   caloTower_et;
  vector<float>   caloTower_eta;
  vector<float>   caloTower_phi;

  // UpgradeTfMuon
  // ****Fix me
 
  // Upgrade
  UShort_t        L1EG_N;
  vector<float>   L1EG_Et;
  vector<float>   L1EG_Eta;
  vector<float>   L1EG_Phi;
  vector<short>   L1EG_IEt;
  vector<short>   L1EG_IEta;
  vector<short>   L1EG_IPhi;
  vector<short>   L1EG_Iso;
  vector<short>   L1EG_Bx;
  vector<short>   L1EG_TowerIPhi;
  vector<short>   L1EG_TowerIEta;
  vector<short>   L1EG_RawEt;
  vector<short>   L1EG_IsoEt;
  vector<short>   L1EG_FootprintEt;
  vector<short>   L1EG_NTT;
  vector<short>   L1EG_Shape;
  vector<short>   L1EG_TowerHoE;

  UShort_t        L1Tau_N;
  vector<float>   L1Tau_Et;
  vector<float>   L1Tau_Eta;
  vector<float>   L1Tau_Phi;
  vector<short>   L1Tau_IEt;
  vector<short>   L1Tau_IEta;
  vector<short>   L1Tau_IPhi;
  vector<short>   L1Tau_Iso;
  vector<short>   L1Tau_Bx;
  vector<short>   L1Tau_TowerIPhi;
  vector<short>   L1Tau_TowerIEta;
  vector<short>   L1Tau_RawEt;
  vector<short>   L1Tau_IsoEt;
  vector<short>   L1Tau_NTT;
  vector<short>   L1Tau_HasEM;
  vector<short>   L1Tau_IsMerged;
  vector<short>   L1Tau_HwQual;

  UShort_t        L1Jet_N;
  vector<float>   L1Jet_Et;
  vector<float>   L1Jet_Eta;
  vector<float>   L1Jet_Phi;
  vector<short>   L1Jet_IEt;
  vector<short>   L1Jet_IEta;
  vector<short>   L1Jet_IPhi;
  vector<short>   L1Jet_Bx;
  vector<short>   L1Jet_TowerIPhi;
  vector<short>   L1Jet_TowerIEta;
  vector<short>   L1Jet_RawEt;
  vector<short>   L1Jet_SeedEt;
  vector<short>   L1Jet_PUEt;
  vector<short>   L1Jet_PUDonutEt0;
  vector<short>   L1Jet_PUDonutEt1;
  vector<short>   L1Jet_PUDonutEt2;
  vector<short>   L1Jet_PUDonutEt3;

  UShort_t        L1Muon_N;
  vector<float>   L1Muon_Et;
  vector<float>   L1Muon_Eta;
  vector<float>   L1Muon_Phi;
  vector<float>   L1Muon_EtaAtVtx;
  vector<float>   L1Muon_PhiAtVtx;
  vector<short>   L1Muon_IEt;
  vector<short>   L1Muon_IEta;
  vector<short>   L1Muon_IPhi;
  vector<short>   L1Muon_IEtaAtVtx;
  vector<short>   L1Muon_IPhiAtVtx;
  vector<short>   L1Muon_IDEta;
  vector<short>   L1Muon_IDPhi;
  vector<short>   L1Muon_Chg;
  vector<unsigned short> L1Muon_Iso;
  vector<unsigned short> L1Muon_Qual;
  vector<unsigned short> L1Muon_TfMuon_Idx;
  vector<short>   L1Muon_Bx;

  UShort_t        L1Sum_N;
  vector<short>   L1Sum_Type;
  vector<float>   L1Sum_Et;
  vector<float>   L1Sum_Phi;
  vector<short>   L1Sum_IEt;
  vector<short>   L1Sum_IPhi;
  vector<float>   L1Sum_Bx;

  // uGT
  Int_t           m_orbitNr;
  Int_t           m_bxNr;
  Int_t           m_bxInEvent;
  Bool_t          m_finalOR;
  Bool_t          m_finalORPreVeto;
  Bool_t          m_finalORVeto;
  Int_t           m_preScColumn;
  vector<bool>    m_algoDecisionInitial;
  vector<bool>    m_algoDecisionPreScaled;
  vector<bool>    m_algoDecisionFinal;

  // HO
  UInt_t          hcalDetId_N;
  UInt_t          hcalQIESample_N;
  vector<int>     hcalDetId_IEta;
  vector<int>     hcalDetId_IPhi;
  vector<int>     hcalQIESample;
  vector<int>     hcalQIESample_Adc;
  vector<int>     hcalQIESample_Dv;
  vector<int>     hcalQIESample_Er;
  
  // Tracks
  vector<float>   *L1Tks_Pt;
  vector<float>   *L1Tks_Eta;
  vector<float>   *L1Tks_Phi;
  vector<float>   *L1Tks_d0;
  vector<float>   *L1Tks_z0;
  vector<float>   *L1Tks_ChiSquared;
  vector<int>     *L1Tks_NStubs;
  vector<float>   *L1Tks_StubPtConsistency;
  vector<float>   *L1Tks_RInv;
  vector<int>     *L1Tks_IsGenuine;
  vector<int>     *L1Tks_IsLoose;
  vector<int>     *L1Tks_IsUnknown;
  vector<int>     *L1Tks_IsCombinatoric;
  vector<int>     *L1Tks_IsFake;

  vector<int>     *L1Tks_TP_PdgId;
  vector<float>   *L1Tks_TP_Pt;
  vector<float>   *L1Tks_TP_Eta;
  vector<float>   *L1Tks_TP_Phi;
  vector<float>   *L1Tks_TP_z0;
  vector<float>   *L1Tks_TP_dxy;

  // Tracking Particles
  vector<float>   *TP_Pt;
  vector<float>   *TP_Eta;
  vector<float>   *TP_Phi;
  vector<float>   *TP_dxy;
  vector<float>   *TP_d0;
  vector<float>   *TP_z0;
  vector<float>   *TP_d0_produced;
  vector<float>   *TP_z0_produced;
  vector<int>     *TP_PdgId;
  vector<int>     *TP_NMatch;
  vector<int>     *TP_NStubs;
  vector<int>     *TP_EventId;
  vector<int>     *TP_Charge;
  
  vector<float>   *TP_Trk_Pt;
  vector<float>   *TP_Trk_Eta;
  vector<float>   *TP_Trk_Phi;
  vector<float>   *TP_Trk_z0;
  vector<float>   *TP_Trk_d0;
  vector<float>   *TP_Trk_ChiSquared;
  vector<int>     *TP_Trk_NStubs;
 
  // =========================================== EMULATOR TREES =============================================

  // Calo Towers Emulator
  //==Calo TPs (Emu)
  Short_t         hcalTPEmu_N;
  vector<short>   hcalTPEmu_ieta;
  vector<short>   hcalTPEmu_iphi;
  vector<short>   hcalTPEmu_Caliphi;
  vector<float>   hcalTPEmu_et;
  vector<short>   hcalTPEmu_compEt;
  vector<short>   hcalTPEmu_fineGrain;

  Short_t         ecalTPEmu_N;
  vector<short>   ecalTPEmu_ieta;
  vector<short>   ecalTPEmu_iphi;
  vector<short>   ecalTPEmu_Caliphi;
  vector<float>   ecalTPEmu_et;
  vector<short>   ecalTPEmu_compEt;
  vector<short>   ecalTPEmu_fineGrain;

  Short_t         ecalEBTPEmu_N;
  vector<short>   ecalEBTPEmu_ieta;
  vector<short>   ecalEBTPEmu_iphi;
  vector<short>   ecalEBTPEmu_Caliphi;
  vector<float>   ecalEBTPEmu_et;
  vector<short>   ecalEBTPEmu_encodEt;
  vector<short>   ecalEBTPEmu_spike;
  vector<short>   ecalEBTPEmu_time;

  //==Calo Towers (Emu)
  Short_t         caloTowerEmu_N;
  vector<short>   caloTowerEmu_ieta;
  vector<short>   caloTowerEmu_iphi;
  vector<short>   caloTowerEmu_iet;
  vector<short>   caloTowerEmu_iem;
  vector<short>   caloTowerEmu_ihad;
  vector<short>   caloTowerEmu_iratio;
  vector<short>   caloTowerEmu_iqual;
  vector<float>   caloTowerEmu_et;
  vector<float>   caloTowerEmu_eta;
  vector<float>   caloTowerEmu_phi;

  //==Calo Clusters (Emu)
  Short_t         caloClusterEmu_N;
  vector<short>   caloClusterEmu_ieta;
  vector<short>   caloClusterEmu_iphi;
  vector<short>   caloClusterEmu_iet;
  vector<short>   caloClusterEmu_iqual;
  vector<float>   caloClusterEmu_et;
  vector<float>   caloClusterEmu_eta;
  vector<float>   caloClusterEmu_phi;

  // Upgrade Emulator
  UShort_t        L1EGEmu_N;
  vector<float>   L1EGEmu_Et;
  vector<float>   L1EGEmu_Eta;
  vector<float>   L1EGEmu_Phi;
  vector<short>   L1EGEmu_IEt;
  vector<short>   L1EGEmu_IEta;
  vector<short>   L1EGEmu_IPhi;
  vector<short>   L1EGEmu_Iso;
  vector<short>   L1EGEmu_Bx;
  vector<short>   L1EGEmu_TowerIPhi;
  vector<short>   L1EGEmu_TowerIEta;
  vector<short>   L1EGEmu_RawEt;
  vector<short>   L1EGEmu_IsoEt;
  vector<short>   L1EGEmu_FootprintEt;
  vector<short>   L1EGEmu_NTT;
  vector<short>   L1EGEmu_Shape;
  vector<short>   L1EGEmu_TowerHoE;

  UShort_t        L1TauEmu_N;
  vector<float>   L1TauEmu_Et;
  vector<float>   L1TauEmu_Eta;
  vector<float>   L1TauEmu_Phi;
  vector<short>   L1TauEmu_IEt;
  vector<short>   L1TauEmu_IEta;
  vector<short>   L1TauEmu_IPhi;
  vector<short>   L1TauEmu_Iso;
  vector<short>   L1TauEmu_Bx;
  vector<short>   L1TauEmu_TowerIPhi;
  vector<short>   L1TauEmu_TowerIEta;
  vector<short>   L1TauEmu_RawEt;
  vector<short>   L1TauEmu_IsoEt;
  vector<short>   L1TauEmu_NTT;
  vector<short>   L1TauEmu_HasEM;
  vector<short>   L1TauEmu_IsMerged;
  vector<short>   L1TauEmu_HwQual;

  UShort_t        L1JetEmu_N;
  vector<float>   L1JetEmu_Et;
  vector<float>   L1JetEmu_Eta;
  vector<float>   L1JetEmu_Phi;
  vector<short>   L1JetEmu_IEt;
  vector<short>   L1JetEmu_IEta;
  vector<short>   L1JetEmu_IPhi;
  vector<short>   L1JetEmu_Bx;
  vector<short>   L1JetEmu_TowerIPhi;
  vector<short>   L1JetEmu_TowerIEta;
  vector<short>   L1JetEmu_RawEt;
  vector<short>   L1JetEmu_SeedEt;
  vector<short>   L1JetEmu_PUEt;
  vector<short>   L1JetEmu_PUDonutEt0;
  vector<short>   L1JetEmu_PUDonutEt1;
  vector<short>   L1JetEmu_PUDonutEt2;
  vector<short>   L1JetEmu_PUDonutEt3;

  UShort_t        L1MuonEmu_N;
  vector<float>   L1MuonEmu_Et;
  vector<float>   L1MuonEmu_Eta;
  vector<float>   L1MuonEmu_Phi;
  vector<float>   L1MuonEmu_EtaAtVtx;
  vector<float>   L1MuonEmu_PhiAtVtx;
  vector<short>   L1MuonEmu_IEt;
  vector<short>   L1MuonEmu_IEta;
  vector<short>   L1MuonEmu_IPhi;
  vector<short>   L1MuonEmu_IEtaAtVtx;
  vector<short>   L1MuonEmu_IPhiAtVtx;
  vector<short>   L1MuonEmu_IDEta;
  vector<short>   L1MuonEmu_IDPhi;
  vector<short>   L1MuonEmu_Chg;
  vector<unsigned short> L1MuonEmu_Iso;
  vector<unsigned short> L1MuonEmu_Qual;
  vector<unsigned short> L1MuonEmu_TfMuon_Idx;
  vector<short>   L1MuonEmu_Bx;

  UShort_t        L1SumEmu_N;
  vector<short>   L1SumEmu_Type;
  vector<float>   L1SumEmu_Et;
  vector<float>   L1SumEmu_Phi;
  vector<short>   L1SumEmu_IEt;
  vector<short>   L1SumEmu_IPhi;
  vector<float>   L1SumEmu_Bx;

  // uGT Emulator
  Int_t           mEmu_orbitNr;
  Int_t           mEmu_bxNr;
  Int_t           mEmu_bxInEvent;
  Bool_t          mEmu_finalOR;
  Bool_t          mEmu_finalORPreVeto;
  Bool_t          mEmu_finalORVeto;
  Int_t           mEmu_preScColumn;
  vector<bool>    mEmu_algoDecisionInitial;
  vector<bool>    mEmu_algoDecisionPreScaled;
  vector<bool>    mEmu_algoDecisionFinal;

 //********* List of branches **********
  

  // Event
  TBranch        *b_Event_run;
  TBranch        *b_Event_event;
  TBranch        *b_Event_lumi; 
  TBranch        *b_Event_bx;  
  TBranch        *b_Event_orbit;  
  TBranch        *b_Event_time;   
  TBranch        *b_Event_nPV;   
  TBranch        *b_Event_nPV_True; 
  TBranch        *b_Event_hlt;  
  TBranch        *b_Event_puWeight;


  //L1PhaseII
  TBranch        *b_EG_N; 
  TBranch        *b_EG_Et; 
  TBranch        *b_EG_Eta;
  TBranch        *b_EG_Phi;
  TBranch        *b_EG_Bx; 
  TBranch        *b_EG_Iso;
  TBranch        *b_EG_zVtx;
  TBranch        *b_EG_HwQual;

  TBranch        *b_tkEM_N;
  TBranch        *b_tkEM_Et;
  TBranch        *b_tkEM_Eta;
  TBranch        *b_tkEM_Phi;
  TBranch        *b_tkEM_Bx;
  TBranch        *b_tkEM_TrkIso;
  TBranch        *b_tkEM_zVtx;
  TBranch        *b_tkEM_HwQual;
  TBranch        *b_tkEM_EGRefPt;
  TBranch        *b_tkEM_EGRefEta;
  TBranch        *b_tkEM_EGRefPhi;


  // Calo Towers
  //==CaloTP
  TBranch        *b_hcalTP_N;  
  TBranch        *b_hcalTP_ieta;
  TBranch        *b_hcalTP_iphi;
  TBranch        *b_hcalTP_Caliphi;
  TBranch        *b_hcalTP_et;   
  TBranch        *b_hcalTP_compEt;
  TBranch        *b_hcalTP_fineGrain; 

  TBranch        *b_ecalTP_N;  
  TBranch        *b_ecalTP_ieta; 
  TBranch        *b_ecalTP_iphi; 
  TBranch        *b_ecalTP_Caliphi; 
  TBranch        *b_ecalTP_et;   
  TBranch        *b_ecalTP_compEt;
  TBranch        *b_ecalTP_fineGrain; 

  TBranch        *b_ecalEBTP_N;  
  TBranch        *b_ecalEBTP_ieta;
  TBranch        *b_ecalEBTP_iphi;
  TBranch        *b_ecalEBTP_Caliphi; 
  TBranch        *b_ecalEBTP_et;   
  TBranch        *b_ecalEBTP_encodEt;
  TBranch        *b_ecalEBTP_spike; 
  TBranch        *b_ecalEBTP_time;  
  
  //==Calo Towers
  TBranch        *b_caloTower_N; 
  TBranch        *b_caloTower_ieta;  
  TBranch        *b_caloTower_iphi;  
  TBranch        *b_caloTower_iet;  
  TBranch        *b_caloTower_iem;  
  TBranch        *b_caloTower_ihad; 
  TBranch        *b_caloTower_iratio;
  TBranch        *b_caloTower_iqual; 
  TBranch        *b_caloTower_et;   
  TBranch        *b_caloTower_eta;  
  TBranch        *b_caloTower_phi;  

  // UpgradeTfMuon
  // ****Fix me

  // Upgrade
  TBranch        *b_L1EG_N;   
  TBranch        *b_L1EG_Et;   
  TBranch        *b_L1EG_Eta;   
  TBranch        *b_L1EG_Phi;   
  TBranch        *b_L1EG_IEt;   
  TBranch        *b_L1EG_IEta;   
  TBranch        *b_L1EG_IPhi;   
  TBranch        *b_L1EG_Iso;   
  TBranch        *b_L1EG_Bx;   
  TBranch        *b_L1EG_TowerIPhi;   
  TBranch        *b_L1EG_TowerIEta;   
  TBranch        *b_L1EG_RawEt;   
  TBranch        *b_L1EG_IsoEt;   
  TBranch        *b_L1EG_FootprintEt;   
  TBranch        *b_L1EG_NTT;   
  TBranch        *b_L1EG_Shape;   
  TBranch        *b_L1EG_TowerHoE;   

  TBranch        *b_L1Tau_N;   
  TBranch        *b_L1Tau_Et;   
  TBranch        *b_L1Tau_Eta;   
  TBranch        *b_L1Tau_Phi;   
  TBranch        *b_L1Tau_IEt;   
  TBranch        *b_L1Tau_IEta;   
  TBranch        *b_L1Tau_IPhi;   
  TBranch        *b_L1Tau_Iso;   
  TBranch        *b_L1Tau_Bx;   
  TBranch        *b_L1Tau_TowerIPhi;   
  TBranch        *b_L1Tau_TowerIEta;   
  TBranch        *b_L1Tau_RawEt;   
  TBranch        *b_L1Tau_IsoEt;   
  TBranch        *b_L1Tau_NTT;   
  TBranch        *b_L1Tau_HasEM;   
  TBranch        *b_L1Tau_IsMerged;   
  TBranch        *b_L1Tau_HwQual;   

  TBranch        *b_L1Jet_N;   
  TBranch        *b_L1Jet_Et;   
  TBranch        *b_L1Jet_Eta;   
  TBranch        *b_L1Jet_Phi;   
  TBranch        *b_L1Jet_IEt;   
  TBranch        *b_L1Jet_IEta;   
  TBranch        *b_L1Jet_IPhi;   
  TBranch        *b_L1Jet_Bx;   
  TBranch        *b_L1Jet_TowerIPhi;   
  TBranch        *b_L1Jet_TowerIEta;   
  TBranch        *b_L1Jet_RawEt;   
  TBranch        *b_L1Jet_SeedEt;   
  TBranch        *b_L1Jet_PUEt;   
  TBranch        *b_L1Jet_PUDonutEt0;   
  TBranch        *b_L1Jet_PUDonutEt1;   
  TBranch        *b_L1Jet_PUDonutEt2;   
  TBranch        *b_L1Jet_PUDonutEt3;   

  TBranch        *b_L1Muon_N;   
  TBranch        *b_L1Muon_Et;   
  TBranch        *b_L1Muon_Eta;   
  TBranch        *b_L1Muon_Phi;   
  TBranch        *b_L1Muon_EtaAtVtx;   
  TBranch        *b_L1Muon_PhiAtVtx;   
  TBranch        *b_L1Muon_IEt;   
  TBranch        *b_L1Muon_IEta;   
  TBranch        *b_L1Muon_IPhi;   
  TBranch        *b_L1Muon_IEtaAtVtx;   
  TBranch        *b_L1Muon_IPhiAtVtx;   
  TBranch        *b_L1Muon_IDEta;   
  TBranch        *b_L1Muon_IDPhi;   
  TBranch        *b_L1Muon_Chg;   
  TBranch        *b_L1Muon_Iso;   
  TBranch        *b_L1Muon_Qual;   
  TBranch        *b_L1Muon_TfMuon_Idx;   
  TBranch        *b_L1Muon_Bx;   

  TBranch        *b_L1Sum_N;   
  TBranch        *b_L1Sum_Type;   
  TBranch        *b_L1Sum_Et;   
  TBranch        *b_L1Sum_Phi;   
  TBranch        *b_L1Sum_IEt;   
  TBranch        *b_L1Sum_IPhi;   
  TBranch        *b_L1Sum_Bx;   

  // uGT
  TBranch        *b_m_orbitNr;   
  TBranch        *b_m_bxNr;   
  TBranch        *b_m_bxInEvent;   
  TBranch        *b_m_finalOR;   
  TBranch        *b_m_finalORPreVeto;   
  TBranch        *b_m_finalORVeto;   
  TBranch        *b_m_preScColumn;   
  TBranch        *b_m_algoDecisionInitial;   
  TBranch        *b_m_algoDecisionPreScaled;   
  TBranch        *b_m_algoDecisionFinal;   
  
  // HO
  TBranch        *b_hcalDetId_N;   
  TBranch        *b_hcalQIESample_N;   
  TBranch        *b_hcalDetId_IEta;   
  TBranch        *b_hcalDetId_IPhi;   
  TBranch        *b_hcalQIESample;   
  TBranch        *b_hcalQIESample_Adc;   
  TBranch        *b_hcalQIESample_Dv;   
  TBranch        *b_hcalQIESample_Er;   

  // Tracks 
  TBranch        *b_L1Tks_Pt;   
  TBranch        *b_L1Tks_Eta;   
  TBranch        *b_L1Tks_Phi;   
  TBranch        *b_L1Tks_d0;   
  TBranch        *b_L1Tks_z0;   
  TBranch        *b_L1Tks_ChiSquared;   
  TBranch        *b_L1Tks_NStubs;   
  TBranch        *b_L1Tks_StubPtConsistency;
  TBranch        *b_L1Tks_RInv;
  TBranch        *b_L1Tks_IsGenuine;   
  TBranch        *b_L1Tks_IsLoose;   
  TBranch        *b_L1Tks_IsUnknown;   
  TBranch        *b_L1Tks_IsCombinatoric;   
  TBranch        *b_L1Tks_IsFake;   

  TBranch        *b_L1Tks_TP_PdgId;   
  TBranch        *b_L1Tks_TP_Pt;   
  TBranch        *b_L1Tks_TP_Eta;   
  TBranch        *b_L1Tks_TP_Phi;   
  TBranch        *b_L1Tks_TP_z0;   
  TBranch        *b_L1Tks_TP_dxy;   

  // Tracking Particles
  TBranch        *b_TP_Pt;   
  TBranch        *b_TP_Eta;   
  TBranch        *b_TP_Phi;   
  TBranch        *b_TP_dxy;   
  TBranch        *b_TP_d0;   
  TBranch        *b_TP_z0;   
  TBranch        *b_TP_d0_produced;   
  TBranch        *b_TP_z0_produced;   
  TBranch        *b_TP_PdgId;   
  TBranch        *b_TP_NMatch;   
  TBranch        *b_TP_NStubs;   
  TBranch        *b_TP_EventId;   
  TBranch        *b_TP_Charge;   

  TBranch        *b_TP_Trk_Pt;   
  TBranch        *b_TP_Trk_Eta;   
  TBranch        *b_TP_Trk_Phi;   
  TBranch        *b_TP_Trk_z0;   
  TBranch        *b_TP_Trk_d0;   
  TBranch        *b_TP_Trk_ChiSquared;   
  TBranch        *b_TP_Trk_NStubs;   

  // =========================================== EMULATOR TREES =============================================

  // Calo Tower Emulator
  //==CaloTP (Emu)
  TBranch        *b_hcalTPEmu_N;  
  TBranch        *b_hcalTPEmu_ieta;
  TBranch        *b_hcalTPEmu_iphi;
  TBranch        *b_hcalTPEmu_Caliphi;
  TBranch        *b_hcalTPEmu_et;   
  TBranch        *b_hcalTPEmu_compEt;
  TBranch        *b_hcalTPEmu_fineGrain; 

  TBranch        *b_ecalTPEmu_N;  
  TBranch        *b_ecalTPEmu_ieta; 
  TBranch        *b_ecalTPEmu_iphi; 
  TBranch        *b_ecalTPEmu_Caliphi; 
  TBranch        *b_ecalTPEmu_et;   
  TBranch        *b_ecalTPEmu_compEt;
  TBranch        *b_ecalTPEmu_fineGrain; 

  TBranch        *b_ecalEBTPEmu_N;  
  TBranch        *b_ecalEBTPEmu_ieta;
  TBranch        *b_ecalEBTPEmu_iphi;
  TBranch        *b_ecalEBTPEmu_Caliphi; 
  TBranch        *b_ecalEBTPEmu_et;   
  TBranch        *b_ecalEBTPEmu_encodEt;
  TBranch        *b_ecalEBTPEmu_spike; 
  TBranch        *b_ecalEBTPEmu_time;  
  
  //
  TBranch        *b_CaloTowerEmu;
  TBranch        *b_CaloClusterEmu;
  //

  //==Calo Towers (Emu)
  TBranch        *b_caloTowerEmu_N; 
  TBranch        *b_caloTowerEmu_ieta;  
  TBranch        *b_caloTowerEmu_iphi;  
  TBranch        *b_caloTowerEmu_iet;  
  TBranch        *b_caloTowerEmu_iem;  
  TBranch        *b_caloTowerEmu_ihad; 
  TBranch        *b_caloTowerEmu_iratio;
  TBranch        *b_caloTowerEmu_iqual; 
  TBranch        *b_caloTowerEmu_et;   
  TBranch        *b_caloTowerEmu_eta;  
  TBranch        *b_caloTowerEmu_phi;  

  //==Calo Clusters (Emu)
  TBranch        *b_caloClusterEmu_N;   
  TBranch        *b_caloClusterEmu_ieta;   
  TBranch        *b_caloClusterEmu_iphi;   
  TBranch        *b_caloClusterEmu_iet;   
  TBranch        *b_caloClusterEmu_iqual;   
  TBranch        *b_caloClusterEmu_et;   
  TBranch        *b_caloClusterEmu_eta;   
  TBranch        *b_caloClusterEmu_phi;   

  // Upgrade Emulator
  TBranch        *b_L1EGEmu_N;  
  TBranch        *b_L1EGEmu_Et;   
  TBranch        *b_L1EGEmu_Eta;   
  TBranch        *b_L1EGEmu_Phi;   
  TBranch        *b_L1EGEmu_IEt;   
  TBranch        *b_L1EGEmu_IEta;   
  TBranch        *b_L1EGEmu_IPhi;   
  TBranch        *b_L1EGEmu_Iso;   
  TBranch        *b_L1EGEmu_Bx;   
  TBranch        *b_L1EGEmu_TowerIPhi;   
  TBranch        *b_L1EGEmu_TowerIEta;   
  TBranch        *b_L1EGEmu_RawEt;   
  TBranch        *b_L1EGEmu_IsoEt;   
  TBranch        *b_L1EGEmu_FootprintEt;   
  TBranch        *b_L1EGEmu_NTT;   
  TBranch        *b_L1EGEmu_Shape;   
  TBranch        *b_L1EGEmu_TowerHoE;   

  TBranch        *b_L1TauEmu_N;   
  TBranch        *b_L1TauEmu_Et;   
  TBranch        *b_L1TauEmu_Eta;   
  TBranch        *b_L1TauEmu_Phi;   
  TBranch        *b_L1TauEmu_IEt;   
  TBranch        *b_L1TauEmu_IEta;   
  TBranch        *b_L1TauEmu_IPhi;   
  TBranch        *b_L1TauEmu_Iso;   
  TBranch        *b_L1TauEmu_Bx;   
  TBranch        *b_L1TauEmu_TowerIPhi;   
  TBranch        *b_L1TauEmu_TowerIEta;   
  TBranch        *b_L1TauEmu_RawEt;   
  TBranch        *b_L1TauEmu_IsoEt;   
  TBranch        *b_L1TauEmu_NTT;   
  TBranch        *b_L1TauEmu_HasEM;   
  TBranch        *b_L1TauEmu_IsMerged;   
  TBranch        *b_L1TauEmu_HwQual;   

  TBranch        *b_L1JetEmu_N;   
  TBranch        *b_L1JetEmu_Et;   
  TBranch        *b_L1JetEmu_Eta;   
  TBranch        *b_L1JetEmu_Phi;   
  TBranch        *b_L1JetEmu_IEt;   
  TBranch        *b_L1JetEmu_IEta;   
  TBranch        *b_L1JetEmu_IPhi;   
  TBranch        *b_L1JetEmu_Bx;   
  TBranch        *b_L1JetEmu_TowerIPhi;   
  TBranch        *b_L1JetEmu_TowerIEta;   
  TBranch        *b_L1JetEmu_RawEt;   
  TBranch        *b_L1JetEmu_SeedEt;   
  TBranch        *b_L1JetEmu_PUEt;   
  TBranch        *b_L1JetEmu_PUDonutEt0;   
  TBranch        *b_L1JetEmu_PUDonutEt1;   
  TBranch        *b_L1JetEmu_PUDonutEt2;   
  TBranch        *b_L1JetEmu_PUDonutEt3;   

  TBranch        *b_L1MuonEmu_N;   
  TBranch        *b_L1MuonEmu_Et;   
  TBranch        *b_L1MuonEmu_Eta;   
  TBranch        *b_L1MuonEmu_Phi;   
  TBranch        *b_L1MuonEmu_EtaAtVtx;   
  TBranch        *b_L1MuonEmu_PhiAtVtx;   
  TBranch        *b_L1MuonEmu_IEt;   
  TBranch        *b_L1MuonEmu_IEta;   
  TBranch        *b_L1MuonEmu_IPhi;   
  TBranch        *b_L1MuonEmu_IEtaAtVtx;   
  TBranch        *b_L1MuonEmu_IPhiAtVtx;   
  TBranch        *b_L1MuonEmu_IDEta;   
  TBranch        *b_L1MuonEmu_IDPhi;   
  TBranch        *b_L1MuonEmu_Chg;   
  TBranch        *b_L1MuonEmu_Iso;   
  TBranch        *b_L1MuonEmu_Qual;   
  TBranch        *b_L1MuonEmu_TfMuon_Idx;   
  TBranch        *b_L1MuonEmu_Bx;   

  TBranch        *b_L1SumEmu_N;   
  TBranch        *b_L1SumEmu_Type;   
  TBranch        *b_L1SumEmu_Et;   
  TBranch        *b_L1SumEmu_Phi;   
  TBranch        *b_L1SumEmu_IEt;   
  TBranch        *b_L1SumEmu_IPhi;   
  TBranch        *b_L1SumEmu_Bx;   

  // uGT Emulator
  TBranch        *b_mEmu_orbitNr;   
  TBranch        *b_mEmu_bxNr;   
  TBranch        *b_mEmu_bxInEvent;   
  TBranch        *b_mEmu_finalOR;   
  TBranch        *b_mEmu_finalORPreVeto;   
  TBranch        *b_mEmu_finalORVeto;   
  TBranch        *b_mEmu_preScColumn;   
  TBranch        *b_mEmu_algoDecisionInitial;   
  TBranch        *b_mEmu_algoDecisionPreScaled;   
  TBranch        *b_mEmu_algoDecisionFinal;   

  virtual void InitReco(TChain *chain);
  
};

void TreeDefinitionReco::InitReco(TChain *chain)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  if (0) cout << "=== TreeDefinitionReco::InitReco()" << endl;
  
  //L1PhaseII
  EG_N     = 0;
  tkEM_N   = 0;

  // Event
  run      = 0;
  event    = 0;
  lumi     = 0;
  bx       = 0;
  orbit    = 0;
  time     = 0;
  nPV      = 0;
  nPV_True = 0;
  puWeight = 0;

  // Calo Towers
  hcalTP_N    = 0;
  ecalTP_N    = 0;
  ecalEBTP_N  = 0;
  caloTower_N = 0;

  // Upgrade
  L1EG_N    = 0;
  L1Tau_N   = 0;
  L1Jet_N   = 0;
  L1Muon_N  = 0;
  L1Sum_N   = 0;

  // uGT
  m_orbitNr        = 0;
  m_bxNr           = 0;
  m_bxInEvent      = 0;
  m_finalOR        = false; 
  m_finalORPreVeto = false; 
  m_finalORVeto    = false; 
  m_preScColumn    = 0;
  
  // HO
  hcalDetId_N     = 0;
  hcalQIESample_N = 0;

  // L1 tracks
  L1Tks_Pt = 0;
  L1Tks_Eta = 0;
  L1Tks_Phi = 0;
  L1Tks_d0 = 0;
  L1Tks_z0 = 0;
  L1Tks_ChiSquared = 0;
  L1Tks_NStubs = 0;
  L1Tks_StubPtConsistency = 0;
  L1Tks_RInv = 0;
  L1Tks_IsGenuine = 0;
  L1Tks_IsLoose = 0;
  L1Tks_IsUnknown = 0;
  L1Tks_IsCombinatoric = 0;
  L1Tks_IsFake = 0;

  L1Tks_TP_PdgId = 0;
  L1Tks_TP_Pt = 0;
  L1Tks_TP_Eta = 0;
  L1Tks_TP_Phi = 0;
  L1Tks_TP_z0 = 0;
  L1Tks_TP_dxy = 0;

  // Tracking Particles
  TP_Pt = 0;
  TP_Eta = 0;
  TP_Phi = 0;
  TP_dxy = 0;
  TP_d0 = 0;
  TP_z0 = 0;
  TP_d0_produced = 0;
  TP_z0_produced = 0;
  TP_PdgId = 0;
  TP_NMatch = 0;
  TP_NStubs = 0;
  TP_EventId = 0;
  TP_Charge = 0;
  
  TP_Trk_Pt = 0;
  TP_Trk_Eta = 0;
  TP_Trk_Phi = 0;
  TP_Trk_z0 = 0;
  TP_Trk_d0 = 0;
  TP_Trk_ChiSquared = 0;
  TP_Trk_NStubs = 0;
 
  // =========================================== EMULATOR TREES =============================================

  // Calo Towers Emulator
  hcalTPEmu_N      = 0;
  ecalTPEmu_N      = 0;
  ecalEBTPEmu_N    = 0;
  caloTowerEmu_N   = 0;
  caloClusterEmu_N = 0;

  // Upgrade Emulator
  L1EGEmu_N    = 0;
  L1TauEmu_N   = 0;
  L1JetEmu_N   = 0;
  L1MuonEmu_N  = 0;
  L1SumEmu_N   = 0;

  // uGT Emulator
  mEmu_orbitNr        = 0;
  mEmu_bxNr           = 0;
  mEmu_bxInEvent      = 0;
  mEmu_finalOR        = false; 
  mEmu_finalORPreVeto = false; 
  mEmu_finalORVeto    = false; 
  mEmu_preScColumn    = 0;

  if (0) cout << "\tSetting branch addresses and branch pointers." << endl;
  if (!chain) return;
  fChain = chain ;
  fCurrent = -1;

  // Event (MainChain)
  if(1)                                                                                                                                                                 
    { 
      fChain->SetBranchAddress("run", &run, &b_Event_run);
      fChain->SetBranchAddress("event", &event, &b_Event_event);
      fChain->SetBranchAddress("lumi", &lumi, &b_Event_lumi);
      fChain->SetBranchAddress("bx", &bx, &b_Event_bx);
      fChain->SetBranchAddress("orbit", &orbit, &b_Event_orbit);
      fChain->SetBranchAddress("time", &time, &b_Event_time);
      fChain->SetBranchAddress("nPV", &nPV, &b_Event_nPV);
      fChain->SetBranchAddress("nPV_True", &nPV_True, &b_Event_nPV_True);
      fChain->SetBranchAddress("hlt", &hlt, &b_Event_hlt);
      fChain->SetBranchAddress("puWeight", &puWeight, &b_Event_puWeight);
    }
  

  //L1PhaseII
  if(doL1PhaseII)
    {
      fL1PhaseII->SetBranchAddress("nEG", &EG_N, &b_EG_N);
      fL1PhaseII->SetBranchAddress("EGEt", &EG_Et, &b_EG_Et);
      fL1PhaseII->SetBranchAddress("EGEta", &EG_Eta, &b_EG_Eta);
      fL1PhaseII->SetBranchAddress("EGPhi", &EG_Phi, &b_EG_Phi);
      fL1PhaseII->SetBranchAddress("EGBx", &EG_Bx, &b_EG_Bx);
      fL1PhaseII->SetBranchAddress("EGIso", &EG_Iso, &b_EG_Iso);
      fL1PhaseII->SetBranchAddress("EGzVtx", &EG_zVtx, &b_EG_zVtx);
      fL1PhaseII->SetBranchAddress("EGHwQual", &EG_HwQual, &b_EG_HwQual);

      fL1PhaseII->SetBranchAddress("nTkEM", &tkEM_N, &b_tkEM_N);
      fL1PhaseII->SetBranchAddress("tkEMEt", &tkEM_Et, &b_tkEM_Et);
      fL1PhaseII->SetBranchAddress("tkEMEta", &tkEM_Eta, &b_tkEM_Eta);
      fL1PhaseII->SetBranchAddress("tkEMPhi", &tkEM_Phi, &b_tkEM_Phi);
      fL1PhaseII->SetBranchAddress("tkEMBx", &tkEM_Bx, &b_tkEM_Bx);
      fL1PhaseII->SetBranchAddress("tkEMTrkIso", &tkEM_TrkIso, &b_tkEM_TrkIso);
      fL1PhaseII->SetBranchAddress("tkEMzVtx", &tkEM_zVtx, &b_tkEM_zVtx);
      fL1PhaseII->SetBranchAddress("tkEMHwQual", &tkEM_HwQual, &b_tkEM_HwQual);
      fL1PhaseII->SetBranchAddress("tkEMEGRefPt", &tkEM_EGRefPt, &b_tkEM_EGRefPt);
      fL1PhaseII->SetBranchAddress("tkEMEGRefEta", &tkEM_EGRefEta, &b_tkEM_EGRefEta);
      fL1PhaseII->SetBranchAddress("tkEMEGRefPhi", &tkEM_EGRefPhi, &b_tkEM_EGRefPhi);

      fChain -> AddFriend(fL1PhaseII);
    }

  /*
  // Calo Towers
  if(doCaloTower)
    {
      // Calo TPs
      fCaloTower->SetBranchAddress("nHCALTP", &hcalTP_N, &b_hcalTP_N);
      fCaloTower->SetBranchAddress("hcalTPieta", &hcalTP_ieta, &b_hcalTP_ieta);
      fCaloTower->SetBranchAddress("hcalTPiphi", &hcalTP_iphi, &b_hcalTP_iphi);
      fCaloTower->SetBranchAddress("hcalTPCaliphi", &hcalTP_Caliphi, &b_hcalTP_Caliphi);
      fCaloTower->SetBranchAddress("hcalTPet", &hcalTP_et, &b_hcalTP_et);
      fCaloTower->SetBranchAddress("hcalTPcompEt", &hcalTP_compEt, &b_hcalTP_compEt);
      fCaloTower->SetBranchAddress("hcalTPfineGrain", &hcalTP_fineGrain, &b_hcalTP_fineGrain);

      fCaloTower->SetBranchAddress("nECALTP", &ecalTP_N, &b_ecalTP_N);
      fCaloTower->SetBranchAddress("ecalTPieta", &ecalTP_ieta, &b_ecalTP_ieta);
      fCaloTower->SetBranchAddress("ecalTPiphi", &ecalTP_iphi, &b_ecalTP_iphi);
      fCaloTower->SetBranchAddress("ecalTPCaliphi", &ecalTP_Caliphi, &b_ecalTP_Caliphi);
      fCaloTower->SetBranchAddress("ecalTPet", &ecalTP_et, &b_ecalTP_et);
      fCaloTower->SetBranchAddress("ecalTPcompEt", &ecalTP_compEt, &b_ecalTP_compEt);
      fCaloTower->SetBranchAddress("ecalTPfineGrain", &ecalTP_fineGrain, &b_ecalTP_fineGrain);

      fCaloTower->SetBranchAddress("nECALEBTP", &ecalEBTP_N, &b_ecalEBTP_N);
      fCaloTower->SetBranchAddress("ecalEBTPieta", &ecalEBTP_ieta, &b_ecalEBTP_ieta);
      fCaloTower->SetBranchAddress("ecalEBTPiphi", &ecalEBTP_iphi, &b_ecalEBTP_iphi);
      fCaloTower->SetBranchAddress("ecalEBTPCaliphi", &ecalEBTP_Caliphi, &b_ecalEBTP_Caliphi);
      fCaloTower->SetBranchAddress("ecalEBTPet", &ecalEBTP_et, &b_ecalEBTP_et);
      fCaloTower->SetBranchAddress("ecalEBTPencodEt", &ecalEBTP_encodEt, &b_ecalEBTP_encodEt);
      fCaloTower->SetBranchAddress("ecalEBTPspike", &ecalEBTP_spike, &b_ecalEBTP_spike);
      fCaloTower->SetBranchAddress("ecalEBTPtime", &ecalEBTP_time, &b_ecalEBTP_time);

      // Calo Towers
      fCaloTower->SetBranchAddress("nTower", &caloTower_N, &b_caloTower_N);
      fCaloTower->SetBranchAddress("ieta", &caloTower_ieta, &b_caloTower_ieta);
      fCaloTower->SetBranchAddress("iphi", &caloTower_iphi, &b_caloTower_iphi);
      fCaloTower->SetBranchAddress("iet", &caloTower_iet, &b_caloTower_iet);
      fCaloTower->SetBranchAddress("iem", &caloTower_iem, &b_caloTower_iem);
      fCaloTower->SetBranchAddress("ihad", &caloTower_ihad, &b_caloTower_ihad);
      fCaloTower->SetBranchAddress("iratio", &caloTower_iratio, &b_caloTower_iratio);
      fCaloTower->SetBranchAddress("iqual", &caloTower_iqual, &b_caloTower_iqual);
      fCaloTower->SetBranchAddress("et", &caloTower_et, &b_caloTower_et);
      fCaloTower->SetBranchAddress("eta", &caloTower_eta, &b_caloTower_eta);
      fCaloTower->SetBranchAddress("phi", &caloTower_phi, &b_caloTower_phi);

      // Add friend (CaloTowerTree)
      fChain -> AddFriend(fCaloTower);
    }

  // Upgrade 
  if(doUpgrade)
    {
      fUpgrade->SetBranchAddress("nEGs", &L1EG_N, &b_L1EG_N);
      fUpgrade->SetBranchAddress("egEt", &L1EG_Et, &b_L1EG_Et);
      fUpgrade->SetBranchAddress("egEta", &L1EG_Eta, &b_L1EG_Eta);
      fUpgrade->SetBranchAddress("egPhi", &L1EG_Phi, &b_L1EG_Phi);
      fUpgrade->SetBranchAddress("egIEt", &L1EG_IEt, &b_L1EG_IEt);
      fUpgrade->SetBranchAddress("egIEta", &L1EG_IEta, &b_L1EG_IEta);
      fUpgrade->SetBranchAddress("egIPhi", &L1EG_IPhi, &b_L1EG_IPhi);
      fUpgrade->SetBranchAddress("egIso", &L1EG_Iso, &b_L1EG_Iso);
      fUpgrade->SetBranchAddress("egBx", &L1EG_Bx, &b_L1EG_Bx);
      fUpgrade->SetBranchAddress("egTowerIPhi", &L1EG_TowerIPhi, &b_L1EG_TowerIPhi);
      fUpgrade->SetBranchAddress("egTowerIEta", &L1EG_TowerIEta, &b_L1EG_TowerIEta);
      fUpgrade->SetBranchAddress("egRawEt", &L1EG_RawEt, &b_L1EG_RawEt);
      fUpgrade->SetBranchAddress("egIsoEt", &L1EG_IsoEt, &b_L1EG_IsoEt);
      fUpgrade->SetBranchAddress("egFootprintEt", &L1EG_FootprintEt, &b_L1EG_FootprintEt);
      fUpgrade->SetBranchAddress("egNTT", &L1EG_NTT, &b_L1EG_NTT);
      fUpgrade->SetBranchAddress("egShape", &L1EG_Shape, &b_L1EG_Shape);
      fUpgrade->SetBranchAddress("egTowerHoE", &L1EG_TowerHoE, &b_L1EG_TowerHoE);

      fUpgrade->SetBranchAddress("nTaus", &L1Tau_N, &b_L1Tau_N);
      fUpgrade->SetBranchAddress("tauEt", &L1Tau_Et, &b_L1Tau_Et);
      fUpgrade->SetBranchAddress("tauEta", &L1Tau_Eta, &b_L1Tau_Eta);
      fUpgrade->SetBranchAddress("tauPhi", &L1Tau_Phi, &b_L1Tau_Phi);
      fUpgrade->SetBranchAddress("tauIEt", &L1Tau_IEt, &b_L1Tau_IEt);
      fUpgrade->SetBranchAddress("tauIEta", &L1Tau_IEta, &b_L1Tau_IEta);
      fUpgrade->SetBranchAddress("tauIPhi", &L1Tau_IPhi, &b_L1Tau_IPhi);
      fUpgrade->SetBranchAddress("tauIso", &L1Tau_Iso, &b_L1Tau_Iso);
      fUpgrade->SetBranchAddress("tauBx", &L1Tau_Bx, &b_L1Tau_Bx);
      fUpgrade->SetBranchAddress("tauTowerIPhi", &L1Tau_TowerIPhi, &b_L1Tau_TowerIPhi);
      fUpgrade->SetBranchAddress("tauTowerIEta", &L1Tau_TowerIEta, &b_L1Tau_TowerIEta);
      fUpgrade->SetBranchAddress("tauRawEt", &L1Tau_RawEt, &b_L1Tau_RawEt);
      fUpgrade->SetBranchAddress("tauIsoEt", &L1Tau_IsoEt, &b_L1Tau_IsoEt);
      fUpgrade->SetBranchAddress("tauNTT", &L1Tau_NTT, &b_L1Tau_NTT);
      fUpgrade->SetBranchAddress("tauHasEM", &L1Tau_HasEM, &b_L1Tau_HasEM);
      fUpgrade->SetBranchAddress("tauIsMerged", &L1Tau_IsMerged, &b_L1Tau_IsMerged);
      fUpgrade->SetBranchAddress("tauHwQual", &L1Tau_HwQual, &b_L1Tau_HwQual);

      fUpgrade->SetBranchAddress("nJets", &L1Jet_N, &b_L1Jet_N);
      fUpgrade->SetBranchAddress("jetEt", &L1Jet_Et, &b_L1Jet_Et);
      fUpgrade->SetBranchAddress("jetEta", &L1Jet_Eta, &b_L1Jet_Eta);
      fUpgrade->SetBranchAddress("jetPhi", &L1Jet_Phi, &b_L1Jet_Phi);
      fUpgrade->SetBranchAddress("jetIEt", &L1Jet_IEt, &b_L1Jet_IEt);
      fUpgrade->SetBranchAddress("jetIEta", &L1Jet_IEta, &b_L1Jet_IEta);
      fUpgrade->SetBranchAddress("jetIPhi", &L1Jet_IPhi, &b_L1Jet_IPhi);
      fUpgrade->SetBranchAddress("jetBx", &L1Jet_Bx, &b_L1Jet_Bx);
      fUpgrade->SetBranchAddress("jetTowerIPhi", &L1Jet_TowerIPhi, &b_L1Jet_TowerIPhi);
      fUpgrade->SetBranchAddress("jetTowerIEta", &L1Jet_TowerIEta, &b_L1Jet_TowerIEta);
      fUpgrade->SetBranchAddress("jetRawEt", &L1Jet_RawEt, &b_L1Jet_RawEt);
      fUpgrade->SetBranchAddress("jetSeedEt", &L1Jet_SeedEt, &b_L1Jet_SeedEt);
      fUpgrade->SetBranchAddress("jetPUEt", &L1Jet_PUEt, &b_L1Jet_PUEt);
      fUpgrade->SetBranchAddress("jetPUDonutEt0", &L1Jet_PUDonutEt0, &b_L1Jet_PUDonutEt0);
      fUpgrade->SetBranchAddress("jetPUDonutEt1", &L1Jet_PUDonutEt1, &b_L1Jet_PUDonutEt1);
      fUpgrade->SetBranchAddress("jetPUDonutEt2", &L1Jet_PUDonutEt2, &b_L1Jet_PUDonutEt2);
      fUpgrade->SetBranchAddress("jetPUDonutEt3", &L1Jet_PUDonutEt3, &b_L1Jet_PUDonutEt3);

      fUpgrade->SetBranchAddress("nMuons", &L1Muon_N, &b_L1Muon_N);
      fUpgrade->SetBranchAddress("muonEt", &L1Muon_Et, &b_L1Muon_Et);
      fUpgrade->SetBranchAddress("muonEta", &L1Muon_Eta, &b_L1Muon_Eta);
      fUpgrade->SetBranchAddress("muonPhi", &L1Muon_Phi, &b_L1Muon_Phi);
      fUpgrade->SetBranchAddress("muonEtaAtVtx", &L1Muon_EtaAtVtx, &b_L1Muon_EtaAtVtx);
      fUpgrade->SetBranchAddress("muonPhiAtVtx", &L1Muon_PhiAtVtx, &b_L1Muon_PhiAtVtx);
      fUpgrade->SetBranchAddress("muonIEt", &L1Muon_IEt, &b_L1Muon_IEt);
      fUpgrade->SetBranchAddress("muonIEta", &L1Muon_IEta, &b_L1Muon_IEta);
      fUpgrade->SetBranchAddress("muonIPhi", &L1Muon_IPhi, &b_L1Muon_IPhi);
      fUpgrade->SetBranchAddress("muonIEtaAtVtx", &L1Muon_IEtaAtVtx, &b_L1Muon_IEtaAtVtx);
      fUpgrade->SetBranchAddress("muonIPhiAtVtx", &L1Muon_IPhiAtVtx, &b_L1Muon_IPhiAtVtx);
      fUpgrade->SetBranchAddress("muonIDEta", &L1Muon_IDEta, &b_L1Muon_IDEta);
      fUpgrade->SetBranchAddress("muonIDPhi", &L1Muon_IDPhi, &b_L1Muon_IDPhi);
      fUpgrade->SetBranchAddress("muonChg", &L1Muon_Chg, &b_L1Muon_Chg);
      fUpgrade->SetBranchAddress("muonIso", &L1Muon_Iso, &b_L1Muon_Iso);
      fUpgrade->SetBranchAddress("muonQual", &L1Muon_Qual, &b_L1Muon_Qual);
      fUpgrade->SetBranchAddress("muonTfMuonIdx", &L1Muon_TfMuon_Idx, &b_L1Muon_TfMuon_Idx);
      fUpgrade->SetBranchAddress("muonBx", &L1Muon_Bx, &b_L1Muon_Bx);

      fUpgrade->SetBranchAddress("nSums", &L1Sum_N, &b_L1Sum_N);
      fUpgrade->SetBranchAddress("sumType", &L1Sum_Type, &b_L1Sum_Type);
      fUpgrade->SetBranchAddress("sumEt", &L1Sum_Et, &b_L1Sum_Et);
      fUpgrade->SetBranchAddress("sumPhi", &L1Sum_Phi, &b_L1Sum_Phi);
      fUpgrade->SetBranchAddress("sumIEt", &L1Sum_IEt, &b_L1Sum_IEt);
      fUpgrade->SetBranchAddress("sumIPhi", &L1Sum_IPhi, &b_L1Sum_IPhi);
      fUpgrade->SetBranchAddress("sumBx", &L1Sum_Bx, &b_L1Sum_Bx);

      // Add friend (UpgradeTree)
      fChain -> AddFriend(fUpgrade);
    }

  // uGT
  if(douGT)
    {
      fuGT->SetBranchAddress("m_orbitNr", &m_orbitNr, &b_m_orbitNr);
      fuGT->SetBranchAddress("m_bxNr", &m_bxNr, &b_m_bxNr);
      fuGT->SetBranchAddress("m_bxInEvent", &m_bxInEvent, &b_m_bxInEvent);
      fuGT->SetBranchAddress("m_finalOR", &m_finalOR, &b_m_finalOR);
      fuGT->SetBranchAddress("m_finalORPreVeto", &m_finalORPreVeto, &b_m_finalORPreVeto);
      fuGT->SetBranchAddress("m_finalORVeto", &m_finalORVeto, &b_m_finalORVeto);
      fuGT->SetBranchAddress("m_preScColumn", &m_preScColumn, &b_m_preScColumn);
      fuGT->SetBranchAddress("m_algoDecisionInitial", &m_algoDecisionInitial, &b_m_algoDecisionInitial);
      fuGT->SetBranchAddress("m_algoDecisionPreScaled", &m_algoDecisionPreScaled, &b_m_algoDecisionPreScaled);
      fuGT->SetBranchAddress("m_algoDecisionFinal", &m_algoDecisionFinal, &b_m_algoDecisionFinal);

      // Add friend (uGTTree) 
      fChain -> AddFriend(fuGT);
    }
  
  // HO
  if(doHO)
    {
      fHO->SetBranchAddress("nHcalDetIds", &hcalDetId_N, &b_hcalDetId_N);
      fHO->SetBranchAddress("nHcalQIESamples", &hcalQIESample_N, &b_hcalQIESample_N);
      fHO->SetBranchAddress("hcalDetIdIEta", &hcalDetId_IEta, &b_hcalDetId_IEta);
      fHO->SetBranchAddress("hcalDetIdIPhi", &hcalDetId_IPhi, &b_hcalDetId_IPhi);
      fHO->SetBranchAddress("hcalQIESample", &hcalQIESample, &b_hcalQIESample);
      fHO->SetBranchAddress("hcalQIESampleAdc", &hcalQIESample_Adc, &b_hcalQIESample_Adc);
      fHO->SetBranchAddress("hcalQIESampleDv", &hcalQIESample_Dv, &b_hcalQIESample_Dv);
      fHO->SetBranchAddress("hcalQIESampleEr", &hcalQIESample_Er, &b_hcalQIESample_Er);

      // Add friend (HOTree)
      fChain -> AddFriend(fHO);
    }
  */
  // Tracks and Tracking Particles
  if(doTracks)
    {
      //Tracks
      fTracks->SetBranchAddress("trk_pt", &L1Tks_Pt, &b_L1Tks_Pt);
      
      fTracks->SetBranchAddress("trk_eta", &L1Tks_Eta, &b_L1Tks_Eta);
      fTracks->SetBranchAddress("trk_phi", &L1Tks_Phi, &b_L1Tks_Phi);
      fTracks->SetBranchAddress("trk_d0", &L1Tks_d0, &b_L1Tks_d0);
      fTracks->SetBranchAddress("trk_z0", &L1Tks_z0, &b_L1Tks_z0);
      fTracks->SetBranchAddress("trk_chi2", &L1Tks_ChiSquared, &b_L1Tks_ChiSquared);
      fTracks->SetBranchAddress("trk_nstub", &L1Tks_NStubs, &b_L1Tks_NStubs);
      fTracks->SetBranchAddress("trk_stubPtConsistency", &L1Tks_StubPtConsistency, &b_L1Tks_StubPtConsistency);
      fTracks->SetBranchAddress("trk_RInv", &L1Tks_RInv, &b_L1Tks_RInv);
      fTracks->SetBranchAddress("trk_genuine", &L1Tks_IsGenuine, &b_L1Tks_IsGenuine);
      fTracks->SetBranchAddress("trk_loose", &L1Tks_IsLoose, &b_L1Tks_IsLoose);
      fTracks->SetBranchAddress("trk_unknown", &L1Tks_IsUnknown, &b_L1Tks_IsUnknown);
      fTracks->SetBranchAddress("trk_combinatoric", &L1Tks_IsCombinatoric, &b_L1Tks_IsCombinatoric);
      fTracks->SetBranchAddress("trk_fake", &L1Tks_IsFake, &b_L1Tks_IsFake);
      fTracks->SetBranchAddress("trk_matchtp_pdgid", &L1Tks_TP_PdgId, &b_L1Tks_TP_PdgId);
      fTracks->SetBranchAddress("trk_matchtp_pt", &L1Tks_TP_Pt, &b_L1Tks_TP_Pt);
      fTracks->SetBranchAddress("trk_matchtp_eta", &L1Tks_TP_Eta, &b_L1Tks_TP_Eta);
      fTracks->SetBranchAddress("trk_matchtp_phi", &L1Tks_TP_Phi, &b_L1Tks_TP_Phi);
      fTracks->SetBranchAddress("trk_matchtp_z0", &L1Tks_TP_z0, &b_L1Tks_TP_z0);
      fTracks->SetBranchAddress("trk_matchtp_dxy", &L1Tks_TP_dxy, &b_L1Tks_TP_dxy);

      fTracks->SetBranchAddress("tp_pt", &TP_Pt, &b_TP_Pt);
      fTracks->SetBranchAddress("tp_eta", &TP_Eta, &b_TP_Eta);
      fTracks->SetBranchAddress("tp_phi", &TP_Phi, &b_TP_Phi);
      fTracks->SetBranchAddress("tp_dxy", &TP_dxy, &b_TP_dxy);
      fTracks->SetBranchAddress("tp_d0", &TP_d0, &b_TP_d0);
      fTracks->SetBranchAddress("tp_z0", &TP_z0, &b_TP_z0);
      fTracks->SetBranchAddress("tp_d0_prod", &TP_d0_produced, &b_TP_d0_produced);
      fTracks->SetBranchAddress("tp_z0_prod", &TP_z0_produced, &b_TP_z0_produced);
      fTracks->SetBranchAddress("tp_pdgid", &TP_PdgId, &b_TP_PdgId);
      fTracks->SetBranchAddress("tp_nmatch", &TP_NMatch, &b_TP_NMatch);
      fTracks->SetBranchAddress("tp_nstub", &TP_NStubs, &b_TP_NStubs);
      fTracks->SetBranchAddress("tp_eventid", &TP_EventId, &b_TP_EventId);
      fTracks->SetBranchAddress("tp_charge", &TP_Charge, &b_TP_Charge);

      fTracks->SetBranchAddress("matchtrk_pt", &TP_Trk_Pt, &b_TP_Trk_Pt);
      fTracks->SetBranchAddress("matchtrk_eta", &TP_Trk_Eta, &b_TP_Trk_Eta);
      fTracks->SetBranchAddress("matchtrk_phi", &TP_Trk_Phi, &b_TP_Trk_Phi);
      fTracks->SetBranchAddress("matchtrk_z0", &TP_Trk_z0, &b_TP_Trk_z0);
      fTracks->SetBranchAddress("matchtrk_d0", &TP_Trk_d0, &b_TP_Trk_d0);
      fTracks->SetBranchAddress("matchtrk_chi2", &TP_Trk_ChiSquared, &b_TP_Trk_ChiSquared);
      fTracks->SetBranchAddress("matchtrk_nstub", &TP_Trk_NStubs, &b_TP_Trk_NStubs);
      
      // Add friend ( TracksTree )
      fChain -> AddFriend(fTracks);
    }
  
  // =========================================== EMULATOR TREES =============================================
  /*
  // Calo Towers Emulator
  if(doCaloTowerEmu)
    {
      // Calo TPs (Emu)
      fCaloTowerEmu->SetBranchAddress("nHCALTP", &hcalTPEmu_N, &b_hcalTPEmu_N);
      fCaloTowerEmu->SetBranchAddress("hcalTPieta", &hcalTPEmu_ieta, &b_hcalTPEmu_ieta);
      fCaloTowerEmu->SetBranchAddress("hcalTPiphi", &hcalTPEmu_iphi, &b_hcalTPEmu_iphi);
      fCaloTowerEmu->SetBranchAddress("hcalTPCaliphi", &hcalTPEmu_Caliphi, &b_hcalTPEmu_Caliphi);
      fCaloTowerEmu->SetBranchAddress("hcalTPet", &hcalTPEmu_et, &b_hcalTPEmu_et);
      fCaloTowerEmu->SetBranchAddress("hcalTPcompEt", &hcalTPEmu_compEt, &b_hcalTPEmu_compEt);
      fCaloTowerEmu->SetBranchAddress("hcalTPfineGrain", &hcalTPEmu_fineGrain, &b_hcalTPEmu_fineGrain);

      fCaloTowerEmu->SetBranchAddress("nECALTP", &ecalTPEmu_N, &b_ecalTPEmu_N);
      fCaloTowerEmu->SetBranchAddress("ecalTPieta", &ecalTPEmu_ieta, &b_ecalTPEmu_ieta);
      fCaloTowerEmu->SetBranchAddress("ecalTPiphi", &ecalTPEmu_iphi, &b_ecalTPEmu_iphi);
      fCaloTowerEmu->SetBranchAddress("ecalTPCaliphi", &ecalTPEmu_Caliphi, &b_ecalTPEmu_Caliphi);
      fCaloTowerEmu->SetBranchAddress("ecalTPet", &ecalTPEmu_et, &b_ecalTPEmu_et);
      fCaloTowerEmu->SetBranchAddress("ecalTPcompEt", &ecalTPEmu_compEt, &b_ecalTPEmu_compEt);
      fCaloTowerEmu->SetBranchAddress("ecalTPfineGrain", &ecalTPEmu_fineGrain, &b_ecalTPEmu_fineGrain);

      fCaloTowerEmu->SetBranchAddress("nECALEBTP", &ecalEBTPEmu_N, &b_ecalEBTPEmu_N);
      fCaloTowerEmu->SetBranchAddress("ecalEBTPieta", &ecalEBTPEmu_ieta, &b_ecalEBTPEmu_ieta);
      fCaloTowerEmu->SetBranchAddress("ecalEBTPiphi", &ecalEBTPEmu_iphi, &b_ecalEBTPEmu_iphi);
      fCaloTowerEmu->SetBranchAddress("ecalEBTPCaliphi", &ecalEBTPEmu_Caliphi, &b_ecalEBTPEmu_Caliphi);
      fCaloTowerEmu->SetBranchAddress("ecalEBTPet", &ecalEBTPEmu_et, &b_ecalEBTPEmu_et);
      fCaloTowerEmu->SetBranchAddress("ecalEBTPencodEt", &ecalEBTPEmu_encodEt, &b_ecalEBTPEmu_encodEt);
      fCaloTowerEmu->SetBranchAddress("ecalEBTPspike", &ecalEBTPEmu_spike, &b_ecalEBTPEmu_spike);
      fCaloTowerEmu->SetBranchAddress("ecalEBTPtime", &ecalEBTPEmu_time, &b_ecalEBTPEmu_time);

      //
      // ---May be needed for taking both CaloTowers and CaloClusters
      // caloTower_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
      // fCaloTowerEmu->SetBranchAddress("L1CaloTower", &caloTower_ );
      //
      
      // Calo Towers (Emu)
      fCaloTowerEmu->SetBranchAddress("nTower", &caloTowerEmu_N, &b_caloTowerEmu_N);
      fCaloTowerEmu->SetBranchAddress("ieta", &caloTowerEmu_ieta, &b_caloTowerEmu_ieta);
      fCaloTowerEmu->SetBranchAddress("iphi", &caloTowerEmu_iphi, &b_caloTowerEmu_iphi);
      fCaloTowerEmu->SetBranchAddress("iet", &caloTowerEmu_iet, &b_caloTowerEmu_iet);
      fCaloTowerEmu->SetBranchAddress("iem", &caloTowerEmu_iem, &b_caloTowerEmu_iem);
      fCaloTowerEmu->SetBranchAddress("ihad", &caloTowerEmu_ihad, &b_caloTowerEmu_ihad);
      fCaloTowerEmu->SetBranchAddress("iratio", &caloTowerEmu_iratio, &b_caloTowerEmu_iratio);
      fCaloTowerEmu->SetBranchAddress("iqual", &caloTowerEmu_iqual, &b_caloTowerEmu_iqual);
      fCaloTowerEmu->SetBranchAddress("et", &caloTowerEmu_et, &b_caloTowerEmu_et);
      fCaloTowerEmu->SetBranchAddress("eta", &caloTowerEmu_eta, &b_caloTowerEmu_eta);
      fCaloTowerEmu->SetBranchAddress("phi", &caloTowerEmu_phi, &b_caloTowerEmu_phi);
      
      // Fix me - Not working
      // Calo Clusters (Emu)
      //fCaloTowerEmu->SetBranchAddress("ieta", &caloClusterEmu_ieta, &b_caloClusterEmu_ieta);
      //fCaloTowerEmu->SetBranchAddress("iphi", &caloClusterEmu_iphi, &b_caloClusterEmu_iphi);
      //fCaloTowerEmu->SetBranchAddress("iet", &caloClusterEmu_iet, &b_caloClusterEmu_iet);
      //fCaloTowerEmu->SetBranchAddress("iqual", &caloClusterEmu_iqual, &b_caloClusterEmu_iqual);
      //fCaloTowerEmu->SetBranchAddress("et", &caloClusterEmu_et, &b_caloClusterEmu_et);
      //fCaloTowerEmu->SetBranchAddress("eta", &caloClusterEmu_eta, &b_caloClusterEmu_eta);
      //fCaloTowerEmu->SetBranchAddress("phi", &caloClusterEmu_phi, &b_caloClusterEmu_phi);

      // Add Friend (CaloTowerEmuTree) 
      fChain -> AddFriend(fCaloTowerEmu);
    }
  */
  // Upgrade Emulator
  if(doUpgradeEmu)
    {
      fUpgradeEmu->SetBranchAddress("nEGs", &L1EGEmu_N, &b_L1EGEmu_N);
      fUpgradeEmu->SetBranchAddress("egEt", &L1EGEmu_Et, &b_L1EGEmu_Et);
      fUpgradeEmu->SetBranchAddress("egEta", &L1EGEmu_Eta, &b_L1EGEmu_Eta);
      fUpgradeEmu->SetBranchAddress("egPhi", &L1EGEmu_Phi, &b_L1EGEmu_Phi);
      fUpgradeEmu->SetBranchAddress("egIEt", &L1EGEmu_IEt, &b_L1EGEmu_IEt);
      fUpgradeEmu->SetBranchAddress("egIEta", &L1EGEmu_IEta, &b_L1EGEmu_IEta);
      fUpgradeEmu->SetBranchAddress("egIPhi", &L1EGEmu_IPhi, &b_L1EGEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("egIso", &L1EGEmu_Iso, &b_L1EGEmu_Iso);
      fUpgradeEmu->SetBranchAddress("egBx", &L1EGEmu_Bx, &b_L1EGEmu_Bx);
      fUpgradeEmu->SetBranchAddress("egTowerIPhi", &L1EGEmu_TowerIPhi, &b_L1EGEmu_TowerIPhi);
      fUpgradeEmu->SetBranchAddress("egTowerIEta", &L1EGEmu_TowerIEta, &b_L1EGEmu_TowerIEta);
      fUpgradeEmu->SetBranchAddress("egRawEt", &L1EGEmu_RawEt, &b_L1EGEmu_RawEt);
      fUpgradeEmu->SetBranchAddress("egIsoEt", &L1EGEmu_IsoEt, &b_L1EGEmu_IsoEt);
      fUpgradeEmu->SetBranchAddress("egFootprintEt", &L1EGEmu_FootprintEt, &b_L1EGEmu_FootprintEt);
      fUpgradeEmu->SetBranchAddress("egNTT", &L1EGEmu_NTT, &b_L1EGEmu_NTT);
      fUpgradeEmu->SetBranchAddress("egShape", &L1EGEmu_Shape, &b_L1EGEmu_Shape);
      fUpgradeEmu->SetBranchAddress("egTowerHoE", &L1EGEmu_TowerHoE, &b_L1EGEmu_TowerHoE);

      fUpgradeEmu->SetBranchAddress("nTaus", &L1TauEmu_N, &b_L1TauEmu_N);
      fUpgradeEmu->SetBranchAddress("tauEt", &L1TauEmu_Et, &b_L1TauEmu_Et);
      fUpgradeEmu->SetBranchAddress("tauEta", &L1TauEmu_Eta, &b_L1TauEmu_Eta);
      fUpgradeEmu->SetBranchAddress("tauPhi", &L1TauEmu_Phi, &b_L1TauEmu_Phi);
      fUpgradeEmu->SetBranchAddress("tauIEt", &L1TauEmu_IEt, &b_L1TauEmu_IEt);
      fUpgradeEmu->SetBranchAddress("tauIEta", &L1TauEmu_IEta, &b_L1TauEmu_IEta);
      fUpgradeEmu->SetBranchAddress("tauIPhi", &L1TauEmu_IPhi, &b_L1TauEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("tauIso", &L1TauEmu_Iso, &b_L1TauEmu_Iso);
      fUpgradeEmu->SetBranchAddress("tauBx", &L1TauEmu_Bx, &b_L1TauEmu_Bx);
      fUpgradeEmu->SetBranchAddress("tauTowerIPhi", &L1TauEmu_TowerIPhi, &b_L1TauEmu_TowerIPhi);
      fUpgradeEmu->SetBranchAddress("tauTowerIEta", &L1TauEmu_TowerIEta, &b_L1TauEmu_TowerIEta);
      fUpgradeEmu->SetBranchAddress("tauRawEt", &L1TauEmu_RawEt, &b_L1TauEmu_RawEt);
      fUpgradeEmu->SetBranchAddress("tauIsoEt", &L1TauEmu_IsoEt, &b_L1TauEmu_IsoEt);
      fUpgradeEmu->SetBranchAddress("tauNTT", &L1TauEmu_NTT, &b_L1TauEmu_NTT);
      fUpgradeEmu->SetBranchAddress("tauHasEM", &L1TauEmu_HasEM, &b_L1TauEmu_HasEM);
      fUpgradeEmu->SetBranchAddress("tauIsMerged", &L1TauEmu_IsMerged, &b_L1TauEmu_IsMerged);
      fUpgradeEmu->SetBranchAddress("tauHwQual", &L1TauEmu_HwQual, &b_L1TauEmu_HwQual);

      fUpgradeEmu->SetBranchAddress("nJets", &L1JetEmu_N, &b_L1JetEmu_N);
      fUpgradeEmu->SetBranchAddress("jetEt", &L1JetEmu_Et, &b_L1JetEmu_Et);
      fUpgradeEmu->SetBranchAddress("jetEta", &L1JetEmu_Eta, &b_L1JetEmu_Eta);
      fUpgradeEmu->SetBranchAddress("jetPhi", &L1JetEmu_Phi, &b_L1JetEmu_Phi);
      fUpgradeEmu->SetBranchAddress("jetIEt", &L1JetEmu_IEt, &b_L1JetEmu_IEt);
      fUpgradeEmu->SetBranchAddress("jetIEta", &L1JetEmu_IEta, &b_L1JetEmu_IEta);
      fUpgradeEmu->SetBranchAddress("jetIPhi", &L1JetEmu_IPhi, &b_L1JetEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("jetBx", &L1JetEmu_Bx, &b_L1JetEmu_Bx);
      fUpgradeEmu->SetBranchAddress("jetTowerIPhi", &L1JetEmu_TowerIPhi, &b_L1JetEmu_TowerIPhi);
      fUpgradeEmu->SetBranchAddress("jetTowerIEta", &L1JetEmu_TowerIEta, &b_L1JetEmu_TowerIEta);
      fUpgradeEmu->SetBranchAddress("jetRawEt", &L1JetEmu_RawEt, &b_L1JetEmu_RawEt);
      fUpgradeEmu->SetBranchAddress("jetSeedEt", &L1JetEmu_SeedEt, &b_L1JetEmu_SeedEt);
      fUpgradeEmu->SetBranchAddress("jetPUEt", &L1JetEmu_PUEt, &b_L1JetEmu_PUEt);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt0", &L1JetEmu_PUDonutEt0, &b_L1JetEmu_PUDonutEt0);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt1", &L1JetEmu_PUDonutEt1, &b_L1JetEmu_PUDonutEt1);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt2", &L1JetEmu_PUDonutEt2, &b_L1JetEmu_PUDonutEt2);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt3", &L1JetEmu_PUDonutEt3, &b_L1JetEmu_PUDonutEt3);

      fUpgradeEmu->SetBranchAddress("nMuons", &L1MuonEmu_N, &b_L1MuonEmu_N);
      fUpgradeEmu->SetBranchAddress("muonEt", &L1MuonEmu_Et, &b_L1MuonEmu_Et);
      fUpgradeEmu->SetBranchAddress("muonEta", &L1MuonEmu_Eta, &b_L1MuonEmu_Eta);
      fUpgradeEmu->SetBranchAddress("muonPhi", &L1MuonEmu_Phi, &b_L1MuonEmu_Phi);
      fUpgradeEmu->SetBranchAddress("muonEtaAtVtx", &L1MuonEmu_EtaAtVtx, &b_L1MuonEmu_EtaAtVtx);
      fUpgradeEmu->SetBranchAddress("muonPhiAtVtx", &L1MuonEmu_PhiAtVtx, &b_L1MuonEmu_PhiAtVtx);
      fUpgradeEmu->SetBranchAddress("muonIEt", &L1MuonEmu_IEt, &b_L1MuonEmu_IEt);
      fUpgradeEmu->SetBranchAddress("muonIEta", &L1MuonEmu_IEta, &b_L1MuonEmu_IEta);
      fUpgradeEmu->SetBranchAddress("muonIPhi", &L1MuonEmu_IPhi, &b_L1MuonEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("muonIEtaAtVtx", &L1MuonEmu_IEtaAtVtx, &b_L1MuonEmu_IEtaAtVtx);
      fUpgradeEmu->SetBranchAddress("muonIPhiAtVtx", &L1MuonEmu_IPhiAtVtx, &b_L1MuonEmu_IPhiAtVtx);
      fUpgradeEmu->SetBranchAddress("muonIDEta", &L1MuonEmu_IDEta, &b_L1MuonEmu_IDEta);
      fUpgradeEmu->SetBranchAddress("muonIDPhi", &L1MuonEmu_IDPhi, &b_L1MuonEmu_IDPhi);
      fUpgradeEmu->SetBranchAddress("muonChg", &L1MuonEmu_Chg, &b_L1MuonEmu_Chg);
      fUpgradeEmu->SetBranchAddress("muonIso", &L1MuonEmu_Iso, &b_L1MuonEmu_Iso);
      fUpgradeEmu->SetBranchAddress("muonQual", &L1MuonEmu_Qual, &b_L1MuonEmu_Qual);
      fUpgradeEmu->SetBranchAddress("muonTfMuonIdx", &L1MuonEmu_TfMuon_Idx, &b_L1MuonEmu_TfMuon_Idx);
      fUpgradeEmu->SetBranchAddress("muonBx", &L1MuonEmu_Bx, &b_L1MuonEmu_Bx);

      fUpgradeEmu->SetBranchAddress("nSums", &L1SumEmu_N, &b_L1SumEmu_N);
      fUpgradeEmu->SetBranchAddress("sumType", &L1SumEmu_Type, &b_L1SumEmu_Type);
      fUpgradeEmu->SetBranchAddress("sumEt", &L1SumEmu_Et, &b_L1SumEmu_Et);
      fUpgradeEmu->SetBranchAddress("sumPhi", &L1SumEmu_Phi, &b_L1SumEmu_Phi);
      fUpgradeEmu->SetBranchAddress("sumIEt", &L1SumEmu_IEt, &b_L1SumEmu_IEt);
      fUpgradeEmu->SetBranchAddress("sumIPhi", &L1SumEmu_IPhi, &b_L1SumEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("sumBx", &L1SumEmu_Bx, &b_L1SumEmu_Bx);

      // Add friend (UpgradeEmuTree)
      fChain -> AddFriend(fUpgradeEmu);
    }
  /*
  // uGT Emulator
  if(douGTEmu)
    {
      fuGTEmu->SetBranchAddress("m_orbitNr", &mEmu_orbitNr, &b_mEmu_orbitNr);
      fuGTEmu->SetBranchAddress("m_bxNr", &mEmu_bxNr, &b_mEmu_bxNr);
      fuGTEmu->SetBranchAddress("m_bxInEvent", &mEmu_bxInEvent, &b_mEmu_bxInEvent);
      fuGTEmu->SetBranchAddress("m_finalOR", &mEmu_finalOR, &b_mEmu_finalOR);
      fuGTEmu->SetBranchAddress("m_finalORPreVeto", &mEmu_finalORPreVeto, &b_mEmu_finalORPreVeto);
      fuGTEmu->SetBranchAddress("m_finalORVeto", &mEmu_finalORVeto, &b_mEmu_finalORVeto);
      fuGTEmu->SetBranchAddress("m_preScColumn", &mEmu_preScColumn, &b_mEmu_preScColumn);
      fuGTEmu->SetBranchAddress("m_algoDecisionInitial", &mEmu_algoDecisionInitial, &b_mEmu_algoDecisionInitial);
      fuGTEmu->SetBranchAddress("m_algoDecisionPreScaled", &mEmu_algoDecisionPreScaled, &b_mEmu_algoDecisionPreScaled);
      fuGTEmu->SetBranchAddress("m_algoDecisionFinal", &mEmu_algoDecisionFinal, &b_mEmu_algoDecisionFinal);

      // Add friend (uGTEmuTree) 
      fChain -> AddFriend(fuGTEmu);
    }
  */
  // Set Make Class for all trees 
  fChain->SetMakeClass(1);
  
  TFriendElement *friendTreeElement = (TFriendElement*)fChain->GetListOfFriends()->First();
  
  while (friendTreeElement){
    friendTreeElement->GetTree()->SetMakeClass(1);
    friendTreeElement = (TFriendElement*)fChain->GetListOfFriends()->After(friendTreeElement);
  }                                                                                                                                                                     

  Notify();
}

#endif // TreeDefinitionReco_h

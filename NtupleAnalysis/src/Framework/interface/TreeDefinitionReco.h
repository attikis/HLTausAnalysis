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
  UShort_t        eg_N;
  vector<float>   eg_Et;
  vector<float>   eg_Eta;
  vector<float>   eg_Phi;
  vector<short>   eg_IEt;
  vector<short>   eg_IEta;
  vector<short>   eg_IPhi;
  vector<short>   eg_Iso;
  vector<short>   eg_Bx;
  vector<short>   eg_TowerIPhi;
  vector<short>   eg_TowerIEta;
  vector<short>   eg_RawEt;
  vector<short>   eg_IsoEt;
  vector<short>   eg_FootprintEt;
  vector<short>   eg_NTT;
  vector<short>   eg_Shape;
  vector<short>   eg_TowerHoE;

  UShort_t        tau_N;
  vector<float>   tau_Et;
  vector<float>   tau_Eta;
  vector<float>   tau_Phi;
  vector<short>   tau_IEt;
  vector<short>   tau_IEta;
  vector<short>   tau_IPhi;
  vector<short>   tau_Iso;
  vector<short>   tau_Bx;
  vector<short>   tau_TowerIPhi;
  vector<short>   tau_TowerIEta;
  vector<short>   tau_RawEt;
  vector<short>   tau_IsoEt;
  vector<short>   tau_NTT;
  vector<short>   tau_HasEM;
  vector<short>   tau_IsMerged;
  vector<short>   tau_HwQual;

  UShort_t        jet_N;
  vector<float>   jet_Et;
  vector<float>   jet_Eta;
  vector<float>   jet_Phi;
  vector<short>   jet_IEt;
  vector<short>   jet_IEta;
  vector<short>   jet_IPhi;
  vector<short>   jet_Bx;
  vector<short>   jet_TowerIPhi;
  vector<short>   jet_TowerIEta;
  vector<short>   jet_RawEt;
  vector<short>   jet_SeedEt;
  vector<short>   jet_PUEt;
  vector<short>   jet_PUDonutEt0;
  vector<short>   jet_PUDonutEt1;
  vector<short>   jet_PUDonutEt2;
  vector<short>   jet_PUDonutEt3;

  UShort_t        muon_N;
  vector<float>   muon_Et;
  vector<float>   muon_Eta;
  vector<float>   muon_Phi;
  vector<float>   muon_EtaAtVtx;
  vector<float>   muon_PhiAtVtx;
  vector<short>   muon_IEt;
  vector<short>   muon_IEta;
  vector<short>   muon_IPhi;
  vector<short>   muon_IEtaAtVtx;
  vector<short>   muon_IPhiAtVtx;
  vector<short>   muon_IDEta;
  vector<short>   muon_IDPhi;
  vector<short>   muon_Chg;
  vector<unsigned short> muon_Iso;
  vector<unsigned short> muon_Qual;
  vector<unsigned short> muon_TfMuon_Idx;
  vector<short>   muon_Bx;

  UShort_t        sum_N;
  vector<short>   sum_Type;
  vector<float>   sum_Et;
  vector<float>   sum_Phi;
  vector<short>   sum_IEt;
  vector<short>   sum_IPhi;
  vector<float>   sum_Bx;



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
  UShort_t        egEmu_N;
  vector<float>   egEmu_Et;
  vector<float>   egEmu_Eta;
  vector<float>   egEmu_Phi;
  vector<short>   egEmu_IEt;
  vector<short>   egEmu_IEta;
  vector<short>   egEmu_IPhi;
  vector<short>   egEmu_Iso;
  vector<short>   egEmu_Bx;
  vector<short>   egEmu_TowerIPhi;
  vector<short>   egEmu_TowerIEta;
  vector<short>   egEmu_RawEt;
  vector<short>   egEmu_IsoEt;
  vector<short>   egEmu_FootprintEt;
  vector<short>   egEmu_NTT;
  vector<short>   egEmu_Shape;
  vector<short>   egEmu_TowerHoE;

  UShort_t        tauEmu_N;
  vector<float>   tauEmu_Et;
  vector<float>   tauEmu_Eta;
  vector<float>   tauEmu_Phi;
  vector<short>   tauEmu_IEt;
  vector<short>   tauEmu_IEta;
  vector<short>   tauEmu_IPhi;
  vector<short>   tauEmu_Iso;
  vector<short>   tauEmu_Bx;
  vector<short>   tauEmu_TowerIPhi;
  vector<short>   tauEmu_TowerIEta;
  vector<short>   tauEmu_RawEt;
  vector<short>   tauEmu_IsoEt;
  vector<short>   tauEmu_NTT;
  vector<short>   tauEmu_HasEM;
  vector<short>   tauEmu_IsMerged;
  vector<short>   tauEmu_HwQual;

  UShort_t        jetEmu_N;
  vector<float>   jetEmu_Et;
  vector<float>   jetEmu_Eta;
  vector<float>   jetEmu_Phi;
  vector<short>   jetEmu_IEt;
  vector<short>   jetEmu_IEta;
  vector<short>   jetEmu_IPhi;
  vector<short>   jetEmu_Bx;
  vector<short>   jetEmu_TowerIPhi;
  vector<short>   jetEmu_TowerIEta;
  vector<short>   jetEmu_RawEt;
  vector<short>   jetEmu_SeedEt;
  vector<short>   jetEmu_PUEt;
  vector<short>   jetEmu_PUDonutEt0;
  vector<short>   jetEmu_PUDonutEt1;
  vector<short>   jetEmu_PUDonutEt2;
  vector<short>   jetEmu_PUDonutEt3;

  UShort_t        muonEmu_N;
  vector<float>   muonEmu_Et;
  vector<float>   muonEmu_Eta;
  vector<float>   muonEmu_Phi;
  vector<float>   muonEmu_EtaAtVtx;
  vector<float>   muonEmu_PhiAtVtx;
  vector<short>   muonEmu_IEt;
  vector<short>   muonEmu_IEta;
  vector<short>   muonEmu_IPhi;
  vector<short>   muonEmu_IEtaAtVtx;
  vector<short>   muonEmu_IPhiAtVtx;
  vector<short>   muonEmu_IDEta;
  vector<short>   muonEmu_IDPhi;
  vector<short>   muonEmu_Chg;
  vector<unsigned short> muonEmu_Iso;
  vector<unsigned short> muonEmu_Qual;
  vector<unsigned short> muonEmu_TfMuon_Idx;
  vector<short>   muonEmu_Bx;

  UShort_t        sumEmu_N;
  vector<short>   sumEmu_Type;
  vector<float>   sumEmu_Et;
  vector<float>   sumEmu_Phi;
  vector<short>   sumEmu_IEt;
  vector<short>   sumEmu_IPhi;
  vector<float>   sumEmu_Bx;


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



  // Upgrade
  TBranch        *b_eg_N;   
  TBranch        *b_eg_Et;   
  TBranch        *b_eg_Eta;   
  TBranch        *b_eg_Phi;   
  TBranch        *b_eg_IEt;   
  TBranch        *b_eg_IEta;   
  TBranch        *b_eg_IPhi;   
  TBranch        *b_eg_Iso;   
  TBranch        *b_eg_Bx;   
  TBranch        *b_eg_TowerIPhi;   
  TBranch        *b_eg_TowerIEta;   
  TBranch        *b_eg_RawEt;   
  TBranch        *b_eg_IsoEt;   
  TBranch        *b_eg_FootprintEt;   
  TBranch        *b_eg_NTT;   
  TBranch        *b_eg_Shape;   
  TBranch        *b_eg_TowerHoE;   
  TBranch        *b_tau_N;   
  TBranch        *b_tau_Et;   
  TBranch        *b_tau_Eta;   
  TBranch        *b_tau_Phi;   
  TBranch        *b_tau_IEt;   
  TBranch        *b_tau_IEta;   
  TBranch        *b_tau_IPhi;   
  TBranch        *b_tau_Iso;   
  TBranch        *b_tau_Bx;   
  TBranch        *b_tau_TowerIPhi;   
  TBranch        *b_tau_TowerIEta;   
  TBranch        *b_tau_RawEt;   
  TBranch        *b_tau_IsoEt;   
  TBranch        *b_tau_NTT;   
  TBranch        *b_tau_HasEM;   
  TBranch        *b_tau_IsMerged;   
  TBranch        *b_tau_HwQual;   
  TBranch        *b_jet_N;   
  TBranch        *b_jet_Et;   
  TBranch        *b_jet_Eta;   
  TBranch        *b_jet_Phi;   
  TBranch        *b_jet_IEt;   
  TBranch        *b_jet_IEta;   
  TBranch        *b_jet_IPhi;   
  TBranch        *b_jet_Bx;   
  TBranch        *b_jet_TowerIPhi;   
  TBranch        *b_jet_TowerIEta;   
  TBranch        *b_jet_RawEt;   
  TBranch        *b_jet_SeedEt;   
  TBranch        *b_jet_PUEt;   
  TBranch        *b_jet_PUDonutEt0;   
  TBranch        *b_jet_PUDonutEt1;   
  TBranch        *b_jet_PUDonutEt2;   
  TBranch        *b_jet_PUDonutEt3;   
  TBranch        *b_muon_N;   
  TBranch        *b_muon_Et;   
  TBranch        *b_muon_Eta;   
  TBranch        *b_muon_Phi;   
  TBranch        *b_muon_EtaAtVtx;   
  TBranch        *b_muon_PhiAtVtx;   
  TBranch        *b_muon_IEt;   
  TBranch        *b_muon_IEta;   
  TBranch        *b_muon_IPhi;   
  TBranch        *b_muon_IEtaAtVtx;   
  TBranch        *b_muon_IPhiAtVtx;   
  TBranch        *b_muon_IDEta;   
  TBranch        *b_muon_IDPhi;   
  TBranch        *b_muon_Chg;   
  TBranch        *b_muon_Iso;   
  TBranch        *b_muon_Qual;   
  TBranch        *b_muon_TfMuon_Idx;   
  TBranch        *b_muon_Bx;   
  TBranch        *b_sum_N;   
  TBranch        *b_sum_Type;   
  TBranch        *b_sum_Et;   
  TBranch        *b_sum_Phi;   
  TBranch        *b_sum_IEt;   
  TBranch        *b_sum_IPhi;   
  TBranch        *b_sum_Bx;   


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
  TBranch        *b_egEmu_N;  
  TBranch        *b_egEmu_Et;   
  TBranch        *b_egEmu_Eta;   
  TBranch        *b_egEmu_Phi;   
  TBranch        *b_egEmu_IEt;   
  TBranch        *b_egEmu_IEta;   
  TBranch        *b_egEmu_IPhi;   
  TBranch        *b_egEmu_Iso;   
  TBranch        *b_egEmu_Bx;   
  TBranch        *b_egEmu_TowerIPhi;   
  TBranch        *b_egEmu_TowerIEta;   
  TBranch        *b_egEmu_RawEt;   
  TBranch        *b_egEmu_IsoEt;   
  TBranch        *b_egEmu_FootprintEt;   
  TBranch        *b_egEmu_NTT;   
  TBranch        *b_egEmu_Shape;   
  TBranch        *b_egEmu_TowerHoE;   
  TBranch        *b_tauEmu_N;   
  TBranch        *b_tauEmu_Et;   
  TBranch        *b_tauEmu_Eta;   
  TBranch        *b_tauEmu_Phi;   
  TBranch        *b_tauEmu_IEt;   
  TBranch        *b_tauEmu_IEta;   
  TBranch        *b_tauEmu_IPhi;   
  TBranch        *b_tauEmu_Iso;   
  TBranch        *b_tauEmu_Bx;   
  TBranch        *b_tauEmu_TowerIPhi;   
  TBranch        *b_tauEmu_TowerIEta;   
  TBranch        *b_tauEmu_RawEt;   
  TBranch        *b_tauEmu_IsoEt;   
  TBranch        *b_tauEmu_NTT;   
  TBranch        *b_tauEmu_HasEM;   
  TBranch        *b_tauEmu_IsMerged;   
  TBranch        *b_tauEmu_HwQual;   
  TBranch        *b_jetEmu_N;   
  TBranch        *b_jetEmu_Et;   
  TBranch        *b_jetEmu_Eta;   
  TBranch        *b_jetEmu_Phi;   
  TBranch        *b_jetEmu_IEt;   
  TBranch        *b_jetEmu_IEta;   
  TBranch        *b_jetEmu_IPhi;   
  TBranch        *b_jetEmu_Bx;   
  TBranch        *b_jetEmu_TowerIPhi;   
  TBranch        *b_jetEmu_TowerIEta;   
  TBranch        *b_jetEmu_RawEt;   
  TBranch        *b_jetEmu_SeedEt;   
  TBranch        *b_jetEmu_PUEt;   
  TBranch        *b_jetEmu_PUDonutEt0;   
  TBranch        *b_jetEmu_PUDonutEt1;   
  TBranch        *b_jetEmu_PUDonutEt2;   
  TBranch        *b_jetEmu_PUDonutEt3;   
  TBranch        *b_muonEmu_N;   
  TBranch        *b_muonEmu_Et;   
  TBranch        *b_muonEmu_Eta;   
  TBranch        *b_muonEmu_Phi;   
  TBranch        *b_muonEmu_EtaAtVtx;   
  TBranch        *b_muonEmu_PhiAtVtx;   
  TBranch        *b_muonEmu_IEt;   
  TBranch        *b_muonEmu_IEta;   
  TBranch        *b_muonEmu_IPhi;   
  TBranch        *b_muonEmu_IEtaAtVtx;   
  TBranch        *b_muonEmu_IPhiAtVtx;   
  TBranch        *b_muonEmu_IDEta;   
  TBranch        *b_muonEmu_IDPhi;   
  TBranch        *b_muonEmu_Chg;   
  TBranch        *b_muonEmu_Iso;   
  TBranch        *b_muonEmu_Qual;   
  TBranch        *b_muonEmu_TfMuon_Idx;   
  TBranch        *b_muonEmu_Bx;   
  TBranch        *b_sumEmu_N;   
  TBranch        *b_sumEmu_Type;   
  TBranch        *b_sumEmu_Et;   
  TBranch        *b_sumEmu_Phi;   
  TBranch        *b_sumEmu_IEt;   
  TBranch        *b_sumEmu_IPhi;   
  TBranch        *b_sumEmu_Bx;   


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

//void TreeDefinitionReco::InitReco(TTree *tree)
void TreeDefinitionReco::InitReco(TChain *chain)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  cout << "=== TreeDefinitionReco::InitReco()" << endl;
  

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
  eg_N    = 0;
  tau_N   = 0;
  jet_N   = 0;
  muon_N  = 0;
  sum_N   = 0;

  // uGT
  m_orbitNr        = 0;
  m_bxNr           = 0;
  m_bxInEvent      = 0;
  m_finalOR        = false; //CHECK (MARINA)
  m_finalORPreVeto = false; //CHECK (MARINA)
  m_finalORVeto    = false; //CHECK (MARINA)
  m_preScColumn    = 0;
  
  
  // HO
  hcalDetId_N     = 0;
  hcalQIESample_N = 0;

  // =========================================== EMULATOR TREES =============================================


  // Calo Towers Emulator
  hcalTPEmu_N      = 0;
  ecalTPEmu_N      = 0;
  ecalEBTPEmu_N    = 0;
  caloTowerEmu_N   = 0;
  caloClusterEmu_N = 0;


  // Upgrade Emulator
  egEmu_N    = 0;
  tauEmu_N   = 0;
  jetEmu_N   = 0;
  muonEmu_N  = 0;
  sumEmu_N   = 0;


  // uGT Emulator
  mEmu_orbitNr        = 0;
  mEmu_bxNr           = 0;
  mEmu_bxInEvent      = 0;
  mEmu_finalOR        = false; //CHECK (MARINA)
  mEmu_finalORPreVeto = false; //CHECK (MARINA)
  mEmu_finalORVeto    = false; //CHECK (MARINA)
  mEmu_preScColumn    = 0;



  cout << "\tSetting branch addresses and branch pointers." << endl;
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
      fUpgrade->SetBranchAddress("nEGs", &eg_N, &b_eg_N);
      fUpgrade->SetBranchAddress("egEt", &eg_Et, &b_eg_Et);
      fUpgrade->SetBranchAddress("egEta", &eg_Eta, &b_eg_Eta);
      fUpgrade->SetBranchAddress("egPhi", &eg_Phi, &b_eg_Phi);
      fUpgrade->SetBranchAddress("egIEt", &eg_IEt, &b_eg_IEt);
      fUpgrade->SetBranchAddress("egIEta", &eg_IEta, &b_eg_IEta);
      fUpgrade->SetBranchAddress("egIPhi", &eg_IPhi, &b_eg_IPhi);
      fUpgrade->SetBranchAddress("egIso", &eg_Iso, &b_eg_Iso);
      fUpgrade->SetBranchAddress("egBx", &eg_Bx, &b_eg_Bx);
      fUpgrade->SetBranchAddress("egTowerIPhi", &eg_TowerIPhi, &b_eg_TowerIPhi);
      fUpgrade->SetBranchAddress("egTowerIEta", &eg_TowerIEta, &b_eg_TowerIEta);
      fUpgrade->SetBranchAddress("egRawEt", &eg_RawEt, &b_eg_RawEt);
      fUpgrade->SetBranchAddress("egIsoEt", &eg_IsoEt, &b_eg_IsoEt);
      fUpgrade->SetBranchAddress("egFootprintEt", &eg_FootprintEt, &b_eg_FootprintEt);
      fUpgrade->SetBranchAddress("egNTT", &eg_NTT, &b_eg_NTT);
      fUpgrade->SetBranchAddress("egShape", &eg_Shape, &b_eg_Shape);
      fUpgrade->SetBranchAddress("egTowerHoE", &eg_TowerHoE, &b_eg_TowerHoE);

      fUpgrade->SetBranchAddress("nTaus", &tau_N, &b_tau_N);
      fUpgrade->SetBranchAddress("tauEt", &tau_Et, &b_tau_Et);
      fUpgrade->SetBranchAddress("tauEta", &tau_Eta, &b_tau_Eta);
      fUpgrade->SetBranchAddress("tauPhi", &tau_Phi, &b_tau_Phi);
      fUpgrade->SetBranchAddress("tauIEt", &tau_IEt, &b_tau_IEt);
      fUpgrade->SetBranchAddress("tauIEta", &tau_IEta, &b_tau_IEta);
      fUpgrade->SetBranchAddress("tauIPhi", &tau_IPhi, &b_tau_IPhi);
      fUpgrade->SetBranchAddress("tauIso", &tau_Iso, &b_tau_Iso);
      fUpgrade->SetBranchAddress("tauBx", &tau_Bx, &b_tau_Bx);
      fUpgrade->SetBranchAddress("tauTowerIPhi", &tau_TowerIPhi, &b_tau_TowerIPhi);
      fUpgrade->SetBranchAddress("tauTowerIEta", &tau_TowerIEta, &b_tau_TowerIEta);
      fUpgrade->SetBranchAddress("tauRawEt", &tau_RawEt, &b_tau_RawEt);
      fUpgrade->SetBranchAddress("tauIsoEt", &tau_IsoEt, &b_tau_IsoEt);
      fUpgrade->SetBranchAddress("tauNTT", &tau_NTT, &b_tau_NTT);
      fUpgrade->SetBranchAddress("tauHasEM", &tau_HasEM, &b_tau_HasEM);
      fUpgrade->SetBranchAddress("tauIsMerged", &tau_IsMerged, &b_tau_IsMerged);
      fUpgrade->SetBranchAddress("tauHwQual", &tau_HwQual, &b_tau_HwQual);

      fUpgrade->SetBranchAddress("nJets", &jet_N, &b_jet_N);
      fUpgrade->SetBranchAddress("jetEt", &jet_Et, &b_jet_Et);
      fUpgrade->SetBranchAddress("jetEta", &jet_Eta, &b_jet_Eta);
      fUpgrade->SetBranchAddress("jetPhi", &jet_Phi, &b_jet_Phi);
      fUpgrade->SetBranchAddress("jetIEt", &jet_IEt, &b_jet_IEt);
      fUpgrade->SetBranchAddress("jetIEta", &jet_IEta, &b_jet_IEta);
      fUpgrade->SetBranchAddress("jetIPhi", &jet_IPhi, &b_jet_IPhi);
      fUpgrade->SetBranchAddress("jetBx", &jet_Bx, &b_jet_Bx);
      fUpgrade->SetBranchAddress("jetTowerIPhi", &jet_TowerIPhi, &b_jet_TowerIPhi);
      fUpgrade->SetBranchAddress("jetTowerIEta", &jet_TowerIEta, &b_jet_TowerIEta);
      fUpgrade->SetBranchAddress("jetRawEt", &jet_RawEt, &b_jet_RawEt);
      fUpgrade->SetBranchAddress("jetSeedEt", &jet_SeedEt, &b_jet_SeedEt);
      fUpgrade->SetBranchAddress("jetPUEt", &jet_PUEt, &b_jet_PUEt);
      fUpgrade->SetBranchAddress("jetPUDonutEt0", &jet_PUDonutEt0, &b_jet_PUDonutEt0);
      fUpgrade->SetBranchAddress("jetPUDonutEt1", &jet_PUDonutEt1, &b_jet_PUDonutEt1);
      fUpgrade->SetBranchAddress("jetPUDonutEt2", &jet_PUDonutEt2, &b_jet_PUDonutEt2);
      fUpgrade->SetBranchAddress("jetPUDonutEt3", &jet_PUDonutEt3, &b_jet_PUDonutEt3);

      fUpgrade->SetBranchAddress("nMuons", &muon_N, &b_muon_N);
      fUpgrade->SetBranchAddress("muonEt", &muon_Et, &b_muon_Et);
      fUpgrade->SetBranchAddress("muonEta", &muon_Eta, &b_muon_Eta);
      fUpgrade->SetBranchAddress("muonPhi", &muon_Phi, &b_muon_Phi);
      fUpgrade->SetBranchAddress("muonEtaAtVtx", &muon_EtaAtVtx, &b_muon_EtaAtVtx);
      fUpgrade->SetBranchAddress("muonPhiAtVtx", &muon_PhiAtVtx, &b_muon_PhiAtVtx);
      fUpgrade->SetBranchAddress("muonIEt", &muon_IEt, &b_muon_IEt);
      fUpgrade->SetBranchAddress("muonIEta", &muon_IEta, &b_muon_IEta);
      fUpgrade->SetBranchAddress("muonIPhi", &muon_IPhi, &b_muon_IPhi);
      fUpgrade->SetBranchAddress("muonIEtaAtVtx", &muon_IEtaAtVtx, &b_muon_IEtaAtVtx);
      fUpgrade->SetBranchAddress("muonIPhiAtVtx", &muon_IPhiAtVtx, &b_muon_IPhiAtVtx);
      fUpgrade->SetBranchAddress("muonIDEta", &muon_IDEta, &b_muon_IDEta);
      fUpgrade->SetBranchAddress("muonIDPhi", &muon_IDPhi, &b_muon_IDPhi);
      fUpgrade->SetBranchAddress("muonChg", &muon_Chg, &b_muon_Chg);
      fUpgrade->SetBranchAddress("muonIso", &muon_Iso, &b_muon_Iso);
      fUpgrade->SetBranchAddress("muonQual", &muon_Qual, &b_muon_Qual);
      fUpgrade->SetBranchAddress("muonTfMuonIdx", &muon_TfMuon_Idx, &b_muon_TfMuon_Idx);
      fUpgrade->SetBranchAddress("muonBx", &muon_Bx, &b_muon_Bx);

      fUpgrade->SetBranchAddress("nSums", &sum_N, &b_sum_N);
      fUpgrade->SetBranchAddress("sumType", &sum_Type, &b_sum_Type);
      fUpgrade->SetBranchAddress("sumEt", &sum_Et, &b_sum_Et);
      fUpgrade->SetBranchAddress("sumPhi", &sum_Phi, &b_sum_Phi);
      fUpgrade->SetBranchAddress("sumIEt", &sum_IEt, &b_sum_IEt);
      fUpgrade->SetBranchAddress("sumIPhi", &sum_IPhi, &b_sum_IPhi);
      fUpgrade->SetBranchAddress("sumBx", &sum_Bx, &b_sum_Bx);
      

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
  
  
    
  // =========================================== EMULATOR TREES =============================================


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


      // Calo Clusters (Emu)
      fCaloTowerEmu->SetBranchAddress("ieta", &caloClusterEmu_ieta, &b_caloClusterEmu_ieta);
      fCaloTowerEmu->SetBranchAddress("iphi", &caloClusterEmu_iphi, &b_caloClusterEmu_iphi);
      fCaloTowerEmu->SetBranchAddress("iet", &caloClusterEmu_iet, &b_caloClusterEmu_iet);
      fCaloTowerEmu->SetBranchAddress("iqual", &caloClusterEmu_iqual, &b_caloClusterEmu_iqual);
      fCaloTowerEmu->SetBranchAddress("et", &caloClusterEmu_et, &b_caloClusterEmu_et);
      fCaloTowerEmu->SetBranchAddress("eta", &caloClusterEmu_eta, &b_caloClusterEmu_eta);
      fCaloTowerEmu->SetBranchAddress("phi", &caloClusterEmu_phi, &b_caloClusterEmu_phi);

      // Add Friend (CaloTowerEmuTree) 
      fChain -> AddFriend(fCaloTowerEmu);

    }



  // Upgrade Emulator
  if(doUpgradeEmu)
    {
      fUpgradeEmu->SetBranchAddress("nEGs", &egEmu_N, &b_egEmu_N);
      fUpgradeEmu->SetBranchAddress("egEt", &egEmu_Et, &b_egEmu_Et);
      fUpgradeEmu->SetBranchAddress("egEta", &egEmu_Eta, &b_egEmu_Eta);
      fUpgradeEmu->SetBranchAddress("egPhi", &egEmu_Phi, &b_egEmu_Phi);
      fUpgradeEmu->SetBranchAddress("egIEt", &egEmu_IEt, &b_egEmu_IEt);
      fUpgradeEmu->SetBranchAddress("egIEta", &egEmu_IEta, &b_egEmu_IEta);
      fUpgradeEmu->SetBranchAddress("egIPhi", &egEmu_IPhi, &b_egEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("egIso", &egEmu_Iso, &b_egEmu_Iso);
      fUpgradeEmu->SetBranchAddress("egBx", &egEmu_Bx, &b_egEmu_Bx);
      fUpgradeEmu->SetBranchAddress("egTowerIPhi", &egEmu_TowerIPhi, &b_egEmu_TowerIPhi);
      fUpgradeEmu->SetBranchAddress("egTowerIEta", &egEmu_TowerIEta, &b_egEmu_TowerIEta);
      fUpgradeEmu->SetBranchAddress("egRawEt", &egEmu_RawEt, &b_egEmu_RawEt);
      fUpgradeEmu->SetBranchAddress("egIsoEt", &egEmu_IsoEt, &b_egEmu_IsoEt);
      fUpgradeEmu->SetBranchAddress("egFootprintEt", &egEmu_FootprintEt, &b_egEmu_FootprintEt);
      fUpgradeEmu->SetBranchAddress("egNTT", &egEmu_NTT, &b_egEmu_NTT);
      fUpgradeEmu->SetBranchAddress("egShape", &egEmu_Shape, &b_egEmu_Shape);
      fUpgradeEmu->SetBranchAddress("egTowerHoE", &egEmu_TowerHoE, &b_egEmu_TowerHoE);

      fUpgradeEmu->SetBranchAddress("nTaus", &tauEmu_N, &b_tauEmu_N);
      fUpgradeEmu->SetBranchAddress("tauEt", &tauEmu_Et, &b_tauEmu_Et);
      fUpgradeEmu->SetBranchAddress("tauEta", &tauEmu_Eta, &b_tauEmu_Eta);
      fUpgradeEmu->SetBranchAddress("tauPhi", &tauEmu_Phi, &b_tauEmu_Phi);
      fUpgradeEmu->SetBranchAddress("tauIEt", &tauEmu_IEt, &b_tauEmu_IEt);
      fUpgradeEmu->SetBranchAddress("tauIEta", &tauEmu_IEta, &b_tauEmu_IEta);
      fUpgradeEmu->SetBranchAddress("tauIPhi", &tauEmu_IPhi, &b_tauEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("tauIso", &tauEmu_Iso, &b_tauEmu_Iso);
      fUpgradeEmu->SetBranchAddress("tauBx", &tauEmu_Bx, &b_tauEmu_Bx);
      fUpgradeEmu->SetBranchAddress("tauTowerIPhi", &tauEmu_TowerIPhi, &b_tauEmu_TowerIPhi);
      fUpgradeEmu->SetBranchAddress("tauTowerIEta", &tauEmu_TowerIEta, &b_tauEmu_TowerIEta);
      fUpgradeEmu->SetBranchAddress("tauRawEt", &tauEmu_RawEt, &b_tauEmu_RawEt);
      fUpgradeEmu->SetBranchAddress("tauIsoEt", &tauEmu_IsoEt, &b_tauEmu_IsoEt);
      fUpgradeEmu->SetBranchAddress("tauNTT", &tauEmu_NTT, &b_tauEmu_NTT);
      fUpgradeEmu->SetBranchAddress("tauHasEM", &tauEmu_HasEM, &b_tauEmu_HasEM);
      fUpgradeEmu->SetBranchAddress("tauIsMerged", &tauEmu_IsMerged, &b_tauEmu_IsMerged);
      fUpgradeEmu->SetBranchAddress("tauHwQual", &tauEmu_HwQual, &b_tauEmu_HwQual);

      fUpgradeEmu->SetBranchAddress("nJets", &jetEmu_N, &b_jetEmu_N);
      fUpgradeEmu->SetBranchAddress("jetEt", &jetEmu_Et, &b_jetEmu_Et);
      fUpgradeEmu->SetBranchAddress("jetEta", &jetEmu_Eta, &b_jetEmu_Eta);
      fUpgradeEmu->SetBranchAddress("jetPhi", &jetEmu_Phi, &b_jetEmu_Phi);
      fUpgradeEmu->SetBranchAddress("jetIEt", &jetEmu_IEt, &b_jetEmu_IEt);
      fUpgradeEmu->SetBranchAddress("jetIEta", &jetEmu_IEta, &b_jetEmu_IEta);
      fUpgradeEmu->SetBranchAddress("jetIPhi", &jetEmu_IPhi, &b_jetEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("jetBx", &jetEmu_Bx, &b_jetEmu_Bx);
      fUpgradeEmu->SetBranchAddress("jetTowerIPhi", &jetEmu_TowerIPhi, &b_jetEmu_TowerIPhi);
      fUpgradeEmu->SetBranchAddress("jetTowerIEta", &jetEmu_TowerIEta, &b_jetEmu_TowerIEta);
      fUpgradeEmu->SetBranchAddress("jetRawEt", &jetEmu_RawEt, &b_jetEmu_RawEt);
      fUpgradeEmu->SetBranchAddress("jetSeedEt", &jetEmu_SeedEt, &b_jetEmu_SeedEt);
      fUpgradeEmu->SetBranchAddress("jetPUEt", &jetEmu_PUEt, &b_jetEmu_PUEt);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt0", &jetEmu_PUDonutEt0, &b_jetEmu_PUDonutEt0);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt1", &jetEmu_PUDonutEt1, &b_jetEmu_PUDonutEt1);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt2", &jetEmu_PUDonutEt2, &b_jetEmu_PUDonutEt2);
      fUpgradeEmu->SetBranchAddress("jetPUDonutEt3", &jetEmu_PUDonutEt3, &b_jetEmu_PUDonutEt3);

      fUpgradeEmu->SetBranchAddress("nMuons", &muonEmu_N, &b_muonEmu_N);
      fUpgradeEmu->SetBranchAddress("muonEt", &muonEmu_Et, &b_muonEmu_Et);
      fUpgradeEmu->SetBranchAddress("muonEta", &muonEmu_Eta, &b_muonEmu_Eta);
      fUpgradeEmu->SetBranchAddress("muonPhi", &muonEmu_Phi, &b_muonEmu_Phi);
      fUpgradeEmu->SetBranchAddress("muonEtaAtVtx", &muonEmu_EtaAtVtx, &b_muonEmu_EtaAtVtx);
      fUpgradeEmu->SetBranchAddress("muonPhiAtVtx", &muonEmu_PhiAtVtx, &b_muonEmu_PhiAtVtx);
      fUpgradeEmu->SetBranchAddress("muonIEt", &muonEmu_IEt, &b_muonEmu_IEt);
      fUpgradeEmu->SetBranchAddress("muonIEta", &muonEmu_IEta, &b_muonEmu_IEta);
      fUpgradeEmu->SetBranchAddress("muonIPhi", &muonEmu_IPhi, &b_muonEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("muonIEtaAtVtx", &muonEmu_IEtaAtVtx, &b_muonEmu_IEtaAtVtx);
      fUpgradeEmu->SetBranchAddress("muonIPhiAtVtx", &muonEmu_IPhiAtVtx, &b_muonEmu_IPhiAtVtx);
      fUpgradeEmu->SetBranchAddress("muonIDEta", &muonEmu_IDEta, &b_muonEmu_IDEta);
      fUpgradeEmu->SetBranchAddress("muonIDPhi", &muonEmu_IDPhi, &b_muonEmu_IDPhi);
      fUpgradeEmu->SetBranchAddress("muonChg", &muonEmu_Chg, &b_muonEmu_Chg);
      fUpgradeEmu->SetBranchAddress("muonIso", &muonEmu_Iso, &b_muonEmu_Iso);
      fUpgradeEmu->SetBranchAddress("muonQual", &muonEmu_Qual, &b_muonEmu_Qual);
      fUpgradeEmu->SetBranchAddress("muonTfMuonIdx", &muonEmu_TfMuon_Idx, &b_muonEmu_TfMuon_Idx);
      fUpgradeEmu->SetBranchAddress("muonBx", &muonEmu_Bx, &b_muonEmu_Bx);

      fUpgradeEmu->SetBranchAddress("nSums", &sumEmu_N, &b_sumEmu_N);
      fUpgradeEmu->SetBranchAddress("sumType", &sumEmu_Type, &b_sumEmu_Type);
      fUpgradeEmu->SetBranchAddress("sumEt", &sumEmu_Et, &b_sumEmu_Et);
      fUpgradeEmu->SetBranchAddress("sumPhi", &sumEmu_Phi, &b_sumEmu_Phi);
      fUpgradeEmu->SetBranchAddress("sumIEt", &sumEmu_IEt, &b_sumEmu_IEt);
      fUpgradeEmu->SetBranchAddress("sumIPhi", &sumEmu_IPhi, &b_sumEmu_IPhi);
      fUpgradeEmu->SetBranchAddress("sumBx", &sumEmu_Bx, &b_sumEmu_Bx);

      // Add friend (UpgradeEmuTree)
      fChain -> AddFriend(fUpgradeEmu);

    }
  



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







  //cout << "**********************SETMAKECLASS-1**********************" << endl;
  
  fChain->SetMakeClass(1);
  
  //TList *friendTreeElement2 = fChain->GetListOfFriends();                                                                                                 
  //friendTreeElement2->ls();
  

  TFriendElement *friendTreeElement = (TFriendElement*)fChain->GetListOfFriends()->First();
  while (friendTreeElement){
    friendTreeElement->GetTree()->SetMakeClass(1);
    friendTreeElement = (TFriendElement*)fChain->GetListOfFriends()->After(friendTreeElement);
  }                                                                                                                                                                      

  //  cout << "**********************SETMAKECLASS-2**********************" << endl;




  
  Notify();
}

#endif // TreeDefinitionReco_h

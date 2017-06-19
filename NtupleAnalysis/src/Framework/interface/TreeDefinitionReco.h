#ifndef TreeDefinitionReco_h
#define TreeDefinitionReco_h

// User
#include "TreeDefinitionBase.h"
#include <TFriendElement.h>

using namespace std;

class TreeDefinitionReco : public virtual TreeDefinitionBase
{
 public:
  
  // TTTracks
  vector<double>  *L1Tks_Pt;
  vector<double>  *L1Tks_Eta;
  vector<double>  *L1Tks_Phi;
  vector<int>     *L1Tks_Charge;
  vector<double>  *L1Tks_POCAx;
  vector<double>  *L1Tks_POCAy;
  vector<double>  *L1Tks_POCAz;
  vector<double>  *L1Tks_ChiSquared;
  vector<double>  *L1Tks_StubPtConsistency;
  vector<double>  *L1Tks_RInv;
  vector<bool>    *L1Tks_IsGenuine;
  vector<bool>    *L1Tks_IsUnknown;
  vector<bool>    *L1Tks_IsCombinatoric;
  vector<bool>    *L1Tks_IsLoose;
  vector<bool>    *L1Tks_IsFake;
  vector<int>     *L1Tks_NStubs;
  vector<int>     *L1Tks_NStubsPS;
  vector<int>     *L1Tks_NStubsBarrel;
  vector<int>     *L1Tks_NStubsEndcap;
  vector<int>     *L1Tks_TP_Index;
  
  // TTPixelTracks
  vector<double> *L1PixTks_Pt;
  vector<double> *L1PixTks_Px;
  vector<double> *L1PixTks_Py;
  vector<double> *L1PixTks_Pz;
  vector<double> *L1PixTks_Eta;
  vector<double> *L1PixTks_Phi;
  vector<int>    *L1PixTks_Charge;
  vector<double> *L1PixTks_POCAx;
  vector<double> *L1PixTks_POCAy;
  vector<double> *L1PixTks_POCAz;
  vector<double> *L1PixTks_ChiSquared;
  vector<double> *L1PixTks_RInv;
  vector<double> *L1PixTks_SigmaRInv;
  vector<double> *L1PixTks_SigmaPhi0;
  vector<double> *L1PixTks_SigmaD0;
  vector<double> *L1PixTks_SigmaT;
  vector<double> *L1PixTks_SigmaZ0;
  vector<int>    *L1PixTks_TTTrackIndex;
  vector<int>    *L1PixTks_NPixHits;
  vector<vector<double> > *L1PixTks_PixHits_X;
  vector<vector<double> > *L1PixTks_PixHits_Y;
  vector<vector<double> > *L1PixTks_PixHits_Z;
  vector<vector<double> > *L1PixTks_PixHits_R;
  vector<vector<double> > *L1PixTks_PixHits_Phi;
  vector<vector<int> >    *L1PixTks_PixHits_Type;
  vector<vector<double> > *L1PixTks_CandPixHits_X;
  vector<vector<double> > *L1PixTks_CandPixHits_Y;
  vector<vector<double> > *L1PixTks_CandPixHits_Z;
  vector<vector<double> > *L1PixTks_CandPixHits_R;
  vector<vector<double> > *L1PixTks_CandPixHits_Phi;
  vector<vector<int> >    *L1PixTks_CandPixHits_Type;

  // Phase-1 L1T Stage2: Jets
  int *nJets;
  vector<float> *jetEt;
  vector<float> *jetEta;
  vector<float> *jetPhi;
  vector<int> *jetIEt;
  vector<int> *jetIEta;
  vector<int> *jetIPhi;
  vector<int> *jetBx;
  vector<int> *jetRawEt;
  vector<int> *jetSeedEt;
  vector<int> *jetTowerIEta;
  vector<int> *jetTowerIPhi;
  vector<int> *jetPUEt;
  vector<int> *jetPUDonutEt0;
  vector<int> *jetPUDonutEt1;
  vector<int> *jetPUDonutEt2;
  vector<int> *jetPUDonutEt3;

  
  // Phase-1 L1T Stage2: Taus
  int           *nTaus;
  vector<float> *tauEt;
  vector<float> *tauEta;
  vector<float> *tauPhi;
  vector<float> *tauIEt;
  vector<int>   *tauIEta;
  vector<int>   *tauIPhi;
  vector<int>   *tauIso;
  vector<int>   *tauBx;
  vector<int>   *tauTowerIPhi;
  vector<int>   *tauTowerIEta;
  vector<int>   *tauRawEt;
  vector<int>   *tauIsoEt;
  vector<int>   *tauNTT;
  vector<int>   *tauHasEM;
  vector<int>   *tauIsMerged;
  vector<int>   *tauHwQual;

  // Phase-1 L1T Stage2: EGamma
  int *nEGs;
  vector<float> *egEt;
  vector<float> *egEta;
  vector<float> *egPhi;
  vector<int> *egIEt;
  vector<int> *egIEta;
  vector<int> *egIPhi;
  vector<int> *egIso;
  vector<int> *egBx;
  vector<int> *egTowerIPhi;
  vector<int> *egTowerIEta;
  vector<int> *egRawEt;
  vector<int> *egIsoEt;
  vector<int> *egFootprintEt;
  vector<int> *egNTT;
  vector<int> *egShape;
  vector<int> *egTowerHoE;

  // Phase-1 L1T Stage2: Muons
  int *nMuons;
  vector<float> *muonEt;
  vector<float> *muonEta;
  vector<float> *muonPhi;
  vector<float> *muonEtaAtVtx;
  vector<float> *muonPhiAtVtx;
  vector<int> *muonIEt;
  vector<int> *muonIEta;
  vector<int> *muonIPhi;
  vector<int> *muonIDEta;
  vector<int> *muonIDPhi;
  vector<int> *muonChg;
  vector<int> *muonIso;
  vector<int> *muonHwQual;
  vector<int> *muonTfMuonIdx;
  vector<int> *muonBx;

  // L1 Sums
  int *nSums;
  vector<int> *sumType;
  vector<float> *sumEt;
  vector<float> *sumPhi;
  vector<int> *sumIEt;
  vector<int> *sumIPhi;
  vector<int> *sumBx;

  // L1TkJets
  vector< double > *L1TkJet_Pt;
  vector< double > *L1TkJet_Eta;
  vector< double > *L1TkJet_Phi;
  vector< double > *L1TkJet_E;
  vector< vector< int > > *L1TkJet_TTTrackIndex; 
  vector< double > *L1TkJet_Vertex;
  

  // TTTracks
  TBranch *b_L1Tks_Pt;
  TBranch *b_L1Tks_Eta;
  TBranch *b_L1Tks_Phi;
  TBranch *b_L1Tks_Charge;
  TBranch *b_L1Tks_POCAx;
  TBranch *b_L1Tks_POCAy;
  TBranch *b_L1Tks_POCAz;
  TBranch *b_L1Tks_ChiSquared;
  TBranch *b_L1Tks_StubPtConsistency;
  TBranch *b_L1Tks_RInv;
  TBranch *b_L1Tks_IsGenuine;
  TBranch *b_L1Tks_IsUnknown;
  TBranch *b_L1Tks_IsCombinatoric;
  TBranch *b_L1Tks_IsLoose;
  TBranch *b_L1Tks_IsFake;
  TBranch *b_L1Tks_NStubs;
  TBranch *b_L1Tks_NStubsPS;
  TBranch *b_L1Tks_NStubsBarrel;
  TBranch *b_L1Tks_NStubsEndcap;
  TBranch *b_L1Tks_TP_Index;

  // TTPixelTracks
  TBranch *b_L1PixTks_Pt;
  TBranch *b_L1PixTks_Px;
  TBranch *b_L1PixTks_Py;
  TBranch *b_L1PixTks_Pz;
  TBranch *b_L1PixTks_Eta;
  TBranch *b_L1PixTks_Phi;
  TBranch *b_L1PixTks_Charge;
  TBranch *b_L1PixTks_NPixHits;
  TBranch *b_L1PixTks_PixHits_X;
  TBranch *b_L1PixTks_PixHits_Y;
  TBranch *b_L1PixTks_PixHits_Z;
  TBranch *b_L1PixTks_PixHits_R;
  TBranch *b_L1PixTks_PixHits_Phi;
  TBranch *b_L1PixTks_PixHits_Type;
  TBranch *b_L1PixTks_CandPixHits_X;
  TBranch *b_L1PixTks_CandPixHits_Y;
  TBranch *b_L1PixTks_CandPixHits_Z;
  TBranch *b_L1PixTks_CandPixHits_R;
  TBranch *b_L1PixTks_CandPixHits_Phi;
  TBranch *b_L1PixTks_CandPixHits_Type;
  TBranch *b_L1PixTks_POCAx;
  TBranch *b_L1PixTks_POCAy;
  TBranch *b_L1PixTks_POCAz;
  TBranch *b_L1PixTks_ChiSquared;
  TBranch *b_L1PixTks_RInv;
  TBranch *b_L1PixTks_SigmaRInv;
  TBranch *b_L1PixTks_SigmaPhi0;
  TBranch *b_L1PixTks_SigmaD0;
  TBranch *b_L1PixTks_SigmaT;
  TBranch *b_L1PixTks_SigmaZ0;
  TBranch *b_L1PixTks_TTTrackIndex;

  // Phase-1 L1T Stage2: Jets
  TBranch *b_nJets;
  TBranch *b_jetEt;
  TBranch *b_jetEta;
  TBranch *b_jetPhi;
  TBranch *b_jetIEt;
  TBranch *b_jetIEta;
  TBranch *b_jetIPhi;
  TBranch *b_jetBx;
  TBranch *b_jetTowerIPhi;
  TBranch *b_jetTowerIEta;
  TBranch *b_jetRawEt;
  TBranch *b_jetSeedEt;
  TBranch *b_jetPUEt;
  TBranch *b_jetPUDonutEt0;
  TBranch *b_jetPUDonutEt1;
  TBranch *b_jetPUDonutEt2;
  TBranch *b_jetPUDonutEt3;

  // Phase-1 L1T Stage2: Taus
  TBranch *b_nTaus;
  TBranch *b_tauEt;
  TBranch *b_tauEta;
  TBranch *b_tauPhi;
  TBranch *b_tauIEt;
  TBranch *b_tauIEta;
  TBranch *b_tauIPhi;
  TBranch *b_tauIso;
  TBranch *b_tauBx;
  TBranch *b_tauTowerIPhi;
  TBranch *b_tauTowerIEta;
  TBranch *b_tauRawEt;
  TBranch *b_tauIsoEt;
  TBranch *b_tauNTT;
  TBranch *b_tauHasEM;
  TBranch *b_tauIsMerged;
  TBranch *b_tauHwQual;

  // Phase-1 L1T Stage2: EG (E-electron, G=Gamma)
  TBranch *b_nEGs;
  TBranch *b_egEt;
  TBranch *b_egEta;
  TBranch *b_egPhi;
  TBranch *b_egIEt;
  TBranch *b_egIEta;
  TBranch *b_egIPhi;
  TBranch *b_egIso;
  TBranch *b_egBx;
  TBranch *b_egTowerIPhi;
  TBranch *b_egTowerIEta;
  TBranch *b_egRawEt;
  TBranch *b_egIsoEt;
  TBranch *b_egFootprintEt;
  TBranch *b_egNTT;
  TBranch *b_egShape;
  TBranch *b_egTowerHoE;

  // Phase-1 L1T Stage2: Muons
  TBranch *b_nMuons;
  TBranch *b_muonEt;
  TBranch *b_muonEta;
  TBranch *b_muonPhi;
  TBranch *b_muonEtaAtVtx;
  TBranch *b_muonPhiAtVtx;
  TBranch *b_muonIEt;
  TBranch *b_muonIEta;
  TBranch *b_muonIPhi;
  TBranch *b_muonIDEta;
  TBranch *b_muonIDPhi;
  TBranch *b_muonChg;
  TBranch *b_muonIso;
  TBranch *b_muonHwQual;
  TBranch *b_muonTfMuonIdx;
  TBranch *b_muonBx;

  // L1 Sums
  TBranch *b_nSums;
  TBranch *b_sumType;
  TBranch *b_sumEt;
  TBranch *b_sumPhi;
  TBranch *b_sumIEt;
  TBranch *b_sumIPhi;
  TBranch *b_sumBx;
  
  // L1TkJets
  TBranch *b_L1TkJet_Pt;
  TBranch *b_L1TkJet_Eta;
  TBranch *b_L1TkJet_Phi;
  TBranch *b_L1TkJet_E;
  TBranch *b_L1TkJet_TTTrackIndex; 
  TBranch *b_L1TkJet_Vertex;
   

  virtual void InitReco(TTree *tree);
};

void TreeDefinitionReco::InitReco(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  cout << "=== TreeDefinitionReco::InitReco()" << endl;
  
  // TTTracks
  L1Tks_Pt                = 0;
  L1Tks_Eta               = 0;
  L1Tks_Phi               = 0;
  L1Tks_Charge            = 0;
  L1Tks_POCAx             = 0;
  L1Tks_POCAy             = 0;
  L1Tks_POCAz             = 0;
  L1Tks_ChiSquared        = 0;
  L1Tks_StubPtConsistency = 0;
  L1Tks_RInv              = 0;
  L1Tks_IsGenuine         = 0;
  L1Tks_IsUnknown         = 0;
  L1Tks_IsCombinatoric    = 0;
  L1Tks_IsLoose           = 0;
  L1Tks_IsFake            = 0;
  L1Tks_NStubs            = 0;
  L1Tks_NStubsPS          = 0;
  L1Tks_NStubsBarrel      = 0;
  L1Tks_NStubsEndcap      = 0;
  L1Tks_TP_Index          = 0;

  // TTPixelTracks
  L1PixTks_Pt                 = 0;
  L1PixTks_Px                 = 0;
  L1PixTks_Py                 = 0;
  L1PixTks_Pz                 = 0;
  L1PixTks_Eta                = 0;
  L1PixTks_Phi                = 0;
  L1PixTks_Charge             = 0;
  L1PixTks_POCAx              = 0;
  L1PixTks_POCAy              = 0;
  L1PixTks_POCAz              = 0;
  L1PixTks_ChiSquared         = 0;
  L1PixTks_RInv               = 0;
  L1PixTks_SigmaRInv          = 0;
  L1PixTks_SigmaPhi0          = 0;
  L1PixTks_SigmaD0            = 0;
  L1PixTks_SigmaT             = 0;
  L1PixTks_SigmaZ0            = 0;
  L1PixTks_TTTrackIndex       = 0;
  L1PixTks_NPixHits           = 0;
  L1PixTks_PixHits_X          = 0;
  L1PixTks_PixHits_Y          = 0;
  L1PixTks_PixHits_Z          = 0;
  L1PixTks_PixHits_R          = 0;
  L1PixTks_PixHits_Phi        = 0;
  L1PixTks_PixHits_Type       = 0;
  L1PixTks_CandPixHits_X      = 0;
  L1PixTks_CandPixHits_Y      = 0;
  L1PixTks_CandPixHits_Z      = 0;
  L1PixTks_CandPixHits_R      = 0;
  L1PixTks_CandPixHits_Phi    = 0;
  L1PixTks_CandPixHits_Type   = 0;

  // Phase-1 L1T Stage2: Jets
  nJets         = 0;
  jetEt         = 0;
  jetEta        = 0;
  jetPhi        = 0;
  jetIEt        = 0;
  jetIEta       = 0;
  jetIPhi       = 0;
  jetBx         = 0;
  jetTowerIPhi  = 0;
  jetTowerIEta  = 0;
  jetRawEt      = 0;
  jetSeedEt     = 0;
  jetPUEt       = 0;
  jetPUDonutEt0 = 0;
  jetPUDonutEt1 = 0;
  jetPUDonutEt2 = 0;
  jetPUDonutEt3 = 0;

  // Phase-1 L1T Stage2: Phase-1 L1T Stage2: Taus
  nTaus         = 0;
  tauEt         = 0;
  tauEta        = 0;
  tauPhi        = 0;
  tauIEt        = 0;
  tauIEta       = 0;
  tauIPhi       = 0;
  tauIso        = 0;
  tauBx         = 0;
  tauTowerIPhi  = 0;
  tauTowerIEta  = 0;
  tauRawEt      = 0;
  tauIsoEt      = 0;
  tauNTT        = 0;
  tauHasEM      = 0;
  tauIsMerged   = 0;
  tauHwQual     = 0;

  // Phase-1 L1T Stage2: EG (E-electron, G=Gamma)
  nEGs          = 0;
  egEt          = 0;
  egEta         = 0;
  egPhi         = 0;
  egIEt         = 0;
  egIEta        = 0;
  egIPhi        = 0;
  egIso         = 0;
  egBx          = 0;
  egTowerIPhi   = 0;
  egTowerIEta   = 0;
  egRawEt       = 0;
  egIsoEt       = 0;
  egFootprintEt = 0;
  egNTT         = 0;
  egShape       = 0;
  egTowerHoE    = 0;

  // Phase-1 L1T Stage2: Muons
  nMuons        = 0;
  muonEt        = 0;
  muonEta       = 0;
  muonPhi       = 0;
  muonEtaAtVtx  = 0;
  muonPhiAtVtx  = 0;
  muonIEt       = 0;
  muonIEta      = 0;
  muonIPhi      = 0;
  muonIDEta     = 0;
  muonIDPhi     = 0;
  muonChg       = 0;
  muonIso       = 0;
  muonHwQual    = 0;
  muonTfMuonIdx = 0;
  muonBx        = 0;

  // L1 Sums
  nSums   = 0;
  sumType = 0;
  sumEt   = 0;
  sumPhi  = 0;
  sumIEt  = 0;
  sumIPhi = 0;
  sumBx   = 0;
  
  // L1TkJets
  L1TkJet_Pt           = 0;
  L1TkJet_Eta          = 0;
  L1TkJet_Phi          = 0;
  L1TkJet_E            = 0;
  L1TkJet_TTTrackIndex = 0;
  L1TkJet_Vertex       = 0;


  cout << "\tSetting branch addresses and branch pointers." << endl;
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;

  fChain->SetMakeClass(1);
  //TFriendElement *friendTreeElement = (TFriendElement*)tree->GetListOfFriends()->Last();
  //friendTreeElement->GetTree()->SetMakeClass(1); //for some reason, friendChain->SetMakeClass(1); does not do the trick

  if(1)
    {
      cout << "\tSetting TTTracks addresses." << endl;
      fChain->SetBranchAddress("L1Tks_Pt"               , &L1Tks_Pt               , &b_L1Tks_Pt);
      fChain->SetBranchAddress("L1Tks_Eta"              , &L1Tks_Eta              , &b_L1Tks_Eta);
      fChain->SetBranchAddress("L1Tks_Phi"              , &L1Tks_Phi              , &b_L1Tks_Phi);
      fChain->SetBranchAddress("L1Tks_Charge"           , &L1Tks_Charge           , &b_L1Tks_Charge);
      fChain->SetBranchAddress("L1Tks_POCAx"            , &L1Tks_POCAx            , &b_L1Tks_POCAx);
      fChain->SetBranchAddress("L1Tks_POCAy"            , &L1Tks_POCAy            , &b_L1Tks_POCAy);
      fChain->SetBranchAddress("L1Tks_POCAz"            , &L1Tks_POCAz            , &b_L1Tks_POCAz);
      fChain->SetBranchAddress("L1Tks_ChiSquared"       , &L1Tks_ChiSquared       , &b_L1Tks_ChiSquared);
      fChain->SetBranchAddress("L1Tks_StubPtConsistency", &L1Tks_StubPtConsistency, &b_L1Tks_StubPtConsistency);
      fChain->SetBranchAddress("L1Tks_RInv"             , &L1Tks_RInv             , &b_L1Tks_RInv);
      fChain->SetBranchAddress("L1Tks_IsGenuine"        , &L1Tks_IsGenuine        , &b_L1Tks_IsGenuine);
      fChain->SetBranchAddress("L1Tks_IsUnknown"        , &L1Tks_IsUnknown        , &b_L1Tks_IsUnknown);
      fChain->SetBranchAddress("L1Tks_IsCombinatoric"   , &L1Tks_IsCombinatoric   , &b_L1Tks_IsCombinatoric);
      fChain->SetBranchAddress("L1Tks_IsLoose"          , &L1Tks_IsLoose          , &b_L1Tks_IsLoose);
      fChain->SetBranchAddress("L1Tks_IsFake"           , &L1Tks_IsFake           , &b_L1Tks_IsFake);
      fChain->SetBranchAddress("L1Tks_NStubs"           , &L1Tks_NStubs           , &b_L1Tks_NStubs);
      fChain->SetBranchAddress("L1Tks_NStubsPS"         , &L1Tks_NStubsPS         , &b_L1Tks_NStubsPS);
      fChain->SetBranchAddress("L1Tks_NStubsBarrel"     , &L1Tks_NStubsBarrel     , &b_L1Tks_NStubsBarrel);
      fChain->SetBranchAddress("L1Tks_NStubsEndcap"     , &L1Tks_NStubsEndcap     , &b_L1Tks_NStubsEndcap);
      fChain->SetBranchAddress("L1Tks_TP_Index"         , &L1Tks_TP_Index         , &b_L1Tks_TP_Index);
    }

  // TTPixelTracks 
  if (0)
    {
      cout << "\tSetting L1 Pixel Tracks addresses." << endl;
      fChain->SetBranchAddress("L1PixTks_Pt"              , &L1PixTks_Pt              , &b_L1PixTks_Pt);
      fChain->SetBranchAddress("L1PixTks_Px"              , &L1PixTks_Px              , &b_L1PixTks_Px);
      fChain->SetBranchAddress("L1PixTks_Py"              , &L1PixTks_Py              , &b_L1PixTks_Py);
      fChain->SetBranchAddress("L1PixTks_Pz"              , &L1PixTks_Pz              , &b_L1PixTks_Pz);
      fChain->SetBranchAddress("L1PixTks_Eta"             , &L1PixTks_Eta             , &b_L1PixTks_Eta);
      fChain->SetBranchAddress("L1PixTks_Phi"             , &L1PixTks_Phi             , &b_L1PixTks_Phi);
      fChain->SetBranchAddress("L1PixTks_Charge"          , &L1PixTks_Charge          , &b_L1PixTks_Charge);
      fChain->SetBranchAddress("L1PixTks_POCAx"           , &L1PixTks_POCAx           , &b_L1PixTks_POCAx);
      fChain->SetBranchAddress("L1PixTks_POCAy"           , &L1PixTks_POCAy           , &b_L1PixTks_POCAy);
      fChain->SetBranchAddress("L1PixTks_POCAz"           , &L1PixTks_POCAz           , &b_L1PixTks_POCAz);
      fChain->SetBranchAddress("L1PixTks_ChiSquared"      , &L1PixTks_ChiSquared      , &b_L1PixTks_ChiSquared);
      fChain->SetBranchAddress("L1PixTks_RInv"            , &L1PixTks_RInv            , &b_L1PixTks_RInv);
      fChain->SetBranchAddress("L1PixTks_SigmaRInv"       , &L1PixTks_SigmaRInv       , &b_L1PixTks_SigmaRInv);
      fChain->SetBranchAddress("L1PixTks_SigmaPhi0"       , &L1PixTks_SigmaPhi0       , &b_L1PixTks_SigmaPhi0);
      fChain->SetBranchAddress("L1PixTks_SigmaD0"         , &L1PixTks_SigmaD0         , &b_L1PixTks_SigmaD0);
      fChain->SetBranchAddress("L1PixTks_SigmaT"          , &L1PixTks_SigmaT          , &b_L1PixTks_SigmaT);
      fChain->SetBranchAddress("L1PixTks_SigmaZ0"         , &L1PixTks_SigmaZ0         , &b_L1PixTks_SigmaZ0);
      fChain->SetBranchAddress("L1PixTks_TTTrackIndex"    , &L1PixTks_TTTrackIndex    , &b_L1PixTks_TTTrackIndex);
      fChain->SetBranchAddress("L1PixTks_NPixHits"        , &L1PixTks_NPixHits        , &b_L1PixTks_NPixHits);
      fChain->SetBranchAddress("L1PixTks_PixHits_X"       , &L1PixTks_PixHits_X       , &b_L1PixTks_PixHits_X);
      fChain->SetBranchAddress("L1PixTks_PixHits_Y"       , &L1PixTks_PixHits_Y       , &b_L1PixTks_PixHits_Y);
      fChain->SetBranchAddress("L1PixTks_PixHits_Z"       , &L1PixTks_PixHits_Z       , &b_L1PixTks_PixHits_Z);
      fChain->SetBranchAddress("L1PixTks_PixHits_R"       , &L1PixTks_PixHits_R       , &b_L1PixTks_PixHits_R);
      fChain->SetBranchAddress("L1PixTks_PixHits_Phi"     , &L1PixTks_PixHits_Phi     , &b_L1PixTks_PixHits_Phi);
      fChain->SetBranchAddress("L1PixTks_PixHits_Type"    , &L1PixTks_PixHits_Type    , &b_L1PixTks_PixHits_Type);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_X"   , &L1PixTks_CandPixHits_X   , &b_L1PixTks_CandPixHits_X);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Y"   , &L1PixTks_CandPixHits_Y   , &b_L1PixTks_CandPixHits_Y);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Z"   , &L1PixTks_CandPixHits_Z   , &b_L1PixTks_CandPixHits_Z);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_R"   , &L1PixTks_CandPixHits_R   , &b_L1PixTks_CandPixHits_R);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Phi" , &L1PixTks_CandPixHits_Phi , &b_L1PixTks_CandPixHits_Phi);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Type", &L1PixTks_CandPixHits_Type, &b_L1PixTks_CandPixHits_Type);
    }

  // Phase-1 L1T Stage2: Jets
  if (1)
    {
      fChain->SetBranchAddress("L1Jet_nJets"     , &nJets        , &b_nJets);
      fChain->SetBranchAddress("L1Jet_Et"        , &jetEt        , &b_jetEt);
      fChain->SetBranchAddress("L1Jet_Eta"       , &jetEta       , &b_jetEta);
      fChain->SetBranchAddress("L1Jet_Phi"       , &jetPhi       , &b_jetPhi);
      fChain->SetBranchAddress("L1Jet_IET"       , &jetIEt       , &b_jetIEt);
      fChain->SetBranchAddress("L1Jet_IEta"      , &jetIEta      , &b_jetIEta);
      fChain->SetBranchAddress("L1Jet_IPhi"      , &jetIPhi      , &b_jetIPhi);
      fChain->SetBranchAddress("L1Jet_Bx"        , &jetBx        , &b_jetBx);
      fChain->SetBranchAddress("L1Jet_RawEt"     , &jetRawEt     , &b_jetRawEt);
      fChain->SetBranchAddress("L1Jet_SeedEt"    , &jetSeedEt    , &b_jetSeedEt);
      // fChain->SetBranchAddress("L1Jet_TowerIEta" , &jetTowerIEta , &b_jetTowerIEta);
      // fChain->SetBranchAddress("L1Jet_TowerIPhi" , &jetTowerIPhi , &b_jetTowerIPhi);
      fChain->SetBranchAddress("L1Jet_PUEt"      , &jetPUEt      , &b_jetPUEt);
      fChain->SetBranchAddress("L1Jet_PUDonutEt0", &jetPUDonutEt0, &b_jetPUDonutEt0);
      fChain->SetBranchAddress("L1Jet_PUDonutEt1", &jetPUDonutEt1, &b_jetPUDonutEt1);
      fChain->SetBranchAddress("L1Jet_PUDonutEt2", &jetPUDonutEt2, &b_jetPUDonutEt2);
      fChain->SetBranchAddress("L1Jet_PUDonutEt3", &jetPUDonutEt3, &b_jetPUDonutEt3);
    }
			   
  // Phase-1 L1T Stage2: Taus
  if (1)
    {
      cout << "\tSetting Phase-1 L1T Stage2 Tau addresses." << endl;
      fChain->SetBranchAddress("L1Tau_nTaus"    , &nTaus        , &b_nTaus);
      fChain->SetBranchAddress("L1Tau_Et"       , &tauEt        , &b_tauEt);
      fChain->SetBranchAddress("L1Tau_Eta"      , &tauEta       , &b_tauEta);
      fChain->SetBranchAddress("L1Tau_Phi"      , &tauPhi       , &b_tauPhi);
      fChain->SetBranchAddress("L1Tau_IET"      , &tauIEt       , &b_tauIEt);
      fChain->SetBranchAddress("L1Tau_IEta"     , &tauIEta      , &b_tauIEta);
      fChain->SetBranchAddress("L1Tau_IPhi"     , &tauIPhi      , &b_tauIPhi);
      fChain->SetBranchAddress("L1Tau_Iso"      , &tauIso       , &b_tauIso);
      fChain->SetBranchAddress("L1Tau_Bx"       , &tauBx        , &b_tauBx);
      // fChain->SetBranchAddress("L1Tau_TowerIPhi", &tauTowerIPhi , &b_tauTowerIPhi);
      // fChain->SetBranchAddress("L1Tau_TowerIEta", &tauTowerIEta , &b_tauTowerIEta);
      fChain->SetBranchAddress("L1Tau_RawEt"    , &tauRawEt     , &b_tauRawEt);
      fChain->SetBranchAddress("L1Tau_IsoEt"    , &tauIsoEt     , &b_tauIsoEt);
      fChain->SetBranchAddress("L1Tau_NTT"      , &tauNTT       , &b_tauNTT);
      fChain->SetBranchAddress("L1Tau_HasEM"    , &tauHasEM     , &b_tauHasEM);
      fChain->SetBranchAddress("L1Tau_IsMerged" , &tauIsMerged  , &b_tauIsMerged);
      fChain->SetBranchAddress("L1Tau_HwQual"   , &tauHwQual    , &b_tauHwQual);
    }

  // L1EG
  if (1)
    {
      cout << "\tSetting Phase-1 L1T Stage2 L1EG addresses." << endl;
      fChain->SetBranchAddress("L1EG_Et"         , &egEt         , &b_egEt);
      fChain->SetBranchAddress("L1EG_Eta"        , &egEta        , &b_egEta);
      fChain->SetBranchAddress("L1EG_Phi"        , &egPhi        , &b_egPhi);
      fChain->SetBranchAddress("L1EG_IET"        , &egIEt        , &b_egIEt);
      fChain->SetBranchAddress("L1EG_IEta"       , &egIEta       , &b_egIEta);
      fChain->SetBranchAddress("L1EG_IPhi"       , &egIPhi       , &b_egIPhi);
      fChain->SetBranchAddress("L1EG_Iso"        , &egIso        , &b_egIso);
      fChain->SetBranchAddress("L1EG_Bx"         , &egBx         , &b_egBx);
      // fChain->SetBranchAddress("L1EG_TowerIPhi"  , &egTowerIPhi  , &b_egTowerIPhi);
      // fChain->SetBranchAddress("L1EG_TowerIEta"  , &egTowerIEta  , &b_egTowerIEta);
      fChain->SetBranchAddress("L1EG_RawEt"      , &egRawEt      , &b_egRawEt);
      fChain->SetBranchAddress("L1EG_IsoEt"      , &egIsoEt      , &b_egIsoEt);
      fChain->SetBranchAddress("L1EG_FootprintEt", &egFootprintEt, &b_egFootprintEt);
      fChain->SetBranchAddress("L1EG_NTT"        , &egNTT        , &b_egNTT);
      fChain->SetBranchAddress("L1EG_Shape"      , &egShape      , &b_egShape);
      fChain->SetBranchAddress("L1EG_TowerHoE"   , &egTowerHoE   , &b_egTowerHoE);
      fChain->SetBranchAddress("L1EG_nEGs"       , &nEGs         , &b_nEGs);
    }
  
  // Muons
  if (0)
    {
      fChain->SetBranchAddress("L1Muon_Et"       , &muonEt       , &b_muonEt);
      fChain->SetBranchAddress("L1Muon_Eta"      , &muonEta      , &b_muonEta);
      fChain->SetBranchAddress("L1Muon_Phi"      , &muonPhi      , &b_muonPhi);
      fChain->SetBranchAddress("L1Muon_EtaAtVtx" , &muonEtaAtVtx , &b_muonEtaAtVtx);
      fChain->SetBranchAddress("L1Muon_PhiAtVtx" , &muonPhiAtVtx , &b_muonPhiAtVtx);
      fChain->SetBranchAddress("L1Muon_IET"      , &muonIEt      , &b_muonIEt);
      fChain->SetBranchAddress("L1Muon_IEta"     , &muonIEta     , &b_muonIEta);
      fChain->SetBranchAddress("L1Muon_IPhi"     , &muonIPhi     , &b_muonIPhi);
      fChain->SetBranchAddress("L1Muon_IDEta"    , &muonIDEta    , &b_muonIDEta);
      fChain->SetBranchAddress("L1Muon_IDPhi"    , &muonIDPhi    , &b_muonIDPhi);
      fChain->SetBranchAddress("L1Muon_Chg"      , &muonChg      , &b_muonChg);
      fChain->SetBranchAddress("L1Muon_Iso"      , &muonIso      , &b_muonIso);
      fChain->SetBranchAddress("L1Muon_HwQual"   , &muonHwQual   , &b_muonHwQual);
      fChain->SetBranchAddress("L1Muon_TfMuonIdx", &muonTfMuonIdx, &b_muonTfMuonIdx);
      fChain->SetBranchAddress("L1Muon_Bx"       , &muonBx       , &b_muonBx);
      fChain->SetBranchAddress("L1Muon_nMuons"   , &nMuons       , &b_nMuons);
    }

  // L1 Sums
  if (1)
    {
      fChain->SetBranchAddress("L1Sum_Type" , &sumType, &b_sumType);
      fChain->SetBranchAddress("L1Sum_Et"   , &sumEt  , &b_sumEt);
      fChain->SetBranchAddress("L1Sum_Phi"  , &sumPhi , &b_sumPhi);
      fChain->SetBranchAddress("L1Sum_IET"  , &sumIEt , &b_sumIEt);
      fChain->SetBranchAddress("L1Sum_IPhi" , &sumIPhi, &b_sumIPhi);
      fChain->SetBranchAddress("L1Sum_Bx"   , &sumBx  , &b_sumBx);
      fChain->SetBranchAddress("L1Sum_nSums", &nSums  , &b_nSums);
    }

  
  // L1TkJets
  if (0)
    {
      cout << "\tSetting L1 Jet addresses." << endl;
      fChain->SetBranchAddress("L1TkJet_Pt"          , &L1TkJet_Pt          , &b_L1TkJet_Pt);
      fChain->SetBranchAddress("L1TkJet_Eta"         , &L1TkJet_Eta         , &b_L1TkJet_Eta);
      fChain->SetBranchAddress("L1TkJet_Phi"         , &L1TkJet_Phi         , &b_L1TkJet_Phi);
      fChain->SetBranchAddress("L1TkJet_E"           , &L1TkJet_E           , &b_L1TkJet_E);
      fChain->SetBranchAddress("L1TkJet_TTTrackIndex", &L1TkJet_TTTrackIndex, &b_L1TkJet_TTTrackIndex);
      fChain->SetBranchAddress("L1TkJet_Vertex"      , &L1TkJet_Vertex      , &b_L1TkJet_Vertex);
    }
  
  Notify();
}

#endif // TreeDefinitionReco_h

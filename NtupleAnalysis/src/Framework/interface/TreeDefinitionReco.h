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
  vector<double>  *L1Tks_Charge;
  vector<double>  *L1Tks_POCAx;
  vector<double>  *L1Tks_POCAy;
  vector<double>  *L1Tks_POCAz;
  vector<double>  *L1Tks_ChiSquared;
  vector<double>  *L1Tks_StubPtConsistency;
  vector<double>  *L1Tks_RInv;
  vector<int>     *L1Tks_IsGenuine;
  vector<int>     *L1Tks_IsUnknown;
  vector<int>     *L1Tks_IsCombinatoric;
  vector<int>     *L1Tks_IsLoose;
  vector<int>     *L1Tks_IsFake;
  vector<int>     *L1Tks_NStubs;
  vector<int>     *L1Tks_NStubsPS;
  vector<int>     *L1Tks_NStubsBarrel;
  vector<int>     *L1Tks_NStubsEndcap;
  vector<int>     *L1Tks_TP_Index;

  // Phase-1 L1T Stage2: Jets
  short *L1Jet_nJets;
  vector<float>   *L1Jet_Et;
  vector<float>   *L1Jet_Eta;
  vector<float>   *L1Jet_Phi;
  vector<short>   *L1Jet_IET;
  vector<short>   *L1Jet_IEta;
  vector<short>   *L1Jet_IPhi;
  vector<short>   *L1Jet_Bx;
  vector<short>   *L1Jet_RawEt;
  vector<short>   *L1Jet_SeedEt;
  vector<short>   *L1Jet_TowerIEta;
  vector<short>   *L1Jet_TowerIPhi;
  vector<short>   *L1Jet_PUEt;
  vector<short>   *L1Jet_PUDonutEt0;
  vector<short>   *L1Jet_PUDonutEt1;
  vector<short>   *L1Jet_PUDonutEt2;
  vector<short>   *L1Jet_PUDonutEt3;
   
  // Phase-1 L1T Stage2: Taus
  short *L1Tau_nTaus;
  vector<float>   *L1Tau_Et;
  vector<float>   *L1Tau_Eta;
  vector<float>   *L1Tau_Phi;
  vector<short>   *L1Tau_IET;
  vector<short>   *L1Tau_IEta;
  vector<short>   *L1Tau_IPhi;
  vector<short>   *L1Tau_Iso;
  vector<short>   *L1Tau_Bx;
  vector<short>   *L1Tau_TowerIEta;
  vector<short>   *L1Tau_TowerIPhi;
  vector<short>   *L1Tau_RawEt;
  vector<short>   *L1Tau_IsoEt;
  vector<short>   *L1Tau_NTT;
  vector<short>   *L1Tau_HasEM;
  vector<short>   *L1Tau_IsMerged;
  vector<short>   *L1Tau_HwQual;

  // Phase-1 L1T Stage2: EGamma
  short *L1EG_nEGs;
  vector<float>   *L1EG_Et;
  vector<float>   *L1EG_Eta;
  vector<float>   *L1EG_Phi;
  vector<short>   *L1EG_IET;
  vector<short>   *L1EG_IEta;
  vector<short>   *L1EG_IPhi;
  vector<short>   *L1EG_Iso;
  vector<short>   *L1EG_Bx;
  vector<short>   *L1EG_TowerIPhi;
  vector<short>   *L1EG_TowerIEta;
  vector<short>   *L1EG_RawEt;
  vector<short>   *L1EG_IsoEt;
  vector<short>   *L1EG_FootprintEt;
  vector<short>   *L1EG_NTT;
  vector<short>   *L1EG_Shape;
  vector<short>   *L1EG_TowerHoE;

  // L1 Sums
  short *L1Sum_nSums;
  vector<short> *L1Sum_Type;
  vector<float> *L1Sum_Et;
  vector<float> *L1Sum_Phi;
  vector<short> *L1Sum_IET;
  vector<short> *L1Sum_IPhi;
  vector<float> *L1Sum_Bx;

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

  // Phase-1 L1T Stage2: Jets
  TBranch *b_L1Jet_nJets;
  TBranch *b_L1Jet_Et;
  TBranch *b_L1Jet_Eta;
  TBranch *b_L1Jet_Phi;
  TBranch *b_L1Jet_IET;
  TBranch *b_L1Jet_IEta;
  TBranch *b_L1Jet_IPhi;
  TBranch *b_L1Jet_Bx;
  TBranch *b_L1Jet_TowerIEta;
  TBranch *b_L1Jet_TowerIPhi;
  TBranch *b_L1Jet_RawEt;
  TBranch *b_L1Jet_SeedEt;
  TBranch *b_L1Jet_PUEt;
  TBranch *b_L1Jet_PUDonutEt0;
  TBranch *b_L1Jet_PUDonutEt1;
  TBranch *b_L1Jet_PUDonutEt2;
  TBranch *b_L1Jet_PUDonutEt3;

  // Phase-1 L1T Stage2: Taus
  TBranch *b_L1Tau_nTaus;
  TBranch *b_L1Tau_Et;
  TBranch *b_L1Tau_Eta;
  TBranch *b_L1Tau_Phi;
  TBranch *b_L1Tau_IET;
  TBranch *b_L1Tau_IEta;
  TBranch *b_L1Tau_IPhi;
  TBranch *b_L1Tau_Iso;
  TBranch *b_L1Tau_Bx;
  TBranch *b_L1Tau_TowerIEta;
  TBranch *b_L1Tau_TowerIPhi;
  TBranch *b_L1Tau_RawEt;
  TBranch *b_L1Tau_IsoEt;
  TBranch *b_L1Tau_NTT;
  TBranch *b_L1Tau_HasEM;
  TBranch *b_L1Tau_IsMerged;
  TBranch *b_L1Tau_HwQual;

  // Phase-1 L1T Stage2: EG (E-electron, G=Gamma)
  TBranch *b_L1EG_nEGs;
  TBranch *b_L1EG_Et;
  TBranch *b_L1EG_Eta;
  TBranch *b_L1EG_Phi;
  TBranch *b_L1EG_IET;
  TBranch *b_L1EG_IEta;
  TBranch *b_L1EG_IPhi;
  TBranch *b_L1EG_Iso;
  TBranch *b_L1EG_Bx;
  TBranch *b_L1EG_TowerIEta;
  TBranch *b_L1EG_TowerIPhi;
  TBranch *b_L1EG_RawEt;
  TBranch *b_L1EG_IsoEt;
  TBranch *b_L1EG_FootprintEt;
  TBranch *b_L1EG_NTT;
  TBranch *b_L1EG_Shape;
  TBranch *b_L1EG_TowerHoE;

  // L1 Sums
  TBranch *b_L1Sum_nSums;
  TBranch *b_L1Sum_Type;
  TBranch *b_L1Sum_Et;
  TBranch *b_L1Sum_Phi;
  TBranch *b_L1Sum_IET;
  TBranch *b_L1Sum_IPhi;
  TBranch *b_L1Sum_Bx;
     

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

  // Phase-1 L1T Stage2: Jets
  L1Jet_nJets      = 0;
  L1Jet_Et         = 0;
  L1Jet_Eta        = 0;
  L1Jet_Phi        = 0;
  L1Jet_IET        = 0;
  L1Jet_IEta       = 0;
  L1Jet_IPhi       = 0;
  L1Jet_Bx         = 0;
  L1Jet_TowerIEta  = 0;
  L1Jet_TowerIPhi  = 0;
  L1Jet_RawEt      = 0;
  L1Jet_SeedEt     = 0;
  L1Jet_PUEt       = 0;
  L1Jet_PUDonutEt0 = 0;
  L1Jet_PUDonutEt1 = 0;
  L1Jet_PUDonutEt2 = 0;
  L1Jet_PUDonutEt3 = 0;

  // Phase-1 L1T Stage2: Phase-1 L1T Stage2: Taus
  L1Tau_nTaus     = 0;
  L1Tau_Et        = 0;
  L1Tau_Eta       = 0;
  L1Tau_Phi       = 0;
  L1Tau_IET       = 0;
  L1Tau_IEta      = 0;
  L1Tau_IPhi      = 0;
  L1Tau_Iso       = 0;
  L1Tau_Bx        = 0;
  L1Tau_TowerIEta = 0;
  L1Tau_TowerIPhi = 0;
  L1Tau_RawEt     = 0;
  L1Tau_IsoEt     = 0;
  L1Tau_NTT       = 0;
  L1Tau_HasEM     = 0;
  L1Tau_IsMerged  = 0;
  L1Tau_HwQual    = 0;

  // Phase-1 L1T Stage2: EG (E-electron, G=Gamma)
  L1EG_nEGs        = 0;
  L1EG_Et          = 0;
  L1EG_Eta         = 0;
  L1EG_Phi         = 0;
  L1EG_IET         = 0;
  L1EG_IEta        = 0;
  L1EG_IPhi        = 0;
  L1EG_Iso         = 0;
  L1EG_Bx          = 0;
  L1EG_TowerIEta   = 0;
  L1EG_TowerIPhi   = 0;
  L1EG_RawEt       = 0;
  L1EG_IsoEt       = 0;
  L1EG_FootprintEt = 0;
  L1EG_NTT         = 0;
  L1EG_Shape       = 0;
  L1EG_TowerHoE    = 0;

  // L1 Sums
  L1Sum_nSums= 0;
  L1Sum_Type = 0;
  L1Sum_Et   = 0;
  L1Sum_Phi  = 0;
  L1Sum_IET  = 0;
  L1Sum_IPhi = 0;
  L1Sum_Bx   = 0;

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

  // Phase-1 L1T Stage2: Jets
  if (1)
    {
      fChain->SetBranchAddress("L1Jet_nJets"     , &L1Jet_nJets     , &b_L1Jet_nJets);
      fChain->SetBranchAddress("L1Jet_Et"        , &L1Jet_Et        , &b_L1Jet_Et);
      fChain->SetBranchAddress("L1Jet_Eta"       , &L1Jet_Eta       , &b_L1Jet_Eta);
      fChain->SetBranchAddress("L1Jet_Phi"       , &L1Jet_Phi       , &b_L1Jet_Phi);
      fChain->SetBranchAddress("L1Jet_IET"       , &L1Jet_IET       , &b_L1Jet_IET);
      fChain->SetBranchAddress("L1Jet_IEta"      , &L1Jet_IEta      , &b_L1Jet_IEta);
      fChain->SetBranchAddress("L1Jet_IPhi"      , &L1Jet_IPhi      , &b_L1Jet_IPhi);
      fChain->SetBranchAddress("L1Jet_Bx"        , &L1Jet_Bx        , &b_L1Jet_Bx);
      fChain->SetBranchAddress("L1Jet_RawEt"     , &L1Jet_RawEt     , &b_L1Jet_RawEt);
      fChain->SetBranchAddress("L1Jet_SeedEt"    , &L1Jet_SeedEt    , &b_L1Jet_SeedEt);
      // fChain->SetBranchAddress("L1Jet_TowerIEta" , &L1JetTower_IEta , &b_L1Jet_TowerIEta);
      // fChain->SetBranchAddress("L1Jet_TowerIPhi" , &L1JetTower_IPhi , &b_L1Jet_TowerIPhi);
      fChain->SetBranchAddress("L1Jet_PUEt"      , &L1Jet_PUEt      , & b_L1Jet_PUEt);
      fChain->SetBranchAddress("L1Jet_PUDonutEt0", &L1Jet_PUDonutEt0, & b_L1Jet_PUDonutEt0);
      fChain->SetBranchAddress("L1Jet_PUDonutEt1", &L1Jet_PUDonutEt1, & b_L1Jet_PUDonutEt1);
      fChain->SetBranchAddress("L1Jet_PUDonutEt2", &L1Jet_PUDonutEt2, & b_L1Jet_PUDonutEt2);
      fChain->SetBranchAddress("L1Jet_PUDonutEt3", &L1Jet_PUDonutEt3, & b_L1Jet_PUDonutEt3);
    }
			   
  // Phase-1 L1T Stage2: Taus
  if (1)
    {
      cout << "\tSetting Phase-1 L1T Stage2 Tau addresses." << endl;
      fChain->SetBranchAddress("L1Tau_nTaus"    , &L1Tau_nTaus     , &b_L1Tau_nTaus);
      fChain->SetBranchAddress("L1Tau_Et"       , &L1Tau_Et        , &b_L1Tau_Et);
      fChain->SetBranchAddress("L1Tau_Eta"      , &L1Tau_Eta       , &b_L1Tau_Eta);
      fChain->SetBranchAddress("L1Tau_Phi"      , &L1Tau_Phi       , &b_L1Tau_Phi);
      fChain->SetBranchAddress("L1Tau_IET"      , &L1Tau_IET       , &b_L1Tau_IET);
      fChain->SetBranchAddress("L1Tau_IEta"     , &L1Tau_IEta      , &b_L1Tau_IEta);
      fChain->SetBranchAddress("L1Tau_IPhi"     , &L1Tau_IPhi      , &b_L1Tau_IPhi);
      fChain->SetBranchAddress("L1Tau_Iso"      , &L1Tau_Iso       , &b_L1Tau_Iso);
      fChain->SetBranchAddress("L1Tau_Bx"       , &L1Tau_Bx        , &b_L1Tau_Bx);
      // fChain->SetBranchAddress("L1Tau_TowerIEta", &L1Tau_TowerIEta , &b_L1Tau_TowerIEta);
      // fChain->SetBranchAddress("L1Tau_TowerIPhi", &L1Tau_TowerIPhi , &b_L1Tau_TowerIPhi);
      fChain->SetBranchAddress("L1Tau_RawEt"    , &L1Tau_RawEt     , &b_L1Tau_RawEt);
      fChain->SetBranchAddress("L1Tau_IsoEt"    , &L1Tau_IsoEt     , &b_L1Tau_IsoEt);
      fChain->SetBranchAddress("L1Tau_NTT"      , &L1Tau_NTT       , &b_L1Tau_NTT);
      fChain->SetBranchAddress("L1Tau_HasEM"    , &L1Tau_HasEM     , &b_L1Tau_HasEM);
      fChain->SetBranchAddress("L1Tau_IsMerged" , &L1Tau_IsMerged  , &b_L1Tau_IsMerged);
      fChain->SetBranchAddress("L1Tau_HwQual"   , &L1Tau_HwQual    , &b_L1Tau_HwQual);
    }

  // L1EG
  if (1)
    {
      cout << "\tSetting Phase-1 L1T Stage2 L1EG addresses." << endl;
      fChain->SetBranchAddress("L1EG_nEGs"       , &L1EG_nEGs       , &b_L1EG_nEGs);
      fChain->SetBranchAddress("L1EG_Et"         , &L1EG_Et         , &b_L1EG_Et);
      fChain->SetBranchAddress("L1EG_Eta"        , &L1EG_Eta        , &b_L1EG_Eta);
      fChain->SetBranchAddress("L1EG_Phi"        , &L1EG_Phi        , &b_L1EG_Phi);
      fChain->SetBranchAddress("L1EG_IET"        , &L1EG_IET        , &b_L1EG_IET);
      fChain->SetBranchAddress("L1EG_IEta"       , &L1EG_IEta       , &b_L1EG_IEta);
      fChain->SetBranchAddress("L1EG_IPhi"       , &L1EG_IPhi       , &b_L1EG_IPhi);
      fChain->SetBranchAddress("L1EG_Iso"        , &L1EG_Iso        , &b_L1EG_Iso);
      fChain->SetBranchAddress("L1EG_Bx"         , &L1EG_Bx         , &b_L1EG_Bx);
      // fChain->SetBranchAddress("L1EG_TowerIEta"  , &L1EG_TowerIEta  , &b_L1EG_TowerIEta);
      // fChain->SetBranchAddress("L1EG_TowerIPhi"  , &L1EG_TowerIPhi  , &b_L1EG_TowerIPhi);
      fChain->SetBranchAddress("L1EG_RawEt"      , &L1EG_RawEt      , &b_L1EG_RawEt);
      fChain->SetBranchAddress("L1EG_IsoEt"      , &L1EG_IsoEt      , &b_L1EG_IsoEt);
      fChain->SetBranchAddress("L1EG_FootprintEt", &L1EG_FootprintEt, &b_L1EG_FootprintEt);
      fChain->SetBranchAddress("L1EG_NTT"        , &L1EG_NTT        , &b_L1EG_NTT);
      fChain->SetBranchAddress("L1EG_Shape"      , &L1EG_Shape      , &b_L1EG_Shape);
      fChain->SetBranchAddress("L1EG_TowerHoE"   , &L1EG_TowerHoE   , &b_L1EG_TowerHoE);
    }
  
  // L1 Sums
  if (1)
    {
      fChain->SetBranchAddress("L1Sum_nSums", &L1Sum_nSums, &b_L1Sum_nSums);
      fChain->SetBranchAddress("L1Sum_Type" , &L1Sum_Type , &b_L1Sum_Type);
      fChain->SetBranchAddress("L1Sum_Et"   , &L1Sum_Et   , &b_L1Sum_Et);
      fChain->SetBranchAddress("L1Sum_Phi"  , &L1Sum_Phi  , &b_L1Sum_Phi);
      fChain->SetBranchAddress("L1Sum_IET"  , &L1Sum_IET  , &b_L1Sum_IET);
      fChain->SetBranchAddress("L1Sum_IPhi" , &L1Sum_IPhi , &b_L1Sum_IPhi);
      fChain->SetBranchAddress("L1Sum_Bx"   , &L1Sum_Bx   , &b_L1Sum_Bx);
    }

  
  Notify();
}

#endif // TreeDefinitionReco_h

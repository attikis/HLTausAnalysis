#ifndef TreeDefinitionReco_h
#define TreeDefinitionReco_h

// User
#include "TreeDefinitionBase.h"

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

  // L1CaloTaus
  vector<double>  *L1CaloTau_E;
  vector<double>  *L1CaloTau_Et;
  vector<double>  *L1CaloTau_Eta;
  vector<double>  *L1CaloTau_Phi;
  vector<int>     *L1CaloTau_Bx;
  vector<int>     *L1CaloTau_Type;

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

  // L1CaloTaus
  TBranch *b_L1CaloTau_E;
  TBranch *b_L1CaloTau_Et;
  TBranch *b_L1CaloTau_Eta;
  TBranch *b_L1CaloTau_Phi;
  TBranch *b_L1CaloTau_Bx;
  TBranch *b_L1CaloTau_Type;

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

  // L1CaloTaus
  L1CaloTau_E             = 0;
  L1CaloTau_Et            = 0;
  L1CaloTau_Eta           = 0;
  L1CaloTau_Phi           = 0;
  L1CaloTau_Bx            = 0;
  L1CaloTau_Type          = 0;

  // L1TkJets
  L1TkJet_Pt              = 0;
  L1TkJet_Eta             = 0;
  L1TkJet_Phi             = 0;
  L1TkJet_E               = 0;
  L1TkJet_TTTrackIndex    = 0;
  L1TkJet_Vertex          = 0;


  cout << "\tSetting branch addresses and branch pointers." << endl;
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1); 

  // TTTracks
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
      fChain->SetBranchAddress("L1Tks_IsGenuine"     , &L1Tks_IsGenuine     , &b_L1Tks_IsGenuine);
      fChain->SetBranchAddress("L1Tks_IsUnknown"     , &L1Tks_IsUnknown     , &b_L1Tks_IsUnknown);
      fChain->SetBranchAddress("L1Tks_IsCombinatoric", &L1Tks_IsCombinatoric, &b_L1Tks_IsCombinatoric);
      fChain->SetBranchAddress("L1Tks_IsLoose"       , &L1Tks_IsLoose       , &b_L1Tks_IsLoose);
      fChain->SetBranchAddress("L1Tks_IsFake"        , &L1Tks_IsFake        , &b_L1Tks_IsFake);
      fChain->SetBranchAddress("L1Tks_NStubs"        , &L1Tks_NStubs        , &b_L1Tks_NStubs);
      fChain->SetBranchAddress("L1Tks_NStubsPS"      , &L1Tks_NStubsPS      , &b_L1Tks_NStubsPS);
      fChain->SetBranchAddress("L1Tks_NStubsBarrel"  , &L1Tks_NStubsBarrel  , &b_L1Tks_NStubsBarrel);
      fChain->SetBranchAddress("L1Tks_NStubsEndcap"  , &L1Tks_NStubsEndcap  , &b_L1Tks_NStubsEndcap);
      fChain->SetBranchAddress("L1Tks_TP_Index"      , &L1Tks_TP_Index      , &b_L1Tks_TP_Index);
    }

  // TTPixelTracks 
  if (0)
    {
      cout << "\tSetting L1 Pixel Tracks addresses." << endl;
      fChain->SetBranchAddress("L1PixTks_Pt"           , &L1PixTks_Pt           , &b_L1PixTks_Pt);
      fChain->SetBranchAddress("L1PixTks_Px"           , &L1PixTks_Px           , &b_L1PixTks_Px);
      fChain->SetBranchAddress("L1PixTks_Py"           , &L1PixTks_Py           , &b_L1PixTks_Py);
      fChain->SetBranchAddress("L1PixTks_Pz"           , &L1PixTks_Pz           , &b_L1PixTks_Pz);
      fChain->SetBranchAddress("L1PixTks_Eta"          , &L1PixTks_Eta          , &b_L1PixTks_Eta);
      fChain->SetBranchAddress("L1PixTks_Phi"          , &L1PixTks_Phi          , &b_L1PixTks_Phi);
      fChain->SetBranchAddress("L1PixTks_Charge"       , &L1PixTks_Charge       , &b_L1PixTks_Charge);
      fChain->SetBranchAddress("L1PixTks_POCAx"        , &L1PixTks_POCAx        , &b_L1PixTks_POCAx);
      fChain->SetBranchAddress("L1PixTks_POCAy"        , &L1PixTks_POCAy        , &b_L1PixTks_POCAy);
      fChain->SetBranchAddress("L1PixTks_POCAz"        , &L1PixTks_POCAz        , &b_L1PixTks_POCAz);
      fChain->SetBranchAddress("L1PixTks_ChiSquared"   , &L1PixTks_ChiSquared   , &b_L1PixTks_ChiSquared);
      fChain->SetBranchAddress("L1PixTks_RInv"         , &L1PixTks_RInv         , &b_L1PixTks_RInv);
      fChain->SetBranchAddress("L1PixTks_SigmaRInv"    , &L1PixTks_SigmaRInv    , &b_L1PixTks_SigmaRInv);
      fChain->SetBranchAddress("L1PixTks_SigmaPhi0"    , &L1PixTks_SigmaPhi0    , &b_L1PixTks_SigmaPhi0);
      fChain->SetBranchAddress("L1PixTks_SigmaD0"      , &L1PixTks_SigmaD0      , &b_L1PixTks_SigmaD0);
      fChain->SetBranchAddress("L1PixTks_SigmaT"       , &L1PixTks_SigmaT       , &b_L1PixTks_SigmaT);
      fChain->SetBranchAddress("L1PixTks_SigmaZ0"      , &L1PixTks_SigmaZ0      , &b_L1PixTks_SigmaZ0);
      fChain->SetBranchAddress("L1PixTks_TTTrackIndex" , &L1PixTks_TTTrackIndex , &b_L1PixTks_TTTrackIndex);
      fChain->SetBranchAddress("L1PixTks_NPixHits"     , &L1PixTks_NPixHits     , &b_L1PixTks_NPixHits);
      fChain->SetBranchAddress("L1PixTks_PixHits_X"    , &L1PixTks_PixHits_X    , &b_L1PixTks_PixHits_X);
      fChain->SetBranchAddress("L1PixTks_PixHits_Y"    , &L1PixTks_PixHits_Y    , &b_L1PixTks_PixHits_Y);
      fChain->SetBranchAddress("L1PixTks_PixHits_Z"    , &L1PixTks_PixHits_Z    , &b_L1PixTks_PixHits_Z);
      fChain->SetBranchAddress("L1PixTks_PixHits_R"    , &L1PixTks_PixHits_R    , &b_L1PixTks_PixHits_R);
      fChain->SetBranchAddress("L1PixTks_PixHits_Phi"  , &L1PixTks_PixHits_Phi  , &b_L1PixTks_PixHits_Phi);
      fChain->SetBranchAddress("L1PixTks_PixHits_Type" , &L1PixTks_PixHits_Type , &b_L1PixTks_PixHits_Type);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_X"    , &L1PixTks_CandPixHits_X    , &b_L1PixTks_CandPixHits_X);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Y"    , &L1PixTks_CandPixHits_Y    , &b_L1PixTks_CandPixHits_Y);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Z"    , &L1PixTks_CandPixHits_Z    , &b_L1PixTks_CandPixHits_Z);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_R"    , &L1PixTks_CandPixHits_R    , &b_L1PixTks_CandPixHits_R);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Phi"  , &L1PixTks_CandPixHits_Phi  , &b_L1PixTks_CandPixHits_Phi);
      fChain->SetBranchAddress("L1PixTks_CandPixHits_Type" , &L1PixTks_CandPixHits_Type , &b_L1PixTks_CandPixHits_Type);
    }

  // L1CaloTaus
  if (0)
    {
      cout << "\tSetting L1 Calo Tau addresses." << endl;
      fChain->SetBranchAddress("L1CaloTau_E"   , &L1CaloTau_E   , &b_L1CaloTau_E);
      fChain->SetBranchAddress("L1CaloTau_Et"  , &L1CaloTau_Et  , &b_L1CaloTau_Et);
      fChain->SetBranchAddress("L1CaloTau_Eta" , &L1CaloTau_Eta , &b_L1CaloTau_Eta);
      fChain->SetBranchAddress("L1CaloTau_Phi" , &L1CaloTau_Phi , &b_L1CaloTau_Phi);
      fChain->SetBranchAddress("L1CaloTau_Bx"  , &L1CaloTau_Bx  , &b_L1CaloTau_Bx);
      fChain->SetBranchAddress("L1CaloTau_Type", &L1CaloTau_Type, &b_L1CaloTau_Type);
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

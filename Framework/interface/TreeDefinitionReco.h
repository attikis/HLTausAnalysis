#ifndef TreeDefinitionReco_h
#define TreeDefinitionReco_h

// User
#include "TreeDefinitionBase.h"

class TreeDefinitionReco : public virtual TreeDefinitionBase
{
 public:
  
  //Event tree
  std::vector<double>  *L1Tks_Pt;
  std::vector<double>  *L1Tks_Px;
  std::vector<double>  *L1Tks_Py;
  std::vector<double>  *L1Tks_Pz;
  std::vector<double>  *L1Tks_Eta;
  std::vector<double>  *L1Tks_Phi;
  std::vector<int>     *L1Tks_Charge;
  std::vector<double>  *L1Tks_POCAx;
  std::vector<double>  *L1Tks_POCAy;
  std::vector<double>  *L1Tks_POCAz;
  std::vector<double>  *L1Tks_ChiSquared;
  std::vector<double>  *L1Tks_StubPtConsistency;
  std::vector<double>  *L1Tks_Sector;
  std::vector<double>  *L1Tks_Wedge;
  std::vector<double>  *L1Tks_RInv;
  std::vector<std::vector<unsigned int> > *L1Tks_Stubs_isPS;
  std::vector<std::vector<unsigned int> > *L1Tks_Stubs_iDisk;
  std::vector<std::vector<unsigned int> > *L1Tks_Stubs_iLayer;
  std::vector<std::vector<unsigned int> > *L1Tks_Stubs_iPhi;
  std::vector<std::vector<unsigned int> > *L1Tks_Stubs_iRing;
  std::vector<std::vector<unsigned int> > *L1Tks_Stubs_iSide;
  std::vector<std::vector<unsigned int> > *L1Tks_Stubs_iZ;
  std::vector<bool> *L1Tks_IsGenuine;
  std::vector<bool> *L1Tks_IsUnknown;
  std::vector<bool> *L1Tks_IsCombinatoric;
  std::vector<double> *L1PixTks_Pt;
  std::vector<double> *L1PixTks_Px;
  std::vector<double> *L1PixTks_Py;
  std::vector<double> *L1PixTks_Pz;
  std::vector<double> *L1PixTks_Eta;
  std::vector<double> *L1PixTks_Phi;
  std::vector<int>    *L1PixTks_Charge;
  std::vector<double> *L1PixTks_POCAx;
  std::vector<double> *L1PixTks_POCAy;
  std::vector<double> *L1PixTks_POCAz;
  std::vector<double> *L1PixTks_ChiSquared;
  // std::vector<double> *L1PixTks_RedChiSquared;
  std::vector<double> *L1PixTks_RInv;
  std::vector<double> *L1PixTks_SigmaRInv;
  std::vector<double> *L1PixTks_SigmaPhi0;
  std::vector<double> *L1PixTks_SigmaD0;
  std::vector<double> *L1PixTks_SigmaT;
  std::vector<double> *L1PixTks_SigmaZ0;
  std::vector<int>    *L1PixTks_TTTrackIndex;
  std::vector<int>    *L1PixTks_NPixHits;
  std::vector<std::vector<double> > *L1PixTks_PixHits_X;
  std::vector<std::vector<double> > *L1PixTks_PixHits_Y;
  std::vector<std::vector<double> > *L1PixTks_PixHits_Z;
  std::vector<std::vector<double> > *L1PixTks_PixHits_R;
  std::vector<std::vector<double> > *L1PixTks_PixHits_Phi;
  std::vector<std::vector<int> >    *L1PixTks_PixHits_Type;
  std::vector<std::vector<double> > *L1PixTks_CandPixHits_X;
  std::vector<std::vector<double> > *L1PixTks_CandPixHits_Y;
  std::vector<std::vector<double> > *L1PixTks_CandPixHits_Z;
  std::vector<std::vector<double> > *L1PixTks_CandPixHits_R;
  std::vector<std::vector<double> > *L1PixTks_CandPixHits_Phi;
  std::vector<std::vector<int> >    *L1PixTks_CandPixHits_Type;

  std::vector<double>  *L1CaloTau_E;
  std::vector<double>  *L1CaloTau_Et;
  std::vector<double>  *L1CaloTau_Eta;
  std::vector<double>  *L1CaloTau_Phi;
  std::vector<int>     *L1CaloTau_Bx;
  std::vector<int>     *L1CaloTau_Type;

  std::vector< double > *L1TkJet_Pt;
  std::vector< double > *L1TkJet_Eta;
  std::vector< double > *L1TkJet_Phi;
  std::vector< double > *L1TkJet_E;
  std::vector< std::vector< int > > *L1TkJet_TTTrackIndex; 
  std::vector< double > *L1TkJet_Vertex;

   

  // List of branches
  TBranch        *b_L1Tks_Pt;
  TBranch        *b_L1Tks_Px;
  TBranch        *b_L1Tks_Py;
  TBranch        *b_L1Tks_Pz;
  TBranch        *b_L1Tks_Eta;
  TBranch        *b_L1Tks_Phi;
  TBranch        *b_L1Tks_Charge;
  TBranch        *b_L1Tks_POCAx;
  TBranch        *b_L1Tks_POCAy;
  TBranch        *b_L1Tks_POCAz;
  TBranch        *b_L1Tks_ChiSquared;
  TBranch        *b_L1Tks_StubPtConsistency;
  TBranch        *b_L1Tks_Sector;
  TBranch        *b_L1Tks_Wedge;
  TBranch        *b_L1Tks_RInv;
  TBranch        *b_L1Tks_Stubs_isPS;
  TBranch        *b_L1Tks_Stubs_iDisk;
  TBranch        *b_L1Tks_Stubs_iLayer;
  TBranch        *b_L1Tks_Stubs_iPhi;
  TBranch        *b_L1Tks_Stubs_iRing;
  TBranch        *b_L1Tks_Stubs_iSide;
  TBranch        *b_L1Tks_Stubs_iZ;
  TBranch        *b_L1Tks_IsGenuine;
  TBranch        *b_L1Tks_IsUnknown;
  TBranch        *b_L1Tks_IsCombinatoric;
  TBranch        *b_L1PixTks_Pt;
  TBranch        *b_L1PixTks_Px;
  TBranch        *b_L1PixTks_Py;
  TBranch        *b_L1PixTks_Pz;
  TBranch        *b_L1PixTks_Eta;
  TBranch        *b_L1PixTks_Phi;
  TBranch        *b_L1PixTks_Charge;
  TBranch        *b_L1PixTks_NPixHits;
  TBranch        *b_L1PixTks_PixHits_X;
  TBranch        *b_L1PixTks_PixHits_Y;
  TBranch        *b_L1PixTks_PixHits_Z;
  TBranch        *b_L1PixTks_PixHits_R;
  TBranch        *b_L1PixTks_PixHits_Phi;
  TBranch        *b_L1PixTks_PixHits_Type;
  TBranch        *b_L1PixTks_CandPixHits_X;
  TBranch        *b_L1PixTks_CandPixHits_Y;
  TBranch        *b_L1PixTks_CandPixHits_Z;
  TBranch        *b_L1PixTks_CandPixHits_R;
  TBranch        *b_L1PixTks_CandPixHits_Phi;
  TBranch        *b_L1PixTks_CandPixHits_Type;

  TBranch        *b_L1PixTks_POCAx;
  TBranch        *b_L1PixTks_POCAy;
  TBranch        *b_L1PixTks_POCAz;
  TBranch        *b_L1PixTks_ChiSquared;
  // TBranch        *b_L1PixTks_RedChiSquared;//new
  TBranch        *b_L1PixTks_RInv;
  TBranch        *b_L1PixTks_SigmaRInv;
  TBranch        *b_L1PixTks_SigmaPhi0;
  TBranch        *b_L1PixTks_SigmaD0;
  TBranch        *b_L1PixTks_SigmaT;
  TBranch        *b_L1PixTks_SigmaZ0;
  TBranch        *b_L1PixTks_TTTrackIndex;
  TBranch        *b_L1CaloTau_E;
  TBranch        *b_L1CaloTau_Et;
  TBranch        *b_L1CaloTau_Eta;
  TBranch        *b_L1CaloTau_Phi;
  TBranch        *b_L1CaloTau_Bx;
  TBranch        *b_L1CaloTau_Type;
  TBranch        *b_L1TkJet_Pt;
  TBranch        *b_L1TkJet_Eta;
  TBranch        *b_L1TkJet_Phi;
  TBranch        *b_L1TkJet_E;
  TBranch        *b_L1TkJet_TTTrackIndex; 
  TBranch        *b_L1TkJet_Vertex;
   

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


  // Set object pointer
  L1Tks_Pt                = 0;
  L1Tks_Px                = 0;
  L1Tks_Py                = 0;
  L1Tks_Pz                = 0;
  L1Tks_Eta               = 0;
  L1Tks_Phi               = 0;
  L1Tks_Charge            = 0;
  L1Tks_POCAx             = 0;
  L1Tks_POCAy             = 0;
  L1Tks_POCAz             = 0;
  L1Tks_ChiSquared        = 0;
  L1Tks_StubPtConsistency = 0;
  L1Tks_Sector            = 0;
  L1Tks_Wedge             = 0;
  L1Tks_RInv              = 0;
  L1Tks_Stubs_isPS        = 0;
  L1Tks_Stubs_iDisk       = 0;
  L1Tks_Stubs_iLayer      = 0;
  L1Tks_Stubs_iPhi        = 0;
  L1Tks_Stubs_iRing       = 0;
  L1Tks_Stubs_iSide       = 0;
  L1Tks_Stubs_iZ          = 0;
  L1Tks_IsGenuine         = 0;
  L1Tks_IsUnknown         = 0;
  L1Tks_IsCombinatoric    = 0;
  L1PixTks_Pt             = 0;
  L1PixTks_Px             = 0;
  L1PixTks_Py             = 0;
  L1PixTks_Pz             = 0;
  L1PixTks_Eta            = 0;
  L1PixTks_Phi            = 0;
  L1PixTks_Charge         = 0;
  L1PixTks_POCAx          = 0;
  L1PixTks_POCAy          = 0;
  L1PixTks_POCAz          = 0;
  L1PixTks_ChiSquared     = 0;
  // L1PixTks_RedChiSquared  = 0;//new
  L1PixTks_RInv           = 0;
  L1PixTks_SigmaRInv      = 0;
  L1PixTks_SigmaPhi0      = 0;
  L1PixTks_SigmaD0        = 0;
  L1PixTks_SigmaT         = 0;
  L1PixTks_SigmaZ0        = 0;
  L1PixTks_TTTrackIndex   = 0;
  L1PixTks_NPixHits       = 0;
  L1PixTks_PixHits_X      = 0;
  L1PixTks_PixHits_Y      = 0;
  L1PixTks_PixHits_Z      = 0;
  L1PixTks_PixHits_R      = 0;
  L1PixTks_PixHits_Phi    = 0;
  L1PixTks_PixHits_Type   = 0;
  L1PixTks_CandPixHits_X      = 0;
  L1PixTks_CandPixHits_Y      = 0;
  L1PixTks_CandPixHits_Z      = 0;
  L1PixTks_CandPixHits_R      = 0;
  L1PixTks_CandPixHits_Phi    = 0;
  L1PixTks_CandPixHits_Type   = 0;
  L1CaloTau_E             = 0;
  L1CaloTau_Et            = 0;
  L1CaloTau_Eta           = 0;
  L1CaloTau_Phi           = 0;
  L1CaloTau_Bx            = 0;
  L1CaloTau_Type          = 0;
  L1TkJet_Pt              = 0;
  L1TkJet_Eta             = 0;
  L1TkJet_Phi             = 0;
  L1TkJet_E               = 0;
  L1TkJet_TTTrackIndex    = 0;
  L1TkJet_Vertex          = 0;


  std::cout << "I N F O ! TreeDefinitionReco::InitReco(...) - Setting branch addresses and branch pointers." << std::endl;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1); 
					
  std::cout << "I N F O ! TreeDefinitionReco::InitReco(...) - Setting L1 Tracks addresses." << std::endl;
  fChain->SetBranchAddress("L1Tks_Pt"               , &L1Tks_Pt               , &b_L1Tks_Pt);
  fChain->SetBranchAddress("L1Tks_Px"               , &L1Tks_Px               , &b_L1Tks_Px);
  fChain->SetBranchAddress("L1Tks_Py"               , &L1Tks_Py               , &b_L1Tks_Py);
  fChain->SetBranchAddress("L1Tks_Pz"               , &L1Tks_Pz               , &b_L1Tks_Pz);
  fChain->SetBranchAddress("L1Tks_Eta"              , &L1Tks_Eta              , &b_L1Tks_Eta);
  fChain->SetBranchAddress("L1Tks_Phi"              , &L1Tks_Phi              , &b_L1Tks_Phi);
  fChain->SetBranchAddress("L1Tks_Charge"           , &L1Tks_Charge           , &b_L1Tks_Charge);
  fChain->SetBranchAddress("L1Tks_POCAx"            , &L1Tks_POCAx            , &b_L1Tks_POCAx);
  fChain->SetBranchAddress("L1Tks_POCAy"            , &L1Tks_POCAy            , &b_L1Tks_POCAy);
  fChain->SetBranchAddress("L1Tks_POCAz"            , &L1Tks_POCAz            , &b_L1Tks_POCAz);
  fChain->SetBranchAddress("L1Tks_ChiSquared"       , &L1Tks_ChiSquared       , &b_L1Tks_ChiSquared);
  fChain->SetBranchAddress("L1Tks_StubPtConsistency", &L1Tks_StubPtConsistency, &b_L1Tks_StubPtConsistency);
  fChain->SetBranchAddress("L1Tks_Sector"           , &L1Tks_Sector           , &b_L1Tks_Sector);
  fChain->SetBranchAddress("L1Tks_Wedge"            , &L1Tks_Wedge            , &b_L1Tks_Wedge);
  fChain->SetBranchAddress("L1Tks_RInv"             , &L1Tks_RInv             , &b_L1Tks_RInv);
  fChain->SetBranchAddress("L1Tks_Stubs_isPS"       , &L1Tks_Stubs_isPS       , &b_L1Tks_Stubs_isPS);
  fChain->SetBranchAddress("L1Tks_Stubs_iDisk"      , &L1Tks_Stubs_iDisk      , &b_L1Tks_Stubs_iDisk);
  fChain->SetBranchAddress("L1Tks_Stubs_iLayer"     , &L1Tks_Stubs_iLayer     , &b_L1Tks_Stubs_iLayer);
  fChain->SetBranchAddress("L1Tks_Stubs_iPhi"       , &L1Tks_Stubs_iPhi       , &b_L1Tks_Stubs_iPhi);
  fChain->SetBranchAddress("L1Tks_Stubs_iRing"      , &L1Tks_Stubs_iRing      , &b_L1Tks_Stubs_iRing);
  fChain->SetBranchAddress("L1Tks_Stubs_iSide"      , &L1Tks_Stubs_iSide      , &b_L1Tks_Stubs_iSide);
  fChain->SetBranchAddress("L1Tks_Stubs_iZ"         , &L1Tks_Stubs_iZ         , &b_L1Tks_Stubs_iZ);
  fChain->SetBranchAddress("L1Tks_IsGenuine"        , &L1Tks_IsGenuine        , &b_L1Tks_IsGenuine);
  fChain->SetBranchAddress("L1Tks_IsUnknown"        , &L1Tks_IsUnknown        , &b_L1Tks_IsUnknown);
  fChain->SetBranchAddress("L1Tks_IsCombinatoric"   , &L1Tks_IsCombinatoric   , &b_L1Tks_IsCombinatoric);

  std::cout << "I N F O ! TreeDefinitionReco::InitReco(...) - Setting L1 Pixel Tracks addresses." << std::endl;
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
  // fChain->SetBranchAddress("L1PixTks_RedChiSquared", &L1PixTks_RedChiSquared, &b_L1PixTks_RedChiSquared);//new
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

  std::cout << "I N F O ! TreeDefinitionReco::InitReco(...) - Setting L1 Calo Tau addresses." << std::endl;
  fChain->SetBranchAddress("L1CaloTau_E"   , &L1CaloTau_E   , &b_L1CaloTau_E);
  fChain->SetBranchAddress("L1CaloTau_Et"  , &L1CaloTau_Et  , &b_L1CaloTau_Et);
  fChain->SetBranchAddress("L1CaloTau_Eta" , &L1CaloTau_Eta , &b_L1CaloTau_Eta);
  fChain->SetBranchAddress("L1CaloTau_Phi" , &L1CaloTau_Phi , &b_L1CaloTau_Phi);
  fChain->SetBranchAddress("L1CaloTau_Bx"  , &L1CaloTau_Bx  , &b_L1CaloTau_Bx);
  fChain->SetBranchAddress("L1CaloTau_Type", &L1CaloTau_Type, &b_L1CaloTau_Type);

  std::cout << "I N F O ! TreeDefinitionReco::InitReco(...) - Setting L1 Jet addresses." << std::endl;
  fChain->SetBranchAddress("L1TkJet_Pt"          , &L1TkJet_Pt          , &b_L1TkJet_Pt);
  fChain->SetBranchAddress("L1TkJet_Eta"         , &L1TkJet_Eta         , &b_L1TkJet_Eta);
  fChain->SetBranchAddress("L1TkJet_Phi"         , &L1TkJet_Phi         , &b_L1TkJet_Phi);
  fChain->SetBranchAddress("L1TkJet_E"           , &L1TkJet_E           , &b_L1TkJet_E);
  fChain->SetBranchAddress("L1TkJet_TTTrackIndex", &L1TkJet_TTTrackIndex, &b_L1TkJet_TTTrackIndex);
  fChain->SetBranchAddress("L1TkJet_Vertex"      , &L1TkJet_Vertex      , &b_L1TkJet_Vertex);

  Notify();
}

#endif // TreeDefinitionReco_h

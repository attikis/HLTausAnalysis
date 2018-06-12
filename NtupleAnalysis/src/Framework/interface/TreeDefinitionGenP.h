#ifndef TreeDefinitionGenP_h
#define TreeDefinitionGenP_h

// User
#include "TreeDefinitionBase.h"
#include <TFriendElement.h>
#include <TList.h>


class TreeDefinitionGenP : public virtual TreeDefinitionBase
{
 public:

  Float_t    weight;
  Float_t    pthat;
  UInt_t   nVtx;
  UInt_t   nMeanPU;

  // GenParticles
  UInt_t              GenP_N;
  std::vector<int>    GenP_PdgId;
  std::vector<int>    GenP_Status;
  std::vector<int>    GenP_Parent;
  std::vector<float>  GenP_Pt;
  std::vector<float>  GenP_Eta;
  std::vector<float>  GenP_Phi;
  std::vector<float>  GenP_E;
  std::vector<int>    GenP_Charge;

  // GenJets 
  Int_t           GenJet_N;
  vector<float>   GenJet_Pt;
  vector<float>   GenJet_Eta;
  vector<float>   GenJet_Phi;
  vector<float>   GenJet_M;

  // **************List of branches***************

  TBranch *b_weight;
  TBranch *b_pthat;
  TBranch *b_nVtx;
  TBranch *b_nMeanPU;

  // GenParticles
  TBranch *b_GenP_N;
  TBranch *b_GenP_PdgId;
  TBranch *b_GenP_Status;
  TBranch *b_GenP_Parent;
  TBranch *b_GenP_Pt;
  TBranch *b_GenP_Eta;
  TBranch *b_GenP_Phi;
  TBranch *b_GenP_E;
  TBranch *b_GenP_Charge;

  // GenJets
  TBranch *b_GenJet_N; 
  TBranch *b_GenJet_Pt;
  TBranch *b_GenJet_Eta;
  TBranch *b_GenJet_Phi;
  TBranch *b_GenJet_M;  

  virtual void InitGenP(TChain *chain);

};

void TreeDefinitionGenP::InitGenP(TChain *chain)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  std::cout << "=== TreeDefinitionGenP::InitGenP()" << std::endl;

  weight  = 0;
  pthat   = 0;
  nVtx    = 0;
  nMeanPU = 0;

  // GenParticles
  GenP_N    = 0;
  
  // GenJets
  GenJet_N  = 0;

  // Set branch addresses and branch pointers
  if (!chain) return;
  fChain = chain;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  // L1 GEN info
  if (doGenerator)
    {
      fGenerator->SetBranchAddress("weight"       , &weight      , &b_weight);
      fGenerator->SetBranchAddress("pthat"        , &pthat       , &b_pthat);
      fGenerator->SetBranchAddress("nVtx"         , &nVtx        , &b_nVtx);
      fGenerator->SetBranchAddress("nMeanPU"      , &nMeanPU     , &b_nMeanPU);

      // GenParticles
      fGenerator->SetBranchAddress("nPart"        , &GenP_N      , &b_GenP_N);
      fGenerator->SetBranchAddress("partId"       , &GenP_PdgId  , &b_GenP_PdgId);
      fGenerator->SetBranchAddress("partStat"     , &GenP_Status , &b_GenP_Status);
      fGenerator->SetBranchAddress("partParent"   , &GenP_Parent , &b_GenP_Parent);
      fGenerator->SetBranchAddress("partPt"       , &GenP_Pt     , &b_GenP_Pt);
      fGenerator->SetBranchAddress("partEta"      , &GenP_Eta    , &b_GenP_Eta);
      fGenerator->SetBranchAddress("partPhi"      , &GenP_Phi    , &b_GenP_Phi);
      fGenerator->SetBranchAddress("partE"        , &GenP_E      , &b_GenP_E);
      fGenerator->SetBranchAddress("partCh"       , &GenP_Charge , &b_GenP_Charge);

      // GenJets
      fGenerator->SetBranchAddress("nJet"         , &GenJet_N    , &b_GenJet_N);
      fGenerator->SetBranchAddress("jetPt"        , &GenJet_Pt   , &b_GenJet_Pt);
      fGenerator->SetBranchAddress("jetEta"       , &GenJet_Eta  , &b_GenJet_Eta);
      fGenerator->SetBranchAddress("jetPhi"       , &GenJet_Phi  , &b_GenJet_Phi);
      fGenerator->SetBranchAddress("jetM"         , &GenJet_M    , &b_GenJet_M);
      
      // Add friend (GeneratorTree)
      fChain -> AddFriend(fGenerator);
    }

  // Set Make Class for Generator tree
  TFriendElement *friendTreeElement = (TFriendElement*)fChain->GetListOfFriends()->Last();
  friendTreeElement->GetTree()->SetMakeClass(1);

  return;

}

#endif  // TreeDefinitionGenP_h

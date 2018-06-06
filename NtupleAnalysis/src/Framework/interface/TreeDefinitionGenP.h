#ifndef TreeDefinitionGenP_h
#define TreeDefinitionGenP_h

// User
#include "TreeDefinitionBase.h"

class TreeDefinitionGenP : public virtual TreeDefinitionBase
{
 public:

  //marina
  float          weight;
  float          pthat;
  unsigned int   nVtx;
  unsigned int   nMeanPU;

  // GenParticles
  unsigned int        GenP_N;
  std::vector<int>    GenP_PdgId;
  std::vector<int>    GenP_Status;
  std::vector<int>    GenP_Parent;
  std::vector<float>  GenP_Pt;
  std::vector<float>  GenP_Eta;
  std::vector<float>  GenP_Phi;
  std::vector<float>  GenP_E;
  std::vector<int>    GenP_Charge;
  

  //marina
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


  virtual void InitGenP(TChain *chain);

};

//void TreeDefinitionGenP::InitGenP(TTree *tree)
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
  
  /*
  // GenParticles
  GenP_N          = 0;
  GenP_Pt         = 0;
  GenP_Eta        = 0;
  GenP_Phi        = 0;
  GenP_E          = 0;
  //GenP_Mass       = 0;
  GenP_Charge     = 0;
  GenP_PdgId      = 0;
  GenP_Status     = 0;
  */


  // Set branch addresses and branch pointers
  if (!chain) return;
  fChain = chain;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  // L1 GEN info
  if (1)
    {
      fChain->SetBranchAddress("weight"       , &weight      , &b_weight);
      fChain->SetBranchAddress("pthat"        , &pthat       , &b_pthat);
      fChain->SetBranchAddress("nVtx"         , &nVtx        , &b_nVtx);
      fChain->SetBranchAddress("nMeanPU"      , &nMeanPU     , &b_nMeanPU);
      fChain->SetBranchAddress("nPart"        , &GenP_N      , &b_GenP_N);
      fChain->SetBranchAddress("partId"       , &GenP_PdgId  , &b_GenP_PdgId);
      fChain->SetBranchAddress("partStat"     , &GenP_Status , &b_GenP_Status);
      fChain->SetBranchAddress("partPt"       , &GenP_Pt     , &b_GenP_Pt);
      fChain->SetBranchAddress("partEta"      , &GenP_Eta    , &b_GenP_Eta);
      fChain->SetBranchAddress("partPhi"      , &GenP_Phi    , &b_GenP_Phi);
      fChain->SetBranchAddress("partE"        , &GenP_E      , &b_GenP_E);
      fChain->SetBranchAddress("partCh"       , &GenP_Charge , &b_GenP_Charge);


    }

  return;

}

#endif  // TreeDefinitionGenP_h

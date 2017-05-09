#ifndef TreeDefinitionGenP_h
#define TreeDefinitionGenP_h

// User
#include "TreeDefinitionBase.h"

class TreeDefinitionGenP : public virtual TreeDefinitionBase
{
 public:

  Int_t           RunNumber;
  Int_t           EvtNumber;
  Double_t        HepMCEvt_VtxX;
  Double_t        HepMCEvt_VtxY;
  Double_t        HepMCEvt_VtxZ;

  // GenParticles
  std::vector<double>  *GenP_Pt;
  std::vector<double>  *GenP_Eta;
  std::vector<double>  *GenP_Phi;
  std::vector<double>  *GenP_Mass;
  std::vector<int>     *GenP_Charge;
  std::vector<int>     *GenP_PdgId;
  std::vector<int>     *GenP_Status;
  std::vector<double>  *GenP_VertexX;
  std::vector<double>  *GenP_VertexY;
  std::vector<double>  *GenP_VertexZ;
  std::vector<std::vector<unsigned short> > *GenP_Mothers;
  std::vector<std::vector<unsigned short> > *GenP_Daughters;

  // Tracking Particles
  std::vector<double>  *TP_Pt;
  std::vector<double>  *TP_Eta;
  std::vector<double>  *TP_Phi;
  std::vector<int>     *TP_Charge;
  std::vector<int>     *TP_PdgId;
  std::vector<double>  *TP_d0_propagated;
  std::vector<double>  *TP_z0_propagated;
  std::vector<int>     *TP_NMatch;
  std::vector<int>     *TP_TTTrackIndex;
  std::vector<int>     *TP_NClusters; 
  std::vector<int>     *TP_NStubs;
  std::vector<int>     *TP_NTracks; 
  std::vector<double>  *TP_d0_produced;
  std::vector<double>  *TP_x0_produced;
  std::vector<double>  *TP_y0_produced;
  std::vector<double>  *TP_z0_produced;
  std::vector<int>     *TP_EventId;
  
  // List of branches
  TBranch *b_RunNumber;
  TBranch *b_EvtNumber;
  TBranch *b_HepMCEvt_VtxX;
  TBranch *b_HepMCEvt_VtxY;
  TBranch *b_HepMCEvt_VtxZ;

  // GenParticles
  TBranch *b_GenP_Pt;
  TBranch *b_GenP_Eta;
  TBranch *b_GenP_Phi;
  TBranch *b_GenP_Mass;
  TBranch *b_GenP_Charge;
  TBranch *b_GenP_PdgId;
  TBranch *b_GenP_Status;
  TBranch *b_GenP_VertexX;
  TBranch *b_GenP_VertexY;
  TBranch *b_GenP_VertexZ;
  TBranch *b_GenP_Mothers;
  TBranch *b_GenP_Daughters;

  // Tracking Particles
  TBranch *b_TP_Pt;
  TBranch *b_TP_Eta;
  TBranch *b_TP_Phi;
  TBranch *b_TP_Charge;
  TBranch *b_TP_PdgId;
  TBranch *b_TP_d0_propagated;
  TBranch *b_TP_z0_propagated;
  TBranch *b_TP_NMatch;
  TBranch *b_TP_TTTrackIndex;
  TBranch *b_TP_NClusters;
  TBranch *b_TP_NStubs;
  TBranch *b_TP_NTracks;
  TBranch *b_TP_d0_produced;
  TBranch *b_TP_x0_produced;
  TBranch *b_TP_y0_produced;
  TBranch *b_TP_z0_produced;
  TBranch *b_TP_EventId;

  virtual void InitGenP(TTree *tree);

};

void TreeDefinitionGenP::InitGenP(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  std::cout << "=== TreeDefinitionGenP::InitGenP()" << std::endl;
   
   
  // GenParticles
  GenP_Pt         = 0;
  GenP_Eta        = 0;
  GenP_Phi        = 0;
  GenP_Mass       = 0;
  GenP_Charge     = 0;
  GenP_PdgId      = 0;
  GenP_Status     = 0;
  GenP_VertexX    = 0;
  GenP_VertexY    = 0;
  GenP_VertexZ    = 0;
  GenP_Mothers    = 0;
  GenP_Daughters  = 0;

  // Tracking Particles
  TP_Pt           = 0;
  TP_Eta          = 0;
  TP_Phi          = 0;
  TP_Charge       = 0;
  TP_PdgId        = 0;
  TP_d0_propagated= 0;
  TP_z0_propagated= 0;
  TP_NMatch       = 0;
  TP_TTTrackIndex = 0;
  TP_NClusters    = 0;
  TP_NStubs       = 0;
  TP_NTracks      = 0;
  TP_d0_produced  = 0;
  TP_x0_produced  = 0;
  TP_y0_produced  = 0;
  TP_z0_produced  = 0;
  TP_x0_produced  = 0;
  TP_EventId      = 0;

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  if (1)
    {
      fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
      fChain->SetBranchAddress("EvtNumber", &EvtNumber, &b_EvtNumber);
    }
   
  if (1)
    {
      std::cout << "\tSetting HepMC addresses." << std::endl;
      fChain->SetBranchAddress("HepMCEvt_VtxX", &HepMCEvt_VtxX, &b_HepMCEvt_VtxX);
      fChain->SetBranchAddress("HepMCEvt_VtxY", &HepMCEvt_VtxY, &b_HepMCEvt_VtxY);
      fChain->SetBranchAddress("HepMCEvt_VtxZ", &HepMCEvt_VtxZ, &b_HepMCEvt_VtxZ);
    }

  // GenParticles
  if (1)
    {
      std::cout << "\tSetting GenParticles addresses." << std::endl;
      fChain->SetBranchAddress("GenP_Pt"       , &GenP_Pt       , &b_GenP_Pt);
      fChain->SetBranchAddress("GenP_Eta"      , &GenP_Eta      , &b_GenP_Eta);
      fChain->SetBranchAddress("GenP_Phi"      , &GenP_Phi      , &b_GenP_Phi);
      fChain->SetBranchAddress("GenP_Mass"     , &GenP_Mass     , &b_GenP_Mass);
      fChain->SetBranchAddress("GenP_Charge"   , &GenP_Charge   , &b_GenP_Charge);
      fChain->SetBranchAddress("GenP_PdgId"    , &GenP_PdgId    , &b_GenP_PdgId);
      fChain->SetBranchAddress("GenP_Status"   , &GenP_Status   , &b_GenP_Status);
      fChain->SetBranchAddress("GenP_VertexX"  , &GenP_VertexX  , &b_GenP_VertexX);
      fChain->SetBranchAddress("GenP_VertexY"  , &GenP_VertexY  , &b_GenP_VertexY);
      fChain->SetBranchAddress("GenP_VertexZ"  , &GenP_VertexZ  , &b_GenP_VertexZ);
      fChain->SetBranchAddress("GenP_Mothers"  , &GenP_Mothers  , &b_GenP_Mothers);
      fChain->SetBranchAddress("GenP_Daughters", &GenP_Daughters, &b_GenP_Daughters);
    }

  // Tracking Particles
  if (1)
    {
      std::cout << "\tSetting Tracking Particles addresses." << std::endl;
      fChain->SetBranchAddress("TP_Pt"           , &TP_Pt           , &b_TP_Pt);
      fChain->SetBranchAddress("TP_Eta"          , &TP_Eta          , &b_TP_Eta);
      fChain->SetBranchAddress("TP_Phi"          , &TP_Phi          , &b_TP_Phi);
      fChain->SetBranchAddress("TP_Charge"       , &TP_Charge       , &b_TP_Charge);
      fChain->SetBranchAddress("TP_PdgId"        , &TP_PdgId        , &b_TP_PdgId);
      fChain->SetBranchAddress("TP_d0_propagated", &TP_d0_propagated, &b_TP_d0_propagated);
      fChain->SetBranchAddress("TP_z0_propagated", &TP_z0_propagated, &b_TP_z0_propagated);
      fChain->SetBranchAddress("TP_NMatch"       , &TP_NMatch       , &b_TP_NMatch);
      fChain->SetBranchAddress("TP_TTTrackIndex" , &TP_TTTrackIndex , &b_TP_TTTrackIndex);
      fChain->SetBranchAddress("TP_NClusters"    , &TP_NClusters    , &b_TP_NClusters);
      fChain->SetBranchAddress("TP_NStubs"       , &TP_NStubs       , &b_TP_NStubs);
      fChain->SetBranchAddress("TP_NTracks"      , &TP_NTracks      , &b_TP_NTracks);
      fChain->SetBranchAddress("TP_x0_produced"  , &TP_x0_produced  , &b_TP_x0_produced);
      fChain->SetBranchAddress("TP_y0_produced"  , &TP_y0_produced  , &b_TP_y0_produced);
      fChain->SetBranchAddress("TP_z0_produced"  , &TP_z0_produced  , &b_TP_z0_produced);
      fChain->SetBranchAddress("TP_EventId"      , &TP_EventId      , &b_TP_EventId);
    }

}

#endif  // TreeDefinitionGenP_h

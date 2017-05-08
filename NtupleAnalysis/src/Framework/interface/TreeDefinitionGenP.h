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

   std::vector<double>  *TP_Pt;
   std::vector<double>  *TP_Px;
   std::vector<double>  *TP_Py;
   std::vector<double>  *TP_Pz;
   std::vector<double>  *TP_Eta;
   std::vector<double>  *TP_Phi;
   std::vector<int>     *TP_NMatch;
   std::vector<int>     *TP_TTTrackIndex;
   std::vector<int>     *TP_Charge;
   std::vector<int>     *TP_PdgId;
   std::vector<double>  *TP_POCAx;
   std::vector<double>  *TP_POCAy;                                                                                                                           
   std::vector<double>  *TP_POCAz;                                                                                                                           
   std::vector<int>     *TP_TTClusters; 
   std::vector<int>     *TP_TTStubs;
   std::vector<int>     *TP_TTTracks; 

  // List of branches
   TBranch *b_RunNumber;
   TBranch *b_EvtNumber;
   TBranch *b_HepMCEvt_VtxX;
   TBranch *b_HepMCEvt_VtxY;
   TBranch *b_HepMCEvt_VtxZ;

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

   TBranch *b_TP_Pt;
   TBranch *b_TP_Px;
   TBranch *b_TP_Py;
   TBranch *b_TP_Pz;
   TBranch *b_TP_Eta;
   TBranch *b_TP_Phi;
   TBranch *b_TP_NMatch;
   TBranch *b_TP_TTTrackIndex;
   TBranch *b_TP_Charge;
   TBranch *b_TP_PdgId;
   TBranch *b_TP_POCAx;
   TBranch *b_TP_POCAy;                                                                                                                           
   TBranch *b_TP_POCAz;                                                                                                                           
   TBranch *b_TP_TTClusters;
   TBranch *b_TP_TTStubs;
   TBranch *b_TP_TTTracks;
  
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
   
  // Set object pointer
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

   TP_Pt           = 0;
   TP_Px           = 0;
   TP_Py           = 0;
   TP_Pz           = 0;
   TP_Eta          = 0;
   TP_Phi          = 0;
   TP_NMatch       = 0;
   TP_TTTrackIndex = 0;
   TP_Charge       = 0;
   TP_PdgId        = 0;
   TP_POCAx        = 0;
   TP_POCAy        = 0;                                                                                                                           
   TP_POCAz        = 0;                                                                                                                           
   TP_TTClusters   = 0;
   TP_TTStubs      = 0;
   TP_TTTracks     = 0;

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
   
   if (1)
     {
       std::cout << "\tSetting GenP addresses." << std::endl;
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

   if (1)
     {
       std::cout << "\tSetting TP addresses." << std::endl;
       fChain->SetBranchAddress("TP_Pt"           , &TP_Pt          , &b_TP_Pt);
       fChain->SetBranchAddress("TP_Px"           , &TP_Px          , &b_TP_Px);
       fChain->SetBranchAddress("TP_Py"           , &TP_Py          , &b_TP_Py);
       fChain->SetBranchAddress("TP_Pz"           , &TP_Pz          , &b_TP_Pz);
       fChain->SetBranchAddress("TP_Eta"          , &TP_Eta         , &b_TP_Eta);
       fChain->SetBranchAddress("TP_Phi"          , &TP_Phi         , &b_TP_Phi);
       fChain->SetBranchAddress("TP_NMatch"       , &TP_NMatch      , &b_TP_NMatch);
       fChain->SetBranchAddress("TP_TTTrackIndex" , &TP_TTTrackIndex, &b_TP_TTTrackIndex);
       fChain->SetBranchAddress("TP_Charge"       , &TP_Charge      , &b_TP_Charge);
       fChain->SetBranchAddress("TP_PdgId"        , &TP_PdgId       , &b_TP_PdgId);
       fChain->SetBranchAddress("TP_x0_produced"  , &TP_POCAx       , &b_TP_POCAx);
       fChain->SetBranchAddress("TP_y0_produced"  , &TP_POCAy       , &b_TP_POCAy);
       fChain->SetBranchAddress("TP_z0_produced"  , &TP_POCAz       , &b_TP_POCAz);
       fChain->SetBranchAddress("TP_TTClusters"   , &TP_TTClusters  , &b_TP_TTClusters);
       fChain->SetBranchAddress("TP_TTStubs"      , &TP_TTStubs     , &b_TP_TTStubs);
       fChain->SetBranchAddress("TP_TTTracks"     , &TP_TTTracks    , &b_TP_TTTracks);
     }

}

#endif  // TreeDefinitionGenP_h

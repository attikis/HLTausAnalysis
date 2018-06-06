#ifndef TreeDefinitionReco_h
#define TreeDefinitionReco_h

// User
#include "TreeDefinitionBase.h"
#include <TFriendElement.h>
#include <TList.h>

using namespace std;

class TreeDefinitionReco : public virtual TreeDefinitionBase
{
 public:
  
  // Event
  UInt_t run;
  //Long64_t event;
  unsigned int lumi;
  //UInt_t bx;

  // Calo Tower 
  short nECALTP;
  vector<short> ecalTP_iEta;
  vector<short> hcalTP_iEta;




  // Event
  TBranch* b_run;
  TBranch* b_lumi;

  // Calo Tower
  TBranch   *b_nECALTP;
  TBranch   *b_ecalTP_iEta;
  TBranch   *b_hcalTP_iEta;


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
  

  run     = 0;
  lumi    = 0;
  nECALTP = 0;


  cout << "\tSetting branch addresses and branch pointers." << endl;
  if (!chain) return;
  fChain = chain ;
  fCurrent = -1;

  
  fChain->SetMakeClass(1);

  //cout << "**********************SETMAKECLASS-1**********************" << endl;
  
  //TList *friendTreeElement2 = fChain->GetListOfFriends();                                                                                                 
  //friendTreeElement2->ls();
  

  TFriendElement *friendTreeElement = (TFriendElement*)fChain->GetListOfFriends()->First();
  while (friendTreeElement){
    friendTreeElement->GetTree()->SetMakeClass(1);
    friendTreeElement = (TFriendElement*)fChain->GetListOfFriends()->After(friendTreeElement);
  }                                                                                                                                                                      

  //cout << "**********************SETMAKECLASS-2**********************" << endl;
  

  if(1)                                                                                                                                                                 
    { 
      fChain->SetBranchAddress("run"                , &run                , &b_run);
      fChain->SetBranchAddress("lumi"               , &lumi               , &b_lumi);
      fChain->SetBranchAddress("nECALTP"            , &nECALTP            , &b_nECALTP);
      fChain->SetBranchAddress("ecalTPieta"         , &ecalTP_iEta        , &b_ecalTP_iEta);
      fChain->SetBranchAddress("hcalTPieta"         , &hcalTP_iEta        , &b_hcalTP_iEta);
    }

  
  Notify();
}

#endif // TreeDefinitionReco_h

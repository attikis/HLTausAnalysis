#ifndef TreeReaderReco_h
#define TreeReaderReco_h

// System
#include <vector>

// User
#include "TreeDefinitionReco.h"
#include "../src/FileOpener.C"

// ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

class TreeReaderReco : public virtual TREEDEFINITIONRECO, public FileOpener 
{
  public: 
    TreeReaderReco() { 
      std::cout << " TreeReaderReco: The contructor needs some arguments \n";
    };
    TreeReaderReco(const std::string SamplePath, const std::string SampleName, TTree *tree = 0);

    virtual ~TreeReaderReco();
    virtual void     Init(TTree *tree);
    virtual Int_t    Cut(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    
};

#endif // TreeReaderReco_h

#ifdef TreeReaderReco_cxx

TreeReaderReco::TreeReaderReco(const std::string SamplePath, const std::string SampleName, TTree *tree)
{
  std::cout << "=== TreeReaderReco::TreeReaderReco()" << std::endl;
  
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
   
    // TChain* chain = new TChain("TkTauFromCaloNTupleMaker/EvtTree");
    TChain* chain = new TChain("Raw2TTreeMaker/EvtTree");  // tree with tracks etc.
    TChain* friendChain = new TChain("l1UpgradeTree/L1UpgradeTree"); // tree with calo taus

    OpenFile(SamplePath, SampleName, chain);
    OpenFile(SamplePath, SampleName, friendChain);    

    // Set the Tree
    tree = chain;
    std::cout << "\tSet tree with name \"" << tree->GetName() << "\" as the principal tree." << std::endl;

    // Add Friend Tree. Then the Tree can access all data of its friend just like its own data
    tree->AddFriend(friendChain);
    std::cout << "\tAdded tree \"" << chain->GetListOfFriends()->Last()->GetName() << "\" as a friend to the principal tree." << std::endl;

  }

  Init(tree);

}

TreeReaderReco::~TreeReaderReco()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Long64_t TreeReaderReco::LoadTree(Long64_t entry)
{
  //Settng the environment to read an entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class())) return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void TreeReaderReco::Init(TTree *tree)
{
  InitReco(tree);
  Notify();
}

Int_t TreeReaderReco::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  std::cout << "entry " << entry << std::endl; // to avoid unused variable warning
  return -1;
}

#endif //TreeReaderReco_cxx

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
#include <TFriendElement.h>

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

    // Creat a new TChain pointer with the TTree name
    TChain* chain = new TChain("Raw2TTreeMaker/EvtTree");

    // Associated the ROOT files to the TChain
    std::cout << "\tGetting ROOT files for dataset " << SampleName << " from dir " << SamplePath << "/" << std::endl;
    OpenFiles(SamplePath, SampleName, chain);

    // Set the Tree
    tree = chain;
    std::cout << "=== TreeReaderReco::TreeReaderReco()" << std::endl;
    std::cout << "\tSet tree with name \"" << tree->GetName() << "\" as the principal tree." << std::endl;

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

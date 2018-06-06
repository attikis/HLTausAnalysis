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
    //TreeReaderReco(const std::string SamplePath, const std::string SampleName, TTree *tree = 0);
    TreeReaderReco(const std::string SamplePath, const std::string SampleName, TChain *chain = 0);

    virtual ~TreeReaderReco();
    //virtual void     Init(TTree *tree);
    virtual void     Init(TChain *chain);
    virtual Int_t    Cut(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    
};

#endif // TreeReaderReco_h

#ifdef TreeReaderReco_cxx

//TreeReaderReco::TreeReaderReco(const std::string SamplePath, const std::string SampleName, TTree *tree
TreeReaderReco::TreeReaderReco(const std::string SamplePath, const std::string SampleName, TChain *chain)
{
  std::cout << "=== TreeReaderReco::TreeReaderReco()" << std::endl;
  
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  //if (tree == 0) {
  if (chain == 0) {
  
    // Creat a new TChain pointer with the TTree name
    //TChain* chain = new TChain("Raw2TTreeMaker/EvtTree");
    TChain* fMainChain    = new TChain("l1EventTree/L1EventTree");
    TChain* fCaloTower        = new TChain("l1CaloTowerTree/L1CaloTowerTree");
    TChain* fUpgradeTfMuon    = new TChain("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree");
    TChain* fUpgrade          = new TChain("l1UpgradeTree/L1UpgradeTree");
    TChain* fuGT              = new TChain("l1uGTTree/L1uGTTree");
    TChain* fHO               = new TChain("l1HOTree/L1HOTree");
    TChain* fUpgradeTfMuonEmu = new TChain("l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree");
    TChain* fCaloTowerEmu     = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
    TChain* fUpgradeEmu       = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    TChain* fuGTEmu           = new TChain("l1uGTEmuTree/L1uGTTree");
    TChain* fGenerator        = new TChain("l1GeneratorTree/L1GenTree");
    
    // Associated the ROOT files to the TChain
    std::cout << "\tGetting ROOT files for dataset " << SampleName << " and adding them to the chain." << std::endl;
    OpenFiles(SamplePath, SampleName, fMainChain);

    std::cout<<"********OPENFILES-1*********"<<std::endl;

    // Fix me: check that trees exist
    OpenFiles(SamplePath, SampleName, fCaloTower        );
    OpenFiles(SamplePath, SampleName, fUpgradeTfMuon    );
    OpenFiles(SamplePath, SampleName, fUpgrade          );
    OpenFiles(SamplePath, SampleName, fuGT              );
    OpenFiles(SamplePath, SampleName, fHO               );
    OpenFiles(SamplePath, SampleName, fUpgradeTfMuonEmu );
    OpenFiles(SamplePath, SampleName, fCaloTowerEmu);
    OpenFiles(SamplePath, SampleName, fUpgradeEmu       );
    OpenFiles(SamplePath, SampleName, fuGTEmu           );
    OpenFiles(SamplePath, SampleName, fGenerator);

    std::cout<<"********OPENFILES-2*********"<<std::endl;

    
    // Set the Tree
    //tree = chain;
    chain = fMainChain;

    std::cout<<"********friends1*********"<<std::endl;

    // Set friendships
    //chain -> AddFriend(fCaloTower);
    chain -> AddFriend(fUpgradeTfMuon);
    chain -> AddFriend(fUpgrade);
    chain -> AddFriend(fuGT);
    chain -> AddFriend(fHO);
    chain -> AddFriend(fUpgradeTfMuonEmu);
    chain -> AddFriend(fCaloTowerEmu);
    chain -> AddFriend(fUpgradeEmu);
    chain -> AddFriend(fuGTEmu);
    chain -> AddFriend(fGenerator);

    std::cout<<"********friends2*********"<<std::endl;

    
    std::cout << "=== TreeReaderReco::TreeReaderReco()" << std::endl;
    //std::cout << "\tSet tree with name \"" << tree->GetName() << "\" as the principal tree." << std::endl;
    std::cout << "\tSet tree with name \"" << chain->GetName() << "\" as the principal tree." << std::endl;

  }

  //Init(tree);
  Init(chain);

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

//void TreeReaderReco::Init(TTree *tree)
void TreeReaderReco::Init(TChain *chain)
{
  InitReco(chain);
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

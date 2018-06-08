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
    virtual bool     CheckTreesExistence(const std::string SamplePath, const std::string SampleName);
};

#endif // TreeReaderReco_h

#ifdef TreeReaderReco_cxx

//TreeReaderReco::TreeReaderReco(const std::string SamplePath, const std::string SampleName, TTree *tree
TreeReaderReco::TreeReaderReco(const std::string SamplePath, const std::string SampleName, TChain *chain)
{
  std::cout << "=== TreeReaderReco::TreeReaderReco()" << std::endl;


  if (!CheckTreesExistence(SamplePath, SampleName)) exit(0);

  
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (chain == 0) {

    // Fix me - Check that tree exists *************************************
    // Creat a new TChain pointer with the TTree name
    TChain* fMainChain    = new TChain("l1EventTree/L1EventTree");
    fCaloTower        = new TChain("l1CaloTowerTree/L1CaloTowerTree");
    fUpgradeTfMuon    = new TChain("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree");
    fUpgrade          = new TChain("l1UpgradeTree/L1UpgradeTree");
    fuGT              = new TChain("l1uGTTree/L1uGTTree");
    fHO               = new TChain("l1HOTree/L1HOTree");
    fUpgradeTfMuonEmu = new TChain("l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree");
    fCaloTowerEmu     = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
    fUpgradeEmu       = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    fuGTEmu           = new TChain("l1uGTEmuTree/L1uGTTree");
    fGenerator        = new TChain("l1GeneratorTree/L1GenTree");
    

    
    // Associated the ROOT files to the TChain
    std::cout << "\tGetting ROOT files for dataset " << SampleName << " and adding them to the chain." << std::endl;
    
    // Add files to the Main Chain (Event)
    OpenFiles(SamplePath, SampleName, fMainChain);
    
    // Add files to chains (if exist)
    if (doCaloTower)           OpenFiles(SamplePath, SampleName, fCaloTower        );
    if (doUpgradeTfMuon)       OpenFiles(SamplePath, SampleName, fUpgradeTfMuon    );
    if (doUpgrade)             OpenFiles(SamplePath, SampleName, fUpgrade          );
    if (douGT)                 OpenFiles(SamplePath, SampleName, fuGT              );
    if (doHO)                  OpenFiles(SamplePath, SampleName, fHO               );
    if (doUpgradeTfMuonEmu)    OpenFiles(SamplePath, SampleName, fUpgradeTfMuonEmu );
    if (doCaloTowerEmu)        OpenFiles(SamplePath, SampleName, fCaloTowerEmu     );
    if (doUpgradeEmu)          OpenFiles(SamplePath, SampleName, fUpgradeEmu       );
    if (douGTEmu)              OpenFiles(SamplePath, SampleName, fuGTEmu           );
    if (doGenerator)           OpenFiles(SamplePath, SampleName, fGenerator        );


    // Set the Tree
    chain = fMainChain;

     
    std::cout << "=== TreeReaderReco::TreeReaderReco()" << std::endl;
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


bool TreeReaderReco::CheckTreesExistence(const std::string SamplePath, const std::string SampleName)
{
  
  // 
  doCaloTower        = true;
  doUpgradeTfMuon    = true;
  doUpgrade          = true;
  douGT              = true;
  doHO               = true;
  doUpgradeTfMuonEmu = true;
  doCaloTowerEmu     = true;
  doUpgradeEmu       = true;
  douGTEmu           = true;
  doGenerator        = true;
  doCheck            = true; // CHECK!
  //


  // Get the file name of the first root file 
  std::string fileName = GetFirstFileName(SamplePath, SampleName);
  
  // Open root file
  TFile* file = TFile::Open(fileName.c_str());


  if (file==0) return false;
  if (file->IsOpen()==0) return false;


  // Take the trees
  TTree *treeEvent            = (TTree*) file->Get("l1EventTree/L1EventTree");
  TTree *treeCaloTower        = (TTree*) file->Get("l1CaloTowerTree/L1CaloTowerTree");
  TTree *treeUpgradeTfMuon    = (TTree*)file->Get("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree");
  TTree *treeUpgrade          = (TTree*)file->Get("l1UpgradeTree/L1UpgradeTree");
  TTree *treeuGT              = (TTree*)file->Get("l1uGTTree/L1uGTTree");
  TTree *treeHO               = (TTree*)file->Get("l1HOTree/L1HOTree");
  TTree *treeUpgradeTfMuonEmu = (TTree*)file->Get("l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree");
  TTree *treeCaloTowerEmu     = (TTree*)file->Get("l1CaloTowerEmuTree/L1CaloTowerTree");
  TTree *treeUpgradeEmu       = (TTree*)file->Get("l1UpgradeEmuTree/L1UpgradeTree");
  TTree *treeuGTEmu           = (TTree*)file->Get("l1uGTEmuTree/L1uGTTree");
  TTree *treeGenerator        = (TTree*)file->Get("l1GeneratorTree/L1GenTree");
  TTree *treeCheck            = (TTree*)file->Get("l1Check/L1CheckTree"); // CHECK!

  // Check if the trees exist
  if (!treeEvent) {
    std::cout<<"L1EventTree not found .... "<<std::endl;
    return false;
  }

  if (!treeCaloTower) {
    std::cout<<"L1CaloTowerTree not found .... "<<std::endl;
    doCaloTower = false;
  }

  if (!treeUpgradeTfMuon) {
    std::cout<<"L1UpgradeTfMuonTree not found .... "<<std::endl;
    doUpgradeTfMuon = false;
  }

  if (!treeUpgrade) {
    std::cout<<"L1UpgradeTree not found .... "<<std::endl;
    doUpgrade = false;
  }

  if (!treeuGT) {
    std::cout<<"L1uGTTree not found .... "<<std::endl;
    douGT = false;
  }

  if (!treeHO) {
    std::cout<<"L1HOTree not found .... "<<std::endl;
    doHO = false;
  }

  if (!treeUpgradeTfMuonEmu) {
    std::cout<<"L1UpgradeTfMuonEmuTree not found .... "<<std::endl;
    doUpgradeTfMuonEmu = false;
  }

  if (!treeCaloTowerEmu) {
    std::cout<<"L1CaloTowerEmuTree not found .... "<<std::endl;
    doCaloTowerEmu = false;
  }

  if (!treeUpgradeEmu) {
    std::cout<<"L1UpgradeEmuTree not found .... "<<std::endl;
    doUpgradeEmu = false;
  }

  if (!treeuGTEmu) {
    std::cout<<"L1uGTEmuTree not found .... "<<std::endl;
    douGTEmu = false;
  }

  if (!treeGenerator) {
    std::cout<<"L1GeneratorTree not found .... "<<std::endl;
    doGenerator = false;
  }

  // CHECK!
  if (!treeCheck) {
    //std::cout<<"TREE EXISTENCE CHECK WORKS ! ! ! ! "<<std::endl;
    doCheck = false;
  }



  return true;
}
#endif //TreeReaderReco_cxx

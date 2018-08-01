
#ifndef TreeDefinitionBase_h
#define TreeDefinitionBase_h

#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloClusterDataFormat.h"
// #include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisCustomGeneratorDataFormat.h"
#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeTfMuonDataFormat.h"
#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisBMTFInputsDataFormat.h"
#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "../../../../../L1Trigger/L1TNtuples/interface/L1AnalysisL1HODataFormat.h"
//#include "../../../../../DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"

// ROOT
#include "TChain.h"

class TreeDefinitionBase 
{
   public:
     TChain  *fChain; //! pointer to the analysed Tree or Chain
     TChain  *fCaloTower;
     TChain  *fUpgradeTfMuon;
     TChain  *fUpgrade;
     TChain  *fuGT;
     TChain  *fHO;
     TChain  *fUpgradeTfMuonEmu;
     TChain  *fCaloTowerEmu;
     TChain  *fUpgradeEmu;
     TChain  *fuGTEmu;
     TChain  *fGenerator;
     TChain  *fTracks;
     Int_t   fCurrent;  //! Current Tree number in a Chain

     bool doCaloTower;
     bool doUpgradeTfMuon;
     bool doUpgrade;
     bool douGT;
     bool doHO;
     bool doUpgradeTfMuonEmu;
     bool doCaloTowerEmu;
     bool doUpgradeEmu;
     bool douGTEmu;
     bool doGenerator;
     bool doTracks;

     // ---May be needed for taking both CaloTowers and CaloClusters
     // L1Analysis::L1AnalysisCaloTPDataFormat      caloTP_;
     // L1Analysis::L1AnalysisL1CaloTowerDataFormat* caloTower_;
     // L1Analysis::L1AnalysisL1CaloClusterDataFormat* caloCluster_;

     virtual Bool_t  Notify();
     virtual void    Show(Long64_t entry = -1);
     virtual Int_t   GetEntry(Long64_t entry);
};

Int_t TreeDefinitionBase::GetEntry(Long64_t entry)
{
  // Read contents of entry
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Bool_t TreeDefinitionBase::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void TreeDefinitionBase::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

#endif // TreeDefinitionBase_h

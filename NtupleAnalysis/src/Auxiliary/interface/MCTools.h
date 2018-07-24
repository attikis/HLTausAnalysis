#ifndef MCTools_h
#define MCTools_h

// System
#include <cmath>
#include <iomanip>

// User
#include "../../Framework/interface/TreeDefinitionGenP.h"
#include "../interface/AuxTools.h"
#include "../interface/Table.h"
// #include "../DataFormat/src/GenParticle.C"

// ROOT
#include <TLorentzVector.h>

class MCTools : public virtual TreeDefinitionGenP
{
 public:
  virtual ~MCTools() { };

  bool IsNeutrino(const int pdgId);
  bool IsBoson(const int pdgId);
  bool IsLepton(const int pdgId);
  bool IsChargedLepton(const int pdgId);
  bool IsQuark(const int pdgId);

  
  void PrintGenp(Int_t Indx, bool bPrintHeaders = true);

  void PrintGenpFullInfo(Int_t Indx, bool bPrintHeaders = true);

  void PrintGenpMinimalInfo(Int_t Indx, bool bPrintHeaders = true);

  void PrintAllMothers(Int_t Indx);
  
  void PrintAllDaughters(Int_t Indx);
  
  void PrintAllDaughtersMinimalInfo(Int_t Indx);

  void PrintGenpDaughters(Int_t Indx, bool bPrintHeaders = true);

  TLorentzVector GetGenpP4(const Int_t Index);

  void PrintTauDecayStatistics(const unsigned int NTaus, 
			       const unsigned int Decay[][7],
			       const unsigned int Nr);

  bool IsHadronicTauDecay(Int_t Indx);
  
  // GenParticle GetGenParticle(unsigned int Index);

  void GetHadronicTauFinalDaugh(Int_t Indx, std::vector<unsigned short>& Daug);

  Int_t GetGenTauDecayMode(std::vector<unsigned short> Daug);
  
  bool RecursivelyLookForMotherId(Int_t Indx,
				  Int_t MoId,
				  const bool posn = false);

  void GetHadronicTauChargedOrNeutralPions(Int_t tauIndex, 
					   Int_t charge,
					   std::vector<unsigned short> &chargedPions);

  void GetHadronicTauNeutralPions(Int_t tauIndex, 
				  std::vector<unsigned short> &neutralPions);

  void GetHadronicTauChargedPions(Int_t tauIndex, 
				  std::vector<unsigned short> &chargedPions);

  double GetLxy(const Int_t iGenP,
		double refPoint_X = 0.0,
		double refPoint_Y = 0.0);
	      
  double GetD0Mag(const Int_t iGenP,
		  const Int_t iMom);
  
  Int_t PosOfMotherId(Int_t indx,
		      Int_t MoId,
		      const bool posn = false);
  
  bool LookForMotherId(Int_t Indx,
		       Int_t MoId,
		       const bool posn = false);
  
  bool IsFinalStateTau(const std::vector<unsigned short>& Daug);
  
  TLorentzVector GetP4(const Int_t iGenP);

  TLorentzVector GetVisibleP4(const std::vector<unsigned short>& Daug);
  
  Int_t GetLdgDaughterIndex(std::vector<unsigned short> Daug, 
			    bool bOnlyChargedDaughters=false);
  
  Double_t GetHadronicTauMaxSignalCone(std::vector<unsigned short> Daug, 
				       bool bOnlyChargedDaughters=true,
				       double minPt=0.0);
  
  Bool_t IsFinalStateHadronicTau(Int_t indx);
  
  Bool_t IsTriggerHadronicTau(const int iGenP, 
			      const bool bApplyTkAcceptanceCut, 
			      const int trueMomId, 
			      const bool bUseAbsMomId,
			      const double visEt);

private:
AuxTools tools;

};

#endif // MCTools_h

#ifndef L1TkTauParticle_h
#define L1TkTauParticle_h

// System
#include <iostream>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../DataFormat/src/L1Jet.C"
#include "../../DataFormat/src/L1Tau.C"
#include "../../DataFormat/interface/TTTrack.h"
#include "../../DataFormat/interface/TTPixelTrack.h"
#include "../../DataFormat/interface/GenParticle.h"

using namespace std;

class L1TkTauParticle{     
 public:
  // Constructors/Destructors
  L1TkTauParticle();

  L1TkTauParticle(double matchCone_dRMin,
		  double matchCone_dRMax,
		  double sigCone_dRMin,
		  double sigCone_dRMax,
		  double isoCone_dRMin,
		  double isoCone_dRMax);

  
  L1TkTauParticle(int caloTau_Index,
		  int matchTk_Index,
		  double matchTk_deltaR,
		  vector<int> sigTks_Index,
		  vector<int> isoTks_Index,
		  double vtxIso,
		  double relIso,
		  double sigCone_minDeltaR, 
		  double sigCone_maxDeltaR, 
		  double isoCone_minDeltaR, 
		  double isoCone_maxDeltaR);

  ~L1TkTauParticle() {};
  
  // Function declaration
  void InitVars_(void);
  int GetCaloTauIndex() const { return caloTau_Index_;}
  int GetMatchTk()   const { return matchTk_Index_;}
  double GetVtxIso()  const { return vtxIso_;}
  double GetRelIso()  const { return relIso_;}
  vector<int> GetSigConeTks()  const { return sigTks_Index_;}
  vector<int> GetIsoConeTks()  const { return isoTks_Index_;} 
  void SetCaloTau(int caloTau_Index) { caloTau_Index_ = caloTau_Index;}
  void SetSignalConeSize(double deltaR_min, double deltaR_max);
  void SetIsolationConeSize(double deltaR_min, double deltaR_max);
  void SetMatchTk(int matchTk_Index) { matchTk_Index_ = matchTk_Index;}
  void SetMatchTkDeltaR(double matchTk_deltaR) { matchTk_deltaR_ = matchTk_deltaR;}
  void SetMatchGenp(int matchGenp_Index, double matchGenp_deltaR);
  void SetVtxIso(double vtxIso) { vtxIso_ = vtxIso;}
  void SetRelIso(double relIso) { relIso_ = relIso;}
  // void PrintProperties(void);
  // NEW
  void SetSigConeTks(vector<int> sigTksIndices) { sigTks_Index_ = sigTksIndices;}
  void SetIsoConeTks(vector<int> isoTksIndices) { isoTks_Index_ = isoTksIndices;}
  double GetMatchConeMin(void) const{ return theMatchCone_dRMin;}
  double GetMatchConeMax(void) const{ return theMatchCone_dRMax;}
  double GetSigConeMin(void) const{ return theSigCone_dRMin;}
  double GetSigConeMax(void) const{ return theSigCone_dRMax;}
  double GetIsoConeMin(void) const{ return theIsoCone_dRMin;}
  double GetIsoConeMax(void) const{ return theIsoCone_dRMax;}
  double GetVtxIsolation(void)const { return theVtxIsolation;}
  TTTrack GetVtxIsolationTrack(void)const {return theVtxIsolationTk;}
  double GetRelIsolation(void)const { return theRelIsolation;}
  double CalculateRelIso(const double deltaZ0_max=999.9, bool bStoreValue=false, bool bInvert_deltaZ0=false, bool bUseIsoCone=false);
  double CalculateVtxIso(bool bStoreValue=false, bool bUseIsoCone=false);
  L1Tau GetCaloTau(void) const{ return theCaloTau;}
  TTTrack GetMatchingTk(void) const{ return theMatchingTk;}
  TTTrack GetSigConeLdgTk(void);
  TTTrack GetIsoConeLdgTk(void);
  TTTrack GetIsoAnnulusLdgTk(void);
  bool HasMatchingTk(void) const{return theMatchingTk.getPt() > 1.0;} 
  double GetMatchingTkDeltaR(void) const{ return theMatchingTk_dR;}  
  vector<TTTrack> GetSigConeTTTracks(void) const { return theSigConeTTTracks;}
  TLorentzVector GetSigConeTTTracksP4(void);
  TLorentzVector GetIsoConeTTTracksP4(void);
  TLorentzVector GetIsoAnnulusTTTracksP4(void);
  vector<TTTrack> GetIsoConeTTTracks(void) const { return theIsoConeTTTracks;}
  vector<TTTrack> GetIsoAnnulusTTTracks(void) const { return theIsoAnnulusTTTracks;}
  vector<TTPixelTrack> GetSigConeTTPixelTracks(void) const { return theIsoConeTTPixelTracks;}
  vector<TTPixelTrack> GetIsoConeTTPixelTracks(void) const { return theSigConeTTPixelTracks;}
  GenParticle GetMatchingGenParticle(void) const { return theMatchingGenParticle;}
  double GetMatchingGenParticleDeltaR(void) const { return theMatchingGenParticle_dR;}
  bool HasMatchingGenParticle(void) const{return GetMatchingGenParticleDeltaR() < 999.9;}
  int GetNProngs(void) const{ return theNProngs;}
  //        
  void SetCaloTau(L1Tau CaloTau){theCaloTau = CaloTau;}
  void SetMatchingTk(TTTrack matchTk){theMatchingTk = matchTk;}
  void SetMatchTkDeltaRNew(double matchTk_dR){theMatchingTk_dR = matchTk_dR;}
  void SetMatchTkDeltaRMin(double matchTk_dRMin){theMatchCone_dRMin = matchTk_dRMin;}
  void SetMatchTkDeltaRMax(double matchTk_dRMax){theMatchCone_dRMax = matchTk_dRMax;}
  void SetVtxIsolationTrack(TTTrack vtxIsoTk){theVtxIsolationTk = vtxIsoTk;}
  void SetVtxIsolation(double isoValue){theVtxIsolation = isoValue;}
  void SetRelIsolation(double isoValue){theRelIsolation = isoValue;}
  void SetSigConeMinDeltaR(double dRMin){theSigCone_dRMin = dRMin;}
  void SetSigConeMaxDeltaR(double dRMax){theSigCone_dRMax = dRMax;}
  void SetIsoConeMinDeltaR(double dRMin){theIsoCone_dRMin = dRMin;}
  void SetIsoConeMaxDeltaR(double dRMax){theIsoCone_dRMax = dRMax;}
  void SetSigConeTracks(vector<TTTrack> sigConeTks){theSigConeTTTracks = sigConeTks;}
  void SetSigConeTracks(vector<TTPixelTrack> sigConeTks){theIsoConeTTPixelTracks = sigConeTks;}
  void SetIsoConeTracks(vector<TTTrack> isoConeTks){theIsoConeTTTracks = isoConeTks;}
  void SetIsoAnnulusTracks(vector<TTTrack> isoAnnulusTks){theIsoAnnulusTTTracks = isoAnnulusTks;}
  void SetIsoConeTracks(vector<TTPixelTrack> isoConeTks){theSigConeTTPixelTracks = isoConeTks;}
  void SetMatchingGenParticle(GenParticle genP){ theMatchingGenParticle = genP;}
  void SetMatchingGenParticleDeltaR(double dR){ theMatchingGenParticle_dR = dR;}
  
  void PrintTTTracks(vector<TTTrack> theTracks,
		     string theTrackType);
  
  void PrintTTPixelTracks(vector<TTPixelTrack> theTracks,
			  string theTrackType);

  void PrintProperties(bool bPrintCaloTau=false,
		       bool bPrintMatchTk=false,
		       bool bPrintSigConeTks=false,
		       bool bPrintIsoConeTks=false,
		       bool bPrintMatchGenParticle=false);
  
  // Variable declaration
  int caloTau_Index_;
  int matchTk_Index_;
  double matchTk_deltaR_;
  int matchGenp_Index_;
  double matchGenp_deltaR_;
  vector<int> sigTks_Index_;
  vector<int> isoTks_Index_;
  double vtxIso_;
  double relIso_;
  double sigCone_minDeltaR_;
  double sigCone_maxDeltaR_;
  double isoCone_minDeltaR_;
  double isoCone_maxDeltaR_;
  // NEW
  L1Tau theCaloTau;
  TTTrack theMatchingTk;
  double theMatchingTk_dR;
  double theMatchCone_dRMin;
  double theMatchCone_dRMax;
  double theSigCone_dRMin;
  double theSigCone_dRMax;
  double theIsoCone_dRMin;
  double theIsoCone_dRMax;
  vector<TTTrack> theSigConeTTTracks;
  vector<TTPixelTrack> theSigConeTTPixelTracks;
  TLorentzVector theSigConeTTTracksP4;
  vector<TTTrack> theIsoConeTTTracks;
  vector<TTTrack> theIsoAnnulusTTTracks;
  vector<TTPixelTrack> theIsoConeTTPixelTracks;
  TLorentzVector theIsoConeTTTracksP4;
  TLorentzVector theIsoAnnulusTTTracksP4;
  double theVtxIsolation;
  double theRelIsolation;
  TTTrack theVtxIsolationTk;
  GenParticle theMatchingGenParticle;
  double theMatchingGenParticle_dR;
  int theNProngs;
  
 private:
  AuxTools auxTools;
  void SetSigConeTTTracksP4_(void);
  void SetIsoConeTTTracksP4_(void);
  void SetIsoAnnulusTTTracksP4_(void);

  
};

#endif

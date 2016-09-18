#ifndef Tracks_h
#define Tracks_h

#include <iostream>

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/MCTools.C"
#include "../utilities/AuxTools.C"
#include "../utilities/HistoTools.C"
#include "../utilities/L1Tracks.C"

class Tracks : public TreeAnalyserMC{
 public:
 Tracks(const string SamplePath,
	const std::string SampleName,
	const std::string text_, 
	const int maxEvents_ = -1, 
	TTree* tree=0) : 
  TreeAnalyserMC("Tracks", SamplePath, SampleName, text_, maxEvents_, tree) {};
  virtual void Loop();
  
 private:
  void Initialise(void);

  void BookHistos(const string tk_WP,
		  const string tk_Type);

  void FillHistos(const unsigned int tk_index,
		  const string tk_Type);

  // Objects/Variables
  AuxTools auxTools_;
  HistoTools histoTools_;
  string track_WP;
  
  // Histograms
  TH1D* hL1Tks_Multiplicity;
  TH1D* hL1Tks_Pt          ;
  TH1D* hL1Tks_Eta         ;
  TH1D* hL1Tks_EtaAbs      ;
  TH1D* hL1Tks_Phi         ;
  TH1D* hL1Tks_Charge      ;
  TH1D* hL1Tks_NPixHits    ;
  TH1D* hL1Tks_NStubs      ;
  TH1D* hL1Tks_NPsStubs    ;
  TH1D* hL1Tks_StubPtCons  ;
  TH1D* hL1Tks_POCAz       ;
  TH1D* hL1Tks_POCAzAbs    ;
  TH1D* hL1Tks_d0          ;
  TH1D* hL1Tks_d0Abs       ;
  TH1D* hL1Tks_d0Sig       ;
  TH1D* hL1Tks_ChiSq       ;
  TH1D* hL1Tks_RedChiSq    ;
  TH1D* hL1Tks_RInv        ;
  TH1D* hL1Tks_SigmaRInv   ;
  TH1D* hL1Tks_SigmaPhi0   ;
  TH1D* hL1Tks_SigmaD0     ;
  TH1D* hL1Tks_SigmaT      ;
  TH1D* hL1Tks_SigmaZ0     ;
  //
  TH1D* hL1Tks_d0_TkInFirstQuad;
  TH1D* hL1Tks_d0_TkInSecondQuad;
  TH1D* hL1Tks_d0_TkInThirdQuad;
  TH1D* hL1Tks_d0_TkInFourthQuad;
  
  TH1D* hL1PixTks_Multiplicity;
  TH1D* hL1PixTks_Pt          ;
  TH1D* hL1PixTks_Eta         ;
  TH1D* hL1PixTks_EtaAbs      ;
  TH1D* hL1PixTks_Phi         ;
  TH1D* hL1PixTks_Charge      ;
  TH1D* hL1PixTks_NPixHits    ;
  TH1D* hL1PixTks_NStubs      ;
  TH1D* hL1PixTks_NPsStubs    ;
  TH1D* hL1PixTks_StubPtCons  ;
  TH1D* hL1PixTks_POCAz       ;
  TH1D* hL1PixTks_POCAzAbs    ;
  TH1D* hL1PixTks_d0          ;
  TH1D* hL1PixTks_d0Abs       ;
  TH1D* hL1PixTks_d0Sig       ;
  TH1D* hL1PixTks_ChiSq       ;
  TH1D* hL1PixTks_RedChiSq    ;
  TH1D* hL1PixTks_RInv        ;
  TH1D* hL1PixTks_SigmaRInv   ;
  TH1D* hL1PixTks_SigmaPhi0   ;
  TH1D* hL1PixTks_SigmaD0     ;
  TH1D* hL1PixTks_SigmaT      ;
  TH1D* hL1PixTks_SigmaZ0     ;

};

#endif // Tracks_h

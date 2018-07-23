#ifndef PixelTracks_cxx
#define PixelTracks_cxx

#include "PixelTracks.h"
#include "../utilities/constants.h"
#include "TFitResult.h"
#include "TF1.h"


//****************************************************************************
void PixelTracks::Loop()
//****************************************************************************
{
  
  if (fChain == 0) return;
  Long64_t nEntries = (MaxEvents == -1) ? fChain->GetEntries() : std::min((int)fChain->GetEntries(), MaxEvents);
  
  ///////////////////////////////////////////////////////////
  /// Initialisations
  ///////////////////////////////////////////////////////////
  Long64_t nbytes = 0, nb = 0;
  L1PixelTrackFit f(3.8112); // Bz in Tesla
  
  // auxTools.PrintPSets(fChain);
  
  ///////////////////////////////////////////////////////////
  // For-loop: Entries
  ///////////////////////////////////////////////////////////
  for (Int_t jentry=0; jentry < nEntries; jentry++){
      
    // Init loop variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    ///////////////////////////////////////////////////////////
    // TTPixelTracks
    ///////////////////////////////////////////////////////////
    const int nPixTks = (int) L1PixTks_Pt->size();
    // For-loop: TTPixelTracks
    for (Int_t pixTk_Index = 0; pixTk_Index < nPixTks; pixTk_Index++) {

      if (0) s->PrintPixTrackProperties(pixTk_Index, false);
      
      // Pixel Hits
      vector<double> pixHits_X    = L1PixTks_PixHits_X->at(pixTk_Index);
      vector<double> pixHits_Y    = L1PixTks_PixHits_Y->at(pixTk_Index);
      vector<double> pixHits_Z    = L1PixTks_PixHits_Z->at(pixTk_Index);
      vector<double> pixHits_R    = L1PixTks_PixHits_R->at(pixTk_Index);
      vector<double> pixHits_Phi  = L1PixTks_PixHits_Phi->at(pixTk_Index);
      vector< int >  pixHits_Type = L1PixTks_PixHits_Type->at(pixTk_Index);
      const int nPixHits = (int) pixHits_X.size();

      // Candidate Pixel Hits
      vector<double> candPixHits_X    = L1PixTks_CandPixHits_X->at(pixTk_Index);
      vector<double> candPixHits_Y    = L1PixTks_CandPixHits_Y->at(pixTk_Index);
      vector<double> candPixHits_Z    = L1PixTks_CandPixHits_Z->at(pixTk_Index);
      vector<double> candPixHits_R    = L1PixTks_CandPixHits_R->at(pixTk_Index);
      vector<double> candPixHits_Phi  = L1PixTks_CandPixHits_Phi->at(pixTk_Index);
      vector< int >  candPixHits_Type = L1PixTks_CandPixHits_Type->at(pixTk_Index);
      const int nCandPixHits = (int) candPixHits_X.size();  
      
      // Corresponding TTTrack
      const int tk_Index  = L1PixTks_TTTrackIndex->at(pixTk_Index);
      double tk_Pt        = L1Tks_Pt ->at(tk_Index);
      double tk_Eta       = L1Tks_Eta->at(tk_Index);
      double tk_Phi       = L1Tks_Phi->at(tk_Index);
      double tk_RInv      = L1Tks_RInv->at(tk_Index);
      double tk_t         = sinh(tk_Eta);
      double tk_d0        = s->GetTrackD0 (tk_Index);
      double tk_z0        = L1Tks_POCAz->at(tk_Index);
      double tk_ChiSq     = L1Tks_ChiSquared->at(tk_Index);
      double tk_RedChiSq  = s->GetTrackRedChiSq(tk_Index);

      
      ///////////////////////////////////////////////////////////
      // Pixel-Hit-Fit on TTTracks
      ///////////////////////////////////////////////////////////
      double sigmarinv, sigmaphi0, sigmad0, sigmat, sigmaz0;
      TTPixelTrack pixFitTk =  f.FitPixelTrack(tk_RInv, tk_Phi, tk_d0, tk_t, tk_z0, candPixHits_X, candPixHits_Y, candPixHits_Z, candPixHits_Type);      
      pixFitTk.PrintProperties();

    } // for (Int_t iTk = 0; iTk < nL1Tks; iTk++) {

    
    auxTools.ProgressBar(jentry, nEntries, 100, 150);
    
  }// For-loop: Entries

  cout << "\n" << endl;
  outFile->cd();
  outFile->Write();

}

#endif // PixelTracks_cxx

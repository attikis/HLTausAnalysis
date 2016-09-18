#ifndef Tracks_cxx
#define Tracks_cxx

#include "Tracks.h"
#include "../utilities/constants.h"

//#define bDebug
//****************************************************************************
void Tracks::Initialise(void)
//****************************************************************************
{

  track_WP = "MediumWP";  
  BookHistos(track_WP, "TTPixelTracks");
  BookHistos(track_WP, "TTTracks");

  return;
}


//****************************************************************************
void Tracks::Loop(void)
//****************************************************************************
{
  
  if (fChain == 0) return;
  
  Long64_t nEntries      = (MaxEvents == -1) ? fChain->GetEntries() : std::min((int)fChain->GetEntries(), MaxEvents);
  Int_t  NPBarDivisions  = 100;
  Int_t PBarWidth        = 150;
  
  std::cout << "Analyzing " << nEntries << " events.\n";

  Initialise();
  // Counters
  unsigned int nTks_Tot           = 0;
  unsigned int nTks_TightWP_Tot   = 0;
  unsigned int nTks_MediumWP_Tot  = 0;
  unsigned int nTks_LooseWP_Tot   = 0;
  unsigned int nTks_vLooseWP_Tot  = 0;
  unsigned int nTks_NoWP_Tot      = 0;
  
  unsigned int nPixTks_Tot           = 0;
  unsigned int nPixTks_TightWP_Tot   = 0;
  unsigned int nPixTks_MediumWP_Tot  = 0;
  unsigned int nPixTks_LooseWP_Tot   = 0;
  unsigned int nPixTks_vLooseWP_Tot  = 0;
  unsigned int nPixTks_NoWP_Tot      = 0;

  ///////////////////////////////////////////////////////////
  // For-loop: Entries
  ///////////////////////////////////////////////////////////
  Long64_t nbytes = 0, nb = 0;
  for (Int_t jentry=0; jentry < nEntries; jentry++){

    // Init loop variables
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    // Sanity checks
    auxTools_.EnsureVectorIsSorted(*L1Tks_Pt    , true);
    auxTools_.EnsureVectorIsSorted(*L1PixTks_Pt , true);
    
    // Counters
    unsigned int nTks           = 0;
    unsigned int nTks_TightWP   = 0;
    unsigned int nTks_MediumWP  = 0;
    unsigned int nTks_LooseWP   = 0;
    unsigned int nTks_vLooseWP  = 0;
    unsigned int nTks_NoWP      = 0;

    unsigned int nPixTks           = 0;
    unsigned int nPixTks_TightWP   = 0;
    unsigned int nPixTks_MediumWP  = 0;
    unsigned int nPixTks_LooseWP   = 0;
    unsigned int nPixTks_vLooseWP  = 0;
    unsigned int nPixTks_NoWP      = 0;

    ///////////////////////////////////////////////////////////
    // TTTracks
    ///////////////////////////////////////////////////////////
    const Int_t nL1Tks = L1Tks_Pt->size();
    Table d0Properties("Entry | Track # | POCA-x | POCA-y | Phi0 | Phi-d0 | d0 | Sgn(d0) | Charge", "Text");
    int i = 0;
    
    // For-loop: L1Tks
    for (Int_t j = 0; j < nL1Tks; j++) {

      // s->PrintTrackProperties(j);
      if (s->SelectTracks(j, track_WP  ) ) FillHistos(j, "TTTracks");
      if (s->SelectTracks(j, "TightWP" ) ) nTks_TightWP++;
      if (s->SelectTracks(j, "MediumWP") ) nTks_MediumWP++;
      if (s->SelectTracks(j, "LooseWP" ) ) nTks_LooseWP++;
      if (s->SelectTracks(j, "vLooseWP") ) nTks_vLooseWP++;
      if (s->SelectTracks(j, "NoWP"    ) ) nTks_NoWP++;
      if (j == 50) break;
      
      string tk_Q = "-";
      double tk_Charge = L1Tks_Charge->at(j);
      if(tk_Charge > 0) tk_Q = "+";
      
      string tk_d0 = "-";
      double tk_dZero = s->GetTrackD0(j);
      if(tk_dZero >= 0) tk_d0 = "+";
      
      string tk_myd0 = "-";
      double tk_mydZero = s->GetTrackD0Sign(j);
      if(tk_mydZero >= 0) tk_myd0 = "+";
      
      // Print d0 Properties     
      // if ( tk_dZero * tk_mydZero < 0 ){
      d0Properties.AddRowColumn(i, auxTools_.ToString(i) );
      d0Properties.AddRowColumn(i, auxTools_.ToString(j) );
      d0Properties.AddRowColumn(i, auxTools_.ToString(L1Tks_POCAx->at(j)  ) );
      d0Properties.AddRowColumn(i, auxTools_.ToString(L1Tks_POCAy->at(j)  ) );
      d0Properties.AddRowColumn(i, auxTools_.ToString(L1Tks_Phi->at(j)*180/PI    ) );
      d0Properties.AddRowColumn(i, auxTools_.ToString(s->GetTrackD0Phi(j)*180/PI ) );
      d0Properties.AddRowColumn(i, tk_d0);
      d0Properties.AddRowColumn(i, tk_myd0);
      d0Properties.AddRowColumn(i, tk_Q );
      i++;
      //}
      
    } // For-loop: L1Tks
    d0Properties.Print();
    s->GetTrackD0Sign(31);
    s->GetTrackD0Sign(37);
    s->GetTrackD0Sign(38);

    ///////////////////////////////////////////////////////////
    // TTPixelTracks
    ///////////////////////////////////////////////////////////
    const Int_t nL1PixTks = L1PixTks_Pt->size();

    // For-loop: L1PixTks
    for (Int_t k = 0; k < nL1PixTks; k++) {
      
      if (s->SelectPixTracks(k, track_WP  ) ) FillHistos(k, "TTPixelTracks");
      if (s->SelectPixTracks(k, "TightWP" ) ) nPixTks_TightWP++;
      if (s->SelectPixTracks(k, "MediumWP") ) nPixTks_MediumWP++;
      if (s->SelectPixTracks(k, "LooseWP" ) ) nPixTks_LooseWP++;
      if (s->SelectPixTracks(k, "vLooseWP") ) nPixTks_vLooseWP++;
      if (s->SelectPixTracks(k, "NoWP"    ) ) nPixTks_NoWP++;
      
    } // For-loop: L1Tks
    
    
    // Fill Multiplicity histograms
    if ( track_WP.compare("TightWP") == 0 ){
      nTks    = nTks_TightWP;
      nPixTks = nPixTks_TightWP;
    }
    else if ( track_WP.compare("MediumWP") == 0 ){
      nTks    = nTks_MediumWP;
      nPixTks = nPixTks_MediumWP;
    }
    else if ( track_WP.compare("LooseWP") == 0 ){
      nTks    = nTks_LooseWP;
      nPixTks = nPixTks_LooseWP;
    }
    else if ( track_WP.compare("vLooseWP") == 0 ){
      nTks    = nTks_vLooseWP;
      nPixTks = nPixTks_vLooseWP;
    }
    else if ( track_WP.compare("NoWP") == 0 ){
      nTks    = nTks_NoWP;
      nPixTks = nPixTks_NoWP;
    }
    else{
      cout << "E R R O R ! Tracks::Loop(...). EXIT" << endl;
      exit(1);
    }

    hL1Tks_Multiplicity   ->Fill(nTks);
    hL1PixTks_Multiplicity->Fill(nPixTks);

    Table tkProperties("Event | Tk-Type | no-WP | vLoose-WP | Loose-WP | Medium-WP | Tight-WP", "Text");
    // Add table rows
    tkProperties.AddRowColumn(0, auxTools_.ToString(jentry) );
    tkProperties.AddRowColumn(0, "TTTracks" );
    tkProperties.AddRowColumn(0, auxTools_.ToString(nTks_NoWP) );
    tkProperties.AddRowColumn(0, auxTools_.ToString(nTks_vLooseWP) );
    tkProperties.AddRowColumn(0, auxTools_.ToString(nTks_LooseWP) );
    tkProperties.AddRowColumn(0, auxTools_.ToString(nTks_MediumWP) );
    tkProperties.AddRowColumn(0, auxTools_.ToString(nTks_TightWP) );
    //
    tkProperties.AddRowColumn(1, auxTools_.ToString(jentry) );
    tkProperties.AddRowColumn(1, "TTPixelTracks" );
    tkProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_NoWP) );
    tkProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_vLooseWP) );
    tkProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_LooseWP) );
    tkProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_MediumWP) );
    tkProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_TightWP) );
    // tkProperties.Print();

    // Update counters
    nTks_Tot           += nTks;
    nTks_TightWP_Tot   += nTks_TightWP;
    nTks_MediumWP_Tot  += nTks_MediumWP;
    nTks_LooseWP_Tot   += nTks_LooseWP;
    nTks_vLooseWP_Tot  += nTks_vLooseWP;
    nTks_NoWP_Tot      += nTks_NoWP;
    
    nPixTks_Tot           += nPixTks;
    nPixTks_TightWP_Tot   += nPixTks_TightWP;
    nPixTks_MediumWP_Tot  += nPixTks_MediumWP;
    nPixTks_LooseWP_Tot   += nPixTks_LooseWP;
    nPixTks_vLooseWP_Tot  += nPixTks_vLooseWP;
    nPixTks_NoWP_Tot      += nPixTks_NoWP;

    
#ifndef bDebug
    auxTools_.ProgressBar(jentry, nEntries, NPBarDivisions, PBarWidth);
#endif
  }// For-loop: Entries
  
    Table tkTotProperties("Tk-Type | no-WP | vLoose-WP | Loose-WP | Medium-WP | Tight-WP", "Text");
    // Add table rows
    tkTotProperties.AddRowColumn(0, "TTTracks" );
    tkTotProperties.AddRowColumn(0, auxTools_.ToString(nTks_NoWP_Tot) );
    tkTotProperties.AddRowColumn(0, auxTools_.ToString(nTks_vLooseWP_Tot) );
    tkTotProperties.AddRowColumn(0, auxTools_.ToString(nTks_LooseWP_Tot) );
    tkTotProperties.AddRowColumn(0, auxTools_.ToString(nTks_MediumWP_Tot) );
    tkTotProperties.AddRowColumn(0, auxTools_.ToString(nTks_TightWP_Tot) );

    //
    tkTotProperties.AddRowColumn(1, "TTPixelTracks" );
    tkTotProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_NoWP_Tot) );
    tkTotProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_vLooseWP_Tot) );
    tkTotProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_LooseWP_Tot) );
    tkTotProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_MediumWP_Tot) );
    tkTotProperties.AddRowColumn(1, auxTools_.ToString(nPixTks_TightWP_Tot) );
    // tkTotProperties.Print();

  // Write to file
  cout << "\n" << endl;
  outFile->cd();
  outFile->Write();

}


//****************************************************************************
void Tracks::BookHistos(const string tk_WP,
			const string tk_Type)
//****************************************************************************
{

  string preName = tk_Type + string("_") + tk_WP + string("_");
  if ( tk_Type.compare("TTPixelTracks") == 0 ) 
    {

      histoTools_.BookHisto_1D(hL1PixTks_Multiplicity , (preName + "Multiplicity").c_str(),   300 ,   -0.5   ,  +299.5    );
      histoTools_.BookHisto_1D(hL1PixTks_Pt           , (preName + "Pt"          ).c_str(),   200 ,   +0.0   ,  +200.0    );
      histoTools_.BookHisto_1D(hL1PixTks_Eta          , (preName + "Eta"         ).c_str(),   600 ,   -3.0   ,    +3.0    );
      histoTools_.BookHisto_1D(hL1PixTks_EtaAbs       , (preName + "EtaAbs"      ).c_str(),   300 ,   +0.0   ,    +3.0    );
      histoTools_.BookHisto_1D(hL1PixTks_Phi          , (preName + "Phi"         ).c_str(),   800 ,   -4.0   ,    +4.0    );
      histoTools_.BookHisto_1D(hL1PixTks_Charge       , (preName + "Charge"      ).c_str(),     7 ,   -3.5   ,    +3.5    );
      histoTools_.BookHisto_1D(hL1PixTks_NPixHits     , (preName + "NPixHits"    ).c_str(),    10 ,   -0.5   ,    +9.5    );
      histoTools_.BookHisto_1D(hL1PixTks_NStubs       , (preName + "NStubs"      ).c_str(),    11 ,   -0.5   ,   +10.5    );
      histoTools_.BookHisto_1D(hL1PixTks_NPsStubs     , (preName + "NPsStubs"    ).c_str(),    11 ,   -0.5   ,   +10.5    );
      histoTools_.BookHisto_1D(hL1PixTks_StubPtCons   , (preName + "StubPtCons"  ).c_str(),   200 ,   +0.0   ,  +200.0    );
      histoTools_.BookHisto_1D(hL1PixTks_POCAz        , (preName + "POCAz"       ).c_str(),   600 ,  -30.0   ,   +30.0    );
      histoTools_.BookHisto_1D(hL1PixTks_POCAzAbs     , (preName + "POCAzAbs"    ).c_str(),   300 ,   +0.0   ,   +30.0    );
      histoTools_.BookHisto_1D(hL1PixTks_d0           , (preName + "d0"          ).c_str(),  2000 ,  -10.0   ,  + 10.0    );
      histoTools_.BookHisto_1D(hL1PixTks_d0Abs        , (preName + "d0Abs"       ).c_str(),  1000 ,   +0.0   ,   +10.0    );
      histoTools_.BookHisto_1D(hL1PixTks_d0Sig        , (preName + "d0Sig"       ).c_str(),   400 , -100.0   ,  +100.0    );
      histoTools_.BookHisto_1D(hL1PixTks_ChiSq        , (preName + "ChiSq"       ).c_str(), 50000 ,   +0.0   ,+50000.0    );
      histoTools_.BookHisto_1D(hL1PixTks_RedChiSq     , (preName + "RedChiSq"    ).c_str(), 50000 ,   +0.0   ,+50000.0    );
      histoTools_.BookHisto_1D(hL1PixTks_RInv         , (preName + "RInv"        ).c_str(), 20000 ,   -1.0   ,    +1.0    );
      histoTools_.BookHisto_1D(hL1PixTks_SigmaRInv    , (preName + "SigmaRInv"   ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1PixTks_SigmaPhi0    , (preName + "SigmaPhi0"   ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1PixTks_SigmaD0      , (preName + "SigmaD0"     ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1PixTks_SigmaT       , (preName + "SigmaT"      ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1PixTks_SigmaZ0      , (preName + "SigmaZ0"     ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      
    }
  else if ( tk_Type.compare("TTTracks") == 0 )
    {

      histoTools_.BookHisto_1D(hL1Tks_Multiplicity , (preName + "Multiplicity").c_str(),   300 ,   -0.5   ,  +299.5    );
      histoTools_.BookHisto_1D(hL1Tks_Pt           , (preName + "Pt"          ).c_str(),   200 ,   +0.0   ,  +200.0    );
      histoTools_.BookHisto_1D(hL1Tks_Eta          , (preName + "Eta"         ).c_str(),   600 ,   -3.0   ,    +3.0    );
      histoTools_.BookHisto_1D(hL1Tks_EtaAbs       , (preName + "EtaAbs"      ).c_str(),   300 ,   +0.0   ,    +3.0    );
      histoTools_.BookHisto_1D(hL1Tks_Phi          , (preName + "Phi"         ).c_str(),   800 ,   -4.0   ,    +4.0    );
      histoTools_.BookHisto_1D(hL1Tks_Charge       , (preName + "Charge"      ).c_str(),     7 ,   -3.5   ,    +3.5    );
      histoTools_.BookHisto_1D(hL1Tks_NPixHits     , (preName + "NPixHits"    ).c_str(),    10 ,   -0.5   ,    +9.5    );
      histoTools_.BookHisto_1D(hL1Tks_NStubs       , (preName + "NStubs"      ).c_str(),    11 ,   -0.5   ,   +10.5    );
      histoTools_.BookHisto_1D(hL1Tks_NPsStubs     , (preName + "NPsStubs"    ).c_str(),    11 ,   -0.5   ,   +10.5    );
      histoTools_.BookHisto_1D(hL1Tks_StubPtCons   , (preName + "StubPtCons"  ).c_str(),   200 ,   +0.0   ,  +200.0    );
      histoTools_.BookHisto_1D(hL1Tks_POCAz        , (preName + "POCAz"       ).c_str(),   600 ,  -30.0   ,   +30.0    );
      histoTools_.BookHisto_1D(hL1Tks_POCAzAbs     , (preName + "POCAzAbs"    ).c_str(),   300 ,   +0.0   ,   +30.0    );
      histoTools_.BookHisto_1D(hL1Tks_d0           , (preName + "d0"          ).c_str(),  2000 ,  -10.0   ,  + 10.0    );
      histoTools_.BookHisto_1D(hL1Tks_d0Abs        , (preName + "d0Abs"       ).c_str(),  1000 ,   +0.0   ,   +10.0    );
      histoTools_.BookHisto_1D(hL1Tks_d0Sig        , (preName + "d0Sig"       ).c_str(),   400 , -100.0   ,  +100.0    );
      histoTools_.BookHisto_1D(hL1Tks_ChiSq        , (preName + "ChiSq"       ).c_str(), 50000 ,   +0.0   ,+50000.0    );
      histoTools_.BookHisto_1D(hL1Tks_RedChiSq     , (preName + "RedChiSq"    ).c_str(), 50000 ,   +0.0   ,+50000.0    );
      histoTools_.BookHisto_1D(hL1Tks_RInv         , (preName + "RInv"        ).c_str(), 20000 ,   -1.0   ,    +1.0    );
      histoTools_.BookHisto_1D(hL1Tks_SigmaRInv    , (preName + "SigmaRInv"   ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1Tks_SigmaPhi0    , (preName + "SigmaPhi0"   ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1Tks_SigmaD0      , (preName + "SigmaD0"     ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1Tks_SigmaT       , (preName + "SigmaT"      ).c_str(),  5000 ,   +0.0   ,    +0.05   );
      histoTools_.BookHisto_1D(hL1Tks_SigmaZ0      , (preName + "SigmaZ0"     ).c_str(),  5000 ,   +0.0   ,    +0.05   );

      histoTools_.BookHisto_1D(hL1Tks_d0_TkInFirstQuad , (preName + "d0_TkInFirstQuad" ).c_str(),  2000 ,  -10.0   ,  + 10.0    );
      histoTools_.BookHisto_1D(hL1Tks_d0_TkInSecondQuad, (preName + "d0_TkInSecondQuad").c_str(),  2000 ,  -10.0   ,  + 10.0    );
      histoTools_.BookHisto_1D(hL1Tks_d0_TkInThirdQuad , (preName + "d0_TkInThirdQuad" ).c_str(),  2000 ,  -10.0   ,  + 10.0    );
      histoTools_.BookHisto_1D(hL1Tks_d0_TkInFourthQuad, (preName + "d0_TkInFourthQuad").c_str(),  2000 ,  -10.0   ,  + 10.0    );
    }
  else
    {
      cout << "E R R O R ! Tracks::BookHistos(...). EXIT" << endl;
      exit(1);
    }


  
  return;
}


//****************************************************************************
void Tracks::FillHistos(const unsigned int tk_index,
			const string tk_Type)		
//****************************************************************************
{

  if ( tk_Type.compare("TTPixelTracks") == 0 ) 
    {
      const unsigned int refTk_index =  L1PixTks_TTTrackIndex->at(tk_index);
      hL1PixTks_Pt          ->Fill( L1PixTks_Pt->at(tk_index)                );
      hL1PixTks_Eta         ->Fill( L1PixTks_Eta->at(tk_index)               );
      hL1PixTks_EtaAbs      ->Fill( abs(L1PixTks_Eta->at(tk_index) )         );
      hL1PixTks_Phi         ->Fill( L1PixTks_Phi->at(tk_index)               );
      hL1PixTks_Charge      ->Fill( L1PixTks_Charge->at(tk_index)            );
      hL1PixTks_NPixHits    ->Fill( L1PixTks_NPixHits->at(tk_index)          );
      hL1PixTks_NStubs      ->Fill( s->GetNumOfStubs(refTk_index)            );
      hL1PixTks_NPsStubs    ->Fill( s->GetNumOfPSStubs(refTk_index)          );
      hL1PixTks_StubPtCons  ->Fill( L1Tks_StubPtConsistency->at(refTk_index) );
      hL1PixTks_POCAz       ->Fill( L1PixTks_POCAz->at(tk_index)             );
      hL1PixTks_POCAzAbs    ->Fill( abs( L1PixTks_POCAz->at(tk_index) )      );
      hL1PixTks_d0          ->Fill( s->GetPixTrackD0(tk_index)               );
      hL1PixTks_d0Abs       ->Fill( abs(s->GetPixTrackD0(tk_index) )         );
      hL1PixTks_d0Sig       ->Fill( s->GetPixTrackD0Sig(tk_index)            );
      hL1PixTks_ChiSq       ->Fill( L1PixTks_ChiSquared->at(tk_index)        );
      hL1PixTks_RedChiSq    ->Fill( s->GetPixTrackRedChiSq(tk_index)         );
      hL1PixTks_RInv        ->Fill( L1PixTks_RInv->at(tk_index)              );
      hL1PixTks_SigmaRInv   ->Fill( L1PixTks_SigmaRInv->at(tk_index)         );
      hL1PixTks_SigmaPhi0   ->Fill( L1PixTks_SigmaPhi0->at(tk_index)         );
      hL1PixTks_SigmaD0     ->Fill( L1PixTks_SigmaD0->at(tk_index)           );
      hL1PixTks_SigmaT      ->Fill( L1PixTks_SigmaT->at(tk_index)            );
      hL1PixTks_SigmaZ0     ->Fill( L1PixTks_SigmaZ0->at(tk_index)           );
    }
  else if ( tk_Type.compare("TTTracks") == 0 )
    {
    hL1Tks_Pt          ->Fill( L1Tks_Pt->at(tk_index)                );
    hL1Tks_Eta         ->Fill( L1Tks_Eta->at(tk_index)               );
    hL1Tks_EtaAbs      ->Fill( abs(L1Tks_Eta->at(tk_index) )         );
    hL1Tks_Phi         ->Fill( L1Tks_Phi->at(tk_index)               );
    hL1Tks_Charge      ->Fill( L1Tks_Charge->at(tk_index)            );
    hL1Tks_NPixHits    ->Fill( 0.0                                   );
    hL1Tks_NStubs      ->Fill( s->GetNumOfStubs(tk_index)            );
    hL1Tks_NPsStubs    ->Fill( s->GetNumOfPSStubs(tk_index)          );
    hL1Tks_StubPtCons  ->Fill( L1Tks_StubPtConsistency->at(tk_index) );
    hL1Tks_POCAz       ->Fill( L1Tks_POCAz->at(tk_index)             );
    hL1Tks_POCAzAbs    ->Fill( abs(L1Tks_POCAz->at(tk_index) )       );
    hL1Tks_d0          ->Fill( s->GetTrackD0(tk_index)               );
    hL1Tks_d0Abs       ->Fill( abs(s->GetTrackD0(tk_index) )         );
    hL1Tks_d0Sig       ->Fill( 0.0                                   );
    hL1Tks_ChiSq       ->Fill( L1Tks_ChiSquared->at(tk_index)        );
    hL1Tks_RedChiSq    ->Fill( L1Tks_RedChiSquared->at(tk_index)     );
    hL1Tks_RInv        ->Fill( L1Tks_RInv->at(tk_index)              );
    hL1Tks_SigmaRInv   ->Fill( 0.0                                   );
    hL1Tks_SigmaPhi0   ->Fill( 0.0                                   );
    hL1Tks_SigmaD0     ->Fill( 0.0                                   );
    hL1Tks_SigmaT      ->Fill( 0.0                                   );
    hL1Tks_SigmaZ0     ->Fill( 0.0                                   );
         
    if (L1Tks_Charge->at(tk_index) > 0){
      
      if (L1Tks_POCAx->at(tk_index) > 0 && L1Tks_POCAy->at(tk_index) > 0) hL1Tks_d0_TkInFirstQuad ->Fill( s->GetTrackD0(tk_index) );
      else if (L1Tks_POCAx->at(tk_index) > 0 && L1Tks_POCAy->at(tk_index) < 0) hL1Tks_d0_TkInFourthQuad->Fill( s->GetTrackD0(tk_index) );
      else if (L1Tks_POCAx->at(tk_index) < 0 && L1Tks_POCAy->at(tk_index) > 0) hL1Tks_d0_TkInSecondQuad->Fill( s->GetTrackD0(tk_index) );
      else if (L1Tks_POCAx->at(tk_index) < 0 && L1Tks_POCAy->at(tk_index) < 0) hL1Tks_d0_TkInThirdQuad ->Fill( s->GetTrackD0(tk_index) );
      else{}
    }

    }
  else
    {
      cout << "E R R O R ! Tracks::FillHistos(...). EXIT" << endl;
      exit(1);
    }
  
  return;
}

#endif // Tracks_cxx

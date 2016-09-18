#ifndef JetProbability_cxx
#define JetProbability_cxx

#include "JetProbability.h"
#include "../utilities/constants.h"

void JetProbability::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = (MaxEvents == -1) ? fChain->GetEntries() : std::min((int)fChain->GetEntries(), MaxEvents);

  std::cout << "Analyzing " << nentries << " events.\n";

  _eta_C             =   0.5; // defining eta points of interest, central and forward.
  _eta_F             =   1.5;
  // Book histograms here
  MakeHisto();

  Int_t nEvts=0;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    if(jentry%1000 == 0)
      std::cout << "Loop over entry " << jentry << "/" << nentries << ".\n";
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    //    cout<< jentry<<"============================================================"<<endl;
    //**************************************************************************************
    // Number of particles/tracks etc in an event.
    //**************************************************************************************
    Size_t nGenp     = GenP_Pt    ->size();
    Size_t nL1Tk     = L1Tks_Pt   ->size();
    Size_t nPixTk    = L1PixTks_Pt->size();
    Size_t nTP       = TP_Pt      ->size();
    Size_t nL1TkJet  = L1TkJet_Pt ->size();
   
    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	if(IsFinalBHadron(iGenp))
	  {

	    vector<Size_t> vDau;
	    GetBHadronChargedImmediateDaughters(iGenp, vDau);
	    //GetFinalDaughter(iGenp, vDau);
	    Double_t BPt  = GenP_Pt->at(iGenp);
	    Double_t BEta = GenP_Eta->at(iGenp);
	    Double_t BPhi = GenP_Phi->at(iGenp);
	    

	    for(Size_t iDau=0; iDau<vDau.size(); iDau++)
	      {
		Int_t    DauIndx   = vDau.at(iDau);
		Int_t    DauCharge = GenP_Charge->at(DauIndx);
		if(DauCharge != 0)// continue;
		  {
		    FillSIPHistogramsD0(DauIndx,BPhi,
		    			hgend0BHadDirect , hGenBHadPhiD0      , hGenBHaddeltaphiD0     , hGenBHadSIP,
		    			hGenMatchedTPB   ,hGenMatchedTPBPhiD0 , hGenMatchedTPBdeltaPhi , hGenMatchedTPBSIP,        
		    			hGenMatchedTkB   ,hGenMatchedTkBPhiD0 , hGenMatchedTkBdeltaPhi , hGenMatchedTkBSIP,
		    			hGenMatchedPixB  ,hGenMatchedPixBPhiD0, hGenMatchedPixBdeltaPhi, hGenMatchedPixBSIP        
		    			);
		    
		    //		    FillSIPHistogramsD0(DauIndx,BPhi);
		    
		    // FillSIPHistogramsDxy(DauIndx,BPhi,
		    // 			hgend0BHadDirect , hGenBHadPhiD0      , hGenBHaddeltaphiD0     , hGenBHadSIP,
		    // 			hGenMatchedTPB   ,hGenMatchedTPBPhiD0 , hGenMatchedTPBdeltaPhi , hGenMatchedTPBSIP,        
		    // 			hGenMatchedTkB   ,hGenMatchedTkBPhiD0 , hGenMatchedTkBdeltaPhi , hGenMatchedTkBSIP,
		    // 			hGenMatchedPixB  ,hGenMatchedPixBPhiD0, hGenMatchedPixBdeltaPhi, hGenMatchedPixBSIP        
		    // 			);
	
		    // TVector3 bHadvec;
		    // bHadvec.SetPtEtaPhi(BPt,BEta,BPhi);
		    // FillSIPHistogramsDotProd(DauIndx , bHadvec,
		    // 			     hgend0BHadDirect , hGenBHadSIP,
		    // 			     hGenMatchedTPB   , hGenMatchedTPBSIP,        
		    // 			     hGenMatchedTkB   , hGenMatchedTkBSIP,
		    // 			     hGenMatchedPixB  , hGenMatchedPixBSIP       
		    // 			     );

		  }//if charged
	      }// end loop on daughters
	    // cout<<<< iGenp<<"\t"<<Gen_ID<<endl;
	  }// if final state BHad

	if(IsLiteHadron(iGenp))
	  {
	   
	    Double_t LPt  = GenP_Pt->at(iGenp);
	    Double_t LEta = GenP_Eta->at(iGenp);
	    Double_t LPhi = GenP_Phi->at(iGenp);
	    if(RecursivelyLookForMotherId(iGenp,92, false))
	      {
		//cout<< iGenp<<"\t\t"<<GenP_PdgId->at(iGenp)<<endl;
		vector<Size_t> vPDau;
		GetFinalDaughter(iGenp, vPDau);
		if(vPDau.size()>1)
		  for(Size_t pDau=0; pDau<vPDau.size(); pDau++)
		    {
		      //		    cout<<GenP_Status->at(iGenp)<<endl;
		      Size_t   DIndx = vPDau.at(pDau);
		      if(!RecursivelyLookForMotherId(DIndx,5, false))
			if(!RecursivelyLookForMotherId(DIndx,4, false))
			  {
			    if(GenP_Charge->at(iGenp)!=0)
			      {
				// FillSIPHistogramsD0(DIndx,LPhi,
				// 		    hgend0UDS        , hGenLHadPhiD0       , hGenLHaddeltaphiD0     , hGenLHadSIP,
				// 		    hGenMatchedTPP   , hGenMatchedTPPPhiD0 , hGenMatchedTPPdeltaPhi , hGenMatchedTPPSIP,        
				// 		    hGenMatchedTkP   , hGenMatchedTkPPhiD0 , hGenMatchedTkPdeltaPhi , hGenMatchedTkPSIP,
				// 		    hGenMatchedPixP  , hGenMatchedPixPPhiD0, hGenMatchedPixPdeltaPhi, hGenMatchedPixPSIP        
				// 		    );
				
				// FillSIPHistogramsDxy(DIndx,LPhi,
				// 		    hgend0UDS        , hGenLHadPhiD0       , hGenLHaddeltaphiD0     , hGenLHadSIP,
				// 		    hGenMatchedTPP   , hGenMatchedTPPPhiD0 , hGenMatchedTPPdeltaPhi , hGenMatchedTPPSIP,        
				// 		    hGenMatchedTkP   , hGenMatchedTkPPhiD0 , hGenMatchedTkPdeltaPhi , hGenMatchedTkPSIP,
				// 		    hGenMatchedPixP  , hGenMatchedPixPPhiD0, hGenMatchedPixPdeltaPhi, hGenMatchedPixPSIP        
				// 		    );

				// TVector3 lHadvec;
				// lHadvec.SetPtEtaPhi(LPt,LEta,LPhi);
				// FillSIPHistogramsDotProd(DIndx , lHadvec,
				// 			 hgend0UDS        , hGenLHadSIP,
				// 			 hGenMatchedTPP   , hGenMatchedTPPSIP,        
				// 			 hGenMatchedTkP   , hGenMatchedTkPSIP,
				// 			 hGenMatchedPixP  , hGenMatchedPixPSIP       
				// 			 );
			      }// if charged
			  }// if not coming from b or c
		    }// loop on daughters
		
	      }// if coming from string
	  }//if is Lite Hadron
	
	
      }// end loop on genp
    
    
    // Update the progress bar
    //tools.ProgressBar(jentry, nentries, 100, 150);
  }// end loop on events
  // Keep this line here!
  outFile->cd();

  // Create and write canvases here


  
  // Uncomment this line to write also the histograms to the file
  outFile->Write();
}
//****************************************************************************************
void JetProbability::MakeHisto()
//****************************************************************************************
{
  hpixd0_C = new TH1D ("hpixd0_C"    ,"hpixd0_C"    , 100,-1,1);
  hpixd0_I = new TH1D ("hpixd0_I"    ,"hpixd0_I"    , 100,-1,1);
  hpixd0_F = new TH1D ("hpixd0_F"    ,"hpixd0_F"    , 100,-1,1);

  hpixeta_C = new TH1D ("hpixeta_C"    ,"hpixeta_C"    , 100,-3,3);
  hpixeta_I = new TH1D ("hpixeta_I"    ,"hpixeta_I"    , 100,-3,3);
  hpixeta_F = new TH1D ("hpixeta_F"    ,"hpixeta_F"    , 100,-3,3);
  
  hjet_pt    = new TH1D ("hjet_pt"    ,"hjet_pt"    , 100,0,200);
  hjet_eta   = new TH1D ("hjet_eta"   ,"hjet_eta"   , 100,-3,3);
  hjet_phi   = new TH1D ("hjet_phi"   ,"hjet_phi"   , 100,-4,4);
  hjet_e     = new TH1D ("hjet_e"     ,"hjet_e"     , 100,0,500);
  hjet_vtx   = new TH1D ("hjet_vtx"   ,"hjet_vtx"   , 100,-25,25);
  hjet_ntk   = new TH1D ("hjet_ntks"  ,"hjet_ntks"  , 100,0,15);
  hjet_multi = new TH1D ("hjet_multi" ,"hjet_multi" , 100,0,15);

  htk_pt     = new TH1D ("htk_pt"    ,"htk_pt"    , 100,0,100);
  htk_eta    = new TH1D ("htk_eta"   ,"htk_eta"   , 100,-3,3);
  htk_phi    = new TH1D ("htk_phi"   ,"htk_phi"   , 100,-4,4);
  htk_pocax  = new TH1D ("htk_pocax" ,"htk_pocax" , 100,-1,1);
  htk_pocay  = new TH1D ("htk_pocay" ,"htk_pocay" , 100,-1,1);
  htk_pocaz  = new TH1D ("htk_pocaz" ,"htk_pocaz" , 100,-25,25);
  htk_chi2   = new TH1D ("htk_chi2"  ,"htk_chi2"  , 100,0,100);
  htk_nstubs = new TH1D ("htk_nstubs","htk_nstubs", 100,0,10);
  htk_d0     = new TH1D ("htk_d0"    ,"htk_d0"    , 100,-1,1);


  hBGen_SignedD0     = new TH1D ("hBGen_SignedD0"   ,"hBGen_SignedD0"   , 100,-1,1);
  hBGen_DeltaPhi     = new TH1D ("hBGen_DeltaPhi"   ,"hBGen_DeltaPhi"   , 100,1,2);
  hBGen_SignedD0Sig  = new TH1D ("hBGen_SignedD0Sig","hBGen_SignedD0Sig", 100,-10,10);
  hsign              = new TH1D ("hsign"            ,"hsign"            , 100,-2 ,2);

  hBGen_SignedD0_C     = new TH1D ("hBGen_SignedD0_C"   ,"hBGen_SignedD0_C"   , 100,-1,1);
  hBGen_DeltaPhi_C     = new TH1D ("hBGen_DeltaPhi_C"   ,"hBGen_DeltaPhi_C"   , 100,1,2);
  hBGen_SignedD0Sig_C  = new TH1D ("hBGen_SignedD0Sig_C","hBGen_SignedD0Sig_C", 100,-10,10);
  hsign_C              = new TH1D ("hsign_C"            ,"hsign_C"            , 100,-2 ,2);

  hBGen_SignedD0_I     = new TH1D ("hBGen_SignedD0_I"   ,"hBGen_SignedD0_I"   , 100,-1,1);
  hBGen_DeltaPhi_I     = new TH1D ("hBGen_DeltaPhi_I"   ,"hBGen_DeltaPhi_I"   , 100,1,2);
  hBGen_SignedD0Sig_I  = new TH1D ("hBGen_SignedD0Sig_I","hBGen_SignedD0Sig_I", 100,-10,10);
  hsign_I              = new TH1D ("hsign_I"            ,"hsign_I"            , 100,-2 ,2);

  hBGen_SignedD0_F     = new TH1D ("hBGen_SignedD0_F"   ,"hBGen_SignedD0_F"   , 100,-1,1);
  hBGen_DeltaPhi_F     = new TH1D ("hBGen_DeltaPhi_F"   ,"hBGen_DeltaPhi_F"   , 100,1,2);
  hBGen_SignedD0Sig_F  = new TH1D ("hBGen_SignedD0Sig_F","hBGen_SignedD0Sig_F", 100,-10,10);
  hsign_F              = new TH1D ("hsign_F"            ,"hsign_F"            , 100,-2 ,2);

  htk_absd0H  = new TH1D ("htk_absd0H"   ,"htk_absd0H"   , 100,0,1);
  htk_absd0P  = new TH1D ("htk_absd0P"   ,"htk_absd0P"   , 100,0,1);


  hAlltk_d0   = new TH1D ("hAlltk_d0"  ,"hAlltk_d0"   , 100,-1,1);
  hAllpix_d0  = new TH1D ("hAllpix_d0" ,"hAllpix_d0"  , 100,-1,1);

  hgend0B     = new TH1D ("hgend0B"   ,"hgend0B"   , 100,0,1);
  hgend0P     = new TH1D ("hgend0P"   ,"hgend0P"   , 100,0,1);
  
  hgend0BHad     = new TH1D ("hgend0BHad"   ,"hgend0BHad"   , 100,0,1);
  hgend0PHad     = new TH1D ("hgend0PHad"   ,"hgend0PHad"   , 100,0,1);

  hBjetTkabsd0  = new TH1D ("hBjetTkabsd0"   ,"hBjetTkabsd0"   , 100,0,1);
  hPjetTkabsd0  = new TH1D ("hPjetTkabsd0"   ,"hPjetTkabsd0"   , 100,0,1);


  hAllTks_d0    = new TH1D ("hAllTks_d0"     ,"hAllTks_d0"     , 100,-1,1);
  hAllTks_dxy   = new TH1D ("hAllTks_dxy"    ,"hAllTks_dxy"    , 100,-1,1);
  hAllTks_phi   = new TH1D ("hAllTks_phi"    ,"hAllTks_phi"    , 100,-4,4);
  hAllTks_atan  = new TH1D ("hAllTks_atan"   ,"hAllTks_atan"   , 100,-4,4);
  
  hAllTks_Deltad0   = new TH1D ("hAllTks_Deltad0"   ,"hAllTks_Deltad0"   , 100,-10,10);
  hAllTks_Deltaphi  = new TH1D ("hAllTks_Deltaphi"  ,"hAllTks_Deltaphi"  , 100,-10,10);

  //--------------------------------------------------------------------------

  hgend0EMu         = new TH1D ("hgend0EMu"        ,"hgend0EMu"         , 100,0,0.5);
  hGenMatchedTkMuE  = new TH1D ("hGenMatchedTkMuE" , "hGenMatchedTkMuE" , 100,0,0.5);
  hGenMatchedPixMuE = new TH1D ("hGenMatchedPixMuE", "hGenMatchedPixMuE", 100,0,0.5);
  hGenMatchedTPMuE  = new TH1D ("hGenMatchedTPMuE" , "hGenMatchedTPMuE" , 100,0,0.5);


  hgend0UDS         = new TH1D ("hgend0UDS"        ,"hgend0UDS"         , 100,0,0.5);
  hGenMatchedTkP    = new TH1D ("hGenMatchedTkP"   , "hGenMatchedTkP"   , 100,0,0.5);
  hGenMatchedPixP   = new TH1D ("hGenMatchedPixP"  , "hGenMatchedPixP"  , 100,0,0.5);
  hGenMatchedTPP    = new TH1D ("hGenMatchedTPP"   , "hGenMatchedTPP"   , 100,0,0.5);


  hgend0BHadDirect  = new TH1D ("hgend0BHadDirect" , "hgend0BHadDirect" , 100,0,0.5);
  hGenMatchedTkB    = new TH1D ("hGenMatchedTkB"   , "hGenMatchedTkB"   , 100,0,0.5);
  hGenMatchedPixB   = new TH1D ("hGenMatchedPixB"  , "hGenMatchedPixB"  , 100,0,0.5); 
  hGenMatchedTPB    = new TH1D ("hGenMatchedTPB"   , "hGenMatchedTPb"   , 100,0,0.5);

  //**********************Histo for B SIP ETC
 
  hGenBHadPhiD0      = new TH1D ("hGenBHadPhiD0"     , "hGenBHadPhiD0"     , 100 , -4  ,4);
  hGenBHaddeltaphiD0 = new TH1D ("hGenBHaddeltaphiD0", "hGenBHaddeltaphiD0", 100 , 0   ,4);
  hGenBHadSIP        = new TH1D ("hGenBHadSIP"       , "hGenBHadSIP"       , 100 , -0.5,0.5);
  
  hGenMatchedTPBPhiD0       = new TH1D ("hGenMatchedTPBPhiD0"      , "hGenMatchedTPBPhiD0"      , 100 , -4,4);
  hGenMatchedTPBdeltaPhi    = new TH1D ("hGenMatchedTPBdeltaPhi"   , "hGenMatchedTPBdeltaPhi"   , 100 , 0,4); 
  hGenMatchedTPBSIP         = new TH1D ("hGenMatchedTPBSIP"        , "hGenMatchedTPBSIP"        , 100 , -0.5,0.5);      

  hGenMatchedTkBPhiD0       = new TH1D ("hGenMatchedTkBPhiD0"      , "hGenMatchedTkBPhiD0"      , 100 , -4,4);
  hGenMatchedTkBdeltaPhi    = new TH1D ("hGenMatchedTkBdeltaPhi"   , "hGenMatchedTkBdeltaPhi"   , 100 , 0,4); 
  hGenMatchedTkBSIP         = new TH1D ("hGenMatchedTkBSIP"        , "hGenMatchedTkBSIP"        , 100 , -0.5,0.5);      

  hGenMatchedPixBPhiD0       = new TH1D ("hGenMatchedPixBPhiD0"      , "hGenMatchedPixBPhiD0"      , 100 , -4,4);
  hGenMatchedPixBdeltaPhi    = new TH1D ("hGenMatchedPixBdeltaPhi"   , "hGenMatchedPixBdeltaPhi"   , 100 , 0,4); 
  hGenMatchedPixBSIP         = new TH1D ("hGenMatchedPixBSIP"        , "hGenMatchedPixBSIP"        , 100 , -0.5,0.5);      

 
  //**********************Histo for B SIP ETC

  //**********************Histo for Prompt SIP ETC
  hGenMatchedTPPPhiD0       = new TH1D ("hGenMatchedTPPPhiD0"      , "hGenMatchedTPPPhiD0"      , 100 , -4,4);
  hGenMatchedTPPdeltaPhi    = new TH1D ("hGenMatchedTPPdeltaPhi"   , "hGenMatchedTPPdeltaPhi"   , 100 , 0,4); 
  hGenMatchedTPPSIP         = new TH1D ("hGenMatchedTPPSIP"        , "hGenMatchedTPPSIP"        , 100 , -0.5,0.5);      

  hGenMatchedTkPPhiD0       = new TH1D ("hGenMatchedTkPPhiD0"      , "hGenMatchedTkPPhiD0"      , 100 , -4,4);
  hGenMatchedTkPdeltaPhi    = new TH1D ("hGenMatchedTkPdeltaPhi"   , "hGenMatchedTkPdeltaPhi"   , 100 , 0,4); 
  hGenMatchedTkPSIP         = new TH1D ("hGenMatchedTkPSIP"        , "hGenMatchedTkPSIP"        , 100 , -0.5,0.5);      

  hGenMatchedPixPPhiD0       = new TH1D ("hGenMatchedPixPPhiD0"      , "hGenMatchedPixPPhiD0"      , 100 , -4,4);
  hGenMatchedPixPdeltaPhi    = new TH1D ("hGenMatchedPixPdeltaPhi"   , "hGenMatchedPixPdeltaPhi"   , 100 , 0,4); 
  hGenMatchedPixPSIP         = new TH1D ("hGenMatchedPixPSIP"        , "hGenMatchedPixPSIP"        , 100 , -0.5,0.5);      

  hGenLHadPhiD0      = new TH1D ("hGenLHadPhiD0", "hGenLHadPhiD0", 100 , -4,4);
  hGenLHaddeltaphiD0 = new TH1D ("hGenLHaddeltaphiD0", "hGenLHaddeltaphiD0", 100 , 0,4);
  hGenLHadSIP        = new TH1D ("hGenLHadSIP", "hGenLHadSIP", 100 , -0.5,0.5);
  //**********************Histo for Prompt SIP ETC

  
  hGenNoPUSIP = new TH1D ("hGenNoPUSIP", "hGenNoPUSIP", 100 , -1,1);
  
  htkDxy = new TH1D ("htkDxy","htkDxy",100,-2,2);
  //__________________________________________________________________
  hgenphiHadDau   = new TH1D ("hgenphiHadDau"  , "hgenphiHadDau", 100 , -4,4);
  //hgenphiD0HadDau = new TH1D ("hgenphiD0HadDau", "hgenphiD0HadDau", 100 , -6,6);
  //hgendeltaphiD0HadDau = new TH1D ("hgendeltaphiD0HadDau", "hgendeltaphiD0HadDau", 100 , -6,6);

  hgend0BHadAll  = new TH1D ("hgend0BHadAll" , "hgend0BHadAll" , 100, 0,1);
  hgend0BHadNoC  = new TH1D ("hgend0BHadNoC" , "hgend0BHadNoC" , 100, 0,1);
  hgend0BHadNoCS = new TH1D ("hgend0BHadNoCS", "hgend0BHadNoCS", 100, 0,1);


  htempgenD0  = new TH1D ("htempgenD0", "htempgenD0", 100, -1,1);
  htemptpD0   = new TH1D ("htemptpD0" , "htemptpD0" , 100, -1,1);



  temp1 = new TH1D("temp1","temp1",100,0,0.5);
  temp2 = new TH1D("temp2","temp2",100,0,0.5);
  temp3 = new TH1D("temp3","temp3",100,0,0.5);




  //****************************************************************************************************
  // Histograms for slices of eta and Pt Pix Tracks
  //****************************************************************************************************

  PixDeltaPhi_C_0203 = new TH1D("PixDeltaPhi_C_0203","PixDeltaPhi_C_0203",100,-4,4);
  PixDeltaPhi_C_0304 = new TH1D("PixDeltaPhi_C_0304","PixDeltaPhi_C_0304",100,-4,4);
  PixDeltaPhi_C_0405 = new TH1D("PixDeltaPhi_C_0405","PixDeltaPhi_C_0405",100,-4,4);
  PixDeltaPhi_C_0506 = new TH1D("PixDeltaPhi_C_0506","PixDeltaPhi_C_0506",100,-4,4);
  PixDeltaPhi_C_0607 = new TH1D("PixDeltaPhi_C_0607","PixDeltaPhi_C_0607",100,-4,4);
  PixDeltaPhi_C_0708 = new TH1D("PixDeltaPhi_C_0708","PixDeltaPhi_C_0708",100,-4,4);
  PixDeltaPhi_C_0809 = new TH1D("PixDeltaPhi_C_0809","PixDeltaPhi_C_0809",100,-4,4);
  PixDeltaPhi_C_0910 = new TH1D("PixDeltaPhi_C_0910","PixDeltaPhi_C_0910",100,-4,4);
  PixDeltaPhi_C_1015 = new TH1D("PixDeltaPhi_C_1015","PixDeltaPhi_C_1015",100,-4,4);
  PixDeltaPhi_C_1520 = new TH1D("PixDeltaPhi_C_1520","PixDeltaPhi_C_1520",100,-4,4);
  PixDeltaPhi_C_2000 = new TH1D("PixDeltaPhi_C_2000","PixDeltaPhi_C_2000",100,-4,4);

  PixSIP_C_0203      = new TH1D ("PixSIP_C_0203", "PixSIP_C_0203", 100 , -0.5,0.5);
  PixSIP_C_0304      = new TH1D ("PixSIP_C_0304", "PixSIP_C_0304", 100 , -0.5,0.5);
  PixSIP_C_0405      = new TH1D ("PixSIP_C_0405", "PixSIP_C_0405", 100 , -0.5,0.5);
  PixSIP_C_0506      = new TH1D ("PixSIP_C_0506", "PixSIP_C_0506", 100 , -0.5,0.5);
  PixSIP_C_0607      = new TH1D ("PixSIP_C_0607", "PixSIP_C_0607", 100 , -0.5,0.5);
  PixSIP_C_0708      = new TH1D ("PixSIP_C_0708", "PixSIP_C_0708", 100 , -0.5,0.5);
  PixSIP_C_0809      = new TH1D ("PixSIP_C_0809", "PixSIP_C_0809", 100 , -0.5,0.5);
  PixSIP_C_0910      = new TH1D ("PixSIP_C_0910", "PixSIP_C_0910", 100 , -0.5,0.5);
  PixSIP_C_1015      = new TH1D ("PixSIP_C_1015", "PixSIP_C_1015", 100 , -0.5,0.5);
  PixSIP_C_1520      = new TH1D ("PixSIP_C_1520", "PixSIP_C_1520", 100 , -0.5,0.5);
  PixSIP_C_2000      = new TH1D ("PixSIP_C_2000", "PixSIP_C_2000", 100 , -0.5,0.5);


  PixDeltaPhi_I_0203 = new TH1D("PixDeltaPhi_I_0203","PixDeltaPhi_I_0203",100,-4,4);
  PixDeltaPhi_I_0304 = new TH1D("PixDeltaPhi_I_0304","PixDeltaPhi_I_0304",100,-4,4);
  PixDeltaPhi_I_0405 = new TH1D("PixDeltaPhi_I_0405","PixDeltaPhi_I_0405",100,-4,4);
  PixDeltaPhi_I_0506 = new TH1D("PixDeltaPhi_I_0506","PixDeltaPhi_I_0506",100,-4,4);
  PixDeltaPhi_I_0607 = new TH1D("PixDeltaPhi_I_0607","PixDeltaPhi_I_0607",100,-4,4);
  PixDeltaPhi_I_0708 = new TH1D("PixDeltaPhi_I_0708","PixDeltaPhi_I_0708",100,-4,4);
  PixDeltaPhi_I_0809 = new TH1D("PixDeltaPhi_I_0809","PixDeltaPhi_I_0809",100,-4,4);
  PixDeltaPhi_I_0910 = new TH1D("PixDeltaPhi_I_0910","PixDeltaPhi_I_0910",100,-4,4);
  PixDeltaPhi_I_1015 = new TH1D("PixDeltaPhi_I_1015","PixDeltaPhi_I_1015",100,-4,4);
  PixDeltaPhi_I_1520 = new TH1D("PixDeltaPhi_I_1520","PixDeltaPhi_I_1520",100,-4,4);
  PixDeltaPhi_I_2000 = new TH1D("PixDeltaPhi_I_2000","PixDeltaPhi_I_2000",100,-4,4);

  PixSIP_I_0203      = new TH1D ("PixSIP_I_0203", "PixSIP_I_0203", 100 , -0.5,0.5);
  PixSIP_I_0304      = new TH1D ("PixSIP_I_0304", "PixSIP_I_0304", 100 , -0.5,0.5);
  PixSIP_I_0405      = new TH1D ("PixSIP_I_0405", "PixSIP_I_0405", 100 , -0.5,0.5);
  PixSIP_I_0506      = new TH1D ("PixSIP_I_0506", "PixSIP_I_0506", 100 , -0.5,0.5);
  PixSIP_I_0607      = new TH1D ("PixSIP_I_0607", "PixSIP_I_0607", 100 , -0.5,0.5);
  PixSIP_I_0708      = new TH1D ("PixSIP_I_0708", "PixSIP_I_0708", 100 , -0.5,0.5);
  PixSIP_I_0809      = new TH1D ("PixSIP_I_0809", "PixSIP_I_0809", 100 , -0.5,0.5);
  PixSIP_I_0910      = new TH1D ("PixSIP_I_0910", "PixSIP_I_0910", 100 , -0.5,0.5);
  PixSIP_I_1015      = new TH1D ("PixSIP_I_1015", "PixSIP_I_1015", 100 , -0.5,0.5);
  PixSIP_I_1520      = new TH1D ("PixSIP_I_1520", "PixSIP_I_1520", 100 , -0.5,0.5);
  PixSIP_I_2000      = new TH1D ("PixSIP_I_2000", "PixSIP_I_2000", 100 , -0.5,0.5);


  PixDeltaPhi_F_0203 = new TH1D("PixDeltaPhi_F_0203","PixDeltaPhi_F_0203",100,-4,4);
  PixDeltaPhi_F_0304 = new TH1D("PixDeltaPhi_F_0304","PixDeltaPhi_F_0304",100,-4,4);
  PixDeltaPhi_F_0405 = new TH1D("PixDeltaPhi_F_0405","PixDeltaPhi_F_0405",100,-4,4);
  PixDeltaPhi_F_0506 = new TH1D("PixDeltaPhi_F_0506","PixDeltaPhi_F_0506",100,-4,4);
  PixDeltaPhi_F_0607 = new TH1D("PixDeltaPhi_F_0607","PixDeltaPhi_F_0607",100,-4,4);
  PixDeltaPhi_F_0708 = new TH1D("PixDeltaPhi_F_0708","PixDeltaPhi_F_0708",100,-4,4);
  PixDeltaPhi_F_0809 = new TH1D("PixDeltaPhi_F_0809","PixDeltaPhi_F_0809",100,-4,4);
  PixDeltaPhi_F_0910 = new TH1D("PixDeltaPhi_F_0910","PixDeltaPhi_F_0910",100,-4,4);
  PixDeltaPhi_F_1015 = new TH1D("PixDeltaPhi_F_1015","PixDeltaPhi_F_1015",100,-4,4);
  PixDeltaPhi_F_1520 = new TH1D("PixDeltaPhi_F_1520","PixDeltaPhi_F_1520",100,-4,4);
  PixDeltaPhi_F_2000 = new TH1D("PixDeltaPhi_F_2000","PixDeltaPhi_F_2000",100,-4,4);

  PixSIP_F_0203      = new TH1D ("PixSIP_F_0203", "PixSIP_F_0203", 100 , -0.5,0.5);
  PixSIP_F_0304      = new TH1D ("PixSIP_F_0304", "PixSIP_F_0304", 100 , -0.5,0.5);
  PixSIP_F_0405      = new TH1D ("PixSIP_F_0405", "PixSIP_F_0405", 100 , -0.5,0.5);
  PixSIP_F_0506      = new TH1D ("PixSIP_F_0506", "PixSIP_F_0506", 100 , -0.5,0.5);
  PixSIP_F_0607      = new TH1D ("PixSIP_F_0607", "PixSIP_F_0607", 100 , -0.5,0.5);
  PixSIP_F_0708      = new TH1D ("PixSIP_F_0708", "PixSIP_F_0708", 100 , -0.5,0.5);
  PixSIP_F_0809      = new TH1D ("PixSIP_F_0809", "PixSIP_F_0809", 100 , -0.5,0.5);
  PixSIP_F_0910      = new TH1D ("PixSIP_F_0910", "PixSIP_F_0910", 100 , -0.5,0.5);
  PixSIP_F_1015      = new TH1D ("PixSIP_F_1015", "PixSIP_F_1015", 100 , -0.5,0.5);
  PixSIP_F_1520      = new TH1D ("PixSIP_F_1520", "PixSIP_F_1520", 100 , -0.5,0.5);
  PixSIP_F_2000      = new TH1D ("PixSIP_F_2000", "PixSIP_F_2000", 100 , -0.5,0.5);

  PixDeltaPhi_C = new TH1D("PixDeltaPhi_C","PixDeltaPhi_C",100  ,-4   ,4  );
  PixSIP_C      = new TH1D("PixSIP_C"     , "PixSIP_C"    , 100 , -0.5,0.5);
  PixDeltaPhi_I = new TH1D("PixDeltaPhi_I","PixDeltaPhi_I",100  ,-4   ,4  );
  PixSIP_I      = new TH1D("PixSIP_I"     , "PixSIP_I"    , 100 , -0.5,0.5);
  PixDeltaPhi_F = new TH1D("PixDeltaPhi_F","PixDeltaPhi_F",100  ,-4   ,4  );
  PixSIP_F      = new TH1D("PixSIP_F"     , "PixSIP_F"    , 100 , -0.5,0.5);




   PixDeltaPhi_T_c2 = new TH1D("PixDeltaPhi_T_c2","PixDeltaPhi_T_c2",100  ,-4   ,4  );
   PixSIP_T_c2      = new TH1D("PixSIP_T_c2"     , "PixSIP_T_c2"    , 100 , -0.5,0.5);

   PixDeltaPhi_M_c2 = new TH1D("PixDeltaPhi_M_c2","PixDeltaPhi_M_c2",100  ,-4   ,4  );
   PixSIP_M_c2      = new TH1D("PixSIP_M_c2"     , "PixSIP_M_c2"    , 100 , -0.5,0.5);

   PixDeltaPhi_L_c2 = new TH1D("PixDeltaPhi_L_c2","PixDeltaPhi_L_c2",100  ,-4   ,4  );
   PixSIP_L_c2      = new TH1D("PixSIP_L_c2"     , "PixSIP_L_c2"    , 100 , -0.5,0.5);

   PixDeltaPhi_VL_c2 = new TH1D("PixDeltaPhi_VL_c2","PixDeltaPhi_VL_c2",100  ,-4   ,4  );
   PixSIP_VL_c2      = new TH1D("PixSIP_VL_c2"     , "PixSIP_VL_c2"    , 100 , -0.5,0.5);


}

//****************************************************************************************
Int_t JetProbability::GetPixTkIndex(Int_t TTTrackIndex)
//****************************************************************************************
{
  Size_t nPix = L1PixTks_Pt->size();
  Int_t PixIndx=-1;
  for(Size_t iPix=0; iPix<nPix; iPix++)
    {
      Int_t tkIndex = L1PixTks_TTTrackIndex->at(iPix);
      if(TTTrackIndex == tkIndex)
	{
	  PixIndx = iPix;
	  break;
	}
    }
  return PixIndx;
}



//****************************************************************************
Bool_t JetProbability::IsWithinEtaRegion(string etaRegion, 
					 Double_t eta)
//****************************************************************************
{

  bool bWithinEtaRegion = false;
  if ( etaRegion.compare("Central") == 0 )           bWithinEtaRegion = (fabs(eta) < _eta_C);
  else if ( etaRegion.compare("Intermediate") == 0 ) bWithinEtaRegion = (fabs(eta) < _eta_F && fabs(eta) >= _eta_C);
  else if ( etaRegion.compare("Forward") == 0 )      bWithinEtaRegion = (fabs(eta) >= _eta_F);
  else{
    cout << "E R R O R ! Tracking::IsWithinEtaRegion(...) - Invalid eta region type \"" << etaRegion << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinEtaRegion;  
}


//****************************************************************************
Bool_t JetProbability::IsWithinPtSlice(string pTSlice, 
				       Double_t pt)
//****************************************************************************
{

  bool bWithinPtSlice = false;
  if ( pTSlice.compare("SL0203") == 0 )         bWithinPtSlice = (pt > 2  && pt <= 3  );
  else if ( pTSlice.compare("SL0304") == 0 )    bWithinPtSlice = (pt > 3  && pt <= 4  );
  else if ( pTSlice.compare("SL0405") == 0 )    bWithinPtSlice = (pt > 4  && pt <= 5  );
  else if ( pTSlice.compare("SL0506") == 0 )    bWithinPtSlice = (pt > 5  && pt <= 6  );
  else if ( pTSlice.compare("SL0607") == 0 )    bWithinPtSlice = (pt > 6  && pt <= 7  );
  else if ( pTSlice.compare("SL0708") == 0 )    bWithinPtSlice = (pt > 7  && pt <= 8  );
  else if ( pTSlice.compare("SL0809") == 0 )    bWithinPtSlice = (pt > 8  && pt <= 9  );
  else if ( pTSlice.compare("SL0910") == 0 )    bWithinPtSlice = (pt > 9  && pt <= 10 );
  else if ( pTSlice.compare("SL1015") == 0 )    bWithinPtSlice = (pt > 10 && pt <= 15 );
  else if ( pTSlice.compare("SL1520") == 0 )    bWithinPtSlice = (pt > 15 && pt <= 20 );
  else if ( pTSlice.compare("SL2000") == 0 )    bWithinPtSlice = (pt > 20);


  else{
    cout << "E R R O R ! Tracking::IsWithinPtSlice(...) - Invalid Pt Slice type \"" << pTSlice << "\". EXIT" << endl;
    exit(1);
  }
  
  return bWithinPtSlice;  
}

//****************************************************************************
Bool_t JetProbability::IsTrackChi2Good(string chi2Sel, 
				       Double_t chi2)
//****************************************************************************
{

  bool bChi2Lessthanselection = false;
  if ( chi2Sel.compare("Tight") == 0 )       bChi2Lessthanselection = (chi2 < 5);
  else if ( chi2Sel.compare("Medium") == 0 ) bChi2Lessthanselection = (chi2 < 10);
  else if ( chi2Sel.compare("Loose") == 0 )  bChi2Lessthanselection = (chi2 < 20);
  else if ( chi2Sel.compare("VLoose") == 0 )  bChi2Lessthanselection = (chi2 < 50);
  else{
    cout << "E R R O R ! Tracking::IsTrackChi2Good(...) - Invalid eta region type \"" << chi2Sel << "\". EXIT" << endl;
    exit(1);
  }
  
  return bChi2Lessthanselection;  
}


// //****************************************************************************
// Bool_t JetProbability::IsJetInBDirection(Int_t JetIndex)
// //****************************************************************************
// {
//   Double_t jet_eta = L1TkJet_Eta    ->at(JetIndex);
//   Double_t jet_phi = L1TkJet_Phi    ->at(JetIndex);
//   Bool_t bIsinBDir = false;
//   //cout<<endl;
//   for(Size_t iGenp=0; iGenp<GenP_Pt->size(); iGenp++)
//     {
//       if(abs(GenP_PdgId->at(iGenp)) != 5) continue;              // selecting b quark from event
//       vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
//       if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that b is used only once
//       Double_t gen_eta     = GenP_Eta->at(iGenp);
//       Double_t gen_phi     = GenP_Phi->at(iGenp);
//       Double_t dr          = auxTools_.DeltaR(jet_eta,jet_phi, 
// 					      gen_eta,gen_phi);
//       //  cout<<JetIndex<<"\tdr\t"<<dr<<endl;
//       if(dr<1.0)
// 	bIsinBDir = true;
//     }
//   return bIsinBDir;
// }

// //****************************************************************************
// Bool_t JetProbability::IsJetInCDirection(Int_t JetIndex)
// //****************************************************************************
// {
//   Double_t jet_eta = L1TkJet_Eta    ->at(JetIndex);
//   Double_t jet_phi = L1TkJet_Phi    ->at(JetIndex);
//   Bool_t bIsinBDir = false;
//   for(Size_t iGenp=0; iGenp<GenP_Pt->size(); iGenp++)
//     {
//       if(abs(GenP_PdgId->at(iGenp)) != 4) continue;              // selecting c quark from event
//       vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
//       if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that c is used only once
//       Double_t gen_eta     = GenP_Eta->at(iGenp);
//       Double_t gen_phi     = GenP_Phi->at(iGenp);
//       Double_t dr          = auxTools_.DeltaR(jet_eta,jet_phi, 
// 					      gen_eta,gen_phi);
//       if(dr<1.0)
// 	bIsinBDir = true;
//     }
//   return bIsinBDir;
// }

// //****************************************************************************
// Bool_t JetProbability::IsJetInUDirection(Int_t JetIndex)
// //****************************************************************************
// {
//   Double_t jet_eta = L1TkJet_Eta    ->at(JetIndex);
//   Double_t jet_phi = L1TkJet_Phi    ->at(JetIndex);
//   Bool_t bIsinBDir = false;
//   //cout<<endl;
//   for(Size_t iGenp=0; iGenp<GenP_Pt->size(); iGenp++)
//     {
//       if(abs(GenP_PdgId->at(iGenp)) != 2) continue;              // selecting b quark from event
//       vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
//       if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that b is used only once
//       Double_t gen_eta     = GenP_Eta->at(iGenp);
//       Double_t gen_phi     = GenP_Phi->at(iGenp);
//       Double_t dr          = auxTools_.DeltaR(jet_eta,jet_phi, 
// 					      gen_eta,gen_phi);
//       // cout<<JetIndex<<"\tdr\t"<<dr<<endl;
//       if(dr<0.5)
// 	bIsinBDir = true;
//     }
//   return bIsinBDir;
// }

// //****************************************************************************
// Bool_t JetProbability::IsJetInDDirection(Int_t JetIndex)
// //****************************************************************************
// {
//   Double_t jet_eta = L1TkJet_Eta    ->at(JetIndex);
//   Double_t jet_phi = L1TkJet_Phi    ->at(JetIndex);
//   Bool_t bIsinBDir = false;
//   //cout<<endl;
//   for(Size_t iGenp=0; iGenp<GenP_Pt->size(); iGenp++)
//     {
//       if(abs(GenP_PdgId->at(iGenp)) != 1) continue;              // selecting b quark from event
//       vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
//       if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that b is used only once
//       Double_t gen_eta     = GenP_Eta->at(iGenp);
//       Double_t gen_phi     = GenP_Phi->at(iGenp);
//       Double_t dr          = auxTools_.DeltaR(jet_eta,jet_phi, 
// 					      gen_eta,gen_phi);
//       //  cout<<JetIndex<<"\tdr\t"<<dr<<endl;
//       if(dr<0.5)
// 	bIsinBDir = true;
//     }
//   return bIsinBDir;
// }

// //****************************************************************************
// Bool_t JetProbability::IsJetHeavyFlavour(Int_t JetIndex)
// //****************************************************************************
// {
//   Double_t jet_eta = L1TkJet_Eta    ->at(JetIndex);
//   Double_t jet_phi = L1TkJet_Phi    ->at(JetIndex);
//   Bool_t bIsHF = false;
//   //cout<<endl;
//   for(Size_t iGenp=0; iGenp<GenP_Pt->size(); iGenp++)
//     {
//       if(abs(GenP_PdgId->at(iGenp)) == 5 || abs(GenP_PdgId->at(iGenp)) == 4)// continue;              // selecting b quark from event
//       {
// 	vector<unsigned short> DauVec = GenP_Daughters->at(iGenp);                                    // getting vector of daughters
// 	if(!IsFinalGenp(iGenp, DauVec)) continue;                                                     // making sure that b is used only once
// 	//cout<<"\tThe quark is\t"<<GenP_PdgId->at(iGenp)<<endl;
// 	Double_t gen_eta     = GenP_Eta->at(iGenp);
// 	Double_t gen_phi     = GenP_Phi->at(iGenp);
// 	Double_t dr          = auxTools_.DeltaR(jet_eta,jet_phi, 
// 					      gen_eta,gen_phi);
// 	//cout<<JetIndex<<"\tdr\t"<<dr<<endl;
// 	if(dr<1.0)
// 	  bIsHF = true;
//       }
//     }
//   return bIsHF;
// }

//****************************************************************************
Bool_t JetProbability::IsFinalGenp (Size_t MotIndx, vector<unsigned short>& Daug)
//****************************************************************************
{
  Int_t MotId = GenP_PdgId->at(MotIndx);
  
  if(Daug.size()==0) return true;
  
  for(Size_t iDau=0; iDau<Daug.size(); iDau++)
    {
      Int_t DauIndx = Daug.at(iDau);
      Int_t DauId   = GenP_PdgId->at(DauIndx);
      
      if(MotId == DauId)
	return false;
    }// end loop on daughters
  return true;
}

// //****************************************************************************
// Int_t JetProbability::FindTrackMatchingQuark(Size_t Index)
// //****************************************************************************
// {
//   Bool_t bIdBMatched = false; // if this module not required, make this variable true.
//   Bool_t bIdCMatched = false;
//   Bool_t bIdSMatched = false;
//   Bool_t bIdUMatched = false;
//   Bool_t bIdDMatched = false;
//   Int_t  quarkID = 0;

//   for(Size_t iTP=0; iTP<TP_Pt->size(); iTP++)
//     {
//       Size_t tp_tkindex = TP_TTTrackIndex->at(iTP);
//       if(tp_tkindex == -1) continue;
//       Double_t tp_eta = TP_Eta->at(iTP);
//       Double_t tp_phi = TP_Phi->at(iTP);
//       Int_t tp_genId = TP_PdgId->at(iTP);
//       //      cout<<"\tTP Track ID\t"<<tp_tkindex<<"\tTrack ID \t"<<Index<<endl;
//       if(tp_tkindex == Index)
// 	{
// 	  for(Size_t iGen =0; iGen<GenP_Pt->size(); iGen++)
// 	    {
// 	      if(GenP_Status->at(iGen)!=1)continue;
// 	      Int_t gen_id = GenP_PdgId->at(iGen);
// 	      if(gen_id != tp_genId)// selecting only if ID match
// 		continue;
// 	      Double_t gen_eta = GenP_Eta->at(iGen);
// 	      Double_t gen_phi = GenP_Phi->at(iGen);
// 	      Double_t temp = auxTools_.DeltaR(tp_eta, tp_phi,
// 					       gen_eta, gen_phi);
// 	      if(temp > 0.0001) // matching TP and Genp on base of dr.... have to improve this part.
// 		continue;
	      
// 	      if(RecursivelyLookForMotherId(iGen,5, false)) 
// 		{
// 		  if(!RecursivelyLookForMotherId(iGen,4, false))
// 		    quarkID =5;
// 		}
// 	      else if(RecursivelyLookForMotherId(iGen,4, false))
// 		{
// 		  quarkID =4;
// 		}
// 	      else
// 		{
// 		  quarkID =123;
// 		}
// 	      // if(IsBHadronProduct(iGen)) quarkID=5;
// 	      // else if(IsCHadronProduct(iGen)) quarkID=4;
// 	      // else if(IsSHadronProduct(iGen)) quarkID=3;
// 	      // else if(IsUHadronProduct(iGen)) quarkID=2;
// 	      // else if(IsHadronProduct(iGen)) quarkID=1;
// 	      //cout<<"\tquark\t"<<quarkID<<endl;

// 	      //bIdMatched =  IsBHadronProduct(iGen);
// 	    }// end loop on Gen
// 	}// end if track matched with index
//     }// end loop on TP

//   return quarkID;
// }

// //**************************************************
// Bool_t JetProbability::IsBHadronProduct(Size_t Indx)
// //**************************************************
// {
//   Size_t nMom = GenP_Mothers->at(Indx).size();
//   if (nMom == 0)return false;
//   Int_t MyId=0;
//   for (Size_t iMom = 0; iMom < nMom; iMom++)
//     {
//       Int_t MomID   = fabs(GenP_PdgId->at(GenP_Mothers->at(Indx).at(iMom)));
//       Int_t Modulus = MomID%10000;
//       if(Modulus>1000)
// 	{
// 	  MyId = Modulus/1000;
// 	}
//       else if(Modulus<1000)
// 	{
// 	  MyId = Modulus/100;
// 	}
//       if(MyId == 5)
// 	return true;
//       if (IsBHadronProduct(GenP_Mothers->at(Indx).at(iMom)))
// 	return true;
      
//     }
//   return false;
// }

// //**************************************************
// Bool_t JetProbability::IsCHadronProduct(Size_t Indx)
// //**************************************************
// {
//   Size_t nMom = GenP_Mothers->at(Indx).size();
//   if (nMom == 0)return false;
//   Int_t MyId=0;
//   for (Size_t iMom = 0; iMom < nMom; iMom++)
//     {
//       Int_t MomID   = fabs(GenP_PdgId->at(GenP_Mothers->at(Indx).at(iMom)));
//       Int_t Modulus = MomID%10000;
//       if(Modulus>1000)
// 	{
// 	  MyId = Modulus/1000;
// 	}
//       else if(Modulus<1000)
// 	{
// 	  MyId = Modulus/100;
// 	}
//       if(MyId == 4)
// 	return true;
//       if (IsCHadronProduct(GenP_Mothers->at(Indx).at(iMom)))
// 	return true;
      
//     }
//   return false;
// }

// //**************************************************
// Bool_t JetProbability::IsSHadronProduct(Size_t Indx)
// //**************************************************
// {
//   Size_t nMom = GenP_Mothers->at(Indx).size();
//   if (nMom == 0)return false;
//   Int_t MyId=0;
//   for (Size_t iMom = 0; iMom < nMom; iMom++)
//     {
//       Int_t MomID   = fabs(GenP_PdgId->at(GenP_Mothers->at(Indx).at(iMom)));
//       Int_t Modulus = MomID%10000;
//       if(Modulus>1000)
// 	{
// 	  MyId = Modulus/1000;
// 	}
//       else if(Modulus<1000)
// 	{
// 	  MyId = Modulus/100;
// 	}
//       if(MyId == 3)
// 	return true;
//       if (IsSHadronProduct(GenP_Mothers->at(Indx).at(iMom)))
// 	return true;
      
//     }
//   return false;
// }

// //**************************************************
// Bool_t JetProbability::IsUHadronProduct(Size_t Indx)
// //**************************************************
// {
//   Size_t nMom = GenP_Mothers->at(Indx).size();
//   if (nMom == 0)return false;
//   Int_t MyId=0;
//   for (Size_t iMom = 0; iMom < nMom; iMom++)
//     {
//       Int_t MomID   = fabs(GenP_PdgId->at(GenP_Mothers->at(Indx).at(iMom)));
//       Int_t Modulus = MomID%10000;
//       if(Modulus>1000)
// 	{
// 	  MyId = Modulus/1000;
// 	}
//       else if(Modulus<1000)
// 	{
// 	  MyId = Modulus/100;
// 	}
//       if(MyId == 2)
// 	return true;
//       if (IsUHadronProduct(GenP_Mothers->at(Indx).at(iMom)))
// 	return true;
      
//     }
//   return false;
// }

// //**************************************************
// Bool_t JetProbability::IsDHadronProduct(Size_t Indx)
// //**************************************************
// {
//   Size_t nMom = GenP_Mothers->at(Indx).size();
//   if (nMom == 0)return false;
//   Int_t MyId=0;
//   for (Size_t iMom = 0; iMom < nMom; iMom++)
//     {
//       Int_t MomID   = fabs(GenP_PdgId->at(GenP_Mothers->at(Indx).at(iMom)));
//       Int_t Modulus = MomID%10000;
//       if(Modulus>1000)
// 	{
// 	  MyId = Modulus/1000;
// 	}
//       else if(Modulus<1000)
// 	{
// 	  MyId = Modulus/100;
// 	}
//       if(MyId == 1)
// 	return true;
//       if (IsDHadronProduct(GenP_Mothers->at(Indx).at(iMom)))
// 	return true;
      
//     }
//   return false;
// }

//**************************************************
void JetProbability::myPrintGenp(Size_t Indx, bool bPrintHeaders)
//**************************************************
{
  unsigned int moth1 = 0;
  unsigned int moth2 = 0;
  unsigned int daug1 = 0;
  unsigned int daug2 = 0;
  unsigned int NDaug;
  static Int_t evOld = 0;

  ofstream myfile;
  myfile.open ("PrintGenp.txt", ios::app);

  if ( (evOld != EvtNumber) || (Indx == 0) ) {
    myfile << std::endl;
    evOld = EvtNumber;

    if (bPrintHeaders)
      {
	myfile << std::setw(10)  << ">>> Run"
	       << std::setw(10)  << "Event"
	       << std::setw(10)  << "Indx"  
	       << std::setw(10) << "Id" 
	       << std::setw(10)  << "Status"
	       << std::setw(10)  << "Mot1"
	       << std::setw(10)  << "Mot2"
	       << std::setw(10)  << "Dau1"
	       << std::setw(10)  << "Dau2"
	       // << std::setw(10)  << "px"
	       // << std::setw(10) << "py"
	       // << std::setw(10) << "pz"
	       << std::setw(10) << "E"
	       << std::setw(10) << "pT"
	       << std::setw(10) << "Eta"
	       << std::setw(10) << "Phi"
	       << std::setw(10) << "M"
	       << std::setw(10) << "Vtx-Z"
	       << std::setw(10) << "D0-Gen"
	       << std::endl;
      }
  }
  if (GenP_Mothers->at(Indx).size() > 0) moth1 = GenP_Mothers->at(Indx).at(0);
  if (GenP_Mothers->at(Indx).size() > 1) moth2 = GenP_Mothers->at(Indx).at(1);
  NDaug = GenP_Daughters->at(Indx).size();
  if (NDaug > 0){
    daug1 = GenP_Daughters->at(Indx).at(0); // First Daughter index
    if (NDaug > 1) {
      daug2 = daug1 + NDaug -1;          // Last Daughter index 
      if (GenP_Daughters->at(Indx).at(NDaug-1) != daug2) daug2 -=1;
    }
  }
  double mass = GenP_Mass->at(Indx);
  double pt   = GenP_Pt->at(Indx);
  double phi  = GenP_Phi->at(Indx);
  double eta  = GenP_Eta->at(Indx);

  TLorentzVector pGen;
  pGen.SetPtEtaPhiM(pt,eta,phi,mass);
  double ene = pGen.E();
  double px  = pGen.Px();
  double py  = pGen.Py();
  double pz  = pGen.Pz();
  myfile << std::fixed;
  myfile << std::setw(10) << std::setprecision(0) << RunNumber
	 << std::setw(10) << std::setprecision(0) << EvtNumber
	 << std::setw(10) << std::setprecision(0) << Indx
	 << std::setw(10) << std::setprecision(0) << GenP_PdgId->at(Indx)
	 << std::setw(10) << std::setprecision(0) << GenP_Status->at(Indx)
	 << std::setw(10) << std::setprecision(0) << moth1
	 << std::setw(10) << std::setprecision(0) << moth2
	 << std::setw(10) << std::setprecision(0) << daug1
	 << std::setw(10) << std::setprecision(0) << daug2
	 // << std::setw(10) << std::setprecision(3) << px
	 // << std::setw(10) << std::setprecision(3) << py
	 // << std::setw(10) << std::setprecision(3) << pz
	 << std::setw(10) << std::setprecision(3) << ene
    	 << std::setw(10) << std::setprecision(3) << pt
	 << std::setw(10) << std::setprecision(3) << eta
	 << std::setw(10) << std::setprecision(3) << phi
	 << std::setw(10) << std::setprecision(3) << mass
	 << std::setw(10) << std::setprecision(3) << GenP_VertexZ->at(Indx)
	 << std::setw(10) << std::setprecision(3) << GetGenpD0(Indx)
	 << std::endl;
}

//**************************************************
void JetProbability::GetFinalDaughter(Size_t Indx, 
				      vector<Size_t> &vFinalDau)
//**************************************************
{
  vector<unsigned short> vDau = GenP_Daughters->at(Indx);
  for(Size_t iDau=0; iDau<vDau.size(); iDau++)
    {
      Int_t DauIndx = vDau.at(iDau);
      Int_t DauStat = GenP_Status->at(DauIndx);
      Bool_t usedDau = false;

      for(Size_t ifd=0; ifd<vFinalDau.size(); ifd++ )
	{
	  if(DauIndx == vFinalDau.at(ifd))
	    usedDau = true;
	}
      if(usedDau) continue;

      if(DauStat == 1)
	{
	  vFinalDau.push_back(DauIndx);
	}
      else if(DauStat > 1)
	{
	  GetFinalDaughter(DauIndx,vFinalDau);
	}
    }
}

// //****************************************************************************
// Int_t JetProbability::MatchTrackToGenp(Size_t Index)
// //****************************************************************************
// {
//   Int_t    genID = GenP_PdgId->at(Index);
//   Double_t genEta= GenP_Eta->at(Index);
//   Double_t genPhi= GenP_Phi->at(Index);
//   Int_t TkIndex = -1;

//   for(Size_t iTP=0; iTP<TP_Pt->size(); iTP++)
//     {
//       Int_t    tpID = TP_PdgId->at(iTP);
//       if(tpID != genID) continue;
      
//       Double_t tpEta= TP_Eta->at(iTP);
//       Double_t tpPhi= TP_Phi->at(iTP);
//       Double_t dR = auxTools_.DeltaR(genEta,genPhi,
// 				    tpEta,tpPhi);
//       if(dR > 0.0001)continue;
//       TkIndex = TP_TTTrackIndex->at(iTP);
      
//     }
//   return TkIndex;
// }


//****************************************************************************
Double_t JetProbability::GetGenpD0(Size_t Index)
//****************************************************************************
{
  Double_t EventVtxX = HepMCEvt_VtxX ;
  Double_t EventVtxY = HepMCEvt_VtxY ;
  Double_t VtxX = GenP_VertexX->at(Index);
  Double_t VtxY = GenP_VertexY->at(Index);
  Double_t Phi  = GenP_Phi    ->at(Index);
  Double_t D0   = -(VtxX-EventVtxX) * TMath::Sin(Phi) + (VtxY-EventVtxY) * TMath::Cos(Phi);
  //Double_t D0   = -VtxX * TMath::Sin(Phi) + VtxY * TMath::Cos(Phi);
  return D0;
}


//****************************************************************************
Double_t JetProbability::GetTPD0(Size_t Index)
//****************************************************************************
{
  Double_t EventVtxX = HepMCEvt_VtxX ;
  Double_t EventVtxY = HepMCEvt_VtxY ;
  Double_t VtxX = TP_POCAx->at(Index);
  Double_t VtxY = TP_POCAy->at(Index);
  Double_t Phi  = TP_Phi  ->at(Index);
  Double_t D0   = -(VtxX-EventVtxX) * TMath::Sin(Phi) + (VtxY-EventVtxY) * TMath::Cos(Phi);
  //Double_t D0   = -VtxX * TMath::Sin(Phi) + VtxY * TMath::Cos(Phi);
  return D0;
}

//****************************************************************************
Double_t JetProbability::GetTPDxy(Size_t Index)
//****************************************************************************
{
  Double_t EventVtxX = HepMCEvt_VtxX ;
  Double_t EventVtxY = HepMCEvt_VtxY ;
  Double_t VtxX      = TP_POCAx->at(Index);
  Double_t VtxY      = TP_POCAy->at(Index);
  Int_t    Charge    = TP_Charge->at(Index);
  Double_t Sgn       = GetSign(Charge); 
  Double_t Dxy       = sqrt(pow(VtxX-EventVtxX,2) + pow(VtxY-EventVtxY,2)) * Sgn;
  //Double_t Dxy       = sqrt(pow(VtxX,2) + pow(VtxY,2)) * Sgn;
  return Dxy;
}

//****************************************************************************
Double_t JetProbability::GetGenpDxy(Size_t Index)
//****************************************************************************
{
  Double_t EventVtxX = HepMCEvt_VtxX ;
  Double_t EventVtxY = HepMCEvt_VtxY ;
  Double_t VtxX      = GenP_VertexX->at(Index);
  Double_t VtxY      = GenP_VertexY->at(Index);
  Int_t    Charge    = GenP_Charge->at(Index);
  Double_t Sgn       = GetSign(Charge); 
  Double_t Dxy       = sqrt(pow(VtxX-EventVtxX,2) + pow(VtxY-EventVtxY,2)) * Sgn;
  //Double_t Dxy       = sqrt(pow(VtxX,2) + pow(VtxY,2)) * Sgn;
  return Dxy;
}

//****************************************************************************
Double_t JetProbability::GetTTTrackDxy(Size_t Index)
//****************************************************************************
{
  Double_t VtxX      = L1Tks_POCAx->at(Index);
  Double_t VtxY      = L1Tks_POCAy->at(Index);
  Double_t RInv      = L1Tks_RInv->at(Index);
  Double_t Sgn       = GetSign(RInv); 
  Double_t Dxy       = sqrt(pow(VtxX,2) + pow(VtxY,2)) * Sgn;
 
  return Dxy;
}

//****************************************************************************
Double_t JetProbability::GetTTPixTrackDxy(Size_t Index)
//****************************************************************************
{
  Double_t VtxX      = L1PixTks_POCAx->at(Index);
  Double_t VtxY      = L1PixTks_POCAy->at(Index);
  Double_t RInv      = L1PixTks_RInv->at(Index);
  Double_t Sgn       = GetSign(RInv); 
  Double_t Dxy       = sqrt(pow(VtxX,2) + pow(VtxY,2)) * Sgn;
  return Dxy;
}



//****************************************************************************
Int_t JetProbability::GetGenMatchingTPIndex(Int_t Index)
//****************************************************************************
{
  Size_t   nTP    = TP_PdgId->size(); 
  Int_t    genId  = GenP_PdgId->at(Index);
  Double_t genEta = GenP_Eta->at(Index);
  Double_t genPhi = GenP_Phi->at(Index);

  Int_t tp_Indx = -1;

  for(Size_t iTP=0; iTP<nTP; iTP++)
    {
      Int_t    tp_Id  = TP_PdgId->at(iTP);
      Double_t tp_Eta = TP_Eta->at(iTP);
      Double_t tp_Phi = TP_Phi->at(iTP);
      Double_t tp_D0  = GetTPD0(iTP);
      if(tp_Id == genId)
	{
	  Double_t dR = auxTools_.DeltaR(tp_Eta,tp_Phi,
					 genEta,genPhi);
	  if(dR<0.0001)
	    {
	      //	      cout<<"matching TP  "<< iTP <<"\tid\t"<< tp_Id<<"\t\t"<<dR<<endl;
	      tp_Indx =iTP;
	    }
	}
    }
  return tp_Indx; 
}

//****************************************************************************
Int_t JetProbability::GetSign(Double_t Quantity)
//****************************************************************************
{
  Int_t Sign =0;
  if (Quantity>0)
    Sign = 1;
  else if (Quantity<0)
    Sign = -1;
  return Sign;
}

//****************************************************************************
Bool_t JetProbability::IsBHadron(Int_t Index)
//****************************************************************************
{
  Int_t Gen_ID = GenP_PdgId->at(Index);
   Int_t Mod = abs(Gen_ID)%10000;
  Int_t Had = 0;
  Bool_t isBHad= false;
  if(Mod > 1000)
    Had=Mod/1000;
  else if(Mod < 1000)
    Had=Mod/100;
  
  if(Had == 5)
    isBHad = true;

  return isBHad;
}

//****************************************************************************
Bool_t JetProbability::IsFinalBHadron(Int_t Index)
//****************************************************************************
{
  Int_t Gen_ID = GenP_PdgId->at(Index);
  Int_t Mod = abs(Gen_ID)%10000;
  Int_t Had = 0;
  Bool_t isFinalBHad = true;
  if(Mod > 1000)
    Had=Mod/1000;
  else if(Mod < 1000)
    Had=Mod/100;
  
  if(Had != 5)
    isFinalBHad = false; // not BHad
  if(Had == 5)
    {
      vector <unsigned short> vBHadDau;
      vBHadDau = GenP_Daughters->at(Index);
      for(Size_t iDau=0; iDau<vBHadDau.size(); iDau++)
	{
	  Size_t DauIndx = vBHadDau.at(iDau);
	  Int_t  DauId   = GenP_PdgId->at(DauIndx);
	  Int_t  DauMod  = abs(DauId)%10000;
	  Int_t  DauHad  = 0;
	  if(DauMod > 1000)
	    DauHad=DauMod/1000;
	  else if(DauMod < 1000)
	    DauHad=DauMod/100;
	  if(DauHad == Had)
	    isFinalBHad = false; // not final BHAd
	}// loop on Dau
    }


  return isFinalBHad;
}

//****************************************************************************
Bool_t JetProbability::IsLiteHadron(Int_t Index)
//****************************************************************************
{
  Int_t Gen_ID = GenP_PdgId->at(Index);
   Int_t Mod = abs(Gen_ID)%10000;
  Int_t Had = 0;
  Bool_t isLHad= false;
  if(Mod > 1000)
    Had=Mod/1000;
  else if(Mod < 1000)
    Had=Mod/100;
  
  if(Had == 1 || Had == 2 || Had == 3)
    isLHad = true;

  return isLHad;
}

//****************************************************************************
void JetProbability::GetBHadronChargedImmediateDaughters(Int_t Index ,
							 vector<Size_t> &vCIDau)
//****************************************************************************
{
  if(!IsFinalBHadron(Index))
    {
      cout<< "Warning:: This is not a Final B Hadron"<<endl;
    }
  vector <unsigned short> vBHadDau;
  vBHadDau = GenP_Daughters->at(Index);
  for(Size_t iDau=0; iDau<vBHadDau.size(); iDau++)
    {
      Int_t DauIndex  = vBHadDau.at(iDau);
      Int_t DauID     = GenP_PdgId->at(DauIndex);
      Int_t DauStatus = GenP_Status->at(DauIndex);
      Int_t DauCharge = GenP_Charge->at(DauIndex);
      if(DauStatus == 1 )
	{
	  if(DauCharge != 0)
	    {
	      vCIDau.push_back(DauIndex);
	    }
	}
      //		cout<< DauIndex<<"\t"<<DauID<<endl;
    }
  

}

//****************************************************************************
Double_t JetProbability::GetPhiD0(Double_t PhiTk, Double_t D0)
//****************************************************************************
{
  Double_t SignD0 = GetSign(D0);
  Double_t PhiD0  = PhiTk + SignD0*(PI/2) ;
  
  if(PhiD0 > PI)
    PhiD0 = PhiD0 - 2*PI;
  else if(PhiD0 <= -PI)
    PhiD0 = PhiD0 + 2*PI;

  return PhiD0;

}
//****************************************************************************
Double_t JetProbability::GetSignofIP(Double_t deltaPhi)
//****************************************************************************
{
  Double_t Sign =0;
  if(deltaPhi<PI/2)
    Sign = 1.0;
  else if (deltaPhi>PI/2)
    Sign = -1.0;
  
  return Sign;
}


//****************************************************************************
void JetProbability::FillSIPHistogramsD0(Int_t DauIndx, Double_t had_Phi,
 					 TH1D* GenAbsD0, TH1D* GenPhiD0, TH1D* GenDeltaPhi,TH1D* GenSIP,
 					 TH1D* TPAbsD0 , TH1D* TPPhiD0 , TH1D* TPDeltaPhi ,TH1D* TPSIP,
 					 TH1D* TkAbsD0 , TH1D* TkPhiD0 , TH1D* TkDeltaPhi ,TH1D* TkSIP,
 					 TH1D* PixAbsD0, TH1D* PixPhiD0, TH1D* PixDeltaPhi,TH1D* PixSIP)
//void JetProbability::FillSIPHistogramsD0(Int_t DauIndx, Double_t had_Phi)
//****************************************************************************
{
  Double_t DauPt    = GenP_Pt->at(DauIndx);
  Double_t DauEta   = GenP_Eta->at(DauIndx);
  Double_t DauPhi   = GenP_Phi->at(DauIndx);
  
  Double_t DauD0    = GetGenpD0(DauIndx);
  Double_t DauD0Phi = GetPhiD0(DauPhi, DauD0);
  Double_t deltaPhi = abs(auxTools_.DeltaPhi(DauD0Phi,had_Phi));
  Double_t DauAbsD0 = abs(DauD0);
  Double_t Sign     = GetSignofIP(deltaPhi);
  Double_t DauSIP   = DauAbsD0 * Sign;
  
  
  if(DauPt >2 && abs(DauEta)<2.5)
    {
      GenAbsD0->Fill(DauAbsD0);
      GenPhiD0->Fill(DauD0Phi);
      GenDeltaPhi->Fill(deltaPhi);
      GenSIP->Fill(DauSIP);
      
      Int_t MatchedTP = GetGenMatchingTPIndex(DauIndx);
      if(MatchedTP > -1)
	{
	  Double_t tp_D0    = GetTPD0(MatchedTP);
	  Double_t tp_Phi   = TP_Phi->at(MatchedTP);
	  Double_t tp_PhiD0 = GetPhiD0(tp_Phi, tp_D0);
	  Double_t tp_dPhi  = abs(auxTools_.DeltaPhi(tp_PhiD0,had_Phi));
	  Double_t tp_AbsD0 = abs(tp_D0);
	  Double_t tp_Sign  = GetSignofIP(tp_dPhi);
	  Double_t tp_SIP   = tp_AbsD0 * tp_Sign;
	  

	  TPAbsD0->Fill(tp_AbsD0);
	  TPPhiD0->Fill(tp_PhiD0);
	  TPDeltaPhi->Fill(tp_dPhi);
	  TPSIP->Fill(tp_SIP);

	  
	  Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
	  if(tp_tkIndx>-1)
	    {
	      Double_t tk_D0    = s->GetTrackD0(tp_tkIndx) * -1;
	      Double_t tk_Phi   = L1Tks_Phi->at(tp_tkIndx);
	      Double_t tk_PhiD0 = GetPhiD0(tk_Phi, tk_D0);
	      Double_t tk_dPhi  = abs(auxTools_.DeltaPhi(tk_PhiD0,had_Phi));
	      Double_t tk_AbsD0 = abs(tk_D0);
	      Double_t tk_Sign  = GetSignofIP(tk_dPhi);
	      Double_t tk_SIP   = tk_AbsD0 * tk_Sign;
	      
	      
	      TkAbsD0   ->Fill(tk_AbsD0);
	      TkPhiD0   ->Fill(tk_PhiD0);
	      TkDeltaPhi->Fill(tk_dPhi);
	      TkSIP     ->Fill(tk_SIP);

	      
	      Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
	      if(tk_pixIndx>-1)
		{
		  Double_t pix_D0    = s->GetPixTrackD0(tk_pixIndx);
		  Double_t pix_Pt    = L1PixTks_Pt->at(tk_pixIndx);
		  Double_t pix_Eta   = L1PixTks_Eta->at(tk_pixIndx);
		  Double_t pix_Phi   = L1PixTks_Phi->at(tk_pixIndx);
		  Double_t pix_Chi2  = L1PixTks_ChiSquared->at(tk_pixIndx);
		  Double_t pix_PhiD0 = GetPhiD0(pix_Phi, pix_D0);
		  Double_t pix_dPhi  = abs(auxTools_.DeltaPhi(pix_PhiD0,had_Phi));
		  Double_t pix_AbsD0 = abs(pix_D0);
		  Double_t pix_Sign  = GetSignofIP(pix_dPhi);
		  Double_t pix_SIP   = pix_AbsD0 * pix_Sign;

		  PixAbsD0   ->Fill(pix_AbsD0);
		  PixPhiD0   ->Fill(pix_PhiD0);
		  PixDeltaPhi->Fill(pix_dPhi);
		  PixSIP     ->Fill(pix_SIP);

		  if(IsTrackChi2Good("Tight",pix_Chi2))
		    {
		      PixDeltaPhi_T_c2->Fill(pix_dPhi);
		      PixSIP_T_c2     ->Fill(pix_SIP);
		    }
		  if(IsTrackChi2Good("Medium",pix_Chi2))
		    {
		      PixDeltaPhi_M_c2->Fill(pix_dPhi);
		      PixSIP_M_c2     ->Fill(pix_SIP);
		    }
		  if(IsTrackChi2Good("Loose",pix_Chi2))
		    {
		      PixDeltaPhi_L_c2->Fill(pix_dPhi);
		      PixSIP_L_c2     ->Fill(pix_SIP);
		    }
		  if(IsTrackChi2Good("VLoose",pix_Chi2))
		    {
		      PixDeltaPhi_VL_c2->Fill(pix_dPhi);
		      PixSIP_VL_c2     ->Fill(pix_SIP);
		    }
		  
		  //**************************************************************************************** Central
		  if(IsWithinEtaRegion("Central",pix_Eta))
		    {
		      PixDeltaPhi_C->Fill(pix_dPhi);
		      PixSIP_C     ->Fill(pix_SIP);
		      if(IsWithinPtSlice("SL0203",pix_Pt))
			{
			  PixDeltaPhi_C_0203->Fill(pix_dPhi);
			  PixSIP_C_0203     ->Fill(pix_SIP);
			}
		      else if(IsWithinPtSlice("SL0304",pix_Pt))
			{
			  PixDeltaPhi_C_0304->Fill(pix_dPhi);
			  PixSIP_C_0304     ->Fill(pix_SIP);
			}
		      else if(IsWithinPtSlice("SL0405",pix_Pt))
			{
			  PixDeltaPhi_C_0405->Fill(pix_dPhi);
			  PixSIP_C_0405     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0506",pix_Pt))
			{
			  PixDeltaPhi_C_0506->Fill(pix_dPhi);
			  PixSIP_C_0506     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0607",pix_Pt))
			{
			  PixDeltaPhi_C_0607->Fill(pix_dPhi);
			  PixSIP_C_0607     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0708",pix_Pt))
			{
			  PixDeltaPhi_C_0708->Fill(pix_dPhi);
			  PixSIP_C_0708     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0809",pix_Pt))
			{
			  PixDeltaPhi_C_0809->Fill(pix_dPhi);
			  PixSIP_C_0809     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0910",pix_Pt))
			{
			  PixDeltaPhi_C_0910->Fill(pix_dPhi);
			  PixSIP_C_0910     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL1015",pix_Pt))
			{
			  PixDeltaPhi_C_1015->Fill(pix_dPhi);
			  PixSIP_C_1015     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL1520",pix_Pt))
			{
			  PixDeltaPhi_C_1520->Fill(pix_dPhi);
			  PixSIP_C_1520     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL2000",pix_Pt))
			{
			  PixDeltaPhi_C_2000->Fill(pix_dPhi);
			  PixSIP_C_2000     ->Fill(pix_SIP);
			} 
		    }// end if eta central
		  //**************************************************************************************** Intermediate
		  else if(IsWithinEtaRegion("Intermediate",pix_Eta))
		    {
		      PixDeltaPhi_I->Fill(pix_dPhi);
		      PixSIP_I     ->Fill(pix_SIP);
		      if(IsWithinPtSlice("SL0203",pix_Pt))
			{
			  PixDeltaPhi_I_0203->Fill(pix_dPhi);
			  PixSIP_I_0203     ->Fill(pix_SIP);
			}
		      else if(IsWithinPtSlice("SL0304",pix_Pt))
			{
			  PixDeltaPhi_I_0304->Fill(pix_dPhi);
			  PixSIP_I_0304     ->Fill(pix_SIP);
			}
		      else if(IsWithinPtSlice("SL0405",pix_Pt))
			{
			  PixDeltaPhi_I_0405->Fill(pix_dPhi);
			  PixSIP_I_0405     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0506",pix_Pt))
			{
			  PixDeltaPhi_I_0506->Fill(pix_dPhi);
			  PixSIP_I_0506     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0607",pix_Pt))
			{
			  PixDeltaPhi_I_0607->Fill(pix_dPhi);
			  PixSIP_I_0607     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0708",pix_Pt))
			{
			  PixDeltaPhi_I_0708->Fill(pix_dPhi);
			  PixSIP_I_0708     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0809",pix_Pt))
			{
			  PixDeltaPhi_I_0809->Fill(pix_dPhi);
			  PixSIP_I_0809     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0910",pix_Pt))
			{
			  PixDeltaPhi_I_0910->Fill(pix_dPhi);
			  PixSIP_I_0910     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL1015",pix_Pt))
			{
			  PixDeltaPhi_I_1015->Fill(pix_dPhi);
			  PixSIP_I_1015     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL1520",pix_Pt))
			{
			  PixDeltaPhi_I_1520->Fill(pix_dPhi);
			  PixSIP_I_1520     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL2000",pix_Pt))
			{
			  PixDeltaPhi_I_2000->Fill(pix_dPhi);
			  PixSIP_I_2000     ->Fill(pix_SIP);
			} 
		    }//end if eta intermediate
		  //**************************************************************************************** Forward
		  else if(IsWithinEtaRegion("Forward",pix_Eta))
		    {
		      PixDeltaPhi_F->Fill(pix_dPhi);
		      PixSIP_F     ->Fill(pix_SIP);
		      if(IsWithinPtSlice("SL0203",pix_Pt))
			{
			  PixDeltaPhi_F_0203->Fill(pix_dPhi);
			  PixSIP_F_0203     ->Fill(pix_SIP);
			}
		      else if(IsWithinPtSlice("SL0304",pix_Pt))
			{
			  PixDeltaPhi_F_0304->Fill(pix_dPhi);
			  PixSIP_F_0304     ->Fill(pix_SIP);
			}
		      else if(IsWithinPtSlice("SL0405",pix_Pt))
			{
			  PixDeltaPhi_F_0405->Fill(pix_dPhi);
			  PixSIP_F_0405     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0506",pix_Pt))
			{
			  PixDeltaPhi_F_0506->Fill(pix_dPhi);
			  PixSIP_F_0506     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0607",pix_Pt))
			{
			  PixDeltaPhi_F_0607->Fill(pix_dPhi);
			  PixSIP_F_0607     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0708",pix_Pt))
			{
			  PixDeltaPhi_F_0708->Fill(pix_dPhi);
			  PixSIP_F_0708     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0809",pix_Pt))
			{
			  PixDeltaPhi_F_0809->Fill(pix_dPhi);
			  PixSIP_F_0809     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL0910",pix_Pt))
			{
			  PixDeltaPhi_F_0910->Fill(pix_dPhi);
			  PixSIP_F_0910     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL1015",pix_Pt))
			{
			  PixDeltaPhi_F_1015->Fill(pix_dPhi);
			  PixSIP_F_1015     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL1520",pix_Pt))
			{
			  PixDeltaPhi_F_1520->Fill(pix_dPhi);
			  PixSIP_F_1520     ->Fill(pix_SIP);
			} 
		      else if(IsWithinPtSlice("SL2000",pix_Pt))
			{
			  PixDeltaPhi_F_2000->Fill(pix_dPhi);
			  PixSIP_F_2000     ->Fill(pix_SIP);
			} 
		    }//end if eta forward

		} // if ttpixtrack found
	    }// if tttrack found
	}//if TP found
    }// pt eta cut on Genp 
  
}


//****************************************************************************
void JetProbability::FillSIPHistogramsDxy(Int_t DauIndx, Double_t had_Phi,
					  TH1D* GenAbsD0, TH1D* GenPhiD0, TH1D* GenDeltaPhi,TH1D* GenSIP,
					  TH1D* TPAbsD0 , TH1D* TPPhiD0 , TH1D* TPDeltaPhi ,TH1D* TPSIP,
					  TH1D* TkAbsD0 , TH1D* TkPhiD0 , TH1D* TkDeltaPhi ,TH1D* TkSIP,
					  TH1D* PixAbsD0, TH1D* PixPhiD0, TH1D* PixDeltaPhi,TH1D* PixSIP)
//****************************************************************************
{
  Double_t DauPt    = GenP_Pt->at(DauIndx);
  Double_t DauEta   = GenP_Eta->at(DauIndx);
  Double_t DauPhi   = GenP_Phi->at(DauIndx);
  
  Double_t DauD0    = GetGenpDxy(DauIndx);
  Double_t DauD0Phi = GetPhiD0(DauPhi, DauD0);
  Double_t deltaPhi = abs(auxTools_.DeltaPhi(DauD0Phi,had_Phi));
  Double_t DauAbsD0 = abs(DauD0);
  Double_t Sign     = GetSignofIP(deltaPhi);
  Double_t DauSIP   = DauAbsD0 * Sign;
  
  
  if(DauPt >2 && abs(DauEta)<2.5)
    {
      GenAbsD0->Fill(DauAbsD0);
      GenPhiD0->Fill(DauD0Phi);
      GenDeltaPhi->Fill(deltaPhi);
      GenSIP->Fill(DauSIP);
      
      Int_t MatchedTP = GetGenMatchingTPIndex(DauIndx);
      if(MatchedTP > -1)
	{
	  Double_t tp_D0    = GetTPDxy(MatchedTP);
	  Double_t tp_Phi   = TP_Phi->at(MatchedTP);
	  Double_t tp_PhiD0 = GetPhiD0(tp_Phi, tp_D0);
	  Double_t tp_dPhi  = abs(auxTools_.DeltaPhi(tp_PhiD0,had_Phi));
	  Double_t tp_AbsD0 = abs(tp_D0);
	  Double_t tp_Sign  = GetSignofIP(tp_dPhi);
	  Double_t tp_SIP   = tp_AbsD0 * tp_Sign;
	  

	  TPAbsD0->Fill(tp_AbsD0);
	  TPPhiD0->Fill(tp_PhiD0);
	  TPDeltaPhi->Fill(tp_dPhi);
	  TPSIP->Fill(tp_SIP);

	  
	  Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
	  if(tp_tkIndx>-1)
	    {
	      Double_t tk_D0    = GetTTTrackDxy(tp_tkIndx);
	      Double_t tk_Phi   = L1Tks_Phi->at(tp_tkIndx);
	      Double_t tk_PhiD0 = GetPhiD0(tk_Phi, tk_D0);
	      Double_t tk_dPhi  = abs(auxTools_.DeltaPhi(tk_PhiD0,had_Phi));
	      Double_t tk_AbsD0 = abs(tk_D0);
	      Double_t tk_Sign  = GetSignofIP(tk_dPhi);
	      Double_t tk_SIP   = tk_AbsD0 * tk_Sign;
	      
	      
	      TkAbsD0   ->Fill(tk_AbsD0);
	      TkPhiD0   ->Fill(tk_PhiD0);
	      TkDeltaPhi->Fill(tk_dPhi);
	      TkSIP     ->Fill(tk_SIP);

	      
	      Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
	      if(tk_pixIndx>-1)
		{
		  Double_t pix_D0    = GetTTPixTrackDxy(tk_pixIndx);
		  Double_t pix_Phi   = L1PixTks_Phi->at(tk_pixIndx);
		  Double_t pix_PhiD0 = GetPhiD0(pix_Phi, pix_D0);
		  Double_t pix_dPhi  = abs(auxTools_.DeltaPhi(pix_PhiD0,had_Phi));
		  Double_t pix_AbsD0 = abs(pix_D0);
		  Double_t pix_Sign  = GetSignofIP(pix_dPhi);
		  Double_t pix_SIP   = pix_AbsD0 * pix_Sign;

		  PixAbsD0   ->Fill(pix_AbsD0);
		  PixPhiD0   ->Fill(pix_PhiD0);
		  PixDeltaPhi->Fill(pix_dPhi);
		  PixSIP     ->Fill(pix_SIP);

		} // if ttpixtrack found
	    }// if tttrack found
	}//if TP found
    }// pt eta cut on Genp 
  
}



//****************************************************************************
void JetProbability::FillSIPHistogramsDotProd(Int_t DauIndx , TVector3 Hadvec,
					      TH1D* GenAbsD0, TH1D* GenSIP,
					      TH1D* TPAbsD0 , TH1D* TPSIP,
					      TH1D* TkAbsD0 , TH1D* TkSIP,
					      TH1D* PixAbsD0, TH1D* PixSIP)
//****************************************************************************
{
  Double_t DauPt    = GenP_Pt->at(DauIndx);
  Double_t DauEta   = GenP_Eta->at(DauIndx);
  Double_t DauPhi   = GenP_Phi->at(DauIndx);
  Double_t DauX     = GenP_VertexX->at(DauIndx)-HepMCEvt_VtxX;
  Double_t DauY     = GenP_VertexX->at(DauIndx)-HepMCEvt_VtxY;
  Double_t DauD0    = GetGenpD0(DauIndx);
  TVector3 D0Vec(DauX,DauY,0);
  Double_t dotProd  = D0Vec.Dot(Hadvec);
  Double_t Sign     = GetSign(dotProd);
  Double_t DauAbsD0 = abs(DauD0);
  Double_t DauSIP   = DauAbsD0 * Sign;
  
  if(DauPt >2 && abs(DauEta)<2.5)
    {
      GenAbsD0->Fill(DauAbsD0);
      GenSIP->Fill(DauSIP);
      
      Int_t MatchedTP = GetGenMatchingTPIndex(DauIndx);
      if(MatchedTP > -1)
	{
	  Double_t tp_D0      = GetTPD0(MatchedTP);
	  Double_t tp_X       = TP_POCAx->at(MatchedTP)-HepMCEvt_VtxX;
	  Double_t tp_Y       = TP_POCAy->at(MatchedTP)-HepMCEvt_VtxY;
	  TVector3 tp_D0Vec(tp_X,tp_Y,0);
	  Double_t tp_dotProd = tp_D0Vec.Dot(Hadvec);
	  Double_t tp_Sign    = GetSign(tp_dotProd);
	  Double_t tp_AbsD0   = abs(tp_D0);
	  Double_t tp_SIP     = tp_AbsD0 * tp_Sign;
	  

	  TPAbsD0->Fill(tp_AbsD0);
	  TPSIP->Fill(tp_SIP);

	  
	  Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
	  if(tp_tkIndx>-1)
	    {
	      Double_t tk_D0    = s->GetTrackD0(tp_tkIndx) * -1;
	      Double_t tk_X   = L1Tks_POCAx->at(tp_tkIndx);
	      Double_t tk_Y   = L1Tks_POCAy->at(tp_tkIndx);
	      TVector3 tk_D0Vec(tk_X,tk_Y,0);
	      Double_t tk_dotProd = tk_D0Vec.Dot(Hadvec);
	      Double_t tk_Sign     = GetSign(tk_dotProd);
	      Double_t tk_AbsD0 = abs(tk_D0);
	      Double_t tk_SIP   = tk_AbsD0 * tk_Sign;
	      
	      TkAbsD0   ->Fill(tk_AbsD0);
	      TkSIP     ->Fill(tk_SIP);

	      
	      Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
	      if(tk_pixIndx>-1)
		{
		  Double_t pix_D0      = s->GetPixTrackD0(tk_pixIndx);		  
		  Double_t pix_X       = L1PixTks_POCAx->at(tk_pixIndx);
		  Double_t pix_Y       = L1PixTks_POCAy->at(tk_pixIndx);
		  TVector3 pix_D0Vec(pix_X,pix_Y,0);
		  Double_t pix_dotProd = pix_D0Vec.Dot(Hadvec);
		  Double_t pix_Sign    = GetSign(pix_dotProd);
		  Double_t pix_AbsD0   = abs(pix_D0);
		  Double_t pix_SIP     = pix_AbsD0 * pix_Sign;


		  PixAbsD0   ->Fill(pix_AbsD0);
		  PixSIP     ->Fill(pix_SIP);

		} // if ttpixtrack found
	    }// if tttrack found
	}//if TP found
    }// pt eta cut on Genp 
  
}



#endif // JetProbability_cxx


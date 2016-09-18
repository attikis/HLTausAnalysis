vector<Size_t> vFinalDau;

    vector<Size_t> promptEMu;
    vector<Size_t> promptUSDdau;
   
    
    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	Int_t Gen_ID = GenP_PdgId->at(iGenp);
	Double_t Gen_D0 = GetGenpD0(iGenp);
	vector<unsigned short> Gen_vDaug  = GenP_Daughters->at(iGenp);
	if(abs(Gen_ID) == 11 || abs(Gen_ID) == 13 ) // e mu in event
	  {
	    if(RecursivelyLookForMotherId(iGenp,24, false)) // look is these are coming from W.
	      {
		if(!RecursivelyLookForMotherId(iGenp,5, false)) // remove leptons from b...
		  {
		    if(!RecursivelyLookForMotherId(iGenp,4, false))// remove lep from  c
		      {
			if(!RecursivelyLookForMotherId(iGenp,15, false))// remove lep from  tau
			  {
			    if(IsFinalGenp (iGenp, Gen_vDaug))  // selecting particle only once 
			      {
				//		cout<< iGenp<<"\t"<<Gen_ID<<"\t\t"<<Gen_D0<<endl;
				hgend0EMu->Fill(abs(Gen_D0));
				promptEMu.push_back(iGenp);
			      } // end if final in list
			  }//
		      } // end if not from c
		  } // end if not from b
	      }// end if from W
	  }// if e mu 

	if(abs(Gen_ID) == 1 || abs(Gen_ID) == 2 || abs(Gen_ID) == 3) // UDS quarks in event.
	  {
	    if(RecursivelyLookForMotherId(iGenp,24, false)) // look is these are coming from W.
	      {
		if(!RecursivelyLookForMotherId(iGenp,5, false)) // remove leptons from b...
		  {
		    if(!RecursivelyLookForMotherId(iGenp,4, false))// remove lep from  c
		      {
			if(IsFinalGenp (iGenp, Gen_vDaug))  // selecting particle only once 
			  {
			    GetFinalDaughter(iGenp,vFinalDau);
			  }
		      }
		  }
	      }
	  }

	//***************************************************
	//HF BHAD DAU
	//***************************************************
	Int_t Mod = abs(Gen_ID)%10000;
	Int_t Had = 0;
	if(Mod > 1000)
	  Had=Mod/1000;
	else if(Mod < 1000)
	  Had=Mod/100;
	
	Bool_t isFinalBhad=true;
	if (Had==5)
	  {
	    Int_t bHadId = GenP_PdgId->at(iGenp);
	    vector <unsigned short> vBHadDau;
	    vBHadDau = GenP_Daughters->at(iGenp);
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
		  isFinalBhad = false;
	      }// loop on Dau
	    if(isFinalBhad)
	      {
		Double_t bHadPhi = GenP_Phi->at(iGenp);
		for(Size_t iDau=0; iDau<vBHadDau.size(); iDau++)
		  {
		    Size_t   DauIndx   = vBHadDau.at(iDau);
		    Int_t    DauId     = GenP_PdgId->at(DauIndx);
		    Int_t    DauStatus = GenP_Status->at(DauIndx);
		    Int_t    DauCharge = GenP_Charge->at(DauIndx);
		    Double_t DauD0     = GetGenpD0(DauIndx);
		    Double_t DauPt     = GenP_Pt->at(DauIndx);
		    Double_t DauEta    = GenP_Eta->at(DauIndx);

		    //***********************************************************************
		    //This part fills the histograms for AbsD0
		    //***********************************************************************
		    if(DauPt>2)
		      {
			if(abs(DauEta)<2.5)
			  {
			    if(DauCharge !=0)
			      {
				if(DauStatus == 1/* && abs(DauId)==13*/ )
				  {
				    //			    cout<<DauIndx<<"\t"<<DauId<<"\t"<<DauStatus<<endl;	    
				    hgend0BHadDirect->Fill(abs(DauD0));

				    Int_t MatchedTP = GetGenMatchingTPIndex(DauIndx);
				    //cout<<"\t\t\tMAtched TP \t\t"<<MatchedTP<<endl;
				    if(MatchedTP > -1)
				      {
					Double_t tp_D0 = GetTPD0(MatchedTP);
					hGenMatchedTPB->Fill(abs(tp_D0));
					Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
					if(tp_tkIndx>-1)
					  {
					    Double_t TTtrackD0 = s->GetTrackD0(tp_tkIndx);
					    // Double_t temp = -L1Tks_POCAx->at(tp_tkIndx)*TMath::Sin(L1Tks_Phi->at(tp_tkIndx)) + L1Tks_POCAy->at(tp_tkIndx)*TMath::Cos(L1Tks_Phi->at(tp_tkIndx)); 
					    // if(abs(TTtrackD0)>0.3)
					    //   cout<<tp_D0 <<"\t"<< TTtrackD0 <<"\t\t"<<temp<<endl;
					    hGenMatchedTkB->Fill(abs(TTtrackD0));
					    Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
					    if(tk_pixIndx>-1)
					      {
						Double_t TTPixtrackD0 = s->GetPixTrackD0(tk_pixIndx);
						hGenMatchedPixB->Fill(abs(TTPixtrackD0));
						// cout<< tk_pixIndx<<"\t"<<L1PixTks_Pt->at(tk_pixIndx)<<"\t"<<L1Tks_Pt->at(tp_tkIndx)<<"\t"<<TP_Pt->at(MatchedTP)<<endl;
					      } // if ttpixtrack found
					  }// if tttrack found
					//  cout<<"\tTP Index\t"<< MatchTP<< "\tTP ID\t"<< TP_PdgId->at(MatchTP)<<endl;
				      }//if TP found
				    
				  }// if status
			      }//if charge
			  }// if eta
		      }//if daupt
		    //***********************************************************************
		    
		    if(DauCharge !=0)
		      {
			if(DauStatus == 1/* && abs(DauId)==13*/ )
			  {
			    
			    Double_t DauPhi = GenP_Phi->at(DauIndx);
			    Double_t cosDauPhi = TMath::Cos(DauPhi);
			    

			    Int_t signPhi    = GetSign(DauPhi); 
			    Int_t signcosPhi = GetSign(cosDauPhi);
			    
			    Double_t DauPhiD0 = DauPhi + (signcosPhi*signPhi*PI/2) ;

			    // if(cosDauPhi>0)
			    //   DauPhiD0 = DauPhi+(PI/2);
			    // else if (cosDauPhi<0)
			    //   DauPhiD0 = DauPhi-(PI/2);
			    
			    // if(DauPhiD0 >  PI)  DauPhiD0 -= PI;
			    // if(DauPhiD0 <= -PI) DauPhiD0 += PI;
			    
			    
			    
			    hgenphiHadDau   ->Fill(DauPhi);
			    
			    hgenphiD0HadDau ->Fill(DauPhiD0);

			    Double_t DeltaPhi = abs(auxTools_.DeltaPhi( bHadPhi,DauPhiD0));
			    //Double_t DeltaPhi = abs(bHadPhi-DauPhiD0);
			    hgendeltaphiD0HadDau ->Fill(DeltaPhi);
			    
			    //	    cout<<"DauPhi\t"<<DauPhi<<"\tD0Phi\t"<<DauPhiD0<<"\tHadPhi\t"<<bHadPhi<<"\tDeltaPhi\t"<<DeltaPhi<<endl;
			    // cout<<DauIndx<<"\t"<<DauId<<"\t"<< DauPhi<<"\t"<<cosDauPhi<<endl;
			  }// is status1
		      }// id charged
		  }// end loop on BHAD dau
		//		cout<<iGenp<<"\t"<< bHadId <<endl;
	      }// end if final BHAD
	  }// end  if BHaf
      }// end loop on Genp

    for(Size_t iFD = 0; iFD<vFinalDau.size(); iFD++)
      {
	Int_t FDindx   = vFinalDau.at(iFD); 
	Int_t FDcharge = GenP_Charge->at(FDindx);

	if(FDcharge != 0)
	  {
	    if(!RecursivelyLookForMotherId(FDindx,5, false)) // remove leptons from b...
	      {
		if(!RecursivelyLookForMotherId(FDindx,4, false)) // remove leptons from c...
		  {
		    Double_t Gen_D0 = GetGenpD0(FDindx);
		    //		    cout<< FDindx<<"\tusd\t"<<GenP_PdgId->at(FDindx)<<"\t\t"<<Gen_D0<<endl;
		    hgend0UDS->Fill(abs(Gen_D0));
		    promptUSDdau.push_back(FDindx);
		  } // remove from c
	      }// remove from b
	  }// only charged
      }// end loop on daugheters of prompt  usd from W

    //    cout<<"promptEMu = "<< promptEMu.size()<<"\tpromptUDS\t"<<promptUSDdau.size()<<endl;
    Size_t npEMu = promptEMu.size();
    Size_t npUSD = promptUSDdau.size();
    for(Size_t ipEMu=0; ipEMu<npEMu; ipEMu++)
      {
	Size_t   pEMuIndex = promptEMu.at(ipEMu);
	
	//	cout<<"from Vector EMu----- "<<pEMuIndex<<"\t"<<GenP_PdgId->at(pEMuIndex)<<endl;
	//	cout<<ipEMu<<endl;
	Int_t MatchedTP = GetGenMatchingTPIndex(pEMuIndex);
	if(MatchedTP > -1)
	  {
	    Double_t tp_D0 = GetTPD0(MatchedTP);
	    hGenMatchedTPMuE->Fill(abs(tp_D0));
	    Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
	    if(tp_tkIndx>-1)
	      {
		Double_t TTtrackD0 = s->GetTrackD0(tp_tkIndx);
		//	Double_t temp = -L1Tks_POCAx->at(tp_tkIndx)*TMath::Sin(L1Tks_Phi->at(tp_tkIndx)) + L1Tks_POCAy->at(tp_tkIndx)*TMath::Cos(L1Tks_Phi->at(tp_tkIndx)); 
		//cout<<tp_D0 <<"\t"<< TTtrackD0 <<"\t\t"<<temp<<endl;
		hGenMatchedTkMuE->Fill(abs(TTtrackD0));
		
		Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
		if(tk_pixIndx>-1)
		  {
		    Double_t TTPixtrackD0 = s->GetPixTrackD0(tk_pixIndx);
		    hGenMatchedPixMuE->Fill(abs(TTPixtrackD0));
		    // cout<< tk_pixIndx<<"\t"<<L1PixTks_Pt->at(tk_pixIndx)<<"\t"<<L1Tks_Pt->at(tp_tkIndx)<<"\t"<<TP_Pt->at(MatchedTP)<<endl;
		  } // if ttpixtrack found
	      }// if tttrack found
	    //  cout<<"\tTP Index\t"<< MatchTP<< "\tTP ID\t"<< TP_PdgId->at(MatchTP)<<endl;
	  }

      }// end loop on EMu Dau
    
    for(Size_t ipUSD=0; ipUSD<npUSD; ipUSD++)
      {
	Size_t   pUSDIndex = promptUSDdau.at(ipUSD);
	
	//	cout<<"from Vector USD----- "<<pUSDIndex<<"\t"<<GenP_PdgId->at(pUSDIndex)<<endl;
	Int_t MatchedTP = GetGenMatchingTPIndex(pUSDIndex);
	//cout<<"\t\t\tMAtched TP \t\t"<<MatchedTP<<endl;
	if(MatchedTP > -1)
	  {
	    Double_t tp_D0 = GetTPD0(MatchedTP);
	    hGenMatchedTPP->Fill(abs(tp_D0));
	    Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
	    if(tp_tkIndx>-1)
	      {
		Double_t TTtrackD0 = s->GetTrackD0(tp_tkIndx);
		// Double_t temp = -L1Tks_POCAx->at(tp_tkIndx)*TMath::Sin(L1Tks_Phi->at(tp_tkIndx)) + L1Tks_POCAy->at(tp_tkIndx)*TMath::Cos(L1Tks_Phi->at(tp_tkIndx)); 
		// if(abs(TTtrackD0)>0.3)
		//   cout<<tp_D0 <<"\t"<< TTtrackD0 <<"\t\t"<<temp<<endl;
		hGenMatchedTkP->Fill(abs(TTtrackD0));
		Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
		if(tk_pixIndx>-1)
		  {
		    Double_t TTPixtrackD0 = s->GetPixTrackD0(tk_pixIndx);
		    hGenMatchedPixP->Fill(abs(TTPixtrackD0));
		    // cout<< tk_pixIndx<<"\t"<<L1PixTks_Pt->at(tk_pixIndx)<<"\t"<<L1Tks_Pt->at(tp_tkIndx)<<"\t"<<TP_Pt->at(MatchedTP)<<endl;
		  } // if ttpixtrack found
	      }// if tttrack found
	    //  cout<<"\tTP Index\t"<< MatchTP<< "\tTP ID\t"<< TP_PdgId->at(MatchTP)<<endl;
	  }//if TP found
	
      }// end loop on dau od USD



//*************************************************************************************************************************************
//*************************************************************************************************************************************//*************************************************************************************************************************************//*************************************************************************************************************************************//*************************************************************************************************************************************//*************************************************************************************************************************************//*************************************************************************************************************************************
// Working code before making methods for filling histos

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
	    TVector3 bHadvec;
	    bHadvec.SetPtEtaPhi(BPt,BEta,BPhi);
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

		    // FillSIPHistogramsDxy(DauIndx,BPhi,
		    // 			hgend0BHadDirect , hGenBHadPhiD0      , hGenBHaddeltaphiD0     , hGenBHadSIP,
		    // 			hGenMatchedTPB   ,hGenMatchedTPBPhiD0 , hGenMatchedTPBdeltaPhi , hGenMatchedTPBSIP,        
		    // 			hGenMatchedTkB   ,hGenMatchedTkBPhiD0 , hGenMatchedTkBdeltaPhi , hGenMatchedTkBSIP,
		    // 			hGenMatchedPixB  ,hGenMatchedPixBPhiD0, hGenMatchedPixBdeltaPhi, hGenMatchedPixBSIP,        
		    // 			);

		    // FillSIPHistogramsDotProd(DauIndx , bHadvec,
		    // 			     hgend0BHadDirect , hGenBHadSIP,
		    // 			     hGenMatchedTPB   , hGenMatchedTPBSIP,        
		    // 			     hGenMatchedTkB   , hGenMatchedTkBSIP,
		    // 			     hGenMatchedPixB  , hGenMatchedPixBSIP,       
		    // 			     );

		    // Double_t DauPt    = GenP_Pt->at(DauIndx);
		    // Double_t DauEta   = GenP_Eta->at(DauIndx);
		    // Double_t DauPhi   = GenP_Phi->at(DauIndx);
		    // Double_t DauX     = GenP_VertexX->at(DauIndx);//-HepMCEvt_VtxX;
		    // Double_t DauY     = GenP_VertexX->at(DauIndx);//-HepMCEvt_VtxY;
		   
		    // Double_t DauD0    = GetGenpDxy(DauIndx);//GetGenpD0(DauIndx);
		    // // Double_t DauD0Phi = GetPhiD0(DauPhi, DauD0);
		    // // Double_t deltaPhi = abs(auxTools_.DeltaPhi(DauD0Phi,BPhi));
		    // // Double_t DauAbsD0 = abs(DauD0);
		    // // Double_t Sign     = GetSignofIP(deltaPhi);
		    // // Double_t DauSIP   = DauAbsD0 * Sign;
		   
		    // TVector3 D0Vec(DauX,DauY,0);
		    // Double_t dotProd = D0Vec.Dot(bHadvec);
		    // Double_t Sign     = GetSign(dotProd);
		    // Double_t DauAbsD0 = abs(DauD0);
		    // Double_t DauSIP   = DauAbsD0 * Sign;
		    
		    // if(DauPt >2 && abs(DauEta)<2.5)
		    //   {

		    // 	hgend0BHadDirect->Fill(DauAbsD0);
		    // 	//hGenBHadPhiD0->Fill(DauD0Phi);
		    // 	//hGenBHaddeltaphiD0->Fill(deltaPhi);
		    // 	hGenBHadSIP->Fill(DauSIP);
			
		    // 	Int_t MatchedTP = GetGenMatchingTPIndex(DauIndx);
		    // 	if(MatchedTP > -1)
		    // 	  {
		    // 	    Double_t tp_D0    = GetTPDxy(MatchedTP);//GetTPD0(MatchedTP);
		    // 	    // Double_t tp_Phi   = TP_Phi->at(MatchedTP);
		    // 	    // Double_t tp_PhiD0 = GetPhiD0(tp_Phi, tp_D0);
		    // 	    // Double_t tp_dPhi  = abs(auxTools_.DeltaPhi(tp_PhiD0,BPhi));
		    // 	    // Double_t tp_AbsD0 = abs(tp_D0);
		    // 	    // Double_t tp_Sign  = GetSignofIP(tp_dPhi);
		    // 	    // Double_t tp_SIP   = tp_AbsD0 * tp_Sign;
			    
		    // 	    Double_t tp_X   = TP_POCAx->at(MatchedTP)-HepMCEvt_VtxX;
		    // 	    Double_t tp_Y   = TP_POCAy->at(MatchedTP)-HepMCEvt_VtxY;
		    // 	    TVector3 tp_D0Vec(tp_X,tp_Y,0);
		    // 	    Double_t tp_dotProd = tp_D0Vec.Dot(bHadvec);
		    // 	    Double_t tp_Sign     = GetSign(tp_dotProd);
		    // 	    Double_t tp_AbsD0 = abs(tp_D0);
		    // 	    Double_t tp_SIP   = tp_AbsD0 * tp_Sign;
			    

		    // 	    // hGenMatchedTPBPhiD0       ->Fill(tp_PhiD0);
		    // 	    // hGenMatchedTPBdeltaPhi ->Fill(tp_dPhi);
		    // 	    hGenMatchedTPBSIP         ->Fill(tp_SIP);
		    // 	    hGenMatchedTPB->Fill(tp_AbsD0);

		    // 	    Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
		    // 	    if(tp_tkIndx>-1)
		    // 	      {
		    // 		Double_t tk_D0    = GetTTTrackDxy(tp_tkIndx);//s->GetTrackD0(tp_tkIndx) * -1;
		    // 		// Double_t tk_Phi   = L1Tks_Phi->at(tp_tkIndx);
		    // 		// Double_t tk_PhiD0 = GetPhiD0(tk_Phi, tk_D0);
		    // 		// Double_t tk_dPhi  = abs(auxTools_.DeltaPhi(tk_PhiD0,BPhi));
		    // 		// Double_t tk_AbsD0 = abs(tk_D0);
		    // 		// Double_t tk_Sign  = GetSignofIP(tk_dPhi);
		    // 		// Double_t tk_SIP   = tk_AbsD0 * tk_Sign;
				
		    // 		Double_t tk_X   = L1Tks_POCAx->at(tp_tkIndx);
		    // 		Double_t tk_Y   = L1Tks_POCAy->at(tp_tkIndx);
		    // 		TVector3 tk_D0Vec(tk_X,tk_Y,0);
		    // 		Double_t tk_dotProd = tk_D0Vec.Dot(bHadvec);
		    // 		Double_t tk_Sign     = GetSign(tk_dotProd);
		    // 		Double_t tk_AbsD0 = abs(tk_D0);
		    // 		Double_t tk_SIP   = tk_AbsD0 * tk_Sign;

		    // 		// hGenMatchedTkBPhiD0       ->Fill(tk_PhiD0);
		    // 		// hGenMatchedTkBdeltaPhi ->Fill(tk_dPhi);
		    // 		hGenMatchedTkBSIP         ->Fill(tk_SIP);
		    // 		hGenMatchedTkB->Fill(tk_AbsD0);
				
		    // 		Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
		    // 		if(tk_pixIndx>-1)
		    // 		  {
		    // 		    Double_t pix_D0    = GetTTPixTrackDxy(tk_pixIndx);// s->GetPixTrackD0(tk_pixIndx);
		    // 		    // Double_t pix_Phi   = L1PixTks_Phi->at(tk_pixIndx);
		    // 		    // Double_t pix_PhiD0 = GetPhiD0(pix_Phi, pix_D0);
		    // 		    // Double_t pix_dPhi  = abs(auxTools_.DeltaPhi(pix_PhiD0,BPhi));
		    // 		    // Double_t pix_AbsD0 = abs(pix_D0);
		    // 		    // Double_t pix_Sign  = GetSignofIP(pix_dPhi);
		    // 		    // Double_t pix_SIP   = pix_AbsD0 * pix_Sign;
				    
		    // 		    Double_t pix_X   = L1PixTks_POCAx->at(tk_pixIndx);
		    // 		    Double_t pix_Y   = L1PixTks_POCAy->at(tk_pixIndx);
		    // 		    TVector3 pix_D0Vec(pix_X,pix_Y,0);
		    // 		    Double_t pix_dotProd = pix_D0Vec.Dot(bHadvec);
		    // 		    Double_t pix_Sign     = GetSign(pix_dotProd);
		    // 		    Double_t pix_AbsD0 = abs(pix_D0);
		    // 		    Double_t pix_SIP   = pix_AbsD0 * pix_Sign;

				    
		    // 		    // hGenMatchedPixBPhiD0       ->Fill(pix_PhiD0);
		    // 		    // hGenMatchedPixBdeltaPhi ->Fill(pix_dPhi);
		    // 		    hGenMatchedPixBSIP         ->Fill(pix_SIP);
		    // 		    hGenMatchedPixB->Fill(pix_AbsD0);
		    // 		  } // if ttpixtrack found
		    // 	      }// if tttrack found
		    // 	  }//if TP found
		    //   }// pt eta cut on Genp 



		  }//if charged
	      }// end loop on daughters
	    // cout<<<< iGenp<<"\t"<<Gen_ID<<endl;
	  }// if final state BHad

	if(IsLiteHadron(iGenp))
	  {
	   
	    Double_t LPt  = GenP_Pt->at(iGenp);
	    Double_t LEta = GenP_Eta->at(iGenp);
	    Double_t LPhi = GenP_Phi->at(iGenp);
	    TVector3 lHadvec;
	    lHadvec.SetPtEtaPhi(LPt,LEta,LPhi);
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
				// Double_t DauPt    = GenP_Pt->at(DIndx);
				// Double_t DauEta   = GenP_Eta->at(DIndx);
				// Double_t DauPhi   = GenP_Phi->at(DIndx);
				// Double_t DauD0    = GetGenpDxy(DIndx);//GetGenpD0(DIndx);
				// Double_t DauD0Phi = GetPhiD0(DauPhi, DauD0);
				// Double_t deltaPhi = abs(auxTools_.DeltaPhi(DauD0Phi,LPhi));
				// Double_t DauAbsD0 = abs(DauD0);
				// Double_t Sign     = GetSignofIP(deltaPhi);
				// Double_t DauSIP   = DauAbsD0 * Sign;

				Double_t DauPt    = GenP_Pt->at(DIndx);
				Double_t DauEta   = GenP_Eta->at(DIndx);
				Double_t DauPhi   = GenP_Phi->at(DIndx);
				Double_t DauX     = GenP_VertexX->at(DIndx)-HepMCEvt_VtxX;
				Double_t DauY     = GenP_VertexX->at(DIndx)-HepMCEvt_VtxY;
				
				Double_t DauD0    = GetGenpDxy(DIndx);//GetGenpD0(DauIndx);
				// Double_t DauD0Phi = GetPhiD0(DauPhi, DauD0);
				// Double_t deltaPhi = abs(auxTools_.DeltaPhi(DauD0Phi,BPhi));
				// Double_t DauAbsD0 = abs(DauD0);
				// Double_t Sign     = GetSignofIP(deltaPhi);
				// Double_t DauSIP   = DauAbsD0 * Sign;
				
				TVector3 D0Vec(DauX,DauY,0);
				Double_t dotProd = D0Vec.Dot(lHadvec);
				Double_t Sign     = GetSign(dotProd);
				Double_t DauAbsD0 = abs(DauD0);
				Double_t DauSIP   = DauAbsD0 * Sign;
				


				if(DauPt >2 && abs(DauEta)<2.5)
				  {
				    // hGenLHadPhiD0->Fill(DauD0Phi);     
				    // hGenLHaddeltaphiD0->Fill(deltaPhi);
				    hGenLHadSIP ->Fill(DauSIP);       
				    hgend0UDS->Fill(DauAbsD0);
				    
				    Int_t MatchedTP = GetGenMatchingTPIndex(DIndx);
				    if(MatchedTP > -1)
				      {
					// Double_t tp_D0    = GetTPDxy(MatchedTP);//GetTPD0(MatchedTP);
					// Double_t tp_Phi   = TP_Phi->at(MatchedTP);
					// Double_t tp_PhiD0 = GetPhiD0(tp_Phi, tp_D0);
					// Double_t tp_dPhi  = abs(auxTools_.DeltaPhi(tp_PhiD0,LPhi));
					// Double_t tp_AbsD0 = abs(tp_D0);
					// Double_t tp_Sign  = GetSignofIP(tp_dPhi);
					// Double_t tp_SIP   = tp_AbsD0 * tp_Sign;
					
					Double_t tp_D0    = GetTPDxy(MatchedTP);//GetTPD0(MatchedTP);
					// Double_t tp_Phi   = TP_Phi->at(MatchedTP);
					// Double_t tp_PhiD0 = GetPhiD0(tp_Phi, tp_D0);
					// Double_t tp_dPhi  = abs(auxTools_.DeltaPhi(tp_PhiD0,BPhi));
					// Double_t tp_AbsD0 = abs(tp_D0);
					// Double_t tp_Sign  = GetSignofIP(tp_dPhi);
					// Double_t tp_SIP   = tp_AbsD0 * tp_Sign;
					
					Double_t tp_X   = TP_POCAx->at(MatchedTP)-HepMCEvt_VtxX;
					Double_t tp_Y   = TP_POCAy->at(MatchedTP)-HepMCEvt_VtxY;
					TVector3 tp_D0Vec(tp_X,tp_Y,0);
					Double_t tp_dotProd = tp_D0Vec.Dot(lHadvec);
					Double_t tp_Sign     = GetSign(tp_dotProd);
					Double_t tp_AbsD0 = abs(tp_D0);
					Double_t tp_SIP   = tp_AbsD0 * tp_Sign;

					// hGenMatchedTPPPhiD0       ->Fill(tp_PhiD0);
					// hGenMatchedTPPdeltaPhi ->Fill(tp_dPhi);
					hGenMatchedTPPSIP         ->Fill(tp_SIP);
					hGenMatchedTPP->Fill(abs(tp_D0));
					// //*************************************************************************************
					// cout<<"*************************************************************************************"<<endl;
					// cout<<GenP_PdgId->at(DIndx)<<"\t\t"<<TP_PdgId->at(MatchedTP)<<endl;
					// cout<<DauPhi<<"\t\t"<<tp_Phi<<endl;
					// cout<<DauD0<<"\t\t"<<tp_D0<<endl;
					// cout<<DauD0Phi<<"\t\t"<<tp_PhiD0<<endl;
					// cout<< deltaPhi <<"\t\t"<<tp_dPhi<<endl;
					// //*************************************************************************************
					
					
					Int_t tp_tkIndx = TP_TTTrackIndex->at(MatchedTP);
					if(tp_tkIndx>-1)
					  {
					    // Double_t tk_D0    = GetTTTrackDxy(tp_tkIndx);//s->GetTrackD0(tp_tkIndx) * -1;
					    // Double_t tk_Phi   = L1Tks_Phi->at(tp_tkIndx);
					    // Double_t tk_PhiD0 = GetPhiD0(tk_Phi, tk_D0);
					    // Double_t tk_dPhi  = abs(auxTools_.DeltaPhi(tk_PhiD0,LPhi));
					    // Double_t tk_AbsD0 = abs(tk_D0);
					    // Double_t tk_Sign  = GetSignofIP(tk_dPhi);
					    // Double_t tk_SIP   = tk_AbsD0 * tk_Sign;
					    
					    Double_t tk_D0    = GetTTTrackDxy(tp_tkIndx);//s->GetTrackD0(tp_tkIndx) * -1;
					    // Double_t tk_Phi   = L1Tks_Phi->at(tp_tkIndx);
					    // Double_t tk_PhiD0 = GetPhiD0(tk_Phi, tk_D0);
					    // Double_t tk_dPhi  = abs(auxTools_.DeltaPhi(tk_PhiD0,BPhi));
					    // Double_t tk_AbsD0 = abs(tk_D0);
					    // Double_t tk_Sign  = GetSignofIP(tk_dPhi);
					    // Double_t tk_SIP   = tk_AbsD0 * tk_Sign;
					    
					    Double_t tk_X   = L1Tks_POCAx->at(tp_tkIndx);
					    Double_t tk_Y   = L1Tks_POCAy->at(tp_tkIndx);
					    TVector3 tk_D0Vec(tk_X,tk_Y,0);
					    Double_t tk_dotProd = tk_D0Vec.Dot(lHadvec);
					    Double_t tk_Sign     = GetSign(tk_dotProd);
					    Double_t tk_AbsD0 = abs(tk_D0);
					    Double_t tk_SIP   = tk_AbsD0 * tk_Sign;
					    
					    // hGenMatchedTkPPhiD0       ->Fill(tk_PhiD0);
					    // hGenMatchedTkPdeltaPhi ->Fill(tk_dPhi);
					    hGenMatchedTkPSIP         ->Fill(tk_SIP);
					    hGenMatchedTkP->Fill(tk_AbsD0);
					    
					    Int_t tk_pixIndx = GetPixTkIndex(tp_tkIndx);
					    if(tk_pixIndx>-1)
					      {
						// Double_t pix_D0    = GetTTPixTrackDxy(tk_pixIndx);// s->GetPixTrackD0(tk_pixIndx);
						// Double_t pix_Phi   = L1PixTks_Phi->at(tk_pixIndx);
						// Double_t pix_PhiD0 = GetPhiD0(pix_Phi, pix_D0);
						// Double_t pix_dPhi  = abs(auxTools_.DeltaPhi(pix_PhiD0,LPhi));
						// Double_t pix_AbsD0 = abs(pix_D0);
						// Double_t pix_Sign  = GetSignofIP(pix_dPhi);
						// Double_t pix_SIP   = pix_AbsD0 * pix_Sign;
						
						Double_t pix_D0    = GetTTPixTrackDxy(tk_pixIndx);// s->GetPixTrackD0(tk_pixIndx);
						// Double_t pix_Phi   = L1PixTks_Phi->at(tk_pixIndx);
						// Double_t pix_PhiD0 = GetPhiD0(pix_Phi, pix_D0);
						// Double_t pix_dPhi  = abs(auxTools_.DeltaPhi(pix_PhiD0,BPhi));
						// Double_t pix_AbsD0 = abs(pix_D0);
						// Double_t pix_Sign  = GetSignofIP(pix_dPhi);
						// Double_t pix_SIP   = pix_AbsD0 * pix_Sign;
						
						Double_t pix_X   = L1PixTks_POCAx->at(tk_pixIndx);
						Double_t pix_Y   = L1PixTks_POCAy->at(tk_pixIndx);
						TVector3 pix_D0Vec(pix_X,pix_Y,0);
						Double_t pix_dotProd = pix_D0Vec.Dot(lHadvec);
						Double_t pix_Sign     = GetSign(pix_dotProd);
						Double_t pix_AbsD0 = abs(pix_D0);
						Double_t pix_SIP   = pix_AbsD0 * pix_Sign;
						
						// hGenMatchedPixPPhiD0       ->Fill(pix_PhiD0);
						// hGenMatchedPixPdeltaPhi ->Fill(pix_dPhi);
						hGenMatchedPixPSIP         ->Fill(pix_SIP);
						hGenMatchedPixP->Fill(pix_AbsD0);
					      } // if ttpixtrack found
					  }// if tttrack found
				      }// if matched TP
				  }//ptEtaCut
			      }// if charged
			  }// if not coming from b or c
		    }// loop on daughters
		
	      }// if coming from string
	  }//if is Lite Hadron
	
	
      }// end loop on genp

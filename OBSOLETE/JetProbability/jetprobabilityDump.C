//JetProbability Plots 11 May.
    //cout<< "Event Number"<< jentry<<endl;

    //****************************************************************************************************************************************************************************
    //Prompt Case
    //****************************************************************************************************************************************************************************
    //**************************************************************************************
    // Loop on GenParticles
    //*****************************************************************.*********************
    vector<Size_t> vFinalDau;
    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	Int_t                  Gen_ID     = GenP_PdgId    ->at(iGenp);
	Int_t                  Gen_Status = GenP_Status   ->at(iGenp);
	vector<unsigned short> Gen_vDaug  = GenP_Daughters->at(iGenp);
	Double_t Gen_D0 = GetGenpD0(iGenp);

	if(abs(Gen_ID) == 11 || abs(Gen_ID) == 13) //looking for electron or muon in an event
	  {
	    if(!RecursivelyLookForMotherId(iGenp,24, false)) continue; // look is these are coming from W.
	    if(RecursivelyLookForMotherId(iGenp,5, false))   continue; // remove leptons from b...
	    if(!IsFinalGenp (iGenp, Gen_vDaug))              continue; // selecting particle only once 
	    hgend0EMu->Fill(abs(Gen_D0));
	    
	    Double_t GenEta = GenP_Eta->at(iGenp);
	    Double_t GenPhi = GenP_Phi->at(iGenp);

	    for(Size_t iTP=0; iTP<nTP; iTP++)
	      {
		Int_t tp_id = TP_PdgId->at(iTP);
		
		Double_t tp_eta = TP_Eta->at(iTP);
		Double_t tp_phi = TP_Phi->at(iTP);
		Double_t tp_x   = TP_POCAx->at(iTP);
		Double_t tp_y   = TP_POCAx->at(iTP);
		Double_t tp_d0 = -tp_x * TMath::Sin(tp_phi) + tp_y * TMath::Cos(tp_phi);
		
		if(Gen_ID == tp_id)
		  {
		    Double_t dtpgenr = auxTools_.DeltaR(GenEta,GenPhi,
							tp_eta,tp_phi);
		    if(dtpgenr < 0.0001)
		      {
			hGenMatchedTPMuE ->Fill(abs(tp_d0));
		      }// dr match

		  }// id match
	      }//loop on TP 

	    Int_t matchTk2 = MatchTrackToGenp(iGenp);
	    if(matchTk2 == -1)
	      continue;
	    Double_t matchL1Tk2_D0 = s->GetTrackD0(matchTk2); 
	    hGenMatchedTkMuE ->Fill(abs(matchL1Tk2_D0));
	    
	    
	    Int_t matchPix =GetPixTkIndex(matchTk2);
	    if(matchPix == -1)
	      continue;
	    Double_t matchPix_D0 = s->GetPixTrackD0(matchPix);
	    hGenMatchedPixMuE ->Fill(abs(matchPix_D0));
	    

	    
	  }// if ID is of lep
	
	if(abs(Gen_ID) == 1 || abs(Gen_ID) == 2 || abs(Gen_ID) == 3) //looking for electron or muon in an event
	  {
	    if(!RecursivelyLookForMotherId(iGenp,24, false)) continue; // look is these are coming from W.
	    if(RecursivelyLookForMotherId(iGenp,5, false))   continue; // remove decay from b...
	    if(RecursivelyLookForMotherId(iGenp,4, false))   continue; // remove decay from c...
	    if(!IsFinalGenp (iGenp, Gen_vDaug))              continue; // selecting particle only once 
	    //cout<<"\tParticle\t"<<iGenp<<"\t\t"<< Gen_ID <<"\t\t"<<Gen_Status<<endl;
	    
	    GetFinalDaughter(iGenp,vFinalDau);
	    
	  }
	
      }// end loop on Genp
    //   cout<<vFinalDau.size()<<endl;
    //**************************************************************************************
    
    //**************************************************************************************
    //Final Daughters of USD
    //**************************************************************************************
    for(Size_t iFD = 0; iFD<vFinalDau.size(); iFD++)
      {
	Int_t FDindx   = vFinalDau.at(iFD); 
	Int_t FDcharge = GenP_Charge->at(FDindx);
	
	if(FDcharge == 0) continue;
	if(RecursivelyLookForMotherId(FDindx,5, false))   continue; // remove decay from b...


	Int_t    FDId  = GenP_PdgId->at(FDindx);
	Double_t FDEta = GenP_Eta->at(FDindx);
	Double_t FDPhi = GenP_Phi->at(FDindx);


	for(Size_t iTP=0; iTP<nTP; iTP++)
	  {
	    Int_t tp_id = TP_PdgId->at(iTP);
	    
	    Double_t tp_eta = TP_Eta->at(iTP);
	    Double_t tp_phi = TP_Phi->at(iTP);
	    Double_t tp_x   = TP_POCAx->at(iTP);
	    Double_t tp_y   = TP_POCAx->at(iTP);
	    Double_t tp_d0 = -tp_x * TMath::Sin(tp_phi) + tp_y * TMath::Cos(tp_phi);
	    
	    if(FDId == tp_id)
	      {
		Double_t dtpgenr = auxTools_.DeltaR(FDEta,FDPhi,
						    tp_eta,tp_phi);
		if(dtpgenr < 0.0001)
		  {
		    //	    cout<<dtpgenr<< "\t"<<  tp_id<< endl;
		    hGenMatchedTPP ->Fill(abs(tp_d0));
		  }// dr match
		
	      }// id match
	  }//loop on TP 

	
	Double_t Gen_D0 = GetGenpD0(FDindx);
	hgend0UDS->Fill(abs(Gen_D0));	
	
	Int_t matchTk = MatchTrackToGenp(FDindx);
	if(matchTk == -1)
	  continue;
	Double_t matchL1Tk_D0 = s->GetTrackD0(matchTk); 
	hGenMatchedTkP ->Fill(abs(matchL1Tk_D0));
	
	
	Int_t matchPix =GetPixTkIndex(matchTk);
	if(matchPix == -1)
	  continue;
	Double_t matchPix_D0 = s->GetPixTrackD0(matchPix);
	hGenMatchedPixP ->Fill(abs(matchPix_D0));
	  
	
	
	
	
      }// end loop on FD
    //**************************************************************************************
    //****************************************************************************************************************************************************************************

    //****************************************************************************************************************************************************************************
    //B Case
    //****************************************************************************************************************************************************************************
    
    // for(Size_t temp=0; temp<nGenp; temp++)
    //   {
    // 	//cout<<EvtNumber<<endl;
    // 	// if(EvtNumber==91144 ||EvtNumber==12507 || EvtNumber==14881 || EvtNumber==41881 || EvtNumber==44876 || EvtNumber==76669 ||EvtNumber==69353)
    // 	//   myPrintGenp(temp,true);
    // 	Int_t id = GenP_PdgId->at(temp);
    // 	Int_t mod = id%10000;
    // 	Int_t bHad=0;
    //   if(mod>1000)
    // 	{
    // 	  bHad = mod/1000;
    // 	}
    //   else if(mod<1000)
    // 	{
    // 	  bHad = mod/100;
    // 	}
    //   if(abs(bHad) != 5)
    // 	continue;
 
    //   //  myPrintGenp(temp,true);
    //   }// end temp loop on Genp
    
    //**************************************************************************************
    // Loop on GenParticles
    //**************************************************************************************
   
    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	//	if(jentry == 12236) myPrintGenp(iGenp,true);
	Int_t                  Gen_ID     = GenP_PdgId    ->at(iGenp);
	Int_t                  Gen_Status = GenP_Status   ->at(iGenp);
	vector<unsigned short> Gen_vDaug  = GenP_Daughters->at(iGenp);
	Double_t Gen_D0  = GetGenpD0(iGenp);
	Double_t Gen_Pt  = GenP_Pt->at(iGenp);
	Double_t Gen_Eta = GenP_Eta->at(iGenp);
	Double_t Gen_Phi = GenP_Phi->at(iGenp);
	
	if(Gen_Pt < 2.0) continue;
	if(Gen_Eta > 2.5) continue;
	
	if(abs(Gen_ID) == 211 || abs(Gen_ID) == 321 || abs(Gen_ID) == 13 ) //looking  muon/pion/Kaon in an event
	  {
	    if (!IsBHadronProduct(iGenp)) continue;
	    if (IsCHadronProduct(iGenp)) continue;
	    if(!IsFinalGenp (iGenp, Gen_vDaug))              continue; // selecting particle only once 
	    //cout<<"\tParticle\t"<<iGenp<<"\t\t"<< Gen_ID <<"\t\t"<<Gen_Status<<endl;
	    Double_t Gen_D0 = GetGenpD0(iGenp);
	    hgend0MuKaPiBHad->Fill(abs(Gen_D0));	


	       for(Size_t iTP=0; iTP<nTP; iTP++)
	      {
	    	Int_t tp_id = TP_PdgId->at(iTP);
		
	    	Double_t tp_eta = TP_Eta->at(iTP);
	    	Double_t tp_phi = TP_Phi->at(iTP);
	    	Double_t tp_x   = TP_POCAx->at(iTP);
	    	Double_t tp_y   = TP_POCAx->at(iTP);
	    	Double_t tp_d0 = -tp_x * TMath::Sin(tp_phi) + tp_y * TMath::Cos(tp_phi);
		
	    	if(Gen_ID == tp_id)
	    	  {
	    	    Double_t dtpgenr = auxTools_.DeltaR(Gen_Eta,Gen_Phi,
	    						tp_eta,tp_phi);
	    	    if(dtpgenr < 0.0001)
	    	      {
	    		//cout<<dtpgenr<< "\t"<<  tp_id<< endl;
	    		hGenMatchedTPB ->Fill(abs(tp_d0));
	    	      }// dr match
		    
	    	  }// id match
	      }//loop on TP 


	    Int_t matchTk = MatchTrackToGenp(iGenp);
	    if(matchTk == -1)
	      continue;
	    Double_t matchL1Tk_D0 = s->GetTrackD0(matchTk); 
	    hGenMatchedTkB ->Fill(abs(matchL1Tk_D0));
	    
	    
	    Int_t matchPix =GetPixTkIndex(matchTk);
	    if(matchPix == -1)
	      continue;
	    Double_t matchPix_D0 = s->GetPixTrackD0(matchPix);
	    hGenMatchedPixB ->Fill(abs(matchPix_D0));
	       

	    
	  }// end if  pion/kaon/muon
      }// end loop on Gen


    //**************************************************************************************
    // Gen Jets SIP
    //**************************************************************************************
    
    //*************************************************
    // Selecting bHadrons in an event
    //*************************************************
    vector <Size_t> bquark;
    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	Int_t Id = GenP_PdgId->at(iGenp);
	Int_t mod = Id%10000;
	Int_t bHad = 0;
	if(mod>1000)  bHad = mod/1000;
	else if(mod<1000)  bHad = mod/100;
	

	if(abs(bHad) != 5) continue;
	
	Bool_t bIsHadDeXi=false;
	vector<unsigned short> DauVec = GenP_Daughters->at(iGenp);
	for(Size_t idau=0; idau<DauVec.size(); idau++)
	  {
	    Size_t dauIndx = DauVec.at(idau);
	    Int_t  dauId  = GenP_PdgId->at(dauIndx);
	    Int_t  dauMod = dauId%10000;
	    Int_t  dauBHad= 0;
	    if(dauMod>1000)  dauBHad = dauMod/1000;
	    else if(dauMod<1000)  dauBHad = dauMod/100;
	 
	    if(dauBHad == bHad){ bIsHadDeXi=true;}
	  }
	if(bIsHadDeXi) continue;
	//	cout<<iGenp<<"\tThis Hadron\t"<<Id << endl<<endl; 
	// if(abs(GenP_PdgId->at(iGenp)) != 5)
	//   continue;
	// vector<unsigned short> DauVec = GenP_Daughters->at(iGenp);
	// if(!IsFinalGenp(iGenp,DauVec))
	//   continue;
	bquark.push_back(iGenp);
	
      }// end loop to collect bquarks
    //*************************************************

    //    if(jentry == 12236) cout<<bquark.size()<<endl;

    //*************************************************
    // Selecting Looking for leading pT dau of bHadrons
    //*************************************************
    vector <Size_t> ldgDau;
    for(Size_t iB=0; iB<bquark.size(); iB++)
      {
	Size_t bIndx= bquark.at(iB);
	Int_t  bId = GenP_PdgId->at(bIndx);

	Int_t MaxPtIndx=-1;
	Double_t MaxPt =0;
	//if(jentry == 12236)cout<<"\tindx\t"<<bIndx<<"\t"<< bId<<endl;
	for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
	  {
	    //	    if(jentry == 12236) cout<< "balay bhai ballay "<<endl;
	    Int_t Gen_Status = GenP_Status ->at(iGenp);
	    Int_t Gen_charge = GenP_Charge ->at(iGenp);
	    if(Gen_Status != 1) continue;
	    if(Gen_charge == 0) continue;
	    if(RecursivelyLookForMotherId(iGenp, bId, true))
	      {
		Int_t Gen_Id = GenP_PdgId->at(iGenp);
		Double_t Gen_Pt = GenP_Pt->at(iGenp);
		//	if(jentry == 12236)cout<<"\tindx\t"<<iGenp<<"\t"<< Gen_Id<<"\t\t"<<Gen_Pt<<"\t\t"<<MaxPt<<"\t"<<MaxPtIndx <<endl;
		
		if(Gen_Pt > MaxPt)
		  {
		    MaxPt = Gen_Pt;
		    MaxPtIndx = iGenp;
		  }
	      }
	  }// emmd loop on Genp

	//	if(jentry == 12236) cout<< "This is fucking index    "<<MaxPtIndx<<endl;
	if(MaxPtIndx <0)
	  {
	    cout<< "None"<<jentry<<endl;
	    continue;
	  }
	ldgDau.push_back(MaxPtIndx);
	//	if(jentry == 12236)cout <<ldgDau.size()<<endl;
      }// end loop on bquarks
    //*************************************************

    
    //    cout<< ldgDau.size()<<endl;
    //*************************************************
    // Making cone aroung leading daughter and selecting particles
    //*************************************************
    for(Size_t temp=0; temp<ldgDau.size(); temp++)
      {
	Int_t ldgIndx = ldgDau.at(temp);
	
	Double_t ldgEta =  GenP_Eta->at(ldgIndx);
	Double_t ldgPhi =  GenP_Phi->at(ldgIndx);
	//	cout<< "Separatot"<<endl;
	vector<Size_t> bGenJet;
	for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
	  {
	    if(GenP_Status->at(iGenp)!=1) continue;
	    if(GenP_Charge->at(iGenp)==0) continue;
	    Double_t eta = GenP_Eta->at(iGenp);
	    Double_t phi = GenP_Phi->at(iGenp);
	    Double_t dR = auxTools_.DeltaR(ldgEta,ldgPhi,
					   eta, phi);
	    if(dR<0.3)
	      {
		bGenJet.push_back(iGenp);
		//	cout<<iGenp<<"\t" <<GenP_PdgId->at(iGenp)<<endl;
	      }

	  }// end loop on Gen
	//*************************************************
	
	//*************************************************
	// Calculating SIP
	//*************************************************
	Double_t JetPhi=0;
	Double_t JetPhiNum=0;
	Double_t JetPhiDen=0;
	for (Size_t iGJ=0; iGJ< bGenJet.size(); iGJ++)
	  {
	    Int_t GenIndx = bGenJet.at(iGJ);
	    Double_t GenPhi = GenP_Phi->at(GenIndx);
	    Double_t GenPt  = GenP_Pt->at(GenIndx);
	    JetPhiNum += GenPhi*GenPt;
	    JetPhiDen += GenPt;
	  }// end loop on GenJet
	JetPhi=JetPhiNum/JetPhiDen;
	//	cout<<endl<<endl;
	for (Size_t iGJ=0; iGJ< bGenJet.size(); iGJ++)
	  {
	    Int_t GenIndx = bGenJet.at(iGJ);
	    Double_t GenPhi = GenP_Phi->at(GenIndx);
	    Double_t DeltaPhi = auxTools_.DeltaPhi(GenPhi,JetPhi) +PI/2;
	    Double_t Sign = 0;
	    if (DeltaPhi<PI/2) Sign=1.0;
	    else if (DeltaPhi>PI/2) Sign=-1.0;
	    
	    Double_t GenD0    = GetGenpD0(GenIndx);
	    Double_t GenAbsD0 = abs(GenD0);
	    Double_t GenSIP   = GenAbsD0 * Sign; 
	    hGenSIP -> Fill(GenSIP); 
	    
	    //	    cout<<"\t"<<GenD0<<"\t"<<GenAbsD0<<"\t"<<GenSIP<<endl;
	    
	  }

	//	cout<<ldgDau.at(temp)<<"\t"<< GenP_Pt->at(ldgDau.at(temp))<<endl;
      }
    //**************************************************************************************


    //**************************************************************************************
    // W/O PU test
    //**************************************************************************************
  
    //*************************************************
    // Selecting bHadrons in an event
    //*************************************************
    vector <Size_t> bList;
    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	Int_t Id = GenP_PdgId->at(iGenp);
	Int_t mod = Id%10000;
	Int_t bHad = 0;
	if(mod>1000)  bHad = mod/1000;
	else if(mod<1000)  bHad = mod/100;
	

	if(abs(bHad) != 5) continue;
	
	Bool_t bIsHadDeXi=false;
	vector<unsigned short> DauVec = GenP_Daughters->at(iGenp);
	for(Size_t idau=0; idau<DauVec.size(); idau++)
	  {
	    Size_t dauIndx = DauVec.at(idau);
	    Int_t  dauId  = GenP_PdgId->at(dauIndx);
	    Int_t  dauMod = dauId%10000;
	    Int_t  dauBHad= 0;
	    if(dauMod>1000)  dauBHad = dauMod/1000;
	    else if(dauMod<1000)  dauBHad = dauMod/100;
	 
	    if(dauBHad == bHad){ bIsHadDeXi=true;}
	  }
	if(bIsHadDeXi) continue;

	bList.push_back(iGenp);
	
      }// end loop to collect bquarks
    //*************************************************
    
    //*************************************************
    // Selecting  dau of bHadrons
    //*************************************************
    
    for(Size_t iB=0; iB<bList.size(); iB++)
      {
	Size_t bIndx= bList.at(iB);
	Int_t  bId = GenP_PdgId->at(bIndx);
	vector <Size_t> dauNoPU;

	for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
	  {
	    //	    if(jentry == 12236) cout<< "balay bhai ballay "<<endl;
	    Int_t Gen_Status = GenP_Status ->at(iGenp);
	    Int_t Gen_charge = GenP_Charge ->at(iGenp);
	    if(Gen_Status != 1) continue;
	    if(Gen_charge == 0) continue;
	    if(RecursivelyLookForMotherId(iGenp, bId, true))
	      {
		dauNoPU.push_back(iGenp);
	      }
	  }// emmd loop on Genp
	//	cout<<dauNoPU.size()<<endl;

	Double_t bHad_phi = GenP_Phi->at(bIndx);
	//	cout<<bIndx<<"\t"<<GenP_PdgId->at(bIndx)<<"\t"<< GenP_Phi->at(bIndx)<<endl;
	
	for(Size_t dau=0; dau<dauNoPU.size(); dau++)
	  {
	    Size_t indexnoPU = dauNoPU.at(dau);
	    Double_t genPhi = GenP_Phi->at(indexnoPU);
	    Double_t DeltaPhi = auxTools_.DeltaPhi(genPhi,bHad_phi) +PI/2;
	    Double_t Sign = 0;
	    if (DeltaPhi<PI/2) Sign=1.0;
	    else if (DeltaPhi>PI/2) Sign=-1.0;


	    Double_t genD0    = GetGenpD0(indexnoPU);
	    Double_t genAbsD0 = abs(genD0);
	    Double_t genSIP   = genAbsD0 * Sign; 
	    hGenNoPUSIP -> Fill(genSIP); 

	    //	    cout<< indexnoPU<<"\t\t"<< GenP_PdgId->at(indexnoPU)<<endl;
	  }
      }// end loop on bquarks
    //*************************************************



    //**************************************************************************************

  hGenMatchedTkB->Scale(1/hGenMatchedTkB->Integral());
  hGenMatchedTkP->Scale(1/hGenMatchedTkP->Integral());

  hGenMatchedTkB->GetYaxis()->SetRangeUser(1E-7,1);
  hGenMatchedTkP->GetYaxis()->SetRangeUser(1E-7,1);

  hGenMatchedPixB->Scale(1/hGenMatchedPixB->Integral());
  hGenMatchedPixP->Scale(1/hGenMatchedPixP->Integral());

  hGenMatchedPixB->GetYaxis()->SetRangeUser(1E-7,1);
  hGenMatchedPixP->GetYaxis()->SetRangeUser(1E-7,1);

  hgend0MuKaPiBHad->Scale(1/hgend0MuKaPiBHad->Integral());
  hgend0EMu->Scale(1/hgend0EMu->Integral());
  hgend0UDS->Scale(1/hgend0EMu->Integral());

  hgend0MuKaPiBHad->GetYaxis()->SetRangeUser(1E-7,1);
  hgend0EMu->Scale(1/hgend0EMu->Integral());
  hgend0UDS->Scale(1/hgend0EMu->Integral());

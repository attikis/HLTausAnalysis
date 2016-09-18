	// //	else if (!IsJetInBDirection(iJet) && !IsJetInCDirection(iJet))
	//   {
	//     for(Size_t ijet_tk=0; ijet_tk<jet_ntk; ijet_tk++)
	//       {
	// 	Int_t tkIndex = L1TkJet_TTTrackIndex->at(iJet).at(ijet_tk);
	// 	Double_t tk_d0    = s->GetTrackD0(tkIndex);//-tk_pocax*TMath::Sin(tk_phi) + tk_pocay*TMath::Cos(tk_phi);
	// 	htk_absd0P-> Fill(abs(tk_d0));
	//       }
	//   }


	// if(IsJetInBDirection(iJet))
	//   {
	//     for(Size_t ijet_tk=0; ijet_tk<jet_ntk; ijet_tk++)
	//       {
	// 	Int_t tkIndex      = L1TkJet_TTTrackIndex->at(iJet).at(ijet_tk);
	// 	Double_t tk_d0     = s->GetTrackD0(tkIndex);
	// 	htk_absd0H->Fill(abs(tk_d0));
	//       }
	//     //	    cout<<iJet<< "\tis in HQ dir \t"<<endl;
	//   }

	// if(IsJetInUDirection(iJet) || IsJetInDDirection(iJet) )
	//   {
	//     for(Size_t ijet_tk=0; ijet_tk<jet_ntk; ijet_tk++)
	//       {
	// 	Int_t tkIndex      = L1TkJet_TTTrackIndex->at(iJet).at(ijet_tk);
	// 	Double_t tk_d0     = s->GetTrackD0(tkIndex);
	// 	htk_absd0P->Fill(abs(tk_d0));
	//       }
	//     //cout<<iJet<< "\tis in lQ dir \t"<<endl;
	//   }

	// if(!IsJetHeavyFlavour(iJet))
	//   {
	//     for(Size_t ijet_tk=0; ijet_tk<jet_ntk; ijet_tk++)
	//       {
	// 	Int_t tkIndex      = L1TkJet_TTTrackIndex->at(iJet).at(ijet_tk);
	// 	Double_t tk_d0     = s->GetTrackD0(tkIndex);
	// 	htk_absd0P->Fill(abs(tk_d0));
	//       }
	//   }

	

	// for(Size_t ijet_pix=0; ijet_pix<njet_pix; ijet_pix++)
	//   {
	//     Size_t iPix = jet_pixTks.at(ijet_pix);
	    
	//     Double_t pix_phi     = L1PixTks_Phi->at(iPix);
	//     Double_t pix_d0      = s->GetPixTrackD0(iPix);
	//     Double_t pix_d0Sig   = s->GetPixTrackD0Sig(iPix);
	//     Double_t pix_dphi    = DeltaPhid0(pix_phi, jet_phi);
	//     Int_t    sign        = Signd0(pix_dphi);
	//     Double_t pix_signedd0    = abs(pix_d0)*sign;
	//     Double_t pix_signedd0Sig = abs(pix_d0Sig)*sign;
	    
	//     if ( IsWithinEtaRegion("Central", jet_eta) )
	//       {
	// 	hsign_C             ->Fill (sign);
	// 	hBGen_SignedD0_C    ->Fill(pix_signedd0);
	// 	hBGen_DeltaPhi_C    ->Fill(pix_dphi);
	// 	hBGen_SignedD0Sig_C ->Fill(pix_signedd0Sig);
	//       }
	    
	//     else if ( IsWithinEtaRegion("Intermediate", jet_eta) )
	//       {
	// 	hsign_I             ->Fill (sign);
	// 	hBGen_SignedD0_I    ->Fill(pix_signedd0);
	// 	hBGen_DeltaPhi_I    ->Fill(pix_dphi);
	// 	hBGen_SignedD0Sig_I ->Fill(pix_signedd0Sig);
	//       }
	//     else if ( IsWithinEtaRegion("Forward", jet_eta) )
	//       {
	// 	hsign_F             ->Fill (sign);
	// 	hBGen_SignedD0_F    ->Fill(pix_signedd0);
	// 	hBGen_DeltaPhi_F    ->Fill(pix_dphi);
	// 	hBGen_SignedD0Sig_F ->Fill(pix_signedd0Sig);
	//       }

	//     hsign             ->Fill (sign);
	//     hBGen_SignedD0    ->Fill(pix_signedd0);
	//     hBGen_DeltaPhi    ->Fill(pix_dphi);
	//     hBGen_SignedD0Sig ->Fill(pix_signedd0Sig);
  
	//   }// loop pix tks in jets


  // for(Size_t itk=0; itk< L1Tks_Pt->size(); itk++)
    //   {
    // 	Double_t tk_d0 = s->GetTrackD0(itk);
    // 	hAlltk_d0->Fill(tk_d0);
    //   }
    // for(Size_t ipix=0; ipix< L1PixTks_Pt->size(); ipix++)
    //   {
    // 	Double_t pix_d0 = s->GetPixTrackD0(ipix);
    // 	hAllpix_d0->Fill(pix_d0);
    //   }



//****************************************************************************************************************************************************************
// This is 9th may Code dump....
//****************************************************************************************************************************************************************

   //    cout<<"=========================================================================================================================\n"<<jentry<<endl<<endl;

    // for(Size_t iTP=0; iTP<TP_Pt->size(); iTP++)
    //   {
    // 	Int_t TP_TKID = TP_TTTrackIndex->at(iTP);
    // 	if(TP_TKID == -1) continue;
    // 	//	cout<<"...................."<<TP_TKID<<endl;
    //   }
    //cout<<"\tEvent Number\t"<<jentry<<endl;
    for(Size_t iGenp=0; iGenp<GenP_PdgId->size(); iGenp++)
      {
    	//if(jentry == 847) myPrintGenp(iGenp, true);
    	Int_t    genID     = GenP_PdgId->at(iGenp);
    	Int_t    genStatus = GenP_Status->at(iGenp); 
    	Int_t    genCharge = GenP_Charge->at(iGenp); 
    	Double_t genVtxX   = GenP_VertexX->at(iGenp);
    	Double_t genVtxY   = GenP_VertexY->at(iGenp);
    	Double_t genVtxZ   = GenP_VertexZ->at(iGenp);
    	Double_t genPHI    = GenP_Phi->at(iGenp);
	Double_t genD0     = -genVtxX*TMath::Sin(genPHI) + genVtxY*TMath::Cos(genPHI);

    	//	Int_t MotIndx = PosOfMotherId(iGenp, 5,false);
    	//if(!RecursivelyLookForMotherId(iGenp,5, false)) continue;
    	//if(RecursivelyLookForMotherId(iGenp,4, false)) continue;
    	//if(genStatus != 1)  continue;
    	//if(genCharge == 0) continue;
    	//	cout<<"\tGenpIndex of last\t"<<iGenp<<endl;


    	//**********************************************************************************************
    	// Correction on D0 from vertex
    	//**********************************************************************************************
    	// Int_t MotIndx = PosOfMotherId(iGenp, 5,false);
    	// //cout<<MotIndx<<endl;
    	// if (MotIndx == 65535)continue;
    	// //cout<<MotIndx<<endl;
    	// Double_t MotVtxX   = GenP_VertexX->at(MotIndx);
    	// Double_t MotVtxY   = GenP_VertexY->at(MotIndx);
    	// Double_t MotVtxZ   = GenP_VertexZ->at(MotIndx);
	
    	// // //	Double_t genD0     = -genVtxX*TMath::Sin(genPHI) + genVtxY*TMath::Cos(genPHI); 
    	// Double_t genD0     = -(genVtxX-MotVtxX)*TMath::Sin(genPHI) + (genVtxY-MotVtxY)*TMath::Cos(genPHI); 
    	// hgend0B->Fill(abs(genD0));
    	//**********************************************************************************************



      	//	myPrintGenp(iGenp,true);
    	//if(abs(GenP_PdgId->at(iGenp)) == 1 ||abs(GenP_PdgId->at(iGenp)) == 2 ) cout<<GenP_PdgId->at(iGenp)<<endl; 
    	//if(abs(GenP_PdgId->at(iGenp)) != 92 ) continue;              // selecting b quark from event
    	if(genStatus != 1) 
    	  continue;
    	if(genCharge == 0) 
    	  continue;
    	//cout<< endl<<endl;
    	// if(!RecursivelyLookForMotherId(iGenp,92, false)) 
    	//   continue;
    	if(IsBHadronProduct(iGenp))
    	  {
    	    if(!IsCHadronProduct(iGenp))
    	      {
    		hgend0BHad->Fill(abs(genD0));
    	      }
    	    //	    cout<<"\t\t"<<iGenp<<endl;
    	  }
    	else// if(IsUHadronProduct(iGenp) || IsDHadronProduct(iGenp))// || IsSHadronProduct(iGenp))
    	  {
    	    if(!IsBHadronProduct(iGenp))
    	      {
    		if(!IsCHadronProduct(iGenp))
    		  {
    		    hgend0PHad->Fill(abs(genD0));
    		  }
    	      }
    	  }
	
    	if(RecursivelyLookForMotherId(iGenp,5, false)) 
    	  {
	
    	    if(!RecursivelyLookForMotherId(iGenp,4, false)) 
    	      {
    		//continue;
    		hgend0B->Fill(abs(genD0));
    	      }
    	  }
    	if(RecursivelyLookForMotherId(iGenp,92, false)) 
    	  {
    	    if(!RecursivelyLookForMotherId(iGenp,5, false)) 
    	      {
    		if(!RecursivelyLookForMotherId(iGenp,4, false)) 
    		  {
    		    hgend0P->Fill(abs(genD0));
    		  }
    	      }
    	  }
	
      }//end loop on Genp



    //****************************************************************************************************
    // Matching Jets with b quark....
    //****************************************************************************************************
    
    vector <Size_t> SelJets;
    vector <Size_t> UsedJets;
    for(Size_t iGenp=0; iGenp<GenP_PdgId->size(); iGenp++)
      {
    	if(abs(GenP_PdgId->at(iGenp)) != 5) continue;              // selecting b quark from event
    	vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
    	if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that b is used only once
	
    	Double_t bEta = GenP_Eta->at(iGenp);
    	Double_t bPhi = GenP_Phi->at(iGenp);
    	Double_t tempdr = 1000;
    	Size_t  jetIndx =-1;
    	for(Size_t iJet=0; iJet<L1TkJet_Pt->size(); iJet++)
    	  {
    	    Bool_t bUsedjet=false;
    	    for(Size_t iused=0; iused<UsedJets.size(); iused++)
    	      {
    		if(UsedJets.at(iused) == iJet)
    		  bUsedjet = true;
    	      }
    	    if(bUsedjet)
    	      continue;


    	    Size_t   jnTks = L1TkJet_TTTrackIndex->at(iJet).size();
    	    if(jnTks<2) continue;

    	    Double_t jEta = L1TkJet_Eta->at(iJet);
    	    Double_t jPhi = L1TkJet_Phi->at(iJet);
    	    Double_t dR = auxTools_.DeltaR(bEta,bPhi,
    					   jEta,jPhi);
    	    if(dR<1)
    	      if(dR < tempdr)
    		{
    		  tempdr = dR;
    		  jetIndx = iJet;
    		}
    	  }// loop on Jets
    	if(jetIndx != -1)
    	  {
    	    SelJets.push_back(jetIndx);
    	    UsedJets.push_back(jetIndx);
    	  }
    	//	cout<< iGenp << endl <<GenP_PdgId->at(iGenp)<<endl << jetIndx<<endl << tempdr<<endl<<endl<<endl;

      }//loop on Genp


    // cout<< "size" <<SelJets.size()<<endl;
    for(Size_t ibjet=0; ibjet<SelJets.size(); ibjet++)
      {
    	Size_t JetINDX = SelJets.at(ibjet);
    	Int_t  JetNTKS = L1TkJet_TTTrackIndex->at(JetINDX).size();
    	for(Size_t ijetTk=0; ijetTk<JetNTKS; ijetTk++)
    	  {
    	    Size_t TkINDX = L1TkJet_TTTrackIndex->at(JetINDX).at(ijetTk);
    	    Double_t TkD0 = s->GetTrackD0(TkINDX);
    	    hBjetTkabsd0 ->Fill(abs(TkD0));
    	  }
    	//	cout<<SelJets.at(a)<<endl<<endl<<endl;
      }
    //****************************************************************************************************
    
    
    //****************************************************************************************************
    // Selecting Jets away from b quark....
    //****************************************************************************************************

    //    cout<<"_____________________________________________________________\n";
    for(Size_t iJet=0; iJet<L1TkJet_Pt->size(); iJet++)
      {
    	Size_t   jnTks = L1TkJet_TTTrackIndex->at(iJet).size();
    	if(jnTks<2) continue;
	
    	Double_t jEta = L1TkJet_Eta->at(iJet);
    	Double_t jPhi = L1TkJet_Phi->at(iJet);
    	Bool_t bisNearB=false;
    	//	cout<< "JetSeparation \n";
    	for(Size_t iGenp=0; iGenp<GenP_PdgId->size(); iGenp++)
    	  {
    	    if(abs(GenP_PdgId->at(iGenp)) != 5) continue;              // selecting b quark from event
    	    vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
    	    if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that b is used only once
	    
    	    Double_t bEta = GenP_Eta->at(iGenp);
    	    Double_t bPhi = GenP_Phi->at(iGenp);
    	    Double_t dR = auxTools_.DeltaR(bEta,bPhi,
    					   jEta,jPhi);
    	    //  cout<<dR<<endl;
    	    if(dR < 2)
    	      bisNearB=true;

    	  }//end loop on genp
    	if(bisNearB)continue;
    	//	cout<<"    Far Tracks      "<< iJet<<endl;
    	Size_t nTksP = L1TkJet_TTTrackIndex->at(iJet).size();
    	for(Size_t itk =0; itk<nTksP; itk++)
    	  {
    	    Size_t tkIndex= L1TkJet_TTTrackIndex->at(iJet).at(itk);
    	    Double_t tkd0 = s->GetTrackD0(tkIndex);
    	    hPjetTkabsd0 ->Fill(abs(tkd0));
    	  }

      }//end loop on Jets
    //****************************************************************************************************



    // 	// if(!RecursivelyLookForMotherId(iGenp,1, false)) 
    // 	//   continue;


    // 	//	cout<<iGenp<<"\t\t"<< genID<<"\t \t"<<genStatus<<endl;

 
    // //	if(RecursivelyLookForMotherId(iGenp,5, false)) continue;
    // 	  // {
    // 	  //   cout<< "yrah b ka hay"<<endl;
    // 	  // }
    // //	else if(RecursivelyLookForMotherId(iGenp,4, false))continue;
    // 	  // {
    // 	  //   cout<< "yrah c ka hay"<<endl;
    // 	  // }
    // 	// cout<<iGenp << "       yeah b k beghair hay"<<endl;
    // 	// vector<unsigned short> DauVec = GenP_Daughters->at(iGenp);
    // 	// for(Size_t tempDau=0; tempDau<DauVec.size(); tempDau++)
    // 	//   {
    // 	//     Int_t DauIndx = DauVec.at(tempDau);
    // 	//     Int_t DauId = GenP_PdgId->at(DauIndx);
    // 	//     Int_t DauStatus = GenP_Status->at(DauIndx);
    // 	//     Int_t DauCharge = GenP_Charge->at(DauIndx);
    // 	//     //if(DauCharge == 0) continue;
	  
    // 	//   }
    
  	
  

    // // vector<Size_t> bquark;
    // // for(Size_t iGenp=0; iGenp<GenP_PdgId->size(); iGenp++)
    // //   {
    // // 	//	myPrintGenp(iGenp,true);
    // // 	//if(abs(GenP_PdgId->at(iGenp)) == 1 ||abs(GenP_PdgId->at(iGenp)) == 2 ) cout<<GenP_PdgId->at(iGenp)<<endl; 
    // // 	if(abs(GenP_PdgId->at(iGenp)) != 5 ) continue;              // selecting b quark from event
    // // 	vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
    // // 	if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that b is used only once	  
    // // 	bquark.push_back(iGenp);
    
    // //   }

    // // vector<Size_t> cquark;
    // // for(Size_t iGenp=0; iGenp<GenP_PdgId->size(); iGenp++)
    // //   {
    // // 	//if(abs(GenP_PdgId->at(iGenp)) == 1 ||abs(GenP_PdgId->at(iGenp)) == 2 ) cout<<GenP_PdgId->at(iGenp)<<endl; 
    // // 	if(abs(GenP_PdgId->at(iGenp)) != 4 ) continue;              // selecting b quark from event
    // // 	vector<unsigned short> DauVec = GenP_Daughters->at(iGenp); // getting vector of daughters
    // // 	if(!IsFinalGenp(iGenp, DauVec)) continue;                  // making sure that b is used only once	  
    // // 	cquark.push_back(iGenp);
    
    // //   }
    // // //    cout<<endl<<endl<<endl<<"\t\tNumber of B's\t\t"<<bquark.size()<<endl; 
				       
    Size_t nJet = L1TkJet_Pt->size();
    
     hjet_multi->Fill(nJet);
    
    // //    cout<<endl<<endl;
    
    // // //**********************************************
    // // // Pion matching jet
    // // //**********************************************
    // // Double_t piEta = GenP_Eta->at(0);
    // // Double_t piPhi = GenP_Phi->at(0);
    // // Double_t pidr = 1000;
    // // Size_t  pijetIndx =-1;
    // // for(Size_t ij=0; ij < nJet; ij++)
    // //   {
    // // 	Double_t jet_eta = L1TkJet_Eta    ->at(ij);
    // // 	Double_t jet_phi = L1TkJet_Phi    ->at(ij);
    // // 	Double_t tempdR = auxTools_.DeltaR(piEta,piPhi,
    // // 					   jet_eta,jet_phi);
    // // 	//if(tempdR<1)
    // // 	if(tempdR < pidr)
    // // 	  {
    // // 	    pidr = tempdR;
    // // 	    pijetIndx = ij;
    // // 	  }
    // //   }
    // // if(pijetIndx == -1) continue;
    // // //**********************************************
    
    // //    cout<<"_________________________Event Separation____________________________________"<<jentry<<endl;
    for(Size_t iJet=0; iJet < nJet; iJet++)
      {
    	//	if(iJet!=pijetIndx) continue;
    	Double_t jet_pt  = L1TkJet_Pt     ->at(iJet);
    	Double_t jet_eta = L1TkJet_Eta    ->at(iJet);
    	Double_t jet_phi = L1TkJet_Phi    ->at(iJet);
    	Double_t jet_e   = L1TkJet_E      ->at(iJet);
    	Double_t jet_vtx = L1TkJet_Vertex ->at(iJet);
    	Size_t   jet_ntk = L1TkJet_TTTrackIndex->at(iJet).size();
    	vector<Int_t> jet_pixTks;
    	GetPixTkVector(iJet,jet_pixTks);
    	Size_t njet_pix = jet_pixTks.size();
    	if(jet_ntk<2) continue;

	


    	//***********************************************************************************************
    	// Validation of L1TkJetObjects
    	//***********************************************************************************************
    	//	if(jet_ntk>0) continue;
	
    	hjet_pt->Fill(jet_pt);
    	hjet_eta->Fill(jet_eta);
    	hjet_phi->Fill(jet_phi);
    	hjet_e->Fill(jet_e);
    	hjet_vtx->Fill(jet_vtx);
    	hjet_ntk->Fill(jet_ntk);
	
    	//	cout<<"\n\n\t Jet Separation\t NTKS  \n________________________________________________________________________________________    "<<jet_ntk<<endl;
	
    	// Bool_t bHF=false;
    	// for(Size_t temp=0; temp<bquark.size(); temp++)
    	//   {
    	//     Int_t genIndx = bquark.at(temp);
    	//     Int_t genId = GenP_PdgId->at(genIndx);
    	//     //  cout<< genId <<endl;
    	//     Double_t geneta=GenP_Eta->at(genIndx);
    	//     Double_t genphi=GenP_Phi->at(genIndx);
    	//     Double_t deltaRr = auxTools_.DeltaR(geneta,genphi,
    	// 					jet_eta,jet_phi);
    	//     if(deltaRr<1.5)
    	//       {
    	// 	bHF=true;
    	// 	break;
    	//       }
    	//   }
    	// for(Size_t temp=0; temp<cquark.size(); temp++)
    	//   {
    	//     Int_t genIndx = cquark.at(temp);
    	//     Int_t genId = GenP_PdgId->at(genIndx);
    	//     //cout<< genId <<endl;
    	//     Double_t geneta=GenP_Eta->at(genIndx);
    	//     Double_t genphi=GenP_Phi->at(genIndx);
    	//     Double_t deltaRr = auxTools_.DeltaR(geneta,genphi,
    	// 					jet_eta,jet_phi);
    	//     if(deltaRr<1.5)
    	//       {
    	// 	bHF=true;
    	// 	break;
    	//       }
    	//   }
	
	
    	Bool_t bHFJet;
    	Bool_t bIsBJet=false;
    	Bool_t bIsCJet=false;
	
    	Bool_t testtttt = IsJetInBDirection(iJet);
    	//	if(IsJetInBDirection(iJet))
    	{
    	  // if(!IsJetInCDirection(iJet))
    	  {
    	    bHFJet=true;
    	    bIsBJet= true;
    	  }
    	}
    	//	else if(IsJetInCDirection(iJet))
    	{
    	  bHFJet=true;
    	  bIsCJet = true;
    	}
    	//else 
    	bHFJet=false;
	
	
    	// if(IsJetInDDirection(iJet))
    	//   cout<<"aho    D"<<endl;
    	// if(IsJetInUDirection(iJet))
    	//   cout<<"aho    U"<<endl;
	
    	//	cout<< "Jet Eta\t"<<jet_eta<<endl;
	
    Bool_t bisHFMatched= false;
    	for(Size_t ijet_tk=0; ijet_tk<jet_ntk; ijet_tk++)
    	  {
    	    Int_t tkIndex      = L1TkJet_TTTrackIndex->at(iJet).at(ijet_tk);
    	    Int_t temp = FindTrackMatchingQuark(tkIndex);
    	    if(temp == 4 || temp == 5)
    	      {
    		bisHFMatched= true;
    		break;
    	      }
    	  }
	
    	for(Size_t ijet_tk=0; ijet_tk<jet_ntk; ijet_tk++)
    	  {
	    
    	    Int_t tkIndex      = L1TkJet_TTTrackIndex->at(iJet).at(ijet_tk);
	    
    	    Double_t tk_pt     = L1Tks_Pt->at(tkIndex);
    	    Double_t tk_eta    = L1Tks_Eta->at(tkIndex);
    	    Double_t tk_phi    = L1Tks_Phi->at(tkIndex);
    	    Double_t tk_pocax  = L1Tks_POCAx->at(tkIndex);
    	    Double_t tk_pocay  = L1Tks_POCAy->at(tkIndex);
    	    Double_t tk_pocaz  = L1Tks_POCAz->at(tkIndex);
    	    Double_t tk_chi2   = L1Tks_ChiSquared->at(tkIndex);
    	    Double_t tk_nstubs = L1Tks_Stubs_iLayer->at(tkIndex).size();
    	    Double_t tk_d0     = s->GetTrackD0(tkIndex);//-tk_pocax*TMath::Sin(tk_phi) + tk_pocay*TMath::Cos(tk_phi);
	   
    	    htk_pt    ->Fill(tk_pt);
    	    htk_eta   ->Fill(tk_eta);
    	    htk_phi   ->Fill(tk_phi);
    	    htk_pocax ->Fill(tk_pocax);
    	    htk_pocay ->Fill(tk_pocay);
    	    htk_pocaz ->Fill(tk_pocaz);
    	    htk_chi2  ->Fill(tk_chi2);
    	    htk_nstubs->Fill(tk_nstubs);
    	    htk_d0    ->Fill(tk_d0);

    	    // Int_t temp = FindTrackMatchingQuark(tkIndex);
    	    //cout<<temp<<endl;
    	    //  cout<< "\tTks Eta\t"<<tk_eta<<endl;
    	    //	    htk_absd0H-> Fill(abs(tk_d0));
    	    // cout<<"..........................................TK INDEX     "<<tkIndex<<endl;
    	    //Int_t temp = FindTrackMatchingQuark(tkIndex);
    	    // if(temp!=0)
    	      //      cout<<"\ttemp in main\t"<< temp<<endl;

    	    //if(bHFJet  && bIsBJet)
    	    if(bisHFMatched) 
    	      htk_absd0H-> Fill(abs(tk_d0));
    	    //if(bHF)
    	    //htk_absd0H-> Fill(abs(tk_d0));
    	    // else if (!bHFJet)
    	    else if(!bisHFMatched)  
    	      htk_absd0P-> Fill(abs(tk_d0));

    	  } // end loop on Jet Tracks
    	//***********************************************************************************************
	

      }// loop on jets


    // //*******************************************************************************************
    // // D0Dxy Phi and atan test
    // //*******************************************************************************************
    // for(Size_t iTks=0; iTks<L1Tks_Pt->size(); iTks++)
    //   {
    // 	Double_t tk_d0    = s->GetTrackD0(iTks);
    // 	Double_t tk_dxy   = s->GetTrackDxy(iTks);
    // 	Double_t tk_phi   = L1Tks_Phi->at(iTks);
    // 	Double_t tk_pocax = L1Tks_POCAx->at(iTks);
    // 	Double_t tk_pocay = L1Tks_POCAy->at(iTks);
    // 	Double_t tk_atan  = atan2(tk_pocay,tk_pocax);
    // 	//GetTrackD0Phi
    // 	Double_t DeltaD0  = tk_d0-tk_dxy;
    // 	Double_t DeltaPhi = tk_phi-tk_atan;
    // 	cout<<DeltaPhi<<endl;
    // 	hAllTks_d0  ->Fill(tk_d0);
    // 	hAllTks_dxy ->Fill(tk_dxy);
    // 	hAllTks_phi ->Fill(tk_phi);
    // 	hAllTks_atan->Fill(tk_atan);

    // 	hAllTks_Deltad0 ->Fill(DeltaD0);
    // 	hAllTks_Deltaphi ->Fill(DeltaPhi);
	
    //   }
    // //*******************************************************************************************



htk_absd0H->Scale(1/htk_absd0H->Integral());
htk_absd0P->Scale(1/htk_absd0P->Integral());

htk_absd0H->GetYaxis()->SetRangeUser(1E-7,0.4);
htk_absd0P->GetYaxis()->SetRangeUser(1E-7,0.4);

hPjetTkabsd0->Scale(1/hPjetTkabsd0->Integral());
hBjetTkabsd0->Scale(1/hBjetTkabsd0->Integral());

hPjetTkabsd0->GetYaxis()->SetRangeUser(1E-7,0.4);
hBjetTkabsd0->GetYaxis()->SetRangeUser(1E-7,0.4);







//************************************************************
//Laborious calculation of D0 phi..
		    // Double_t DaucosPhi = TMath::Cos(DauPhi);
		    
		    // Int_t    SigncosPhi= GetSign(DaucosPhi);
		    // Int_t    SignBPhi  = GetSign(BPhi);
		    
		//		cout<<SignD0<<endl;
		
		// if(DaucosPhi>0)
		//   {
		//     if(abs(DauPhi)>abs(BPhi))
		//       {
		// 	DauD0Phi = -(PI/2 - abs(DauPhi))*SignBPhi;
		//       }
		//     else if (abs(DauPhi)<abs(BPhi))
		//       {
		// 	DauD0Phi = (PI/2 + abs(DauPhi))*SignBPhi;
		//       }
		//   }
		// else if (DaucosPhi<0)
		//   {
		//     if(abs(DauPhi)>abs(BPhi))
		//       {
		// 	DauD0Phi = (abs(DauPhi) - PI/2)*SignBPhi;
		//       }
		//     else if (abs(DauPhi)<abs(BPhi))
		//       {
		// 	DauD0Phi = -(3*PI/2 - abs(DauPhi))*SignBPhi;
		//       }
		//   }


		
		//		cout<<"\tDau\t"<< DauIndx<< "\t"<<DauID<<"\t"<<DauPhi<<"\t\t"<<DauD0Phi<<endl;


		// if(CIDauPt<2)continue;
		// if(abs(CIDauEta)>2.5) continue;
		// //cout << "\t" << CIDauIndx <<"\t"<<CIDauID<<endl;
		// Double_t CIDauD0 = GetGenpD0(CIDauIndx);
		// hgend0BHadDirect->Fill(abs(CIDauD0));

//other stuff
	    Int_t    Gen_ID  = GenP_PdgId->at(iGenp);
	    Double_t BPt  = GenP_Pt->at(iGenp);
	    Double_t BEta = GenP_Eta->at(iGenp);
	    
	    TVector3 bvec;
	    bvec.SetPtEtaPhi(BPt,BEta,BPhi);

	    //    cout<< iGenp<< "\t"<<Gen_ID<<"\t"<<BPhi<<endl;
	    //cout<<vIDau.size()<<endl;


    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	if(IsFinalBHadron(iGenp))
	  {
	    vector<Size_t> vDau;
	    //	    GetBHadronChargedImmediateDaughters(iGenp, vDau);
	    GetFinalDaughter(iGenp, vDau);
	    Double_t BPhi = GenP_Phi->at(iGenp);

	    for(Size_t iDau=0; iDau<vDau.size(); iDau++)
	      {
		Int_t    DauIndx   = vDau.at(iDau);
		Int_t    DauID     = GenP_PdgId->at(DauIndx);
		Int_t    DauStatus = GenP_Status->at(DauIndx);
		Int_t    DauCharge = GenP_Charge->at(DauIndx);
		Double_t DauPt     = GenP_Pt->at(DauIndx);
		Double_t DauEta    = GenP_Eta->at(DauIndx);
		if(DauCharge != 0)// continue;
		  {
		    Double_t daux = GenP_VertexX->at(DauIndx);
		    Double_t dauy = GenP_VertexY->at(DauIndx);
		    TVector3 d0vec;
		    d0vec.SetXYZ(daux,dauy,0);

		    Double_t dotProd = bvec.Dot(d0vec);
		    
		    //	cout<<dotProd<<endl;
		    
		    Double_t DauPhi    = GenP_Phi->at(DauIndx);
		    Double_t DauD0  = GetGenpD0(DauIndx);
		    //Double_t SignD0 = GetSign(DauD0);
		    
		    Double_t SignD0;
		    
		    if(DauD0>0)
		      SignD0 = 1;
		    else if(DauD0<0)
		      SignD0 = -1;
		    
		    
		    
		    Double_t DauD0Phi = DauPhi + SignD0*(PI/2) ;
		    
		    if(DauD0Phi > PI)
		      DauD0Phi = DauD0Phi - 2*PI;
		    else if(DauD0Phi <= -PI)
		      DauD0Phi = DauD0Phi + 2*PI;
		    
		    
		    hgenphiD0HadDau->Fill(DauD0Phi);
		    Double_t deltaPhi = abs(auxTools_.DeltaPhi(DauD0Phi,BPhi));
		    hgendeltaphiD0HadDau->Fill(deltaPhi);
		    
		    Double_t DauAbsD0 = abs(GetGenpD0(DauIndx));
		    Double_t Sign =0;
		    if(deltaPhi<PI/2)
		      Sign = 1.0;
		    else if (deltaPhi>PI/2)
		      Sign = -1.0;
		    
		    // if(dotProd>0)
		    //   Sign = 1.0;
		    // else if (dotProd<0)
		    //   Sign = -1.0;
		    
		    Double_t DauSIP =  DauAbsD0 * Sign;
		    
		    hGenSIP->Fill(DauSIP);


	      }
	    // cout<<<< iGenp<<"\t"<<Gen_ID<<endl;
	  }// if final state BHad
	
      }// end loop on genp


 // for(Size_t iPix=0; iPix<nPixTk; iPix++)
    //   {
    // 	Double_t pix_d0 = s->GetPixTrackD0(iPix);
    // 	Double_t pix_eta = L1PixTks_Eta->at(iPix);
    // 	if ( IsWithinEtaRegion("Central", pix_eta) )
    // 	  {
    // 	    hpixeta_C->Fill(pix_eta);
    // 	    hpixd0_C->Fill(pix_d0);
    // 	  }
    // 	else if ( IsWithinEtaRegion("Intermediate", pix_eta) )
    // 	  {
    // 	    hpixeta_I->Fill(pix_eta);
    // 	    hpixd0_I->Fill(pix_d0);
    // 	  }
    // 	else if ( IsWithinEtaRegion("Forward", pix_eta) )
    // 	  {
    // 	    hpixeta_F->Fill(pix_eta);
    // 	    hpixd0_F->Fill(pix_d0);
    // 	  }
	
    // 	htk_d0->Fill(pix_d0);
    //   }
    



	//*********************************************************************************************************************

    for(Size_t iGenp=0; iGenp<nGenp; iGenp++)
      {
	if(IsFinalBHadron(iGenp))
	  {
	    //	    cout<<"\tThe Hadron\t"<<iGenp<<endl;
	    vector<Size_t> vDau;
	    GetBHadronChargedImmediateDaughters(iGenp, vDau);
	    // GetFinalDaughter(iGenp, vDau);
	    //cout<<"\tDau Num\t" << vDau.size()<<endl;
	    Double_t BPhi = GenP_Phi->at(iGenp);
	    for(Size_t iDau=0; iDau<vDau.size(); iDau++)
	      {
		Int_t    DauIndx   = vDau.at(iDau);
		Int_t    DauCharge = GenP_Charge->at(DauIndx);
		if(DauCharge != 0)// continue;
		  {
		    
		    Double_t DauPt    = GenP_Pt->at(DauIndx);
		    Double_t DauEta   = GenP_Eta->at(DauIndx);
		    Double_t DauPhi   = GenP_Phi->at(DauIndx);
		    Double_t DauD0    = GetGenpD0(DauIndx);
		    Double_t DauD0Phi = GetPhiD0(DauPhi, DauD0);
		    Double_t deltaPhi = abs(auxTools_.DeltaPhi(DauD0Phi,BPhi));
		    Double_t DauAbsD0 = abs(DauD0);
		    Double_t Sign     = GetSignofIP(deltaPhi);
		    Double_t DauSIP   = DauAbsD0 * Sign;
		    
		    if(DauPt >2 && abs(DauEta)<2.5)
		      {
			Int_t MatchedTP = GetGenMatchingTPIndex(DauIndx);
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
				    // hGenMatchedTkB->Fill(abs(TTtrackD0));
				    // hGenMatchedTPB->Fill(abs(tp_D0));
				    // hgend0BHadDirect->Fill(abs(DauD0));
				    // cout<<"==========================================================================================================="<<endl;
				    // cout<<DauIndx<<"\t"<<GenP_Pt->at(DauIndx)<<"\t"<<TP_Pt->at(MatchedTP)<<"\t"<<L1Tks_Pt->at(tp_tkIndx)<<"\t"<<L1PixTks_Pt->at(tk_pixIndx)<<endl;
				    // cout<<DauIndx<<"\t"<<GenP_Eta->at(DauIndx)<<"\t"<<TP_Eta->at(MatchedTP)<<"\t"<<L1Tks_Eta->at(tp_tkIndx)<<"\t"<<L1PixTks_Eta->at(tk_pixIndx)<<endl;
				    // cout<<DauIndx<<"\t"<<GenP_Phi->at(DauIndx)<<"\t"<<TP_Phi->at(MatchedTP)<<"\t"<<L1Tks_Phi->at(tp_tkIndx)<<"\t"<<L1PixTks_Phi->at(tk_pixIndx)<<endl;
				    // cout<<DauIndx<<"\t"<<DauD0<<"\t"<<tp_D0<<"\t"<<TTtrackD0<<"\t"<<TTPixtrackD0<<endl;
				    // cout<<"==========================================================================================================="<<endl;
				    // cout<< tk_pixIndx<<"\t"<<L1PixTks_Pt->at(tk_pixIndx)<<"\t"<<L1Tks_Pt->at(tp_tkIndx)<<"\t"<<TP_Pt->at(MatchedTP)<<endl;
				  } // if ttpixtrack found
			      }// if tttrack found
			    //  cout<<"\tTP Index\t"<< MatchTP<< "\tTP ID\t"<< TP_PdgId->at(MatchTP)<<endl;
			  }//if TP found
			hgend0BHadDirect->Fill(abs(DauD0));
		      }// pt eta cut on Genp 


		    hgenphiD0HadDau->Fill(DauD0Phi);
		    hgendeltaphiD0HadDau->Fill(deltaPhi);
		    hGenSIP->Fill(DauSIP);

		  }//if charged
	      }// end loop on daughters
	    // cout<<<< iGenp<<"\t"<<Gen_ID<<endl;
	  }// if final state BHad
	
      }// end loop on genp

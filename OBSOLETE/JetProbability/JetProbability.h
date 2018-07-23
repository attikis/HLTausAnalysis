#ifndef JetProbability_h
#define JetProbability_h

#include <iostream>

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/AuxTools.C"
#include "../utilities/Table.C"
#include "../utilities/Datasets.C"
#include "../utilities/L1Tracks.C"



class JetProbability : public TreeAnalyserMC
{
  public:
    JetProbability(const std::string SamplePath,
	const std::string SampleName,
	const std::string text_, 
	const int maxEvents_ = -1, 
	TTree* tree=0) : 
    TreeAnalyserMC("JetProbability", SamplePath, SampleName, text_, maxEvents_, tree) {
      auxTools_.StopwatchStart();
      mcSample = SampleName;
    };

    ~JetProbability() {};  

    virtual void Loop();

  private:
    
    std::string mcSample;
    AuxTools auxTools_;
    Datasets datasets_;

    // Add your private variables/methods here
    void  MakeHisto(void);
    Int_t GetPixTkIndex(Int_t TTTrackIndex);
    /* void  GetPixTkVector(Size_t jetIndx,  */
    /* 			 vector<Int_t> &pixTks); */
    
    /* Double_t DeltaPhid0(Double_t tkPhi, */
    /* 			Double_t jetPhi); */
    /* Int_t Signd0(Double_t dPhi); */
    Bool_t IsWithinEtaRegion(string etaRegion, 
			     Double_t eta);

    Bool_t IsWithinPtSlice(string pTSlice, 
			   Double_t pt);
    
    Bool_t IsTrackChi2Good(string chi2Sel, 
			   Double_t chi2);

    /* Bool_t IsJetInBDirection(Int_t JetIndex); */
    /* Bool_t IsJetInCDirection(Int_t JetIndex); */
    /* Bool_t IsJetInUDirection(Int_t JetIndex); */
    /* Bool_t IsJetInDDirection(Int_t JetIndex); */
    Bool_t IsFinalGenp (Size_t MotIndx,  
     			vector<unsigned short>& Daug); 
    /* Bool_t IsJetHeavyFlavour(Int_t JetIndex); */

    /* Int_t FindTrackMatchingQuark(Size_t Index); */
    
    /* Bool_t IsBHadronProduct(Size_t Indx); */
    /* Bool_t IsCHadronProduct(Size_t Indx); */
    /* Bool_t IsSHadronProduct(Size_t Indx); */
    /* Bool_t IsDHadronProduct(Size_t Indx); */
    /* Bool_t IsUHadronProduct(Size_t Indx); */
    
    void   myPrintGenp(Size_t Indx, bool bPrintHeaders);
  
    void GetFinalDaughter(Size_t Indx, 
			  vector<Size_t> &vFinalDau);
    
    /* Int_t MatchTrackToGenp(Size_t Index); */

    Double_t GetGenpD0(Size_t Index);
    Double_t GetTPD0(Size_t Index);

    Int_t GetGenMatchingTPIndex(Int_t Index);

    Int_t GetSign(Double_t Quantity);

    Bool_t IsBHadron(Int_t Index);
    
    Bool_t IsLiteHadron(Int_t Index);

    Bool_t IsFinalBHadron(Int_t Index);
    

    void GetBHadronChargedImmediateDaughters(Int_t Index ,
					     vector<Size_t> &vCIDau);

    Double_t GetPhiD0(Double_t PhiTk, Double_t D0);
    Double_t GetSignofIP(Double_t deltaPhi);


    Double_t GetTPDxy(Size_t Index);
    Double_t GetGenpDxy(Size_t Index);
    Double_t GetTTTrackDxy(Size_t Index);
    Double_t GetTTPixTrackDxy(Size_t Index);


    void FillSIPHistogramsD0(Int_t DauIndx, Double_t had_Phi,
    			     TH1D* GenAbsD0, TH1D* GenPhiD0, TH1D* GenDeltaPhi,TH1D* GenSIP,
    			     TH1D* TPAbsD0 , TH1D* TPPhiD0 , TH1D* TPDeltaPhi ,TH1D* TPSIP,
    			     TH1D* TkAbsD0 , TH1D* TkPhiD0 , TH1D* TkDeltaPhi ,TH1D* TkSIP,
    			     TH1D* PixAbsD0, TH1D* PixPhiD0, TH1D* PixDeltaPhi,TH1D* PixSIP);
    //void FillSIPHistogramsD0(Int_t DauIndx, Double_t had_Phi); 
    
    void FillSIPHistogramsDxy(Int_t DauIndx, Double_t had_Phi,
			      TH1D* GenAbsD0, TH1D* GenPhiD0, TH1D* GenDeltaPhi,TH1D* GenSIP,
			      TH1D* TPAbsD0 , TH1D* TPPhiD0 , TH1D* TPDeltaPhi ,TH1D* TPSIP,
			      TH1D* TkAbsD0 , TH1D* TkPhiD0 , TH1D* TkDeltaPhi ,TH1D* TkSIP,
			      TH1D* PixAbsD0, TH1D* PixPhiD0, TH1D* PixDeltaPhi,TH1D* PixSIP);


    void FillSIPHistogramsDotProd(Int_t DauIndx , TVector3 Hadvec,
				  TH1D* GenAbsD0, TH1D* GenSIP,
				  TH1D* TPAbsD0 , TH1D* TPSIP,
				  TH1D* TkAbsD0 , TH1D* TkSIP,
				  TH1D* PixAbsD0, TH1D* PixSIP);

    Double_t _eta_C;
    Double_t _eta_F;


    //Histograms

    TH1D* hpixd0_C;
    TH1D* hpixd0_I;
    TH1D* hpixd0_F;
    
    TH1D* hpixeta_C;
    TH1D* hpixeta_I;
    TH1D* hpixeta_F;

    TH1D* hjet_pt;    
    TH1D* hjet_eta;   
    TH1D* hjet_phi;   
    TH1D* hjet_e;     
    TH1D* hjet_vtx;   
    TH1D* hjet_ntk;   
    TH1D* hjet_multi; 
    
    TH1D* htk_pt;    
    TH1D* htk_eta; 
    TH1D* htk_phi;   
    TH1D* htk_pocax; 
    TH1D* htk_pocay; 
    TH1D* htk_pocaz;
    TH1D* htk_chi2;  
    TH1D* htk_nstubs;
    TH1D* htk_d0;    


    TH1D* hBGen_SignedD0;
    TH1D* hBGen_DeltaPhi;
    TH1D* hBGen_SignedD0Sig;
    TH1D* hsign;

    TH1D* hBGen_SignedD0_C;
    TH1D* hBGen_DeltaPhi_C;
    TH1D* hBGen_SignedD0Sig_C;
    TH1D* hsign_C;

    TH1D* hBGen_SignedD0_I;
    TH1D* hBGen_DeltaPhi_I;
    TH1D* hBGen_SignedD0Sig_I;
    TH1D* hsign_I;

    TH1D* hBGen_SignedD0_F;
    TH1D* hBGen_DeltaPhi_F;
    TH1D* hBGen_SignedD0Sig_F;
    TH1D* hsign_F;

    TH1D* htk_absd0H;
    TH1D* htk_absd0P;

    TH1D* hAlltk_d0;
    TH1D* hAllpix_d0;
    
    TH1D* hgend0B;
    TH1D* hgend0P;

    TH1D* hgend0BHad;
    TH1D* hgend0PHad;
    
    TH1D* hBjetTkabsd0;
    TH1D* hPjetTkabsd0;


    TH1D* hAllTks_d0;  
    TH1D* hAllTks_dxy;
    TH1D* hAllTks_phi; 
    TH1D* hAllTks_atan;
    
    TH1D* hAllTks_Deltad0;
    TH1D* hAllTks_Deltaphi;

    TH1D* hgend0EMu;
    TH1D* hGenMatchedTkMuE;
    TH1D* hGenMatchedPixMuE;
    TH1D* hGenMatchedTPMuE;

    TH1D* hgend0UDS;
    TH1D* hGenMatchedTkP;
    TH1D* hGenMatchedPixP;
    TH1D* hGenMatchedTPP;

    TH1D* hgend0BHadDirect;
    TH1D* hGenMatchedTkB;
    TH1D* hGenMatchedPixB;    
    TH1D* hGenMatchedTPB;

    //*********************************8
    TH1D* hGenBHadPhiD0;      
    TH1D* hGenBHaddeltaphiD0; 
    TH1D* hGenBHadSIP;
    
    TH1D* hGenMatchedTPBPhiD0;      
    TH1D* hGenMatchedTPBdeltaPhi;
    TH1D* hGenMatchedTPBSIP;        
    
    TH1D* hGenMatchedTkBPhiD0;      
    TH1D* hGenMatchedTkBdeltaPhi;
    TH1D* hGenMatchedTkBSIP;        
    
    TH1D* hGenMatchedPixBPhiD0;     
    TH1D* hGenMatchedPixBdeltaPhi;
    TH1D* hGenMatchedPixBSIP;        
    //*********************************8

    //*********************************8
    TH1D* hGenLHadPhiD0;      
    TH1D* hGenLHaddeltaphiD0; 
    TH1D* hGenLHadSIP;        

    TH1D* hGenMatchedTPPPhiD0;      
    TH1D* hGenMatchedTPPdeltaPhi;
    TH1D* hGenMatchedTPPSIP;        
    
    TH1D* hGenMatchedTkPPhiD0;      
    TH1D* hGenMatchedTkPdeltaPhi;
    TH1D* hGenMatchedTkPSIP;        
    
    TH1D* hGenMatchedPixPPhiD0;     
    TH1D* hGenMatchedPixPdeltaPhi;
    TH1D* hGenMatchedPixPSIP;       

    //*********************************8

    TH1D* hGenNoPUSIP; 
    TH1D* hgenphiHadDau;



    TH1D* hgend0BHadAll;
    TH1D* hgend0BHadNoC;
    TH1D* hgend0BHadNoCS;
    

    TH1D* htempgenD0;
    TH1D* htemptpD0;

    TH1D* htkDxy;


    TH1D* temp1;
    TH1D* temp2;
    TH1D* temp3;

    //*********************************************************************************Histo for SLices of Eta and PT
    
    TH1D* PixDeltaPhi_C_0203;
    TH1D* PixDeltaPhi_C_0304;
    TH1D* PixDeltaPhi_C_0405;
    TH1D* PixDeltaPhi_C_0506;
    TH1D* PixDeltaPhi_C_0607;
    TH1D* PixDeltaPhi_C_0708;
    TH1D* PixDeltaPhi_C_0809;
    TH1D* PixDeltaPhi_C_0910;
    TH1D* PixDeltaPhi_C_1015;
    TH1D* PixDeltaPhi_C_1520;
    TH1D* PixDeltaPhi_C_2000;
    
    TH1D* PixSIP_C_0203;
    TH1D* PixSIP_C_0304;
    TH1D* PixSIP_C_0405;
    TH1D* PixSIP_C_0506;
    TH1D* PixSIP_C_0607;
    TH1D* PixSIP_C_0708;
    TH1D* PixSIP_C_0809;
    TH1D* PixSIP_C_0910;
    TH1D* PixSIP_C_1015;
    TH1D* PixSIP_C_1520;
    TH1D* PixSIP_C_2000;

    TH1D* PixDeltaPhi_I_0203;
    TH1D* PixDeltaPhi_I_0304;
    TH1D* PixDeltaPhi_I_0405;
    TH1D* PixDeltaPhi_I_0506;
    TH1D* PixDeltaPhi_I_0607;
    TH1D* PixDeltaPhi_I_0708;
    TH1D* PixDeltaPhi_I_0809;
    TH1D* PixDeltaPhi_I_0910;
    TH1D* PixDeltaPhi_I_1015;
    TH1D* PixDeltaPhi_I_1520;
    TH1D* PixDeltaPhi_I_2000;
    
    TH1D* PixSIP_I_0203;
    TH1D* PixSIP_I_0304;
    TH1D* PixSIP_I_0405;
    TH1D* PixSIP_I_0506;
    TH1D* PixSIP_I_0607;
    TH1D* PixSIP_I_0708;
    TH1D* PixSIP_I_0809;
    TH1D* PixSIP_I_0910;
    TH1D* PixSIP_I_1015;
    TH1D* PixSIP_I_1520;
    TH1D* PixSIP_I_2000;

    TH1D* PixDeltaPhi_F_0203;
    TH1D* PixDeltaPhi_F_0304;
    TH1D* PixDeltaPhi_F_0405;
    TH1D* PixDeltaPhi_F_0506;
    TH1D* PixDeltaPhi_F_0607;
    TH1D* PixDeltaPhi_F_0708;
    TH1D* PixDeltaPhi_F_0809;
    TH1D* PixDeltaPhi_F_0910;
    TH1D* PixDeltaPhi_F_1015;
    TH1D* PixDeltaPhi_F_1520;
    TH1D* PixDeltaPhi_F_2000;
    
    TH1D* PixSIP_F_0203;
    TH1D* PixSIP_F_0304;
    TH1D* PixSIP_F_0405;
    TH1D* PixSIP_F_0506;
    TH1D* PixSIP_F_0607;
    TH1D* PixSIP_F_0708;
    TH1D* PixSIP_F_0809;
    TH1D* PixSIP_F_0910;
    TH1D* PixSIP_F_1015;
    TH1D* PixSIP_F_1520;
    TH1D* PixSIP_F_2000;


    TH1D* PixDeltaPhi_C; 
    TH1D* PixSIP_C;      
    
    TH1D* PixDeltaPhi_I; 
    TH1D* PixSIP_I;      
    
    TH1D* PixDeltaPhi_F; 
    TH1D* PixSIP_F;      

    
    TH1D* PixDeltaPhi_T_c2;
    TH1D* PixSIP_T_c2;     

    TH1D* PixDeltaPhi_M_c2;
    TH1D* PixSIP_M_c2;     

    TH1D* PixDeltaPhi_L_c2;
    TH1D* PixSIP_L_c2;     
    
    TH1D* PixDeltaPhi_VL_c2;
    TH1D* PixSIP_VL_c2;     


};

#endif // JetProbability_h

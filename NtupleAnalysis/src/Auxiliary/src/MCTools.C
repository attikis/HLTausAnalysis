#ifndef MCTools_cxx
#define MCTools_cxx

// User
#include "../interface/MCTools.h"

//#define DEBUG


// ****************************************************************************
bool MCTools::IsNeutrino(const int pdgId){
// ****************************************************************************
  //
  // Description:
  // Returns true if genParticle is neutrino (v_e, v_mu, v_tau), else false.
  //

  if( (abs(pdgId) == 12) || (abs(pdgId) == 14) || (abs(pdgId) == 16) ) return true;
  else return false;
}


// ****************************************************************************
bool MCTools::IsBoson(const int pdgId){
// ****************************************************************************
  //
  // Description:
  // Returns true if genParticle is boson, else false.
  //

  if( (abs(pdgId) == 21) || (abs(pdgId) == 22) || (abs(pdgId) == 23) ||
      (abs(pdgId) == 24) || (abs(pdgId) == 25) || (abs(pdgId) == 32) ||
      (abs(pdgId) == 33) || (abs(pdgId) == 34) || (abs(pdgId) == 35) ||
      (abs(pdgId) == 36) || (abs(pdgId) == 37) ) return true;
  else return false;
}


// ****************************************************************************
bool MCTools::IsLepton(const int pdgId){
// ****************************************************************************  
  //
  // Description:
  // Returns true if genParticle is a lepton (e, mu, tau & associated neutrinos), else false.
  //

  if( (abs(pdgId) == 11) || (abs(pdgId) == 12)  ||
      (abs(pdgId) == 13) || (abs(pdgId) == 14)  ||
      (abs(pdgId) == 15) || (abs(pdgId) == 16) ) return true;
  else return false;
}


// ****************************************************************************
bool MCTools::IsChargedLepton(const int pdgId){
// ****************************************************************************
  //
  // Description:
  // Returns true if genParticle is a charged lepton (e, mu, tau), else false.
  //

  if( (abs(pdgId) == 11) || (abs(pdgId) == 13) || (abs(pdgId) == 15) ) return true;
  else return false;
}


// ****************************************************************************
bool MCTools::IsQuark(const int pdgId){
// ****************************************************************************

  //
  // Description:
  // Returns true if genParticle is a quark (u, d, c, s, t, b), else false.
  //

  if( (abs(pdgId) == 1) || (abs(pdgId) == 2)  ||
      (abs(pdgId) == 3) || (abs(pdgId) == 4)  ||
      (abs(pdgId) == 5) || (abs(pdgId) == 6) ) return true;
  else return false;
}


/* marina1
// ****************************************************************************
bool MCTools::RecursivelyLookForMotherId(Int_t Indx,
					 Int_t MoId,
					 const bool posn)
// ****************************************************************************
{
  unsigned short nMoms = GenP_Mothers->at(Indx).size();
  if (nMoms == 0)return false;
  Int_t MyId;

#ifdef DEBUG
  std::cout << "Looking for mother #" << MoId << " of particle "<< Indx << std::endl;
#endif

  for (unsigned short i = 0; i < nMoms; i++){
    if (!posn) {
      MyId = fabs(GenP_PdgId->at(GenP_Mothers->at(Indx).at(i)));
      MoId = fabs(MoId);
    }
    else {
      MyId = GenP_PdgId->at(GenP_Mothers->at(Indx).at(i));
    } 
    if (MyId == MoId)
      return true;
    if (RecursivelyLookForMotherId(GenP_Mothers->at(Indx).at(i), MoId, posn) )
      return true;
  }
  return false;
}


// *****************************************************************
Int_t MCTools::PosOfMotherId(Int_t Indx,
			   Int_t MoId,
			   const bool posn)
// *****************************************************************
{
  unsigned short nMoms = GenP_Mothers->at(Indx).size();
  if (nMoms == 0)return 65535;
  Int_t MyId;
  Int_t Position;
  for (unsigned short i = 0; i < nMoms; i++){
    Position = GenP_Mothers->at(Indx).at(i);
    if (!posn) {
      MyId = fabs(GenP_PdgId->at(Position));
      MoId = fabs(MoId);
    }
    else {
      MyId = GenP_PdgId->at(Position);
    } 
    Int_t status = GenP_Status->at(Position);
    if (MyId == MoId) {
      if (status == 3) {
	if(GenP_Daughters->at(Position).size()>2)
	  {
	    Position = GenP_Daughters->at(Position).at(2);
	    return Position;
	  }
	else
	  {
	    return 65535;
	  }
      }
      else {
	return Position;
      }
    }
    else {
      Position = PosOfMotherId(Position,MoId,posn);
      return Position;
    }
  }
  return 65535;
}

// *****************************************************************
bool MCTools::LookForMotherId(Int_t Indx,
			      Int_t MoId,
			      const bool posn)
// *****************************************************************
{
  unsigned short nMoms = GenP_Mothers->at(Indx).size();
  if (nMoms == 0)
    {
      std::cout << "Look  No Mothers\n"; 
      return false;
    }
  Int_t MyId;
  for (unsigned short i = 0; i < nMoms; i++) {
    if (!posn) {
      MyId = fabs(GenP_PdgId->at(GenP_Mothers->at(Indx).at(i)));
      MoId = fabs(MoId);
    }
    else {
      MyId = GenP_PdgId->at(GenP_Mothers->at(Indx).at(i));
    } 
    if (MyId == MoId) return true;
  }
  return false;
}


// ****************************************************************************
TLorentzVector MCTools::GetP4(const Int_t iGenP)
// ****************************************************************************
{

  // Sanity check on 4-momentum
  TLorentzVector p4;
  p4.SetPtEtaPhiM(GenP_Pt->at(iGenP), GenP_Eta->at(iGenP), GenP_Phi->at(iGenP), GenP_Mass->at(iGenP));
  return p4;

}





// ****************************************************************************
TLorentzVector MCTools::GetVisibleP4(const std::vector<unsigned short>& Daug)
// ****************************************************************************
{

  
  //  Returns the 4-vector sum of all visible daughters. If one would use this
  //  to calculate the visible eT from a hadronic tau then the eT calculated would 
  //  correcpond to the vector-sum eT of all decay products (and NOT the scalar eT sum)
  //  In fact, the eT calculated using this 4-vector would be identical to the eT of 
  //  the tau itself or that of the intermediate resonance (e.g. rho, alpha_1)
 

  TLorentzVector p4;
  if (Daug.size() == 0) return p4;

  // For-loop: Daughters
  for (unsigned short i = 0; i< Daug.size(); i++){

    const unsigned short index = Daug.at(i);
    const Int_t id             = abs( GenP_PdgId->at(index) );

    // Skip invisible daughters (neutrinos)
    if( (id == 12)  || (id == 14)  || (id == 16) ){ continue; }

    // Get properties
    Double_t pt   = GenP_Pt  ->at(index);
    Double_t eta  = GenP_Eta ->at(index);
    Double_t phi  = GenP_Phi ->at(index);
    Double_t mass = GenP_Mass->at(index);

    TLorentzVector tmp;
    tmp.SetPtEtaPhiM(pt, eta, phi, mass);
    p4 += tmp;

  } // For-loop: Daughters

  return p4;
}

// ****************************************************************************
Int_t MCTools::GetLdgDaughterIndex(std::vector<unsigned short> Daug, bool bOnlyChargedDaughters)
// ****************************************************************************
{

  // Declarations
  Double_t ldgPt   = -1.0;
  Int_t ldgPtIndex = -1;
  if (Daug.size() == 0) return ldgPtIndex;

  // For-loop: Daughters
  for (size_t i = 0; i< Daug.size(); i++){
    
    const Int_t Indx = Daug.at(i);
    const Int_t id   = fabs( GenP_PdgId->at(Indx) );

    // Skip neutrinos
    if( (id == 12)  || (id == 14)  || (id == 16) ) continue;

    // Get pt and charge
    Double_t pt     = GenP_Pt->at(Indx);
    Double_t charge = GenP_Charge->at(Indx);
    
    if(bOnlyChargedDaughters && fabs(charge) < 1 ) continue;

    // Find leading daughter index
    if (pt > ldgPt ){
      ldgPt = pt;
      ldgPtIndex = Indx;
    }
    
  }  // For-loop: Daughters

  return ldgPtIndex;
}


// ****************************************************************************
Double_t MCTools::GetHadronicTauMaxSignalCone(std::vector<unsigned short> Daug, 
					      bool bOnlyChargedDaughters, 
					      double minPt)
// ****************************************************************************
{

  if (Daug.size() <= 1) return -1;

  // Get Ldg Charged Track index and direction
  const Int_t ldgIndex  = GetLdgDaughterIndex(Daug, true);
  const Double_t ldgEta = GenP_Eta->at(ldgIndex);
  const Double_t ldgPhi = GenP_Phi->at(ldgIndex);
  Double_t deltaRMax    = -1.0;
  Int_t indexMax        = -1.0;
  Double_t dEta         = -1.0;
  Double_t dPhi         = -1.0;

  // For-loop: Daughters
  for (size_t i = 0; i< Daug.size(); i++){

    const Int_t index = Daug.at(i);
    // Skip self
    if (index == ldgIndex) continue;

    // Skip neutrinos
    Int_t Id = abs(GenP_PdgId->at(index));
    if( (Id == 12)  || (Id == 14)  || (Id == 16) ){ continue; }

    // Get Daughter properties
    Double_t pt     = GenP_Pt ->at(index);
    Double_t eta    = GenP_Eta->at(index);
    Double_t phi    = GenP_Phi->at(index);
    Double_t charge = GenP_Charge->at(index);

    // Consider only charged daughters above minPt
    if (pt < minPt) continue;
    if (bOnlyChargedDaughters){ if( fabs(charge) < 1 ) continue; }

    Double_t deltaR = tools.DeltaR(ldgEta, ldgPhi, eta, phi);
    if (deltaR > deltaRMax){
      deltaRMax = deltaR;
      indexMax  = index;
      dEta  = tools.DeltaEta(ldgEta, eta);
      dPhi = tools.DeltaPhi(ldgPhi, phi);
    }
    
  } // For-loop: Daughters

  return deltaRMax;
}


// ****************************************************************************
void MCTools::GetHadronicTauChargedOrNeutralPions(Int_t tauIndex, 
						  Int_t charge,
						  std::vector<unsigned short> &chargedPions)
// ****************************************************************************
{
  
  if (GenP_Daughters->at(tauIndex).size() == 0) return;

  // Get the pi+/-,pi0, K+/-, K0,K0L,KOS,eta,omegas and gammas
  std::vector<unsigned short> hTau_Dau;
  GetHadronicTauFinalDaughters(tauIndex, hTau_Dau);
  
  // For-loop: Daughters
  for (unsigned short i = 0; i< hTau_Dau.size(); i++){

    // Get Daughter properties
    Int_t daughter_index     = hTau_Dau.at(i);
    Double_t daughter_charge = GenP_Charge->at(daughter_index);

    // Keep only the pi+/-, K+/-, omegas
    if( fabs(daughter_charge) != charge ) continue;
    
    // Save to container
    chargedPions.push_back(daughter_index);
    
  } // For-loop: Daughters
  
  return;
}


// ****************************************************************************
void MCTools::GetHadronicTauChargedPions(Int_t tauIndex, 
					 std::vector<unsigned short> &chargedPions)
// ****************************************************************************
{
  
  GetHadronicTauChargedOrNeutralPions(tauIndex, 1, chargedPions);
  return;
}


// ****************************************************************************
void MCTools::GetHadronicTauNeutralPions(Int_t tauIndex, 
					 std::vector<unsigned short> &neutralPions)
// ****************************************************************************
{
  
  GetHadronicTauChargedOrNeutralPions(tauIndex, 0, neutralPions);
  return;
}


// ****************************************************************************
void MCTools::GetHadronicTauFinalDaughters(Int_t Indx,
					   std::vector<unsigned short>& Daug)
// ****************************************************************************
{

  std::cout << "GenP_PdgId->size() = " << GenP_PdgId->size() << std::endl;
  std::cout << "GenP_PdgId->at("<<Indx<<") = " << GenP_PdgId->at(Indx) << std::endl;
  std::cout << "GenP_Daughters->at(0).size() = " << GenP_Daughters->at(0).size() << std::endl;
  
  if (GenP_Daughters->at(Indx).size() == 0) return;
  for (unsigned short i = 0; i< GenP_Daughters->at(Indx).size(); i++)
    {

    Int_t IndxDau = GenP_Daughters->at(Indx).at(i);
    Int_t IdDau = fabs(GenP_PdgId->at(IndxDau));
    Int_t IdMo  = fabs(GenP_PdgId->at(Indx));
    // Keep only the pi+/-,pi0, K+/-, 
    // K0,K0L,KOS,eta,omegas and gammas from tau->tau+gamma transition
    if ( (IdDau == 111 || IdDau == 211 || IdDau == 321 ||   //pi0,pi+/-,K+/-
          IdDau == 130 || IdDau == 310 || IdDau == 311 ||   //K0L,K0S,K0
	  IdDau == 211 || IdDau == 223)                     //eta and omega
       ||(IdDau == 22 && IdMo==15 && Daug.size()!=0 ) ){  //Avoid leptonic decay
      // Because of the mixing of Generator and Detector simulation particles
      // in the list check if the particle has a parent already in the list
      // If it does then it comes from a hadronic Int_teraction with the
      // detector material and it is not part of the tau decay

    Int_t ifound = 0;      
      for (unsigned short j=0; j< Daug.size(); j++){
	if (RecursivelyLookForMotherId(IndxDau, GenP_PdgId->at(Daug[j]),true)) {
	  ifound += 1;
	}
      }
      if (ifound == 0) Daug.push_back(IndxDau);
    }

    GetHadronicTauFinalDaughters(IndxDau,Daug);
  }
  return;
}

// // **************************************************************************
// bool MCTools::IsHadronicTauDecay(Int_t Indx)
// // **************************************************************************
// {
//   std::vector<unsigned short> Daug_tmp;
//   GetHadronicTauFinalDaughters(Indx, Daug_tmp);
//   if (Daug_tmp.size() != 0) {return true;}
//   else {return false;}

// }


// **************************************************************************
bool MCTools::IsFinalStateTau(const std::vector<unsigned short>& Daug)
// **************************************************************************
{
  
  for (unsigned short i = 0; i < Daug.size(); i++){
    
    const unsigned short Indx = Daug[i];
    Int_t id = GenP_PdgId->at(Indx);
    if ( TMath::Abs(id) == 15) return false;
  
  }
  return true;

}


// *************************************************************************
Bool_t MCTools::IsFinalStateHadronicTau(Int_t indx)
// *************************************************************************
 {

   // Only consider taus
   bool bIsTau = (TMath::Abs( GenP_PdgId->at(indx) ) == 15);
   if (!bIsTau) return false;

   // Only consider final state taus (do not decay to self)
   bool bIsFinalStateTau  = IsFinalStateTau( GenP_Daughters->at(indx) );
   if (!bIsFinalStateTau) return false;

   // // Only consider tau-jets
   // bool bIsHadronicTauDecay = IsHadronicTauDecay(indx);
   // if (!bIsHadronicTauDecay) return false;

   return true;
 }



// ****************************************************************************
Bool_t MCTools::IsTriggerHadronicTau(const int iGenP, 
				     const bool bApplyTkAcceptanceCut, 
				     const int trueMomId, 
				     const bool bUseAbsMomId,
				     const double visEt)
// ****************************************************************************
{
    
  // Only consider final state hadronic taus
  if ( !IsFinalStateHadronicTau(iGenP) ) return false;

  // Get hTau properties
  std::vector<unsigned short> hTau_Products;
  GetHadronicTauFinalDaughters(iGenP, hTau_Products);
  
  // Acceptance cuts using visible 4-momenta
  TLorentzVector hTau_visP4 = GetVisibleP4(hTau_Products);
  if (hTau_visP4.Et() <= visEt) return false;

  const bool bIsOutsideEtaAcceptance = fabs(hTau_visP4.Eta()) > 2.3;
  if (bApplyTkAcceptanceCut * bIsOutsideEtaAcceptance) return false;

  // Only allow defined hadronic decay modes (reject decays with leptons e or mu)
  Int_t hTau_decayMode = GetGenTauDecayMode(hTau_Products);
  if (hTau_decayMode == -1) return false;

  // Use the default value 0 to disable looking for a true mother
  if (trueMomId == 0) return true;
  else return RecursivelyLookForMotherId(iGenP, trueMomId, bUseAbsMomId);

}


// *************************************************************************
Int_t MCTools::GetGenTauDecayMode(std::vector<unsigned short> Daug)
// *************************************************************************
{
  unsigned int nPipm = 0;
  unsigned int nPi0s = 0;
  unsigned int nKpm  = 0;
  unsigned int nCharged=0;
  unsigned int nOtherCharged = 0; // Other than pions, Kaons.
  unsigned int nOtherNeutral = 0; // K_S, K_L, Sigma, etc..
  unsigned int nLeptons = 0; // e, mu

  if (Daug.size() <= 0) return -1;

  // For-loop: Daughters
  for (unsigned short j = 0; j < Daug.size(); j++) {
    Int_t partId   = fabs(GenP_PdgId->at(Daug.at(j)));
    float charge = fabs(GenP_Charge->at(Daug.at(j)));
    
    // Leptonic decays
    if ( fabs(partId) == 11 || fabs(partId) == 13 ) nLeptons++;
    
    if (charge > 0) {
      nCharged++;
      if (partId == 211) {
	nPipm++;
      }
      else if (partId == 321){
	nKpm++;
      }
      else {
	nOtherCharged++;
      }
    }
    else {
      if (partId == 111){
	nPi0s ++;
      }
      else{
	nOtherNeutral++;
      }
    }
  }  // For-loop: Daughters
  
  // Determine return value
  if (nLeptons !=0 ) return -1;
  
  if (nOtherNeutral !=0 ) {
    return ((nCharged/2)*10 + 8);      // it returns 8, 18, 28 for 1, 3, 5 prong
  }
  else{
    return ((nCharged/2)*10 + nPi0s); // it returns a value depending on nProngs+nPi0s (i.e. 3prongs+3pi0 =13)
  }
}


// ******************************************************************************
void MCTools::PrintTauDecayStatistics(const unsigned int NTaus,
				      const unsigned int Decay[][7],
				      const unsigned int Nrows)
// ******************************************************************************
{
  double HadDecaysFrac = 0.6476; // The total tau->hadrons from PDG 
                          // The number is derived from total (1) 
                          // the leptonic decays: 0.3524 (e and mu) 
  std::cout << "=========================================" << std::endl;
  std::cout << "    Statistics of Hadronic Tau Decays    " << std::endl;
  std::cout << " ** Adjusted to 64.76% as expected (PDG) " << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << " Number of hadronic tau decays: " << std::setw(6) << NTaus 
	    << std::endl;
  std::cout << "=========================================" << std::endl;
  for (unsigned int i = 0; i < Nrows; i++){
    double fraction1 = Decay[i][6]*100*HadDecaysFrac/NTaus;
    std::cout << std::setw(1) << (i*2+1) << " prong decays:  " 
	      << std::setw(5) << Decay[i][6]
	      << "--->"  << std::setw(6) << std::setprecision(3) << fraction1 << "%"
	      << std::endl;  
    for (unsigned int j = 0; j < 5; j++) {
      double fraction2 = Decay[i][j]*100*HadDecaysFrac/NTaus;
      std::cout << "    " << std::setw(1) << (i*2+1) << "-prong"<< std::setw(1)
		<< j << "Pi0s:" << std::setw(5) << Decay[i][j] << "--->"
		<< std::setw(6) << std::setprecision(3) << fraction2 << "%"
		<< std::endl;
    }
    double fraction2 = Decay[i][5]*100*HadDecaysFrac/NTaus;
    std::cout << "    " << std::setw(1) << (i*2+1) << "-prong" 
	      << "+OtherNeutrals:" << std::setw(5) << Decay[i][5] << "--->"
	      << std::setw(6) << std::setprecision(3) << fraction2 << "%"
	      << std::endl;
  }
  std::cout << "=========================================" << std::endl;
}

// ********************************************************
void MCTools::PrintAllDaughters(Int_t Indx)
// ********************************************************
{

  if (GenP_Daughters->at(Indx).size() == 0) return;

  for (unsigned int i=0; i< GenP_Daughters->at(Indx).size(); i++)
    {
      Int_t DauIndx = GenP_Daughters->at(Indx).at(i);
      std::cout << std::fixed;
      std::cout << "(" << std::setw(4) << Indx << ")" << std::setw(6) 
		<< GenP_PdgId->at(Indx)<< " ---> " 
		<< "(" << std::setw(4) << DauIndx << ")" << std::setw(6)
		<< GenP_PdgId->at(DauIndx)
		<< std::endl;
      // Print all daughters of the daughter itself
      PrintAllDaughters(DauIndx);
    }
  return;
}


// ********************************************************
void MCTools::PrintAllDaughtersMinimalInfo(Int_t Indx)
// ********************************************************
{

  if (GenP_Daughters->at(Indx).size() == 0) return;

  for (unsigned int i=0; i< GenP_Daughters->at(Indx).size(); i++)
    {
      Int_t DauIndx = GenP_Daughters->at(Indx).at(i);
      PrintGenpMinimalInfo(DauIndx, false);
      PrintAllDaughtersMinimalInfo(DauIndx);
    }
  return;
}

// ********************************************************
void MCTools::PrintGenpDaughters(Int_t Indx, bool bPrintHeaders)
// ********************************************************
{

  if (GenP_Daughters->at(Indx).size() == 0) return;

  for (unsigned int i=0; i< GenP_Daughters->at(Indx).size(); i++)
    {
      Int_t DauIndx = GenP_Daughters->at(Indx).at(i);
      PrintGenp(DauIndx, bPrintHeaders);
    }
  return;
}

// ******************************************************
void MCTools::PrintAllMothers(Int_t Indx)
// ******************************************************
{

  if (GenP_Mothers->at(Indx).size() == 0) return;
  for (unsigned int i = 0; i < GenP_Mothers->at(Indx).size(); i++){
    Int_t MoIndx = GenP_Mothers->at(Indx).at(i);
    std::cout << std::fixed;
    std::cout << "(" << std::setw(4) << Indx << ")" << std::setw(4) 
	      << GenP_PdgId->at(Indx)<< " <--- " << "(" << std::setw(3) 
	      << MoIndx << ")" << std::setw(4) << GenP_PdgId->at(MoIndx)
	      << std::endl;
    PrintAllMothers(MoIndx);
  }
  return;
}

// **************************************************
void MCTools::PrintGenp(Int_t Indx, bool bPrintHeaders)
// **************************************************
{
  unsigned int moth1 = 0;
  unsigned int moth2 = 0;
  unsigned int daug1 = 0;
  unsigned int daug2 = 0;
  unsigned int NDaug;
  static Int_t evOld = 0;

  if ( (evOld != EvtNumber) || (Indx == 0) ) {
    std::cout << std::endl;
    evOld = EvtNumber;

    if (bPrintHeaders)
      {
	std::cout << std::setw(10)  << ">>> Run"
		  << std::setw(10)  << "Event"
		  << std::setw(10)  << "Indx"  
		  << std::setw(10) << "Id" 
		  << std::setw(10)  << "Status"
		  << std::setw(10)  << "Mot1"
		  << std::setw(10)  << "Mot2"
		  << std::setw(10)  << "Dau1"
		  << std::setw(10)  << "Dau2"
		  << std::setw(10)  << "px"
		  << std::setw(10) << "py"
		  << std::setw(10) << "pz"
		  << std::setw(10) << "E"
		  << std::setw(10) << "Eta"
		  << std::setw(10) << "Phi"
		  << std::setw(10) << "M"
		  << std::setw(10) << "Vtx-Z"
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
  std::cout << std::fixed;
  std::cout << std::setw(10) << std::setprecision(0) << RunNumber
	    << std::setw(10) << std::setprecision(0) << EvtNumber
	    << std::setw(10) << std::setprecision(0) << Indx
	    << std::setw(10) << std::setprecision(0) << GenP_PdgId->at(Indx)
	    << std::setw(10) << std::setprecision(0) << GenP_Status->at(Indx)
	    << std::setw(10) << std::setprecision(0) << moth1
	    << std::setw(10) << std::setprecision(0) << moth2
	    << std::setw(10) << std::setprecision(0) << daug1
	    << std::setw(10) << std::setprecision(0) << daug2
	    << std::setw(10) << std::setprecision(3) << px
	    << std::setw(10) << std::setprecision(3) << py
	    << std::setw(10) << std::setprecision(3) << pz
	    << std::setw(10) << std::setprecision(3) << ene
	    << std::setw(10) << std::setprecision(3) << eta
	    << std::setw(10) << std::setprecision(3) << phi
	    << std::setw(10) << std::setprecision(3) << mass
	    << std::setw(10) << std::setprecision(3) << GenP_VertexZ->at(Indx)
	    << std::endl;
}


// **************************************************
void MCTools::PrintGenpFullInfo(Int_t Indx, bool bPrintHeaders)
// **************************************************
{

  // std::cout << "*** GenpFullInfo (" << Indx << ") " << std::endl;
  PrintGenp(Indx, bPrintHeaders);
  std::cout << "    Mothers: " << std::endl;
  PrintAllMothers(Indx);
  std::cout << "   Daughters: " << std::endl;
  PrintAllDaughters(Indx);
  
  return;
}


// **************************************************
void MCTools::PrintGenpMinimalInfo(Int_t Indx, bool bPrintHeaders)
// **************************************************
{

  double pt    = GenP_Pt->at(Indx);
  double eta   = GenP_Eta->at(Indx);
  double phi   = GenP_Phi->at(Indx);
  double mass  = GenP_Mass->at(Indx);
  double PdgId = GenP_PdgId->at(Indx);

  TLorentzVector pGen;
  pGen.SetPtEtaPhiM(pt,eta,phi,mass);
  double E   = pGen.E();
  double Et  = pGen.Et();
  double Eta = pGen.Eta();
  double Phi = pGen.Phi();

  std::cout << std::fixed;
  if (bPrintHeaders)
    {
      std::cout << std::setw(10) << ">>> Run"
		<< std::setw(10) << "Event"
		<< std::setw(10) << "Indx"  
		<< std::setw(10) << "Energy"
		<< std::setw(10) << "Et"
		<< std::setw(10) << "Eta"
		<< std::setw(10) << "Phi"
		<< std::setw(10) << "Id" 
		<< std::endl;
    }

  std::cout << std::setw(10) << std::setprecision(0) << RunNumber
	    << std::setw(10) << std::setprecision(0) << EvtNumber
	    << std::setw(10) << std::setprecision(0) << Indx
	    << std::setw(10) << std::setprecision(3) << E
	    << std::setw(10) << std::setprecision(3) << Et
	    << std::setw(10) << std::setprecision(3) << Eta
	    << std::setw(10) << std::setprecision(3) << Phi
	    << std::setw(10) << std::setprecision(0) << PdgId
	    << std::endl;
  
  return;
}


// ****************************************************************************
TLorentzVector MCTools::GetGenpP4(const Int_t Index)
// ****************************************************************************
{

  // Construct the 4-momentum of the gen-particles
  Double_t pt   = GenP_Pt  ->at(Index);
  Double_t eta  = GenP_Eta ->at(Index);
  Double_t phi  = GenP_Phi ->at(Index);
  Double_t mass = GenP_Mass->at(Index);
    
  TLorentzVector p4;
  p4.SetPtEtaPhiM(pt, eta, phi, mass);
  return p4;
}



// ****************************************************************************
double MCTools::GetLxy(const Int_t iGenP,
		       double refPoint_X,
		       double refPoint_Y)
// ****************************************************************************
{

  // The distance traversed by a long-lived particle is Lxy (decay length).
  // Assuming the origin (0, 0) as the Primary Vertex (the point where the
  // particle was produced) and  (vtx_X, vtx_Y) as the Secondary Vertex (the
  // point where the long-lived particle reaches before decaying), then Lxy is
  // obtained from Pythagoras theorem.
  // Particles that decay promptly should have Lxy very close to zero.
  double refP_X = refPoint_X;
  double refP_Y = refPoint_Y;
  double vtx_X  = GenP_VertexX->at(iGenP);
  double vtx_Y  = GenP_VertexY->at(iGenP);
  double Lxy    = sqrt(pow( (vtx_X - refP_X), 2) + pow( (vtx_Y - refP_Y), 2));

  return Lxy;
}


// ****************************************************************************
double MCTools::GetD0Mag(const Int_t iGenP,
			 const Int_t iMom)
// ****************************************************************************
{

  // The distance traversed by a long-lived particle is Lxy (decay length).
  // Particles that decay promptly should have Lxy very close to zero.
  // Simple trigonometry will reveal that the tangent of the angle between the
  // long-lived particles' decay product and it's own direction will give:
  // sin( |phi_mom - phi_Phi| ) = d0/Lxy
  // Then d0 can be simply obtained by multiplying the tangent of the azimuthal
  // angle difference with the decay length Lxy:
  double Lxy      = GetLxy(iGenP);
  double phi      = GenP_Phi->at(iGenP);
  double phi_mom  = GenP_Phi->at(iMom);
  double deltaPhi = fabs(phi - phi_mom);
  double d0Mag    = sin(phi)*Lxy;

  return d0Mag;
}
*/ //marina1

#endif  // MCTools_cxx


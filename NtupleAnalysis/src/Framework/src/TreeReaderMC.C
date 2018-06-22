#ifndef TreeReaderMC_cxx
#define TreeReaderMC_cxx

// User
#include "../interface/TreeReaderMC.h"

// ROOT
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//============================================================================
void TreeReaderMC::InitVars_()
//============================================================================
{

  cfg_DEBUG = false;
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::InitVars_()" << std::endl;
  
  return;
}


//============================================================================
void TreeReaderMC::GetHadronicTauFinalDaughters(GenParticle p,
					    vector<unsigned short> &Daug)
//============================================================================
{

  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetHadronicTauFinalDaughters()" << std::endl;
  //
  // Description:
  // Keep only the detector-related particles (tracker/calo)
  //  pi+/-, pi0, K+/-, K0, K0L, KOS.
  // Avoid leptonic decays.
  //
  // Additional info:
  // We keep the K+/-(pdgId=321) and K0L(pdgId=130) and we don't get their decay
  // products (pi+/-) because typically they live long enought to reach the ECAL/HCAL
  // c*tau (K+/-) =  3.713 m
  // c*tau (K0L)  = 15.33  m
  // c*tau (K0S)  =  0.026842 m
  //
  
  // TreeReaderMC check
  if (p.daughtersIndex().size() < 1) return;

  // Define the allowed final daugther pdgIds
  std::vector<unsigned int> pdgIds;
  pdgIds.push_back(111); // pi0
  pdgIds.push_back(211); // pi+
  pdgIds.push_back(16);  // nu_tau
  pdgIds.push_back(130); // K0_L
  pdgIds.push_back(310); // K0_S
  pdgIds.push_back(311); // K0
  pdgIds.push_back(321); // K+
  // pdgIds.push_back(223);   // omega0
  // pdgIds.push_back(221);   // eta0
  // pdgIds.push_back(213);   // rho(770)+
  // pdgIds.push_back(20213); // alpha1(1260)+ 
    

  // For-loop: All daughters
  for (unsigned int i = 0; i < p.daughtersIndex().size(); i++) 
    {

      unsigned int index = p.daughtersIndex().at(i);
      GenParticle d      = GetGenParticle(index);
      int genD_pdgId     = abs(d.pdgId());


      // Check if particle is of interest
      if ( std::find(pdgIds.begin(), pdgIds.end(), genD_pdgId) == pdgIds.end() )
	{
	  GetHadronicTauFinalDaughters(d, Daug);
      	}
      else
	{
	  Daug.push_back( d.index() );
	}
    }
  
  return;

}


//============================================================================
vector<GenParticle> TreeReaderMC::GetGenParticles(bool bPrintList)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetGenParticles()" << std::endl;
    
  vector<GenParticle> theGenParticles;
  Table info("Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | Mothers | Daughters |", "Text"); //LaTeX or Text    
  
  // For-loop: All GenParticles
  for (Size_t genP_index = 0; genP_index < GenP_Pt.size(); genP_index++)
    {
      GenParticle p = GetGenParticle(genP_index);
      SetGenParticleMomsAndDads(p);
      theGenParticles.push_back(p);

      info.AddRowColumn(genP_index, auxTools_.ToString(p.index())  );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.pdgId())  );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.status()) );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.charge()) );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.pt())     );
      if (p.pdgId() != 2212) info.AddRowColumn(genP_index, auxTools_.ToString(p.eta())    );
      else info.AddRowColumn(genP_index, "infty");
      info.AddRowColumn(genP_index, auxTools_.ToString(p.phi())    );
      info.AddRowColumn(genP_index, auxTools_.ToString(p.energy()) );
      info.AddRowColumn(genP_index, "(" + auxTools_.ToString(p.vx())  + ", " + auxTools_.ToString(p.vy())  + ", " + auxTools_.ToString(p.vz())  + ", " + ")");
      if (p.mothersIndex().size() < 4)
	{
	  info.AddRowColumn(genP_index, auxTools_.ConvertIntVectorToString(p.mothersIndex() ) );
	}
      else info.AddRowColumn(genP_index, ".. Too many .." );
      if (p.daughtersIndex().size() < 6)
	{
	  info.AddRowColumn(genP_index, auxTools_.ConvertIntVectorToString(p.daughtersIndex()) );
	}
      else info.AddRowColumn(genP_index, ".. Too many .." );

      // p.PrintProperties();
      // p.PrintDaughters();
      // p.PrintFinalDaughters();
      // p.PrintMothers();
    }

 
  if (bPrintList) info.Print();
  
  return theGenParticles;
}


//============================================================================
vector<GenParticle> TreeReaderMC::GetGenParticles(int pdgId, bool isLastCopy)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetGenParticles()" << std::endl;
  
  // Initialise variables
  vector<GenParticle> myGenParticles;

  // For-lool: All GenParticles
  for (unsigned int genP_index = 0; genP_index < GenP_Pt.size(); genP_index++)
    {

      // Only examine particles of specific pdgId
      if ( abs(GenP_PdgId.at(genP_index) ) != pdgId) continue;

      // Get the genParticle
      GenParticle p = GetGenParticle(genP_index);
      SetGenParticleMomsAndDads(p);
		
      if (!isLastCopy) myGenParticles.push_back(p);
      else
	{
	  
	  // Determine if it's a last copy
	  bool save = true;
	  bool bDecaysToSelf = false;
	  for (unsigned int i = 0; i < p.daughtersIndex().size(); i++) 
	    {

	      GenParticle d = GetGenParticle(p.daughtersIndex().at(i) );
	      if (abs(d.pdgId()) == pdgId)
		{
		  bDecaysToSelf = true;
		  break;
		}
	    }

	  if (bDecaysToSelf) continue;	  
	  myGenParticles.push_back(p);
	  
	  // Print properties?
	  if (0) p.PrintDaughters();
	  
	}
    }

  return myGenParticles;
}


//============================================================================
GenParticle TreeReaderMC::GetGenParticle(unsigned int Index)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetGenParticle()" << std::endl;

  // Get the GenParticle properties
  double Pt      = GenP_Pt.at(Index);
  double Eta     = GenP_Eta.at(Index);
  double Phi     = GenP_Phi.at(Index);
  double Mass    = GenP_Mass.at(Index);
  int Charge     = GenP_Charge.at(Index);
  int PdgId      = GenP_PdgId.at(Index);
  int Status     = GenP_Status.at(Index);
  double VertexX = GenP_VertexX.at(Index);
  double VertexY = GenP_VertexY.at(Index);
  double VertexZ = GenP_VertexZ.at(Index);
  vector<unsigned short> MothersIndex   = GenP_Mothers.at(Index);
  vector<unsigned short> DaughtersIndex = GenP_Daughters.at(Index);

  // Construct the GenParticle
  GenParticle theGenParticle(Index, Pt, Eta, Phi, Mass, Charge, PdgId, Status, VertexX, VertexY, VertexZ, MothersIndex, DaughtersIndex);
  
  if (cfg_DEBUG*0) theGenParticle.PrintProperties();
  
  return theGenParticle;
}



//============================================================================
void TreeReaderMC::SetGenParticleMomsAndDads(GenParticle &p)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::SetGenParticleMomsAndDads()" << std::endl;
  
  // Variable declaration
  vector<GenParticle> Mothers;
  vector<GenParticle> Daughters;

  // For-loops:
  for (unsigned int i = 0; i < p.mothersIndex().size();   i++) Mothers.push_back  ( GetGenParticle( p.mothersIndex().at(i)   ) );
  for (unsigned int j = 0; j < p.daughtersIndex().size(); j++) Daughters.push_back( GetGenParticle( p.daughtersIndex().at(j) ) );
  p.SetMothers(Mothers);
  p.SetDaughters(Daughters);
  SetGenParticleFinalDaughters(p);  
  
  return;
}


//============================================================================
vector<GenParticle> TreeReaderMC::GetHadronicGenTaus(vector<GenParticle> GenTaus,
						       double visEt,
						       double visEta)
//============================================================================
{
  // Sanity check
  vector<GenParticle> genP_hadGenTaus;
  if (GenTaus.size() < 1 ) return genP_hadGenTaus;
    
  // For-loop: GenTaus
  for (vector<GenParticle>::iterator tau = GenTaus.begin(); tau != GenTaus.end(); tau++)
    {

      // If no hadronic decay products found (pi+/-, pi0, K+/-, K0, K0L), skip this tau
      bool bIsHadronicTau = (tau->finalDaughtersCharged().size() > 0);

      if (0) tau->PrintFinalDaughters();
      
      // Apply acceptance cuts
      bool bPassVisEt  = ( tau->p4vis().Et() >= visEt);
      bool bPassVisEta = ( abs(tau->p4vis().Eta()) <= visEta );


      if (!(bPassVisEt)) continue;
      if (!(bPassVisEta)) continue;
      if (!bIsHadronicTau) continue;
      
      // Save this hadronic generator tau
      genP_hadGenTaus.push_back(*tau);
    }
  
  return genP_hadGenTaus;
}      


//============================================================================
void TreeReaderMC::SetGenParticleFinalDaughters(GenParticle &p)
//============================================================================
{

  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::SetGenParticleFinalDaughters()" << std::endl;
  
  // Variable declaration
  vector<GenParticle> finalDaughters;
  vector<unsigned short> finalDaughtersIndex;

  // For taus, get hadronic decay products (pi+/-, pi0, K+/-, K0, K0L)
  GetHadronicTauFinalDaughters(p, finalDaughtersIndex);

  // For-loop: All daughter indices
  for (unsigned int i = 0; i < finalDaughtersIndex.size(); i++)
    {
      GenParticle d = GetGenParticle( finalDaughtersIndex.at(i) );
      SetGenParticleMomsAndDads(d);
      finalDaughters.push_back(d);
    }
  
  // Need to manually set the final daughters (not the direct ones)
  p.SetFinalDaughters(finalDaughters);

  return;
}



//============================================================================
vector<TTTrack> TreeReaderMC::GetTTTracks(const double minPt,
				      const double minEta,
				      const double maxEta,
				      const double maxChiSqRed,
				      const unsigned int minStubs,
					  //const unsigned int minStubsPS,
					  //const unsigned int maxStubsPS,
				      const unsigned nFitParams,
				      bool bPrintList)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetTTTracks()" << std::endl;

  vector<TTTrack> theTTTracks;

  // For-loop: All TTTracsk
  for (Size_t iTk = 0; iTk < L1Tks_Pt->size(); iTk++)
    { 
      if (cfg_DEBUG*0) cout << "Getting TTTrack " << iTk << endl;
      TTTrack tk = GetTTTrack(iTk, nFitParams);

      // Apply selection criteria
      if (tk.getPt() < minPt) continue;
      if (abs(tk.getEta()) > maxEta) continue;
      if (abs(tk.getEta()) < minEta) continue;
      if (tk.getChi2Red() > maxChiSqRed) continue;
      if (tk.getNumOfStubs() < minStubs) continue;
      double z0 = tk.getZ0();
      double d0 = tk.getD0();
      theTTTracks.push_back( tk );

    }

  if (bPrintList) PrintTTTrackCollection(theTTTracks);
  
  return theTTTracks;
}


//============================================================================
TTTrack TreeReaderMC::GetTTTrack(unsigned int Index,
			     const unsigned int nFitParams)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetTTTrack()" << std::endl;
  
  // Initialise variables
  TVector3 p;  

  if (cfg_DEBUG*0) cout << "Getting track properties" << endl;

  // Get the track properties
  float pt  = L1Tks_Pt->at(Index);
  float eta = L1Tks_Eta->at(Index);
  float phi = L1Tks_Phi->at(Index);
  float d0  = L1Tks_d0->at(Index);
  float z0  = L1Tks_z0->at(Index);
  p.SetPtEtaPhi(pt, eta, phi);
  //float aRInv        = L1Tks_RInv->at(Index);
  float aChi2         = L1Tks_ChiSquared->at(Index);
  bool isGenuine      = L1Tks_IsGenuine->at(Index);
  bool isUnknown      = L1Tks_IsUnknown->at(Index);
  bool isCombinatoric = L1Tks_IsCombinatoric->at(Index);
  bool isLoose        = L1Tks_IsLoose->at(Index);
  bool isFake         = L1Tks_IsFake->at(Index);
  int  nStubs         = L1Tks_NStubs->at(Index);
  int  matchTP_pdgId  = L1Tks_TP_PdgId->at(Index);
  float  matchTP_pt   = L1Tks_TP_Pt   ->at(Index);
  float  matchTP_eta  = L1Tks_TP_Eta  ->at(Index);
  float  matchTP_phi  = L1Tks_TP_Phi  ->at(Index);
  float  matchTP_z0   = L1Tks_TP_z0   ->at(Index);
  float  matchTP_dxy  = L1Tks_TP_dxy  ->at(Index);
  
  if (cfg_DEBUG*0)  cout << "Constructing the track TTTrack" << endl;
  
  // Construct the TTTrack
  TTTrack theTTTrack(Index, p, d0, z0, aChi2, isGenuine, isUnknown, isCombinatoric, isLoose, isFake, nStubs, matchTP_pdgId, matchTP_pt, matchTP_eta, matchTP_phi, matchTP_z0, matchTP_dxy, nFitParams);

  if (cfg_DEBUG*0) theTTTrack.PrintProperties();

  /* //marina
  // Get the uniquely matched TP and print its properties
  if (cfg_DEBUG*0 && matchTP_index >= 0) {
    TrackingParticle theTP;
    theTP = GetTrackingParticle(matchTP_index);
    // theTTTrack.SetTP(theTP); // cannot implement. cyclic
    theTP.PrintProperties();
  if (cfg_DEBUG*0)  cout << "" << endl;
    }
  */ //marina
  return theTTTrack;
}



//============================================================================
vector<TrackingParticle> TreeReaderMC::GetTrackingParticles(bool bPrintList)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetTrackingParticles()" << std::endl;
  
  vector<TrackingParticle> theTrackingParticles;

  // For-loop: Tracking Particles
  for (Size_t iTP = 0; iTP < TP_Pt->size(); iTP++)
    {
      TrackingParticle p = GetTrackingParticle(iTP);
      theTrackingParticles.push_back(p);
      
    }
  if (bPrintList) PrintTrackingParticleCollection(theTrackingParticles);
  
  return theTrackingParticles;
}



//============================================================================
TrackingParticle TreeReaderMC::GetTrackingParticle(unsigned int Index)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::GetTrackingParticle()" << std::endl;
  
  // Get the track properties
  float pt            = TP_Pt->at(Index);
  float eta           = TP_Eta->at(Index);
  float phi           = TP_Phi->at(Index);
  int charge          = TP_Charge->at(Index);
  int pdgId           = TP_PdgId->at(Index);
  float d0_propagated = TP_d0->at(Index);
  float z0_propagated = TP_z0->at(Index);
  int nMatch          = TP_NMatch->at(Index);
  float ttTrackPt     = TP_Trk_Pt->at(Index);
  float ttTrackEta    = TP_Trk_Eta->at(Index);
  float ttTrackPhi    = TP_Trk_Phi->at(Index);
  float ttTrackZ0     = TP_Trk_z0->at(Index);
  float ttTrackD0     = TP_Trk_d0->at(Index);
  float ttTrackChi2   = TP_Trk_ChiSquared->at(Index);
  int ttTrackNstubs   = TP_Trk_NStubs->at(Index);
  int nStubs          = TP_NStubs->at(Index);
  float z0_produced   = TP_z0_produced->at(Index);
  float d0_produced   = TP_d0_produced->at(Index);
  float dxy           = TP_dxy->at(Index);
  int eventId         = TP_EventId->at(Index);  

  // Construct higher-level variables
  TVector3 p;
  p.SetPtEtaPhi(pt, eta, phi);

  if (cfg_DEBUG*0) cout << "Constructing TP" << endl;

  // Construct the TP
  TrackingParticle theTrackingParticle(Index, p, d0_propagated, z0_propagated, d0_produced, z0_produced, dxy, charge, pdgId, nMatch, nStubs, ttTrackPt, ttTrackEta, ttTrackPhi, ttTrackZ0, ttTrackD0, ttTrackChi2, ttTrackNstubs, eventId);  
  if (cfg_DEBUG*0*0) theTrackingParticle.PrintProperties();

  // Get the uniquely matched TTTrack
  // TTTrack theTrack;
  // if (ttTrackIndex >= 0) theTrack= GetTTTrack(ttTrackIndex, 4);
  // theTrackingParticle.SetTTTTrack(theTrack); // cannot use! cyclic problems 
  // if (cfg_DEBUG*0) theTrack.PrintProperties(); // cannot use! cyclic problems
  
  return theTrackingParticle;
}


//============================================================================
void TreeReaderMC::PrintTrackingParticleCollection(vector<TrackingParticle> collection)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::PrintTrackingParticleCollection()" << std::endl;
  
  Table info("Index | Pt | Eta | Phi | PdgId | Q | z0 prop | d0 prop | z0 prod | d0 prod | dxy | NMatch | NStubs | TTTrkPt | TTTrkEta | TTTrkPhi | TTTrkZ0 | TTTrkD0 | TTTrkChi2 | TTTrkNStubs | Event-Id", "Text");

  // For-loop: Tracking Particles
  int row=0;
  for (vector<TrackingParticle>::iterator p = collection.begin(); p != collection.end(); p++)
    {
      // Construct table
      info.AddRowColumn(row, auxTools_.ToString( p->index()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Perp(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Eta() , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Phi() , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getPdgId()) );
      info.AddRowColumn(row, p->getQ());
      info.AddRowColumn(row, auxTools_.ToString( p->getZ0propagated()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getD0propagated()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getZ0produced()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getD0produced()   , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getDxy()   , 3) );
      //info.AddRowColumn(row, auxTools_.ToString( p->getD0Phi(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getNMatch())       );
      info.AddRowColumn(row, auxTools_.ToString( p->getNStubs())       );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackPt()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackEta()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackPhi()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackZ0()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackD0()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackChi2()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getTTTrackNStubs()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getEventId())     );
      row++;
    }
  

  info.Print();
  
  return;    
}


//============================================================================
void TreeReaderMC::PrintTTTrackCollection(vector<TTTrack> collection)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::PrintTTTrackCollection()" << std::endl;

  Table info("Index | Pt | Eta | Phi | z0 | d0 | Chi2 | DOF | Chi2Red | Stubs | Genuine | Unknown | Comb.", "Text");
      
  // For-loop: All TTTracsk
  int row=0;
  for (vector<TTTrack>::iterator p = collection.begin(); p != collection.end(); p++)    
    {
      // Construct table
      info.AddRowColumn(row, auxTools_.ToString( p->index()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Perp(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Eta() , 3) );  
      info.AddRowColumn(row, auxTools_.ToString( p->getMomentum().Phi() , 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getZ0(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getD0(), 3) );  
      //info.AddRowColumn(row, p->getQ());
      info.AddRowColumn(row, auxTools_.ToString( p->getChi2(),3 ));
      info.AddRowColumn(row, auxTools_.ToString( p->getDOF()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getChi2Red(), 3) );
      info.AddRowColumn(row, auxTools_.ToString( p->getNumOfStubs()) );// + " (" + auxTools_.ToString(p->getNumOfStubsPS()) + ")");
      info.AddRowColumn(row, auxTools_.ToString( p->getIsGenuine()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getIsUnknown()) );
      info.AddRowColumn(row, auxTools_.ToString( p->getIsCombinatoric()) );
      row++;
    }

  info.Print();

  return;
}

//============================================================================
void TreeReaderMC::PrintGenParticleCollection(vector<GenParticle> collection)
//============================================================================
{
  if (cfg_DEBUG*0) std::cout << "=== TreeReaderMC::PrintGenParticleCollection()" << std::endl;
  
  
  for (vector<GenParticle>::iterator p = collection.begin(); p != collection.end(); p++)
    {
      if (cfg_DEBUG*0) std::cout << "\nGenParticle:" << endl;
      p->PrintProperties();

      // For more information
      if (cfg_DEBUG*0){

	std::cout << "\nDaughters:" << endl;
	p->PrintDaughters();
	std::cout << "\nFinal Daughters:" << endl;
	p->PrintFinalDaughters();
	std::cout << "\nFinal Daughters (Charged):" << endl;
	p->PrintFinalDaughtersCharged();
	std::cout << "\nFinal Daughters (Neutral):" << endl;
	p->PrintFinalDaughtersNeutral();
	std::cout << "\nMothers:" << endl;
	p->PrintMothers();
      }
      
    }
  
  return;
}


//============================================================================
vector<L1Tau> TreeReaderMC::GetL1Taus(bool bPrintList)
//============================================================================
{
  vector<L1Tau> theL1Taus;
  L1Tau theL1Tau;

  // For-loop: All L1Taus
  for (Size_t iCalo = 0; iCalo < L1TauEmu_Et.size(); iCalo++)
    { 
      theL1Tau = GetL1Tau(iCalo);
      //if (bPrintList) theL1Tau.PrintProperties(false);
      theL1Taus.push_back( theL1Tau );
    }
  
  if (bPrintList) PrintL1TauCollection(theL1Taus); 
  return theL1Taus;
}


//============================================================================
L1Tau TreeReaderMC::GetL1Tau(unsigned int Index)
//============================================================================
{

  if (0)
    {
      std::cout << " 1: " << L1TauEmu_Et.at(Index) << std::endl;
      std::cout << " 2: " << L1TauEmu_Eta.at(Index) << std::endl;
      std::cout << " 3: " << L1TauEmu_Phi.at(Index) << std::endl;
      std::cout << " 4: " << L1TauEmu_IEt.at(Index) << std::endl;
      std::cout << " 5: " << L1TauEmu_IEta.at(Index) << std::endl;
      std::cout << " 6: " << L1TauEmu_IPhi.at(Index) << std::endl;
      std::cout << " 7: " << L1TauEmu_Iso.at(Index) << std::endl;
      std::cout << " 8: " << L1TauEmu_Bx.at(Index) << std::endl;
      std::cout << " 9: " << L1TauEmu_TowerIPhi.at(Index) << std::endl;
      std::cout << "10: " << L1TauEmu_TowerIEta.at(Index) << std::endl;
      std::cout << "11: " << L1TauEmu_RawEt.at(Index) << std::endl;
      std::cout << "12: " << L1TauEmu_IsoEt.at(Index) << std::endl;
      std::cout << "13: " << L1TauEmu_NTT.at(Index) << std::endl;
      std::cout << "14: " << L1TauEmu_HasEM.at(Index) << std::endl;
      std::cout << "15: " << L1TauEmu_IsMerged.at(Index) << std::endl;
      std::cout << "16: " << L1TauEmu_HwQual.at(Index) << std::endl;
    }
  
  L1Tau theL1Tau(Index,
			 L1TauEmu_Et.at(Index),
			 L1TauEmu_Eta.at(Index),
			 L1TauEmu_Phi.at(Index),
			 L1TauEmu_IEt.at(Index),
			 L1TauEmu_IEta.at(Index),
			 L1TauEmu_IPhi.at(Index),
			 L1TauEmu_Iso.at(Index),
		         L1TauEmu_Bx.at(Index),
			 L1TauEmu_TowerIPhi.at(Index),
			 L1TauEmu_TowerIEta.at(Index),
			 L1TauEmu_RawEt.at(Index),
			 L1TauEmu_IsoEt.at(Index),
			 L1TauEmu_NTT.at(Index),
			 L1TauEmu_HasEM.at(Index),
			 L1TauEmu_IsMerged.at(Index),
			 L1TauEmu_HwQual.at(Index)
			 );

  return theL1Tau;
}



//============================================================================
vector<L1Jet> TreeReaderMC::GetL1Jets(bool bPrintList)
//============================================================================
{

  vector<L1Jet> theL1Jets;
  L1Jet theL1Jet;

  // For-loop: All L1 Jets
  for (Size_t i = 0; i < L1JetEmu_Et.size(); i++)
    {
      theL1Jet = GetL1Jet(i);
      theL1Jets.push_back(theL1Jet);
    }
  
  if (bPrintList) PrintL1JetCollection(theL1Jets); 
  return theL1Jets;
}


//============================================================================
L1Jet TreeReaderMC::GetL1Jet(unsigned int Index)
//============================================================================
{

  L1Jet theL1Jet(Index,
		 L1JetEmu_Et.at(Index),
		 L1JetEmu_Eta.at(Index),
		 L1JetEmu_Phi.at(Index),
		 L1JetEmu_IEt.at(Index),
		 L1JetEmu_IEta.at(Index),
		 L1JetEmu_IPhi.at(Index),
		 L1JetEmu_Bx.at(Index),
		 L1JetEmu_RawEt.at(Index),
		 L1JetEmu_SeedEt.at(Index),
		 L1JetEmu_TowerIEta.at(Index),
		 L1JetEmu_TowerIPhi.at(Index),
		 L1JetEmu_PUEt.at(Index),
		 L1JetEmu_PUDonutEt0.at(Index),
		 L1JetEmu_PUDonutEt1.at(Index),
		 L1JetEmu_PUDonutEt2.at(Index),
		 L1JetEmu_PUDonutEt3.at(Index));

  return theL1Jet;
}


//============================================================================
vector<L1EG> TreeReaderMC::GetL1EGs(bool bPrintList)
//============================================================================
{
  vector<L1EG> theL1EGs;
  L1EG theL1EG;

  // For-loop: All L1 EG
  for (Size_t i = 0; i < L1EGEmu_Et.size(); i++)
    {
      theL1EG = GetL1EG(i);
      theL1EGs.push_back(theL1EG);
    }
  
  if (bPrintList) PrintL1EGCollection(theL1EGs); 
  return theL1EGs;
}

//============================================================================
L1EG TreeReaderMC::GetL1EG(unsigned int Index)
//============================================================================
{
  
  L1EG theL1EG(Index,
	       L1EGEmu_Et.at(Index),
	       L1EGEmu_Eta.at(Index),
	       L1EGEmu_Phi.at(Index),
	       L1EGEmu_IEt.at(Index),
	       L1EGEmu_IEta.at(Index),
	       L1EGEmu_IPhi.at(Index),
	       L1EGEmu_Iso.at(Index),
	       L1EGEmu_Bx.at(Index),
	       L1EGEmu_TowerIPhi.at(Index),
	       L1EGEmu_TowerIEta.at(Index),
	       L1EGEmu_RawEt.at(Index),
	       L1EGEmu_IsoEt.at(Index),
	       L1EGEmu_FootprintEt.at(Index),
	       L1EGEmu_NTT.at(Index),
	       L1EGEmu_Shape.at(Index),
	       L1EGEmu_TowerHoE.at(Index));
    
    return theL1EG;
}

//============================================================================
L1Sum TreeReaderMC::GetL1Sum(unsigned int Index)
//============================================================================
{

  L1Sum theL1Sum(Index,
		 L1SumEmu_Et.at(Index),
		 L1SumEmu_Phi.at(Index),
		 L1SumEmu_IEt.at(Index),
		 L1SumEmu_IPhi.at(Index),
		 L1SumEmu_Type.at(Index),
		 L1SumEmu_Bx.at(Index));

  return theL1Sum;
}


//============================================================================
vector<L1Sum> TreeReaderMC::GetL1Sums(bool bPrintList)
//============================================================================
{

  vector<L1Sum> theL1Sums;
  L1Sum theL1Sum;

  // For-loop: All L1 EG
  for (Size_t i = 0; i < L1SumEmu_Et.size(); i++)
    {
      theL1Sum = GetL1Sum(i);
      theL1Sums.push_back(theL1Sum);
    }
  
  if (bPrintList) PrintL1SumCollection(theL1Sums); 
  return theL1Sums;
}

//============================================================================
void TreeReaderMC::PrintL1TauCollection(vector<L1Tau> collection)
//============================================================================
{
  
  Table info("Index | Et | Eta | Phi | IET | IEta | IPhi | Iso | Bx | TowerIPhi | TowerIEta | RawEt | IsoEt | NTT | HasEM | IsMerged | HwQual | Type", "Text");

  // For-loop: All L1Taus
  int row=0;
  for (auto p = collection.begin(); p != collection.end(); p++)
  {
    // Construct table
    info.AddRowColumn(row, auxTools_.ToString( p->getIndex(), 1)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getEt(), 4)        );
    info.AddRowColumn(row, auxTools_.ToString( p->getEta(), 4)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getPhi(), 4)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getIEt() , 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getIEta(), 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getIPhi(), 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getIso() , 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getBx()  , 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getTowerIPhi(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getTowerIEta(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getRawEt(), 3)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getIsoEt(), 3)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getNTT(), 3)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getHasEM(), 3)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getIsMerged(), 3)  );
    info.AddRowColumn(row, auxTools_.ToString( p->getHwQual(), 3)    );
    row++;
    
  }

  if (collection.size()>0) info.Print();
  return;
}

//============================================================================
void TreeReaderMC::PrintL1JetCollection(vector<L1Jet> collection)
//============================================================================
{
  
  Table info("Index | Et | Eta | Phi | IET | IEta | IPhi | Bx | Raw Et | Seed Et | TowerIPhi | TowerIEta | PU Et | PUDonutEt0 | PUDonutEt1 | PUDonutEt2 | PUDonutEt3", "Text");

  // For-loop: All L1Taus
  int row=0;
  for (auto p = collection.begin(); p != collection.end(); p++, row++)
  {
    // Construct table
    info.AddRowColumn(row, auxTools_.ToString( p->getIndex() )        );
    info.AddRowColumn(row, auxTools_.ToString( p->getEt(), 3)         );
    info.AddRowColumn(row, auxTools_.ToString( p->getEta(), 3)        );
    info.AddRowColumn(row, auxTools_.ToString( p->getPhi(), 3)        );
    info.AddRowColumn(row, auxTools_.ToString( p->getIEt() , 3)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getIEta(), 3)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getIPhi(), 3)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getBx() , 3)        );
    info.AddRowColumn(row, auxTools_.ToString( p->getRawEt(), 3)      );
    info.AddRowColumn(row, auxTools_.ToString( p->getSeedEt(), 3)     );
    info.AddRowColumn(row, auxTools_.ToString( p->getTowerIPhi(), 3)  );
    info.AddRowColumn(row, auxTools_.ToString( p->getTowerIEta(), 3)  );
    info.AddRowColumn(row, auxTools_.ToString( p->getPUEt(), 3)       );
    info.AddRowColumn(row, auxTools_.ToString( p->getPUDonutEt0(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getPUDonutEt1(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getPUDonutEt2(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getPUDonutEt3(), 3) );
  }

  info.Print();
  return;
}


//============================================================================
void TreeReaderMC::PrintL1EGCollection(vector<L1EG> collection)
//============================================================================
{
  
  Table info("Index | Et | Eta | Phi | IET | IEta | IPhi | Iso | Bx | TowerIPhi | TowerIEta | Raw Et | Iso Et| Footprint Et | NTT | Shape | TowerHoE", "Text");
  
  // For-loop: All L1EG
  int row=0;
  for (auto p = collection.begin(); p != collection.end(); p++, row++)
  {
    // Construct table
  info.AddRowColumn(row, auxTools_.ToString( p->getIndex() ) );
  info.AddRowColumn(row, auxTools_.ToString( p->getEt()         , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getEta()        , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getPhi()        , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getIEt()        , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getIEta()       , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getIPhi()       , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getIso()        , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getBx()         , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getTowerIPhi()  , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getTowerIEta()  , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getRawEt()      , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getIsoEt()      , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getFootprintEt(), 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getNTT()        , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getShape()      , 3) );
  info.AddRowColumn(row, auxTools_.ToString( p->getTowerHoE()   , 3) );
  
  }

  info.Print();
  return;
}


//============================================================================
void TreeReaderMC::PrintL1SumCollection(vector<L1Sum> collection)
//============================================================================
{
  
  Table info("Index | Et | Phi | IEt | IPhi | Type | Bx", "Text");
  
  // For-loop: All L1Sum
  int row=0;
  for (auto p = collection.begin(); p != collection.end(); p++, row++)
  {
    // Construct table
    info.AddRowColumn(row, auxTools_.ToString( p->getIndex() ) );
    info.AddRowColumn(row, auxTools_.ToString( p->getEt()  , 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getPhi() , 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getIEt() , 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getIPhi(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getType(), 3) );
    info.AddRowColumn(row, auxTools_.ToString( p->getBx()  , 3) );
  }

  info.Print();
  return;
}


//============================================================================
void TreeReaderMC::PrintL1TkTauParticleCollection(vector<L1TkTauParticle> collection)
//============================================================================
{
  
  Table info("Match-Cone | Sig-Cone | Iso-Cone | Calo-Et | Calo-Eta | Tk-Pt | Tk-Eta | Tk-dR | Gen-Pt | Gen-Eta | Gen-dR | Sig-Tks | Iso-Tks | VtxIso | RelIso", "Text");
  
  // For-loop: All L1TkTauParticles
  int row=0;
  for (vector<L1TkTauParticle>::iterator p = collection.begin(); p != collection.end(); p++)
  {
    // Construct table
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchConeMin(), 2) + " < dR < " + auxTools_.ToString( p->GetMatchConeMax(), 2) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetSigConeMin()  , 2) + " < dR < " + auxTools_.ToString( p->GetSigConeMax()  , 2) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetIsoConeMin()  , 2) + " < dR < " + auxTools_.ToString( p->GetIsoConeMax()  , 2) ); 
    info.AddRowColumn(0, auxTools_.ToString( p->GetCaloTau().et()       ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetCaloTau().eta()      ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingTk().getPt() ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingTk().getEta(),3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingTkDeltaR()   ,3  ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingGenParticle().pt() , 3) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingGenParticle().eta(), 3) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetMatchingGenParticleDeltaR(), 3) );  
    info.AddRowColumn(0, auxTools_.ToString( p->GetSigConeTTTracks().size() ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetIsoConeTTTracks().size() ) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetVtxIsolation(), 3) );
    info.AddRowColumn(0, auxTools_.ToString( p->GetRelIsolation(), 3) );
    row++;
  }

  if (collection.size()>0) info.Print();
  return;
}

#endif //TreeReaderMC_cxx 

#ifndef GenParticle_cxx
#define GenParticle_cxx

// User
#include "../interface/GenParticle.h"

//****************************************************************************
GenParticle::GenParticle()
//****************************************************************************
{

}

//****************************************************************************
GenParticle::GenParticle(unsigned short Index,
			 double Pt,
			 double Eta,
			 double Phi,
			 double Mass,
			 int Charge,
			 int PdgId,
			 int Status,
			 double VertexX,
			 double VertexY,
			 double VertexZ,
			 vector<unsigned short> MothersIndex,
			 vector<unsigned short> DaughtersIndex)
//****************************************************************************
{

  theIndex     = Index;
  thePt        = Pt;
  theEta       = Eta;
  thePhi       = Phi;
  theMass      = Mass;
  theCharge    = Charge;
  thePdgId     = PdgId;
  theStatus    = Status;
  theVertexX   = VertexX;
  theVertexY   = VertexY;
  theVertexZ   = VertexZ;
  theMothersIndex   = MothersIndex;
  theDaughtersIndex = DaughtersIndex;
  theP4.SetPtEtaPhiM(thePt, theEta, thePhi, theMass);
  theVertex.SetXYZ(VertexX, VertexY, VertexZ);
  if (0) PrintProperties();
  
}


//****************************************************************************
GenParticle::~GenParticle()
//****************************************************************************
{

}



//****************************************************************************
void GenParticle::PrintProperties(bool bPrintTitleRows)
//****************************************************************************
{
  
  Table info("Index | Pt | Eta | Phi | Mass | Charge | PdgId | Status | vX | vY | vZ | mothers | daughters |", "Text");
  
  info.AddRowColumn(0, auxTools.ToString( index() , 1) );
  info.AddRowColumn(0, auxTools.ToString( pt()    , 3) );
  info.AddRowColumn(0, auxTools.ToString( eta()   , 3) );
  info.AddRowColumn(0, auxTools.ToString( phi()   , 3) );
  info.AddRowColumn(0, auxTools.ToString( mass()  , 3) );
  info.AddRowColumn(0, auxTools.ToString( charge(), 1) );
  info.AddRowColumn(0, auxTools.ToString( pdgId() , 3) );
  info.AddRowColumn(0, auxTools.ToString( status(), 1) );
  info.AddRowColumn(0, auxTools.ToString( vx()    , 3) );
  info.AddRowColumn(0, auxTools.ToString( vy()    , 3) );
  info.AddRowColumn(0, auxTools.ToString( vz()    , 3) );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(mothersIndex()) );
  info.AddRowColumn(0, auxTools.ConvertIntVectorToString(daughtersIndex()) );
  info.Print(bPrintTitleRows);

  return;
}


//****************************************************************************
void GenParticle::PrintDaughters(bool bPrintTitleRows)
//****************************************************************************
{
  PrintDaughters_(theDaughters, bPrintTitleRows);
  return;
}


//****************************************************************************
void GenParticle::PrintFinalDaughters(bool bPrintTitleRows)
//****************************************************************************
{
  PrintDaughters_(theFinalDaughters, bPrintTitleRows);
  return;
}

//****************************************************************************
void GenParticle::PrintFinalDaughtersCharged(bool bPrintTitleRows)
//****************************************************************************
{
  PrintDaughters_(theFinalDaughtersCharged, bPrintTitleRows);
  return;
}

//****************************************************************************
void GenParticle::PrintFinalDaughtersNeutral(bool bPrintTitleRows)
//****************************************************************************
{
  PrintDaughters_(theFinalDaughtersNeutral, bPrintTitleRows);
  return;
}


//****************************************************************************
void GenParticle::PrintDaughters_(vector<GenParticle> daughters, bool bPrintTitleRows)
//****************************************************************************
{

  Table info("Index | Pt | Eta | Phi | Mass | Charge | PdgId | Status | vX | vY | vZ | mothers | daughters |", "Text");
  int row = 0;  
  for (vector<GenParticle>::iterator d = daughters.begin(); d != daughters.end(); d++)
    {
      // d->PrintProperties(bPrintTitleRows);
      info.AddRowColumn(row, auxTools.ToString( d->index() , 1) );
      info.AddRowColumn(row, auxTools.ToString( d->pt()    , 3) );
      info.AddRowColumn(row, auxTools.ToString( d->eta()   , 3) );
      info.AddRowColumn(row, auxTools.ToString( d->phi()   , 3) );
      info.AddRowColumn(row, auxTools.ToString( d->mass()  , 3) );
      info.AddRowColumn(row, auxTools.ToString( d->charge(), 1) );
      info.AddRowColumn(row, auxTools.ToString( d->pdgId() , 3) );
      info.AddRowColumn(row, auxTools.ToString( d->status(), 1) );
      info.AddRowColumn(row, auxTools.ToString( d->vx()    , 3) );
      info.AddRowColumn(row, auxTools.ToString( d->vy()    , 3) );
      info.AddRowColumn(row, auxTools.ToString( d->vz()    , 3) );
      info.AddRowColumn(row, auxTools.ConvertIntVectorToString(d->mothersIndex()) );
      info.AddRowColumn(row, auxTools.ConvertIntVectorToString(d->daughtersIndex()) );
      row++;
    }
  info.Print(bPrintTitleRows);
  
  return;
}


//****************************************************************************
void GenParticle::PrintMothers(bool bPrintTitleRows)
//****************************************************************************
{
  
  for (vector<GenParticle>::iterator m = theMothers.begin(); m != theMothers.end(); m++)
    {
      m->PrintProperties(bPrintTitleRows);
    }

  return;
}


//****************************************************************************
vector<GenParticle> GenParticle::mothers(void)
//****************************************************************************
{

  // if (theMothers.size() < 1)
  //   {
  //     std::cout << "=== GenParticle::mothers() - theMothers.size() = " << theMothers.size() << ". Have you called GenParticle::SetMothers()?" << std::endl;
  //   }
  return theMothers;
}


//****************************************************************************
vector<GenParticle> GenParticle::daughters(void)
//****************************************************************************
{

  // if (theDaughters.size() < 1)
  //   {
  //     std::cout << "=== GenParticle::daughters() - theDaughters.size() = " << theDaughters.size() << ". Have you called GenParticle::SetDaughters()?" << std::endl;
  //   }
  return theDaughters;
}



//****************************************************************************
void GenParticle::SetFinalDaughtersCharged(void)
//****************************************************************************
{

  for (vector<GenParticle>::iterator d = theFinalDaughters.begin(); d != theFinalDaughters.end(); d++)
    {
      if (d->charge() == 0) continue;
      theFinalDaughtersCharged.push_back(*d);
    }
  
  return;
}


//****************************************************************************
void GenParticle::SetFinalDaughtersNeutral(void)
//****************************************************************************
{
  std::vector<unsigned int> nuIds;
  nuIds.push_back(12);  // nu_e
  nuIds.push_back(14);  // nu_mu
  nuIds.push_back(16);  // nu_tau
  
  for (vector<GenParticle>::iterator d = theFinalDaughters.begin(); d != theFinalDaughters.end(); d++)
    {
      if (d->charge() != 0) continue;
      if ( std::find(nuIds.begin(), nuIds.end(), abs(d->pdgId()) ) != nuIds.end() ) continue;
      theFinalDaughtersNeutral.push_back(*d);
    }
  
  return;
}


//****************************************************************************
TLorentzVector GenParticle::p4vis(void)
//****************************************************************************
{

  TLorentzVector p4vis;
  std::vector<unsigned int> nuIds;
  nuIds.push_back(12);  // nu_e
  nuIds.push_back(14);  // nu_mu
  nuIds.push_back(16);  // nu_tau


  for (vector<GenParticle>::iterator d = theFinalDaughters.begin(); d != theFinalDaughters.end(); d++)
    {

      if ( std::find(nuIds.begin(), nuIds.end(), abs(d->pdgId()) ) != nuIds.end() ) continue;

      p4vis += d->p4();
    }

  return p4vis;
  
}

#endif

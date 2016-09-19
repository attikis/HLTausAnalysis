#ifndef L1TkPrimaryVertex_cxx
#define L1TkPrimaryVertex_cxx

// System
#include <iostream>
#include <algorithm>
#include <cmath>

// User
#include "../../Auxiliary/interface/constants.h"
#include "../interface/L1TkPrimaryVertex.h"

// ROOT
#include <TLorentzVector.h>

#ifdef USING_MC
#include "../../Framework/interface/TreeDefinitionGenP.h"
#endif

//#define DEBUG
//****************************************************************************
double L1TkPrimaryVertex::GetPrimaryVertexZ(string algo, 
					    vector<int> &PVTks_index) 
//****************************************************************************
{

  double pv_z = -9999.9;
  
  // To be called from main, selects which PV-algorithm is called
  if(algo.compare("slowPT")==0)       pv_z = SlowVtxProducer(algo, PVTks_index);
  else if(algo.compare("slowPT2")==0) pv_z = SlowVtxProducer(algo, PVTks_index);
  else if(algo.compare("fast")==0)    pv_z = FastVtxProducer(algo, PVTks_index);
  else if(algo.compare("fastPT")==0)  pv_z = FastVtxProducer(algo, PVTks_index);
  else if(algo.compare("fastPT2")==0) pv_z = FastVtxProducer(algo, PVTks_index);
  else if(algo.compare("mc")==0)      pv_z = s->t->HepMCEvt_VtxZ;
  else if(algo.compare("default")==0) pv_z = 0.0;
  else{
    cout << "E R R O R ! L1TkPrimaryVertex::TkBasedPVProducer(...) - Unknown PV-producer algorithm \"" << algo << "\". EXIT" << endl;
    exit(1);
  }

  return pv_z;
}


//****************************************************************************
double L1TkPrimaryVertex::SlowVtxProducer(string algo, 
				       vector<int> &PVTks_index)
//****************************************************************************
{

  // Algorithm loops over all tracks for each bin and calculates sum
  Double_t sum_ptMax =  -999.9; 
  Double_t pv_z      = +1000.0;
  Double_t zmax      =   +25.0;
  Int_t    nBins     = zmax * 10 * 2; //  bin size (in mm)
  
  // For-loop: z-axis bins
  for(Int_t iBin = 0; iBin <= nBins; iBin++)
    {
      
      // Select one bin (1 mm accross) at a time and calculate the sum pT of associated tracks
      Double_t zCoordinate = -zmax * 10 + iBin;  // convert to mm 
      zCoordinate          = zCoordinate / 10.0; // convert to cm
      Double_t delta_z     = 0.1;                // bin-witdh to look for tracks  

      vector<int> tmpPvTks_index;
      Double_t tmp_sum_pt = VtxProducerSumPtBin(zCoordinate, delta_z, algo, tmpPvTks_index);

      // Look for the peak
      if( tmp_sum_pt > sum_ptMax)
	{
	  
	  // Update PV z-coordinate and associated sum_pT  
	  sum_ptMax = tmp_sum_pt;
	  pv_z   = zCoordinate;
	  // Update the vector with the tracks associated with selected PV 
	  PVTks_index.erase( PVTks_index.begin(), PVTks_index.end() );
	  PVTks_index.insert(PVTks_index.begin(), tmpPvTks_index.begin(), tmpPvTks_index.end());

	}

    }//For-loop: Bins

  return pv_z;
}


//****************************************************************************
 double L1TkPrimaryVertex::VtxProducerSumPtBin(const double zCoordinate, 
					    const double delta_z, 
					    const string algo,
					    vector<int> &pvTks_index)
//****************************************************************************
{

  double sum_pt         = 0.0;
  unsigned int nTks     = s->t->L1Tks_Pt->size();  
  const double tk_ptMax = 50.0;

  // For-loop: All Tracks
  for(unsigned int iTk = 0; iTk < nTks; iTk++)
    {

      // Does the track pass quality criteria?
      if ( !s->SelectTracks(iTk, "VtxProducer") ) continue;
      
      const double tk_pocaz = s->t->L1Tks_POCAz->at(iTk);
      double tk_pt    = s->t->L1Tks_Pt->at(iTk);

      // Is this track within the zCoordinate of the bin of interest?
      if( abs(tk_pocaz - zCoordinate) > delta_z) continue;
      
      // Saturate the track pT       
      if(tk_pt > tk_ptMax) tk_pt = tk_ptMax;
      
      // Save index of all tracks coming from the calculate PV
      pvTks_index.push_back(iTk);

      // Calculate the sum_pT according to algorithm selection
      if(algo.compare("slowPT") == 0) sum_pt += tk_pt;
      else if(algo.compare("slowPT2") == 0) sum_pt +=(tk_pt*tk_pt);
      else{}

    } // For-loop: All Tracks

  return sum_pt;

}


//****************************************************************************
double L1TkPrimaryVertex::FastVtxProducer(string algo,
				       vector<int> &PVTks_index)
//****************************************************************************
{

  // Declare variables
  vector <int> selTks_index;
  double zcenter = 0.0;
  const double tk_ptMax = 50.0;

  // Create histogram to store POC-z of all tracks
  const double xMin  = - 25.0;
  const double xMax  = + 25.0;
  const int  nBins   = +500;
  TH1D *htemp;
  htemp = new TH1D("htemp", "htemp", nBins, xMin, xMax);
  
  // For-loop: All Tks
  for(unsigned int iTk = 0; iTk < s->t->L1Tks_Pt->size(); iTk++)
    {

      // Does the track pass quality criteria?
      if ( !s->SelectTracks(iTk, "VtxProducer") ) continue;

      // Get track variables
      double weight   = 0.0;
      double tk_pt    = s->t->L1Tks_Pt->at(iTk);
      double tk_pocaz = s->t->L1Tks_POCAz->at(iTk);

      // Saturate the track pT
      if(tk_pt > tk_ptMax) tk_pt = tk_ptMax;

      // Weight histogram entries according to algorithm selection
      if(algo.compare("fast")==0)   weight = 1.0;
      else if(algo.compare("fastPT")==0) weight = tk_pt;
      else if(algo.compare("fastPT2")==0) weight = tk_pt*tk_pt;
      else{
	cout << "E R R O R ! L1TkPrimaryVertex::FastVtxProducer(...) - Unknown algorithm type \"" << algo << "\". EXIT" << endl;
	exit(1);
      }

      // Fill the histogram with the POCA-z, weighted accordingly
      htemp->Fill(tk_pocaz, weight);

      // Save index of all selected tracks
      selTks_index.push_back(iTk);
    }
  
  // Get the PV z-coordinate by determining the peak of the histo
  Double_t pv_z = -1000.0;
  pv_z = PeakFinder(htemp, zcenter);

  // For-loop: All Selected tracks
  for(unsigned int k=0; k < selTks_index.size();k++)
    {

      Int_t tkindx    = selTks_index.at(k);
      Double_t tmptkz = s->t->L1Tks_POCAz->at(tkindx);
      
      // Select only tracks withing 3 mm of PV z-coordinate center (3 bins * 1mm wide each)
      if(tmptkz >= zcenter-0.15 && tmptkz <= zcenter+0.15) PVTks_index.push_back(tkindx);

    }// For-loop: All Selected tracks

  // Delete the no longer needed histogram and return the calculated PV-z
  delete htemp;

  return pv_z;
  
}

//****************************************************************************
double L1TkPrimaryVertex::PeakFinder(TH1D *htemp,
				  Double_t &center)
//****************************************************************************
{
  
  // Declare variables
  double sum_max = -1000.0;
  double pv_z    = -1000.0;
  const unsigned int nBins = htemp->GetNbinsX();

  // For-loop: All Bins in groups of 3
  for(unsigned int i=2; i <= nBins-1; i++)
    {

      // Get bin contents
      double a0 = htemp->GetBinContent(i-1);
      double a1 = htemp->GetBinContent(i);
      double a2 = htemp->GetBinContent(i+1);

      // Sum bin contents 
      double threeBinSum = a0+a1+a2;

      // Find 3-bin weigthed maximum
      if(threeBinSum > sum_max)
	{

	  // Update maximum max
	  sum_max = threeBinSum;
	  
	  // Calculate the pv_z as the weighted average
	  double z0 = htemp -> GetBinCenter(i-1);
	  double z1 = htemp -> GetBinCenter(i);
	  double z2 = htemp -> GetBinCenter(i+1);
	  pv_z = (a0*z0 + a1*z1 + a2*z2)/threeBinSum;
	  center = z1;
	}
    } // For-loop: All Bins in groups of 3

  return pv_z;

}

#endif


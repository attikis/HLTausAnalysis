#ifndef AuxTools_cxx
#define AuxTools_cxx

// System
#include <string>
#include <sstream>
#include <iostream>

// User
#include "../interface/AuxTools.h"


// #define myDEBUG

//****************************************************************************
void AuxTools::ProgressBar(Long64_t entry,
			   Long64_t total, 
			   Int_t resolution,
			   Int_t barWidth)
//****************************************************************************  
{  

  // Correct for loop index
  entry = entry + 1;

  // Sanity check (otherwise crashes)
  if (total < resolution) resolution = total;

  // Only update "resolution" times.
  if ( entry % (total/resolution) != 0 ) return;

  // Calculuate the ratio of complete-to-incomplete.
  Double_t ratio      = (Double_t) entry / (Double_t) total;
  Double_t percentage = 100.0*ratio;
  Int_t completed     = ratio*barWidth;

  // Show the percentage completed
  cout << "Progress: " << setprecision(3) << percentage << " % " << " (" << entry << "/" << total << ")";
  
  // Show the progress bar
  // for (Int_t i=0; i < completed; i++) cout << "=";
  for (Int_t j=completed; j < barWidth; j++) cout << " ";
  
  // Move to the first column and flush
  cout << "\r";
  fflush(stdout);
  
  return;
}


//****************************************************************************
long long AuxTools::nCr(int n,
			int r)
//****************************************************************************
{

  if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++) {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;  
}
 

//****************************************************************************
template<class TYPE> float AuxTools::Sgn(TYPE myNumber)
//****************************************************************************
{

  if (myNumber >= 0) return +1.0;
  else return -1.0;
  
}
//****************************************************************************



//****************************************************************************
template <typename T>  std::string AuxTools::ToString(const T a_value, const int n)
//****************************************************************************
{
  // comment: equivalent to "to_string" but 
  // it is only available in new gcc compiler versions
  std::ostringstream out;
  out << std::setprecision(n) << a_value;

  return out.str();
}


//****************************************************************************  
Double_t AuxTools::DeltaPhi(const Double_t phi1, 
			    const Double_t phi2)
//****************************************************************************  
{
  // See: https://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_4_2/doc/html/d1/d92/DataFormats_2Math_2interface_2deltaPhi_8h_source.html
  Double_t result = phi1 - phi2;
  while (result > PI) result -= 2*PI;
  while (result <= -PI) result += 2*PI; 

  return result;
}


//****************************************************************************  
Double_t AuxTools::DeltaEta(const Double_t eta1, 
			    const Double_t eta2)
//****************************************************************************  
{
  // See: https://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_4_2/doc/html/d1/d92/DataFormats_2Math_2interface_2deltaPhi_8h_source.html
  Double_t deltaEta = fabs ( eta1 - eta2 );
  return deltaEta;
}


//****************************************************************************  
Double_t AuxTools::DeltaR(const Double_t eta1, 
			  const Double_t phi1, 
			  const Double_t eta2, 
			  const Double_t phi2)
//****************************************************************************  
{
  // See: https://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_3_9/doc/html/d5/d6b/DataFormats_2Math_2interface_2deltaR_8h_source.html
  Double_t deltaEta = DeltaEta(eta1, eta2);
  Double_t deltaPhi = DeltaPhi(phi1, phi2);
  Double_t deltaR   = sqrt( pow(deltaPhi, 2) + pow(deltaEta, 2) );
  return deltaR;
}


//****************************************************************************  
TLorentzVector AuxTools::GetTLorentzVector( Double_t pt, 
					    Double_t eta, 
					    Double_t phi, 
					    Double_t e)
//****************************************************************************  
{

  TLorentzVector v4;      // initialized by (0., 0., 0., 0.) 
  v4.SetPtEtaPhiE(pt,eta,phi,e);

  return v4;
}


//****************************************************************************  
TVector3 AuxTools::GetTVector3(Double_t px, 
			       Double_t py, 
			       Double_t pz)
//****************************************************************************  
{

  TVector3 v3(px , py , pz);

  return v3;
}


//****************************************************************************  
TVector2 AuxTools::GetTVector2(Double_t eta, 
			       Double_t phi)
//****************************************************************************  
{

  TVector2 v2(eta , phi);

  return v2;
}


//****************************************************************************  
template<class TYPE> void AuxTools::EnsureVectorIsSorted(const vector<TYPE> myVector, 
							 Bool_t bDescendingOrder)
//****************************************************************************  
{
  if (myVector.size() < 2) return;
  
  // Get value of first element
  Double_t firstVal = myVector.front();  // last elemet is: myVector.back()
  
  // For-loop: Vector Elements
  for( Size_t i = 0; i < myVector.size(); i++){

    if(bDescendingOrder){ 
      if (myVector[i] > firstVal){
	cout << "E R R O R ! AuxTools::EnsureVectorIsSorted(...) - [0] = " << myVector[0] << ", ["<< i << "] = " << myVector[i]  << endl;
	PrintVector(myVector);
	exit(0);
      }
    }
    else{ 
      if (myVector[i] < firstVal){
	cout << "E R R O R ! AuxTools::EnsureVectorIsSorted(...) - [0] = " << myVector[0] << ", ["<< i << "] = " << myVector[i]  << endl;
	PrintVector(myVector);
	exit(0);
      }
    }

  } // For-loop: Vector Elements

  return;

}

//****************************************************************************  
template<class TYPE> Bool_t AuxTools::VectorIsSorted(const vector<TYPE> myVector, 
						     Bool_t bDescendingOrder)
//****************************************************************************  
{
  if (myVector.size() < 2) return true;
  
  // Get value of first element
  Double_t firstVal = myVector.front();  // last elemet is: myVector.back()
  
  // For-loop: Vector Elements
  for( Size_t i = 0; i < myVector.size(); i++){
    
    if(bDescendingOrder){ if (myVector[i] > firstVal) return false; }
    else{ if (myVector[i] < firstVal) return false; }
    
  }// For-loop: Vector Elements
  
  return true;

}

//****************************************************************************  
template<class TYPE> void AuxTools::PrintVector(const vector<TYPE> myVector, 
						string title)
//****************************************************************************  
{
  
  if (myVector.size() < 1) return;

  cout << title << endl;
  Table table("Index | Value ", "Text", "c c");
  for(int index = 0; index < (int) myVector.size(); index++)
    {       
      table.AddRowColumn(index, ToString( index) );
      table.AddRowColumn(index, ToString(myVector[index]) );
    }

  table.Print();
  return;
}

//****************************************************************************  
template<class TYPE> string AuxTools::ConvertIntVectorToString(const vector<TYPE> myVector)
//****************************************************************************  
{
  
  if (myVector.size() < 1) return "";
  
  string text = "";
  for(int index = 0; index < (int) myVector.size(); index++)
    {       

      TYPE value = myVector[index];
      if (value == 999999) value = 0;
      if(index == (int) myVector.size()-1) text += ToString( value );
      else text += ToString( value ) + ","; 

    }

  return text;
}


//****************************************************************************  
double AuxTools::Divide(int numerator, 
			int denominator)
//****************************************************************************  
{

  if (denominator <= 0){
    cout << "W A R N I N G ! AuxTools::Divide(...) - "
              << "denominator has illegal value \"" << denominator;
    cout << "Returning 0.";
    // cout << "Exiting \n";
    // exit(1);
    return 0.0;
  }
  else return double(numerator)/double(denominator);
}


//****************************************************************************  
void AuxTools::Efficiency(int nPass, 
			  int nTotal, 
			  const string errType, 
			  double &eff, 
			  double &err )
//****************************************************************************  
{
  
  if (nTotal == 0)
    {
      eff = 0.0;
      err = 0.0;
      return;
    }
  
  
  eff = Divide(nPass, nTotal);
  if (errType.compare("binomial") == 0){
    // double errSquared = eff*(1-eff)/nTotal;
    // err = TMath::Sqrt(errSquared);
    err = (1.0/nTotal) * sqrt(nPass * (1.0 - nPass/nTotal) ); //Louise
      }
  else{
    cout << "W A R N I N G ! AuxTools::Efficiency(...) - "
	      << "Invalid error type \"" << errType << "\" selected."
	      << "Only the \"binomial\" error type is supported at the moment.";
    cout << "Exiting \n";
    exit(1);
  }
  return;
}


//****************************************************************************  
void AuxTools::StopwatchStart(void)
//****************************************************************************  
{

  stopwatch_start = clock();
  return;

}  



//****************************************************************************  
void AuxTools::PrintPSets(TTree *t)
//****************************************************************************  
{
  
  TFile *f = t->GetCurrentFile();
  TDirectory* named = (TDirectory*) f->Get("TkTauFromCaloNTupleMaker/configInfo");
  cout << "\nI N F O ! AuxTools::PrintPSets(...) - Printing paremeters used for generating these ROOT files" << endl;
  named->Get("parameterSet")->Print();
  cout << "\n" << endl;
  
  return;

}  



//****************************************************************************  
void AuxTools::StopwatchStop(const int myPrecision, 
			     const string myUnits,
			     const string myTitle)
//****************************************************************************  
{
  
  double units = 0.0;
  if (myUnits.compare("seconds") == 0)
    {
      units = 1.0;
    }
  else if (myUnits.compare("minutes") == 0)
    {
      units = 60.0;
    }
  else if (myUnits.compare("hours") == 0)
    {
      units = 3600.0;
    }
  else
    {
      cout << "E R R O R ! AuxTools::StopwatchEnd(...) - Invalid unit of time \"" << myUnits << "\" provided. ";
      cout << "Please provide one of the following:\n\"seconds\", \"minutes\", \"hours\". EXIT." << endl;
      exit(1);
    }

  stopwatch_stop = clock();
  double elapsed_time_secs  = double(stopwatch_stop - stopwatch_start) / CLOCKS_PER_SEC;
  double elapsed_time_units = elapsed_time_secs/double(units);
  
  cout << myTitle << ": " << setprecision(myPrecision) << elapsed_time_units << " (" << myUnits << ")" << endl;
  return;

}  


//****************************************************************************  
char* AuxTools::AppendCharToCharArray(char* array,
				      char a)
//****************************************************************************  
{
    size_t len = strlen(array);
    char* ret = new char[len+2];

    strcpy(ret, array);    
    ret[len] = a;
    ret[len+1] = '\0';

    return ret;
}


//****************************************************************************
void AuxTools::ReplaceString(string &myString,
			     string oldSubstring,
			     string newSubstring)
//****************************************************************************
{
  size_t index = 0;
  while (true) {
    
    // Locate the substring to replace
    index = myString.find(oldSubstring, index);
    if (index == std::string::npos) break;

    // Make the replacement
    myString.replace(index, oldSubstring.length(), newSubstring);

    // Advance index forward so the next iteration doesn't pick it up as well
    index += oldSubstring.length();
  }

  return;
}


#endif // AuxTools_cc

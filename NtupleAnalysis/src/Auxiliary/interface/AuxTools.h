#ifndef AuxTools_h
#define AuxTools_h

// System
#include <cmath>
#include <iomanip>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include "constants.h"
#include <ctime>

// User
#include "../src/Table.C"

using namespace std;
using namespace constants;

class AuxTools{

 public:
  AuxTools() { };
  virtual ~AuxTools() { };

  void ProgressBar(Long64_t entry,
		   Long64_t total, 
		   Int_t resolution, 
		   Int_t barWidth=0);
  
  template <typename T> std::string ToString(const T a_value, const int n=3);

  template<class TYPE> float Sgn(TYPE myNumber);

  long long nCr(int n, int r);


  Double_t DeltaEta(const Double_t eta1, 
		    const Double_t eta2);

  Double_t DeltaPhi(const Double_t phi1, 
		    const Double_t phi2);

  Double_t DeltaR(const Double_t eta1, 
		  const Double_t phi1, 
		  const Double_t eta2, 
		  const Double_t phi2);

  TLorentzVector GetTLorentzVector( Double_t pt, 
				    Double_t eta, 
				    Double_t phi, 
				    Double_t e);

  TVector3 GetTVector3(Double_t px, 
		       Double_t py, 
		       Double_t pz);

  TVector2 GetTVector2(Double_t eta, 
		       Double_t phi);
  
  template<class TYPE> void EnsureVectorIsSorted(const vector<TYPE> myVector, 
						 Bool_t bDescendingOrder);

  template<class TYPE> Bool_t VectorIsSorted(const vector<TYPE> myVector, 
					     Bool_t bDescendingOrder);

  template<class TYPE> void PrintVector(const vector<TYPE> myVector, 
					string title="");

  void PrintPSets(TTree *t);
    
  template<class TYPE> string ConvertIntVectorToString(const vector<TYPE> myVector);

  double Divide(int numerator, 
		int denominator);

  void Efficiency(int nPass, 
		  int nTotal, 
		  const string errType, 
		  double &eff, 
		  double &err);

  void StopwatchStart(void);

  void StopwatchStop(const int myPrecision = 5, 
		     const string myUnits = "seconds");

  char* AppendCharToCharArray(char* array,
			      char a);

  void ReplaceString(string &myText,
		     string oldSubstring,
		     string newSubstring);
  
 private:
  clock_t stopwatch_start;
  clock_t stopwatch_stop;
};

#endif

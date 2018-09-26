#ifndef L1PixelTrackFit_h
#define L1PixelTrackFit_h

// System
#include <iostream>
#include <vector>
#include <algorithm>

// User
#include "../../Auxiliary/src/AuxTools.C"
#include "../../DataFormat/src/TTPixelTrack.C"

// ROOT
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"

using namespace std;

class L1PixelTrackFit{     
 public:
  // Constructors/Destructors
  L1PixelTrackFit(double Bz=3.8112); // for CMS Tracker: Bz=3.8112
  ~L1PixelTrackFit();
  
 // Function declaration
 TTPixelTrack FitPixelTrack(double rinv,
 			    double phi0,
			    double d0,
			    double t,
			    double z0,
			    vector<double> candPixHits_X,
			    vector<double> candPixHits_Y,
			    vector<double> candPixHits_Z,
			    vector< int >  candPixHits_Type);

 void Multifit(double rinv,
	       double phi0,
	       double d0,
	       double t,
	       double z0,
	       std::vector<ROOT::Math::XYZVector> hitL1,
	       std::vector<ROOT::Math::XYZVector> hitL2,
	       std::vector<ROOT::Math::XYZVector> hitL3,
	       std::vector<ROOT::Math::XYZVector> hitL4,
	       std::vector<ROOT::Math::XYZVector> hitD1,
	       std::vector<ROOT::Math::XYZVector> hitD2,
	       std::vector<ROOT::Math::XYZVector> hitD3,		
	       bool& success,
	       double& invrfinal,
	       double& phi0final,
	       double& d0final,
	       double& tfinal,
	       double& z0final,
	       double& chisqfinal,
	       int& nhit,		
	       double &sigmarinv,
	       double &sigmaphi0,
	       double &sigmad0,
	       double &sigmat,
	       double &sigmaz0,
	       std::vector<ROOT::Math::XYZVector> &pixHits,
	       unsigned int cfg_NPixelHits_Min=3);
	       

  void TrackFit(double rinv,
		double phi0,
		double d0,
		double t,
		double z0,
		double &rinvfit,
		double &phi0fit,
		double &d0fit,
		double &tfit,
		double &z0fit,
		double &chisqfit,
		double &sigmarinv,
		double &sigmaphi0,
		double &sigmad0,
		double &sigmat,
		double &sigmaz0,
		std::vector<ROOT::Math::XYZVector> fithits,
		std::vector<bool> fitbarrel);
  
  void CalculateDerivatives(double rinv,
			    double phi0,
			    double t,
			    double z0,
			    std::vector<ROOT::Math::XYZVector> fithits,
			    std::vector<bool> fitbarrel,  
			    double D[5][8],
			    double MinvDt[5][8]);
  
  void LinearTrackFit(double rinv,
		      double phi0,
		      double d0,
		      double t,
		      double z0,
		      double &rinvfit,
		      double &phi0fit,
		      double &d0fit,
		      double &tfit,
		      double &z0fit,
		      double &chisqfit,
		      double &sigmarinv,
		      double &sigmaphi0,
		      double &sigmad0,
		      double &sigmat,
		      double &sigmaz0, 
		      std::vector<ROOT::Math::XYZVector> fithits,
		      std::vector<bool> fitbarrel,  
		      double D[5][8],
		      double MinvDt[5][8]);
  
  void Invert(double M[5][10],
	      unsigned int n);
   
  void PrintProperties(void);

  // Additional below this point
  void SaveHitPermutationCandidate(int i0temp,
				   int i1temp,
				   int i2temp,
				   int i3temp,
				   double chiSq,
				   double chiSqDof,
				   double rinv,
				   double phi, 
				   double d0,
				   double t,
				   double z0,
				   string remarks="");
  
  void PrintHitPermutationCandidate(bool bSort=true);

  void SortHitPermutationCandidates(void);

  void ClearCandInfo(void);
  
 private:
 // Function declaration
 void _InitVars(void);
 
 void _SetCandPixHitContainers(vector<double> candPixHits_X,
			       vector<double> candPixHits_Y,
			       vector<double> candPixHits_Z,
			       vector< int >  candPixHits_Type);

 // Variable declaration
 AuxTools auxTools;
 struct SortAscendingFabs{ bool operator() (double a, double b) const { return fabs(a) > fabs(b); } };
 struct SortDescendingFabs{ bool operator() (double a, double b) const { return fabs(a) < fabs(b); } };
 struct SortAscending{ bool operator() (double a, double b) const { return a > b; } }; 
 struct SortDescending{ bool operator() (double a, double b) const { return a < b; } }; 

 struct SortDescendingTarget{
   const std::vector<double>& target;  
 SortDescendingTarget(const std::vector<double>& target): target(target) {}
   bool operator()(int a, int b) const { return target[a] < target[b]; }
 }; // ignore complier warning
 
 struct SortClosestToOneTarget{
   const std::vector<double>& target;  
 SortClosestToOneTarget(const std::vector<double>& target): target(target) {}
   bool operator()(int a, int b) const { return fabs(target[a]-1) < fabs(target[b]-1); }
 }; // ignore complier warning
 

 struct SortAscendingTarget{
   const std::vector<double>& target;  
 SortAscendingTarget(const std::vector<double>& target): target(target) {}
   bool operator()(int a, int b) const { return target[a] > target[b]; }
 }; // ignore complier warning

 double mMagneticFieldStrength; // Bz at (0,0,0). Assume uniform B-field inslide solenoid.
 std::vector<ROOT::Math::XYZVector> pixHitL1;
 std::vector<ROOT::Math::XYZVector> pixHitL2;
 std::vector<ROOT::Math::XYZVector> pixHitL3;
 std::vector<ROOT::Math::XYZVector> pixHitL4;
 std::vector<ROOT::Math::XYZVector> pixHitD1;
 std::vector<ROOT::Math::XYZVector> pixHitD2;
 std::vector<ROOT::Math::XYZVector> pixHitD3;
 double rinvfit;
 double phi0fit;
 double d0fit;
 double tfit;
 double z0fit;
 double chisqfit;
 double sigmarinv;
 double sigmaphi0;
 double sigmad0;
 double sigmat;
 double sigmaz0;
 std::vector<ROOT::Math::XYZVector> pixHits;
 std::vector<ROOT::Math::XYZVector> candidatePixHits;
 // Additional below this point
 std::vector<int> cand_i0;
 std::vector<int> cand_i1;
 std::vector<int> cand_i2;
 std::vector<int> cand_i3;
 std::vector<int> cand_nHits;
 std::vector<double> cand_chiSq;
 std::vector<double> cand_chiSqDof;
 std::vector<double> cand_rinv;
 std::vector<double> cand_phi; 
 std::vector<double> cand_d0;
 std::vector<double> cand_t;
 std::vector<double> cand_z0;
 std::vector<string> cand_remarks;
 bool cand_sorted;
 
};

#endif

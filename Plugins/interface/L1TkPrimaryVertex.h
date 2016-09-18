#ifndef L1TkPrimaryVertex_h
#define L1TkPrimaryVertex_h

// User
#include "../../Auxiliary/src/MCTools.C"
#include "../../Auxiliary/src/AuxTools.C"
#include "../../Auxiliary/src/L1Tracks.C"
#include "../../Auxiliary/src/Table.C"
#include "../../Auxiliary/src/Datasets.C"

using namespace std;

class L1TkPrimaryVertex
{
 public:
  // Constructors/Destructors
  L1TkPrimaryVertex(L1Tracks *s_) { s = s_; };
  ~L1TkPrimaryVertex() {};
  
  // Member Functions
  double GetPrimaryVertexZ(string algo, 
			   vector<int> &PVTks_index);  
    
 private:
  // Member Functions    
  double SlowVtxProducer(string algo, 
			 vector<int> &PVTks_index);
  
  double VtxProducerSumPtBin(const double zCoordinate, 
			     const double delta_z, 
			     const string algo,
			     vector<int> &pvTks_index);
  
  double FastVtxProducer(string algo,
			 vector<int> &PVTks_index);
  
  double PeakFinder(TH1D *htemp,
		    Double_t &center);

  // Variables
  L1Tracks* s;

};

#endif

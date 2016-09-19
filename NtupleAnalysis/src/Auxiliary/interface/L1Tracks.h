#ifndef L1Tracks_h
#define L1Tracks_h

// User
#include "../../Framework/interface/TreeDefinitionReco.h"
#include "../../DataFormat/src/TTPixelTrack.C"
#include "../../Plugins/src/L1PixelTrackFit.C"

#ifdef USING_MC
#include "../../Framework/interface/TreeDefinitionGenP.h" 
#define TREEANALYSER TreeAnalyserMC
#include "../../Framework/interface/TreeAnalyserMC.h"
#else
#include "../../Framework/interface/TreeDefinitionReco.h"
#define TREEANALYSER TreeAnalyserReco
#include "../../Framework/interface/TreeAnalyserReco.h"
#endif // USING_MC


using namespace std;

class L1Tracks
{
 public:
  // Constructors/Destructors
  L1Tracks(TREEANALYSER* t_) { t = t_; };

  // Member Functions
  bool SelectTracks(int iTrack, 
		    const string selectionType,
		    const int nFitParams=5);

  bool SelectPixTracks(int iTrack, 
		       const string selectionType);

  bool SelectPixTracksRefit(TTPixelTrack t,
			    int iTrack,
			    const string selectionType);

  void GetPixTrackAllHits(vector<TVector3> &pixHits_all,
			  vector<int> &pixHits_PixTkRefIndex);
  
  void GetPixTrackSharedHits(int iTrack,
			     vector<int>    &pixTk_PixHits_PixTkRefIndex,
			     vector<double> &pixTk_PixHits_X,
			     vector<double> &pixTk_PixHits_Y,
			     vector<double> &pixTk_PixHits_Z,
			     vector<double> &pixTk_PixHits_R,
			     vector<double> &pixTk_PixHits_Phi,
			     vector<int> &pixTk_PixHits_Type,
			     vector<TVector3> pixHits_XYZ,
			     vector<int> pixHits_TTPixelTrackIndex);
 
  bool HasMustPixelHitInLayerAndDisk(int iTrack,
				     vector<int> pixTk_Barrel_Type,
				     vector<int> pixTk_Endcap_Type,
				     double endcapEtaBoundary=1.5);

  bool HasMustPixelHitPattern(int pixTk_PixHits_Pattern,
			      vector<int> pixTk_PixHits_AllowedPatterns);
      
  bool HasMustPixelHitInLayeOrDisk(int iTrack,
				   vector<int> pixTk_HitTypes);
  
  double GetTrackCharge(const int iTrack);

  double GetPixTrackCharge(const int iTrack);

  double GetTrackPt(const int iTrack);

  double GetPixTrackPt(const int iTrack);
  
  double GetTrackD0(const int iTrack);

  double GetTrackD0Mag(const int iTrack);

  double GetTrackD0Sign(const int iTrack);

  double GetTrackD0Phi(const int iTrack);
  
  int GetPixelIndexOfTrack(const int tk_Index);

  double GetPixTrackD0(const int iTrack);

  double GetPixTrackD0Sig(const int iTrack);

  double GetTrackPOCAzSig(const int iTrack);

  double GetPixTrackPOCAzSig(const int iTrack);

  double GetTrackSigmaPt(const int iTrack);

  double GetPixTrackSigmaPt(const int iTrack);

  double GetTrackRedChiSq(const int iTrack, 
			  const int nFitParams=5);

  double GetPixTrackRedChiSq(const int iTrack);
  
  double GetTrackStubPtConsistency(const int iTrack);
      
  void PrintTrackProperties(const int iTrack,
			    const int nFitParams=5);

  void PrintPixTrackProperties(const int iPixTrack,
			       bool bPrintRefTTTrack=true);

  void PrintPixTrackRefitProperties(TTPixelTrack tk,
				    int pixTk_index,
				    bool bPrintRefTTPixelTrack=true);
    
  unsigned int GetNumOfStubs(const int iTrack);

  unsigned int GetNumOfPSStubs(const int iTrack);

  unsigned int GetNumOfBarrelStubs(const int iTrack);
    
  unsigned int GetNumOfEndcapStubs(const int iTrack);

  TLorentzVector GetP4FromTracks(const vector<int> tks_Index);

  TLorentzVector GetP4FromPixTracks(const vector<int> tks_Index);

  TLorentzVector GetP4FromPixTracksRefit(const vector<int> tks_Index);

  int GetChargeFromTracks(const vector<int> tks_Index);

  int GetChargeFromPixTracks(const vector<int> tks_Index);

  int GetDOF(int iTrack,
	     int nFitParams=5);

  double GetRedChiSquared(int iTrack,
			  int nFitParams);

  // Variables
  TREEANALYSER* t;

 private:
  // Variables
  MCTools mcTools;
  AuxTools auxTools;

};

#endif //ObjectSelect_h

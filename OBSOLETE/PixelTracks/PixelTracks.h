#ifndef PixelTracks_h
#define PixelTracks_h

#include <iostream>

#include "../utilities/TreeAnalyserMC.C"
#include "../utilities/MCTools.C"
#include "../utilities/AuxTools.C"
#include "../utilities/HistoTools.C"
#include "../utilities/L1Tracks.C"
#include "../utilities/TrackingParticles.C"
#include "../utilities/L1PixelTrackFit.C"
#include "../utilities/TTPixelTrack.C"

class PixelTracks : public TreeAnalyserMC{
 public:
 PixelTracks(const string SamplePath,
	       const std::string SampleName,
	       const std::string text_, 
	       const int maxEvents_ = -1, 
	       TTree* tree=0) : 
  TreeAnalyserMC("PixelTracks", SamplePath, SampleName, text_, maxEvents_, tree) {};
  virtual void Loop();
  
 private:

  // Objects/Variables
  AuxTools auxTools;

};

#endif // PixelTracks_h

#ifndef TreeAnalyserReco_cxx
#define TreeAnalyserReco_cxx

// System
#include <typeinfo>

#define TREEDEFINITIONRECO TreeDefinitionReco
// User
#include "../interface/TreeAnalyserReco.h"
#include "../../Auxiliary/src/L1Tracks.C"
#include "../../Plugins/src/L1PixelTrackFit.C"

void TreeAnalyserReco::InitSelector()
{
  s = new L1Tracks(this);
  f  = new L1PixelTrackFit(this);
}

#endif // TreeAnalyserReco_cxx

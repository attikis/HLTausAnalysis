#ifndef Datasets_h
#define Datasets_h

// System
#include <algorithm>

// User
#include "../interface/Table.h"
#include "../interface/AuxTools.h"

using namespace std;

class Datasets{
 public:
  // Default and overloaded constructors & destructors
  Datasets() {}; 
  Datasets(string alias,
	   string datasetPath, 
	   string datasetPathExt, 
	   string productionType,
	   string cmssw,
	   string geometry,
	   int nPileUp, 
	   int nEvents,
	   int McTauMomPdgId,
	   int nMcTaus);

  ~Datasets() {};

  // Member Functions
  void PrintSamples(const string mcProduction, 
		    const string textFormat = "Text", 
		    bool bCompact=false);

  const string GetDatasetPathFromAlias(const string datasetAlias);
  
  Datasets GetDataset(const string datasetAlias);

  // Variables
  AuxTools auxTools_;

  vector<Datasets> datasets_TTI2023Updg14D;
  vector<Datasets> datasets_PhaseIISpring17D;
  vector<Datasets> datasets_PhaseIIFall17D;

  string alias_;
  string process_;
  string productionType_;
  string cmssw_;
  string geometry_;
  string datasetPath_;
  string datasetPathExt_;
  int nPileUp_;
  int nEvents_;
  int nMcTaus_;
  int McTauMomPdgId_;
  vector<string> mcProductions;

  private:
  // Member Functions
  void CreateMcProductions_(void);

  void PrintDatasetsVector_(vector<Datasets> myDatasets, 
			    const string textFormat = "Text", 
			    bool bCompact=false);

  void IsValidMcProduction_(const string mcProduction,
			    vector<Datasets> &myDatasets);

  bool IsValidDatasetAlias(const string datasetName);

  bool IsValidDatasetName(const string datasetName);

};

#endif

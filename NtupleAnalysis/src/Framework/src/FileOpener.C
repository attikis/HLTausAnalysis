#ifndef FileOpener_cxx 
#define FileOpener_cxx

// System
#include <fstream>

// User 
#include "../interface/FileOpener.h"

// ROOT
#include "TChain.h"

TChain* FileOpener::OpenFile(const std::string SamplePath, const std::string SampleName, TChain* mychain)
{

  // Get the dataset attributes
  const string datasetPath = datasets_.GetDatasetPathFromAlias(SampleName);
  
  // Ensure file exists
  std::string FullFileName     = SamplePath + "/" + datasetPath + "/res/output-" + datasetPath + ".root"; 
  std::string FullFileName_alt = SamplePath + "/" + datasetPath + "/results/output-" + datasetPath + ".root";
  std::string name = FullFileName;
  std::ifstream GetFile(FullFileName.c_str());
  if (!GetFile.good() )
    {
      cout << "=== FileOpener::OpenFile()\n\tFile \"" << FullFileName << "\" does not exist. Trying alternative paths ..." << endl;
      name = FullFileName_alt;

      std::ifstream GetFile_alt(FullFileName_alt.c_str());
      if (!GetFile_alt.good() )
	{
	  std::ifstream GetFile_alt(FullFileName_alt.c_str());
	cout << "=== FileOpener::OpenFile()\n\tFile \"" << FullFileName_alt << "\" also does not exist. EXIT" << endl;
	exit(1);
	}
    }
  
  /// Inform user of ROOT file added to the TChain
  cout << "=== FileOpener::OpenFile()\n\tAdding file " << name << " to the chain" << endl;
  mychain -> Add(name.c_str());

  return mychain; 
}

#endif // FileOpener_cxx

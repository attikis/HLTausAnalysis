#ifndef FileOpener_cxx 
#define FileOpener_cxx

// System
#include <fstream>
#include <vector>

// User 
#include "../interface/FileOpener.h"

// ROOT
#include "TChain.h"
#include "TSystem.h"

TChain* FileOpener::OpenFiles(const string multicrabPath, const string dataset, TChain* mychain)
{
  cout << "=== FileOpener::OpenFiles()" << endl;

  // Get the dataset attributes
  string rootFileName;
  const string datasetPath = datasets_.GetDatasetPathFromAlias(dataset);
  const string fullPath    = multicrabPath + "/" + datasetPath + "/results/";
  const string filePrefix  = "raw2TTree_";
  vector<string> dirs      = GetListOfFiles(fullPath, filePrefix);
  
  // Sanity check
  cout << "\tFound " << dirs.size() << " files of type " << filePrefix << "*.root under " << fullPath << endl;
  if (dirs.size() < 1) exit(1);
  else
    {
      // cout << "\tAdding " << dirs.size() << " files to the TChain" << endl;
    }
  
  // For-loop: All ROOT files
  for (vector<string>::iterator f = dirs.begin(); f != dirs.end(); f++)
    {
      rootFileName = *f;
      ifstream file(rootFileName);
      
      if (!file.good() )
	{
	  cout << "\tFile \"" << rootFileName << "\" does not exist" << endl;
   	  exit(1);
   	}
      else
	{
	  // cout << "\tAdding file " << rootFileName << " to the chain" << endl;
	  mychain -> Add(rootFileName.c_str());
	}
    }
  
  return mychain; 
}


vector<string> FileOpener::GetListOfFiles(const string fullPath, const string filePrefix)
{

  // std::cout << "=== FileOpener::GetListOfFiles()" << std::endl;
  
  // Declare variables
  const char* entry;
  std::vector<string> dirs;
  const char* ext = ".root";
  TString str;  
  char* dir  = gSystem->ExpandPathName(fullPath.c_str()); // cast a "std::string" to a "char*" with c_str()
  void* dirp = gSystem->OpenDirectory(dir);
  vector<string> fileNames;


  // Get all files in the dir fullPath
  while((entry = (char*)gSystem->GetDirEntry(dirp)))
    {
      str = entry;
      if(str.EndsWith(ext)) fileNames.push_back((string)str);
    }
 
  // For-loop: All ROOT files 
  for (auto f = fileNames.begin(); f != fileNames.end(); f++)
    {
      // cout << "\tProcessing file " << *f << endl;

      // Skip all not-ROOT files 
      if (!EndsWith(*f, ".root")) continue;

      // Look for files containing the prefix
      size_t found = f->find(filePrefix);

      // If prefix found in string add it for use as input ROOT file
      if (found!=std::string::npos) dirs.push_back(fullPath + "/" + *f);
      else std::cout << "\tSkipping " << *f << endl;
    }

  // cout << "\tFound " << dirs.size() << " files of type " << filePrefix << "*.root under directory " << fullPath << endl;
  return dirs;
}

bool FileOpener::EndsWith(const std::string &str, const std::string &suffix)
{
  // cout << "=== FileOpener::EndsWith()" << endl;
  bool sizeOk   = str.size() >= suffix.size();
  if (!sizeOk) return false;

  // cout << "\tstr = " << str << ", suffix  = " << suffix << endl;
  bool suffixOk = str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
  return suffixOk;
}

#endif // FileOpener_cxx

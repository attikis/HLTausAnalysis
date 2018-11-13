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

  if (0) cout << "=== FileOpener::OpenFile()" << endl;

  // Get the dataset attributes
  string rootFileName;
  const string datasetPath = datasets_.GetDatasetPathFromAlias(dataset);
  const string fullPath    = multicrabPath + "/" + datasetPath + "/results/";
  vector<string> dirs      = GetListOfFiles(fullPath);

  // Sanity check
  if (0) {
    if (dirs.size() < 1)
      {
	cout << "\tFound " << dirs.size() << " ROOT files under directory " << fullPath << "! EXIT" << endl;
	exit(1);      
      }
    else if (0) cout << "\tFound " << dirs.size() << " ROOT files under directory " << fullPath << endl;
  }

  // For-loop: All ROOT files
  for (vector<string>::iterator f = dirs.begin(); f != dirs.end(); f++)
    {     
      rootFileName = *f;
      if (0) cout << "\tAdding file " << rootFileName << " to the chain" << endl;
      
      ifstream file(rootFileName);
      
      if (!file.good() )
	{
   	  if(0) cout << "\tFile \"" << rootFileName << "\" also does not exist. EXIT" << endl;
   	  exit(1);
   	}
      else mychain -> Add(rootFileName.c_str());
    }
  
  return mychain; 
}


string FileOpener::GetFirstFileName(const string multicrabPath, const string dataset)
{

  if (0) cout << "=== FileOpener::GetFirstFileName()" << endl;

  // Get the dataset attributes
  string rootFileName;
  const string datasetPath = datasets_.GetDatasetPathFromAlias(dataset);
  const string fullPath    = multicrabPath + "/" + datasetPath + "/results/";
  vector<string> dirs      = GetListOfFiles(fullPath);

  // Sanity check
  if (dirs.size() < 1)
    {
      cout << "\tFound " << dirs.size() << " ROOT files under directory " << fullPath << "! EXIT" << endl;
      exit(1);      
    }
  else if(0) cout << "\tFound " << dirs.size() << " ROOT files under directory " << fullPath << endl;

  vector<string>::iterator f = dirs.begin();
  rootFileName = *f;

  ifstream file(rootFileName);
  
  if (!file.good() )
    {
      if(0) cout << "\tFile \"" << rootFileName << "\" also does not exist. EXIT" << endl;
      exit(1);
    }
  
  return rootFileName; 
}


vector<string> FileOpener::GetListOfFiles(const string fullPath)
{
  
  // Declare variables
  const char* entry;
  vector<string> dirs;
  vector<string> dirs_new;
  
  Int_t n = 0;
  const char* ext = ".root";
  
  
  char* dir  = gSystem->ExpandPathName(fullPath.c_str()); // cast a "std::string" to a "char*" with c_str()
  void* dirp = gSystem->OpenDirectory(dir);
  const char* filename[1000];
  TString str;
  
  // Get all files in the dir fullPath
  while((entry = (char*)gSystem->GetDirEntry(dirp)))
    {
      str = entry;
      if(str.EndsWith(ext)) filename[n++] = gSystem->ConcatFileName(dir, entry);
    }
  
  // For-loop: All files
  for (Int_t i = 0; i < n; i++)
    {
      // Printf("file -> %s", filename[i]);
      string fileName = string(filename[i]);
      size_t found    = fileName.find("L1Ntuple");
      if (found!=std::string::npos) dirs.push_back(fileName);
      else {}//cout << "\tSkipping " << fileName << endl;
      
    }


  return dirs;
}


#endif // FileOpener_cxx

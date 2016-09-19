#ifndef Datasets_cxx
#define Datasets_cxx

// User
#include "../interface/Datasets.h"

//****************************************************************************
Datasets::Datasets(string alias,
		   string datasetPath, 
		   string datasetPathExt, 
		   string productionType,
		   string cmssw,
		   string geometry,
		   int nPileUp, 
		   int nEvents,
		   int McTauMomPdgId,
		   int nMcTaus)
//****************************************************************************
{
  
  alias_          = alias;
  productionType_ = productionType;
  cmssw_          = cmssw;
  geometry_       = geometry;
  datasetPath_    = datasetPath;
  datasetPathExt_ = datasetPathExt;
  nPileUp_        = nPileUp;
  nEvents_        = nEvents;
  nMcTaus_        = nMcTaus;
  McTauMomPdgId_  = McTauMomPdgId;

}


//****************************************************************************
void Datasets::PrintSamples(const string mcProduction, 
			    const string textFormat, 
			    bool bCompact)
//****************************************************************************
{
  
  // Declarations
  vector<Datasets> myDatasets;

  // First create the tables for the known MC Productions
  if (datasets_TTI2023Updg14D.size() < 1 ) CreateMcProductions_();

  // Check that the mcProduction variable is valid
  IsValidMcProduction_(mcProduction, myDatasets);

  // Print full or compact infromation table
  PrintDatasetsVector_(myDatasets, textFormat, bCompact);

  return;
}


//****************************************************************************
void Datasets::IsValidMcProduction_(const string mcProduction, 
				    vector<Datasets> &myDatasets)
//****************************************************************************
{

  // Declare currently supported MC-Productions
  mcProductions.push_back("2023TTIUpg14");

  // Check if the selected MC Prodution is valid
  for( int i = 0 ; i < (int) mcProductions.size(); i++){

    if ( mcProduction.compare(mcProductions.at(i)) == 0) {
      myDatasets = datasets_TTI2023Updg14D;
      return;
    }

  }

  // If this is reached then the MC-Production is invalid  
  cout << "\nE R R O R ! Datasets::PrintSamples(...) - Unknown MC Production \"" << mcProduction << "\"." << endl;
  cout << "EXIT" << endl;
  exit(1);
      
  return;
}


//****************************************************************************
void Datasets::IsValidDatasetName(const string datasetName)
//****************************************************************************
{


  // First create the tables for the known MC Productions
  if (datasets_TTI2023Updg14D.size() < 1 ) CreateMcProductions_();

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){
    if ( datasetName.compare( datasets_TTI2023Updg14D.at(iD).datasetPath_) == 0) return;
  }
  
  // If this is reached then the dataset sample is invalid
  cout << "\nE R R O R ! Datasets::IsValidDataset(...) - Unknown dataset name \"" << datasetName << "\". See below all available dataset names." << endl;
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){
    cout << "\"" << datasets_TTI2023Updg14D.at(iD).datasetPath_ << "\"" << endl;
  }
  cout << "EXIT" << endl;
  exit(1);
  
  return;
}


//****************************************************************************
void Datasets::IsValidDatasetAlias(const string datasetName)
//****************************************************************************
{


  // First create the tables for the known MC Productions
  if (datasets_TTI2023Updg14D.size() < 1 ) CreateMcProductions_();

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){
    if ( datasetName.compare( datasets_TTI2023Updg14D.at(iD).alias_) == 0) return;
  }
  
  // If this is reached then the dataset sample is invalid
  cout << "\nE R R O R ! Datasets::IsValidDatasetAlias(...) - Unknown dataset alias \"" << datasetName << "\". See below all available dataset aliases. EXIT" << endl;
  PrintDatasetsVector_(datasets_TTI2023Updg14D, "Text", true);
  cout << "EXIT" << endl;
  exit(1);
  
  return;
}


//****************************************************************************
void Datasets::PrintDatasetsVector_(const vector<Datasets> myDatasets, 
				    const string textFormat, 
				    bool bCompact)
//****************************************************************************
{

  if (myDatasets.size() < 1) return;
  if (bCompact)
    {
      Table mcSamples("Alias | Dataset | Production | CMSSW | PU | Events", textFormat, "l l l c c c");
      for (int i = 0; i < (int) myDatasets.size(); i++){

	mcSamples.AddRowColumn(i, myDatasets.at(i).alias_);
	mcSamples.AddRowColumn(i, myDatasets.at(i).datasetPath_);
	mcSamples.AddRowColumn(i, myDatasets.at(i).productionType_);
	mcSamples.AddRowColumn(i, myDatasets.at(i).cmssw_);
	mcSamples.AddRowColumn(i, auxTools_.ToString( myDatasets.at(i).nPileUp_) );
	mcSamples.AddRowColumn(i, auxTools_.ToString( myDatasets.at(i).nEvents_) );

      }
      mcSamples.Print();
    }
  else
    {
      Table mcSamples("Alias | Dataset | Dataset (CMS DAS) | Production | CMSSW | Geometry | PU | Events", textFormat, "l l l l c c c c");
      for (int i = 0; i < (int) myDatasets.size(); i++){
	
	mcSamples.AddRowColumn(i, myDatasets.at(i).alias_);
	mcSamples.AddRowColumn(i, myDatasets.at(i).datasetPath_);
	mcSamples.AddRowColumn(i, myDatasets.at(i).datasetPathExt_); 	// mcSamples.AddRowColumn(i, "X");
	mcSamples.AddRowColumn(i, myDatasets.at(i).productionType_);
	mcSamples.AddRowColumn(i, myDatasets.at(i).cmssw_);
	mcSamples.AddRowColumn(i, myDatasets.at(i).geometry_);
	mcSamples.AddRowColumn(i, auxTools_.ToString( myDatasets.at(i).nPileUp_) );
	mcSamples.AddRowColumn(i, auxTools_.ToString( myDatasets.at(i).nEvents_) );
	
      }
      mcSamples.Print();
    }
      
  return;
}


//****************************************************************************
void Datasets::CreateMcProductions_(void)
//****************************************************************************
{

  // TTI2023Upgd14D
  string path     = "TTI2023Upg14D-PU140bx25";
  string pathExt1 = "/" + path + "_PH2_1K_FB_V3-v2/GEN-SIM-DIGI-RAW";
  string pathExt2 = "/" + path + "_PH2_1K_FB_V3-v1/GEN-SIM-DIGI-RAW";
  string cmssw    = "6_2_0_SLHC12";
  string geometry = "Extended2023TTI";
  string PP       = "Private";
  string CP       = "Central";

  // Private Production
  Datasets SingleMuMinus_TTI2023_Pt_2_10_NoPU ( "SingleMuMinus_Pt_2_10_NoPU", "SingleMuMinus_E2023TTI_Pt_2_10_NoPU", "SingleMuMinus_TTI2023_Pt_2_10_NoPU", PP, cmssw, geometry,   0, 20000, 0, 0);
  Datasets SingleMuPlus_TTI2023_Pt_2_10_NoPU  ( "SingleMuPlus_Pt_2_10_NoPU" , "SingleMuPlus_E2023TTI_Pt_2_10_NoPU" , "SingleMuPlus_TTI2023_Pt_2_10_NoPU" , PP, cmssw, geometry,   0, 20000, 0, 0);
  Datasets SingleMuon_E2023TTI_NoPU           ( "SingleMuon_NoPU"           , "SingleMuon_E2023TTI_NoPU"           , "SingleMuon_E2023TTI_NoPU"          , PP, cmssw, geometry,   0, 20000, 0, 0);
  Datasets SingleMuPlus_E2023TTI_NoPU         ( "SingleMuPlus_NoPU"         , "SingleMuPlus_E2023TTI_NoPU"         , "SingleMuPlus_E2023TTI_NoPU"        , PP, cmssw, geometry,   0, 20000, 0, 0);
  Datasets SingleMuon_E2023TTI_PU140          ( "SingleMuon_PU140"          , "SingleMuon_E2023TTI_PU140"          , "SingleMuon_E2023TTI_PU140"         , PP, cmssw, geometry, 140, 20000, 0, 0);
  // Central Production
  Datasets MinBias_TTI2023       ( "MinBias"       , "Neutrino_Pt2to20_gun_"                   + path, "Neutrino_Pt2to20_gun_"                   + pathExt1, CP, cmssw, geometry, 140, 284300,  0, 0);
  Datasets PiPlus_TTI2023        ( "PiPlus"        , "SinglePionPlusFlatPt0p2To50_"            + path, "SinglePionPlusFlatPt0p2To50_"            + pathExt1, CP, cmssw, geometry, 140,  50000,  0, 0);
  Datasets PiMinus_TTI2023       ( "PiMinus"       , "SinglePionMinusFlatPt0p2To50_"           + path, "SinglePionMinusFlatPt0p2To50_"           + pathExt1, CP, cmssw, geometry, 140,  50000,  0, 0);
  Datasets SingleTauGun1p_TTI2023( "SingleTauGun1p", "SingleTauOneProngFlatPt10To100_"         + path, "SingleTauOneProngFlatPt10To100_"         + pathExt1, CP, cmssw, geometry, 140,  50000,  0, 1);
  Datasets DiTauGun3p_TTI2023    ( "DiTauGun3p"    , "TauThreeProngsEnriched_"                 + path, "TauThreeProngsEnriched_"                 + pathExt1, CP, cmssw, geometry, 140,  50000,  0, 2);
  Datasets VBF_HToTauTau_TTI2023 ( "VBF"           , "VBF_HToTauTau_125_14TeV_powheg_pythia6_" + path, "VBF_HToTauTau_125_14TeV_powheg_pythia6_" + pathExt2, CP, cmssw, geometry, 140,  24977, 25, 2);
  Datasets TTbar_TTI2023         ( "TTBar"         , "PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV_"  + path, "PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV_"  + pathExt1, CP, cmssw, geometry, 140,  99535, 24, 2);
  Datasets HPlus200_TTI2023      ( "HPlus200"      , "PYTHIA_Tauola_TB_ChargedHiggs200_14TeV_" + path, "PYTHIA_Tauola_TB_ChargedHiggs200_14TeV_" + pathExt1, CP, cmssw, geometry, 140,  50000, 37, 1);
  Datasets HPlus160_TTI2023( "HPlus160" , "PYTHIA6_Tauola_TTbar_ChargedHiggs160_taunu_14TeV_"  + path, "PYTHIA6_Tauola_TTbar_ChargedHiggs160_taunu_14TeV_"+pathExt1, CP, cmssw, geometry, 140, 49680, 37, 1);
  Datasets SingleElectron_TTI2023( "SingleElectron", "SingleElectronFlatPt0p2To50_" + path, "SingleElectronFlatPt0p2To50_" + pathExt1, CP, cmssw, geometry, 140,  50000, 0, 0);
  Datasets SinglePositron_TTI2023( "SinglePositron", "SinglePositronFlatPt0p2To50_" + path, "SinglePositronFlatPt0p2To50_" + pathExt1, CP, cmssw, geometry, 140,  50000, 0, 0);
  Datasets SingleMuPlus_TTI2023  ( "SingleMuPlus"  , "SingleMuPlusFlatPt0p2To150_"  + path, "SingleMuPlusFlatPt0p2To150_"  + pathExt1, CP, cmssw, geometry, 140, 100000, 0, 0);
  Datasets SingleMuMinus_TTI2023 ( "SingleMuMinus" , "SingleMuMinusFlatPt0p2To150_" + path, "SingleMuMinusFlatPt0p2To150_" + pathExt1, CP, cmssw, geometry, 140, 168000, 0, 0);
  Datasets SinglePhoton_TTI2023  ( "SinglePhoton"  , "SinglePhotonFlatPt5To75_"     + path, "SinglePhotonFlatPt5To75_"     + pathExt1, CP, cmssw, geometry, 140, 100000, 0, 0);  
  // 
  datasets_TTI2023Updg14D.push_back(SingleMuMinus_TTI2023_Pt_2_10_NoPU);
  datasets_TTI2023Updg14D.push_back(SingleMuPlus_TTI2023_Pt_2_10_NoPU);
  datasets_TTI2023Updg14D.push_back(SingleMuon_E2023TTI_NoPU);
  datasets_TTI2023Updg14D.push_back(SingleMuPlus_E2023TTI_NoPU);
  datasets_TTI2023Updg14D.push_back(SingleMuon_E2023TTI_PU140);
  datasets_TTI2023Updg14D.push_back(MinBias_TTI2023);
  datasets_TTI2023Updg14D.push_back(PiPlus_TTI2023);
  datasets_TTI2023Updg14D.push_back(PiMinus_TTI2023);
  datasets_TTI2023Updg14D.push_back(SingleTauGun1p_TTI2023);
  datasets_TTI2023Updg14D.push_back(DiTauGun3p_TTI2023);
  datasets_TTI2023Updg14D.push_back(VBF_HToTauTau_TTI2023);
  datasets_TTI2023Updg14D.push_back(TTbar_TTI2023);
  datasets_TTI2023Updg14D.push_back(HPlus200_TTI2023);
  datasets_TTI2023Updg14D.push_back(HPlus160_TTI2023);
  datasets_TTI2023Updg14D.push_back(SingleElectron_TTI2023);
  datasets_TTI2023Updg14D.push_back(SinglePositron_TTI2023);
  datasets_TTI2023Updg14D.push_back(SingleMuPlus_TTI2023);
  datasets_TTI2023Updg14D.push_back(SingleMuMinus_TTI2023);
  datasets_TTI2023Updg14D.push_back(SinglePhoton_TTI2023);

  return;
}


//****************************************************************************
Datasets Datasets::GetDataset(const string datasetAlias)
//****************************************************************************
{

  Datasets d;
  bool bSuccess = false;

  // First create the tables for the known MC Productions
  if (datasets_TTI2023Updg14D.size() < 1 ) CreateMcProductions_();

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){

    if ( datasetAlias.compare( datasets_TTI2023Updg14D.at(iD).alias_) != 0) continue;
    bSuccess = true;
    d = datasets_TTI2023Updg14D.at(iD);
  }
  
  if(!bSuccess){
    cout << "\nE R R O R ! Datasets::GetDataset(...) - Unexpected error! Could not find dataset alias \"" << datasetAlias << "\"." << endl;
    exit(1);
  }

  return d;
}


//****************************************************************************
const string Datasets::GetDatasetPathFromAlias(const string datasetAlias)
//****************************************************************************
{

  string datasetPath = "none";

  // First create the tables for the known MC Productions
  if (datasets_TTI2023Updg14D.size() < 1 ) CreateMcProductions_();

  // Ensure that the supplied sample name is valids
  IsValidDatasetAlias(datasetAlias);  

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){

    if ( datasetAlias.compare( datasets_TTI2023Updg14D.at(iD).alias_) == 0){
      datasetPath = datasets_TTI2023Updg14D.at(iD).datasetPath_;
      return datasetPath;
    }

  }
  
  // If this is reached then the dataset sample is invalid
  cout << "\nE R R O R ! Datasets::GetDatasetPathFromAlias(...) - Unexpected error! Could not find dataset path for dataset alias \"" << datasetAlias << "\"." << endl;
  cout << "EXIT" << endl;
  exit(1);

  return datasetPath;
}

#endif

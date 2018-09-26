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
  if (datasetPathExt == "") datasetPathExt_ = datasetPath;
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
  // mcProductions.push_back("2023TTIUpg14");     // obsolete
  // mcProductions.push_back("PhaseIISpring17D"); // obsolete 
  mcProductions.push_back("PhaseIIFall17D");

  // Check if the selected MC Prodution is valid
  for( int i = 0 ; i < (int) mcProductions.size(); i++){

    if ( mcProduction.compare(mcProductions.at(i)) == 0) {
      myDatasets = datasets_TTI2023Updg14D;
      return;
    }

  }

  // If this is reached then the MC-Production is invalid  
  cout << "=== Datasets::PrintSamples() - Unknown MC Production \"" << mcProduction << "\"." << endl;
  cout << "EXIT" << endl;
  exit(1);
      
  return;
}


//****************************************************************************
bool Datasets::IsValidDatasetName(const string datasetName)
//****************************************************************************
{

  // First create the tables for the known MC Productions
  if (datasets_TTI2023Updg14D.size() < 1 ) CreateMcProductions_();

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){
    if ( datasetName.compare( datasets_TTI2023Updg14D.at(iD).datasetPath_) == 0) return true;
  }
  
  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_PhaseIISpring17D.size(); iD++)
    {
      if ( datasetName.compare( datasets_PhaseIISpring17D.at(iD).datasetPath_) == 0) return true;
    }

  // Check if the selected MC Prodution is valid                                                                                                                        
  for( int iD = 0 ; iD < (int) datasets_PhaseIIFall17D.size(); iD++)
    {
      if ( datasetName.compare( datasets_PhaseIIFall17D.at(iD).datasetPath_) == 0) return true;
    }

  // If this is reached then the dataset sample is invalid
  cout << "=== Datasets::IsValidDataset() - Unknown dataset name \"" << datasetName << "\". See below all available dataset names." << endl;
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){
    cout << "\"" << datasets_TTI2023Updg14D.at(iD).datasetPath_ << "\"" << endl;
  }
  cout << "EXIT" << endl;
  exit(1);
  
  return true;
}


//****************************************************************************
bool Datasets::IsValidDatasetAlias(const string datasetName)
//****************************************************************************
{

  // First create the tables for the known MC Productions
  // if (datasets_TTI2023Updg14D.size() < 1 )   CreateMcProductions_(); // obsolete
  // if (datasets_PhaseIISpring17D.size() < 1 ) CreateMcProductions_(); // obsolete
  if (datasets_PhaseIIFall17D.size() < 1 ) CreateMcProductions_();


  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++)
    {
      if ( datasetName.compare( datasets_TTI2023Updg14D.at(iD).alias_) == 0) return true;
    }

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_PhaseIISpring17D.size(); iD++)
    {
      if ( datasetName.compare( datasets_PhaseIISpring17D.at(iD).alias_) == 0) return true;
    }

  // Check if the selected MC Prodution is valid                                                                                                                        
  for( int iD = 0 ; iD < (int) datasets_PhaseIIFall17D.size(); iD++)
    {
      if ( datasetName.compare( datasets_PhaseIIFall17D.at(iD).alias_) == 0) return true;
    }

  // If this is reached then the dataset sample is invalid
  cout << "=== Datasets::IsValidDatasetAlias() - Unknown dataset alias \"" << datasetName << "\". See below all available dataset aliases. EXIT" << endl;

  // TTI2023Updg14D 
  if (0) PrintDatasetsVector_(datasets_TTI2023Updg14D, "Text", true);
  cout << "" << endl;

  // PhaseIISpring17D
  if (0) PrintDatasetsVector_(datasets_PhaseIISpring17D, "Text", true);

  // PhaseIIFall17D
  if (1) PrintDatasetsVector_(datasets_PhaseIIFall17D, "Text", true);
  cout << "EXIT" << endl;
  exit(1);
  
  return true;
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
  string cmssw    = "9_3_7";
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


  // PhaseIISpring17D
  string path_NoPU  = "/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW";
  string path_PU140 = "/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW";
  string path_PU200 = "/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW";
  string path_NoPU_pilot  = "/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW";
  string path_PU140_pilot = "/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW";
  string path_PU200_pilot = "/PhaseIISpring17D-PU200_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW";
  cmssw    = "9_0_0";
  geometry = "Phase2";
  
  Datasets SingleNeutrinoPU140_PhaseIISpring17D("SingleNeutrino_PU140", "SingleNeutrino_PU140", "/SingleNeutrino" + path_PU140, CP, cmssw, geometry, 140, 500000, 0, 0);
  Datasets SingleNeutrinoPU200_PhaseIISpring17D("SingleNeutrino_PU200", "SingleNeutrino_PU200", "/SingleNeutrino" + path_PU200, CP, cmssw, geometry, 140, 500000, 0, 0);
  
  Datasets SinglePion0NoPU_PhaseIISpring17D("SinglePion0_NoPU"  , "SinglePion0_NoPU" , "/SinglePion0_FlatPt-8to100" + path_NoPU , CP, cmssw, geometry,   0, 500000, 0, 0);
  Datasets SinglePion0PU140_PhaseIISpring17D("SinglePion0_PU140", "SinglePion0_PU140", "/SinglePion0_FlatPt-8to100" + path_PU140, CP, cmssw, geometry, 140, 500000, 0, 0);
  Datasets SinglePion0PU200_PhaseIISpring17D("SinglePion0_PU200", "SinglePion0_PU200", "/SinglePion0_FlatPt-8to100" + path_PU200, CP, cmssw, geometry, 200, 499400, 0, 0);
  
  Datasets SinglePionNoPU_PhaseIISpring17D("SinglePion_NoPU"  , "SinglePion_NoPU" , "/SinglePion_FlatPt-8to100" + path_NoPU , CP, cmssw, geometry,   0, 492588, 0, 0);
  Datasets SinglePionPU140_PhaseIISpring17D("SinglePion_PU140", "SinglePion_PU140", "/SinglePion_FlatPt-8to100" + path_PU140, CP, cmssw, geometry, 140, 499850, 0, 0);
  Datasets SinglePionPU200_PhaseIISpring17D("SinglePion_PU200", "SinglePion_PU200", "/SinglePion_FlatPt-8to100" + path_PU200, CP, cmssw, geometry, 200, 498400, 0, 0);
  
  Datasets SingleTauNoPU_PhaseIISpring17D("SingleTau_NoPU"  , "SingleTau_NoPU" , "/SingleTau_FlatPt-8to150" + path_NoPU , CP, cmssw, geometry, 0, 242578, 0, 1);
  Datasets SingleTauPU140_PhaseIISpring17D("SingleTau_PU140", "SingleTau_PU140", "/SingleTau_FlatPt-8to150" + path_PU140, CP, cmssw, geometry, 0, 244605, 0, 1);
  Datasets SingleTauPU200_PhaseIISpring17D("SingleTau_PU200", "SingleTau_PU200", "/SingleTau_FlatPt-8to150" + path_PU200, CP, cmssw, geometry, 0, 245355, 0, 1);
  
  Datasets TTNoPU_PhaseIISpring17D("TT_TuneCUETP8M1_14TeV_NoPU"  , "TT_TuneCUETP8M1_14TeV_NoPU", "/TT_TuneCUETP8M1_14TeV-powheg-pythia8/" + path_NoPU , CP, cmssw, geometry, 0, 5000, 24, 2);
  Datasets TTPU140_PhaseIISpring17D("TT_TuneCUETP8M1_14TeV_PU140", "TT_TuneCUETP8M1_14TeV_PU140", "/TT_TuneCUETP8M1_14TeV-powheg-pythia8/" + path_PU140, CP, cmssw, geometry, 0, 5000, 24, 2);
  Datasets TTPU200_PhaseIISpring17D("TT_TuneCUETP8M1_14TeV_PU200", "TT_TuneCUETP8M1_14TeV_PU200", "/TT_TuneCUETP8M1_14TeV-powheg-pythia8/" + path_PU200, CP, cmssw, geometry, 0, 5000, 24, 2);
    
  datasets_PhaseIISpring17D.push_back(SingleNeutrinoPU140_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SingleNeutrinoPU200_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SinglePion0NoPU_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SinglePion0PU140_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SinglePion0PU200_PhaseIISpring17D);  
  datasets_PhaseIISpring17D.push_back(SinglePionNoPU_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SinglePionPU140_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SinglePionPU200_PhaseIISpring17D);  
  datasets_PhaseIISpring17D.push_back(SingleTauNoPU_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SingleTauPU140_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(SingleTauPU200_PhaseIISpring17D);  
  datasets_PhaseIISpring17D.push_back(TTNoPU_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(TTPU140_PhaseIISpring17D);
  datasets_PhaseIISpring17D.push_back(TTPU200_PhaseIISpring17D);

  // PhaseIIFall17D                                                                                                                                                     
  string path_NoPU_v1  = "/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW";
  string path_NoPU_v2  = "/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW"; 
  string path_PU140_v1 = "/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW";
  string path_PU140_v2 = "/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW"; 
  string path_PU200_v1 = "/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW";
  string path_PU200_v2 = "/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW"; 
  cmssw    = "10_1_5";
  geometry = "Phase2";

  Datasets SinglePionNoPU_PhaseIIFall17D("SinglePion_FlatPt_2to100_NoPU", "SinglePion_FlatPt_2to100_NoPU", "/SinglePion_FlatPt-2to100" + path_NoPU_v1 , CP, cmssw, geometry, 140, 500000, 0, 0);

  Datasets TTNoPU_PhaseIIFall17D("TT_TuneCUETP8M2T4_14TeV_L1TnoPU", "TT_TuneCUETP8M2T4_14TeV_powheg_pythia8_PhaseIIFall17D_L1TnoPU_93X", "/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8" + path_NoPU_v2, CP, cmssw, geometry, 0, 99328, 24, 2);
  Datasets TTPU140_PhaseIIFall17D("TT_TuneCUETP8M2T4_14TeV_L1TPU140", "TT_TuneCUETP8M2T4_14TeV_powheg_pythia8_PhaseIIFall17D_L1TPU140_93X", "/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8" + path_PU140_v2, CP, cmssw, geometry, 140, 99328, 24, 2);

  Datasets TTPU200_PhaseIIFall17D("TT_TuneCUETP8M2T4_14TeV_L1TPU200", "TT_TuneCUETP8M2T4_14TeV_powheg_pythia8_PhaseIIFall17D_L1TPU200_93X", "/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8" + path_PU200_v2, CP, cmssw, geometry, 200, 99328, 24, 2);

  Datasets SingleTauNoPU_2023D17("SingleTau_L1TnoPU", "RelValSingleTauFlatPt2To100_pythia8_93X_upgrade2023_realistic_v5_2023D17noPU_93X", "/RelValSingleTauFlatPt2To100_pythia8/CMSSW_9_3_7-93X_upgrade2023_realistic_v5_2023D17noPU-v2/GEN-SIM-DIGI-RAW", CP, cmssw, geometry, 0, 9000, 0, 2); // not 1 tau but 2!
  Datasets SingleTauPU200_2023D17("SingleTau_L1TPU200", "RelValSingleTauFlatPt2To100_pythia8_PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200_93X", "/RelValSingleTauFlatPt2To100_pythia8/CMSSW_9_3_7-PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/GEN-SIM-DIGI-RAW", CP, cmssw, geometry, 0, 9000, 0, 2); // not 1 tau but 2!

  Datasets SingleNeutrinoPU140_PhaseIIFall17D("SingleNeutrino_L1TPU140", "SingleNeutrino_PhaseIIFall17D_L1TPU140_93X", "/SingleNeutrino" + path_PU140_v1, CP, cmssw, geometry, 140, 500000, 0, 0);
  Datasets SingleNeutrinoPU200_PhaseIIFall17D("SingleNeutrino_L1TPU200", "SingleNeutrino_PhaseIIFall17D_L1TPU200_93X", "/SingleNeutrino" + path_PU200_v1, CP, cmssw, geometry, 200, 498400, 0, 0);

  Datasets GluGluHToTauTauM125NoPU_PhaseIIFall17D("GluGluHToTauTau_14TeV_L1TnoPU", "GluGluHToTauTau_M125_14TeV_powheg_pythia8_PhaseIIFall17D_L1TnoPU_93X", "/GluGluHToTauTau_M125_14TeV_powheg_pythia8" + path_NoPU_v1 , CP, cmssw, geometry, 0, 100000, 25 , 2);
  Datasets GluGluHToTauTauM125PU140_PhaseIIFall17D("GluGluHToTauTau_14TeV_L1TPU140", "GluGluHToTauTau_M125_14TeV_powheg_pythia8_PhaseIIFall17D_L1TPU140_93X", "/GluGluHToTauTau_M125_14TeV_powheg_pythia8" + path_PU140_v1 ,CP, cmssw, geometry, 0, 100000, 25 , 2);
  Datasets GluGluHToTauTauM125PU200_PhaseIIFall17D("GluGluHToTauTau_14TeV_L1TPU200", "GluGluHToTauTau_M125_14TeV_powheg_pythia8_PhaseIIFall17D_L1TPU200_93X", "/GluGluHToTauTau_M125_14TeV_powheg_pythia8" + path_PU140_v1 ,CP, cmssw, geometry, 0, 100000, 25 , 2);

  Datasets HPlus200NoPU_PhaseIIFall17D ( "ChargedHiggs200_14TeV_L1TnoPU"  , "PYTHIA_Tauola_TB_ChargedHiggs200_14TeV_PhaseIIFall17D_L1TnoPU_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV" + path_NoPU_v2, CP, cmssw, geometry, 0,  100000, 37, 1);
  Datasets HPlus200PU140_PhaseIIFall17D ( "ChargedHiggs200_14TeV_L1TPU140"  , "PYTHIA_Tauola_TB_ChargedHiggs200_14TeV_PhaseIIFall17D_L1TPU140_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV" + path_PU140_v2, CP, cmssw, geometry, 140, 100000, 37, 1);
  Datasets HPlus200PU200_PhaseIIFall17D ( "ChargedHiggs200_14TeV_L1TPU200"  , "PYTHIA_Tauola_TB_ChargedHiggs200_14TeV_PhaseIIFall17D_L1TPU200_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV" + path_PU200_v2, CP, cmssw, geometry, 200, 100000, 37, 1);

  Datasets HPlus500NoPU_PhaseIIFall17D ( "ChargedHiggs500_14TeV_L1TnoPU"  , "PYTHIA_Tauola_TB_ChargedHiggs500_14TeV_PhaseIIFall17D_L1TnoPU_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV" + path_NoPU_v2, CP, cmssw, geometry, 0,  100000, 37, 1);
  Datasets HPlus500PU140_PhaseIIFall17D ( "ChargedHiggs500_14TeV_L1TPU140"  , "PYTHIA_Tauola_TB_ChargedHiggs500_14TeV_PhaseIIFall17D_L1TPU140_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV" + path_PU140_v2, CP, cmssw, geometry, 140, 100000, 37, 1);
  Datasets HPlus500PU200_PhaseIIFall17D ( "ChargedHiggs500_14TeV_L1TPU200"  , "PYTHIA_Tauola_TB_ChargedHiggs500_14TeV_PhaseIIFall17D_L1TPU200_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV" + path_PU200_v2, CP, cmssw, geometry, 200, 99000, 37, 1);

  Datasets HPlus1000NoPU_PhaseIIFall17D ( "ChargedHiggs1000_14TeV_L1TnoPU"  , "PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV_PhaseIIFall17D_L1TnoPU_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV" + path_NoPU_v2, CP, cmssw, geometry, 0,  100000, 37, 1);
  Datasets HPlus1000PU140_PhaseIIFall17D ( "ChargedHiggs1000_14TeV_L1TPU140"  , "PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV_PhaseIIFall17D_L1TPU140_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV" + path_PU140_v2, CP, cmssw, geometry, 140, 100000, 37, 1);
  Datasets HPlus1000PU200_PhaseIIFall17D ( "ChargedHiggs1000_14TeV_L1TPU200"  , "PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV_PhaseIIFall17D_L1TPU200_93X" , "/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV" + path_PU200_v2, CP, cmssw, geometry, 200, 100000, 37, 1);

  Datasets SingleElectronNoPU_PhaseIIFall17D("SingleE_L1TnoPU", "SingleE_FlatPt_2to100_PhaseIIFall17D_L1TnoPU_93X","/SingleE_FlatPt-2to100" + path_NoPU_v1, CP, cmssw, geometry, 0, 243571, 0, 0);
  Datasets SingleElectronPU140_PhaseIIFall17D("SingleE_L1TPU140", "SingleE_FlatPt_2to100_PhaseIIFall17D_L1TPU140_93X","/SingleE_FlatPt-2to100" + path_NoPU_v1, CP, cmssw, geometry, 0, 250000, 0, 0);
  Datasets SingleElectronPU200_PhaseIIFall17D("SingleE_L1TPU200", "SingleE_FlatPt_2to100_PhaseIIFall17D_L1TPU200_93X","/SingleE_FlatPt-2to100" + path_NoPU_v1, CP, cmssw, geometry, 0, 248000, 0, 0);

  // datasets_PhaseIIFall17D.push_back(SinglePionNoPU_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(TTNoPU_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(TTPU140_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(TTPU200_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(SingleTauNoPU_2023D17);
  datasets_PhaseIIFall17D.push_back(SingleTauPU200_2023D17);
  datasets_PhaseIIFall17D.push_back(SingleNeutrinoPU140_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(SingleNeutrinoPU200_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(GluGluHToTauTauM125NoPU_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(GluGluHToTauTauM125PU140_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(GluGluHToTauTauM125PU200_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus200NoPU_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus200PU140_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus200PU200_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus500NoPU_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus500PU140_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus500PU200_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus1000NoPU_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus1000PU140_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(HPlus1000PU200_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(SingleElectronNoPU_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(SingleElectronPU140_PhaseIIFall17D);
  datasets_PhaseIIFall17D.push_back(SingleElectronPU200_PhaseIIFall17D);

  return;
}


//****************************************************************************
Datasets Datasets::GetDataset(const string datasetAlias)
//****************************************************************************
{

  Datasets d;
  bool bSuccess = false;

  // First create the tables for the known MC Productions
  // if (datasets_TTI2023Updg14D.size() < 1 ) CreateMcProductions_();   // obsolete
  // if (datasets_PhaseIISpring17D.size() < 1 ) CreateMcProductions_(); // obsolete
  if (datasets_PhaseIIFall17D.size() < 1 ) CreateMcProductions_();


  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_TTI2023Updg14D.size(); iD++){

    if ( datasetAlias.compare( datasets_TTI2023Updg14D.at(iD).alias_) != 0) continue;
    bSuccess = true;
    d = datasets_TTI2023Updg14D.at(iD);
  }

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_PhaseIISpring17D.size(); iD++){

    if ( datasetAlias.compare( datasets_PhaseIISpring17D.at(iD).alias_) != 0) continue;
    bSuccess = true;
    d = datasets_PhaseIISpring17D.at(iD);
  }
  
  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_PhaseIIFall17D.size(); iD++){

    bool foundAlias = datasetAlias.compare( datasets_PhaseIIFall17D.at(iD).alias_);
    bool foundDset  = datasetAlias.compare( datasets_PhaseIIFall17D.at(iD).datasetPath_);

    // Sanity check
    if ( foundAlias != 0 && foundDset !=0) continue;
    bSuccess = true;
    d = datasets_PhaseIIFall17D.at(iD);
  }
    
  if(!bSuccess){
    cout << "=== Datasets::GetDataset() - Unexpected error! Could not find dataset alias \"" << datasetAlias << "\"." << endl;
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
  // if (datasets_TTI2023Updg14D.size() < 1 )   CreateMcProductions_(); // obsolete
  // if (datasets_PhaseIISpring17D.size() < 1 ) CreateMcProductions_(); // obsolete
  if (datasets_PhaseIIFall17D.size() < 1 ) CreateMcProductions_();

  // Ensure that the supplied sample name is valids
  if (!IsValidDatasetName(datasetAlias)) bool isValidAlias = IsValidDatasetAlias(datasetAlias);

  // Check if the selected MC Prodution is valid
  for( int iD = 0 ; iD < (int) datasets_PhaseIIFall17D.size(); iD++){

    // if ( datasetAlias.compare( datasets_TTI2023Updg14D.at(iD).alias_) == 0)
    //   {
    //   datasetPath = datasets_TTI2023Updg14D.at(iD).datasetPath_;
    //   return datasetPath;
    // }

    // if ( datasetAlias.compare( datasets_PhaseIISpring17D.at(iD).alias_) == 0)
    //   {
    // 	datasetPath = datasets_PhaseIISpring17D.at(iD).datasetPath_;
    // 	return datasetPath;
    //   }

    if ( datasetAlias.compare( datasets_PhaseIIFall17D.at(iD).alias_) == 0)
      {
	datasetPath = datasets_PhaseIIFall17D.at(iD).datasetPath_;
	return datasetPath;
      }
    
    if  ( datasetAlias.compare( datasets_PhaseIIFall17D.at(iD).datasetPath_) == 0)
      {
	datasetPath = datasets_PhaseIIFall17D.at(iD).datasetPath_;
	return datasetPath;
      }

  }
  
  // If this is reached then the dataset sample is invalid
  cout << "=== Datasets::GetDatasetPathFromAlias() - Unexpected error! Could not find dataset path for dataset alias \"" << datasetAlias << "\"." << endl;
  cout << "EXIT" << endl;
  exit(1);

  return datasetPath;
}

#endif

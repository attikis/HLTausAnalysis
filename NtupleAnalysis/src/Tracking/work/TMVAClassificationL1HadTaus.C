// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/TMVAGui.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassificationTopRec( TString myMethodList = "", TString fout = "TMVA_TopRec.root")

{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library

   TMVA::Tools::Instance();

   // to get access to the GUI and all tmva macros
   TString thisdir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   gROOT->SetMacroPath(thisdir + ":" + gROOT->GetMacroPath());
   gROOT->ProcessLine(".L TMVAGui.C");

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 1;
   Use["CutsGA"]          = 1;
   Use["CutsSA"]          = 1;
   // 
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // 
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   std::cout<<"Output File:"<<fout<<std::endl;
   TString outfileName( fout );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   
   //Check the scatterplots of the anticorrelated variables.

   //List of TopRec input variables
  if(1){
    // TBranch *IsFatHiggsJet_T;
    // TBranch *IsFatTopJet_T;
    // TBranch *IsFatWJet_T;
    // TBranch *IsFatPartTopJet_T;
    // TBranch *weight_T;
    factory->AddVariable( "seedPt", "p_{T}", "GeV", 'F' );
    factory->AddVariable( "seedChi2", "#chi^{2}", "", 'F' );
    factory->AddVariable( "seedStubs", "N_{stubs}", "", 'I' );

    // factory->AddVariable( "FatJet_Pt", "p_{T}", "GeV", 'F' );
    // factory->AddVariable( "FatJet_Eta", "#eta", "", 'F');
    // factory->AddVariable( "FatJet_Phi", "#phi", "", 'F');
    // factory->AddVariable( "FatJet_Mass", "mass", "GeV", 'F');
    // factory->AddVariable( "FatJet_CSV", "CSV", "", 'F');
    // factory->AddVariable( "FatJet_tau1", "#tau_{1}", "", 'F');
    // factory->AddVariable( "FatJet_tau2","#tau_{2}", "", 'F');
    // factory->AddVariable( "FatJet_tau3", "#tau_{3}", "", 'F');
    // factory->AddVariable( "FatJet_tau21", "#tau_{21}", "", 'F');
    // factory->AddVariable( "FatJet_tau32", "#tau_{32}", "", 'F');
    // factory->AddVariable( "FatJet_PrunedMass", "pruned mass", "GeV", 'F');
    // factory->AddVariable( "FatJet_SDMass", "sd mass",  "GeV", 'F');
    // factory->AddVariable( "FatJet_corrPrunedMass", "corr pruned mass", "GeV", 'F');
    // factory->AddVariable( "FatJet_sdsubjet1_CSV", "sd subjet1 CSV", "", 'F');
    // factory->AddVariable( "FatJet_sdsubjet1_Eta", "sd subjet1 #eta", "", 'F');
    // factory->AddVariable( "FatJet_sdsubjet1_Mass", "sd subjet1 mass", "GeV", 'F');
    // factory->AddVariable( "FatJet_sdsubjet1_Phi", "sd subjet1 #phi", "", 'F');
    // factory->AddVariable( "FatJet_sdsubjet1_Pt", "sd subjet1 p_{T}", "GeV", 'F');
    // factory->AddVariable( "FatJet_sdsubjet2_CSV", "sd subjet2 CSV", "", 'F');
    // factory->AddVariable( "FatJet_sdsubjet2_Eta", "sd subjet2 #eta", "", 'F');
    // factory->AddVariable( "FatJet_sdsubjet2_Mass", "sd subjet2 mass", "GeV", 'F');
    // factory->AddVariable( "FatJet_sdsubjet2_Phi", "sd subjet2 #phi", "", 'F');
    // factory->AddVariable( "FatJet_sdsubjet2_Pt", "sd subjet2 p_{T}", "GeV", 'F');
    // factory->AddVariable( "FatJet_hadronFlavour", "hadron flavour", "", 'I');
    // factory->AddVariable( "FatJet_partonFlavour",  "parton flavour", 'I');    
    // factory->AddVariable( "FatJet_nsoftdropSubjets", "N SD subjets", "", 'I');
    // factory->AddVariable( "FatJet_nDaughters", "n daughters", "", 'I');
  }


   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
  // TString pseudomulticrab = "BoostedTopRecoTree_181029_143902_ANSelections_PlotsAndTrees";
  // TString fnameS_TT = pseudomulticrab+"/TT/res/histograms-TT.root";
  // TString fnameB_TT = pseudomulticrab+"/TT/res/histograms-TT.root";
  TString fnameS_TauGun = "histograms-RelValSingleTauFlatPt2To100_pythia8_PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200_93X.root";
  TString fnameB_NeuGun = "histograms-SingleNeutrino_PhaseIIFall17D_L1TPU200_93X.root";
  
  TFile *inputS = TFile::Open(fnameS_TT);
  TFile *inputB = TFile::Open(fnameB_TT);
  
  std::cout << "--- TMVAClassification       : Using input files: " << inputS->GetName() << ", "<<inputB->GetName()<<std::endl;

      // --- Register the training and test trees

   TTree *signal     = (TTree*)inputS->Get("tree");
   TTree *background = (TTree*)inputB->Get("tree");
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,               signalWeight     ); 
   factory->AddBackgroundTree( background,           backgroundWeight ); 

   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   //factory->SetSignalWeightExpression    ("eventWeight");
   //factory->SetBackgroundWeightExpression("eventWeight");                                  //Set weights!

   // Apply additional cuts on the signal and background samples (can be different)

    // TCut mucuts = "LdgTrijetLdgJetBDisc > 0.8 && LdgTrijetLdgJetBDisc<1.01 &&   LdgTrijetSubldgJetBDisc > 0.8 && LdgTrijetSubldgJetBDisc < 1.01 && SubldgTrijetLdgJetBDisc > 0.8 && SubldgTrijetLdgJetBDisc < 1.01 && SubldgTrijetSubldgJetBDisc > 0.8 && SubldgTrijetSubldgJetBDisc < 1.01";

    // TCut mucutb = "LdgTrijetLdgJetBDisc > 0.8 && LdgTrijetLdgJetBDisc<1.01 &&   LdgTrijetSubldgJetBDisc > 0.8 && LdgTrijetSubldgJetBDisc < 1.01 && SubldgTrijetLdgJetBDisc > 0.8 && SubldgTrijetLdgJetBDisc < 1.01 && SubldgTrijetSubldgJetBDisc > 0.8 && SubldgTrijetSubldgJetBDisc < 1.01";

   //TCut mycuts = "IsFatTopJet == 1";
   //TCut mycutb = "IsFatTopJet == 0";
    // Tell the factory how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   // To also specify the number of testing events, use:
     factory->PrepareTrainingAndTestTree( mycuts, mycutb,
					  //"nTrain_Signal=20000:nTrain_Background=20000:NTest_Signal=20000:nTest_Background=20000:SplitMode=Random:NormMode=NumEvents:!V" );
					  "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:NegWeightTreatment=IgnoreNegWeightsInTraining" );


   if (Use["BDT"])  // Adaptive Boost
     factory->BookMethod( TMVA::Types::kBDT, "BDT",
			  "H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:IgnoreNegWeightsInTraining=True" );



   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "H:!V:NTrees=400:MinNodeSize=2.5%:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:IgnoreNegWeightsInTraining=True" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "H:!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate:IgnoreNegWeightsInTraining=True" );

   // if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
   //    factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
   //                         "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:IgnoreNegWeightsInTraining=True" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:IgnoreNegWeightsInTraining=True" );


   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   //   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}

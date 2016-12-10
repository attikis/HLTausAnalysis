#!/usr/bin/env python
'''

Usage:
cd HLTausAnalysis
source setup.csh
source /Users/attikis/ROOT/v5-34-00-patches/bin/thisroot.csh
./plotL1TkTaus.py -m results/test/

Comments:

'''
#================================================================================================
# All imported modules
#================================================================================================
# System
import os
import sys
import numpy
import math
from optparse import OptionParser

# User
import HLTausAnalysis.NtupleAnalysis.tools.plotter as m_plotter
import HLTausAnalysis.NtupleAnalysis.tools.histos as m_histos
from HLTausAnalysis.NtupleAnalysis.tools.l1tktauAux import *
import HLTausAnalysis.NtupleAnalysis.tools.styles as m_styles
import HLTausAnalysis.NtupleAnalysis.tools.aux as m_aux

# ROOT
import ROOT

#================================================================================================
# Options here
#================================================================================================
analysis         = "CaloPlusTracks"
bDoL1TkTau       = True
bDoL1TkTauExtra  = True
bDoMatchTk       = True
bDoSigTks        = True
bDoIsoTks        = True
bDoEfficiencies  = False # Seg-Faults
datasetList      = ["VBF", "MinBias"] #["MinBias", "VBF", "HPlus160"]
saveFormats      = ["png"]
savePath         = "plots/"

#================================================================================================
# Function Definition
#================================================================================================
def CreateDatasetDict(inputPath, analysis, outputExt):
    '''
    '''
    datasetPaths= {}
    for dataset in GetDatasetsList():
        datasetPaths[dataset] = inputPath + analysis + "_Histograms_" + dataset + outputExt + ".root"
    return datasetPaths


def GetDatasetsList():
    '''
    '''
    datasets = ["MinBias", "VBF", "PiPlus", "PiMinus", "SingleTauGun1p",
                "DiTauGun3p", "TTbar", "HPlus160", "HPlus200", "SingleElectron",
                "SinglePositron", "SingleMuPlus", "SingleMuMinus", "SinglePhoton",
                "SingleMuon_NoPU", "SingleMuMinus_Pt_2_10_NoPU", "SingleMuPlus_Pt_2_10_NoPU",
                "SingleMuPlus_NoPU"]    
    return datasets


def PrintPSet(datasetPath, dataset):
    '''
    '''
    p = m_plotter.Plotter( Verbose=False, BatchMode=parseOpts.batchMode )
    p.GetDatasets("TTI2023Upg14D", datasetPath, [dataset])
    p.PrintPSet("TkTauFromCaloNTupleMaker/configInfo/parameterSet")    
    return


def DoPlots(hList, datasetPaths, datasetList, saveExt=""):
    '''
    '''
    p = m_plotter.Plotter( Verbose=False, BatchMode=parseOpts.batchMode )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
    p.SetBoolUseDatasetAsLegEntry(True)
    p.AddHisto(hList)
    p.EnableColourPalette(False)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = False)
    #p.GetTLegend().SetHeader( "" )
    p.SaveHistos(True, savePath, saveFormats, saveExt)
    return


def DoEfficiency(hList, datasetPaths, datasetList, cutDirList, saveExt=""):
    '''
    '''
    for cutDir in cutDirList:
        
        p = m_plotter.Plotter( Verbose=False, BatchMode=parseOpts.batchMode )
        for dataset in datasetList:
            p.AddDataset(dataset, datasetPaths[dataset])
        p.SetBoolUseDatasetAsLegEntry(True)
        p.AddHisto(hList)
        p.EnableColourPalette(False)
        p.DrawEfficiency(cutDir, "binomial")
        p.GetTLegend().SetHeader( "L1TkTau" )
        p.SaveHistos(True, savePath, saveFormats, saveExt)        
    return


def main(opts):

    datasetPaths = CreateDatasetDict(opts.mcrab, analysis, "")

    if opts.verbose:
        PrintPSet("/Users/attikis/disk/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v7/test/", "VBF")

        
    if bDoL1TkTau:
        if opts.verbose:
            print "=== Doing L1TkTau"        
        # DoPlots( hL1TkTau_Multiplicity, datasetPaths, datasetList)        
        DoPlots( hL1TkTau_CHF        , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_NHF        , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_NHFAbs     , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_InvMass    , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_Rtau       , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_RelIso     , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_VtxIso     , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_VtxIsoAbs  , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_InvMassIncl, datasetPaths, datasetList) 
        DoPlots( hL1TkTau_ResolutionCaloEt , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_ResolutionCaloEta, datasetPaths, datasetList) 
        DoPlots( hL1TkTau_ResolutionCaloPhi, datasetPaths, datasetList) 
  
    if bDoL1TkTauExtra:
        if opts.verbose:
            print "=== Doing L1TkTauExtra"
        DoPlots( hL1TkTau_SigConeRMin, datasetPaths, datasetList) 
        DoPlots( hL1TkTau_SigConeRMax, datasetPaths, datasetList) 
        DoPlots( hL1TkTau_IsoConeRMin, datasetPaths, datasetList) 
        DoPlots( hL1TkTau_IsoConeRMax, datasetPaths, datasetList) 
        DoPlots( hL1TkTau_DeltaRGenP , datasetPaths, datasetList) 
        DoPlots( hL1TkTau_Charge     , datasetPaths, datasetList) 

        
    if bDoEfficiencies:
        if opts.verbose:
            print "=== Doing Efficiencies"
        DoEfficiency( hL1TkTau_CHF                  , datasetPaths, datasetList, ["<"] )
        DoEfficiency( hL1TkTau_NHF                  , datasetPaths, datasetList, ["<", ">"] )
        DoEfficiency( hL1TkTau_NHFAbs               , datasetPaths, datasetList, ["<", ">"] )
        DoEfficiency( hL1TkTau_InvMass              , datasetPaths, datasetList, ["<"] )
        DoEfficiency( hL1TkTau_Rtau                 , datasetPaths, datasetList, [">", "<"] )
        DoEfficiency( hL1TkTau_RelIso               , datasetPaths, datasetList, ["<"] )
        DoEfficiency( hL1TkTau_VtxIsoAbs            , datasetPaths, datasetList, [">"] )

        
    if bDoMatchTk:
        if opts.verbose:
            print "=== Doing MatchTks"
        DoPlots( hL1TkTau_MatchTk_PtRel         , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_DeltaR        , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_Pt            , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_d0Abs         , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_ChiSquared    , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_RedChiSquared , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_NStubs        , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_NPsStubs      , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_NBarrelStubs  , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_NEndcapStubs  , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_Eta           , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_POCAz         , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_d0            , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_IsGenuine     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_IsUnknown     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_IsCombinatoric, datasetPaths, datasetList)
        DoPlots( hL1TkTau_MatchTk_PtMinusCaloEt , datasetPaths, datasetList)
        
    if bDoSigTks:
        if opts.verbose:
            print "=== Doing SigTks"
        DoPlots( hL1TkTau_NSigTks          , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_Pt        , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_Eta       , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_d0        , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_d0Abs     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_d0Sig     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_d0SigAbs  , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_PtRel     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_POCAz     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_DeltaPOCAz, datasetPaths, datasetList)
        DoPlots( hL1TkTau_SigTks_StubPtCons, datasetPaths, datasetList)
        
        
    if bDoIsoTks:
        if opts.verbose:
            print "=== Doing IsoTks"
        DoPlots( hL1TkTau_NIsoTks          , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_Pt        , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_Eta       , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_d0        , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_d0Abs     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_d0Sig     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_d0SigAbs  , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_PtRel     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_POCAz     , datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_DeltaPOCAz, datasetPaths, datasetList)
        DoPlots( hL1TkTau_IsoTks_StubPtCons, datasetPaths, datasetList)
    return

        
#================================================================================================
# Main
#================================================================================================
if __name__ == "__main__":    

    parser = OptionParser(usage="Usage: %prog [options]" , add_help_option=False,conflict_handler="resolve")
    parser.add_option("-m", "--mcrab"    , dest="mcrab"    , action="store", help="Path to the multicrab directory for input")
    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=True, help="Enables batch mode (canvas creation does NOT generates a window)")
    parser.add_option("-v", "--verbose"  , dest="verbose"  , action="store_true", default=False, help="Enables verbose mode (for debugging purposes)")
    (parseOpts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if parseOpts.mcrab == None:
        print "Not enough arguments passed to script execution. Printing docstring & EXIT."
        print __doc__
        sys.exit(0)
    else:
        pass

    # Program execution
    main(parseOpts)

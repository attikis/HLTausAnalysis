#!/usr/bin/env python
'''

Usage:
cd HLTausAnalysis
source setup.csh
source /Users/attikis/ROOT/v5-34-00-patches/bin/thisroot.csh
./plotTrigger.py -m results/test/

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
from HLTausAnalysis.NtupleAnalysis.tools.tauTriggerAux import *
import HLTausAnalysis.NtupleAnalysis.tools.styles as m_styles
import HLTausAnalysis.NtupleAnalysis.tools.aux as m_aux
import HLTausAnalysis.NtupleAnalysis.tools.datasets as m_datasets

# ROOT
import ROOT

#================================================================================================
# Options here
#================================================================================================
analysis         = "CaloPlusTracks"
bTurnOns         = False
bSingleTau       = True
bDiTau           = False
bDiTau_Indist    = False
bDiTau_Dist_Calo = False
bDiTau_Dist_Tk   = False
datasetList      = ["VBF"] #["MinBias", "VBF", "HPlus160"]
saveFormats      = ["png"]
savePath         = "plots/"

#================================================================================================
# Function Definition
#================================================================================================
def CreateDatasetDict(inputPath, analysis, outputExt):
    '''
    '''
    # inputPath = os.getcwd() + inputPath
    if not inputPath.endswith("/"):
        inputPath = inputPath + "/"
        
    datasetPaths= {}
    for dataset in GetDatasetsList():
        rFile = inputPath + analysis + "_Histograms_" + dataset + outputExt + ".root" 
        datasetPaths[dataset] = rFile
        #print "[%s] = %s" % (dataset, rFile)
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


def PrintTreeCreationParameters():
    
    inputPath = "/Users/attikis/my_work/cms/lxplus/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v2/SigTk5_R01to04_TTI2023Upg14D_DefaultTrackingSequence_08Aug2014_164642/"
    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.GetDatasets("TTI2023Upg14D", inputPath, ["minbias"])
    p.PrintPSet("TkTauFromCaloAnalyzer/configInfo/parameterSet")
    return


def DoPlots(hList, datasetPaths, datasetList, saveExt="", bLegHeader=None):

    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
    p.AddHisto(hList)
    p.SetBoolUseDatasetAsLegEntry(False)
    p.EnableColourPalette(True)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = False)
    if bLegHeader != None:
        p.GetTLegend().SetHeader(bLegHeader)
    p.SaveHistos(True, savePath, saveFormats, saveExt)
    return


def DoROC(rateToEffMap, datasetPaths, signalDataset, rocSaveName, bSaveAuxHistos=False):
    
    p0 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    RateHistoList, EffHistoList = p0.GetROCHistoLists(rateToEffMap)

    ### Rates
    p1 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p1.SetBoolUseDatasetAsLegEntry(False)
    p1.AddDataset("MinBias", datasetPaths["MinBias"])
    p1.AddHisto(RateHistoList)
    p1.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p1.SaveHistos(bSaveAuxHistos, savePath, saveFormats)

    ### Efficiencies
    p2 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p2.SetBoolUseDatasetAsLegEntry(False)
    p2.AddDataset(signalDataset, datasetPaths[signalDataset])
    p2.AddHisto(EffHistoList)
    p2.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p2.SaveHistos(bSaveAuxHistos, savePath, saveFormats)
    
    ### ROCs
    p3 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p3.SetBoolUseDatasetAsLegEntry(False)
    p3.ConvertToROC( p1.GetHistos(), p2.GetHistos(), decimals = 3, bDrawZValues=True, **ROC)
    p3.DrawMultigraph(rocSaveName, **ROC)
    p3.SetTLegendHeader(signalDataset, "")
    p3.SaveHistos(True, savePath, saveFormats)
    return


def main(opts):

    datasetPaths = CreateDatasetDict(opts.mcrab, analysis, "")

    if opts.verbose:
        PrintPSet("/Users/attikis/disk/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v7/test/", "VBF")

    # DoPlots( hHepMCEvt_VtxZ   , datasetList)
    # DoPlots( hHepMCEvt_VtxX_VtxY   , datasetList)

    if bTurnOns:
        if opts.verbose:
            print "=== Doing Turn-Ons"
        DoPlots( SingleTau_TurnOn      , datasetPaths, datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        DoPlots( SingleTau_TurnOn_50KHz, datasetPaths, datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        DoPlots( DiTau_TurnOn          , datasetPaths, datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        DoPlots( DiTau_TurnOn_50KHz    , datasetPaths, datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )

        
    if bSingleTau:
        if opts.verbose:
            print "=== Doing SingleTau"
        DoPlots( SingleTau_Rate  , datasetPaths, ["MinBias"], "")
        DoPlots( SingleTau_Eff   , datasetPaths, datasetList, "")
        DoROC  ( SingleTau_ROCs  , datasetPaths, datasetList[0], "SingleTau_ROCs" + "")


    if bDiTau:
        if opts.verbose:
            print "=== Doing DiTau"
        #for h in DiTau_Rate_CaloIso:
        #    DoPlots( h, ["MinBias"])
        #for h in DiTau_Eff_CaloIso:
        #    DoPlots( h , datasetList)
        DoPlots( DiTau_Rate   , ["MinBias"], "")
        DoPlots( DiTau_Eff    , datasetPaths, datasetList, "")
        DoROC  ( DiTau_ROCs_TP, datasetList[0], "DiTau_ROCs" + "")


    if bDiTau_Indist:
        if opts.verbose:
            print "=== Doing DiTau (Indistinguishable)"

        ### Turn-On Curves
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        ### ROCs
        DoPlots( DiTau_Rate  , ["MinBias"])
        DoPlots( DiTau_Eff   , datasetList)
        DoROC  ( DiTau_ROCs  , datasetList[0], "DiTau_ROCs_Indist")
    

    if bDiTau_Dist_Calo:
        if opts.verbose:
            print "=== Doing DiTau (Distinguishable)"
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        for h in DiTau_Rate_CaloIso:
            DoPlots( h, ["MinBias"])
        for h in DiTau_Eff_CaloIso:
            DoPlots( h , datasetList)
        DoROC  ( DiTau_ROCs_CaloIso, datasetList[0], "DiTau_ROCs_CaloIso")


    if bDiTau_Dist_Tk:
        if opts.verbose:
            print "=== Doing DiTau (Distinguishable, Tk-Other)"

        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        for h in DiTau_Rate_TkIso:
            DoPlots( h, ["MinBias"])
        for h in DiTau_Eff_TkIso:
            DoPlots( h , datasetList)
        DoROC  ( DiTau_ROCs_TkIso, datasetList[0], "DiTau_ROCs_TkIso")

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

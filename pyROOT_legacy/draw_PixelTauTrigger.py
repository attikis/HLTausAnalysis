#!/usr/bin/env python
###############################################################
# All imported modules
###############################################################
import ROOT
import os
import sys
import numpy
import math
import tools.plotter as m_plotter
import tools.histos as m_histos
from tools.tauTriggerAux import *
import tools.datasets as m_datasets
import tools.styles as m_styles
import tools.aux as m_aux

###############################################################
# Options
###############################################################
bUseTTPixelTks   = False
bTurnOns         = False
bSingleTau       = False
bDiTau           = True
bDiTau_Indist    = False
bDiTau_Dist_Calo = False
bDiTau_Dist_Tk   = False
bColourPalette   = True

###############################################################
# General Settings
###############################################################
datasetList      = ["VBF"]
saveFormats      = ["png", "pdf"]
if bUseTTPixelTks:
    tkType = "pixTks"
else:
    tkType = "ttTks"

#savePath     = "/Users/attikis/Desktop/HLTaus/TauTrigger/%s/" % (tkType)
savePath      = "/Users/attikis/postdoc/HLTaus/2016/CMSUpgradeWeek_13Sept2016/figures/TauTrigger/" + tkType + "/"

###############################################################
# Helper Functions
###############################################################
def CreateDatasetDict(inputPath, outputExt):
    '''
    '''
    datasetPaths= {}
    for dataset in GetDatasetsList():
        datasetPaths[dataset] = inputPath + dataset + outputExt + ".root"
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


###############################################################
# Main
###############################################################
def PrintTreePSet():
    '''
    '''
    baseDir   = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v3/"
    inputPath = baseDir + "AllTracks_CaloTausCorr_FitParam4_TTI2023Upg14D_17Feb2015_165244/"
    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.GetDatasets("TTI2023Upg14D", inputPath, ["minbias"])
    p.PrintPSet("TkTauFromCaloAnalyzer/configInfo/parameterSet")
    return


def DoPlots(hList, datasetList, saveExt="", bLegHeader=None):
    '''
    '''
    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
    p.AddHisto(hList)
    p.SetBoolUseDatasetAsLegEntry(False)
    p.EnableColourPalette(True)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = False)
    if bLegHeader != None:
        p.GetTLegend().SetHeader(bLegHeader)
    else:
        p.SetTLegendHeader(dataset, "" )
    p.SaveHistos(True, savePath, saveFormats, saveExt)
    return


def DoROC(rateToEffMap, signalDataset, rocSaveName, tkType=None, bSaveAuxHistos=False):
    '''
    '''
    p0 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    RateHistoList, EffHistoList = p0.GetROCHistoLists(rateToEffMap)

    # Rates
    p1 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p1.SetBoolUseDatasetAsLegEntry(False)
    p1.AddDataset("MinBias", datasetPaths["MinBias"])
    p1.AddHisto(RateHistoList)
    p1.EnableColourPalette(bColourPalette)
    p1.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    if tkType != None:
        p1.GetTLegend().SetHeader(tkType)
    else:
        p1.SetTLegendHeader("MinBias", "" )        
    p1.SaveHistos(bSaveAuxHistos, savePath, saveFormats)

    # Efficiencies
    p2 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p2.SetBoolUseDatasetAsLegEntry(False)
    p2.AddDataset(signalDataset, datasetPaths[signalDataset])
    p2.AddHisto(EffHistoList)
    p2.EnableColourPalette(bColourPalette)
    p2.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    if tkType != None:
        p2.GetTLegend().SetHeader(tkType)
    else:
        p2.SetTLegendHeader(signalDataset, "" )        
    p2.SaveHistos(bSaveAuxHistos, savePath, saveFormats)
    
    # ROCs
    p3 = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p3.SetBoolUseDatasetAsLegEntry(False)
    p3.ConvertToROC( p1.GetHistos(), p2.GetHistos(), decimals = 3, bDrawZValues=True, **ROC)
    p3.DrawMultigraph(rocSaveName, **ROC)
    p3.SetTLegendHeader(signalDataset, "")
    p3.SaveHistos(True, savePath, saveFormats)
    return

###############################################################
# Main Functions
###############################################################
if __name__ == "__main__":

    if bUseTTPixelTks:
        datasetPaths = CreateDatasetDict("Macros/PixelTauTrigger/results/pixTks/Default/PixelTauTrigger_Histograms_", "")
    else:
        datasetPaths = CreateDatasetDict("Macros/PixelTauTrigger/results/ttTks/PixelTauTrigger_Histograms_", "")

    # Turn-On Curves
    if bTurnOns:
        DoPlots( SingleTau_TurnOn      , datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        DoPlots( DiTau_TurnOn          , datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        # DoPlots( SingleTau_TurnOn_50KHz, datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        # DoPlots( DiTau_TurnOn_50KHz    , datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )

    # SingleTau
    if bSingleTau:
        DoPlots( SingleTau_Rate  , ["MinBias"], "", "MinBias")
        DoPlots( SingleTau_Eff   , datasetList, "")
        DoROC  ( SingleTau_ROCs  , datasetList[0], "SingleTau_ROCs", None, False)


    # DiTau
    if bDiTau:
        DoPlots( DiTau_Rate   , ["MinBias"], "", "MinBias")
        DoPlots( DiTau_Eff    , datasetList, "")
        DoROC  ( DiTau_ROCs_TP, datasetList[0], "DiTau_ROCs")

    # DiTau: Indistinguishable
    if bDiTau_Indist:
        ### Turn-On Curves
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        ### ROCs
        DoPlots( DiTau_Rate  , ["MinBias"], "", "MinBias")
        DoPlots( DiTau_Eff   , datasetList)
        DoROC  ( DiTau_ROCs  , datasetList[0], "DiTau_ROCs_Indist")
    
    # DiTau: Distinguishable (Calo-Other)
    if bDiTau_Dist_Calo:
        ### Turn-On Curves
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        ### ROCs
        for h in DiTau_Rate_CaloIso:
            DoPlots( h, ["MinBias"], "", "MinBias")
        for h in DiTau_Eff_CaloIso:
            DoPlots( h , datasetList)
        DoROC  ( DiTau_ROCs_CaloIso, datasetList[0], "DiTau_ROCs_CaloIso")

    # DiTau: Distinguishable (Tk-Other)
    if bDiTau_Dist_Tk:
        # Turn-On Curves
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        # ROCs
        for h in DiTau_Rate_TkIso:
            DoPlots( h, ["MinBias"], "", "MinBias")
        for h in DiTau_Eff_TkIso:
            DoPlots( h , datasetList)
        DoROC  ( DiTau_ROCs_TkIso, datasetList[0], "DiTau_ROCs_TkIso")    
 
###############################################################

#!/usr/bin/env python
###############################################################
### All imported modules
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
### Options here
###############################################################
bTurnOns         = True
bSingleTau       = True
bDiTau           = True
bDiTau_Indist    = False
bDiTau_Dist_Calo = False
bDiTau_Dist_Tk   = False
bUseTTPixelTks   = True
nFitParams       = "5FitParams/"
datasetList      = ["VBF"]
saveFormats      = ["png"]#, "pdf"]

if bUseTTPixelTks:
    tkType = "pixTks"
else:
    tkType = "ttTks"

#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/5FitParams/" + tkType + "/optimisation/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/" + nFitParams + tkType + "/default/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/vtxIso0p2cm/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/tkm_vLooseWP/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/trial5/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/trial6/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/finalTrial/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/" + tkType + "/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/"

#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/5FitParams/" + tkType + "/optimisation/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/5FitParams/" + tkType + "/optimisation/new/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/5FitParams/" + tkType + "/optimisation/finalTrial/"
#inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/TauTrigger/results/5FitParams/" + tkType + "/default/"

inputPath   = "/Users/attikis/my_work/cms/lxplus/hltaus/pyROOT/Macros/PixelTauTrigger/results/" + tkType + "/Default/"
savePath    = "/Users/attikis/Desktop/HLTaus/"
#savePath    = "/Users/attikis/talks/Phase2_Trigger_Workshop_17March2015/figures/5FitParams/" + tkType + "/" + datasetList[0] + "/"

datasetPaths   = {}
datasetPaths["MinBias"]           = inputPath + "TauTrigger_Histograms_MinBias.root"       
datasetPaths["PiPlus"]            = inputPath + "TauTrigger_Histograms_PiPlus.root"        
datasetPaths["PiMinus"]           = inputPath + "TauTrigger_Histograms_PiMinus.root"        
datasetPaths["SingleTauGun1p"]    = inputPath + "TauTrigger_Histograms_SingleTauGun1p.root"
datasetPaths["DiTauGun3p"]        = inputPath + "TauTrigger_Histograms_DiTauGun3p.root"    
datasetPaths["VBF"]               = inputPath + "TauTrigger_Histograms_VBF.root"           
datasetPaths["TTbar"]             = inputPath + "TauTrigger_Histograms_TTBar.root"
datasetPaths["HPlus160"]          = inputPath + "TauTrigger_Histograms_HPlus160.root"
datasetPaths["HPlus200"]          = inputPath + "TauTrigger_Histograms_HPlus200.root"

###############################################################
### Main
###############################################################
def PrintTreeCreationParameters():
    
    inputPath = "/Users/attikis/my_work/cms/lxplus/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v2/SigTk5_R01to04_TTI2023Upg14D_DefaultTrackingSequence_08Aug2014_164642/"
    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.GetDatasets("TTI2023Upg14D", inputPath, ["minbias"])
    p.PrintPSet("TkTauFromCaloAnalyzer/configInfo/parameterSet")
    
    return

###############################################################
def DoPlots(hList, datasetList, saveExt="", bLegHeader=None):

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

###############################################################
def DoROC(rateToEffMap, signalDataset, rocSaveName, bSaveAuxHistos=False):
    
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

###############################################################
if __name__ == "__main__":

    ### Miscellaneous
    #DoPlots( hHepMCEvt_VtxZ   , datasetList)
    #DoPlots( hHepMCEvt_VtxX_VtxY   , datasetList)

    if bTurnOns:
        DoPlots( SingleTau_TurnOn      , datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        DoPlots( SingleTau_TurnOn_50KHz, datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        DoPlots( DiTau_TurnOn          , datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )
        DoPlots( DiTau_TurnOn_50KHz    , datasetList, "",  m_datasets.DatasetClass().ConvertDatasetToLatex(datasetList[0]) )

    ### SingleTau
    if bSingleTau:
        DoPlots( SingleTau_Rate  , ["MinBias"], "")
        DoPlots( SingleTau_Eff   , datasetList, "")
        DoROC  ( SingleTau_ROCs  , datasetList[0], "SingleTau_ROCs" + "")


    ### DiTau
    if bDiTau:
        #for h in DiTau_Rate_CaloIso:
        #    DoPlots( h, ["MinBias"])
        #for h in DiTau_Eff_CaloIso:
        #    DoPlots( h , datasetList)
        DoPlots( DiTau_Rate   , ["MinBias"], "")
        DoPlots( DiTau_Eff    , datasetList, "")
        DoROC  ( DiTau_ROCs_TP, datasetList[0], "DiTau_ROCs" + "")

    ### DiTau: Indistinguishable
    if bDiTau_Indist:
        ### Turn-On Curves
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        ### ROCs
        DoPlots( DiTau_Rate  , ["MinBias"])
        DoPlots( DiTau_Eff   , datasetList)
        DoROC  ( DiTau_ROCs  , datasetList[0], "DiTau_ROCs_Indist")
    
    ### DiTau: Distinguishable (Calo-Other)
    if bDiTau_Dist_Calo:
        ### Turn-On Curves
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        ### ROCs
        for h in DiTau_Rate_CaloIso:
            DoPlots( h, ["MinBias"])
        for h in DiTau_Eff_CaloIso:
            DoPlots( h , datasetList)
        DoROC  ( DiTau_ROCs_CaloIso, datasetList[0], "DiTau_ROCs_CaloIso")

    ### DiTau: Distinguishable (Tk-Other)
    if bDiTau_Dist_Tk:
        ### Turn-On Curves
        DoPlots( DiTau_TurnOn      , datasetList)
        DoPlots( DiTau_TurnOn_50KHz, datasetList)
        ### ROCs
        for h in DiTau_Rate_TkIso:
            DoPlots( h, ["MinBias"])
        for h in DiTau_Eff_TkIso:
            DoPlots( h , datasetList)
        DoROC  ( DiTau_ROCs_TkIso, datasetList[0], "DiTau_ROCs_TkIso")    
 
###############################################################

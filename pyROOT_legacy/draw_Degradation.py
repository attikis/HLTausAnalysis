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
import tools.styles as m_styles
import tools.aux as m_aux

###############################################################
### Options here
###############################################################
bBatchMode     = True
bVerbose       = False
bSavePlots     = True
inputPath_def  = "/Users/attikis/talks/TauTrigger_26Sep2014/figures/Optimised/SingleTau_Histograms_%s.root"
inputPath_3s   = "/Users/attikis/talks/TauTrigger_26Sep2014/figures/Optimised/Degraded/3sigma/SingleTau_Histograms_%s.root"
inputPath_5s   = "/Users/attikis/talks/TauTrigger_26Sep2014/figures/Optimised/Degraded/5sigma/SingleTau_Histograms_%s.root"
inputPath_10s  = "/Users/attikis/talks/TauTrigger_26Sep2014/figures/Optimised/Degraded/10sigma/SingleTau_Histograms_%s.root"
histoFolder    = ""
#mySavePath     = "/Users/attikis/talks/TauTrigger_26Sep2014/figures/Degraded/"
mySavePath     = ""
mySaveFormats  = ["png", "pdf"]
dataset        = "VBF" #VBF, HPlus160

############################################################### 
### Histogram Options
############################################################### 
yMin       = 1E-05
yMax       = 1E0
yMinRatio  = 0.5
yMaxRatio  = 1.5
bRatio     = False
bInvRatio  = False
ratioLabel = "ratio"
normFactor = None

Rate   = {"xLabel": "E_{T}"             , "xUnits": "GeV", "xMin": 0.0 , "xMax": +100.0, "binWidthX": None, "xCutLines": []  , "xCutBoxes": [[0.0, 15.0, ROOT.kBlack]], "gridX": True, "logX": False,
          "yLabel": "Rate (kHz) / %0.0f", "yUnits": "kHz", "yMin": 1E+1, "yMax": 1E+03   , "binWidthY": None, "yCutLines": [50], "yCutBoxes": []                        , "gridY": True, "logY": True ,
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.68, "yLegMax": 0.82 }

EvtEff = {"xLabel": "E_{T}"             , "xUnits": "GeV", "xMin": 0.0 , "xMax": +100.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [[0.0, 15.0, ROOT.kBlack]], "gridX": True, "logX": False,
          "yLabel": "Efficiency / %0.0f", "yUnits": ""   , "yMin": 0.0 , "yMax": +1.05 , "binWidthY": None, "yCutLines": [], "yCutBoxes": []                        , "gridY": True, "logY": False,
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.68, "yLegMax": 0.82 }

Purity = {"xLabel": "E_{T}"         , "xUnits": "GeV", "xMin": 0.0, "xMax": +100.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [[0.0, 15.0, ROOT.kBlack]], "gridX": True, "logX": False,
          "yLabel": "Purity / %0.0f", "yUnits": ""   , "yMin": 0.0, "yMax": +1.07 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
          "ratioLabel": "Ratio", "ratio": bRatio  , "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.52, "yLegMax": 0.66 }

ROC    = {"xLabel": "Efficiency"      , "xUnits": ""   , "xMin": 0.00, "xMax": +0.8 , "binWidthX": None, "xCutLines": []     , "xCutBoxes": [], "gridX": True, "logX": False,
          "yLabel": "Background Rate" , "yUnits": "kHz", "yMin": 1E+0, "yMax": +1E+03, "binWidthY": None, "yCutLines": [50.0], "yCutBoxes": [], "gridY": True, "logY": True,
          #"legOptions": "FP", "drawOptions": "APE2", "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.20, "yLegMax": 0.46}
          "legOptions": "FP", "drawOptions": "AP", "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.20, "yLegMax": 0.46}

############################################################### 
### Histograms
############################################################### 
hCalo_Rate         = m_histos.TH1orTH2( histoFolder, "Calo_Rate"        , "Calo"             , "SingleTau_Rate"           , **Rate )
hTkConfirmed_Rate  = m_histos.TH1orTH2( histoFolder, "TkConfirmed_Rate" , "+Tk"              , "SingleTau_Rate_Tk"        , **Rate )
hVtxIso1p00_Rate   = m_histos.TH1orTH2( histoFolder, "VtxIso1p00_Rate"  , "+VtxIso (1.00 cm)", "SingleTau_Rate_VtxIso1p00", **Rate )
hRelIso0p10_Rate   = m_histos.TH1orTH2( histoFolder, "RelIso0p10_Rate"  , "+RelIso (0.10)"   , "SingleTau_Rate_RelIso0p10", **Rate )

hCalo_EvtEff         = m_histos.TH1orTH2( histoFolder, "Calo_EvtEff"        , "Calo"             , "SingleTau_PerEvtEff"           , **EvtEff )
hTkConfirmed_EvtEff  = m_histos.TH1orTH2( histoFolder, "TkConfirmed_EvtEff" , "+Tk"              , "SingleTau_PerEvtEff_Tk"        , **EvtEff )
hVtxIso1p00_EvtEff   = m_histos.TH1orTH2( histoFolder, "VtxIso1p00_EvtEff"  , "+VtxIso (1.00 cm)", "SingleTau_PerEvtEff_VtxIso1p00", **EvtEff )
hRelIso0p10_EvtEff   = m_histos.TH1orTH2( histoFolder, "RelIso0p10_EvtEff"  , "+RelIso (0.10)"   , "SingleTau_PerEvtEff_RelIso0p10", **EvtEff )

hCalo_Purity         = m_histos.TH1orTH2( histoFolder, "Calo_Purity"        , "Calo"             , "SingleTau_PerTauPurity"           , **Purity )
hTkConfirmed_Purity  = m_histos.TH1orTH2( histoFolder, "TkConfirmed_Purity" , "+Tk"              , "SingleTau_PerTauPurity_Tk"        , **Purity )
hVtxIso1p00_Purity   = m_histos.TH1orTH2( histoFolder, "VtxIso1p00_Purity"  , "+VtxIso (1.00 cm)", "SingleTau_PerTauPurity_VtxIso1p00", **Purity )
hRelIso0p10_Purity   = m_histos.TH1orTH2( histoFolder, "RelIso0p10_Purity"  , "+RelIso (0.10)"   , "SingleTau_PerTauPurity_RelIso0p10", **Purity )

############################################################### 
### Histo Lists
############################################################### 
RateList = []
RateList.append(hCalo_Rate)
RateList.append(hTkConfirmed_Rate)
RateList.append(hVtxIso1p00_Rate)
RateList.append(hRelIso0p10_Rate)

PerEvtEff = []
PerEvtEff.append(hCalo_EvtEff)
PerEvtEff.append(hTkConfirmed_EvtEff)
PerEvtEff.append(hVtxIso1p00_EvtEff)
PerEvtEff.append(hRelIso0p10_EvtEff)

ROCHistoMap_Calo = {}
ROCHistoMap_Calo[hCalo_Rate] = hCalo_EvtEff
ROCHistoMap_Tk = {}
ROCHistoMap_Tk[hTkConfirmed_Rate] = hTkConfirmed_EvtEff
ROCHistoMap_VtxIso = {}
ROCHistoMap_VtxIso[hVtxIso1p00_Rate]  = hVtxIso1p00_EvtEff
ROCHistoMap_RelIso = {}
ROCHistoMap_RelIso[hRelIso0p10_Rate]  = hRelIso0p10_EvtEff

###############################################################
### Main
###############################################################
def DoPlots(hList, dataset, saveNameExt):

    p = m_plotter.Plotter( bVerbose, bBatchMode )
    if dataset == "HPlus160":
        datasetFileExt = "liteHiggs"
    elif dataset == "HPlus200":
        datasetFileExt = "heavyHiggs"
    else:
        datasetFileExt = dataset
    p.AddDataset(dataset + ":default" , inputPath_def %  (datasetFileExt) )
    p.AddDataset(dataset + ":3#sigma" , inputPath_3s  %  (datasetFileExt) )
    p.AddDataset(dataset + ":5#sigma" , inputPath_5s  %  (datasetFileExt) )
    p.AddDataset(dataset + ":10#sigma", inputPath_10s %  (datasetFileExt) )
    p.SetBoolUseDatasetAsLegEntry(True)
    p.AddHisto(hList)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p.SetTLegendHeader(dataset, "" )
    #p.AppendToCanvasName("_" + saveNameExt)
    p.SaveHistos(True, mySavePath + dataset + "/", mySaveFormats)

    return

###############################################################
def DoROC(rateToEffMap, signalDataset, saveNameExt):
    
    p0 = m_plotter.Plotter( bVerbose, bBatchMode )
    RateHistoList, EffHistoList = p0.GetROCHistoLists(rateToEffMap)

    ### Background Rate
    p1 = m_plotter.Plotter( bVerbose, bBatchMode )
    p1.AddDataset("MinBias:default" , inputPath_def % ("nugun") )
    p1.AddDataset("MinBias:3#sigma" , inputPath_3s  %  ("nugun") )
    p1.AddDataset("MinBias:5#sigma" , inputPath_5s  %  ("nugun") )
    p1.AddDataset("MinBias:10#sigma", inputPath_10s %  ("nugun") )
    p1.AddHisto(RateHistoList)
    p1.SetBoolUseDatasetAsLegEntry(True)
    p1.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p1.SetTLegendHeader("MinBias", "" ) #saveNameExt
    p1.SaveHistos(True, mySavePath + signalDataset + "/", mySaveFormats)

    ### Per-Event Efficiencies
    p2 = m_plotter.Plotter( bVerbose, bBatchMode )
    if signalDataset == "HPlus160":
        signalDatasetFileExt = "liteHiggs"
    elif signalDataset == "HPlus200":
        signalDatasetFileExt = "heavyHiggs"
    else:
        signalDatasetFileExt = signalDataset

    p2.AddDataset(signalDataset + ":default" , inputPath_def %  ( signalDatasetFileExt ) )
    p2.AddDataset(signalDataset + ":3#sigma" , inputPath_3s  %  ( signalDatasetFileExt ) )
    p2.AddDataset(signalDataset + ":5#sigma" , inputPath_5s  %  ( signalDatasetFileExt ) )
    p2.AddDataset(signalDataset + ":10#sigma", inputPath_10s %  ( signalDatasetFileExt ) )
    p2.SetBoolUseDatasetAsLegEntry(True)
    p2.AddHisto(EffHistoList)
    p2.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p2.SetTLegendHeader(signalDataset, "" ) #saveNameExt
    p2.SaveHistos(True, mySavePath + signalDataset + "/", mySaveFormats)
    
    ### ROC Curves
    p3 = m_plotter.Plotter( bVerbose, bBatchMode )
    p3.AddDataset(signalDataset + ":default" , inputPath_def %  ( signalDatasetFileExt ) )
    p3.AddDataset(signalDataset + ":3#sigma" , inputPath_3s  %  ( signalDatasetFileExt ) )
    p3.AddDataset(signalDataset + ":5#sigma" , inputPath_5s  %  ( signalDatasetFileExt ) )
    p3.AddDataset(signalDataset + ":10#sigma", inputPath_10s %  ( signalDatasetFileExt ) )
    p3.SetBoolUseDatasetAsLegEntry(True)
    p3.ConvertToROC( p1.GetHistos() , p2.GetHistos(), decimals = 3)
    p3.DrawMultigraph("SingleTau_ROC_" + saveNameExt, **ROC)
    p3.SetTLegendHeader(signalDataset, "" ) #saveNameExt
    p3.SaveHistos(bSavePlots, mySavePath + signalDataset + "/", mySaveFormats)
    return

###############################################################
if __name__ == "__main__":

    #DoPlots( hTkConfirmed_Purity, dataset, "Tk")
    #DoPlots( hVtxIso1p00_Purity , dataset, "VtxIso")
    #DoPlots( hRelIso0p10_Purity , dataset, "RelIso")
    
    DoROC( ROCHistoMap_Tk, dataset, "Tk")
    DoROC( ROCHistoMap_VtxIso, dataset, "VtxIso")
    DoROC( ROCHistoMap_RelIso, dataset, "RelIso")

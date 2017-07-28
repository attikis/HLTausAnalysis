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
histoFolder    = ""
inputPath_def  = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/Validation/Validation_Histograms_%s.root"
inputPath_3s   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/Validation/3sigma/Validation_Histograms_%s.root"
inputPath_5s   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/Validation/5sigma/Validation_Histograms_%s.root"
inputPath_10s  = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/Validation/10sigma/Validation_Histograms_%s.root"
mySavePath     = "/Users/attikis/talks/TauTrigger_26Sep2014/figures/Degraded/"
mySaveFormats  = ["png", "pdf"]
datasetList    = ["MinBias"]

############################################################### 
### Histogram Options
############################################################### 
yMin       = 1E-04
yMax       = 1E0
yMinRatio  = 0.0
yMaxRatio  = 5.0
bRatio     = False
bInvRatio  = False
ratioLabel = "ratio"
normFactor = "One"


Pt         = {"xLabel": "p_{T}"         , "xUnits": "GeV/c", "xMin": 0.00 , "xMax": 20.0, "binWidthX": 2.0 , "xCutLines": [20], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
             "yLabel": "Entries / %0.0f", "yUnits": ""     , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
             "ratioLabel": "1/ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": 0.5, "yMaxRatio": 1.5, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
              "xLegMin": 0.7, "xLegMax": 0.85, "yLegMin": 0.68, "yLegMax": 0.82, "styleType": "random"}

POCAz      = {"xLabel": "z_{vtx}"       , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": 1.0  , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.1f", "yUnits": ""    , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
              "xLegMin": 0.7, "xLegMax": 0.85, "yLegMin": 0.68, "yLegMax": 0.82, "styleType": "random"}

Eta        = {"xLabel": "#eta"            , "xUnits": ""     , "xMin": -2.6 , "xMax": +2.6 , "binWidthX": 0.25, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin , "yMax": yMax, "binWidthY": None , "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
              "xLegMin": 0.7, "xLegMax": 0.85, "yLegMin": 0.28, "yLegMax": 0.42, "styleType": "random"}

NStubs     = {"xLabel": "Stubs Multiplicity", "xUnits": "", "xMin": 0.5 , "xMax": 8.5, "binWidthX": 1.0 , "xCutLines": [4, 5], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": yMin, "yMax": yMax, "binWidthY": None, "yCutLines": []    , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" ,
              "xLegMin": 0.2, "xLegMax": 0.35, "yLegMin": 0.68, "yLegMax": 0.82, "styleType": "random"}

StubPtCons = {"xLabel": "p_{T}^{stub} constistency", "xUnits": "", "xMin": 0.0  , "xMax": 200.0, "binWidthX": 10.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": yMin , "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
              "xLegMin": 0.7, "xLegMax": 0.85, "yLegMin": 0.68, "yLegMax": 0.82, "styleType": "random"}

############################################################### 
### Group Histogram (of same binning)
############################################################### 
hL1Tks_Pt                 = m_histos.TH1orTH2( histoFolder, "L1Tks_Pt"               , "Default", "L1Tks_Pt"               , **Pt )
hL1Tks_Eta                = m_histos.TH1orTH2( histoFolder, "L1Tks_Eta"              , "Default", "L1Tks_Eta"              , **Eta )
hL1Tks_POCAz              = m_histos.TH1orTH2( histoFolder, "L1Tks_POCAz"            , "Default", "L1Tks_POCAz"            , **POCAz )
hL1Tks_StubPtConsistency  = m_histos.TH1orTH2( histoFolder, "L1Tks_StubPtConsistency", "Default", "L1Tks_StubPtConsistency", **StubPtCons )
hL1Tks_StubMultiplicity   = m_histos.TH1orTH2( histoFolder, "L1Tks_StubMultiplicity" , "Default", "L1Tks_StubMultiplicity" , **NStubs )

###############################################################
### Main
###############################################################
def DoPlots(histoObjectList, dataset):
    '''
    '''

    p = m_plotter.Plotter(bVerbose, bBatchMode)
    p.EnableColourPalette(True)
    p.AddDataset( dataset + ":default" , inputPath_def %  ("nugun") )
    p.AddDataset( dataset + ":3#sigma" , inputPath_3s  %  ("nugun") )
    p.AddDataset( dataset + ":5#sigma" , inputPath_5s  %  ("nugun") )
    p.AddDataset( dataset + ":10#sigma", inputPath_10s %  ("nugun") )
    p.AddHisto( histoObjectList )
    p.SetBoolUseDatasetAsLegEntry(True)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p.GetTLegend().SetHeader( dataset )
    p.SaveHistos(bSavePlots, mySavePath, mySaveFormats)

    return

###############################################################
if __name__ == "__main__":

    DoPlots( [hL1Tks_Pt]   , datasetList[0] )
    DoPlots( [hL1Tks_Eta]  , datasetList[0] )
    DoPlots( [hL1Tks_POCAz], datasetList[0] )
    DoPlots( [hL1Tks_StubPtConsistency], datasetList[0] )
    DoPlots( [hL1Tks_StubMultiplicity], datasetList[0] )

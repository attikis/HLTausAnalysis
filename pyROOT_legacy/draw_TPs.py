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
pi = 4*math.atan(1)

###############################################################
### Options here
###############################################################
saveFormats   = ["png", "pdf"]
datasetList   = ["VBF"]
inputPath     = "Macros/TauTrigger/Validation_Histograms_"
savePath      = ""

datasetPaths  = {}
datasetPaths["MinBias"]         = inputPath + "MinBias.root"
datasetPaths["VBF"]             = inputPath + "VBF.root"
datasetPaths["PiPlus"]          = inputPath + "PiPlus.root"
datasetPaths["PiMinus"]         = inputPath + "PiMinus.root"
datasetPaths["SingleTauGun1p"]  = inputPath + "SingleTauGun1p.root"
datasetPaths["DiTauGun3p"]      = inputPath + "DiTauGun3p.root"
datasetPaths["TTbar"]           = inputPath + "TTBar.root"
datasetPaths["HPlus160"]        = inputPath + "HPlus160.root"
datasetPaths["HPlus200"]        = inputPath + "HPlus200.root"


############################################################### 
### Histogram Options
############################################################### 
NTks = {
    "xLabel": "Multiplicity"    , "xUnits": "" , "xMin": -0.5 , "xMax": 299.5, "binWidthX": 5.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": 1e-05, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Pt = {
    "xLabel": "p_{T}"           , "xUnits": "GeV/c", "xMin": 0.00 , "xMax": 100.0, "binWidthX": 2.0 , "xCutLines": [15], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": 1E-05, "yMax": +1.0 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Eta = {
    "xLabel": "#eta"            , "xUnits": ""     , "xMin": -3.5 , "xMax": +3.5 , "binWidthX": 0.1 , "xCutLines": [-3.0, 0, +3.0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1E-03, "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Phi = {
    "xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.4  , "xMax": +3.4 , "binWidthX": 0.1, "xCutLines": [-pi, 0, +pi], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1E-03, "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

POCAz = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": 0.1, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1E-05, "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Charge = {
    "xLabel": "Charge (e)"     , "xUnits": "", "xMin": -3.5 , "xMax": 3.5 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-05, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.26, "xLegMax": 0.45, "yLegMin": 0.72, "yLegMax": 0.82 }

d0 = {
    "xLabel": "d_{0}"          , "xUnits": "cm", "xMin": -1.0 , "xMax": +1.0, "binWidthX": 0.05, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1E-05, "yMax": +1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

PdgId = {
    "xLabel": "PdgId"          , "xUnits": "", "xMin": 0.0, "xMax": +250.0, "binWidthX": 1.0 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": "", "yMin": 1.0, "yMax": None  , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }


############################################################### 
### Group Histogram (of same binning)
############################################################### 
histoFolder    = ""
hTPs_Multiplicity = m_histos.TH1orTH2( histoFolder, "TPs_Multiplicity", "TPs", None, **NTks )
hTPs_Pt           = m_histos.TH1orTH2( histoFolder, "TPs_Pt"          , "TPs", None, **Pt )
hTPs_Eta          = m_histos.TH1orTH2( histoFolder, "TPs_Eta"         , "TPs", None, **Eta )
hTPs_Phi          = m_histos.TH1orTH2( histoFolder, "TPs_Phi"         , "TPs", None, **Phi )
hTPs_Charge       = m_histos.TH1orTH2( histoFolder, "TPs_Charge"      , "TPs", None, **Charge )
hTPs_POCAz        = m_histos.TH1orTH2( histoFolder, "TPs_POCAz"       , "TPs", None, **POCAz )
hTPs_d0           = m_histos.TH1orTH2( histoFolder, "TPs_d0"          , "TPs", None, **d0 )
hTPs_PdgId        = m_histos.TH1orTH2( histoFolder, "TPs_PdgId"       , "TPs", None, **PdgId )
hTPs_NMatch       = m_histos.TH1orTH2( histoFolder, "TPs_NMatch"      , "TPs", None, **NTks )
# std::vector<std::vector<int>> TP_TTTrackIndex;

###############################################################
### Main
###############################################################
def DoPlots(hList, datasetList, bColourPalette=False, saveExt=""):

    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
        p.AddHisto(hList)
        p.EnableColourPalette(bColourPalette)
        p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = False)
        if len(hList) > 1:
            p.SetTLegendHeader(datasetList[0], "")
        p.SaveHistos(True, savePath, saveFormats, saveExt)

    return

###############################################################
if __name__ == "__main__":

    DoPlots( [hTPs_Multiplicity], datasetList, True )
    DoPlots( [hTPs_Pt          ], datasetList, True )
    DoPlots( [hTPs_Eta         ], datasetList, True )
    DoPlots( [hTPs_Phi         ], datasetList, True )
    DoPlots( [hTPs_Charge      ], datasetList, True )
    DoPlots( [hTPs_POCAz       ], datasetList, True )
    DoPlots( [hTPs_d0          ], datasetList, True )
    DoPlots( [hTPs_PdgId       ], datasetList, True )
    DoPlots( [hTPs_NMatch      ], datasetList, True )

    

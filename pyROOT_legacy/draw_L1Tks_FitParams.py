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
bDo1DPlots    = True
bDo2DPlots    = False

datasetList   = ["PiPlus"]
inputPath_4   = "Macros/TauTrigger/results/validation/4FitParams/Validation_Histograms_"
inputPath_5   = "Macros/TauTrigger/results/validation/5FitParams/Validation_Histograms_"
mySavePath    = "/Users/attikis/talks/TauTrigger_20Feb2015/figures/4vs5FitParams/" + datasetList[0] + "/"
#mySavePath    = ""
mySaveFormats = ["png", "pdf"]

datasetPaths_4  = {}
datasetPaths_4["MinBias"]         = inputPath_4 + "MinBias.root"
datasetPaths_4["VBF"]             = inputPath_4 + "VBF.root"
datasetPaths_4["PiPlus"]          = inputPath_4 + "PiPlus.root"
datasetPaths_4["PiMinus"]         = inputPath_4 + "PiMinus.root"
datasetPaths_4["SingleTauGun1p"]  = inputPath_4 + "SingleTauGun1p.root"
datasetPaths_4["DiTauGun3p"]      = inputPath_4 + "DiTauGun3p.root"
datasetPaths_4["TTbar"]           = inputPath_4 + "TTBar.root"
datasetPaths_4["HPlus160"]        = inputPath_4 + "HPlus160.root"
datasetPaths_4["HPlus200"]        = inputPath_4 + "HPlus200.root"

datasetPaths_5  = {}
datasetPaths_5["MinBias"]         = inputPath_5 + "MinBias.root"
datasetPaths_5["VBF"]             = inputPath_5 + "VBF.root"
datasetPaths_5["PiPlus"]          = inputPath_5 + "PiPlus.root"
datasetPaths_5["PiMinus"]         = inputPath_5 + "PiMinus.root"
datasetPaths_5["SingleTauGun1p"]  = inputPath_5 + "SingleTauGun1p.root"
datasetPaths_5["DiTauGun3p"]      = inputPath_5 + "DiTauGun3p.root"
datasetPaths_5["TTbar"]           = inputPath_5 + "TTBar.root"
datasetPaths_5["HPlus160"]        = inputPath_5 + "HPlus160.root"
datasetPaths_5["HPlus200"]        = inputPath_5 + "HPlus200.root"

############################################################### 
### Histogram Options
############################################################### 
yMin       = 1E-04
yMax       = 1E+00
yMinRatio  = 0.0
yMaxRatio  = 2.2
bRatio     = True
bInvRatio  = False
ratioLabel = "ratio"
normFactor = "One"

NTks = {
    "xLabel": "Multiplicity"    , "xUnits": "" , "xMin": -0.5 , "xMax": 299.5, "binWidthX": 5.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin , "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Pt = {
    "xLabel": "p_{T}"           , "xUnits": "GeV/c", "xMin": 0.00 , "xMax": 150.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": 1E-05, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Eta = {
    "xLabel": "#eta"            , "xUnits": ""     , "xMin": -3.5 , "xMax": +3.5 , "binWidthX": 0.1 , "xCutLines": [-3.0, 0, +3.0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1E-03, "yMax": 2E-1 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Phi = {
    "xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.4 , "xMax": +3.4 , "binWidthX": 0.1, "xCutLines": [-pi, 0, +pi], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f" , "yUnits": ""       , "yMin": 5E-3 , "yMax": 1E-1 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

POCAz = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": 0.5 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1E-05 , "yMax": 5E-01 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

POCAzSig = {
    "xLabel": "z_{0} / #sigma(z_{0})", "xUnits": "", "xMin": -900.0, "xMax": +900.0, "binWidthX": 50.0 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f"      , "yUnits": "", "yMin": +5E-03, "yMax": 1e-01 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True ,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85 }

ChiSq  = {
    "xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +150.0, "binWidthX": 2.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": ""    , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

RedChiSq = {
    "xLabel": "#chi^{2}/d.o.f.", "xUnits": "", "xMin": -0.5, "xMax": +50.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Charge = {
    "xLabel": "Charge (e)"     , "xUnits": "", "xMin": -1.5 , "xMax": 1.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin , "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.26, "xLegMax": 0.45, "yLegMin": 0.72, "yLegMax": 0.82 }

d0 = {
    "xLabel": "d_{0}"   , "xUnits": "cm", "xMin": -5.0 , "xMax": +5.0 , "binWidthX": 0.05 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": "", "yMin": 1E-05, "yMax": yMax, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }


D0Sig = {
    "xLabel": "d_{0}/#sigma(d_{0})", "xUnits": "", "xMin": -20.0 , "xMax": +20.0, "binWidthX": 0.50, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f"    , "yUnits": "", "yMin": 1E-05, "yMax": 1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85 }

d0Abs = {
    "xLabel": "|d_{0}|" , "xUnits": "cm", "xMin": +0.0 , "xMax": +5.0, "binWidthX": 0.05 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": "", "yMin": 1E-05, "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

NStubs = {
    "xLabel": "Stubs Multiplicity", "xUnits": "", "xMin": -0.5, "xMax": 12.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": 1e-4, "yMax":  2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85 }

NPsStubs = {
    "xLabel": "PS Stubs Multiplicity", "xUnits": "", "xMin": -0.5, "xMax": 9.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"      , "yUnits": "", "yMin": 1e-4, "yMax": 2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85 }

NBStubs = {
    "xLabel": "Barrel Stubs Multiplicity", "xUnits": "", "xMin": -0.5, "xMax": 7.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-3, "yMax": 2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.22, "yLegMax": 0.32 }

NEStubs = {
    "xLabel": "Endcap Stubs Multiplicity", "xUnits": "", "xMin": -0.5, "xMax": 6.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-3, "yMax": 2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.22, "yLegMax": 0.32 }

IsGenuine = {
    "xLabel": "is Genuine"     , "xUnits": "", "xMin": -0.5, "xMax": 1.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-2, "yMax": 1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85}

IsUnknown = {
    "xLabel": "is Unknown"     , "xUnits": "", "xMin": -0.5, "xMax": 1.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-2, "yMax": 1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85}

IsCombinatoric = {
    "xLabel": "is Combinatoric", "xUnits": "", "xMin": -0.5, "xMax": 1.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-2, "yMax": 1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85}

StubPtCons = {
    "xLabel": "p_{T}^{stub} constistency", "xUnits": "", "xMin": 0.0  , "xMax": 200.0, "binWidthX": 5.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-5 , "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.85 }

ChiSq_Z0  = {
    "xLabel": "#chi^{2} / %0.0f", "xUnits": "" , "xMin":  +1.0, "xMax": +50000.0, "binWidthX": 1.0, "gridX": False, "logX": True,
    "yLabel": "z_{0} / %0.2f"  , "yUnits": "cm", "yMin": -30.0, "yMax":    +30.0, "binWidthY": 0.5 , "gridY": True, "logY": False,
    "zLabel": "Entries"        , "zUnits": ""  , "zMin":  +1.0, "zMax":     None, "zCutLines": [], "zCutBoxes": [], "logZ": True, "normaliseTo": None, 
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

ChiSq_Pt  = {
    "xLabel": "#chi^{2} / %0.0f", "xUnits": ""     , "xMin": +1.0, "xMax": +50000.0, "binWidthX": 1.0, "gridX": False, "logX": True,
    "yLabel": "p_{T} / %0.0f"   , "yUnits": "GeV/c", "yMin": +0.0, "yMax":   +100.0, "binWidthY": 2.0, "gridY": False, "logY": False,
    "zLabel": "Entries"         , "zUnits": ""     , "zMin": +1.0, "zMax": None, "zCutLines": [], "zCutBoxes": [], "logZ": True, "normaliseTo": None,
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

ChiSq_Eta  = {
    "xLabel": "#chi^{2} / %0.0f", "xUnits": "", "xMin": +1.0, "xMax": +50000.0, "binWidthX": 1.0, "gridX": False, "logX": True,
    "yLabel": "#eta / %0.2f "   , "yUnits": "", "yMin": -3.0, "yMax":     +3.0, "binWidthY": 0.1, "gridY": True , "logY": False, 
    "zLabel": "Entries"         , "zUnits": "", "zMin": +1.0, "zMax":     None, "zCutLines": [], "zCutBoxes": [], "logZ": True, "normaliseTo": None,
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

NStubs_Eta = {
    "xLabel": "Stubs / % 0.0f", "XUnits": "", "xMin": +0.0, "xMax": +10.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "#eta / %0.2f " , "yUnits": "", "yMin": -3.0, "yMax":  +3.0, "binWidthY": 0.1 , "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "zLabel": "Entries"       , "zUnits": "", "zMin": +1.0, "zMax":  None, "zCutLines": [], "zCutBoxes": [], "logZ": True, "normaliseTo": None,
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

############################################################### 
### Group Histogram (of same binning)
############################################################### 
folder    = ""
hL1Tks_Multiplicity      = m_histos.TH1orTH2( folder, "L1Tks_Multiplicity"     , "TTTracks", None, **NTks           )
hL1Tks_Pt                = m_histos.TH1orTH2( folder, "L1Tks_Pt"               , "TTTracks", None, **Pt             )
hL1Tks_Eta               = m_histos.TH1orTH2( folder, "L1Tks_Eta"              , "TTTracks", None, **Eta            )
hL1Tks_Phi               = m_histos.TH1orTH2( folder, "L1Tks_Phi"              , "TTTracks", None, **Phi            )
hL1Tks_Charge            = m_histos.TH1orTH2( folder, "L1Tks_Charge"           , "TTTracks", None, **Charge         )
hL1Tks_POCAz             = m_histos.TH1orTH2( folder, "L1Tks_POCAz"            , "TTTracks", None, **POCAz          )
hL1Tks_d0                = m_histos.TH1orTH2( folder, "L1Tks_d0"               , "TTTracks", None, **d0             )
hL1Tks_d0Abs             = m_histos.TH1orTH2( folder, "L1Tks_d0Abs"            , "TTTracks", None, **d0Abs          )
hL1Tks_ChiSquared        = m_histos.TH1orTH2( folder, "L1Tks_ChiSquared"       , "TTTracks", None, **ChiSq          )
hL1Tks_RedChiSquared     = m_histos.TH1orTH2( folder, "L1Tks_RedChiSquared"    , "TTTracks", None, **RedChiSq       )
hL1Tks_NStubs            = m_histos.TH1orTH2( folder, "L1Tks_NStubs"           , "TTTracks", None, **NStubs         )
hL1Tks_NPsStubs          = m_histos.TH1orTH2( folder, "L1Tks_NPsStubs"         , "TTTracks", None, **NPsStubs       )
hL1Tks_NBarrelStubs      = m_histos.TH1orTH2( folder, "L1Tks_NBarrelStubs"     , "TTTracks", None, **NBStubs        )
hL1Tks_NEndcapStubs      = m_histos.TH1orTH2( folder, "L1Tks_NEndcapStubs"     , "TTTracks", None, **NEStubs        )
hL1Tks_StubPtConsistency = m_histos.TH1orTH2( folder, "L1Tks_StubPtConsistency", "TTTracks", None, **StubPtCons     )
hL1Tks_IsGenuine         = m_histos.TH1orTH2( folder, "L1Tks_IsGenuine"        , "TTTracks", None, **IsGenuine      )
hL1Tks_IsUnknown         = m_histos.TH1orTH2( folder, "L1Tks_IsUnknown"        , "TTTracks", None, **IsUnknown      )
hL1Tks_IsCombinatoric    = m_histos.TH1orTH2( folder, "L1Tks_IsCombinatoric"   , "TTTracks", None, **IsCombinatoric )

### 2D
hL1Tks_ChiSquared_Vs_POCAz = m_histos.TH1orTH2( folder, "L1Tks_ChiSquared_Vs_POCAz", "TTTracks", None, **ChiSq_Z0   )
hL1Tks_ChiSquared_Vs_Pt    = m_histos.TH1orTH2( folder, "L1Tks_ChiSquared_Vs_Pt"   , "TTTracks", None, **ChiSq_Pt   )
hL1Tks_ChiSquared_Vs_Eta   = m_histos.TH1orTH2( folder, "L1Tks_ChiSquared_Vs_Eta"  , "TTTracks", None, **ChiSq_Eta  )
hL1Tks_NStubs_Vs_Eta       = m_histos.TH1orTH2( folder, "L1Tks_NStubs_Vs_Eta"      , "TTTracks", None, **NStubs_Eta )


###############################################################
### Main
###############################################################
def DoPlots(hList, datasetList, datasetPaths, bColourPalette=False, saveExt=""):

    p = m_plotter.Plotter(False, True)
    for dataset in datasetList:
        p.AddDataset( dataset, datasetPaths[dataset] )
        p.AddHisto(hList)
        p.SetBoolUseDatasetAsLegEntry(True)
        p.EnableColourPalette(bColourPalette)
        p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = False)
        p.GetTLegend().SetHeader( "" )
        p.SaveHistos(True, mySavePath, mySaveFormats, saveExt)
    return

###############################################################
def DoPlotsCompare(hList, datasetList, bColourPalette=False):

    p = m_plotter.Plotter(False, True)
    for dataset in datasetList:
        p.AddDataset( dataset + ":4-Parameters"   , datasetPaths_4[dataset] )
        p.AddDataset( dataset + ":5-Parameters", datasetPaths_5[dataset] )
        p.AddHisto(hList)
        p.SetBoolUseDatasetAsLegEntry(True)
        p.EnableColourPalette(bColourPalette)
        p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = False)
        p.GetTLegend().SetHeader( "" )
        p.SaveHistos(True, mySavePath, mySaveFormats)

    return

###############################################################
if __name__ == "__main__":

    if bDo1DPlots:
        DoPlotsCompare( [hL1Tks_Multiplicity  ], datasetList, True )
        DoPlotsCompare( [hL1Tks_Pt            ], datasetList, True )
        DoPlotsCompare( [hL1Tks_Eta           ], datasetList, True )
        DoPlotsCompare( [hL1Tks_Phi           ], datasetList, True )
        DoPlotsCompare( [hL1Tks_Charge        ], datasetList, True )
        DoPlotsCompare( [hL1Tks_POCAz         ], datasetList, True )
        DoPlotsCompare( [hL1Tks_d0            ], datasetList, True )
        DoPlotsCompare( [hL1Tks_d0Abs         ], datasetList, True )
        DoPlotsCompare( [hL1Tks_ChiSquared    ], datasetList, True )
        DoPlotsCompare( [hL1Tks_RedChiSquared ], datasetList, True )
        DoPlotsCompare( [hL1Tks_NStubs        ], datasetList, True )
        DoPlotsCompare( [hL1Tks_NPsStubs      ], datasetList, True )
        DoPlotsCompare( [hL1Tks_NBarrelStubs  ], datasetList, True )
        DoPlotsCompare( [hL1Tks_NEndcapStubs  ], datasetList, True )
        DoPlotsCompare( [hL1Tks_IsGenuine     ], datasetList, True )
        DoPlotsCompare( [hL1Tks_IsUnknown     ], datasetList, True )
        DoPlotsCompare( [hL1Tks_IsCombinatoric], datasetList, True )

    ### 2D
    if bDo2DPlots:
        DoPlots( [hL1Tks_ChiSquared_Vs_POCAz], datasetList, datasetPaths_4, True, saveExt="_4FitParams")
        DoPlots( [hL1Tks_ChiSquared_Vs_Pt   ], datasetList, datasetPaths_4, True, saveExt="_4FitParams")
        DoPlots( [hL1Tks_ChiSquared_Vs_Eta  ], datasetList, datasetPaths_4, True, saveExt="_4FitParams")
        DoPlots( [hL1Tks_NStubs_Vs_Eta      ], datasetList, datasetPaths_4, True, saveExt="_4FitParams")

        DoPlots( [hL1Tks_ChiSquared_Vs_POCAz], datasetList, datasetPaths_5, True, saveExt="_5FitParams")
        DoPlots( [hL1Tks_ChiSquared_Vs_Pt   ], datasetList, datasetPaths_5, True, saveExt="_5FitParams")
        DoPlots( [hL1Tks_ChiSquared_Vs_Eta  ], datasetList, datasetPaths_5, True, saveExt="_5FitParams")
        DoPlots( [hL1Tks_NStubs_Vs_Eta      ], datasetList, datasetPaths_5, True, saveExt="_5FitParams")

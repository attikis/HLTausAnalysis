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
import tools.datasets as m_datasets
import tools.aux as m_aux
pi = 4*math.atan(1)

###############################################################
### Options here
###############################################################
bDoTPs                = False #no
bDoMatchedTPs         = False #no
bDoEff                = False #yes (talk)
bDoMatchedTks         = False #yes
bDoMatchedTks_Res     = True #yes (talk)

bDoResVsPt_Pt         = False #maybe
bDoResVsPt_PtRel      = False #no
bDoResVsPt_Eta        = False #no
bDoResVsPt_Phi        = False #maybe
bDoResVsPt_Z0         = False #no
bDoResVsPt_d0         = False #no

bDoResVsEta_Pt        = False #no
bDoResVsEta_Eta       = False #no
bDoResVsEta_Phi       = False #no
bDoResVsEta_Z0        = False #no
bDoResVsEta_d0        = False #no

bDo2D                 = False #yes

bDo2DResVsPt_Pt       = True #yes (talk)
bDo2DResVsPt_PtRel    = True #yes
bDo2DMeanResVsPt_Pt   = True #yes
bDo2DResVsPt_Z0       = True #yes (talk)
bDo2DResVsPt_d0       = True #yes (talk)
bDo2DResVsPt_Phi      = True #yes (talk)
bDo2DResVsPt_Eta      = True #yes (talk)

bDo2DResVsEta_Pt      = True #yes (talk)
bDo2DResVsEta_PtRel   = True #yes (talk)
bDo2DMeanResVsEta_Eta = False #no
bDo2DResVsEta_Z0      = True #yes (talk)
bDo2DResVsEta_d0      = True #yes (talk)
bDo2DResVsEta_Phi     = True #yes (talk)
bDo2DResVsEta_Eta     = True #yes (talk)

nFitParams    = "5"
bDoPixelTks   = True
#inputPath     = "Macros/Tracking/results/" + nFitParams + "FitParams/Tracking_Histograms_"
#inputPath     = "Macros/Tracking/results/" + nFitParams + "FitParams/0p5xWindowsSF/Tracking_Histograms_"
inputPath     = "Macros/Tracking/results/" + nFitParams + "FitParams/1xWindowsSF/Tracking_Histograms_"
#inputPath     = "Macros/Tracking/results/" + nFitParams + "FitParams/2xWindowsSF/Tracking_Histograms_"
datasetList   = ["PiMinus"] #["SingleMuPlus", "PiPlus", "SinglePositron"]

if bDoPixelTks:
    ext       = "_Pixel"
    savePath  = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_12June2015/figures/" + nFitParams + "FitParams/pixTks/" + datasetList[0] + "/"
else:
    ext       = ""
    savePath  = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_12June2015/figures/" + nFitParams + "FitParams/ttTks/" + datasetList[0] + "/"
savePath      = ""
saveFormats   = ["png"]

datasetPaths  = {}
datasetPaths["MinBias"]                   = inputPath + "MinBias" + ext + ".root"
datasetPaths["VBF"]                       = inputPath + "VBF" + ext + ".root"
datasetPaths["PiPlus"]                    = inputPath + "PiPlus" + ext + ".root"
datasetPaths["PiMinus"]                   = inputPath + "PiMinus" + ext + ".root"
datasetPaths["SingleTauGun1p"]            = inputPath + "SingleTauGun1p" + ext + ".root"
datasetPaths["DiTauGun3p"]                = inputPath + "DiTauGun3p" + ext + ".root"
datasetPaths["TTbar"]                     = inputPath + "TTBar" + ext + ".root"
datasetPaths["HPlus160"]                  = inputPath + "HPlus160" + ext + ".root"
datasetPaths["HPlus200"]                  = inputPath + "HPlus200" + ext + ".root"
datasetPaths["SingleElectron"]            = inputPath + "SingleElectron" + ext + ".root"
datasetPaths["SinglePositron"]            = inputPath + "SinglePositron" + ext + ".root"
datasetPaths["SingleMuPlus"  ]            = inputPath + "SingleMuPlus" + ext + ".root"
datasetPaths["SingleMuMinus" ]            = inputPath + "SingleMuMinus" + ext + ".root"
datasetPaths["SinglePhoton"  ]            = inputPath + "SinglePhoton" + ext + ".root"  
# PRIVATE PRODUCTION
datasetPaths["SingleMuon_NoPU"]           = inputPath + "SingleMuon_NoPU" + ext + ".root"
datasetPaths["SingleMuon_PU140"]          = inputPath + "SingleMuon_PU140" + ext + ".root"
datasetPaths["SingleMuPlus_Pt_2_10_NoPU"] = inputPath + "SingleMuPlus_Pt_2_10_NoPU" + ext + ".root"
datasetPaths["SingleMuMinus_Pt_2_10_NoPU"]= inputPath + "SingleMuMinus_Pt_2_10_NoPU" + ext + ".root"



############################################################### 
### Definitions=
############################################################### 
C = "|#eta| < 0.8"
I = "0.8 < |#eta| < 1.6"
F = "|#eta| #geq 1.6"
L = "p_{T} < 5 GeVc^{-1}"
M = "5 < p_{T} < 15 GeVc^{-1}"
H = "p_{T} #geq 15 GeVc^{-1}"
EtaLines  = [-1.6, -0.8, +0.8, +1.6]
#EtaRegions= []
EtaRegions= [[-2.5, -1.6, ROOT.kRed+1], [+2.5, +1.6, ROOT.kRed+1], [-1.6, -0.8, ROOT.kYellow-4], [+0.8, +1.6, ROOT.kYellow-4], [-0.8, +0.8, ROOT.kGreen+1] ]

############################################################### 
### Histogram Options
############################################################### 
pTmax = 50.0 #100.0

Pt = {
    "xLabel": "p_{T}"           , "xUnits": "GeVc^{-1}", "xMin": 0.00 , "xMax": pTmax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": 1E-04, "yMax": +1e0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtL = {
    "xLabel": "p_{T}"           , "xUnits": "GeVc^{-1}", "xMin": 0.00 , "xMax": 5.0  , "binWidthX": 0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": 1E-04, "yMax": +1e0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


Eta = {
    "xLabel": "#eta"           , "xUnits": ""     , "xMin": -2.6 , "xMax": +2.6 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""     , "yMin": 1E-04, "yMax": +1e0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


Phi = {
    "xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.2 , "xMax": +3.2 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1E-04, "yMax": +1e0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


z0 = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -25.0, "xMax": +25.0  , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1E-04, "yMax": +1e0   , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


d0 = {
    # "xLabel": "d_{0}"          , "xUnits": "#mum", "xMin": -0.02, "xMax": +0.02, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "xLabel": "d_{0}"          , "xUnits": "#mum", "xMin": -50, "xMax": +50 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1E-04, "yMax": +2.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


NStubs = {
    "xLabel": "Stubs"          , "xUnits": "", "xMin": +0.0, "xMax": 15.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-4, "yMax": +2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.25, "yLegMax": 0.35
}

ChiSq  = {
    "xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +150.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.1f" , "yUnits": ""    , "yMin": 1E-04, "yMax": +1e0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

RedChiSq = {
    "xLabel": "#chi^{2}_{N}"   , "xUnits": "", "xMin": +0.00, "xMax": +10.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.1f", "yUnits": "", "yMin": 1e-04, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtRes = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual", "xUnits": "GeVc^{-1}", "xMin": -5.0 , "xMax": +5.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.3f"          , "yUnits": ""     , "yMin": 1E-04, "yMax": +1e0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    #"ratioLabel": "Ratio", "ratio": True, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtResAlt = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual", "xUnits": "GeVc^{-1}", "xMin": -5.0 , "xMax": +5.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.2f"          , "yUnits": ""     , "yMin": 1E-04, "yMax": +1e0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    #"ratioLabel": "Ratio", "ratio": True, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtResRel = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual/ p_{T}", "xUnits": "", "xMin": -2.0 , "xMax": +2.0, "binWidthX": 0.05, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.3f"      , "yUnits": "", "yMin": 1E-04, "yMax": +1e0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


EtaRes = { #residual = (L1 - Sim)
    "xLabel": "#eta residual", "xUnits": "GeVc^{-1}", "xMin": -0.01 , "xMax": 0.01, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.4f"          , "yUnits": ""     , "yMin": 1E-04, "yMax": +1e0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PhiRes = { #residual = (L1 - Sim)
    "xLabel": "#phi residual"  , "xUnits": "rads", "xMin": -0.005, "xMax": +0.005, "binWidthX": 0.0002, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.4f", "yUnits": ""    , "yMin": 1E-04 , "yMax": +1e0  , "binWidthY": None  , "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


Z0Res = { #residual =  (L1 - Sim)
    "xLabel": "z_{0} residual", "xUnits": "cm", "xMin": -1.0 , "xMax": +1.0, "binWidthX": None, "xCutLines": [-0.2, 0, +0.2], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    #"xLabel": "z_{0} residual", "xUnits": "cm", "xMin": -1.0 , "xMax": +1.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.3f"          , "yUnits": ""  , "yMin": 1E-04, "yMax": +1e0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

d0Res = { #residual =  (L1 - Sim)
    "xLabel": "d_{0} residual", "xUnits": "cm", "xMin": -0.5 , "xMax": +0.5, "binWidthX": None, "xCutLines": [-0.2, 0, +0.2], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    #"xLabel": "z_{0} residual", "xUnits": "cm", "xMin": -1.0 , "xMax": +1.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.3f"          , "yUnits": ""  , "yMin": 1E-04, "yMax": +1e0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


ChiSqEta = {
    "xLabel": "#eta / %0.2f"    , "xUnits": "", "xMin": -2.6 , "xMax": +2.6, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "log #chi^{2} / %0.2f", "yUnits": "", "yMin": +0.0 , "yMax": +5.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

RedChiSqEta = {
    "xLabel": "#eta / %0.2f"        , "xUnits": ""    , "xMin": -2.6 , "xMax": +2.6, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "log #chi^{2}_{N} / %0.2f", "yUnits": ""    , "yMin": +0.0 , "yMax": +5.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


dZ0Eta = {
    "xLabel": "#eta / %0.2f"             , "xUnits": ""  , "xMin": -2.6, "xMax": +2.6, "binWidthX": None, "xCutLines": [], "xCutBoxes": EtaRegions,  "gridX": True, "logX": False, "logXRatio": False,
    #"yLabel": "|#Deltaz_{0}| (L1 - Sim) / % 0.2f", "yUnits": "cm", "yMin":  0.0, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "yLabel": "|z_{0} residual| / % 0.2f", "yUnits": "cm", "yMin":  0.0, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dd0Eta = {
    "xLabel": "#eta / %0.2f"             , "xUnits": ""  , "xMin": -2.6, "xMax": +2.6, "binWidthX": None, "xCutLines": [], "xCutBoxes": EtaRegions,  "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "|d_{0} residual| / % 0.2f", "yUnits": "cm", "yMin":  0.0, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dEtaEta = {
    "xLabel": "#eta / %0.2f"                     , "xUnits": "", "xMin": -2.6, "xMax": +2.6  , "binWidthX": None, "xCutLines": [] , "gridX": True, "logX": False, "logXRatio": False,
    #"yLabel": "|#Delta#eta| (L1 - Sim) / % 0.2f" , "yUnits": "", "yMin":  0.0, "yMax": +0.015, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "yLabel": "|#eta residual| / % 0.4f" , "yUnits": "", "yMin":  0.0, "yMax": +0.015, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dPhiEta = {
    "xLabel": "#eta / %0.2f"                     , "xUnits": ""    , "xMin": -2.6, "xMax": +2.6  , "binWidthX": None, "xCutLines": [] , "gridX": True, "logX": False, "logXRatio": False,
    #"yLabel": "|#Delta#phi| (L1 - Sim) / % 0.2f" , "yUnits": "rads", "yMin":  0.0, "yMax": +0.010, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "yLabel": "|#phi residual| / % 0.4f" , "yUnits": "rads", "yMin":  0.0, "yMax": +0.010, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


dPtRelEta = {
    "xLabel": "#eta / %0.2f"                            , "xUnits": "", "xMin": -2.6, "xMax": +2.6, "binWidthX": None, "xCutLines": [] , "gridX": True, "logX": False, "logXRatio": False,
    #"yLabel": "|#Deltap_{T}/p_{T}| (L1 - Sim) / % 0.2f" , "yUnits": "", "yMin":  0.0, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "yLabel": "|p_{T} residual/p_{T}| / % 0.2f" , "yUnits": "", "yMin":  0.0, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dPtEta = {
    "xLabel": "#eta / %0.2f"              , "xUnits": ""     , "xMin": -2.6, "xMax": +2.6, "binWidthX": None, "xCutLines": [] , "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "|p_{T} residual| / % 0.2f" , "yUnits": "GeVc^{-1}", "yMin":  0.0, "yMax": +5.0, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, "logYRatio": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


PtPtRes = {
    "xLabel": "p_{T}"                                , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": pTmax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "p_{T} resolution (GeVc^{-1}) / % 0.2f", "yUnits": ""         , "yMin": 0.0 , "yMax":   2.5, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PtEtaRes = {
    "xLabel": "p_{T}"                   , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": pTmax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "#eta resolution / % 0.2f", "yUnits": ""     , "yMin": 0.0 , "yMax":  +0.0055, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PtPhiRes = {
    "xLabel": "p_{T}"                          , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": pTmax  , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#phi resolution (rads) / % 0.2f", "yUnits": ""         , "yMin": 0.00, "yMax":   0.0035, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtZ0Res = {
    "xLabel": "p_{T}"                    , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": pTmax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "z_{0} resolution / % 0.2f", "yUnits": ""         , "yMin": 0.0 , "yMax":   0.5, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

Ptd0Res = {
    "xLabel": "p_{T}"                    , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": pTmax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "d_{0} resolution / % 0.2f", "yUnits": ""     , "yMin": 0.0 , "yMax":   0.12, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


MeanPtPtRes = {
    "xLabel": "p_{T}"                         , "xUnits": "GeVc^{-1}", "xMin":  0.0, "xMax": pTmax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Mean(p_{T} residual)  / % 0.0f", "yUnits": ""         , "yMin": -0.6, "yMax":  +1.6, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "L", "legOptions": "LP", "logYRatio": False,  
    #"ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, #don't know why
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtPtResRel = {
    "xLabel": "p_{T}"                            , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": pTmax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "p_{T} resolution / p_{T} / % 0.2f", "yUnits": ""         , "yMin": 0.00, "yMax": 0.20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaPtRes = {
    "xLabel": "#eta"                     , "xUnits": "", "xMin": -2.6 , "xMax": +2.6 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "p_{T} resolution (GeVc^{-1}) / % 0.2f", "yUnits": "", "yMin": 0.0 , "yMax":  None , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaEtaRes = {
    "xLabel": "#eta"                    , "xUnits": ""     , "xMin": -2.6 , "xMax": +2.6 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "#eta resolution / % 0.2f", "yUnits": ""     , "yMin": 0.0 , "yMax":  +0.005, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaPhiRes = {
    "xLabel": "#eta"                           , "xUnits": "", "xMin": -2.6 , "xMax": +2.6 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "#phi resolution (rads) / % 0.2f", "yUnits": "", "yMin": 0.0 , "yMax":  0.004, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaZ0Res = {
    "xLabel": "#eta"                          , "xUnits": "", "xMin": -2.6 , "xMax": +2.6, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "z_{0} resolution (cm) / % 0.2f", "yUnits": "", "yMin": 0.0 , "yMax":  +0.7, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

Etad0Res = {
    "xLabel": "#eta"                          , "xUnits": "", "xMin": -2.6 , "xMax": +2.6, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "d_{0} resolution (cm) / % 0.2f", "yUnits": "", "yMin": 0.0 , "yMax":  +0.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

MeanEtaEtaRes = {
    "xLabel": "#eta"                          , "xUnits": "" , "xMin": -2.600, "xMax": +2.600, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Mean(#eta residual)  / % 0.2f" , "yUnits": "" , "yMin": -0.001, "yMax": +0.001, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaPtResRel = {
    "xLabel": "#eta"                             , "xUnits": "" , "xMin": -2.6 , "xMax": +2.6, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "p_{T} resolution / p_{T} / % 0.2f", "yUnits": "" , "yMin": 0.0 , "yMax":  +0.1, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


EffPt = {
    "xLabel": "p_{T}"              , "xUnits": "GeVc^{-1}", "xMin": 0.0 , "xMax": +50.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Efficiency / %0.1f" , "yUnits": ""     , "yMin": 0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}

EffPtL = {
    "xLabel": "p_{T}"              , "xUnits": "GeVc^{-1}", "xMin": 0.0, "xMax": +5.0  , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Efficiency / %0.1f" , "yUnits": ""     , "yMin": 0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", "logYRatio": False, 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


EffEta = {
    #"xLabel": "#eta"               , "xUnits": "" , "xMin": -2.6, "xMax": +2.6  , "binWidthX": None, "xCutLines": [], "xCutBoxes": EtaRegions, "gridX": True, "logX": False, "logXRatio": False, 
    "xLabel": "#eta"               , "xUnits": "" , "xMin": -2.6, "xMax": +2.6  , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Efficiency / %0.2f" , "yUnits": "" , "yMin": 0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


EffPhi = {
    "xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.2 , "xMax": +3.2 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Efficiency / %0.2f" , "yUnits": ""  , "yMin": 0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01 , "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


Effz0 = {
    "xLabel": "z_{0}"              , "xUnits": "cm", "xMin": -25.0 , "xMax": +25.0 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Efficiency / %0.2f" , "yUnits": ""  , "yMin":  +0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


Effd0 = {
    "xLabel": "d_{0}"             , "xUnits": "#mum", "xMin": -50  , "xMax": +50   , "binWidthX": None , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logYRatio": False,
    "yLabel": "Efficiency / %0.2f", "yUnits": ""    , "yMin": +0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.20, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


############################################################### 
### Group Histogram (of same binning)
############################################################### 
hFolder = ""

### TPs
tp_pt     = m_histos.TH1orTH2( hFolder, "tp_pt"    , "TP"      , None, **Pt)
tp_pt_L   = m_histos.TH1orTH2( hFolder, "tp_pt_L"  , "TP, " + L, None, **PtL)
tp_pt_C   = m_histos.TH1orTH2( hFolder, "tp_pt_C"  , "TP, " + C, None, **Pt)
tp_pt_I   = m_histos.TH1orTH2( hFolder, "tp_pt_I"  , "TP, " + I, None, **Pt)
tp_pt_F   = m_histos.TH1orTH2( hFolder, "tp_pt_F"  , "TP, " + F, None, **Pt)
tp_eta    = m_histos.TH1orTH2( hFolder, "tp_eta"   , "TP"      , None, **Eta)
tp_eta_L  = m_histos.TH1orTH2( hFolder, "tp_eta_L" , "TP, " + L, None, **Eta)
tp_eta_M  = m_histos.TH1orTH2( hFolder, "tp_eta_M" , "TP, " + M, None, **Eta)
tp_eta_H  = m_histos.TH1orTH2( hFolder, "tp_eta_H" , "TP, " + H, None, **Eta)
tp_phi    = m_histos.TH1orTH2( hFolder, "tp_phi"   , "TP"      , None, **Phi)
tp_z0     = m_histos.TH1orTH2( hFolder, "tp_z0"    , "TP"      , None, **z0)
tp_d0     = m_histos.TH1orTH2( hFolder, "tp_d0"    , "TP"      , None, **d0)

### TPs (TTTrack-Matched)
match_tp_pt    = m_histos.TH1orTH2( hFolder, "match_tp_pt"   , "Matched TP"      , None, **Pt )
match_tp_pt_L  = m_histos.TH1orTH2( hFolder, "match_tp_pt_L" , "Matched TP, " + L, None, **PtL)
match_tp_pt_C  = m_histos.TH1orTH2( hFolder, "match_tp_pt_C" , "Matched TP, " + C, None, **Pt )
match_tp_pt_I  = m_histos.TH1orTH2( hFolder, "match_tp_pt_I" , "Matched TP, " + I, None, **Pt )
match_tp_pt_F  = m_histos.TH1orTH2( hFolder, "match_tp_pt_F" , "Matched TP, " + F, None, **Pt )
match_tp_eta   = m_histos.TH1orTH2( hFolder, "match_tp_eta"  , "Matched TP"      , None, **Eta)
match_tp_eta_L = m_histos.TH1orTH2( hFolder, "match_tp_eta_L", "Matched TP, " + L, None, **Eta)
match_tp_eta_M = m_histos.TH1orTH2( hFolder, "match_tp_eta_M", "Matched TP, " + M, None, **Eta)
match_tp_eta_H = m_histos.TH1orTH2( hFolder, "match_tp_eta_H", "Matched TP, " + H, None, **Eta)
match_tp_phi   = m_histos.TH1orTH2( hFolder, "match_tp_phi"  , "Matched TP"      , None, **Phi)
match_tp_z0    = m_histos.TH1orTH2( hFolder, "match_tp_z0"   , "Matched TP"      , None, **z0 ) 
match_tp_d0    = m_histos.TH1orTH2( hFolder, "match_tp_d0"   , "Matched TP"      , None, **d0 )

### Efficiencies
eff_pt    = m_histos.TH1orTH2( hFolder, "eff_pt"   , "overall", None, **EffPt )
eff_pt_L  = m_histos.TH1orTH2( hFolder, "eff_pt_L" , "" + L   , None, **EffPtL)
eff_pt_C  = m_histos.TH1orTH2( hFolder, "eff_pt_C" , "" + C   , None, **EffPt)
eff_pt_I  = m_histos.TH1orTH2( hFolder, "eff_pt_I" , "" + I   , None, **EffPt)
eff_pt_F  = m_histos.TH1orTH2( hFolder, "eff_pt_F" , "" + F   , None, **EffPt)
eff_eta   = m_histos.TH1orTH2( hFolder, "eff_eta"  , "overall", None, **EffEta)
eff_eta_L = m_histos.TH1orTH2( hFolder, "eff_eta_L", "" + L   , None, **EffEta)
eff_eta_M = m_histos.TH1orTH2( hFolder, "eff_eta_M", "" + M   , None, **EffEta)
eff_eta_H = m_histos.TH1orTH2( hFolder, "eff_eta_H", "" + H   , None, **EffEta)
eff_phi   = m_histos.TH1orTH2( hFolder, "eff_phi"  , "overall", None, **EffPhi)
eff_z0    = m_histos.TH1orTH2( hFolder, "eff_z0"   , "overall", None, **Effz0 )
eff_d0    = m_histos.TH1orTH2( hFolder, "eff_d0"   , "overall", None, **Effd0 )

### TTTrack (TP-Matched)
label = "" #"Tk, "
match_trk_nstub   = m_histos.TH1orTH2( hFolder, "match_trk_nstub"  , "Tk"     , None, **NStubs)
match_trk_nstub_C = m_histos.TH1orTH2( hFolder, "match_trk_nstub_C", label + C, None, **NStubs)
match_trk_nstub_I = m_histos.TH1orTH2( hFolder, "match_trk_nstub_I", label + I, None, **NStubs)
match_trk_nstub_F = m_histos.TH1orTH2( hFolder, "match_trk_nstub_F", label + F, None, **NStubs)

### TTTrack (TP-Matched): ChiSq histograms (last bin is an overflow bin)
match_trk_chi2     = m_histos.TH1orTH2( hFolder, "match_trk_chi2"    , "Tk"                  , None, **ChiSq)
match_trk_chi2_L   = m_histos.TH1orTH2( hFolder, "match_trk_chi2_L"  , label + L             , None, **ChiSq)
match_trk_chi2_M   = m_histos.TH1orTH2( hFolder, "match_trk_chi2_M"  , label + M             , None, **ChiSq)
match_trk_chi2_H   = m_histos.TH1orTH2( hFolder, "match_trk_chi2_H"  , label + H             , None, **ChiSq)
match_trk_chi2_C_L = m_histos.TH1orTH2( hFolder, "match_trk_chi2_C_L", label + C + " , "  + L, None, **ChiSq)
match_trk_chi2_I_L = m_histos.TH1orTH2( hFolder, "match_trk_chi2_I_L", label + I + " , "  + L, None, **ChiSq)
match_trk_chi2_F_L = m_histos.TH1orTH2( hFolder, "match_trk_chi2_F_L", label + F + " , "  + L, None, **ChiSq)
match_trk_chi2_C_M = m_histos.TH1orTH2( hFolder, "match_trk_chi2_C_M", label + C + " , "  + M, None, **ChiSq)
match_trk_chi2_I_M = m_histos.TH1orTH2( hFolder, "match_trk_chi2_I_M", label + I + " , "  + M, None, **ChiSq)
match_trk_chi2_F_M = m_histos.TH1orTH2( hFolder, "match_trk_chi2_F_M", label + F + " , "  + M, None, **ChiSq)
match_trk_chi2_C_H = m_histos.TH1orTH2( hFolder, "match_trk_chi2_C_H", label + C + " , "  + H, None, **ChiSq)
match_trk_chi2_I_H = m_histos.TH1orTH2( hFolder, "match_trk_chi2_I_H", label + I + " , "  + H, None, **ChiSq)
match_trk_chi2_F_H = m_histos.TH1orTH2( hFolder, "match_trk_chi2_F_H", label + F + " , "  + H, None, **ChiSq)

### TTTrack (TP-Matched): RedChiSq histograms (lastbin is an overflow bin)
match_trk_chi2_dof     = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof"    , "Tk"                  , None, **RedChiSq)
match_trk_chi2_dof_L   = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_L"  , label + L             , None, **RedChiSq)
match_trk_chi2_dof_M   = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_M"  , label + M             , None, **RedChiSq)
match_trk_chi2_dof_H   = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_H"  , label + H             , None, **RedChiSq)
match_trk_chi2_dof_C_L = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_C_L", label + C + " , "  + L, None, **RedChiSq)
match_trk_chi2_dof_I_L = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_I_L", label + I + " , "  + L, None, **RedChiSq)
match_trk_chi2_dof_F_L = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_F_L", label + F + " , "  + L, None, **RedChiSq)
match_trk_chi2_dof_C_M = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_C_M", label + C + " , "  + M, None, **RedChiSq)
match_trk_chi2_dof_I_M = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_I_M", label + I + " , "  + M, None, **RedChiSq)
match_trk_chi2_dof_F_M = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_F_M", label + F + " , "  + M, None, **RedChiSq)
match_trk_chi2_dof_C_H = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_C_H", label + C + " , "  + H, None, **RedChiSq)
match_trk_chi2_dof_I_H = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_I_H", label + I + " , "  + H, None, **RedChiSq)
match_trk_chi2_dof_F_H = m_histos.TH1orTH2( hFolder, "match_trk_chi2_dof_F_H", label + F + " , "  + H, None, **RedChiSq)

### TTTrack (TP-Matched): Resolution histograms
res_pt      = m_histos.TH1orTH2( hFolder, "res_pt"     , "Tk"     , None, **PtRes   )
res_pt_C    = m_histos.TH1orTH2( hFolder, "res_pt_C"   , label + C, None, **PtResAlt)
res_pt_I    = m_histos.TH1orTH2( hFolder, "res_pt_I"   , label + I, None, **PtResAlt)
res_pt_F    = m_histos.TH1orTH2( hFolder, "res_pt_F"   , label + F, None, **PtResAlt)

res_ptRel   = m_histos.TH1orTH2( hFolder, "res_ptRel"  , label    , None, **PtResRel)
res_ptRel_C = m_histos.TH1orTH2( hFolder, "res_ptRel_C", label + C, None, **PtResRel)
res_ptRel_I = m_histos.TH1orTH2( hFolder, "res_ptRel_I", label + I, None, **PtResRel)
res_ptRel_F = m_histos.TH1orTH2( hFolder, "res_ptRel_F", label + F, None, **PtResRel)

res_eta     = m_histos.TH1orTH2( hFolder, "res_eta"    , label    , None, **EtaRes  )
res_eta_C   = m_histos.TH1orTH2( hFolder, "res_eta_C"  , label + C, None, **EtaRes  )
res_eta_I   = m_histos.TH1orTH2( hFolder, "res_eta_I"  , label + I, None, **EtaRes  )
res_eta_F   = m_histos.TH1orTH2( hFolder, "res_eta_F"  , label + F, None, **EtaRes  )

res_phi     = m_histos.TH1orTH2( hFolder, "res_phi"    , label    , None, **PhiRes  )
res_phi_C   = m_histos.TH1orTH2( hFolder, "res_phi_C"  , label + C, None, **PhiRes  )
res_phi_I   = m_histos.TH1orTH2( hFolder, "res_phi_I"  , label + I, None, **PhiRes  )
res_phi_F   = m_histos.TH1orTH2( hFolder, "res_phi_F"  , label + F, None, **PhiRes  )

res_z0      = m_histos.TH1orTH2( hFolder, "res_z0"     , label    , None, **Z0Res   )
res_z0_C    = m_histos.TH1orTH2( hFolder, "res_z0_C"   , label + C, None, **Z0Res   )
res_z0_I    = m_histos.TH1orTH2( hFolder, "res_z0_I"   , label + I, None, **Z0Res   )
res_z0_F    = m_histos.TH1orTH2( hFolder, "res_z0_F"   , label + F, None, **Z0Res   )

res_d0      = m_histos.TH1orTH2( hFolder, "res_d0"     , label    , None, **d0Res   )
res_d0_C    = m_histos.TH1orTH2( hFolder, "res_d0_C"   , label + C, None, **d0Res   )
res_d0_I    = m_histos.TH1orTH2( hFolder, "res_d0_I"   , label + I, None, **d0Res   )
res_d0_F    = m_histos.TH1orTH2( hFolder, "res_d0_F"   , label + F, None, **d0Res   )

### Resolution vs. pt histograms  
ptRange = ["0-5","5-10", "10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50","50-55", "55-60","60-65","65-70","70-75","75-80","80-85","85-90","90-95","95-100"]
h_resVsPt_pt      = []
h_resVsPt_pt_C    = []
h_resVsPt_pt_I    = []
h_resVsPt_pt_F    = []

h_resVsPt_ptRel   = []
h_resVsPt_ptRel_C = []
h_resVsPt_ptRel_I = []
h_resVsPt_ptRel_F = []

h_resVsPt_eta     = []
h_resVsPt_eta_C   = []
h_resVsPt_eta_I   = []
h_resVsPt_eta_F   = []

h_resVsPt_phi     = []
h_resVsPt_phi_C   = []
h_resVsPt_phi_I   = []
h_resVsPt_phi_F   = []

h_resVsPt_z0      = []
h_resVsPt_z0_C    = []
h_resVsPt_z0_I    = []
h_resVsPt_z0_F    = []

h_resVsPt_d0      = []
h_resVsPt_d0_C    = []
h_resVsPt_d0_I    = []
h_resVsPt_d0_F    = []

for i in range(0, len(ptRange)):
    resVsPt_pt      = m_histos.TH1orTH2( hFolder, "resVsPt_pt_"      + ptRange[i], ptRange[i] + " GeV"      , None, **PtRes)
    resVsPt_pt_C    = m_histos.TH1orTH2( hFolder, "resVsPt_pt_C_"    + ptRange[i], ptRange[i] + " GeV, " + C, None, **PtRes)
    resVsPt_pt_I    = m_histos.TH1orTH2( hFolder, "resVsPt_pt_I_"    + ptRange[i], ptRange[i] + " GeV, " + I, None, **PtRes)
    resVsPt_pt_F    = m_histos.TH1orTH2( hFolder, "resVsPt_pt_F_"    + ptRange[i], ptRange[i] + " GeV, " + F, None, **PtRes)

    resVsPt_ptRel   = m_histos.TH1orTH2( hFolder, "resVsPt_ptRel_"   + ptRange[i], ptRange[i] + " GeV"      , None, **PtResRel)
    resVsPt_ptRel_C = m_histos.TH1orTH2( hFolder, "resVsPt_ptRel_C_" + ptRange[i], ptRange[i] + " GeV, " + C, None, **PtResRel)
    resVsPt_ptRel_I = m_histos.TH1orTH2( hFolder, "resVsPt_ptRel_I_" + ptRange[i], ptRange[i] + " GeV, " + I, None, **PtResRel)
    resVsPt_ptRel_F = m_histos.TH1orTH2( hFolder, "resVsPt_ptRel_F_" + ptRange[i], ptRange[i] + " GeV, " + F, None, **PtResRel)

    resVsPt_eta     = m_histos.TH1orTH2( hFolder, "resVsPt_eta_"     + ptRange[i], ptRange[i] + " GeV"      , None, **EtaRes)
    resVsPt_eta_C   = m_histos.TH1orTH2( hFolder, "resVsPt_eta_C_"   + ptRange[i], ptRange[i] + " GeV, " + C, None, **EtaRes)
    resVsPt_eta_I   = m_histos.TH1orTH2( hFolder, "resVsPt_eta_I_"   + ptRange[i], ptRange[i] + " GeV, " + I, None, **EtaRes)
    resVsPt_eta_F   = m_histos.TH1orTH2( hFolder, "resVsPt_eta_F_"   + ptRange[i], ptRange[i] + " GeV, " + F, None, **EtaRes)

    resVsPt_phi     = m_histos.TH1orTH2( hFolder, "resVsPt_phi_"     + ptRange[i], ptRange[i] + " GeV"      , None, **PhiRes)
    resVsPt_phi_C   = m_histos.TH1orTH2( hFolder, "resVsPt_phi_C_"   + ptRange[i], ptRange[i] + " GeV, " + C, None, **PhiRes)
    resVsPt_phi_I   = m_histos.TH1orTH2( hFolder, "resVsPt_phi_I_"   + ptRange[i], ptRange[i] + " GeV, " + I, None, **PhiRes)
    resVsPt_phi_F   = m_histos.TH1orTH2( hFolder, "resVsPt_phi_F_"   + ptRange[i], ptRange[i] + " GeV, " + F, None, **PhiRes)

    resVsPt_z0      = m_histos.TH1orTH2( hFolder, "resVsPt_z0_"      + ptRange[i], ptRange[i] + " GeV"      , None, **Z0Res)
    resVsPt_z0_C    = m_histos.TH1orTH2( hFolder, "resVsPt_z0_C_"    + ptRange[i], ptRange[i] + " GeV, " + C, None, **Z0Res)
    resVsPt_z0_I    = m_histos.TH1orTH2( hFolder, "resVsPt_z0_I_"    + ptRange[i], ptRange[i] + " GeV, " + I, None, **Z0Res)
    resVsPt_z0_F    = m_histos.TH1orTH2( hFolder, "resVsPt_z0_F_"    + ptRange[i], ptRange[i] + " GeV, " + F, None, **Z0Res)

    resVsPt_d0      = m_histos.TH1orTH2( hFolder, "resVsPt_d0_"      + ptRange[i], ptRange[i] + " GeV"      , None, **d0Res)
    resVsPt_d0_C    = m_histos.TH1orTH2( hFolder, "resVsPt_d0_C_"    + ptRange[i], ptRange[i] + " GeV, " + C, None, **d0Res)
    resVsPt_d0_I    = m_histos.TH1orTH2( hFolder, "resVsPt_d0_I_"    + ptRange[i], ptRange[i] + " GeV, " + I, None, **d0Res)
    resVsPt_d0_F    = m_histos.TH1orTH2( hFolder, "resVsPt_d0_F_"    + ptRange[i], ptRange[i] + " GeV, " + F, None, **d0Res)


    ### Append to lists
    h_resVsPt_pt.append(resVsPt_pt)
    h_resVsPt_pt_C.append(resVsPt_pt_C)
    h_resVsPt_pt_I.append(resVsPt_pt_I)
    h_resVsPt_pt_F.append(resVsPt_pt_F)

    h_resVsPt_ptRel.append(resVsPt_ptRel)
    h_resVsPt_ptRel_C.append(resVsPt_ptRel_C)
    h_resVsPt_ptRel_I.append(resVsPt_ptRel_I)
    h_resVsPt_ptRel_F.append(resVsPt_ptRel_F)

    h_resVsPt_eta.append(resVsPt_eta)
    h_resVsPt_eta_C.append(resVsPt_eta_C)
    h_resVsPt_eta_I.append(resVsPt_eta_I)
    h_resVsPt_eta_F.append(resVsPt_eta_F)

    h_resVsPt_phi.append(resVsPt_phi)
    h_resVsPt_phi_C.append(resVsPt_phi_C)
    h_resVsPt_phi_I.append(resVsPt_phi_I)
    h_resVsPt_phi_F.append(resVsPt_phi_F)

    h_resVsPt_z0.append(resVsPt_z0)
    h_resVsPt_z0_C.append(resVsPt_z0_C)
    h_resVsPt_z0_I.append(resVsPt_z0_I)
    h_resVsPt_z0_F.append(resVsPt_z0_F)

    h_resVsPt_d0.append(resVsPt_d0)
    h_resVsPt_d0_C.append(resVsPt_d0_C)
    h_resVsPt_d0_I.append(resVsPt_d0_I)
    h_resVsPt_d0_F.append(resVsPt_d0_F)


### resolution vs. eta histograms
etaRange     = ["0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0", "1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2.0", "2.1","2.2","2.3","2.4","2.5"]
etaRangePlus = ["0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0", "1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2.0", "2.1","2.2","2.3","2.4","2.5","2.6"]
h_resVsEta_pt      = []
h_resVsEta_ptRel   = []
h_resVsEta_ptRel_L = []
h_resVsEta_ptRel_M = []
h_resVsEta_ptRel_H = []

h_resVsEta_eta     = []
h_resVsEta_eta_L   = []
h_resVsEta_eta_M   = []
h_resVsEta_eta_H   = []

h_resVsEta_phi     = []
h_resVsEta_phi_L   = []
h_resVsEta_phi_M   = []
h_resVsEta_phi_H   = []

h_resVsEta_z0      = []
h_resVsEta_z0_L    = []
h_resVsEta_z0_M    = []
h_resVsEta_z0_H    = []

h_resVsEta_d0      = []
h_resVsEta_d0_L    = []
h_resVsEta_d0_M    = []
h_resVsEta_d0_H    = []

for i in range(0, len(etaRange)):    
    resVsEta_pt      = m_histos.TH1orTH2( hFolder, "resVsEta_pt_"      + etaRange[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **PtRes)
    resVsEta_ptRel   = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel_"   + etaRange[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **PtResRel)
    resVsEta_ptRel_L = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel_L_" + etaRange[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtResRel)
    resVsEta_ptRel_M = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel_M_" + etaRange[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtResRel)
    resVsEta_ptRel_H = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel_H_" + etaRange[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtResRel)

    resVsEta_eta     = m_histos.TH1orTH2( hFolder, "resVsEta_eta_"     + etaRange[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **EtaRes)
    resVsEta_eta_L   = m_histos.TH1orTH2( hFolder, "resVsEta_eta_L_"   + etaRange[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **EtaRes)
    resVsEta_eta_M   = m_histos.TH1orTH2( hFolder, "resVsEta_eta_M_"   + etaRange[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **EtaRes)
    resVsEta_eta_H   = m_histos.TH1orTH2( hFolder, "resVsEta_eta_H_"   + etaRange[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **EtaRes)

    resVsEta_phi     = m_histos.TH1orTH2( hFolder, "resVsEta_phi_"     + etaRange[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **PhiRes)
    resVsEta_phi_L   = m_histos.TH1orTH2( hFolder, "resVsEta_phi_L_"   + etaRange[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PhiRes)
    resVsEta_phi_M   = m_histos.TH1orTH2( hFolder, "resVsEta_phi_M_"   + etaRange[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PhiRes)
    resVsEta_phi_H   = m_histos.TH1orTH2( hFolder, "resVsEta_phi_H_"   + etaRange[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PhiRes)

    resVsEta_z0      = m_histos.TH1orTH2( hFolder, "resVsEta_z0_"      + etaRange[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **Z0Res)
    resVsEta_z0_L    = m_histos.TH1orTH2( hFolder, "resVsEta_z0_L_"    + etaRange[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **Z0Res)
    resVsEta_z0_M    = m_histos.TH1orTH2( hFolder, "resVsEta_z0_M_"    + etaRange[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **Z0Res)
    resVsEta_z0_H    = m_histos.TH1orTH2( hFolder, "resVsEta_z0_H_"    + etaRange[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **Z0Res)

    resVsEta_d0      = m_histos.TH1orTH2( hFolder, "resVsEta_d0_"      + etaRange[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **d0Res)
    resVsEta_d0_L    = m_histos.TH1orTH2( hFolder, "resVsEta_d0_L_"    + etaRange[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **d0Res)
    resVsEta_d0_M    = m_histos.TH1orTH2( hFolder, "resVsEta_d0_M_"    + etaRange[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **d0Res)
    resVsEta_d0_H    = m_histos.TH1orTH2( hFolder, "resVsEta_d0_H_"    + etaRange[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **d0Res)

    ### Append to lists
    h_resVsEta_pt.append(resVsEta_pt)
    h_resVsEta_ptRel.append(resVsEta_ptRel)
    h_resVsEta_ptRel_L.append(resVsEta_ptRel_L)
    h_resVsEta_ptRel_M.append(resVsEta_ptRel_M)
    h_resVsEta_ptRel_H.append(resVsEta_ptRel_H)

    h_resVsEta_eta.append(resVsEta_eta)
    h_resVsEta_eta_L.append(resVsEta_eta_L)
    h_resVsEta_eta_M.append(resVsEta_eta_M)
    h_resVsEta_eta_H.append(resVsEta_eta_H)

    h_resVsEta_phi.append(resVsEta_phi)
    h_resVsEta_phi_L.append(resVsEta_phi_L)
    h_resVsEta_phi_M.append(resVsEta_phi_M)
    h_resVsEta_phi_H.append(resVsEta_phi_H)

    h_resVsEta_z0.append(resVsEta_z0)
    h_resVsEta_z0_L.append(resVsEta_z0_L)
    h_resVsEta_z0_M.append(resVsEta_z0_M)
    h_resVsEta_z0_H.append(resVsEta_z0_H)

    h_resVsEta_d0.append(resVsEta_d0)
    h_resVsEta_d0_L.append(resVsEta_d0_L)
    h_resVsEta_d0_M.append(resVsEta_d0_M)
    h_resVsEta_d0_H.append(resVsEta_d0_H)


### 2D histograms
logChiSq_eta     = m_histos.TH1orTH2( hFolder, "2d_logchi2_eta"    , "", None, **ChiSqEta)
logChiSq_dof_eta = m_histos.TH1orTH2( hFolder, "2d_logchi2_dof_eta", "", None, **RedChiSqEta)
dz0_eta          = m_histos.TH1orTH2( hFolder, "2d_dz0_eta"        , "", None, **dZ0Eta)
dd0_eta          = m_histos.TH1orTH2( hFolder, "2d_dd0_eta"        , "", None, **dd0Eta)
deta_eta         = m_histos.TH1orTH2( hFolder, "2d_deta_eta"       , "", None, **dEtaEta)
dphi_eta         = m_histos.TH1orTH2( hFolder, "2d_dphi_eta"       , "", None, **dPhiEta)
dpt_eta          = m_histos.TH1orTH2( hFolder, "2d_dpt_eta"        , "", None, **dPtEta)
dptRel_eta       = m_histos.TH1orTH2( hFolder, "2d_dptRel_eta"     , "", None, **dPtRelEta)


### resolution vs. pT histograms. NB: Filled with RMS and RMSError
resVsPt_pt   = m_histos.TH1orTH2( hFolder, "resVsPt_pt"  , "", None, **PtPtRes)
resVsPt_pt_C = m_histos.TH1orTH2( hFolder, "resVsPt_pt_C", C , None, **PtPtRes)
resVsPt_pt_I = m_histos.TH1orTH2( hFolder, "resVsPt_pt_I", I , None, **PtPtRes)
resVsPt_pt_F = m_histos.TH1orTH2( hFolder, "resVsPt_pt_F", F , None, **PtPtRes)

resVsPt_ptRel   = m_histos.TH1orTH2( hFolder,  "resVsPt_ptRel"  ,  "", None, **PtPtResRel)
resVsPt_ptRel_C = m_histos.TH1orTH2( hFolder,  "resVsPt_ptRel_C",  C , None, **PtPtResRel)
resVsPt_ptRel_I = m_histos.TH1orTH2( hFolder,  "resVsPt_ptRel_I",  I , None, **PtPtResRel)
resVsPt_ptRel_F = m_histos.TH1orTH2( hFolder,  "resVsPt_ptRel_F",  F , None, **PtPtResRel)

mresVsPt_pt    = m_histos.TH1orTH2( hFolder, "mresVsPt_pt"  ,  "", None, **MeanPtPtRes)
mresVsPt_pt_C  = m_histos.TH1orTH2( hFolder, "mresVsPt_pt_C",  C , None, **MeanPtPtRes)
mresVsPt_pt_I  = m_histos.TH1orTH2( hFolder, "mresVsPt_pt_I",  I , None, **MeanPtPtRes)
mresVsPt_pt_F  = m_histos.TH1orTH2( hFolder, "mresVsPt_pt_F",  F , None, **MeanPtPtRes)

resVsPt_eta   = m_histos.TH1orTH2( hFolder, "resVsPt_eta"  ,  "", None, **PtEtaRes)
resVsPt_eta_C = m_histos.TH1orTH2( hFolder, "resVsPt_eta_C",  C, None, **PtEtaRes)
resVsPt_eta_I = m_histos.TH1orTH2( hFolder, "resVsPt_eta_I",  I, None, **PtEtaRes)
resVsPt_eta_F = m_histos.TH1orTH2( hFolder, "resVsPt_eta_F",  F, None, **PtEtaRes)

resVsPt_phi   = m_histos.TH1orTH2( hFolder, "resVsPt_phi"  ,  "", None, **PtPhiRes)
resVsPt_phi_C = m_histos.TH1orTH2( hFolder, "resVsPt_phi_C",  C , None, **PtPhiRes)
resVsPt_phi_I = m_histos.TH1orTH2( hFolder, "resVsPt_phi_I",  I , None, **PtPhiRes)
resVsPt_phi_F = m_histos.TH1orTH2( hFolder, "resVsPt_phi_F",  F , None, **PtPhiRes)

resVsPt_z0    = m_histos.TH1orTH2( hFolder, "resVsPt_z0"  ,  "", None, **PtZ0Res)
resVsPt_z0_C  = m_histos.TH1orTH2( hFolder, "resVsPt_z0_C",  C , None, **PtZ0Res)
resVsPt_z0_I  = m_histos.TH1orTH2( hFolder, "resVsPt_z0_I",  I , None, **PtZ0Res)
resVsPt_z0_F  = m_histos.TH1orTH2( hFolder, "resVsPt_z0_F",  F , None, **PtZ0Res)

resVsPt_d0    = m_histos.TH1orTH2( hFolder, "resVsPt_d0"  ,  "", None, **Ptd0Res)
resVsPt_d0_C  = m_histos.TH1orTH2( hFolder, "resVsPt_d0_C",  C , None, **Ptd0Res)
resVsPt_d0_I  = m_histos.TH1orTH2( hFolder, "resVsPt_d0_I",  I , None, **Ptd0Res)
resVsPt_d0_F  = m_histos.TH1orTH2( hFolder, "resVsPt_d0_F",  F , None, **Ptd0Res)


### resolution vs. eta histograms
mresVsEta_eta   = m_histos.TH1orTH2( hFolder, "mresVsEta_eta"  , "", None, **MeanEtaEtaRes)
mresVsEta_eta_L = m_histos.TH1orTH2( hFolder, "mresVsEta_eta_L", L , None, **MeanEtaEtaRes)
mresVsEta_eta_M = m_histos.TH1orTH2( hFolder, "mresVsEta_eta_M", M , None, **MeanEtaEtaRes)
mresVsEta_eta_H = m_histos.TH1orTH2( hFolder, "mresVsEta_eta_H", H , None, **MeanEtaEtaRes)

resVsEta_pt   = m_histos.TH1orTH2( hFolder, "resVsEta_pt"  , "", None, **EtaPtRes)
resVsEta_pt_L = m_histos.TH1orTH2( hFolder, "resVsEta_pt_L", L , None, **EtaPtRes)
resVsEta_pt_M = m_histos.TH1orTH2( hFolder, "resVsEta_pt_M", M , None, **EtaPtRes)
resVsEta_pt_H = m_histos.TH1orTH2( hFolder, "resVsEta_pt_H", H , None, **EtaPtRes)

resVsEta_ptRel   = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel"  , "", None, **EtaPtResRel)
resVsEta_ptRel_L = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel_L", L , None, **EtaPtResRel)
resVsEta_ptRel_M = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel_M", M , None, **EtaPtResRel)
resVsEta_ptRel_H = m_histos.TH1orTH2( hFolder, "resVsEta_ptRel_H", H , None, **EtaPtResRel)

resVsEta_eta   = m_histos.TH1orTH2( hFolder, "resVsEta_eta"  , "", None, **EtaEtaRes)
resVsEta_eta_L = m_histos.TH1orTH2( hFolder, "resVsEta_eta_L", L , None, **EtaEtaRes)
resVsEta_eta_M = m_histos.TH1orTH2( hFolder, "resVsEta_eta_M", M , None, **EtaEtaRes)
resVsEta_eta_H = m_histos.TH1orTH2( hFolder, "resVsEta_eta_H", H , None, **EtaEtaRes)

resVsEta_phi   = m_histos.TH1orTH2( hFolder, "resVsEta_phi"  , "", None, **EtaPhiRes)
resVsEta_phi_L = m_histos.TH1orTH2( hFolder, "resVsEta_phi_L", L , None, **EtaPhiRes)
resVsEta_phi_M = m_histos.TH1orTH2( hFolder, "resVsEta_phi_M", M , None, **EtaPhiRes)
resVsEta_phi_H = m_histos.TH1orTH2( hFolder, "resVsEta_phi_H", H , None, **EtaPhiRes)

resVsEta_z0   = m_histos.TH1orTH2( hFolder, "resVsEta_z0"  , "", None, **EtaZ0Res)
resVsEta_z0_L = m_histos.TH1orTH2( hFolder, "resVsEta_z0_L", L , None, **EtaZ0Res)
resVsEta_z0_M = m_histos.TH1orTH2( hFolder, "resVsEta_z0_M", M , None, **EtaZ0Res)
resVsEta_z0_H = m_histos.TH1orTH2( hFolder, "resVsEta_z0_H", H , None, **EtaZ0Res)

resVsEta_d0   = m_histos.TH1orTH2( hFolder, "resVsEta_d0"  , "", None, **Etad0Res)
resVsEta_d0_L = m_histos.TH1orTH2( hFolder, "resVsEta_d0_L", L , None, **Etad0Res)
resVsEta_d0_M = m_histos.TH1orTH2( hFolder, "resVsEta_d0_M", M , None, **Etad0Res)
resVsEta_d0_H = m_histos.TH1orTH2( hFolder, "resVsEta_d0_H", H , None, **Etad0Res)




###############################################################
### Main
###############################################################
def DoPlots(hList, datasetList, bColourPalette=False, saveExt=""):

    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.SetBoolUseDatasetAsLegEntry(not bColourPalette)
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
    
    if(bDoTPs):
        DoPlots( [tp_pt]   , datasetList, True )
        DoPlots( [tp_pt_L] , datasetList, True )
        DoPlots( [tp_pt_C] , datasetList, True )
        DoPlots( [tp_pt_I] , datasetList, True )
        DoPlots( [tp_pt_F] , datasetList, True )
        DoPlots( [tp_eta]  , datasetList, True )
        DoPlots( [tp_eta_L], datasetList, True )
        DoPlots( [tp_eta_M], datasetList, True )
        DoPlots( [tp_eta_H], datasetList, True )
        DoPlots( [tp_phi]  , datasetList, True )
        DoPlots( [tp_z0]   , datasetList, True )
        DoPlots( [tp_d0]   , datasetList, True )

    if(bDoMatchedTPs):
        DoPlots( [match_tp_pt   ], datasetList, True )
        DoPlots( [match_tp_pt_L ], datasetList, True )
        DoPlots( [match_tp_pt_C ], datasetList, True )
        DoPlots( [match_tp_pt_I ], datasetList, True )
        DoPlots( [match_tp_pt_F ], datasetList, True )
        DoPlots( [match_tp_eta  ], datasetList, True )
        DoPlots( [match_tp_eta_L], datasetList, True )
        DoPlots( [match_tp_eta_M], datasetList, True )
        DoPlots( [match_tp_eta_H], datasetList, True )
        DoPlots( [match_tp_phi  ], datasetList, True )
        DoPlots( [match_tp_z0   ], datasetList, True )
        DoPlots( [match_tp_d0   ], datasetList, True )

    if(bDoEff):
        bColourPalette = True
        if len(datasetList) > 1:
            bColourPalette = False 
        DoPlots( [eff_pt]   , datasetList, bColourPalette)
        DoPlots( [eff_pt_L] , datasetList, bColourPalette)
        #DoPlots( [eff_pt_C] , datasetList, bColourPalette)
        #DoPlots( [eff_pt_I] , datasetList, bColourPalette)
        #DoPlots( [eff_pt_F] , datasetList, bColourPalette)
        DoPlots( [eff_eta ] , datasetList, bColourPalette)
        #DoPlots( [eff_eta_L], datasetList, bColourPalette)
        #DoPlots( [eff_eta_M], datasetList, bColourPalette)
        #DoPlots( [eff_eta_H], datasetList, bColourPalette)
        DoPlots( [eff_phi]  , datasetList, bColourPalette)
        DoPlots( [eff_z0 ]  , datasetList, bColourPalette)
        DoPlots( [eff_d0]   , datasetList, bColourPalette)
        if (len(datasetList) < 2):
            DoPlots( [eff_pt    , eff_pt_C, eff_pt_I, eff_pt_F]   , datasetList, True, "_CIF")
            DoPlots( [eff_eta   , eff_eta_L, eff_eta_M, eff_eta_H], datasetList, True, "_LMH")

    if(bDoMatchedTks):
        ### TTTrack (TP-Matched)
        #DoPlots( [match_trk_nstub]  , datasetList, True )
        #DoPlots( [match_trk_nstub_C], datasetList, True )
        #DoPlots( [match_trk_nstub_I], datasetList, True )
        #DoPlots( [match_trk_nstub_F], datasetList, True )
        DoPlots( [match_trk_nstub_C, match_trk_nstub_I, match_trk_nstub_F], datasetList, True, "IF")

        ### TTTrack (TP-Matched): ChiSq histograms (last bin is an overflow bin)
        DoPlots( [match_trk_chi2_L    , match_trk_chi2_M    , match_trk_chi2_H]    , datasetList, True, "MH")
        DoPlots( [match_trk_chi2_dof_L, match_trk_chi2_dof_M, match_trk_chi2_dof_H], datasetList, True, "MH")


        ### TTTrack (TP-Matched): ChiSq histograms (last bin is an overflow bin)
        #DoPlots( [match_trk_chi2    ], datasetList, True)
        #DoPlots( [match_trk_chi2_C_L], datasetList, True)
        #DoPlots( [match_trk_chi2_I_L], datasetList, True)
        #DoPlots( [match_trk_chi2_F_L], datasetList, True)
        DoPlots( [match_trk_chi2_C_L, match_trk_chi2_I_L, match_trk_chi2_F_L], datasetList, True, "_IL_FL")

        #DoPlots( [match_trk_chi2_C_M], datasetList, True)
        #DoPlots( [match_trk_chi2_I_M], datasetList, True)
        #DoPlots( [match_trk_chi2_F_M], datasetList, True)
        #DoPlots( [match_trk_chi2_C_H], datasetList, True)
        #DoPlots( [match_trk_chi2_I_H], datasetList, True)
        #DoPlots( [match_trk_chi2_F_H], datasetList, True)
        DoPlots( [match_trk_chi2_C_M, match_trk_chi2_I_M, match_trk_chi2_F_M], datasetList, True, "_IM_FM")

        ### TTTrack (TP-Matched): RedChiSq histograms (lastbin is an overflow bin)
        #DoPlots( [match_trk_chi2_dof    ], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_C_L], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_I_L], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_F_L], datasetList, True)
        DoPlots( [match_trk_chi2_dof_C_L, match_trk_chi2_dof_I_L, match_trk_chi2_dof_F_L], datasetList, True, "_IL_FL")

        #DoPlots( [match_trk_chi2_dof_C_M], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_I_M], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_F_M], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_C_H], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_I_H], datasetList, True)
        #DoPlots( [match_trk_chi2_dof_F_H], datasetList, True)
        DoPlots( [match_trk_chi2_dof_C_M, match_trk_chi2_dof_I_M, match_trk_chi2_dof_F_M], datasetList, True, "_ML_FM")

    if(bDoMatchedTks_Res):
        ### TTTrack (TP-Matched): Resolution histograms
        #DoPlots( [res_pt   ], datasetList, True)
        #DoPlots( [res_pt_C ], datasetList, True)
        #DoPlots( [res_pt_I ], datasetList, True)
        #DoPlots( [res_pt_F ], datasetList, True)
        DoPlots( [res_pt_C, res_pt_I, res_pt_F ], datasetList, True, "IF")

        #DoPlots( [res_ptRel   ], datasetList, True)
        #DoPlots( [res_ptRel_C ], datasetList, True)
        #DoPlots( [res_ptRel_I ], datasetList, True)
        #DoPlots( [res_ptRel_F ], datasetList, True)
        DoPlots( [res_ptRel_C, res_ptRel_I, res_ptRel_F ], datasetList, True, "IF")

        #DoPlots( [res_eta   ], datasetList, True)
        #DoPlots( [res_eta_C ], datasetList, True)
        #DoPlots( [res_eta_I ], datasetList, True)
        #DoPlots( [res_eta_F ], datasetList, True)
        DoPlots( [res_eta_C, res_eta_I, res_eta_F ], datasetList, True, "IF")

        #DoPlots( [res_phi   ], datasetList, True)
        #DoPlots( [res_phi_C ], datasetList, True)
        #DoPlots( [res_phi_I ], datasetList, True)
        #DoPlots( [res_phi_F ], datasetList, True)
        DoPlots( [res_phi_C, res_phi_I, res_phi_F ], datasetList, True, "IF")

        #DoPlots( [res_z0   ], datasetList, True)
        #DoPlots( [res_z0_C ], datasetList, True)
        #DoPlots( [res_z0_I ], datasetList, True)
        #DoPlots( [res_z0_F ], datasetList, True)
        DoPlots( [res_z0_C, res_z0_I, res_z0_F ], datasetList, True, "IF")

        #DoPlots( [res_d0   ], datasetList, True)
        #DoPlots( [res_d0_C ], datasetList, True)
        #DoPlots( [res_d0_I ], datasetList, True)
        #DoPlots( [res_d0_F ], datasetList, True)
        DoPlots( [res_d0_C, res_d0_I, res_d0_F ], datasetList, True, "IF")


    if(bDoResVsPt_Pt):
        for i in range(0, 20):
            #DoPlots( [ h_resVsPt_pt[i]  ], datasetList, True)
            #DoPlots( [ h_resVsPt_pt_C[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_pt_I[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_pt_F[i]], datasetList, True)
            DoPlots( [ h_resVsPt_pt_C[i] , h_resVsPt_pt_I[i], h_resVsPt_pt_F[i]], datasetList, True, "_IF")

    if(bDoResVsPt_PtRel):
        for i in range(0, 20):
            #DoPlots( [ h_resVsPt_ptRel[i]  ], datasetList, True)
            #DoPlots( [ h_resVsPt_ptRel_C[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_ptRel_I[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_ptRel_F[i]], datasetList, True)
            DoPlots( [ h_resVsPt_ptRel_C[i] , h_resVsPt_ptRel_I[i], h_resVsPt_ptRel_F[i]], datasetList, True, "_IF")

    if(bDoResVsPt_Z0):
        for i in range(0, 20):
            #DoPlots( [ h_resVsPt_z0[i]  ], datasetList, True)
            #DoPlots( [ h_resVsPt_z0_C[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_z0_I[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_z0_F[i]], datasetList, True)
            DoPlots( [ h_resVsPt_z0_C[i], h_resVsPt_z0_I[i], h_resVsPt_z0_F[i]], datasetList, True)

    if(bDoResVsPt_d0):
        for i in range(0, 20):
            #DoPlots( [ h_resVsPt_d0[i]  ], datasetList, True)
            #DoPlots( [ h_resVsPt_d0_C[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_d0_I[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_d0_F[i]], datasetList, True)
            DoPlots( [ h_resVsPt_d0_C[i], h_resVsPt_d0_I[i], h_resVsPt_d0_F[i]], datasetList, True)

    if(bDoResVsPt_Phi):
        for i in range(0, 20):
            #DoPlots( [ h_resVsPt_phi[i]  ], datasetList, True)
            #DoPlots( [ h_resVsPt_phi_C[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_phi_I[i]], datasetList, True)
            #DoPlots( [ h_resVsPt_phi_F[i]], datasetList, True)
            DoPlots( [ h_resVsPt_phi_C[i], h_resVsPt_phi_I[i], h_resVsPt_phi_F[i]], datasetList, True)

    if(bDoResVsPt_Eta):
        for i in range(0, 20):
            DoPlots( [ h_resVsPt_eta[i]]  , datasetList, True)
            DoPlots( [ h_resVsPt_eta_C[i]], datasetList, True)
            DoPlots( [ h_resVsPt_eta_I[i]], datasetList, True)
            DoPlots( [ h_resVsPt_eta_F[i]], datasetList, True)
            DoPlots( [ h_resVsPt_eta_C[i], h_resVsPt_eta_I[i], h_resVsPt_eta_F[i]], datasetList, True)

    if(bDoResVsEta_Eta):
        for i in range(0, 25):
            #DoPlots( [ h_resVsEta_eta[i] ]     , datasetList, True)
            #DoPlots( [ h_resVsEta_eta_L[i] ]   , datasetList, True)
            #DoPlots( [ h_resVsEta_eta_M[i] ]   , datasetList, True)
            #DoPlots( [ h_resVsEta_eta_H[i] ]   , datasetList, True)
            DoPlots( [ h_resVsEta_eta_L[i], h_resVsEta_eta_M[i], h_resVsEta_eta_H[i] ], datasetList, True, "MH")

    if(bDoResVsEta_Z0):
        for i in range(0, 25):
            #DoPlots( [ h_resVsEta_z0[i] ]      , datasetList, True)
            #DoPlots( [ h_resVsEta_z0_L[i] ]    , datasetList, True)
            #DoPlots( [ h_resVsEta_z0_M[i] ]    , datasetList, True)
            #DoPlots( [ h_resVsEta_z0_H[i] ]    , datasetList, True)
            DoPlots( [ h_resVsEta_z0_L[i], h_resVsEta_z0_M[i], h_resVsEta_z0_H[i] ], datasetList, True, "MH")

    if(bDoResVsEta_d0):
        for i in range(0, 25):
            #DoPlots( [ h_resVsEta_d0[i] ]      , datasetList, True)
            #DoPlots( [ h_resVsEta_d0_L[i] ]    , datasetList, True)
            #DoPlots( [ h_resVsEta_d0_M[i] ]    , datasetList, True)
            #DoPlots( [ h_resVsEta_d0_H[i] ]    , datasetList, True)
            DoPlots( [ h_resVsEta_d0_L[i], h_resVsEta_d0_M[i], h_resVsEta_d0_H[i] ], datasetList, True, "MH")

    if(bDoResVsEta_Phi):
        for i in range(0, 25):
            #DoPlots( [ h_resVsEta_phi[i] ]     , datasetList, True)
            #DoPlots( [ h_resVsEta_phi_L[i] ]   , datasetList, True)
            #DoPlots( [ h_resVsEta_phi_M[i] ]   , datasetList, True)
            #DoPlots( [ h_resVsEta_phi_H[i] ]   , datasetList, True)
            DoPlots( [ h_resVsEta_phi_L[i], h_resVsEta_phi_M[i], h_resVsEta_phi_H[i] ], datasetList, True, "MH")


    if(bDo2D):
        DoPlots( [ logChiSq_eta ]     , datasetList, True)
        DoPlots( [ logChiSq_dof_eta ] , datasetList, True)
        DoPlots( [ dz0_eta ]          , datasetList, True)
        DoPlots( [ deta_eta ]         , datasetList, True)
        DoPlots( [ dphi_eta ]         , datasetList, True)
        DoPlots( [ dpt_eta ]          , datasetList, True)
        DoPlots( [ dptRel_eta ]       , datasetList, True)


    if(bDo2DResVsPt_Pt):
        #DoPlots( [ resVsPt_pt  ]   , datasetList, True)
        #DoPlots( [ resVsPt_pt_C]   , datasetList, True)
        #DoPlots( [ resVsPt_pt_I]   , datasetList, True)
        #DoPlots( [ resVsPt_pt_F]   , datasetList, True)
        DoPlots( [ resVsPt_pt_C, resVsPt_pt_I, resVsPt_pt_F], datasetList, True, "IF")

    if(bDo2DResVsPt_PtRel):
        #DoPlots( [ resVsPt_ptRel  ]   , datasetList, True)
        #DoPlots( [ resVsPt_ptRel_C]   , datasetList, True)
        #DoPlots( [ resVsPt_ptRel_I]   , datasetList, True)
        #DoPlots( [ resVsPt_ptRel_F]   , datasetList, True)
        DoPlots( [ resVsPt_ptRel_C, resVsPt_ptRel_I, resVsPt_ptRel_F], datasetList, True, "IF")

    if(bDo2DMeanResVsPt_Pt):        
        #DoPlots( [ mresVsPt_pt  ]   , datasetList, True)
        #DoPlots( [ mresVsPt_pt_C]   , datasetList, True)
        #DoPlots( [ mresVsPt_pt_I]   , datasetList, True)
        #DoPlots( [ mresVsPt_pt_F]   , datasetList, True)
        DoPlots( [ mresVsPt_pt_C, mresVsPt_pt_I, mresVsPt_pt_F], datasetList, True, "IF")
        #DoPlots( [ mresVsPt_pt_C, mresVsPt_pt_I, mresVsPt_pt_F], ["PiPlus", "PiMinus"], False, "IF")

    if(bDo2DMeanResVsEta_Eta):
        #DoPlots( [ mresVsEta_eta  ], datasetList, True)
        #DoPlots( [ mresVsEta_eta_L], datasetList, True)
        #DoPlots( [ mresVsEta_eta_M], datasetList, True)
        #DoPlots( [ mresVsEta_eta_H], datasetList, True)
        DoPlots( [ mresVsEta_eta_L, mresVsEta_eta_M, mresVsEta_eta_H], datasetList, True, "MH")
        #DoPlots( [ mresVsEta_eta_L, mresVsEta_eta_M, mresVsEta_eta_H], ["PiPlus", "PiMinus"], True, "MH")

    if(bDo2DResVsPt_Eta):
        #DoPlots( [ resVsPt_eta  ], datasetList, True)
        #DoPlots( [ resVsPt_eta_C], datasetList, True)
        #DoPlots( [ resVsPt_eta_I], datasetList, True)
        #DoPlots( [ resVsPt_eta_F], datasetList, True)
        DoPlots( [ resVsPt_eta_C, resVsPt_eta_I, resVsPt_eta_F], datasetList, True, "IF")

    if(bDo2DResVsPt_Phi):
        #DoPlots( [ resVsPt_phi  ]   , datasetList, True)
        #DoPlots( [ resVsPt_phi_C]   , datasetList, True)
        #DoPlots( [ resVsPt_phi_I]   , datasetList, True)
        #DoPlots( [ resVsPt_phi_F]   , datasetList, True)
        DoPlots( [ resVsPt_phi_C, resVsPt_phi_I, resVsPt_phi_F], datasetList, True, "IF")

    if(bDo2DResVsPt_Z0):
        #DoPlots( [ resVsPt_z0  ]   , datasetList, True)
        #DoPlots( [ resVsPt_z0_C]   , datasetList, True)
        #DoPlots( [ resVsPt_z0_I]   , datasetList, True)
        #DoPlots( [ resVsPt_z0_F]   , datasetList, True)
        DoPlots( [ resVsPt_z0_C, resVsPt_z0_I, resVsPt_z0_F], datasetList, True, "IF")

    if(bDo2DResVsPt_d0):
        #DoPlots( [ resVsPt_d0  ]   , datasetList, True)
        #DoPlots( [ resVsPt_d0_C]   , datasetList, True)
        #DoPlots( [ resVsPt_d0_I]   , datasetList, True)
        #DoPlots( [ resVsPt_d0_F]   , datasetList, True)
        DoPlots( [ resVsPt_d0_C, resVsPt_d0_I, resVsPt_d0_F], datasetList, True, "IF")

                    
    if(bDo2DResVsEta_Pt):
        #DoPlots( [ resVsEta_pt  ], datasetList, True)
        #DoPlots( [ resVsEta_pt_L], datasetList, True)
        #DoPlots( [ resVsEta_pt_M], datasetList, True)
        #DoPlots( [ resVsEta_pt_H], datasetList, True)
        DoPlots( [ resVsEta_pt_L, resVsEta_pt_M, resVsEta_pt_H], datasetList, True, "MH")

    if(bDo2DResVsEta_PtRel):
        #DoPlots( [ resVsEta_ptRel  ], datasetList, True)
        #DoPlots( [ resVsEta_ptRel_L], datasetList, True)
        #DoPlots( [ resVsEta_ptRel_M], datasetList, True)
        #DoPlots( [ resVsEta_ptRel_H], datasetList, True)
        DoPlots( [ resVsEta_ptRel_L, resVsEta_ptRel_M, resVsEta_ptRel_H], datasetList, True, "MH")

    if(bDo2DResVsEta_Eta):
        #DoPlots( [ resVsEta_eta  ], datasetList, True)
        #DoPlots( [ resVsEta_eta_L], datasetList, True)
        #DoPlots( [ resVsEta_eta_M], datasetList, True)
        #DoPlots( [ resVsEta_eta_H], datasetList, True)
        DoPlots( [ resVsEta_eta_L, resVsEta_eta_M, resVsEta_eta_H], datasetList, True, "MH")

    if(bDo2DResVsEta_Phi):
        #DoPlots( [ resVsEta_phi  ], datasetList, True)
        #DoPlots( [ resVsEta_phi_L], datasetList, True)
        #DoPlots( [ resVsEta_phi_M], datasetList, True)
        #DoPlots( [ resVsEta_phi_H], datasetList, True)
        DoPlots( [ resVsEta_phi_L, resVsEta_phi_M, resVsEta_phi_H], datasetList, True, "MH")

    if(bDo2DResVsEta_Z0):
        #DoPlots( [ resVsEta_z0  ], datasetList, True)
        #DoPlots( [ resVsEta_z0_L], datasetList, True)
        #DoPlots( [ resVsEta_z0_M], datasetList, True)
        #DoPlots( [ resVsEta_z0_H], datasetList, True)
        DoPlots( [ resVsEta_z0_L, resVsEta_z0_M, resVsEta_z0_H], datasetList, True, "MH")

    if(bDo2DResVsEta_d0):
        #DoPlots( [ resVsEta_d0  ], datasetList, True)
        #DoPlots( [ resVsEta_d0_L], datasetList, True)
        #DoPlots( [ resVsEta_d0_M], datasetList, True)
        #DoPlots( [ resVsEta_d0_H], datasetList, True)
        DoPlots( [ resVsEta_d0_L, resVsEta_d0_M, resVsEta_d0_H], datasetList, True, "MH")

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
### Scatter residual plots
bDo2D = False

### Tracks and Tracking Particles
bTPs                     = False
bMatchedTPs              = False
bMatchedTks              = False  #talk (1 dataset)
bMatchedTks_EtaRange     = False  #talk (1 dataset)

### Tracks
bTkChiSq                 = False
bTkChiSq_EtaRange        = False
bTkRedChiSq              = False
bTkRedChiSq_EtaRange     = False

### Inclusive Residuals & Residuals
bResolutions             = False #talk
bResiduals_EtaRange      = False #talk (datasets)
bResiduals_CIF           = False #talk (1 dataset)
            
### Tracking Efficiency
bTkEfficiency            = False #talk
bTkEfficiency_Vs_Pt      = False 
bTkEfficiency_Vs_Eta     = False 
bTkEfficiency_Vs_d0z0phi = False #talk

### pt Residuals
pt_Resolution_Vs_Pt                    = False #talk (datasets)
pt_Residual_Vs_Pt_PtRange              = False #resolution fits
pt_Resolution_Vs_Pt_EtaRange           = False #talk (datasets)
pt_Residual_Vs_Pt_EtaRange             = False #resolution fits
pt_Resolution_Vs_Pt_CIF                = False #talk (1 dataset)
pt_Residual_Vs_Pt_CIF_PtRange          = False #resolution fits 
pt_Resolution_Vs_Eta                   = False #talk (datasets)
pt_Residual_Vs_Eta_EtaRange            = False #resolution fits
pt_Resolution_Vs_Eta_PtRange           = False #talk (datasets)
pt_Residual_Vs_Eta_PtRange_EtaRange    = False #resolution fits
pt_Resolution_Vs_Eta_LMH               = False #talk (1 dataset)
pt_Residual_Vs_Eta_LMH_EtaRange        = False #resolution fits (avoid)
        
### ptRel Residuals
ptRel_Resolution_Vs_Pt                    = False #talk (datasets)
ptRel_Residual_Vs_Pt_PtRange              = False #resolution fits
ptRel_Resolution_Vs_Pt_EtaRange           = False #talk (datasets)
ptRel_Residual_Vs_Pt_EtaRange             = False #resolution fits
ptRel_Resolution_Vs_Pt_CIF                = False #talk (1 dataset)
ptRel_Residual_Vs_Pt_CIF_PtRange          = False #resolution fits 
ptRel_Resolution_Vs_Eta                   = False #talk (datasets)
ptRel_Residual_Vs_Eta_EtaRange            = False #resolution fits
ptRel_Resolution_Vs_Eta_PtRange           = False #talk (datasets)
ptRel_Residual_Vs_Eta_PtRange_EtaRange    = False #resolution fits
ptRel_Resolution_Vs_Eta_LMH               = False #talk (1 dataset)
ptRel_Residual_Vs_Eta_LMH_EtaRange        = False #resolution fits (avoid)

### phi Residuals
phi_Resolution_Vs_Pt                    = False #talk (datasets)
phi_Residual_Vs_Pt_PtRange              = False #resolution fits (1 dataset)
phi_Resolution_Vs_Pt_EtaRange           = False #talk (datasets)
phi_Residual_Vs_Pt_EtaRange             = False #resolution fits (1 dataset)
phi_Resolution_Vs_Pt_CIF                = False #talk (1 dataset)
phi_Residual_Vs_Pt_CIF_PtRange          = False #resolution fits (avoid)
phi_Resolution_Vs_Eta                   = False #talk (datasets)
phi_Residual_Vs_Eta_EtaRange            = False #resolution fits (1 dataset)
phi_Resolution_Vs_Eta_PtRange           = False #talk (datasets)
phi_Residual_Vs_Eta_PtRange_EtaRange    = False #resolution fits (1 dataset)
phi_Resolution_Vs_Eta_LMH               = False #talk (1 dataset)
phi_Residual_Vs_Eta_LMH_EtaRange        = False #resolution fits (avoid)

### eta Residuals
eta_Resolution_Vs_Pt                    = True  #talk (datasets)
eta_Residual_Vs_Pt_PtRange              = False #resolution fits (1 dataset)
eta_Resolution_Vs_Pt_EtaRange           = False #talk (datasets)
eta_Residual_Vs_Pt_EtaRange             = False #resolution fits (1 dataset)
eta_Resolution_Vs_Pt_CIF                = False #talk (1 dataset)
eta_Residual_Vs_Pt_CIF_PtRange          = False #resolution fits (avoid)
eta_Resolution_Vs_Eta                   = False #talk (datasets)
eta_Residual_Vs_Eta_EtaRange            = False #resolution fits (1 dataset)
eta_Resolution_Vs_Eta_PtRange           = True  #talk (datasets)
eta_Residual_Vs_Eta_PtRange_EtaRange    = False #resolution fits (1 dataset)
eta_Resolution_Vs_Eta_LMH               = False #talk (1 dataset)
eta_Residual_Vs_Eta_LMH_EtaRange        = False #resolution fits (avoid)

### z0 Residuals
z0_Resolution_Vs_Pt                    = True  #talk (datasets)
z0_Residual_Vs_Pt_PtRange              = False #resolution fits (1 dataset)
z0_Resolution_Vs_Pt_EtaRange           = False #talk (datasets)
z0_Residual_Vs_Pt_EtaRange             = False #resolution fits (1 dataset)
z0_Resolution_Vs_Pt_CIF                = False #talk (1 dataset)
z0_Residual_Vs_Pt_CIF_PtRange          = False #resolution fits (avoid)
z0_Resolution_Vs_Eta                   = True  #talk (datasets)
z0_Residual_Vs_Eta_EtaRange            = False #resolution fits (1 dataset)
z0_Resolution_Vs_Eta_PtRange           = False #talk (datasets)
z0_Residual_Vs_Eta_PtRange_EtaRange    = False #resolution fits (1 dataset)
z0_Resolution_Vs_Eta_LMH               = False #talk (1 dataset)
z0_Residual_Vs_Eta_LMH_EtaRange        = False #resolution fits (avoid)

### d0 Residuals
d0_Resolution_Vs_Pt                    = True  #talk (datasets)
d0_Residual_Vs_Pt_PtRange              = False #resolution fits (1 dataset)
d0_Resolution_Vs_Pt_EtaRange           = False #talk (datasets)
d0_Residual_Vs_Pt_EtaRange             = False #resolution fits (1 dataset)
d0_Resolution_Vs_Pt_CIF                = False #talk (1 dataset)
d0_Residual_Vs_Pt_CIF_PtRange          = False #resolution fits  (avoid)
d0_Resolution_Vs_Eta                   = True  #talk (datasets)
d0_Residual_Vs_Eta_EtaRange            = False #resolution fits (1 dataset)
d0_Resolution_Vs_Eta_PtRange           = False #talk (datasets)
d0_Residual_Vs_Eta_PtRange_EtaRange    = False #resolution fits (1 dataset)
d0_Resolution_Vs_Eta_LMH               = False #talk (1 dataset)
d0_Residual_Vs_Eta_LMH_EtaRange        = False #resolution fits (avoid)

###############################################################
### General Settings
###############################################################
nFitParams    = "5"
bDoPixelTks   = True
bRatio        = True
if bDoPixelTks:
    ext       = "_Pixel"
    tkTypeDir = "pixTks"
else:
    ext       = ""
    tkTypeDir = "ttTks"

datasetList   = ["PiMinus"] #SingleMuMinus
savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_03July2015/figures/" + nFitParams + "FitParams/" + tkTypeDir + "/" + datasetList[0] + "/"
#savePath      = ""
saveFormats   = ["png", "pdf"]


inputPath_0p25x     = "Macros/Tracking/results/" + nFitParams + "FitParams/0p25xWindowsSF/Tracking_Histograms_"
inputPath_0p5x      = "Macros/Tracking/results/" + nFitParams + "FitParams/0p5xWindowsSF/Tracking_Histograms_"
inputPath_1p0x      = "Macros/Tracking/results/" + nFitParams + "FitParams/1p0xWindowsSF/Tracking_Histograms_"
#inputPath_1p0x      = "Macros/Tracking/Tracking_Histograms_"
inputPath_1p5x      = "Macros/Tracking/results/" + nFitParams + "FitParams/1p5xWindowsSF/Tracking_Histograms_"
inputPath_2p0x      = "Macros/Tracking/results/" + nFitParams + "FitParams/2p0xWindowsSF/Tracking_Histograms_"
inputPath_4p0x      = "Macros/Tracking/results/" + nFitParams + "FitParams/4p0xWindowsSF/Tracking_Histograms_"
inputPath_0p5xB2p0E = "Macros/Tracking/results/" + nFitParams + "FitParams/0p5xBarrel_2p0xEndcap_WindowsSF/Tracking_Histograms_"

datasetPaths_0p25x  = {}
datasetPaths_0p25x["MinBias"]                 = inputPath_0p25x + "MinBias" + ext + ".root"
datasetPaths_0p25x["VBF"]                     = inputPath_0p25x + "VBF" + ext + ".root"
datasetPaths_0p25x["PiPlus"]                  = inputPath_0p25x + "PiPlus" + ext + ".root"
datasetPaths_0p25x["PiMinus"]                 = inputPath_0p25x + "PiMinus" + ext + ".root"
datasetPaths_0p25x["SingleTauGun1p"]          = inputPath_0p25x + "SingleTauGun1p" + ext + ".root"
datasetPaths_0p25x["DiTauGun3p"]              = inputPath_0p25x + "DiTauGun3p" + ext + ".root"
datasetPaths_0p25x["TTbar"]                   = inputPath_0p25x + "TTBar" + ext + ".root"
datasetPaths_0p25x["HPlus160"]                = inputPath_0p25x + "HPlus160" + ext + ".root"
datasetPaths_0p25x["HPlus200"]                = inputPath_0p25x + "HPlus200" + ext + ".root"
datasetPaths_0p25x["SingleElectron"]          = inputPath_0p25x + "SingleElectron" + ext + ".root"
datasetPaths_0p25x["SinglePositron"]          = inputPath_0p25x + "SinglePositron" + ext + ".root"
datasetPaths_0p25x["SingleMuPlus"  ]          = inputPath_0p25x + "SingleMuPlus" + ext + ".root"
datasetPaths_0p25x["SingleMuMinus" ]          = inputPath_0p25x + "SingleMuMinus" + ext + ".root"
datasetPaths_0p25x["SinglePhoton"  ]          = inputPath_0p25x + "SinglePhoton" + ext + ".root"  

datasetPaths_0p5x  = {}
datasetPaths_0p5x["MinBias"]                 = inputPath_0p5x + "MinBias" + ext + ".root"
datasetPaths_0p5x["VBF"]                     = inputPath_0p5x + "VBF" + ext + ".root"
datasetPaths_0p5x["PiPlus"]                  = inputPath_0p5x + "PiPlus" + ext + ".root"
datasetPaths_0p5x["PiMinus"]                 = inputPath_0p5x + "PiMinus" + ext + ".root"
datasetPaths_0p5x["SingleTauGun1p"]          = inputPath_0p5x + "SingleTauGun1p" + ext + ".root"
datasetPaths_0p5x["DiTauGun3p"]              = inputPath_0p5x + "DiTauGun3p" + ext + ".root"
datasetPaths_0p5x["TTbar"]                   = inputPath_0p5x + "TTBar" + ext + ".root"
datasetPaths_0p5x["HPlus160"]                = inputPath_0p5x + "HPlus160" + ext + ".root"
datasetPaths_0p5x["HPlus200"]                = inputPath_0p5x + "HPlus200" + ext + ".root"
datasetPaths_0p5x["SingleElectron"]          = inputPath_0p5x + "SingleElectron" + ext + ".root"
datasetPaths_0p5x["SinglePositron"]          = inputPath_0p5x + "SinglePositron" + ext + ".root"
datasetPaths_0p5x["SingleMuPlus"  ]          = inputPath_0p5x + "SingleMuPlus" + ext + ".root"
datasetPaths_0p5x["SingleMuMinus" ]          = inputPath_0p5x + "SingleMuMinus" + ext + ".root"
datasetPaths_0p5x["SinglePhoton"  ]          = inputPath_0p5x + "SinglePhoton" + ext + ".root"  

datasetPaths_1p0x  = {}
datasetPaths_1p0x["MinBias"]                   = inputPath_1p0x + "MinBias" + ext + ".root"
datasetPaths_1p0x["VBF"]                       = inputPath_1p0x + "VBF" + ext + ".root"
datasetPaths_1p0x["PiPlus"]                    = inputPath_1p0x + "PiPlus" + ext + ".root"
datasetPaths_1p0x["PiMinus"]                   = inputPath_1p0x + "PiMinus" + ext + ".root"
datasetPaths_1p0x["SingleTauGun1p"]            = inputPath_1p0x + "SingleTauGun1p" + ext + ".root"
datasetPaths_1p0x["DiTauGun3p"]                = inputPath_1p0x + "DiTauGun3p" + ext + ".root"
datasetPaths_1p0x["TTbar"]                     = inputPath_1p0x + "TTBar" + ext + ".root"
datasetPaths_1p0x["HPlus160"]                  = inputPath_1p0x + "HPlus160" + ext + ".root"
datasetPaths_1p0x["HPlus200"]                  = inputPath_1p0x + "HPlus200" + ext + ".root"
datasetPaths_1p0x["SingleElectron"]            = inputPath_1p0x + "SingleElectron" + ext + ".root"
datasetPaths_1p0x["SinglePositron"]            = inputPath_1p0x + "SinglePositron" + ext + ".root"
datasetPaths_1p0x["SingleMuPlus"  ]            = inputPath_1p0x + "SingleMuPlus" + ext + ".root"
datasetPaths_1p0x["SingleMuMinus" ]            = inputPath_1p0x + "SingleMuMinus" + ext + ".root"
datasetPaths_1p0x["SinglePhoton"  ]            = inputPath_1p0x + "SinglePhoton" + ext + ".root"  

datasetPaths_1p5x  = {}
datasetPaths_1p5x["MinBias"]                   = inputPath_1p5x + "MinBias" + ext + ".root"
datasetPaths_1p5x["VBF"]                       = inputPath_1p5x + "VBF" + ext + ".root"
datasetPaths_1p5x["PiPlus"]                    = inputPath_1p5x + "PiPlus" + ext + ".root"
datasetPaths_1p5x["PiMinus"]                   = inputPath_1p5x + "PiMinus" + ext + ".root"
datasetPaths_1p5x["SingleTauGun1p"]            = inputPath_1p5x + "SingleTauGun1p" + ext + ".root"
datasetPaths_1p5x["DiTauGun3p"]                = inputPath_1p5x + "DiTauGun3p" + ext + ".root"
datasetPaths_1p5x["TTbar"]                     = inputPath_1p5x + "TTBar" + ext + ".root"
datasetPaths_1p5x["HPlus160"]                  = inputPath_1p5x + "HPlus160" + ext + ".root"
datasetPaths_1p5x["HPlus200"]                  = inputPath_1p5x + "HPlus200" + ext + ".root"
datasetPaths_1p5x["SingleElectron"]            = inputPath_1p5x + "SingleElectron" + ext + ".root"
datasetPaths_1p5x["SinglePositron"]            = inputPath_1p5x + "SinglePositron" + ext + ".root"
datasetPaths_1p5x["SingleMuPlus"  ]            = inputPath_1p5x + "SingleMuPlus" + ext + ".root"
datasetPaths_1p5x["SingleMuMinus" ]            = inputPath_1p5x + "SingleMuMinus" + ext + ".root"
datasetPaths_1p5x["SinglePhoton"  ]            = inputPath_1p5x + "SinglePhoton" + ext + ".root"  

datasetPaths_2p0x  = {}
datasetPaths_2p0x["MinBias"]                   = inputPath_2p0x + "MinBias" + ext + ".root"
datasetPaths_2p0x["VBF"]                       = inputPath_2p0x + "VBF" + ext + ".root"
datasetPaths_2p0x["PiPlus"]                    = inputPath_2p0x + "PiPlus" + ext + ".root"
datasetPaths_2p0x["PiMinus"]                   = inputPath_2p0x + "PiMinus" + ext + ".root"
datasetPaths_2p0x["SingleTauGun1p"]            = inputPath_2p0x + "SingleTauGun1p" + ext + ".root"
datasetPaths_2p0x["DiTauGun3p"]                = inputPath_2p0x + "DiTauGun3p" + ext + ".root"
datasetPaths_2p0x["TTbar"]                     = inputPath_2p0x + "TTBar" + ext + ".root"
datasetPaths_2p0x["HPlus160"]                  = inputPath_2p0x + "HPlus160" + ext + ".root"
datasetPaths_2p0x["HPlus200"]                  = inputPath_2p0x + "HPlus200" + ext + ".root"
datasetPaths_2p0x["SingleElectron"]            = inputPath_2p0x + "SingleElectron" + ext + ".root"
datasetPaths_2p0x["SinglePositron"]            = inputPath_2p0x + "SinglePositron" + ext + ".root"
datasetPaths_2p0x["SingleMuPlus"  ]            = inputPath_2p0x + "SingleMuPlus" + ext + ".root"
datasetPaths_2p0x["SingleMuMinus" ]            = inputPath_2p0x + "SingleMuMinus" + ext + ".root"
datasetPaths_2p0x["SinglePhoton"  ]            = inputPath_2p0x + "SinglePhoton" + ext + ".root"  

datasetPaths_4p0x  = {}
datasetPaths_4p0x["MinBias"]                   = inputPath_4p0x + "MinBias" + ext + ".root"
datasetPaths_4p0x["VBF"]                       = inputPath_4p0x + "VBF" + ext + ".root"
datasetPaths_4p0x["PiPlus"]                    = inputPath_4p0x + "PiPlus" + ext + ".root"
datasetPaths_4p0x["PiMinus"]                   = inputPath_4p0x + "PiMinus" + ext + ".root"
datasetPaths_4p0x["SingleTauGun1p"]            = inputPath_4p0x + "SingleTauGun1p" + ext + ".root"
datasetPaths_4p0x["DiTauGun3p"]                = inputPath_4p0x + "DiTauGun3p" + ext + ".root"
datasetPaths_4p0x["TTbar"]                     = inputPath_4p0x + "TTBar" + ext + ".root"
datasetPaths_4p0x["HPlus160"]                  = inputPath_4p0x + "HPlus160" + ext + ".root"
datasetPaths_4p0x["HPlus200"]                  = inputPath_4p0x + "HPlus200" + ext + ".root"
datasetPaths_4p0x["SingleElectron"]            = inputPath_4p0x + "SingleElectron" + ext + ".root"
datasetPaths_4p0x["SinglePositron"]            = inputPath_4p0x + "SinglePositron" + ext + ".root"
datasetPaths_4p0x["SingleMuPlus"  ]            = inputPath_4p0x + "SingleMuPlus" + ext + ".root"
datasetPaths_4p0x["SingleMuMinus" ]            = inputPath_4p0x + "SingleMuMinus" + ext + ".root"
datasetPaths_4p0x["SinglePhoton"  ]            = inputPath_4p0x + "SinglePhoton" + ext + ".root"  

datasetPaths_0p5xB2p0E  = {}
datasetPaths_0p5xB2p0E["MinBias"]              = inputPath_0p5xB2p0E + "MinBias" + ext + ".root"
datasetPaths_0p5xB2p0E["VBF"]                  = inputPath_0p5xB2p0E + "VBF" + ext + ".root"
datasetPaths_0p5xB2p0E["PiPlus"]               = inputPath_0p5xB2p0E + "PiPlus" + ext + ".root"
datasetPaths_0p5xB2p0E["PiMinus"]              = inputPath_0p5xB2p0E + "PiMinus" + ext + ".root"
datasetPaths_0p5xB2p0E["SingleTauGun1p"]       = inputPath_0p5xB2p0E + "SingleTauGun1p" + ext + ".root"
datasetPaths_0p5xB2p0E["DiTauGun3p"]           = inputPath_0p5xB2p0E + "DiTauGun3p" + ext + ".root"
datasetPaths_0p5xB2p0E["TTbar"]                = inputPath_0p5xB2p0E + "TTBar" + ext + ".root"
datasetPaths_0p5xB2p0E["HPlus160"]             = inputPath_0p5xB2p0E + "HPlus160" + ext + ".root"
datasetPaths_0p5xB2p0E["HPlus200"]             = inputPath_0p5xB2p0E + "HPlus200" + ext + ".root"
datasetPaths_0p5xB2p0E["SingleElectron"]       = inputPath_0p5xB2p0E + "SingleElectron" + ext + ".root"
datasetPaths_0p5xB2p0E["SinglePositron"]       = inputPath_0p5xB2p0E + "SinglePositron" + ext + ".root"
datasetPaths_0p5xB2p0E["SingleMuPlus"  ]       = inputPath_0p5xB2p0E + "SingleMuPlus" + ext + ".root"
datasetPaths_0p5xB2p0E["SingleMuMinus" ]       = inputPath_0p5xB2p0E + "SingleMuMinus" + ext + ".root"
datasetPaths_0p5xB2p0E["SinglePhoton"  ]       = inputPath_0p5xB2p0E + "SinglePhoton" + ext + ".root"  

###############################################################
### Definitions
###############################################################
yMin          =  0.0
ptMax         = 50.0
etaMax        =  2.4 #2.5
nPt_Range     = int(ptMax/5.0)
nEta_Range    = 12 #25
C = "|#eta| #leq 0.8"
I = "0.8 < |#eta| #leq 1.6"
F = "|#eta| > 1.6"
L = "p_{T} #leq 6 GeV"
M = "6 < p_{T} #leq 16 GeV"
H = "p_{T} > 16 GeV"
EtaLines   = [-1.6, -0.8, +0.8, +1.6]
EtaRange = [[-etaMax, -1.6, ROOT.kRed+1], [+etaMax, +1.6, ROOT.kRed+1], [-1.6, -0.8, ROOT.kYellow-4], [+0.8, +1.6, ROOT.kYellow-4], [-0.8, +0.8, ROOT.kGreen+1] ]

############################################################### 
### Histogram Options
############################################################### 
Pt = {
    "xLabel": "p_{T}"           , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": ""         , "yMin": 1E00, "yMax": None , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtL = {
    "xLabel": "p_{T}"           , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": +5.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f" , "yUnits": ""         , "yMin": 1E00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


Eta = {
    "xLabel": "#eta"           , "xUnits": ""     , "xMin": -etaMax , "xMax": +etaMax, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""     , "yMin": +1e00, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


Phi = {
    "xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.2 , "xMax": +3.2, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": +1E00, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

z0 = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -25.0, "xMax": +25.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": +1E00, "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


d0 = {
    "xLabel": "d_{0}"          , "xUnits": "#mum", "xMin": -40 , "xMax": +40 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": ""    , "yMin": +1e0, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


NStubs = {
    "xLabel": "Stubs"          , "xUnits": "", "xMin": +3.0, "xMax": 10.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-4, "yMax": +2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.25, "yLegMax": 0.35
}


ChiSq  = {
    "xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +150.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f" , "yUnits": ""    , "yMin": 1E-04, "yMax": +1e0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False,  "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


RedChiSq = {
    "xLabel": "#chi^{2}_{N}"   , "xUnits": "", "xMin": +0.00, "xMax": +10.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f", "yUnits": "", "yMin": 1e-04, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.9, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}



PtRes = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual" , "xUnits": "GeVc^{-1}", "xMin": -5.0, "xMax": +5.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""         , "yMin": +1e0, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PtResAlt = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual" , "xUnits": "GeVc^{-1}", "xMin": -5.0, "xMax": +5.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f", "yUnits": ""         , "yMin": +1e0, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True ,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PtResRel = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual / p_{T}", "xUnits": "", "xMin": -0.12, "xMax": +0.12, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.3f"       , "yUnits": "", "yMin": +1e+0, "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


EtaRes = { #residual = (L1 - Sim)
    "xLabel": "#eta residual"   , "xUnits": "", "xMin": -0.003, "xMax": 0.003, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.5f" , "yUnits": "", "yMin": yMin  , "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PhiRes = { #residual = (L1 - Sim)
    "xLabel": "#phi residual"  , "xUnits": "rads", "xMin": -0.01, "xMax": +0.01, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.4f", "yUnits": ""    , "yMin": +1E00 , "yMax": None  , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True ,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


Z0Res = { #residual =  (L1 - Sim)
    "xLabel": "z_{0} residual" , "xUnits": "cm", "xMin": -0.04, "xMax": +0.04, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.3f", "yUnits": ""  , "yMin": yMin, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


d0Res = { #residual =  (L1 - Sim)
    "xLabel": "d_{0} residual" , "xUnits": "cm", "xMin": -0.035, "xMax": +0.035, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.4f", "yUnits": ""  , "yMin": +yMin, "yMax": 0.08 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


ChiSqEta = {
    "xLabel": "#eta / %0.2f"        , "xUnits": "", "xMin": -etaMax , "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "log #chi^{2} / %0.2f", "yUnits": "", "yMin": +0.0 , "yMax": +5.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

RedChiSqEta = {
    "xLabel": "#eta / %0.2f"            , "xUnits": "", "xMin": -etaMax , "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "log #chi^{2}_{N} / %0.2f", "yUnits": "", "yMin": +0.0 , "yMax": +5.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


dZ0Eta = {
    "xLabel": "#eta / %0.2f"             , "xUnits": ""  , "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutBoxes": EtaRange, "gridX": True, "logX": False,
    "yLabel": "|z_{0} residual| / % 0.2f", "yUnits": "cm", "yMin":  0.0   , "yMax": +1.2   , "binWidthY": None, "yCutBoxes": []      , "gridY": True, "logY": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False, "drawOptions": "COLZ", 
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dd0Eta = {
    "xLabel": "#eta / %0.2f"             , "xUnits": ""  , "xMin": -etaMax, "xMax": +etaMax, "xCutLines": [], "xCutBoxes": EtaRange,  "gridX": True, "logX": False,
    "yLabel": "|d_{0} residual| / % 0.2f", "yUnits": "cm", "yMin":  0.0   , "yMax": +1.2   , "yCutLines": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dEtaEta = {
    "xLabel": "#eta / %0.2f"             , "xUnits": "", "xMin": -etaMax, "xMax": +etaMax  , "binWidthX": None, "xCutLines": [] , "gridX": True, "logX": False, 
    "yLabel": "|#eta residual| / % 0.4f" , "yUnits": "", "yMin":  0.0, "yMax": +0.015, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dPhiEta = {
    "xLabel": "#eta / %0.2f"                     , "xUnits": ""    , "xMin": -etaMax, "xMax": +etaMax  , "binWidthX": None, "xCutLines": [] , "gridX": True, "logX": False, 
    "yLabel": "|#phi residual| / % 0.4f" , "yUnits": "rads", "yMin":  0.0, "yMax": +0.010, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


dPtRelEta = {
    "xLabel": "#eta / %0.2f"                    , "xUnits": "", "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutLines": [] , "gridX": True, "logX": False, 
    "yLabel": "|p_{T} residual/p_{T}| / % 0.2f" , "yUnits": "", "yMin":  0.0, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}

dPtEta = {
    "xLabel": "#eta / %0.2f"              , "xUnits": ""         , "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "|p_{T} residual| / % 0.2f" , "yUnits": "GeVc^{-1}", "yMin":  0.0, "yMax": +5.0, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


PtPtRes = {
    "xLabel": "p_{T}"                    , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax , "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "p_{T} resolution / % 0.0f", "yUnits": ""         , "yMin": 0.0 , "yMax": etaMax, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PtEtaRes = {
    "xLabel": "p_{T}"                   , "xUnits": "GeVc^{-1}", "xMin": 0.0, "xMax": ptMax  , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#eta resolution / % 0.2f", "yUnits": ""         , "yMin": 0.0, "yMax": +0.0055, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PtPhiRes = {
    "xLabel": "p_{T}"                   , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#phi resolution / % 0.0f", "yUnits": ""         , "yMin": 0.00, "yMax": 0.005, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtZ0Res = {
    "xLabel": "p_{T}"                         , "xUnits": "GeVc^{-1}", "xMin": 0.0, "xMax": ptMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "z_{0} resolution (cm) / % 0.0f", "yUnits": ""         , "yMin": 0.0, "yMax": +0.3, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

Ptd0Res = {
    "xLabel": "p_{T}"                         , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax, "binWidthX": None, "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "d_{0} resolution (cm) / % 0.0f", "yUnits": ""         , "yMin": 0.0 , "yMax": 0.08 , "binWidthY": None, "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


PtPtResRel = {
    "xLabel": "p_{T}"                            , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False,
    "yLabel": "p_{T} resolution / p_{T} / % 0.1f", "yUnits": ""         , "yMin": 0.00, "yMax": None , "binWidthY": None, "yCutLines": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaPtRes = {
    "xLabel": "#eta"                     , "xUnits": "", "xMin": -etaMax , "xMax": +etaMax , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "p_{T} resolution / % 0.1f", "yUnits": "", "yMin": 0.0 , "yMax":  None , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaEtaRes = {
    "xLabel": "#eta"                    , "xUnits": "", "xMin": 0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#eta resolution / % 0.2f", "yUnits": "", "yMin": 0.0, "yMax": 0.0055, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaPhiRes = {
    "xLabel": "#eta"                    , "xUnits": "", "xMin": -etaMax, "xMax": +etaMax , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#phi resolution / % 0.2f", "yUnits": "", "yMin": 0.0 , "yMax":  0.005, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaZ0Res = {
    "xLabel": "#eta"                          , "xUnits": "", "xMin": 0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "z_{0} resolution (cm) / % 0.2f", "yUnits": "", "yMin": 0.0, "yMax": +0.3, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

Etad0Res = {
    "xLabel": "#eta"                          , "xUnits": "", "xMin": +0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "d_{0} resolution (cm) / % 0.2f", "yUnits": "", "yMin": +0.0, "yMax": +0.08  , "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

EtaPtResRel = {
    "xLabel": "#eta"                             , "xUnits": "" , "xMin": +0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "p_{T} resolution / p_{T} / % 0.1f", "yUnits": "" , "yMin": +0.0, "yMax": +0.1, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


EffPt = {
    "xLabel": "p_{T}"              , "xUnits": "GeVc^{-1}", "xMin": 0.0 , "xMax": +50.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Efficiency / %0.0f" , "yUnits": ""         , "yMin": 0.00, "yMax":  +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False ,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}

EffPtL = {
    "xLabel": "p_{T}"              , "xUnits": "GeVc^{-1}", "xMin": 0.0, "xMax": +5.0  , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Efficiency / %0.1f" , "yUnits": ""     , "yMin": 0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


EffEta = {
    "xLabel": "#eta"               , "xUnits": "" , "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.2f" , "yUnits": "" , "yMin":  0.6, "yMax": +1.1, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio":  0.5, "yMaxRatio": 1.6 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


EffPhi = {
    "xLabel": "#phi"              , "xUnits": "rads", "xMin": -3.2, "xMax": +3.2, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.2f", "yUnits": ""    , "yMin": 0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2 , "yMaxRatio": 1.8 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


Effz0 = {
    "xLabel": "z_{0}"              , "xUnits": "cm", "xMin": -25.0 , "xMax": +25.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.2f" , "yUnits": ""  , "yMin":  +0.00, "yMax":  +1.2, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False ,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.4, "yMaxRatio": 1.6, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


Effd0 = {
    "xLabel": "d_{0}"             , "xUnits": "#mum", "xMin": -10  , "xMax": +10 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Efficiency / %0.2f", "yUnits": ""    , "yMin": +0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2, "yMaxRatio": 1.6, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.16, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


############################################################### 
### Group Histogram (of same binning)
############################################################### 
### TPs
tp_pt     = m_histos.TH1orTH2( "", "tp_pt"    , "TP"      , None, **Pt)
tp_pt_L   = m_histos.TH1orTH2( "", "tp_pt_L"  , "TP, " + L, None, **PtL)
tp_pt_C   = m_histos.TH1orTH2( "", "tp_pt_C"  , "TP, " + C, None, **Pt)
tp_pt_I   = m_histos.TH1orTH2( "", "tp_pt_I"  , "TP, " + I, None, **Pt)
tp_pt_F   = m_histos.TH1orTH2( "", "tp_pt_F"  , "TP, " + F, None, **Pt)
tp_eta    = m_histos.TH1orTH2( "", "tp_eta"   , "TP"      , None, **Eta)
tp_eta_L  = m_histos.TH1orTH2( "", "tp_eta_L" , "TP, " + L, None, **Eta)
tp_eta_M  = m_histos.TH1orTH2( "", "tp_eta_M" , "TP, " + M, None, **Eta)
tp_eta_H  = m_histos.TH1orTH2( "", "tp_eta_H" , "TP, " + H, None, **Eta)
tp_phi    = m_histos.TH1orTH2( "", "tp_phi"   , "TP"      , None, **Phi)
tp_z0     = m_histos.TH1orTH2( "", "tp_z0"    , "TP"      , None, **z0)
tp_d0     = m_histos.TH1orTH2( "", "tp_d0"    , "TP"      , None, **d0)

### TPs (TTTrack-Matched)
match_tp_pt    = m_histos.TH1orTH2( "", "match_tp_pt"   , "Matched TP"      , None, **Pt )
match_tp_pt_L  = m_histos.TH1orTH2( "", "match_tp_pt_L" , "Matched TP, " + L, None, **PtL)
match_tp_pt_C  = m_histos.TH1orTH2( "", "match_tp_pt_C" , "Matched TP, " + C, None, **Pt )
match_tp_pt_I  = m_histos.TH1orTH2( "", "match_tp_pt_I" , "Matched TP, " + I, None, **Pt )
match_tp_pt_F  = m_histos.TH1orTH2( "", "match_tp_pt_F" , "Matched TP, " + F, None, **Pt )
match_tp_eta   = m_histos.TH1orTH2( "", "match_tp_eta"  , "Matched TP"      , None, **Eta)
match_tp_eta_L = m_histos.TH1orTH2( "", "match_tp_eta_L", "Matched TP, " + L, None, **Eta)
match_tp_eta_M = m_histos.TH1orTH2( "", "match_tp_eta_M", "Matched TP, " + M, None, **Eta)
match_tp_eta_H = m_histos.TH1orTH2( "", "match_tp_eta_H", "Matched TP, " + H, None, **Eta)
match_tp_phi   = m_histos.TH1orTH2( "", "match_tp_phi"  , "Matched TP"      , None, **Phi)
match_tp_z0    = m_histos.TH1orTH2( "", "match_tp_z0"   , "Matched TP"      , None, **z0 ) 
match_tp_d0    = m_histos.TH1orTH2( "", "match_tp_d0"   , "Matched TP"      , None, **d0 )

### Efficiencies
eff_pt    = m_histos.TH1orTH2( "", "eff_pt"   , "overall", None, **EffPt )
eff_pt_L  = m_histos.TH1orTH2( "", "eff_pt_L" , "" + L   , None, **EffPtL)
eff_pt_C  = m_histos.TH1orTH2( "", "eff_pt_C" , "" + C   , None, **EffPt)
eff_pt_I  = m_histos.TH1orTH2( "", "eff_pt_I" , "" + I   , None, **EffPt)
eff_pt_F  = m_histos.TH1orTH2( "", "eff_pt_F" , "" + F   , None, **EffPt)
eff_eta   = m_histos.TH1orTH2( "", "eff_eta"  , "overall", None, **EffEta)
eff_eta_L = m_histos.TH1orTH2( "", "eff_eta_L", "" + L   , None, **EffEta)
eff_eta_M = m_histos.TH1orTH2( "", "eff_eta_M", "" + M   , None, **EffEta)
eff_eta_H = m_histos.TH1orTH2( "", "eff_eta_H", "" + H   , None, **EffEta)
eff_phi   = m_histos.TH1orTH2( "", "eff_phi"  , "overall", None, **EffPhi)
eff_z0    = m_histos.TH1orTH2( "", "eff_z0"   , "overall", None, **Effz0 )
eff_d0    = m_histos.TH1orTH2( "", "eff_d0"   , "overall", None, **Effd0 )

### TTTrack (TP-Matched)
label = "" #"Tk, "
match_trk_nstub   = m_histos.TH1orTH2( "", "match_trk_nstub"  , "inclusive"    , None, **NStubs)
match_trk_nstub_C = m_histos.TH1orTH2( "", "match_trk_nstub_C", label + C, None, **NStubs)
match_trk_nstub_I = m_histos.TH1orTH2( "", "match_trk_nstub_I", label + I, None, **NStubs)
match_trk_nstub_F = m_histos.TH1orTH2( "", "match_trk_nstub_F", label + F, None, **NStubs)

### TTTrack (TP-Matched): ChiSq histograms (last bin is an overflow bin)
match_trk_chi2     = m_histos.TH1orTH2( "", "match_trk_chi2"    , "inclusive"           , None, **ChiSq)
match_trk_chi2_L   = m_histos.TH1orTH2( "", "match_trk_chi2_L"  , label + L             , None, **ChiSq)
match_trk_chi2_M   = m_histos.TH1orTH2( "", "match_trk_chi2_M"  , label + M             , None, **ChiSq)
match_trk_chi2_H   = m_histos.TH1orTH2( "", "match_trk_chi2_H"  , label + H             , None, **ChiSq)
match_trk_chi2_C_L = m_histos.TH1orTH2( "", "match_trk_chi2_C_L", label + C + " , "  + L, None, **ChiSq)
match_trk_chi2_I_L = m_histos.TH1orTH2( "", "match_trk_chi2_I_L", label + I + " , "  + L, None, **ChiSq)
match_trk_chi2_F_L = m_histos.TH1orTH2( "", "match_trk_chi2_F_L", label + F + " , "  + L, None, **ChiSq)
match_trk_chi2_C_M = m_histos.TH1orTH2( "", "match_trk_chi2_C_M", label + C + " , "  + M, None, **ChiSq)
match_trk_chi2_I_M = m_histos.TH1orTH2( "", "match_trk_chi2_I_M", label + I + " , "  + M, None, **ChiSq)
match_trk_chi2_F_M = m_histos.TH1orTH2( "", "match_trk_chi2_F_M", label + F + " , "  + M, None, **ChiSq)
match_trk_chi2_C_H = m_histos.TH1orTH2( "", "match_trk_chi2_C_H", label + C + " , "  + H, None, **ChiSq)
match_trk_chi2_I_H = m_histos.TH1orTH2( "", "match_trk_chi2_I_H", label + I + " , "  + H, None, **ChiSq)
match_trk_chi2_F_H = m_histos.TH1orTH2( "", "match_trk_chi2_F_H", label + F + " , "  + H, None, **ChiSq)

### TTTrack (TP-Matched): RedChiSq histograms (lastbin is an overflow bin)
match_trk_chi2_dof     = m_histos.TH1orTH2( "", "match_trk_chi2_dof"    , "inclusive"           , None, **RedChiSq)
match_trk_chi2_dof_L   = m_histos.TH1orTH2( "", "match_trk_chi2_dof_L"  , label + L             , None, **RedChiSq)
match_trk_chi2_dof_M   = m_histos.TH1orTH2( "", "match_trk_chi2_dof_M"  , label + M             , None, **RedChiSq)
match_trk_chi2_dof_H   = m_histos.TH1orTH2( "", "match_trk_chi2_dof_H"  , label + H             , None, **RedChiSq)
match_trk_chi2_dof_C_L = m_histos.TH1orTH2( "", "match_trk_chi2_dof_C_L", label + C + " , "  + L, None, **RedChiSq)
match_trk_chi2_dof_I_L = m_histos.TH1orTH2( "", "match_trk_chi2_dof_I_L", label + I + " , "  + L, None, **RedChiSq)
match_trk_chi2_dof_F_L = m_histos.TH1orTH2( "", "match_trk_chi2_dof_F_L", label + F + " , "  + L, None, **RedChiSq)
match_trk_chi2_dof_C_M = m_histos.TH1orTH2( "", "match_trk_chi2_dof_C_M", label + C + " , "  + M, None, **RedChiSq)
match_trk_chi2_dof_I_M = m_histos.TH1orTH2( "", "match_trk_chi2_dof_I_M", label + I + " , "  + M, None, **RedChiSq)
match_trk_chi2_dof_F_M = m_histos.TH1orTH2( "", "match_trk_chi2_dof_F_M", label + F + " , "  + M, None, **RedChiSq)
match_trk_chi2_dof_C_H = m_histos.TH1orTH2( "", "match_trk_chi2_dof_C_H", label + C + " , "  + H, None, **RedChiSq)
match_trk_chi2_dof_I_H = m_histos.TH1orTH2( "", "match_trk_chi2_dof_I_H", label + I + " , "  + H, None, **RedChiSq)
match_trk_chi2_dof_F_H = m_histos.TH1orTH2( "", "match_trk_chi2_dof_F_H", label + F + " , "  + H, None, **RedChiSq)

### TTTrack (TP-Matched): Resolution histograms
res_pt      = m_histos.TH1orTH2( "", "res_pt"     , "inclusive", None, **PtRes   )
res_pt_C    = m_histos.TH1orTH2( "", "res_pt_C"   , label + C  , None, **PtResAlt)
res_pt_I    = m_histos.TH1orTH2( "", "res_pt_I"   , label + I  , None, **PtResAlt)
res_pt_F    = m_histos.TH1orTH2( "", "res_pt_F"   , label + F  , None, **PtResAlt)

res_ptRel   = m_histos.TH1orTH2( "", "res_ptRel"  , "inclusive", None, **PtResRel)
res_ptRel_C = m_histos.TH1orTH2( "", "res_ptRel_C", label + C  , None, **PtResRel)
res_ptRel_I = m_histos.TH1orTH2( "", "res_ptRel_I", label + I  , None, **PtResRel)
res_ptRel_F = m_histos.TH1orTH2( "", "res_ptRel_F", label + F  , None, **PtResRel)

res_eta     = m_histos.TH1orTH2( "", "res_eta"    , "inclusive", None, **EtaRes  )
res_eta_C   = m_histos.TH1orTH2( "", "res_eta_C"  , label + C  , None, **EtaRes  )
res_eta_I   = m_histos.TH1orTH2( "", "res_eta_I"  , label + I  , None, **EtaRes  )
res_eta_F   = m_histos.TH1orTH2( "", "res_eta_F"  , label + F  , None, **EtaRes  )

res_phi     = m_histos.TH1orTH2( "", "res_phi"    , "inclusive", None, **PhiRes  )
res_phi_C   = m_histos.TH1orTH2( "", "res_phi_C"  , label + C  , None, **PhiRes  )
res_phi_I   = m_histos.TH1orTH2( "", "res_phi_I"  , label + I  , None, **PhiRes  )
res_phi_F   = m_histos.TH1orTH2( "", "res_phi_F"  , label + F  , None, **PhiRes  )

res_z0      = m_histos.TH1orTH2( "", "res_z0"     , "inclusive", None, **Z0Res   )
res_z0_C    = m_histos.TH1orTH2( "", "res_z0_C"   , label + C  , None, **Z0Res   )
res_z0_I    = m_histos.TH1orTH2( "", "res_z0_I"   , label + I  , None, **Z0Res   )
res_z0_F    = m_histos.TH1orTH2( "", "res_z0_F"   , label + F  , None, **Z0Res   )

res_d0      = m_histos.TH1orTH2( "", "res_d0"     , "inclusive", None, **d0Res   )
res_d0_C    = m_histos.TH1orTH2( "", "res_d0_C"   , label + C  , None, **d0Res   )
res_d0_I    = m_histos.TH1orTH2( "", "res_d0_I"   , label + I  , None, **d0Res   )
res_d0_F    = m_histos.TH1orTH2( "", "res_d0_F"   , label + F  , None, **d0Res   )

### Resolution vs. pt histograms  
ptRange = ["0-5","5-10", "10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50"] #"50-55", "55-60","60-65","65-70","70-75","75-80","80-85","85-90","90-95","95-100"]
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
    resVsPt_pt      = m_histos.TH1orTH2( "", "resVsPt_pt_"      + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **PtRes)
    resVsPt_pt_C    = m_histos.TH1orTH2( "", "resVsPt_pt_C_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **PtRes)
    resVsPt_pt_I    = m_histos.TH1orTH2( "", "resVsPt_pt_I_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **PtRes)
    resVsPt_pt_F    = m_histos.TH1orTH2( "", "resVsPt_pt_F_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **PtRes)

    resVsPt_ptRel   = m_histos.TH1orTH2( "", "resVsPt_ptRel_"   + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **PtResRel)
    resVsPt_ptRel_C = m_histos.TH1orTH2( "", "resVsPt_ptRel_C_" + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **PtResRel)
    resVsPt_ptRel_I = m_histos.TH1orTH2( "", "resVsPt_ptRel_I_" + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **PtResRel)
    resVsPt_ptRel_F = m_histos.TH1orTH2( "", "resVsPt_ptRel_F_" + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **PtResRel)

    resVsPt_eta     = m_histos.TH1orTH2( "", "resVsPt_eta_"     + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **EtaRes)
    resVsPt_eta_C   = m_histos.TH1orTH2( "", "resVsPt_eta_C_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **EtaRes)
    resVsPt_eta_I   = m_histos.TH1orTH2( "", "resVsPt_eta_I_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **EtaRes)
    resVsPt_eta_F   = m_histos.TH1orTH2( "", "resVsPt_eta_F_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **EtaRes)

    resVsPt_phi     = m_histos.TH1orTH2( "", "resVsPt_phi_"     + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **PhiRes)
    resVsPt_phi_C   = m_histos.TH1orTH2( "", "resVsPt_phi_C_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **PhiRes)
    resVsPt_phi_I   = m_histos.TH1orTH2( "", "resVsPt_phi_I_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **PhiRes)
    resVsPt_phi_F   = m_histos.TH1orTH2( "", "resVsPt_phi_F_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **PhiRes)

    resVsPt_z0      = m_histos.TH1orTH2( "", "resVsPt_z0_"      + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **Z0Res)
    resVsPt_z0_C    = m_histos.TH1orTH2( "", "resVsPt_z0_C_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **Z0Res)
    resVsPt_z0_I    = m_histos.TH1orTH2( "", "resVsPt_z0_I_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **Z0Res)
    resVsPt_z0_F    = m_histos.TH1orTH2( "", "resVsPt_z0_F_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **Z0Res)

    resVsPt_d0      = m_histos.TH1orTH2( "", "resVsPt_d0_"      + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **d0Res)
    resVsPt_d0_C    = m_histos.TH1orTH2( "", "resVsPt_d0_C_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **d0Res)
    resVsPt_d0_I    = m_histos.TH1orTH2( "", "resVsPt_d0_I_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **d0Res)
    resVsPt_d0_F    = m_histos.TH1orTH2( "", "resVsPt_d0_F_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **d0Res)


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
etaRange     = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2"]
etaRangePlus = ["0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4"]
h_resVsEta_pt      = []
h_resVsEta_pt_L    = []
h_resVsEta_pt_M    = []
h_resVsEta_pt_H    = []

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
    resVsEta_pt      = m_histos.TH1orTH2( "", "h_resVsEta_pt_"   + etaRangePlus[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **PtRes)
    resVsEta_pt_L    = m_histos.TH1orTH2( "", "h_resVsEta_pt_L_" + etaRangePlus[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtRes)
    resVsEta_pt_M    = m_histos.TH1orTH2( "", "h_resVsEta_pt_M_" + etaRangePlus[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtRes)
    resVsEta_pt_H    = m_histos.TH1orTH2( "", "h_resVsEta_pt_H_" + etaRangePlus[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtRes)

    resVsEta_ptRel   = m_histos.TH1orTH2( "", "h_resVsEta_ptRel_"   + etaRangePlus[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **PtResRel)
    resVsEta_ptRel_L = m_histos.TH1orTH2( "", "h_resVsEta_ptRel_L_" + etaRangePlus[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtResRel)
    resVsEta_ptRel_M = m_histos.TH1orTH2( "", "h_resVsEta_ptRel_M_" + etaRangePlus[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtResRel)
    resVsEta_ptRel_H = m_histos.TH1orTH2( "", "h_resVsEta_ptRel_H_" + etaRangePlus[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PtResRel)

    resVsEta_eta     = m_histos.TH1orTH2( "", "h_resVsEta_eta_"     + etaRangePlus[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **EtaRes)
    resVsEta_eta_L   = m_histos.TH1orTH2( "", "h_resVsEta_eta_L_"   + etaRangePlus[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **EtaRes)
    resVsEta_eta_M   = m_histos.TH1orTH2( "", "h_resVsEta_eta_M_"   + etaRangePlus[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **EtaRes)
    resVsEta_eta_H   = m_histos.TH1orTH2( "", "h_resVsEta_eta_H_"   + etaRangePlus[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **EtaRes)

    resVsEta_phi     = m_histos.TH1orTH2( "", "h_resVsEta_phi_"     + etaRangePlus[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **PhiRes)
    resVsEta_phi_L   = m_histos.TH1orTH2( "", "h_resVsEta_phi_L_"   + etaRangePlus[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PhiRes)
    resVsEta_phi_M   = m_histos.TH1orTH2( "", "h_resVsEta_phi_M_"   + etaRangePlus[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PhiRes)
    resVsEta_phi_H   = m_histos.TH1orTH2( "", "h_resVsEta_phi_H_"   + etaRangePlus[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **PhiRes)

    resVsEta_z0      = m_histos.TH1orTH2( "", "h_resVsEta_z0_"   + etaRangePlus[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **Z0Res)
    resVsEta_z0_L    = m_histos.TH1orTH2( "", "h_resVsEta_z0_L_" + etaRangePlus[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **Z0Res)
    resVsEta_z0_M    = m_histos.TH1orTH2( "", "h_resVsEta_z0_M_" + etaRangePlus[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **Z0Res)
    resVsEta_z0_H    = m_histos.TH1orTH2( "", "h_resVsEta_z0_H_" + etaRangePlus[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **Z0Res)

    resVsEta_d0      = m_histos.TH1orTH2( "", "h_resVsEta_d0_"   + etaRangePlus[i], etaRange[i] + " < |#eta| < " + etaRangePlus[i]              , None, **d0Res)
    resVsEta_d0_L    = m_histos.TH1orTH2( "", "h_resVsEta_d0_L_" + etaRangePlus[i], L  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **d0Res)
    resVsEta_d0_M    = m_histos.TH1orTH2( "", "h_resVsEta_d0_M_" + etaRangePlus[i], M  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **d0Res)
    resVsEta_d0_H    = m_histos.TH1orTH2( "", "h_resVsEta_d0_H_" + etaRangePlus[i], H  + " , " +  etaRange[i] + " < |#eta| < " + etaRangePlus[i], None, **d0Res)

    ### Append to lists
    h_resVsEta_pt.append(resVsEta_pt)
    h_resVsEta_pt_L.append(resVsEta_pt_L)
    h_resVsEta_pt_M.append(resVsEta_pt_M)
    h_resVsEta_pt_H.append(resVsEta_pt_H)

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
logChiSq_eta     = m_histos.TH1orTH2( "", "2d_logchi2_eta"    , "", None, **ChiSqEta)
logChiSq_dof_eta = m_histos.TH1orTH2( "", "2d_logchi2_dof_eta", "", None, **RedChiSqEta)
dz0_eta          = m_histos.TH1orTH2( "", "2d_dz0_eta"        , "", None, **dZ0Eta)
dd0_eta          = m_histos.TH1orTH2( "", "2d_dd0_eta"        , "", None, **dd0Eta)
deta_eta         = m_histos.TH1orTH2( "", "2d_deta_eta"       , "", None, **dEtaEta)
dphi_eta         = m_histos.TH1orTH2( "", "2d_dphi_eta"       , "", None, **dPhiEta)
dpt_eta          = m_histos.TH1orTH2( "", "2d_dpt_eta"        , "", None, **dPtEta)
dptRel_eta       = m_histos.TH1orTH2( "", "2d_dptRel_eta"     , "", None, **dPtRelEta)


### resolution vs. pT histograms. NB: Filled with RMS and RMSError
resVsPt_pt   = m_histos.TH1orTH2( "", "resVsPt_pt"  , "inclusive", None, **PtPtRes)
resVsPt_pt_C = m_histos.TH1orTH2( "", "resVsPt_pt_C", C          , None, **PtPtRes)
resVsPt_pt_I = m_histos.TH1orTH2( "", "resVsPt_pt_I", I          , None, **PtPtRes)
resVsPt_pt_F = m_histos.TH1orTH2( "", "resVsPt_pt_F", F          , None, **PtPtRes)

resVsPt_ptRel   = m_histos.TH1orTH2( "",  "resVsPt_ptRel"  , "inclusive", None, **PtPtResRel)
resVsPt_ptRel_C = m_histos.TH1orTH2( "",  "resVsPt_ptRel_C",  C         , None, **PtPtResRel)
resVsPt_ptRel_I = m_histos.TH1orTH2( "",  "resVsPt_ptRel_I",  I         , None, **PtPtResRel)
resVsPt_ptRel_F = m_histos.TH1orTH2( "",  "resVsPt_ptRel_F",  F         , None, **PtPtResRel)

resVsPt_eta   = m_histos.TH1orTH2( "", "resVsPt_eta"  , "inclusive", None, **PtEtaRes)
resVsPt_eta_C = m_histos.TH1orTH2( "", "resVsPt_eta_C",  C         , None, **PtEtaRes)
resVsPt_eta_I = m_histos.TH1orTH2( "", "resVsPt_eta_I",  I         , None, **PtEtaRes)
resVsPt_eta_F = m_histos.TH1orTH2( "", "resVsPt_eta_F",  F         , None, **PtEtaRes)

resVsPt_phi   = m_histos.TH1orTH2( "", "resVsPt_phi"  ,  "inclusive", None, **PtPhiRes)
resVsPt_phi_C = m_histos.TH1orTH2( "", "resVsPt_phi_C",  C          , None, **PtPhiRes)
resVsPt_phi_I = m_histos.TH1orTH2( "", "resVsPt_phi_I",  I          , None, **PtPhiRes)
resVsPt_phi_F = m_histos.TH1orTH2( "", "resVsPt_phi_F",  F          , None, **PtPhiRes)

resVsPt_z0    = m_histos.TH1orTH2( "", "resVsPt_z0"  ,  "inclusive", None, **PtZ0Res)
resVsPt_z0_C  = m_histos.TH1orTH2( "", "resVsPt_z0_C",  C          , None, **PtZ0Res)
resVsPt_z0_I  = m_histos.TH1orTH2( "", "resVsPt_z0_I",  I          , None, **PtZ0Res)
resVsPt_z0_F  = m_histos.TH1orTH2( "", "resVsPt_z0_F",  F          , None, **PtZ0Res)

resVsPt_d0    = m_histos.TH1orTH2( "", "resVsPt_d0"  ,  "inclusive", None, **Ptd0Res)
resVsPt_d0_C  = m_histos.TH1orTH2( "", "resVsPt_d0_C",  C          , None, **Ptd0Res)
resVsPt_d0_I  = m_histos.TH1orTH2( "", "resVsPt_d0_I",  I          , None, **Ptd0Res)
resVsPt_d0_F  = m_histos.TH1orTH2( "", "resVsPt_d0_F",  F          , None, **Ptd0Res)


### resolution vs. eta histograms
resVsEta_pt   = m_histos.TH1orTH2( "", "resVsEta_pt"  , "inclusive", None, **EtaPtRes)
resVsEta_pt_L = m_histos.TH1orTH2( "", "resVsEta_pt_L", L          , None, **EtaPtRes)
resVsEta_pt_M = m_histos.TH1orTH2( "", "resVsEta_pt_M", M          , None, **EtaPtRes)
resVsEta_pt_H = m_histos.TH1orTH2( "", "resVsEta_pt_H", H          , None, **EtaPtRes)

resVsEta_ptRel   = m_histos.TH1orTH2( "", "resVsEta_ptRel"  , "inclusive", None, **EtaPtResRel)
resVsEta_ptRel_L = m_histos.TH1orTH2( "", "resVsEta_ptRel_L", L          , None, **EtaPtResRel)
resVsEta_ptRel_M = m_histos.TH1orTH2( "", "resVsEta_ptRel_M", M          , None, **EtaPtResRel)
resVsEta_ptRel_H = m_histos.TH1orTH2( "", "resVsEta_ptRel_H", H          , None, **EtaPtResRel)

resVsEta_eta   = m_histos.TH1orTH2( "", "resVsEta_eta"  , "inclusive", None, **EtaEtaRes)
resVsEta_eta_L = m_histos.TH1orTH2( "", "resVsEta_eta_L", L          , None, **EtaEtaRes)
resVsEta_eta_M = m_histos.TH1orTH2( "", "resVsEta_eta_M", M          , None, **EtaEtaRes)
resVsEta_eta_H = m_histos.TH1orTH2( "", "resVsEta_eta_H", H          , None, **EtaEtaRes)

resVsEta_phi   = m_histos.TH1orTH2( "", "resVsEta_phi"  , "inclusive", None, **EtaPhiRes)
resVsEta_phi_L = m_histos.TH1orTH2( "", "resVsEta_phi_L", L          , None, **EtaPhiRes)
resVsEta_phi_M = m_histos.TH1orTH2( "", "resVsEta_phi_M", M          , None, **EtaPhiRes)
resVsEta_phi_H = m_histos.TH1orTH2( "", "resVsEta_phi_H", H          , None, **EtaPhiRes)

resVsEta_z0   = m_histos.TH1orTH2( "", "resVsEta_z0"  , "inclusive", None, **EtaZ0Res)
resVsEta_z0_L = m_histos.TH1orTH2( "", "resVsEta_z0_L", L          , None, **EtaZ0Res)
resVsEta_z0_M = m_histos.TH1orTH2( "", "resVsEta_z0_M", M          , None, **EtaZ0Res)
resVsEta_z0_H = m_histos.TH1orTH2( "", "resVsEta_z0_H", H          , None, **EtaZ0Res)

resVsEta_d0   = m_histos.TH1orTH2( "", "resVsEta_d0"  , "inclusive", None, **Etad0Res)
resVsEta_d0_L = m_histos.TH1orTH2( "", "resVsEta_d0_L", L          , None, **Etad0Res)
resVsEta_d0_M = m_histos.TH1orTH2( "", "resVsEta_d0_M", M          , None, **Etad0Res)
resVsEta_d0_H = m_histos.TH1orTH2( "", "resVsEta_d0_H", H          , None, **Etad0Res)


###############################################################
### Main
###############################################################
def DoPlots(hList, dataset, bColourPalette=False, saveExt=""):
    
    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.SetBoolUseDatasetAsLegEntry(bColourPalette)

    p.AddDataset(dataset + ":1x-BE"         , datasetPaths_1p0x[dataset]  )
    ### p.AddDataset(dataset + ":0.25x-BE"      , datasetPaths_0p25x[dataset] )
    p.AddDataset(dataset + ":0.5x-BE"       , datasetPaths_0p5x[dataset]  )
    p.AddDataset(dataset + ":1.5x-BE "      , datasetPaths_1p5x[dataset]  )
    p.AddDataset(dataset + ":2x-BE"         , datasetPaths_2p0x[dataset]  )
    ### p.AddDataset(dataset + ":4x-BE "        , datasetPaths_4p0x[dataset]  )
    # p.AddDataset(dataset + ":0.5x-B, 2x-E", datasetPaths_0p5xB2p0E[dataset]  )
    
    if (len(p.DatasetToRootFileMap.keys()) > 1):
        p.SetBoolUseDatasetAsLegEntry(True)
    else:
        p.SetBoolUseDatasetAsLegEntry(False)
        
    p.EnableColourPalette(bColourPalette)
    p.AddHisto(hList)
    p.SetupStatsBox()
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = True)
    p.SetTLegendHeader(dataset, "" )
    #p.AppendToCanvasName("_" + saveNameExt)
    p.SaveHistos(True, savePath, saveFormats, saveExt)

    return

###############################################################
if __name__ == "__main__":
    bColourPalette = True
    
    if(bTPs):
        DoPlots( [tp_pt]   , datasetList[0], bColourPalette)
        DoPlots( [tp_pt_L] , datasetList[0], bColourPalette )
        DoPlots( [tp_pt_C] , datasetList[0], bColourPalette )
        DoPlots( [tp_pt_I] , datasetList[0], bColourPalette )
        DoPlots( [tp_pt_F] , datasetList[0], bColourPalette )
        DoPlots( [tp_eta]  , datasetList[0], bColourPalette )
        DoPlots( [tp_eta_L], datasetList[0], bColourPalette )
        DoPlots( [tp_eta_M], datasetList[0], bColourPalette )
        DoPlots( [tp_eta_H], datasetList[0], bColourPalette )
        DoPlots( [tp_phi]  , datasetList[0], bColourPalette )
        DoPlots( [tp_z0]   , datasetList[0], bColourPalette )
        DoPlots( [tp_d0]   , datasetList[0], bColourPalette )

        
    if(bMatchedTPs):
        DoPlots( [match_tp_pt   ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_pt_L ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_pt_C ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_pt_I ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_pt_F ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_eta  ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_eta_L], datasetList[0], bColourPalette)
        DoPlots( [match_tp_eta_M], datasetList[0], bColourPalette)
        DoPlots( [match_tp_eta_H], datasetList[0], bColourPalette)
        DoPlots( [match_tp_phi  ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_z0   ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_d0   ], datasetList[0], bColourPalette)


    if(bMatchedTks):
        DoPlots( [match_trk_nstub]  , datasetList[0], bColourPalette)        
        DoPlots( [match_trk_chi2_L    , match_trk_chi2_M    , match_trk_chi2_H]    , datasetList[0], bColourPalette, "MH")
        DoPlots( [match_trk_chi2_dof_L, match_trk_chi2_dof_M, match_trk_chi2_dof_H], datasetList[0], bColourPalette, "MH")
        DoPlots( [match_trk_nstub_C   , match_trk_nstub_I   , match_trk_nstub_F]   , datasetList[0], bColourPalette, "IF")
        if(bMatchedTks_EtaRange):
            DoPlots( [match_trk_nstub_C], datasetList[0], bColourPalette)
            DoPlots( [match_trk_nstub_I], datasetList[0], bColourPalette)
            DoPlots( [match_trk_nstub_F], datasetList[0], bColourPalette)

        
    if(bResolutions):
        DoPlots( [res_pt   ], datasetList[0], bColourPalette)
        DoPlots( [res_ptRel], datasetList[0], bColourPalette)
        DoPlots( [res_eta  ], datasetList[0], bColourPalette)
        DoPlots( [res_phi  ], datasetList[0], bColourPalette)
        DoPlots( [res_z0   ], datasetList[0], bColourPalette)
        DoPlots( [res_d0   ], datasetList[0], bColourPalette)

        
    if(bResiduals_EtaRange):
        DoPlots( [res_pt_C ]   , datasetList[0], bColourPalette)
        DoPlots( [res_pt_I ]   , datasetList[0], bColourPalette)
        DoPlots( [res_pt_F ]   , datasetList[0], bColourPalette)
        DoPlots( [res_ptRel_C ], datasetList[0], bColourPalette)
        DoPlots( [res_ptRel_I ], datasetList[0], bColourPalette)
        DoPlots( [res_ptRel_F ], datasetList[0], bColourPalette)
        DoPlots( [res_eta_C ]  , datasetList[0], bColourPalette)
        DoPlots( [res_eta_I ]  , datasetList[0], bColourPalette)
        DoPlots( [res_eta_F ]  , datasetList[0], bColourPalette)
        DoPlots( [res_phi_C ]  , datasetList[0], bColourPalette)
        DoPlots( [res_phi_I ]  , datasetList[0], bColourPalette)
        DoPlots( [res_phi_F ]  , datasetList[0], bColourPalette)
        DoPlots( [res_z0_C ]   , datasetList[0], bColourPalette)
        DoPlots( [res_z0_I ]   , datasetList[0], bColourPalette)
        DoPlots( [res_z0_F ]   , datasetList[0], bColourPalette)
        DoPlots( [res_d0_C ]   , datasetList[0], bColourPalette)
        DoPlots( [res_d0_I ]   , datasetList[0], bColourPalette)
        DoPlots( [res_d0_F ]   , datasetList[0], bColourPalette)
            

    if(bResiduals_CIF):        
        DoPlots( [res_pt_C   , res_pt_I   , res_pt_F    ], datasetList[0], bColourPalette, "IF")
        DoPlots( [res_ptRel_C, res_ptRel_I, res_ptRel_F ], datasetList[0], bColourPalette, "IF")
        DoPlots( [res_eta_C  , res_eta_I  , res_eta_F   ], datasetList[0], bColourPalette, "IF")
        DoPlots( [res_phi_C  , res_phi_I  , res_phi_F   ], datasetList[0], bColourPalette, "IF")
        DoPlots( [res_z0_C   , res_z0_I   , res_z0_F    ], datasetList[0], bColourPalette, "IF")
        DoPlots( [res_d0_C   , res_d0_I   , res_d0_F    ], datasetList[0], bColourPalette, "IF")
        
        
    if(bTkChiSq):
        DoPlots( [match_trk_chi2], datasetList[0], bColourPalette)
        DoPlots( [match_trk_chi2_C_L, match_trk_chi2_I_L, match_trk_chi2_F_L], datasetList[0], True, "_IL_FL")
        DoPlots( [match_trk_chi2_C_M, match_trk_chi2_I_M, match_trk_chi2_F_M], datasetList[0], True, "_IM_FM")
    if(bTkChiSq_EtaRange):
       DoPlots( [match_trk_chi2_C_L], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_I_L], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_F_L], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_C_M], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_I_M], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_F_M], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_C_H], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_I_H], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_F_H], datasetList[0], bColourPalette)
    if(bTkRedChiSq):        
       DoPlots( [match_trk_chi2_dof], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_C_L, match_trk_chi2_dof_I_L, match_trk_chi2_dof_F_L], datasetList[0], True, "_IL_FL")
       DoPlots( [match_trk_chi2_dof_C_M, match_trk_chi2_dof_I_M, match_trk_chi2_dof_F_M], datasetList[0], True, "_ML_FM")
    if(bTkRedChiSq_EtaRange):        
       DoPlots( [match_trk_chi2_dof_C_L], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_I_L], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_F_L], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_C_M], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_I_M], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_F_M], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_C_H], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_I_H], datasetList[0], bColourPalette)
       DoPlots( [match_trk_chi2_dof_F_H], datasetList[0], bColourPalette)


    if(bTkEfficiency):
       DoPlots( [eff_pt]   , datasetList[0], bColourPalette)
       DoPlots( [eff_pt_L] , datasetList[0], bColourPalette)
       DoPlots( [eff_eta ] , datasetList[0], bColourPalette)
    if (bTkEfficiency_Vs_Pt):
        DoPlots( [eff_pt_C] , datasetList[0], bColourPalette)
        DoPlots( [eff_pt_I] , datasetList[0], bColourPalette)
        DoPlots( [eff_pt_F] , datasetList[0], bColourPalette)
        DoPlots( [eff_pt    , eff_pt_C, eff_pt_I, eff_pt_F]   , datasetList[0], True, "_CIF")
    if (bTkEfficiency_Vs_Eta):
        DoPlots( [eff_eta_L], datasetList[0], bColourPalette)
        DoPlots( [eff_eta_M], datasetList[0], bColourPalette)
        DoPlots( [eff_eta_H], datasetList[0], bColourPalette)
        DoPlots( [eff_eta   , eff_eta_L, eff_eta_M, eff_eta_H], datasetList[0], True, "_LMH")
    if (bTkEfficiency_Vs_d0z0phi):
        DoPlots( [eff_phi]  , datasetList[0], bColourPalette)
        DoPlots( [eff_z0 ]  , datasetList[0], bColourPalette)
        DoPlots( [eff_d0]   , datasetList[0], bColourPalette)
                            

    if(bDo2D):
        DoPlots( [ logChiSq_eta ]     , datasetList[0], bColourPalette)
        DoPlots( [ logChiSq_dof_eta ] , datasetList[0], bColourPalette)
        DoPlots( [ dz0_eta ]          , datasetList[0], bColourPalette)
        DoPlots( [ deta_eta ]         , datasetList[0], bColourPalette)
        DoPlots( [ dphi_eta ]         , datasetList[0], bColourPalette)
        DoPlots( [ dpt_eta ]          , datasetList[0], bColourPalette)
        DoPlots( [ dptRel_eta ]       , datasetList[0], bColourPalette)


    ### pt residuals
    if(pt_Resolution_Vs_Pt):
        DoPlots( [ resVsPt_pt  ], datasetList[0], bColourPalette)
    if(pt_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resVsPt_pt_C], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_pt_I], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_pt_F], datasetList[0], bColourPalette)
    if(pt_Resolution_Vs_Pt_CIF):
        DoPlots( [ resVsPt_pt_C, resVsPt_pt_I, resVsPt_pt_F], datasetList[0], True, "IF")

    if(pt_Residual_Vs_Pt_PtRange or pt_Residual_Vs_Pt_EtaRange or pt_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(pt_Residual_Vs_Pt_PtRange):
                DoPlots( h_resVsPt_pt[i], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_resVsPt_pt_C[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_pt_I[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_pt_F[i]], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_resVsPt_pt_C[i], h_resVsPt_pt_I[i], h_resVsPt_pt_F[i]], datasetList[0], bColourPalette)
                
    if(pt_Resolution_Vs_Eta):
        DoPlots( [ resVsEta_pt  ], datasetList[0], bColourPalette)
    if (pt_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resVsEta_pt_L ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_pt_M ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_pt_H ], datasetList[0], bColourPalette)
    if (pt_Resolution_Vs_Eta_LMH):
        DoPlots( [ resVsEta_pt_L, resVsEta_pt_M, resVsEta_pt_H ], datasetList[0], True, "MH")

    if(pt_Residual_Vs_Eta_EtaRange or pt_Residual_Vs_Eta_PtRange_EtaRange or pt_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(pt_Residual_Vs_Eta_EtaRange):
                DoPlots( h_resVsEta_pt[i], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_resVsEta_pt_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_pt_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_pt_H[i] ], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_resVsEta_pt_L[i], h_resVsEta_pt_M[i], h_resVsEta_pt_H[i] ], datasetList[0], True, "MH")

                
    ### ptRel residuals
    if(ptRel_Resolution_Vs_Pt):
        DoPlots( [ resVsPt_ptRel  ], datasetList[0], bColourPalette)
    if(ptRel_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resVsPt_ptRel_C], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_ptRel_I], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_ptRel_F], datasetList[0], bColourPalette)
    if(ptRel_Resolution_Vs_Pt_CIF):
        DoPlots( [ resVsPt_ptRel_C, resVsPt_ptRel_I, resVsPt_ptRel_F], datasetList[0], True, "IF")

    if(ptRel_Residual_Vs_Pt_PtRange or ptRel_Residual_Vs_Pt_EtaRange or ptRel_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(ptRel_Residual_Vs_Pt_PtRange):
                DoPlots( h_resVsPt_ptRel[i], datasetList[0], bColourPalette)
            if (ptRel_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_resVsPt_ptRel_C[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_ptRel_I[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_ptRel_F[i]], datasetList[0], bColourPalette)
            if (ptRel_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_resVsPt_ptRel_C[i], h_resVsPt_ptRel_I[i], h_resVsPt_ptRel_F[i]], datasetList[0], bColourPalette)
                
    if(ptRel_Resolution_Vs_Eta):
        DoPlots( [ resVsEta_ptRel  ], datasetList[0], bColourPalette)
    if (ptRel_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resVsEta_ptRel_L ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_ptRel_M ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_ptRel_H ], datasetList[0], bColourPalette)
    if (ptRel_Resolution_Vs_Eta_LMH):
        DoPlots( [ resVsEta_ptRel_L, resVsEta_ptRel_M, resVsEta_ptRel_H ], datasetList[0], True, "MH")

    if(ptRel_Residual_Vs_Eta_EtaRange or ptRel_Residual_Vs_Eta_PtRange_EtaRange or ptRel_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(ptRel_Residual_Vs_Eta_EtaRange):
                DoPlots( h_resVsEta_ptRel[i], datasetList[0], bColourPalette)
            if (ptRel_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_resVsEta_ptRel_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_ptRel_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_ptRel_H[i] ], datasetList[0], bColourPalette)
            if (ptRel_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_resVsEta_ptRel_L[i], h_resVsEta_ptRel_M[i], h_resVsEta_ptRel_H[i] ], datasetList[0], True, "MH")


    ### eta residuals
    if(eta_Resolution_Vs_Pt):
        DoPlots( [ resVsPt_eta  ], datasetList[0], bColourPalette)
    if(eta_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resVsPt_eta_C], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_eta_I], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_eta_F], datasetList[0], bColourPalette)
    if(eta_Resolution_Vs_Pt_CIF):
        DoPlots( [ resVsPt_eta_C, resVsPt_eta_I, resVsPt_eta_F], datasetList[0], True, "IF")

    if(eta_Residual_Vs_Pt_PtRange or eta_Residual_Vs_Pt_EtaRange or eta_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(eta_Residual_Vs_Pt_PtRange):
                DoPlots( h_resVsPt_eta[i], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_resVsPt_eta_C[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_eta_I[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_eta_F[i]], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_resVsPt_eta_C[i], h_resVsPt_eta_I[i], h_resVsPt_eta_F[i]], datasetList[0], bColourPalette)
                
    if(eta_Resolution_Vs_Eta):
        DoPlots( [ resVsEta_eta  ], datasetList[0], bColourPalette)
    if (eta_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resVsEta_eta_L ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_eta_M ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_eta_H ], datasetList[0], bColourPalette)
    if (eta_Resolution_Vs_Eta_LMH):
        DoPlots( [ resVsEta_eta_L, resVsEta_eta_M, resVsEta_eta_H ], datasetList[0], True, "MH")

    if(eta_Residual_Vs_Eta_EtaRange or eta_Residual_Vs_Eta_PtRange_EtaRange or eta_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(eta_Residual_Vs_Eta_EtaRange):
                DoPlots( h_resVsEta_eta[i], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_resVsEta_eta_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_eta_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_eta_H[i] ], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_resVsEta_eta_L[i], h_resVsEta_eta_M[i], h_resVsEta_eta_H[i] ], datasetList[0], True, "MH")

        
    ### phi residuals
    if(phi_Resolution_Vs_Pt):
        DoPlots( [ resVsPt_phi  ], datasetList[0], bColourPalette)
    if(phi_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resVsPt_phi_C], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_phi_I], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_phi_F], datasetList[0], bColourPalette)
    if(phi_Resolution_Vs_Pt_CIF):
        DoPlots( [ resVsPt_phi_C, resVsPt_phi_I, resVsPt_phi_F], datasetList[0], True, "IF")

    if(phi_Residual_Vs_Pt_PtRange or phi_Residual_Vs_Pt_EtaRange or phi_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(phi_Residual_Vs_Pt_PtRange):
                DoPlots( h_resVsPt_phi[i], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_resVsPt_phi_C[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_phi_I[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_phi_F[i]], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_resVsPt_phi_C[i], h_resVsPt_phi_I[i], h_resVsPt_phi_F[i]], datasetList[0], bColourPalette)
                
    if(phi_Resolution_Vs_Eta):
        DoPlots( [ resVsEta_phi  ], datasetList[0], bColourPalette)
    if (phi_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resVsEta_phi_L ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_phi_M ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_phi_H ], datasetList[0], bColourPalette)
    if (phi_Resolution_Vs_Eta_LMH):
        DoPlots( [ resVsEta_phi_L, resVsEta_phi_M, resVsEta_phi_H ], datasetList[0], True, "MH")

    if(phi_Residual_Vs_Eta_EtaRange or phi_Residual_Vs_Eta_PtRange_EtaRange or phi_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(phi_Residual_Vs_Eta_EtaRange):
                DoPlots( h_resVsEta_phi[i], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_resVsEta_phi_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_phi_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_phi_H[i] ], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_resVsEta_phi_L[i], h_resVsEta_phi_M[i], h_resVsEta_phi_H[i] ], datasetList[0], True, "MH")
                
                
    ### z0 residuals
    if(z0_Resolution_Vs_Pt):
        DoPlots( [ resVsPt_z0  ], datasetList[0], bColourPalette)
    if(z0_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resVsPt_z0_C ], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_z0_I ], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_z0_F ], datasetList[0], bColourPalette)
    if(z0_Resolution_Vs_Pt_CIF):
        DoPlots( [ resVsPt_z0_C, resVsPt_z0_I, resVsPt_z0_F], datasetList[0], True, "IF")

    if(z0_Residual_Vs_Pt_PtRange or z0_Residual_Vs_Pt_EtaRange or z0_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(z0_Residual_Vs_Pt_PtRange):
                DoPlots( h_resVsPt_z0[i], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_resVsPt_z0_C[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_z0_I[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_z0_F[i] ], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_resVsPt_z0_C[i], h_resVsPt_z0_I[i], h_resVsPt_z0_F[i]], datasetList[0], bColourPalette)
                
    if(z0_Resolution_Vs_Eta):
        DoPlots( [ resVsEta_z0  ], datasetList[0], bColourPalette)
    if (z0_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resVsEta_z0_L ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_z0_M ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_z0_H ], datasetList[0], bColourPalette)
    if (z0_Resolution_Vs_Eta_LMH):
        DoPlots( [ resVsEta_z0_L, resVsEta_z0_M, resVsEta_z0_H ], datasetList[0], True, "MH")

    if(z0_Residual_Vs_Eta_EtaRange or z0_Residual_Vs_Eta_PtRange_EtaRange or z0_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(z0_Residual_Vs_Eta_EtaRange):
                DoPlots( h_resVsEta_z0[i], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_resVsEta_z0_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_z0_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_z0_H[i] ], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_resVsEta_z0_L[i], h_resVsEta_z0_M[i], h_resVsEta_z0_H[i] ], datasetList[0], True, "MH")

                
    ### d0 residuals
    if(d0_Resolution_Vs_Pt):
        DoPlots( [ resVsPt_d0  ], datasetList[0], bColourPalette)
    if(d0_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resVsPt_d0_C ], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_d0_I ], datasetList[0], bColourPalette)
        DoPlots( [ resVsPt_d0_F ], datasetList[0], bColourPalette)
    if(d0_Resolution_Vs_Pt_CIF):
        DoPlots( [ resVsPt_d0_C, resVsPt_d0_I, resVsPt_d0_F], datasetList[0], True, "IF")

    if(d0_Residual_Vs_Pt_PtRange or d0_Residual_Vs_Pt_EtaRange or d0_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(d0_Residual_Vs_Pt_PtRange):
                DoPlots( h_resVsPt_d0[i], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_resVsPt_d0_C[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_d0_I[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsPt_d0_F[i] ], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_resVsPt_d0_C[i], h_resVsPt_d0_I[i], h_resVsPt_d0_F[i]], datasetList[0], bColourPalette)
                
    if(d0_Resolution_Vs_Eta):
        DoPlots( [ resVsEta_d0  ], datasetList[0], bColourPalette)
    if (d0_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resVsEta_d0_L ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_d0_M ], datasetList[0], bColourPalette)
        DoPlots( [ resVsEta_d0_H ], datasetList[0], bColourPalette)
    if (d0_Resolution_Vs_Eta_LMH):
        DoPlots( [ resVsEta_d0_L, resVsEta_d0_M, resVsEta_d0_H ], datasetList[0], True, "MH")

    if(d0_Residual_Vs_Eta_EtaRange or d0_Residual_Vs_Eta_PtRange_EtaRange or d0_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(d0_Residual_Vs_Eta_EtaRange):
                DoPlots( h_resVsEta_d0[i], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_resVsEta_d0_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_d0_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_resVsEta_d0_H[i] ], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_resVsEta_d0_L[i], h_resVsEta_d0_M[i], h_resVsEta_d0_H[i] ], datasetList[0], True, "MH")


############################################################################################################################################

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
bDo2D                      = False # 1 dataset
bDo_Pixel_FitHits          = False # datasets
bDo_Pixel_FitHits_2D       = False # 1 dataset
bDo_CandPixel_FitHits      = False # datasets
bDo_CandPixel_FitHits_2D   = False # 1 dataset
bDo_Pixel_SharedFitHits    = False
bDo_Pixel_SharedFitHits_2D = False


### Tracks and Tracking Particles
bTPs                     = False
bMatchedTPs              = False
bMatchedTks              = False


### Tracks 
bTkChiSq                 = False # datasets
bTkChiSq_Slices          = False
bTkChiSq_EtaRange        = False
bTkRedChiSq              = False
bTkRedChiSq_Slices       = False 
bTkRedChiSq_EtaRange     = False 


### Inclusive Residuals & Residuals
bResiduals               = False 
bResiduals_EtaRange      = False
bResiduals_CIF           = False
            

### Tracking Efficiency
bTkEfficiency            = False  # datasets
bTkEfficiency_Vs_Pt      = False
bTkEfficiency_Vs_Pt_CIF  = False
bTkEfficiency_Vs_Eta     = False
bTkEfficiency_Vs_Eta_LMH = False
bTkEfficiency_Vs_d0z0phi = False 


### pt Residuals
pt_Residual_Vs_Pt_PtRange              = False
pt_Residual_Vs_Pt_EtaRange             = False
pt_Residual_Vs_Eta_EtaRange            = False
pt_Residual_Vs_Eta_PtRange_EtaRange    = False
pt_Residual_Vs_Pt_CIF_PtRange          = False 
pt_Residual_Vs_Eta_LMH_EtaRange        = False 
### pt Resolutions
pt_Resolution_Vs_Pt                    = False # datasets (opt)
pt_Resolution_Vs_Eta                   = False # datasets (opt)
pt_Resolution_Vs_Pt_EtaRange           = False # datasets (opt)
pt_Resolution_Vs_Eta_PtRange           = False # datasets (opt)
pt_Resolution_Vs_Pt_CIF                = False 
pt_Resolution_Vs_Eta_LMH               = False 


### phi Residuals
phi_Residual_Vs_Pt_PtRange             = False 
phi_Residual_Vs_Pt_EtaRange            = False 
phi_Residual_Vs_Eta_EtaRange           = False 
phi_Residual_Vs_Eta_PtRange_EtaRange   = False 
phi_Residual_Vs_Pt_CIF_PtRange         = False 
phi_Residual_Vs_Eta_LMH_EtaRange       = False 
### phi Resolutions
phi_Resolution_Vs_Pt                   = False
phi_Resolution_Vs_Eta                  = False
phi_Resolution_Vs_Pt_EtaRange          = False
phi_Resolution_Vs_Eta_PtRange          = False
phi_Resolution_Vs_Pt_CIF               = False
phi_Resolution_Vs_Eta_LMH              = False


### eta Residuals
eta_Residual_Vs_Pt_PtRange             = False # 1 dataset
eta_Residual_Vs_Pt_EtaRange            = False # 1 dataset
eta_Residual_Vs_Eta_EtaRange           = False # 1 dataset
eta_Residual_Vs_Eta_PtRange_EtaRange   = False # 1 dataset
eta_Residual_Vs_Pt_CIF_PtRange         = False 
eta_Residual_Vs_Eta_LMH_EtaRange       = False 
### eta Resolutions
eta_Resolution_Vs_Pt                   = True  # datasets
eta_Resolution_Vs_Eta                  = True  # datasets
eta_Resolution_Vs_Pt_EtaRange          = False # datasets (opt)
eta_Resolution_Vs_Eta_PtRange          = False # datasets (opt)
eta_Resolution_Vs_Pt_CIF               = False 
eta_Resolution_Vs_Eta_LMH              = False 


### z0 Residuals
z0_Residual_Vs_Pt_PtRange             = False # 1 dataset
z0_Residual_Vs_Pt_EtaRange            = False # 1 dataset
z0_Residual_Vs_Eta_EtaRange           = False # 1 dataset
z0_Residual_Vs_Eta_PtRange_EtaRange   = False # 1 dataset
z0_Residual_Vs_Pt_CIF_PtRange         = False 
z0_Residual_Vs_Eta_LMH_EtaRange       = False
### z0 Resolutions
z0_Resolution_Vs_Pt                   = True  # datasets
z0_Resolution_Vs_Eta                  = True  # datasets
z0_Resolution_Vs_Pt_EtaRange          = False # datasets (opt)
z0_Resolution_Vs_Eta_PtRange          = False # datasets (opt)
z0_Resolution_Vs_Pt_CIF               = False 
z0_Resolution_Vs_Eta_LMH              = False


### d0 Residuals
d0_Residual_Vs_Pt_PtRange             = False # 1 dataset
d0_Residual_Vs_Pt_EtaRange            = False # 1 dataset
d0_Residual_Vs_Eta_EtaRange           = False # 1 dataset
d0_Residual_Vs_Eta_PtRange_EtaRange   = False # 1 dataset
d0_Residual_Vs_Pt_CIF_PtRange         = False 
d0_Residual_Vs_Eta_LMH_EtaRange       = False 
### d0 Resolutions
d0_Resolution_Vs_Pt                   = True  # datasets
d0_Resolution_Vs_Eta                  = True  # datasets
d0_Resolution_Vs_Pt_EtaRange          = False # datasets (opt)
d0_Resolution_Vs_Eta_PtRange          = False # datasets (opt)
d0_Resolution_Vs_Pt_CIF               = False 
d0_Resolution_Vs_Eta_LMH              = False 


###############################################################
### Helper Functions
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
### General Settings
###############################################################
bDoPixelTks   = True
bRatio        = True
saveFormats   = ["png", "pdf"]
ext           = ""
tkTypeDir     = "5FitParams/pixTks"
datasetList   = ["SingleMuon_NoPU"] #SingleMuPlus_NoPU , PiPlus, SingleMuon_NoPU
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/AllDatasets/"
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/Default/"
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/Priority4/"
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/AtLeast4Not/"
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/AtLeast4/"
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/3PixHitsCand/"
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/3PixHitsCand_EtaLE1p0/"
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/3PixHitsCand_EtaGE1p0_LE1p6/"
savePath      = ""
#savePath      = "/Users/attikis/talks/post_doc.git/HLTaus/2015/L1Tracks_13November2015/figures/" + tkTypeDir + "/" + datasetList[0] + "/Default/HitPatterns/"


###############################################################
### Save Definitions
###############################################################r
### Algorithms
datasetPaths_Test        = CreateDatasetDict("Macros/PixelRefitting/PixelRefitting_Histograms_", ext)
datasetPaths_Default     = CreateDatasetDict("Macros/PixelRefitting/results/Default/PixelRefitting_Histograms_", ext)
datasetPaths_Priority4   = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/PixelRefitting_Histograms_", ext)
datasetPaths_AtLeast4Not = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/PixelRefitting_Histograms_", ext)
datasetPaths_AtLeast4    = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4/PixelRefitting_Histograms_", ext)

### Algos & Windows
datasetPaths_Default_2xWindow     = CreateDatasetDict("Macros/PixelRefitting/results/rphi6mm_rz10mm/Default/PixelRefitting_Histograms_"    , ext)
datasetPaths_AtLeast4Not_2xWindow = CreateDatasetDict("Macros/PixelRefitting/results/rphi6mm_rz10mm/AtLeast4Not/PixelRefitting_Histograms_", ext)
datasetPaths_Priority4_2xWindow   = CreateDatasetDict("Macros/PixelRefitting/results/rphi6mm_rz10mm/Priority4/PixelRefitting_Histograms_"  , ext)

### 1 Lost Hit Investigation
datasetPaths_3PixHitsCand                = CreateDatasetDict("Macros/PixelRefitting/results/Default/3PixHitsCand/PixelRefitting_Histograms_" , ext)
datasetPaths_3PixHitsCand_EtaLE1p0       = CreateDatasetDict("Macros/PixelRefitting/results/Default/3PixHitsCand_EtaLE1p0/PixelRefitting_Histograms_" , ext)
datasetPaths_3PixHitsCand_EtaGE1p0_LE1p6 = CreateDatasetDict("Macros/PixelRefitting/results/Default/3PixHitsCand_EtaGE1p0_LE1p6/PixelRefitting_Histograms_" , ext)
datasetPaths_3PixHitsCand_EtaGE1p6_LE2p5 = CreateDatasetDict("Macros/PixelRefitting/results/Default/3PixHitsCand_EtaGE1p6_LE2p5/PixelRefitting_Histograms_" , ext)

### Default Algo
datasetPaths_Default_HP_7             = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_7/PixelRefitting_Histograms_" , ext)
datasetPaths_Default_HP_11            = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_11/PixelRefitting_Histograms_" , ext)
datasetPaths_Default_HP_13            = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_13/PixelRefitting_Histograms_" , ext)
datasetPaths_Default_HP_14            = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_14/PixelRefitting_Histograms_", ext)
datasetPaths_Default_HP_15            = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_15/PixelRefitting_Histograms_" , ext)
datasetPaths_Default_HP_23            = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_23/PixelRefitting_Histograms_" , ext)
datasetPaths_Default_HP_51            = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_51/PixelRefitting_Histograms_" , ext)
datasetPaths_Default_HP_112           = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_112/PixelRefitting_Histograms_", ext)
datasetPaths_Default_HP_113           = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_113/PixelRefitting_Histograms_", ext)
datasetPaths_Default_HP_7or11or13     = CreateDatasetDict("Macros/PixelRefitting/results/Default/HitPattern_7or11or13/PixelRefitting_Histograms_" , ext)
datasetPaths_Default_HP_7or11or13or14 = CreateDatasetDict("Macros/PixelRefitting/results/Default//HitPattern_7or11or13or14/PixelRefitting_Histograms_" , ext)

### Priority4 Algo
datasetPaths_Priority4_HP_7             = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_7/PixelRefitting_Histograms_" , ext)
datasetPaths_Priority4_HP_11            = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_11/PixelRefitting_Histograms_" , ext)
datasetPaths_Priority4_HP_13            = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_13/PixelRefitting_Histograms_" , ext)
datasetPaths_Priority4_HP_14            = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_14/PixelRefitting_Histograms_", ext)
datasetPaths_Priority4_HP_15            = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_15/PixelRefitting_Histograms_" , ext)
datasetPaths_Priority4_HP_23            = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_23/PixelRefitting_Histograms_" , ext)
datasetPaths_Priority4_HP_51            = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_51/PixelRefitting_Histograms_" , ext)
datasetPaths_Priority4_HP_112           = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_112/PixelRefitting_Histograms_", ext)
datasetPaths_Priority4_HP_113           = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_113/PixelRefitting_Histograms_", ext)
datasetPaths_Priority4_HP_7or11or13     = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_7or11or13/PixelRefitting_Histograms_" , ext)
datasetPaths_Priority4_HP_7or11or13or14 = CreateDatasetDict("Macros/PixelRefitting/results/Priority4/HitPattern_7or11or13or14/PixelRefitting_Histograms_" , ext)


### AtLeast4Not Algo
datasetPaths_AtLeast4Not_HP_7             = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_7/PixelRefitting_Histograms_" , ext)
datasetPaths_AtLeast4Not_HP_11            = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_11/PixelRefitting_Histograms_" , ext)
datasetPaths_AtLeast4Not_HP_13            = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_13/PixelRefitting_Histograms_" , ext)
datasetPaths_AtLeast4Not_HP_14            = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_14/PixelRefitting_Histograms_", ext)
datasetPaths_AtLeast4Not_HP_15            = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_15/PixelRefitting_Histograms_" , ext)
datasetPaths_AtLeast4Not_HP_23            = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_23/PixelRefitting_Histograms_" , ext)
datasetPaths_AtLeast4Not_HP_51            = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_51/PixelRefitting_Histograms_" , ext)
datasetPaths_AtLeast4Not_HP_112           = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_112/PixelRefitting_Histograms_", ext)
datasetPaths_AtLeast4Not_HP_113           = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_113/PixelRefitting_Histograms_", ext)
datasetPaths_AtLeast4Not_HP_7or11or13     = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_7or11or13/PixelRefitting_Histograms_" , ext)
datasetPaths_AtLeast4Not_HP_7or11or13or14 = CreateDatasetDict("Macros/PixelRefitting/results/AtLeast4Not/HitPattern_7or11or13or14/PixelRefitting_Histograms_" , ext)


###############################################################
### Histo Definitions
###############################################################
yMin          =  0.0
ptMax         = 50.0
etaMax        =  2.4 #2.5
nPt_Range     = int(ptMax/5.0)
nEta_Range    = 12 #25 or 12
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
    "yLabel": "Entries / %0.0f" , "yUnits": ""         , "yMin": 1E00, "yMax": None , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True
    , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

PtL = {
    "xLabel": "p_{T}"           , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": +5.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f" , "yUnits": ""         , "yMin": 1E00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


Eta = {
    "xLabel": "#eta"           , "xUnits": ""     , "xMin": -etaMax , "xMax": +etaMax, "binWidthX": None, "xCutLines": [0], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""     , "yMin": +1e00, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90, "xCutBoxes": [[-1.0, -1.6, ROOT.kBlue], [+1.0, +1.6, ROOT.kBlue]]
    #"xCutBoxes": [[-1.0, -1.6, ROOT.kBlue], [-1.0, +1.0, ROOT.kRed], [+1.0, +1.6, ROOT.kBlue]]
    #"xCutBoxes": [[-1.0, -1.6, ROOT.kBlue], [+1.0, +1.6, ROOT.kBlue]]
    #"xCutBoxes": [[-1.0, +1.0, ROOT.kRed]]
}


Phi = {
    "xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.2 , "xMax": +3.2, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": +1E00, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

z0 = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -25.0, "xMax": +25.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": +1E00, "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


d0 = {
    "xLabel": "d_{0}"          , "xUnits": "#mum", "xMin": -40 , "xMax": +40 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": ""    , "yMin": +1e0, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}

dxy = {
    "xLabel": "d_{xy} = (x_{0}^{2} + y_{0}^{2})^{#frac{1}{2}}", "xUnits": "#mum", "xMin": 0, "xMax": +50 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True,
    "yLabel": "Entries / %0.2f", "yUnits": "", "yMin": +1e0, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.80, "yLegMax": 0.90
}


ChiSq  = {
    "xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +150.0, "binWidthX": 2.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f" , "yUnits": ""    , "yMin": 1E-04, "yMax": None , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.15, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.86, "yLegMin": 0.62, "yLegMax": 0.82
}


RedChiSq = {
    "xLabel": "#chi^{2}_{N}"   , "xUnits": "", "xMin": +0.00, "xMax": +10.0, "binWidthX": 0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f", "yUnits": "", "yMin": 1e-04, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 2.15, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.86, "yLegMin": 0.62, "yLegMax": 0.82
}



PtRes = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual" , "xUnits": "GeVc^{-1}", "xMin": -5.0, "xMax": +5.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""         , "yMin": +1e0, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PtResAlt = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual" , "xUnits": "GeVc^{-1}", "xMin": -5.0, "xMax": +5.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.1f", "yUnits": ""         , "yMin": +1e0, "yMax": None, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True ,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PtResRel = { #residual = (L1 - Sim)
    "xLabel": "p_{T} residual / p_{T}", "xUnits": "", "xMin": -0.12, "xMax": +0.12, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.3f"       , "yUnits": "", "yMin": +1e+0, "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


EtaRes = { #residual = (L1 - Sim)
    "xLabel": "#eta residual"   , "xUnits": "", "xMin": -0.004, "xMax": 0.004, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.4f" , "yUnits": "", "yMin": yMin  , "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


EtaResN = { #residual = (L1 - Sim)
    "xLabel": "#eta residual"   , "xUnits": "", "xMin": -0.004, "xMax": 0.004, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.4f" , "yUnits": "", "yMin": 1e-4  , "yMax": 1e0  , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PhiRes = { #residual = (L1 - Sim)
    "xLabel": "#phi residual"  , "xUnits": "rads", "xMin": -0.01, "xMax": +0.01, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.4f", "yUnits": ""    , "yMin": +1E00 , "yMax": None  , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True ,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PhiResN = { #residual = (L1 - Sim)
    "xLabel": "#phi residual"  , "xUnits": "rads", "xMin": -0.01, "xMax": +0.01, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.4f", "yUnits": ""    , "yMin": 1e-4 , "yMax": 1e0  , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True ,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


Z0Res = { #residual =  (L1 - Sim)
    "xLabel": "z_{0} residual" , "xUnits": "cm", "xMin": -0.02, "xMax": +0.02, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.3f", "yUnits": ""  , "yMin": +yMin, "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


Z0ResN = { #residual =  (L1 - Sim)
    "xLabel": "z_{0} residual" , "xUnits": "cm", "xMin": -0.02, "xMax": +0.02, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.3f", "yUnits": ""  , "yMin": 1e-4, "yMax": 1e0   , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


d0Res = { #residual =  (L1 - Sim)
    "xLabel": "d_{0} residual" , "xUnits": "cm", "xMin": -0.02, "xMax": +0.02, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.3f", "yUnits": ""  , "yMin": +yMin, "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

d0ResN = { #residual =  (L1 - Sim)
    "xLabel": "d_{0} residual" , "xUnits": "cm", "xMin": -0.02, "xMax": +0.02, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.3f", "yUnits": ""  , "yMin": 1e-4, "yMax": 1e0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
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


dPtEta = {
    "xLabel": "#eta / %0.2f"              , "xUnits": ""         , "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "|p_{T} residual| / % 0.2f" , "yUnits": "GeVc^{-1}", "yMin":  0.0, "yMax": +5.0, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 0, "zMax": None , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


PtPtRes = {
    "xLabel": "p_{T}"                    , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax , "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "p_{T} resolution / % 0.0f", "yUnits": ""         , "yMin": 0.0 , "yMax": etaMax, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "logYRatio": False, "logXRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PtEtaRes = {
    "xLabel": "p_{T}"                   , "xUnits": "GeVc^{-1}", "xMin": 0.0, "xMax": ptMax , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    #"yLabel": "#eta resolution / % 0.1f", "yUnits": ""         , "yMin": 0.0, "yMax": None  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "yLabel": "#eta resolution / % 0.1f", "yUnits": ""         , "yMin": 0.0, "yMax": +3e-3, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": True, "yMinRatio": 0.0 , "yMaxRatio": 1.75 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PtPhiRes = {
    "xLabel": "p_{T}"                   , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#phi resolution / % 0.1f", "yUnits": ""         , "yMin": 0.00, "yMax": 0.005, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PtZ0Res = {
    "xLabel": "p_{T}"                         , "xUnits": "GeVc^{-1}", "xMin": 0.0, "xMax": ptMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "z_{0} resolution (cm) / % 0.1f", "yUnits": ""         , "yMin": 0.0, "yMax": None, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False,
    #"yLabel": "z_{0} resolution (cm) / % 0.1f", "yUnits": ""         , "yMin": 0.0, "yMax": +0.02, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": True, "yMinRatio": 0.0 , "yMaxRatio": 1.75 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

Ptd0Res = {
    "xLabel": "p_{T}"                         , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax, "binWidthX": None, "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "d_{0} resolution (cm) / % 0.1f", "yUnits": ""         , "yMin": 0.0 , "yMax": None , "binWidthY": None, "yCutBoxes": [], "gridY": True, "logY": False,
    #"yLabel": "d_{0} resolution (cm) / % 0.1f", "yUnits": ""         , "yMin": 0.0 , "yMax": 9e-3 , "binWidthY": None, "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": True, "yMinRatio": 0.0 , "yMaxRatio": 1.75 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PtPtResRel = {
    "xLabel": "p_{T}"                            , "xUnits": "GeVc^{-1}", "xMin": 0.00, "xMax": ptMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False,
    "yLabel": "p_{T} resolution / p_{T} / % 0.1f", "yUnits": ""         , "yMin": 0.00, "yMax": None , "binWidthY": None, "yCutLines": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

EtaPtRes = {
    "xLabel": "#eta"                     , "xUnits": "", "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "p_{T} resolution / % 0.1f", "yUnits": "", "yMin": 0.0    , "yMax": None   , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    #"yLabel": "p_{T} resolution / % 0.1f", "yUnits": "", "yMin": 0.0    , "yMax":  +0.4  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


EtaEtaRes = {
    "xLabel": "#eta"                    , "xUnits": "", "xMin": 0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#eta resolution / % 0.1f", "yUnits": "", "yMin": 0.0, "yMax": None   , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    #"yLabel": "#eta resolution / % 0.1f", "yUnits": "", "yMin": 0.0, "yMax": 3e-3, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": True, "yMinRatio": 0.0 , "yMaxRatio": 1.75, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

EtaPhiRes = {
    "xLabel": "#eta"                    , "xUnits": "", "xMin": -etaMax, "xMax": +etaMax , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "#phi resolution / % 0.2f", "yUnits": "", "yMin": 0.0 , "yMax":  None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    #"yLabel": "#phi resolution / % 0.2f", "yUnits": "", "yMin": 0.0 , "yMax":  0.005, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

EtaZ0Res = {
    "xLabel": "#eta"                          , "xUnits": "", "xMin": 0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "z_{0} resolution (cm) / % 0.3f", "yUnits": "", "yMin": 0.0, "yMax": None   , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    #"yLabel": "z_{0} resolution (cm) / % 0.3f", "yUnits": "", "yMin": 0.0, "yMax": +0.02, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": True, "yMinRatio": 0.0 , "yMaxRatio": 1.75 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

Etad0Res = {
    "xLabel": "#eta"                          , "xUnits": "", "xMin": +0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "d_{0} resolution (cm) / % 0.2f", "yUnits": "", "yMin": +0.0, "yMax": None   , "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False,
    #"yLabel": "d_{0} resolution (cm) / % 0.2f", "yUnits": "", "yMin": +0.0, "yMax": +0.012  , "binWidthY": None, "yCutLines": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": True, "yMinRatio": 0.0 , "yMaxRatio": 1.75 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "logXRatio": False, "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

EtaPtResRel = {
    "xLabel": "#eta"                             , "xUnits": "" , "xMin": +0.0, "xMax": +etaMax, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "p_{T} resolution / p_{T} / % 0.1f", "yUnits": "" , "yMin": +0.0, "yMax": +0.1, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


EffPt = {
    "xLabel": "p_{T}"              , "xUnits": "GeVc^{-1}", "xMin": 0.0 , "xMax": +50.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Efficiency / %0.0f" , "yUnits": ""         , "yMin": 0.00, "yMax":  +1.15, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False ,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 1.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.66, "xLegMax": 0.90, "yLegMin": 0.25, "yLegMax": 0.40
}

EffPtL = {
    "xLabel": "p_{T}"              , "xUnits": "GeVc^{-1}", "xMin": 0.0, "xMax": +5.0  , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Efficiency / %0.1f" , "yUnits": ""     , "yMin": 0.00, "yMax": +1.15, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 1.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "logYRatio": False, "xLegMin": 0.66, "xLegMax": 0.90, "yLegMin": 0.25, "yLegMax": 0.40
}


EffEta = {
    #"xLabel": "#eta", "xUnits": "" , "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutLines": [-1.1, -1.5, +1.1, +1.5], "xCutBoxes": [], "gridX": True, "logX": False,
    "xLabel": "#eta", "xUnits": "" , "xMin": -etaMax, "xMax": +etaMax, "binWidthX": None, "xCutBoxes": [[-1.1, -1.5, ROOT.kBlack], [+1.1, +1.5, ROOT.kBlack]],
    "yLabel": "Efficiency / %0.2f" , "yUnits": "" , "yMin":  0.0, "yMax": +1.15, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio":  0.0, "yMaxRatio": 1.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "logYRatio": False, "xLegMin": 0.66, "xLegMax": 0.90, "yLegMin": 0.25, "yLegMax": 0.40
}


EffPhi = {
    "xLabel": "#phi"              , "xUnits": "rads", "xMin": -3.2, "xMax": +3.2, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.2f", "yUnits": ""    , "yMin": 0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0 , "yMaxRatio": 2.15 , "normaliseTo": None, "drawOptions": "P", "legOptions": "LP", 
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


Effz0 = {
    "xLabel": "z_{0}"              , "xUnits": "cm", "xMin": -25.0 , "xMax": +25.0, "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.2f" , "yUnits": ""  , "yMin":  +0.00, "yMax":  +1.2, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False ,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.4, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


Effd0 = {
    "xLabel": "d_{0}"             , "xUnits": "#mum", "xMin": -10  , "xMax": +10 , "binWidthX": None, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Efficiency / %0.2f", "yUnits": ""    , "yMin": +0.00, "yMax": +1.2, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.30, "yLegMax": 0.42
}


NPixHits = {
    "xLabel": "Number of Hits"    , "xUnits": "", "xMin": -0.5  , "xMax": +5.5, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": 40e3, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    #"ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PixHitsRho = {
    "xLabel": "Hit radius"     , "xUnits": "cm", "xMin": 0.0, "xMax": +20, "binWidthX": 0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1.0, "yMax": 1e6, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PixHitsZ = {
    "xLabel": "Hit z-coordinate", "xUnits": "cm", "xMin": -60, "xMax": +60 , "binWidthX": 1.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PixHitsType = {
    "xLabel": "Hit Type"       , "xUnits": "", "xMin": -3.0 , "xMax": +4.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}


PixHitsPattern = {
    "xLabel": "Hit Pattern"    , "xUnits": "", "xMin": -0.5  , "xMax": 120.5, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    #"ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80, "xCutBoxes": []
}

PixHitsPatternHL = {
    "xLabel": "Hit Pattern"    , "xUnits": "", "xMin": -0.5  , "xMax": 120.5, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    #"yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": 0.6, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    #"ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80,
    "xCutBoxes": [
        [  7.0,   7.0, ROOT.kBlue],
#        [ 15.0,  15.0, ROOT.kRed],
#        [ 23.0,  23.0, ROOT.kRed],
#        [ 51.0,  51.0, ROOT.kRed],
#        [113.0, 113.0, ROOT.kRed],
#        [  7.0,  14.0, ROOT.kBlue],
#        [112.0, 112.0, ROOT.kBlue],
    ]
}


PixHitsPatternVsEta = {
    "xLabel": "Hit Pattern / %0.0f" , "xUnits": "", "xMin": -0.5  , "xMax": None, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False,
    "yLabel": "#eta / %0.2f", "yUnits": "", "yMin": -etaMax , "yMax": +etaMax , "binWidthY": None, "yCutLines": [0.8, 1.6], "yCutBoxes": [], "gridY": True, "logY": False, 
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": None, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}

NPixHitsCand = { #xenios
    #"xLabel": "Number of Candidate Hits"    , "xUnits": "", "xMin": -0.5  , "xMax": +15.5, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    #"xLabel": "Number of Candidate Hits"    , "xUnits": "", "xMin": -0.5  , "xMax": +29.5, "binWidthX": None, "xCutLines": [4], "xCutBoxes": [], "gridX": True, "logX": False, 
    #"yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    # "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.2, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLabel": "Number of Candidate Hits"    , "xUnits": "", "xMin": -0.5  , "xMax": +14.5, "binWidthX": None, "xCutLines": [4], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": 1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.2, "yMaxRatio": 2.15, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.62, "xLegMax": 0.90, "yLegMin": 0.70, "yLegMax": 0.80
}


PixHitsRhoCand = {
    "xLabel": "Candidate Hit radius"     , "xUnits": "cm", "xMin": 0.0, "xMax": +20, "binWidthX": 0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1.0, "yMax": 1e6, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PixHitsZCand = {
    "xLabel": "Candidate Hit z-coordinate", "xUnits": "cm", "xMin": -60, "xMax": +60 , "binWidthX": 1.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PixHitsTypeCand = { #xenios
    "xLabel": "Candidate Hit Type"          , "xUnits": "", "xMin": -3.0  , "xMax": +4.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    #"ratioLabel": "Ratio", "ratio": bRatio, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

EtaVsN = {
    "xLabel": "#eta / %0.2f", "xUnits": "", "xMin": -etaMax , "xMax": +etaMax , "binWidthX": None, "xCutLines": [-1.4, +1.4], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Number of Hits / %0.1f"    , "yUnits": "", "yMin": -0.5  , "yMax": +5.5, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, 
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": None, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}

EtaVsNCand = {
    "xLabel": "#eta / %0.2f", "xUnits": "", "xMin": -etaMax , "xMax": +etaMax , "binWidthX": None, "xCutLines": [-1.4, +1.4], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Number of Candidate Hits / %0.1f", "yUnits": "", "yMin": -0.5  , "yMax": +29.5, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, 
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": None, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}


PixHitsZVsRho = {
    "xLabel": "Hit z-coordinate / %0.2f", "xUnits": "cm", "xMin": -60, "xMax": +60 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit radius / %0.2f"      , "yUnits": "cm", "yMin": 0  , "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    #"zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 3e1, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True,
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 1e2, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}

PixHitsZVsRhoNorm = {
    "xLabel": "Hit z-coordinate / %0.2f", "xUnits": "cm", "xMin": -60, "xMax": +60 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit radius / %0.2f"      , "yUnits": "cm", "yMin": 0  , "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 5e-5, "zMax": 1e-03, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": "One", 
}

PixHitsXVsY = {
    "xLabel": "Hit x-coordinate / %0.2f", "xUnits": "cm", "xMin": -20, "xMax": +20 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit y-coordinate / %0.2f", "yUnits": "cm", "yMin": -20, "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 1e2, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}

PixHitsXVsYNorm = {
    "xLabel": "Hit x-coordinate / %0.2f", "xUnits": "cm", "xMin": -20, "xMax": +20 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit y-coordinate / %0.2f", "yUnits": "cm", "yMin": -20, "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 1e-6, "zMax": 5e-03, "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": "One", 
}

PixHitsZVsRhoCand = {
    "xLabel": "Hit z-coordinate / %0.2f", "xUnits": "cm", "xMin": -60, "xMax": +60 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit radius / %0.2f"      , "yUnits": "cm", "yMin": 0  , "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    #"zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 3e1 , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True,
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 1e3 , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}

PixHitsXVsYCand = {
    "xLabel": "Hit x-coordinate / %0.2f", "xUnits": "cm", "xMin": -20, "xMax": +20 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit y-coordinate / %0.2f", "yUnits": "cm", "yMin": -20, "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 1e3 , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}


PixHitsRhoShared = {
    "xLabel": "Hit radius"     , "xUnits": "cm", "xMin": 0.0, "xMax": +20, "binWidthX": 0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1.0, "yMax": 1e6, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PixHitsZShared = {
    "xLabel": "Hit z-coordinate", "xUnits": "cm", "xMin": -60, "xMax": +60 , "binWidthX": 1.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PixHitsTypeShared = {
    "xLabel": "Hit Type"          , "xUnits": "", "xMin": -3.0  , "xMax": +4.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +0.00, "yMax": None, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, 
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 2.15, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.18, "xLegMax": 0.4, "yLegMin": 0.78, "yLegMax": 0.90
}

PixHitsZVsRhoShared = {
    "xLabel": "Hit z-coordinate / %0.2f", "xUnits": "cm", "xMin": -60, "xMax": +60 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit radius / %0.2f"      , "yUnits": "cm", "yMin": 0  , "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 100 , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
}

PixHitsXVsYShared = {
    "xLabel": "Hit x-coordinate / %0.2f", "xUnits": "cm", "xMin": -20, "xMax": +20 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Hit y-coordinate / %0.2f", "yUnits": "cm", "yMin": -20, "yMax": +20 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False,
    "zLabel": "Entries", "zUnits": "", "zMin": 1.0, "zMax": 100 , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": True, 
    "drawOptions": "COLZ", "logYRatio": False, "logXRatio": False, "xLegMin": 0.68, "xLegMax": 0.90, "yLegMin": 0.78, "yLegMax": 0.80, "normaliseTo": None, 
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
tp_dxy    = m_histos.TH1orTH2( "", "tp_dxy"   , "TP"      , None, **dxy)

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
efficiency_pt    = m_histos.TH1orTH2( "", "efficiency_pt"   , "overall", None, **EffPt )
efficiency_pt_L  = m_histos.TH1orTH2( "", "efficiency_pt_L" , "" + L   , None, **EffPtL)
efficiency_pt_C  = m_histos.TH1orTH2( "", "efficiency_pt_C" , "" + C   , None, **EffPt)
efficiency_pt_I  = m_histos.TH1orTH2( "", "efficiency_pt_I" , "" + I   , None, **EffPt)
efficiency_pt_F  = m_histos.TH1orTH2( "", "efficiency_pt_F" , "" + F   , None, **EffPt)
efficiency_eta   = m_histos.TH1orTH2( "", "efficiency_eta"  , "overall", None, **EffEta)
efficiency_eta_L = m_histos.TH1orTH2( "", "efficiency_eta_L", "" + L   , None, **EffEta)
efficiency_eta_M = m_histos.TH1orTH2( "", "efficiency_eta_M", "" + M   , None, **EffEta)
efficiency_eta_H = m_histos.TH1orTH2( "", "efficiency_eta_H", "" + H   , None, **EffEta)
efficiency_phi   = m_histos.TH1orTH2( "", "efficiency_phi"  , "overall", None, **EffPhi)
efficiency_z0    = m_histos.TH1orTH2( "", "efficiency_z0"   , "overall", None, **Effz0 )
efficiency_d0    = m_histos.TH1orTH2( "", "efficiency_d0"   , "overall", None, **Effd0 )

### TTTrack (TP-Matched): ChiSq histograms (last bin is an overflow bin)
label = ""
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
residual_pt      = m_histos.TH1orTH2( "", "residual_pt"     , "inclusive", None, **PtRes   )
residual_pt_C    = m_histos.TH1orTH2( "", "residual_pt_C"   , label + C  , None, **PtResAlt)
residual_pt_I    = m_histos.TH1orTH2( "", "residual_pt_I"   , label + I  , None, **PtResAlt)
residual_pt_F    = m_histos.TH1orTH2( "", "residual_pt_F"   , label + F  , None, **PtResAlt)

residual_eta     = m_histos.TH1orTH2( "", "residual_eta"    , "inclusive", None, **EtaResN  )
residual_eta_C   = m_histos.TH1orTH2( "", "residual_eta_C"  , label + C  , None, **EtaResN  )
residual_eta_I   = m_histos.TH1orTH2( "", "residual_eta_I"  , label + I  , None, **EtaResN  )
residual_eta_F   = m_histos.TH1orTH2( "", "residual_eta_F"  , label + F  , None, **EtaResN  )

residual_phi     = m_histos.TH1orTH2( "", "residual_phi"    , "inclusive", None, **PhiResN  )
residual_phi_C   = m_histos.TH1orTH2( "", "residual_phi_C"  , label + C  , None, **PhiResN  )
residual_phi_I   = m_histos.TH1orTH2( "", "residual_phi_I"  , label + I  , None, **PhiResN  )
residual_phi_F   = m_histos.TH1orTH2( "", "residual_phi_F"  , label + F  , None, **PhiResN  )

residual_z0      = m_histos.TH1orTH2( "", "residual_z0"     , "inclusive", None, **Z0ResN   )
residual_z0_C    = m_histos.TH1orTH2( "", "residual_z0_C"   , label + C  , None, **Z0ResN   )
residual_z0_I    = m_histos.TH1orTH2( "", "residual_z0_I"   , label + I  , None, **Z0ResN   )
residual_z0_F    = m_histos.TH1orTH2( "", "residual_z0_F"   , label + F  , None, **Z0ResN   )

residual_d0      = m_histos.TH1orTH2( "", "residual_d0"     , "inclusive", None, **d0ResN  )
residual_d0_C    = m_histos.TH1orTH2( "", "residual_d0_C"   , label + C  , None, **d0ResN  )
residual_d0_I    = m_histos.TH1orTH2( "", "residual_d0_I"   , label + I  , None, **d0ResN  )
residual_d0_F    = m_histos.TH1orTH2( "", "residual_d0_F"   , label + F  , None, **d0ResN  )

### Pixel Hits
label = "hits"
pixTk_pixHits_N             = m_histos.TH1orTH2( "", "pixTk_pixHits_N"            , label, None, **NPixHits    )
pixTk_pixHits_Rho           = m_histos.TH1orTH2( "", "pixTk_pixHits_Rho"          , label, None, **PixHitsRho  )
pixTk_pixHits_Z             = m_histos.TH1orTH2( "", "pixTk_pixHits_Z"            , label, None, **PixHitsZ    )
pixTk_pixHits_Type          = m_histos.TH1orTH2( "", "pixTk_pixHits_Type"         , label, None, **PixHitsType )
pixTk_pixHits_Pattern       = m_histos.TH1orTH2( "", "pixTk_pixHits_Pattern"      , label, None, **PixHitsPattern )
pixTk_pixHits_PatternHL     = m_histos.TH1orTH2( "", "pixTk_pixHits_Pattern"      , label, None, **PixHitsPatternHL)
pixTk_pixHits_ZVsRho        = m_histos.TH1orTH2( "", "pixTk_pixHits_ZVsRho"       , label, None, **PixHitsZVsRho )
pixTk_pixHits_ZVsRho_Norm   = m_histos.TH1orTH2( "", "pixTk_pixHits_ZVsRho"       , label, None, **PixHitsZVsRho )
pixTk_pixHits_ZVsRho_C      = m_histos.TH1orTH2( "", "pixTk_pixHits_ZVsRho_C"     , label, None, **PixHitsZVsRho )
pixTk_pixHits_ZVsRho_I      = m_histos.TH1orTH2( "", "pixTk_pixHits_ZVsRho_I"     , label, None, **PixHitsZVsRho )
pixTk_pixHits_ZVsRho_F      = m_histos.TH1orTH2( "", "pixTk_pixHits_ZVsRho_F"     , label, None, **PixHitsZVsRho )
pixTk_pixHits_ZVsRho_EtaGE2 = m_histos.TH1orTH2( "", "pixTk_pixHits_ZVsRho_EtaGE2", label, None, **PixHitsZVsRho )
pixTk_pixHits_XVsY          = m_histos.TH1orTH2( "", "pixTk_pixHits_XVsY"         , label, None, **PixHitsXVsY )
pixTk_pixHits_XVsY_Norm     = m_histos.TH1orTH2( "", "pixTk_pixHits_XVsY"         , label, None, **PixHitsXVsY )
pixTk_pixHits_XVsY_C        = m_histos.TH1orTH2( "", "pixTk_pixHits_XVsY_C"       , label, None, **PixHitsXVsY )
pixTk_pixHits_XVsY_I        = m_histos.TH1orTH2( "", "pixTk_pixHits_XVsY_I"       , label, None, **PixHitsXVsY )
pixTk_pixHits_XVsY_F        = m_histos.TH1orTH2( "", "pixTk_pixHits_XVsY_F"       , label, None, **PixHitsXVsY )
pixTk_pixHits_XVsY_EtaGE2   = m_histos.TH1orTH2( "", "pixTk_pixHits_XVsY_EtaGE2"  , label, None, **PixHitsXVsY )
pixTk_pixHits_PatternVsEta  = m_histos.TH1orTH2( "", "pixTk_pixHits_PatternVsEta" , label, None, **PixHitsPatternVsEta )

### Shared Pixel Hits
label = "shared hits"
pixTk_sharedPixHits_N             = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_N"            , label, None, **NPixHits    )
pixTk_sharedPixHits_Rho           = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_Rho"          , label, None, **PixHitsRhoShared    )
pixTk_sharedPixHits_Z             = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_Z"            , label, None, **PixHitsZShared      )
pixTk_sharedPixHits_Type          = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_Type"         , label, None, **PixHitsTypeShared   )
pixTk_sharedPixHits_ZVsRho        = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_ZVsRho"       , label, None, **PixHitsZVsRhoShared )
pixTk_sharedPixHits_ZVsRho_C      = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_ZVsRho_C"     , label, None, **PixHitsZVsRhoShared )
pixTk_sharedPixHits_ZVsRho_I      = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_ZVsRho_I"     , label, None, **PixHitsZVsRhoShared )
pixTk_sharedPixHits_ZVsRho_F      = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_ZVsRho_F"     , label, None, **PixHitsZVsRhoShared )
pixTk_sharedPixHits_ZVsRho_EtaGE2 = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_ZVsRho_EtaGE2", label, None, **PixHitsZVsRhoShared )
pixTk_sharedPixHits_XVsY          = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_XVsY"         , label, None, **PixHitsXVsYShared )
pixTk_sharedPixHits_XVsY_C        = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_XVsY_C"       , label, None, **PixHitsXVsYShared )
pixTk_sharedPixHits_XVsY_I        = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_XVsY_I"       , label, None, **PixHitsXVsYShared )
pixTk_sharedPixHits_XVsY_F        = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_XVsY_F"       , label, None, **PixHitsXVsYShared )
pixTk_sharedPixHits_XVsY_EtaGE2   = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_XVsY_EtaGE2"  , label, None, **PixHitsXVsYShared )

### Candidate Pixel Hits
label = "candidate hits"
pixTk_candPixHits_N             = m_histos.TH1orTH2( "", "pixTk_candPixHits_N"            , label, None, **NPixHitsCand      )
pixTk_candPixHits_Rho           = m_histos.TH1orTH2( "", "pixTk_candPixHits_Rho"          , label, None, **PixHitsRhoCand    )
pixTk_candPixHits_Z             = m_histos.TH1orTH2( "", "pixTk_candPixHits_Z"            , label, None, **PixHitsZCand      )
pixTk_candPixHits_Type          = m_histos.TH1orTH2( "", "pixTk_candPixHits_Type"         , label, None, **PixHitsTypeCand   )
pixTk_candPixHits_ZVsRho        = m_histos.TH1orTH2( "", "pixTk_candPixHits_ZVsRho"       , label, None, **PixHitsZVsRhoCand )
pixTk_candPixHits_ZVsRho_C      = m_histos.TH1orTH2( "", "pixTk_candPixHits_ZVsRho_C"     , label, None, **PixHitsZVsRhoCand )
pixTk_candPixHits_ZVsRho_I      = m_histos.TH1orTH2( "", "pixTk_candPixHits_ZVsRho_I"     , label, None, **PixHitsZVsRhoCand )
pixTk_candPixHits_ZVsRho_F      = m_histos.TH1orTH2( "", "pixTk_candPixHits_ZVsRho_F"     , label, None, **PixHitsZVsRhoCand )
pixTk_candPixHits_ZVsRho_EtaGE2 = m_histos.TH1orTH2( "", "pixTk_candPixHits_ZVsRho_EtaGE2", label, None, **PixHitsZVsRhoCand )
pixTk_candPixHits_XVsY          = m_histos.TH1orTH2( "", "pixTk_candPixHits_XVsY"         , label, None, **PixHitsXVsYCand )
pixTk_candPixHits_XVsY_C        = m_histos.TH1orTH2( "", "pixTk_candPixHits_XVsY_C"       , label, None, **PixHitsXVsYCand )
pixTk_candPixHits_XVsY_I        = m_histos.TH1orTH2( "", "pixTk_candPixHits_XVsY_I"       , label, None, **PixHitsXVsYCand )
pixTk_candPixHits_XVsY_F        = m_histos.TH1orTH2( "", "pixTk_candPixHits_XVsY_F"       , label, None, **PixHitsXVsYCand )
pixTk_candPixHits_XVsY_EtaGE2   = m_histos.TH1orTH2( "", "pixTk_candPixHits_XVsY_EtaGE2"  , label, None, **PixHitsXVsYCand )

### Resolution vs. pt histograms  
ptRange = ["0-5","5-10", "10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50"] #"50-55", "55-60","60-65","65-70","70-75","75-80","80-85","85-90","90-95","95-100"]
h_residualVsPt_pt      = []
h_residualVsPt_pt_C    = []
h_residualVsPt_pt_I    = []
h_residualVsPt_pt_F    = []

h_residualVsPt_eta     = []
h_residualVsPt_eta_C   = []
h_residualVsPt_eta_I   = []
h_residualVsPt_eta_F   = []

h_residualVsPt_phi     = []
h_residualVsPt_phi_C   = []
h_residualVsPt_phi_I   = []
h_residualVsPt_phi_F   = []

h_residualVsPt_z0      = []
h_residualVsPt_z0_C    = []
h_residualVsPt_z0_I    = []
h_residualVsPt_z0_F    = []

h_residualVsPt_d0      = []
h_residualVsPt_d0_C    = []
h_residualVsPt_d0_I    = []
h_residualVsPt_d0_F    = []

for i in range(0, len(ptRange)):
    residualVsPt_pt      = m_histos.TH1orTH2( "", "residualVsPt_pt_"      + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **PtRes)
    residualVsPt_pt_C    = m_histos.TH1orTH2( "", "residualVsPt_pt_C_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **PtRes)
    residualVsPt_pt_I    = m_histos.TH1orTH2( "", "residualVsPt_pt_I_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **PtRes)
    residualVsPt_pt_F    = m_histos.TH1orTH2( "", "residualVsPt_pt_F_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **PtRes)

    residualVsPt_eta     = m_histos.TH1orTH2( "", "residualVsPt_eta_"     + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **EtaRes)
    residualVsPt_eta_C   = m_histos.TH1orTH2( "", "residualVsPt_eta_C_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **EtaRes)
    residualVsPt_eta_I   = m_histos.TH1orTH2( "", "residualVsPt_eta_I_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **EtaRes)
    residualVsPt_eta_F   = m_histos.TH1orTH2( "", "residualVsPt_eta_F_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **EtaRes)

    residualVsPt_phi     = m_histos.TH1orTH2( "", "residualVsPt_phi_"     + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **PhiRes)
    residualVsPt_phi_C   = m_histos.TH1orTH2( "", "residualVsPt_phi_C_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **PhiRes)
    residualVsPt_phi_I   = m_histos.TH1orTH2( "", "residualVsPt_phi_I_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **PhiRes)
    residualVsPt_phi_F   = m_histos.TH1orTH2( "", "residualVsPt_phi_F_"   + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **PhiRes)

    residualVsPt_z0      = m_histos.TH1orTH2( "", "residualVsPt_z0_"      + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **Z0Res)
    residualVsPt_z0_C    = m_histos.TH1orTH2( "", "residualVsPt_z0_C_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **Z0Res)
    residualVsPt_z0_I    = m_histos.TH1orTH2( "", "residualVsPt_z0_I_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **Z0Res)
    residualVsPt_z0_F    = m_histos.TH1orTH2( "", "residualVsPt_z0_F_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **Z0Res)

    residualVsPt_d0      = m_histos.TH1orTH2( "", "residualVsPt_d0_"      + ptRange[i], ptRange[i] + " GeVc^{-1}"      , None, **d0Res)
    residualVsPt_d0_C    = m_histos.TH1orTH2( "", "residualVsPt_d0_C_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + C, None, **d0Res)
    residualVsPt_d0_I    = m_histos.TH1orTH2( "", "residualVsPt_d0_I_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + I, None, **d0Res)
    residualVsPt_d0_F    = m_histos.TH1orTH2( "", "residualVsPt_d0_F_"    + ptRange[i], ptRange[i] + " GeVc^{-1}, " + F, None, **d0Res)


    ### Append to lists
    h_residualVsPt_pt.append(residualVsPt_pt)
    h_residualVsPt_pt_C.append(residualVsPt_pt_C)
    h_residualVsPt_pt_I.append(residualVsPt_pt_I)
    h_residualVsPt_pt_F.append(residualVsPt_pt_F)

    h_residualVsPt_eta.append(residualVsPt_eta)
    h_residualVsPt_eta_C.append(residualVsPt_eta_C)
    h_residualVsPt_eta_I.append(residualVsPt_eta_I)
    h_residualVsPt_eta_F.append(residualVsPt_eta_F)

    h_residualVsPt_phi.append(residualVsPt_phi)
    h_residualVsPt_phi_C.append(residualVsPt_phi_C)
    h_residualVsPt_phi_I.append(residualVsPt_phi_I)
    h_residualVsPt_phi_F.append(residualVsPt_phi_F)

    h_residualVsPt_z0.append(residualVsPt_z0)
    h_residualVsPt_z0_C.append(residualVsPt_z0_C)
    h_residualVsPt_z0_I.append(residualVsPt_z0_I)
    h_residualVsPt_z0_F.append(residualVsPt_z0_F)

    h_residualVsPt_d0.append(residualVsPt_d0)
    h_residualVsPt_d0_C.append(residualVsPt_d0_C)
    h_residualVsPt_d0_I.append(residualVsPt_d0_I)
    h_residualVsPt_d0_F.append(residualVsPt_d0_F)


### resolution vs. eta histograms
etaRange     = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2"]
etaRangePlus = ["0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.4"]
h_residualVsEta_pt      = []
h_residualVsEta_pt_L    = []
h_residualVsEta_pt_M    = []
h_residualVsEta_pt_H    = []

h_residualVsEta_eta     = []
h_residualVsEta_eta_L   = []
h_residualVsEta_eta_M   = []
h_residualVsEta_eta_H   = []

h_residualVsEta_phi     = []
h_residualVsEta_phi_L   = []
h_residualVsEta_phi_M   = []
h_residualVsEta_phi_H   = []

h_residualVsEta_z0      = []
h_residualVsEta_z0_L    = []
h_residualVsEta_z0_M    = []
h_residualVsEta_z0_H    = []

h_residualVsEta_d0      = []
h_residualVsEta_d0_L    = []
h_residualVsEta_d0_M    = []
h_residualVsEta_d0_H    = []

for i in range(0, len(etaRange)):    
    residualVsEta_pt      = m_histos.TH1orTH2( "", "residualVsEta_pt_"   + etaRangePlus[i], etaRange[i] + " #leq |#eta| < " + etaRangePlus[i]              , None, **PtRes)
    residualVsEta_pt_L    = m_histos.TH1orTH2( "", "residualVsEta_pt_L_" + etaRangePlus[i], L  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **PtRes)
    residualVsEta_pt_M    = m_histos.TH1orTH2( "", "residualVsEta_pt_M_" + etaRangePlus[i], M  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **PtRes)
    residualVsEta_pt_H    = m_histos.TH1orTH2( "", "residualVsEta_pt_H_" + etaRangePlus[i], H  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **PtRes)

    residualVsEta_eta     = m_histos.TH1orTH2( "", "residualVsEta_eta_"     + etaRangePlus[i], etaRange[i] + " #leq |#eta| < " + etaRangePlus[i]              , None, **EtaRes)
    residualVsEta_eta_L   = m_histos.TH1orTH2( "", "residualVsEta_eta_L_"   + etaRangePlus[i], L  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **EtaRes)
    residualVsEta_eta_M   = m_histos.TH1orTH2( "", "residualVsEta_eta_M_"   + etaRangePlus[i], M  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **EtaRes)
    residualVsEta_eta_H   = m_histos.TH1orTH2( "", "residualVsEta_eta_H_"   + etaRangePlus[i], H  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **EtaRes)

    residualVsEta_phi     = m_histos.TH1orTH2( "", "residualVsEta_phi_"     + etaRangePlus[i], etaRange[i] + " #leq |#eta| < " + etaRangePlus[i]              , None, **PhiRes)
    residualVsEta_phi_L   = m_histos.TH1orTH2( "", "residualVsEta_phi_L_"   + etaRangePlus[i], L  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **PhiRes)
    residualVsEta_phi_M   = m_histos.TH1orTH2( "", "residualVsEta_phi_M_"   + etaRangePlus[i], M  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **PhiRes)
    residualVsEta_phi_H   = m_histos.TH1orTH2( "", "residualVsEta_phi_H_"   + etaRangePlus[i], H  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **PhiRes)

    residualVsEta_z0      = m_histos.TH1orTH2( "", "residualVsEta_z0_"   + etaRangePlus[i], etaRange[i] + " #leq |#eta| < " + etaRangePlus[i]              , None, **Z0Res)
    residualVsEta_z0_L    = m_histos.TH1orTH2( "", "residualVsEta_z0_L_" + etaRangePlus[i], L  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **Z0Res)
    residualVsEta_z0_M    = m_histos.TH1orTH2( "", "residualVsEta_z0_M_" + etaRangePlus[i], M  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **Z0Res)
    residualVsEta_z0_H    = m_histos.TH1orTH2( "", "residualVsEta_z0_H_" + etaRangePlus[i], H  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **Z0Res)

    residualVsEta_d0      = m_histos.TH1orTH2( "", "residualVsEta_d0_"   + etaRangePlus[i], etaRange[i] + " #leq |#eta| < " + etaRangePlus[i]              , None, **d0Res)
    residualVsEta_d0_L    = m_histos.TH1orTH2( "", "residualVsEta_d0_L_" + etaRangePlus[i], L  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **d0Res)
    residualVsEta_d0_M    = m_histos.TH1orTH2( "", "residualVsEta_d0_M_" + etaRangePlus[i], M  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **d0Res)
    residualVsEta_d0_H    = m_histos.TH1orTH2( "", "residualVsEta_d0_H_" + etaRangePlus[i], H  + " , " +  etaRange[i] + " #leq |#eta| < " + etaRangePlus[i], None, **d0Res)

    ### Append to lists
    h_residualVsEta_pt.append(residualVsEta_pt)
    h_residualVsEta_pt_L.append(residualVsEta_pt_L)
    h_residualVsEta_pt_M.append(residualVsEta_pt_M)
    h_residualVsEta_pt_H.append(residualVsEta_pt_H)

    h_residualVsEta_eta.append(residualVsEta_eta)
    h_residualVsEta_eta_L.append(residualVsEta_eta_L)
    h_residualVsEta_eta_M.append(residualVsEta_eta_M)
    h_residualVsEta_eta_H.append(residualVsEta_eta_H)

    h_residualVsEta_phi.append(residualVsEta_phi)
    h_residualVsEta_phi_L.append(residualVsEta_phi_L)
    h_residualVsEta_phi_M.append(residualVsEta_phi_M)
    h_residualVsEta_phi_H.append(residualVsEta_phi_H)

    h_residualVsEta_z0.append(residualVsEta_z0)
    h_residualVsEta_z0_L.append(residualVsEta_z0_L)
    h_residualVsEta_z0_M.append(residualVsEta_z0_M)
    h_residualVsEta_z0_H.append(residualVsEta_z0_H)

    h_residualVsEta_d0.append(residualVsEta_d0)
    h_residualVsEta_d0_L.append(residualVsEta_d0_L)
    h_residualVsEta_d0_M.append(residualVsEta_d0_M)
    h_residualVsEta_d0_H.append(residualVsEta_d0_H)


### 2D histograms
logChiSq_eta     = m_histos.TH1orTH2( "", "2d_logchi2_eta"    , "", None, **ChiSqEta)
logChiSq_dof_eta = m_histos.TH1orTH2( "", "2d_logchi2_dof_eta", "", None, **RedChiSqEta)
dz0_eta          = m_histos.TH1orTH2( "", "2d_dz0_eta"        , "", None, **dZ0Eta)
dd0_eta          = m_histos.TH1orTH2( "", "2d_dd0_eta"        , "", None, **dd0Eta)
deta_eta         = m_histos.TH1orTH2( "", "2d_deta_eta"       , "", None, **dEtaEta)
dphi_eta         = m_histos.TH1orTH2( "", "2d_dphi_eta"       , "", None, **dPhiEta)
dpt_eta          = m_histos.TH1orTH2( "", "2d_dpt_eta"        , "", None, **dPtEta)
pixTk_pixHits_EtaVsN       = m_histos.TH1orTH2( "", "pixTk_pixHits_EtaVsN"      , "", None, **EtaVsN)
pixTk_candPixHits_EtaVsN   = m_histos.TH1orTH2( "", "pixTk_candPixHits_EtaVsN"  , "", None, **EtaVsNCand)
pixTk_sharedPixHits_EtaVsN = m_histos.TH1orTH2( "", "pixTk_sharedPixHits_EtaVsN", "", None, **EtaVsN)

### Resolution Vs pT histograms
resolutionVsPt_pt      = m_histos.TH1orTH2( "", "resolutionVsPt_pt"  , "inclusive", None, **PtPtRes)
resolutionVsPt_pt_C    = m_histos.TH1orTH2( "", "resolutionVsPt_pt_C", C          , None, **PtPtRes)
resolutionVsPt_pt_I    = m_histos.TH1orTH2( "", "resolutionVsPt_pt_I", I          , None, **PtPtRes)
resolutionVsPt_pt_F    = m_histos.TH1orTH2( "", "resolutionVsPt_pt_F", F          , None, **PtPtRes)

resolutionVsPt_eta     = m_histos.TH1orTH2( "", "resolutionVsPt_eta"  , "inclusive", None, **PtEtaRes)
resolutionVsPt_eta_C   = m_histos.TH1orTH2( "", "resolutionVsPt_eta_C",  C         , None, **PtEtaRes)
resolutionVsPt_eta_I   = m_histos.TH1orTH2( "", "resolutionVsPt_eta_I",  I         , None, **PtEtaRes)
resolutionVsPt_eta_F   = m_histos.TH1orTH2( "", "resolutionVsPt_eta_F",  F         , None, **PtEtaRes)

resolutionVsPt_phi     = m_histos.TH1orTH2( "", "resolutionVsPt_phi"  ,  "inclusive", None, **PtPhiRes)
resolutionVsPt_phi_C   = m_histos.TH1orTH2( "", "resolutionVsPt_phi_C",  C          , None, **PtPhiRes)
resolutionVsPt_phi_I   = m_histos.TH1orTH2( "", "resolutionVsPt_phi_I",  I          , None, **PtPhiRes)
resolutionVsPt_phi_F   = m_histos.TH1orTH2( "", "resolutionVsPt_phi_F",  F          , None, **PtPhiRes)

resolutionVsPt_z0      = m_histos.TH1orTH2( "", "resolutionVsPt_z0"  ,  "inclusive", None, **PtZ0Res)
resolutionVsPt_z0_C    = m_histos.TH1orTH2( "", "resolutionVsPt_z0_C",  C          , None, **PtZ0Res)
resolutionVsPt_z0_I    = m_histos.TH1orTH2( "", "resolutionVsPt_z0_I",  I          , None, **PtZ0Res)
resolutionVsPt_z0_F    = m_histos.TH1orTH2( "", "resolutionVsPt_z0_F",  F          , None, **PtZ0Res)

resolutionVsPt_d0      = m_histos.TH1orTH2( "", "resolutionVsPt_d0"  ,  "inclusive", None, **Ptd0Res)
resolutionVsPt_d0_C    = m_histos.TH1orTH2( "", "resolutionVsPt_d0_C",  C          , None, **Ptd0Res)
resolutionVsPt_d0_I    = m_histos.TH1orTH2( "", "resolutionVsPt_d0_I",  I          , None, **Ptd0Res)
resolutionVsPt_d0_F    = m_histos.TH1orTH2( "", "resolutionVsPt_d0_F",  F          , None, **Ptd0Res)


### Resolution Vs eta histograms
resolutionVsEta_pt      = m_histos.TH1orTH2( "", "resolutionVsEta_pt"  , "inclusive", None, **EtaPtRes)
resolutionVsEta_pt_L    = m_histos.TH1orTH2( "", "resolutionVsEta_pt_L", L          , None, **EtaPtRes)
resolutionVsEta_pt_M    = m_histos.TH1orTH2( "", "resolutionVsEta_pt_M", M          , None, **EtaPtRes)
resolutionVsEta_pt_H    = m_histos.TH1orTH2( "", "resolutionVsEta_pt_H", H          , None, **EtaPtRes)

resolutionVsEta_eta     = m_histos.TH1orTH2( "", "resolutionVsEta_eta"  , "inclusive", None, **EtaEtaRes)
resolutionVsEta_eta_L   = m_histos.TH1orTH2( "", "resolutionVsEta_eta_L", L          , None, **EtaEtaRes)
resolutionVsEta_eta_M   = m_histos.TH1orTH2( "", "resolutionVsEta_eta_M", M          , None, **EtaEtaRes)
resolutionVsEta_eta_H   = m_histos.TH1orTH2( "", "resolutionVsEta_eta_H", H          , None, **EtaEtaRes)

resolutionVsEta_phi     = m_histos.TH1orTH2( "", "resolutionVsEta_phi"  , "inclusive", None, **EtaPhiRes)
resolutionVsEta_phi_L   = m_histos.TH1orTH2( "", "resolutionVsEta_phi_L", L          , None, **EtaPhiRes)
resolutionVsEta_phi_M   = m_histos.TH1orTH2( "", "resolutionVsEta_phi_M", M          , None, **EtaPhiRes)
resolutionVsEta_phi_H   = m_histos.TH1orTH2( "", "resolutionVsEta_phi_H", H          , None, **EtaPhiRes)

resolutionVsEta_z0      = m_histos.TH1orTH2( "", "resolutionVsEta_z0"  , "inclusive", None, **EtaZ0Res)
resolutionVsEta_z0_L    = m_histos.TH1orTH2( "", "resolutionVsEta_z0_L", L          , None, **EtaZ0Res)
resolutionVsEta_z0_M    = m_histos.TH1orTH2( "", "resolutionVsEta_z0_M", M          , None, **EtaZ0Res)
resolutionVsEta_z0_H    = m_histos.TH1orTH2( "", "resolutionVsEta_z0_H", H          , None, **EtaZ0Res)

resolutionVsEta_d0      = m_histos.TH1orTH2( "", "resolutionVsEta_d0"  , "inclusive", None, **Etad0Res)
resolutionVsEta_d0_L    = m_histos.TH1orTH2( "", "resolutionVsEta_d0_L", L          , None, **Etad0Res)
resolutionVsEta_d0_M    = m_histos.TH1orTH2( "", "resolutionVsEta_d0_M", M          , None, **Etad0Res)
resolutionVsEta_d0_H    = m_histos.TH1orTH2( "", "resolutionVsEta_d0_H", H          , None, **Etad0Res)


###############################################################
### Main
###############################################################
def DoPlots(hList, dataset, bColourPalette=False, saveExt=""):

    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.SetBoolUseDatasetAsLegEntry(bColourPalette)

    ### Testing
    p.AddDataset(dataset + ":test"           , datasetPaths_Test[dataset])
    ### Algos
    #p.AddDataset(dataset + ":Default"        , datasetPaths_Default  [dataset])
    #p.AddDataset(dataset + ":4 Layers"       , datasetPaths_AtLeast4[dataset])
    #p.AddDataset(dataset + ":4 (3) Layers"   , datasetPaths_Priority4[dataset])
    #p.AddDataset(dataset + ":4 (3) Recovered", datasetPaths_AtLeast4Not[dataset])

    ### 1 Lost Hit
    #p.AddDataset(dataset + ":1 Lost Hit", datasetPaths_3PixHitsCand[dataset])
    #p.AddDataset(dataset + ":1 Lost Hit, |#eta|<1.0", datasetPaths_3PixHitsCand_EtaLE1p0[dataset])
    #p.AddDataset(dataset + ":1 Lost Hit, 1.0 #leq |#eta| #leq 1.6", datasetPaths_3PixHitsCand_EtaGE1p0_LE1p6[dataset])
    #p.AddDataset(dataset + ":1 Lost Hit, 1.6 #leq |#eta| #leq 2.5", datasetPaths_3PixHitsCand_EtaGE1p6_LE2p5[dataset])

    ### Algos & 2xWindow
    #p.AddDataset(dataset + ":Default, 2x"        , datasetPaths_Default_2xWindow  [dataset])
    #p.AddDataset(dataset + ":4 (3) Layers, 2x"   , datasetPaths_Priority4_2xWindow[dataset])
    #p.AddDataset(dataset + ":4 (3) Recovered, 2x", datasetPaths_AtLeast4Not_2xWindow[dataset])

    ### HP, Default Algo
    #p.AddDataset(dataset + ":Default HP   7"          , datasetPaths_Default_HP_7[dataset])
    #p.AddDataset(dataset + ":Default HP  11"          , datasetPaths_Default_HP_11[dataset])
    #p.AddDataset(dataset + ":Default HP  13"          , datasetPaths_Default_HP_13[dataset])
    #p.AddDataset(dataset + ":Default HP  14"          , datasetPaths_Default_HP_14[dataset])
    #p.AddDataset(dataset + ":Default HP  15"          , datasetPaths_Default_HP_15[dataset])
    #p.AddDataset(dataset + ":Default HP  23"          , datasetPaths_Default_HP_23[dataset])
    #p.AddDataset(dataset + ":Default HP  51"          , datasetPaths_Default_HP_51[dataset])
    #p.AddDataset(dataset + ":Default HP 112"          , datasetPaths_Default_HP_112[dataset])
    #p.AddDataset(dataset + ":Default HP 113"          , datasetPaths_Default_HP_113[dataset])
    #p.AddDataset(dataset + ":Default HP 7, 11, 13"    , datasetPaths_Default_HP_7or11or13[dataset])
    #p.AddDataset(dataset + ":Default HP 7, 11, 13, 14", datasetPaths_Default_HP_7or11or13or14[dataset])

    ### HP, Priority 4
    #p.AddDataset(dataset + ":4 (3) HP   7"          , datasetPaths_Priority4_HP_7[dataset])
    #p.AddDataset(dataset + ":4 (3) HP  11"          , datasetPaths_Priority4_HP_11[dataset])
    #p.AddDataset(dataset + ":4 (3) HP  13"          , datasetPaths_Priority4_HP_13[dataset])
    #p.AddDataset(dataset + ":4 (3) HP  14"          , datasetPaths_Priority4_HP_14[dataset])
    #p.AddDataset(dataset + ":4 (3) HP  15"          , datasetPaths_Priority4_HP_15[dataset])
    #p.AddDataset(dataset + ":4 (3) HP  23"          , datasetPaths_Priority4_HP_23[dataset])
    #p.AddDataset(dataset + ":4 (3) HP  51"          , datasetPaths_Priority4_HP_51[dataset])
    #p.AddDataset(dataset + ":4 (3) HP 112"          , datasetPaths_Priority4_HP_112[dataset])
    #p.AddDataset(dataset + ":4 (3) HP 113"          , datasetPaths_Priority4_HP_113[dataset])
    #p.AddDataset(dataset + ":4 (3) HP 7, 11, 13"    , datasetPaths_Priority4_HP_7or11or13[dataset])
    #p.AddDataset(dataset + ":4 (3) HP 7, 11, 13, 14", datasetPaths_Priority4_HP_7or11or13or14[dataset])

    ### HP, AtLeast4Not
    #p.AddDataset(dataset + ":4 (3)R HP   7"          , datasetPaths_AtLeast4Not_HP_7[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP  11"          , datasetPaths_AtLeast4Not_HP_11[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP  13"          , datasetPaths_AtLeast4Not_HP_13[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP  14"          , datasetPaths_AtLeast4Not_HP_14[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP  15"          , datasetPaths_AtLeast4Not_HP_15[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP  23"          , datasetPaths_AtLeast4Not_HP_23[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP  51"          , datasetPaths_AtLeast4Not_HP_51[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP 112"          , datasetPaths_AtLeast4Not_HP_112[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP 113"          , datasetPaths_AtLeast4Not_HP_113[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP 7, 11, 13"    , datasetPaths_AtLeast4Not_HP_7or11or13[dataset])
    #p.AddDataset(dataset + ":4 (3)R HP 7, 11, 13, 14", datasetPaths_AtLeast4Not_HP_7or11or13or14[dataset])

    
    if (len(p.DatasetToRootFileMap.keys()) > 1):
        p.SetBoolUseDatasetAsLegEntry(True)
    else:
        p.SetBoolUseDatasetAsLegEntry(False)
        
    p.EnableColourPalette(bColourPalette)
    p.AddHisto(hList)
    #p.SetupStatsBox(-1, -1, -1, -1, 000000000)
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
        #DoPlots( [tp_pt_L] , datasetList[0], bColourPalette )
        #DoPlots( [tp_pt_C] , datasetList[0], bColourPalette )
        #DoPlots( [tp_pt_I] , datasetList[0], bColourPalette )
        #DoPlots( [tp_pt_F] , datasetList[0], bColourPalette )
        DoPlots( [tp_eta]  , datasetList[0], bColourPalette )
        #DoPlots( [tp_eta_L], datasetList[0], bColourPalette )
        #DoPlots( [tp_eta_M], datasetList[0], bColourPalette )
        #DoPlots( [tp_eta_H], datasetList[0], bColourPalette )
        #DoPlots( [tp_phi]  , datasetList[0], bColourPalette )
        DoPlots( [tp_z0]   , datasetList[0], bColourPalette )
        DoPlots( [tp_d0]   , datasetList[0], bColourPalette )
        DoPlots( [tp_dxy]  , datasetList[0], bColourPalette )

        
    if(bMatchedTPs):
        DoPlots( [match_tp_pt   ], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_pt_L ], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_pt_C ], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_pt_I ], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_pt_F ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_eta  ], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_eta_L], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_eta_M], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_eta_H], datasetList[0], bColourPalette)
        #DoPlots( [match_tp_phi  ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_z0   ], datasetList[0], bColourPalette)
        DoPlots( [match_tp_d0   ], datasetList[0], bColourPalette)


    if(bMatchedTks):
        DoPlots( [match_trk_chi2_L    , match_trk_chi2_M    , match_trk_chi2_H]    , datasetList[0], bColourPalette, "MH")
        DoPlots( [match_trk_chi2_dof_L, match_trk_chi2_dof_M, match_trk_chi2_dof_H], datasetList[0], bColourPalette, "MH")

        
    if(bResiduals):
        DoPlots( [residual_pt   ], datasetList[0], bColourPalette)
        DoPlots( [residual_eta  ], datasetList[0], bColourPalette)
        DoPlots( [residual_phi  ], datasetList[0], bColourPalette)
        DoPlots( [residual_z0   ], datasetList[0], bColourPalette)
        DoPlots( [residual_d0   ], datasetList[0], bColourPalette)

        
    if(bResiduals_EtaRange):
        DoPlots( [residual_pt_C ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_pt_I ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_pt_F ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_eta_C ]  , datasetList[0], bColourPalette)
        DoPlots( [residual_eta_I ]  , datasetList[0], bColourPalette)
        DoPlots( [residual_eta_F ]  , datasetList[0], bColourPalette)
        DoPlots( [residual_phi_C ]  , datasetList[0], bColourPalette)
        DoPlots( [residual_phi_I ]  , datasetList[0], bColourPalette)
        DoPlots( [residual_phi_F ]  , datasetList[0], bColourPalette)
        DoPlots( [residual_z0_C ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_z0_I ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_z0_F ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_d0_C ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_d0_I ]   , datasetList[0], bColourPalette)
        DoPlots( [residual_d0_F ]   , datasetList[0], bColourPalette)
            

    if(bResiduals_CIF):        
        DoPlots( [residual_pt_C   , residual_pt_I   , residual_pt_F    ], datasetList[0], bColourPalette, "IF")
        DoPlots( [residual_eta_C  , residual_eta_I  , residual_eta_F   ], datasetList[0], bColourPalette, "IF")
        DoPlots( [residual_phi_C  , residual_phi_I  , residual_phi_F   ], datasetList[0], bColourPalette, "IF")
        DoPlots( [residual_z0_C   , residual_z0_I   , residual_z0_F    ], datasetList[0], bColourPalette, "IF")
        DoPlots( [residual_d0_C   , residual_d0_I   , residual_d0_F    ], datasetList[0], bColourPalette, "IF")
        
        
    if(bTkChiSq):
        DoPlots( [match_trk_chi2], datasetList[0], bColourPalette)
    if(bTkChiSq_Slices):        
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
    if(bTkRedChiSq_Slices):
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
       DoPlots( [efficiency_pt]   , datasetList[0], bColourPalette)
       DoPlots( [efficiency_eta ] , datasetList[0], bColourPalette)
    if (bTkEfficiency_Vs_Pt):
        DoPlots( [efficiency_pt_C] , datasetList[0], bColourPalette)
        DoPlots( [efficiency_pt_I] , datasetList[0], bColourPalette)
        DoPlots( [efficiency_pt_F] , datasetList[0], bColourPalette)
        DoPlots( [efficiency_pt_L] , datasetList[0], bColourPalette)
    if (bTkEfficiency_Vs_Pt_CIF):
        DoPlots( [efficiency_pt    , efficiency_pt_C, efficiency_pt_I, efficiency_pt_F]   , datasetList[0], True, "_CIF")
    if (bTkEfficiency_Vs_Eta_LMH):
        DoPlots( [efficiency_eta   , efficiency_eta_L, efficiency_eta_M, efficiency_eta_H], datasetList[0], True, "_LMH")
    if (bTkEfficiency_Vs_Eta):
        DoPlots( [efficiency_eta_L], datasetList[0], bColourPalette)
        DoPlots( [efficiency_eta_M], datasetList[0], bColourPalette)
        DoPlots( [efficiency_eta_H], datasetList[0], bColourPalette)
    if (bTkEfficiency_Vs_d0z0phi):
        DoPlots( [efficiency_phi]  , datasetList[0], bColourPalette)
        DoPlots( [efficiency_z0 ]  , datasetList[0], bColourPalette)
        DoPlots( [efficiency_d0]   , datasetList[0], bColourPalette)
                            

    if(bDo2D):
        #DoPlots( [ logChiSq_eta ]     , datasetList[0], bColourPalette)
        #DoPlots( [ logChiSq_dof_eta ] , datasetList[0], bColourPalette)
        #DoPlots( [ dz0_eta ]          , datasetList[0], bColourPalette)
        #DoPlots( [ deta_eta ]         , datasetList[0], bColourPalette)
        #DoPlots( [ dphi_eta ]         , datasetList[0], bColourPalette)
        #DoPlots( [ dpt_eta ]          , datasetList[0], bColourPalette)
        DoPlots( [ pixTk_pixHits_EtaVsN      ], datasetList[0], bColourPalette)
        DoPlots( [ pixTk_candPixHits_EtaVsN  ], datasetList[0], bColourPalette)
        DoPlots( [ pixTk_sharedPixHits_EtaVsN], datasetList[0], bColourPalette)


    ### Pt Resolutions Vs Pt
    if(pt_Resolution_Vs_Pt):
        DoPlots( [ resolutionVsPt_pt  ], datasetList[0], bColourPalette)
    if(pt_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resolutionVsPt_pt_C], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_pt_I], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_pt_F], datasetList[0], bColourPalette)
    if(pt_Resolution_Vs_Pt_CIF):
        DoPlots( [ resolutionVsPt_pt_C, resolutionVsPt_pt_I, resolutionVsPt_pt_F], datasetList[0], True, "IF")
    ### Pt Residuals Vs Pt
    if(pt_Residual_Vs_Pt_PtRange or pt_Residual_Vs_Pt_EtaRange or pt_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(pt_Residual_Vs_Pt_PtRange):
                DoPlots( h_residualVsPt_pt[i], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_residualVsPt_pt_C[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_pt_I[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_pt_F[i]], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_residualVsPt_pt_C[i], h_residualVsPt_pt_I[i], h_residualVsPt_pt_F[i]], datasetList[0], bColourPalette)
    ### Pt Resolutions Vs Eta
    if(pt_Resolution_Vs_Eta):
        DoPlots( [ resolutionVsEta_pt  ], datasetList[0], bColourPalette)
    if (pt_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resolutionVsEta_pt_L ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_pt_M ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_pt_H ], datasetList[0], bColourPalette)
    if (pt_Resolution_Vs_Eta_LMH):
        DoPlots( [ resolutionVsEta_pt_L, resolutionVsEta_pt_M, resolutionVsEta_pt_H ], datasetList[0], True, "MH")
    ### Pt Residuals Vs Eta
    if(pt_Residual_Vs_Eta_EtaRange or pt_Residual_Vs_Eta_PtRange_EtaRange or pt_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(pt_Residual_Vs_Eta_EtaRange):
                DoPlots( h_residualVsEta_pt[i], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_residualVsEta_pt_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_pt_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_pt_H[i] ], datasetList[0], bColourPalette)
            if (pt_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_residualVsEta_pt_H[i], h_residualVsEta_pt_M[i], h_residualVsEta_pt_L[i] ], datasetList[0], True, "MH")
            
    ### Eta Resolutions Vs Pt
    if(eta_Resolution_Vs_Pt):
        DoPlots( [ resolutionVsPt_eta  ], datasetList[0], bColourPalette)
    if(eta_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resolutionVsPt_eta_C], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_eta_I], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_eta_F], datasetList[0], bColourPalette)
    if(eta_Resolution_Vs_Pt_CIF):
        DoPlots( [ resolutionVsPt_eta_C, resolutionVsPt_eta_I, resolutionVsPt_eta_F], datasetList[0], True, "IF")

    ### Eta Residuals Vs Pt
    if(eta_Residual_Vs_Pt_PtRange or eta_Residual_Vs_Pt_EtaRange or eta_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(eta_Residual_Vs_Pt_PtRange):
                DoPlots( h_residualVsPt_eta[i], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_residualVsPt_eta_C[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_eta_I[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_eta_F[i]], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_residualVsPt_eta_C[i], h_residualVsPt_eta_I[i], h_residualVsPt_eta_F[i]], datasetList[0], bColourPalette)
    ### Eta Resolutions Vs Eta
    if(eta_Resolution_Vs_Eta):
        DoPlots( [ resolutionVsEta_eta  ], datasetList[0], bColourPalette)
    if (eta_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resolutionVsEta_eta_L ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_eta_M ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_eta_H ], datasetList[0], bColourPalette)
    if (eta_Resolution_Vs_Eta_LMH):
        DoPlots( [ resolutionVsEta_eta_H, resolutionVsEta_eta_M, resolutionVsEta_eta_L ], datasetList[0], True, "ML")
    ### Eta Residuals Vs Eta
    if(eta_Residual_Vs_Eta_EtaRange or eta_Residual_Vs_Eta_PtRange_EtaRange or eta_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(eta_Residual_Vs_Eta_EtaRange):
                DoPlots( h_residualVsEta_eta[i], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_residualVsEta_eta_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_eta_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_eta_H[i] ], datasetList[0], bColourPalette)
            if (eta_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_residualVsEta_eta_H[i], h_residualVsEta_eta_M[i], h_residualVsEta_eta_L[i] ], datasetList[0], True, "ML")

        
    ### Phi Resolutions Vs Pt
    if(phi_Resolution_Vs_Pt):
        DoPlots( [ resolutionVsPt_phi  ], datasetList[0], bColourPalette)
    if(phi_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resolutionVsPt_phi_C], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_phi_I], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_phi_F], datasetList[0], bColourPalette)
    if(phi_Resolution_Vs_Pt_CIF):
        DoPlots( [ resolutionVsPt_phi_C, resolutionVsPt_phi_I, resolutionVsPt_phi_F], datasetList[0], True, "IF")
    ### Phi Residuals Vs Pt
    if(phi_Residual_Vs_Pt_PtRange or phi_Residual_Vs_Pt_EtaRange or phi_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(phi_Residual_Vs_Pt_PtRange):
                DoPlots( h_residualVsPt_phi[i], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_residualVsPt_phi_C[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_phi_I[i]], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_phi_F[i]], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_residualVsPt_phi_C[i], h_residualVsPt_phi_I[i], h_residualVsPt_phi_F[i]], datasetList[0], bColourPalette)
    ### Phi Resolutions Vs Eta
    if(phi_Resolution_Vs_Eta):
        DoPlots( [ resolutionVsEta_phi  ], datasetList[0], bColourPalette)
    if (phi_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resolutionVsEta_phi_L ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_phi_M ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_phi_H ], datasetList[0], bColourPalette)
    if (phi_Resolution_Vs_Eta_LMH):
        DoPlots( [ resolutionVsEta_phi_H, resolutionVsEta_phi_M, resolutionVsEta_phi_L ], datasetList[0], True, "ML")
    ### Phi Residuals Vs Eta
    if(phi_Residual_Vs_Eta_EtaRange or phi_Residual_Vs_Eta_PtRange_EtaRange or phi_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(phi_Residual_Vs_Eta_EtaRange):
                DoPlots( h_residualVsEta_phi[i], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_residualVsEta_phi_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_phi_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_phi_H[i] ], datasetList[0], bColourPalette)
            if (phi_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_residualVsEta_phi_H[i], h_residualVsEta_phi_M[i], h_residualVsEta_phi_L[i] ], datasetList[0], True, "ML")
                
                
    ### z0 Resolutions Vs Pt
    if(z0_Resolution_Vs_Pt):
        DoPlots( [ resolutionVsPt_z0  ], datasetList[0], bColourPalette)
    if(z0_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resolutionVsPt_z0_C ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_z0_I ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_z0_F ], datasetList[0], bColourPalette)
    if(z0_Resolution_Vs_Pt_CIF):
        DoPlots( [ resolutionVsPt_z0_C, resolutionVsPt_z0_I, resolutionVsPt_z0_F], datasetList[0], True, "IF")
    ### z0 Residuals Vs Pt
    if(z0_Residual_Vs_Pt_PtRange or z0_Residual_Vs_Pt_EtaRange or z0_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(z0_Residual_Vs_Pt_PtRange):
                DoPlots( h_residualVsPt_z0[i], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_residualVsPt_z0_C[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_z0_I[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_z0_F[i] ], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_residualVsPt_z0_C[i], h_residualVsPt_z0_I[i], h_residualVsPt_z0_F[i]], datasetList[0], bColourPalette)
    ### z0 Resolutions Vs Eta
    if(z0_Resolution_Vs_Eta):
        DoPlots( [ resolutionVsEta_z0  ], datasetList[0], bColourPalette)
    if (z0_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resolutionVsEta_z0_L ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_z0_M ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_z0_H ], datasetList[0], bColourPalette)
    if (z0_Resolution_Vs_Eta_LMH):
        DoPlots( [ resolutionVsEta_z0_H, resolutionVsEta_z0_M, resolutionVsEta_z0_L ], datasetList[0], True, "ML")
    ### z0 Residuals Vs Eta
    if(z0_Residual_Vs_Eta_EtaRange or z0_Residual_Vs_Eta_PtRange_EtaRange or z0_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(z0_Residual_Vs_Eta_EtaRange):
                DoPlots( h_residualVsEta_z0[i], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_residualVsEta_z0_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_z0_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_z0_H[i] ], datasetList[0], bColourPalette)
            if (z0_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_residualVsEta_z0_H[i], h_residualVsEta_z0_M[i], h_residualVsEta_z0_L[i] ], datasetList[0], True, "ML")

                
    ### d0 Resolutions Vs Pt
    if(d0_Resolution_Vs_Pt):
        DoPlots( [ resolutionVsPt_d0  ], datasetList[0], bColourPalette)
    if(d0_Resolution_Vs_Pt_EtaRange):
        DoPlots( [ resolutionVsPt_d0_C ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_d0_I ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsPt_d0_F ], datasetList[0], bColourPalette)
    if(d0_Resolution_Vs_Pt_CIF):
        DoPlots( [ resolutionVsPt_d0_C, resolutionVsPt_d0_I, resolutionVsPt_d0_F], datasetList[0], True, "IF")
    ### d0 Residuals Vs Pt
    if(d0_Residual_Vs_Pt_PtRange or d0_Residual_Vs_Pt_EtaRange or d0_Residual_Vs_Pt_CIF_PtRange):
        for i in range(0, nPt_Range):
            if(d0_Residual_Vs_Pt_PtRange):
                DoPlots( h_residualVsPt_d0[i], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Pt_EtaRange):
                DoPlots( [ h_residualVsPt_d0_C[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_d0_I[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsPt_d0_F[i] ], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Pt_CIF_PtRange):
                DoPlots( [ h_residualVsPt_d0_C[i], h_residualVsPt_d0_I[i], h_residualVsPt_d0_F[i]], datasetList[0], bColourPalette)
    ### d0 Resolutions Vs Eta
    if(d0_Resolution_Vs_Eta):
        DoPlots( [ resolutionVsEta_d0  ], datasetList[0], bColourPalette)
    if (d0_Resolution_Vs_Eta_PtRange):
        DoPlots( [ resolutionVsEta_d0_L ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_d0_M ], datasetList[0], bColourPalette)
        DoPlots( [ resolutionVsEta_d0_H ], datasetList[0], bColourPalette)
    if (d0_Resolution_Vs_Eta_LMH):
        DoPlots( [ resolutionVsEta_d0_H, resolutionVsEta_d0_M, resolutionVsEta_d0_L ], datasetList[0], True, "ML")
    ### d0 Residuals Vs Eta
    if(d0_Residual_Vs_Eta_EtaRange or d0_Residual_Vs_Eta_PtRange_EtaRange or d0_Residual_Vs_Eta_LMH_EtaRange):
        for i in range(0, nEta_Range):
            if(d0_Residual_Vs_Eta_EtaRange):
                DoPlots( h_residualVsEta_d0[i], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Eta_PtRange_EtaRange):
                DoPlots( [ h_residualVsEta_d0_L[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_d0_M[i] ], datasetList[0], bColourPalette)
                DoPlots( [ h_residualVsEta_d0_H[i] ], datasetList[0], bColourPalette)
            if (d0_Residual_Vs_Eta_LMH_EtaRange):
                DoPlots( [ h_residualVsEta_d0_H[i], h_residualVsEta_d0_M[i], h_residualVsEta_d0_L[i] ], datasetList[0], True, "ML")

    ### Pixel Hits (best-fit)
    if(bDo_Pixel_FitHits):
        DoPlots( [ pixTk_pixHits_N        ], datasetList[0], True)
        #DoPlots( [ pixTk_pixHits_Rho      ], datasetList[0], True)
        #DoPlots( [ pixTk_pixHits_Z        ], datasetList[0], True)
        DoPlots( [ pixTk_pixHits_Type     ], datasetList[0], True)
        DoPlots( [ pixTk_pixHits_Pattern  ], datasetList[0], True)
        DoPlots( [ pixTk_pixHits_PatternHL], datasetList[0], True, "HL")

    if(bDo_Pixel_FitHits_2D):
        #DoPlots( [ pixTk_pixHits_ZVsRho ]       , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_ZVsRho_Norm ]  , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_ZVsRho_C ]     , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_ZVsRho_I ]     , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_ZVsRho_F ]     , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_ZVsRho_EtaGE2 ], datasetList[0], True)
        #DoPlots( [ pixTk_pixHits_XVsY ]         , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_XVsY_Norm ]    , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_XVsY_C ]       , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_XVsY_I ]       , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_XVsY_F ]       , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_XVsY_EtaGE2 ]  , datasetList[0], True)
        DoPlots( [ pixTk_pixHits_PatternVsEta]  , datasetList[0], True)

    ### Pixel Hits (shared)
    if(bDo_Pixel_SharedFitHits):
        DoPlots( [ pixTk_sharedPixHits_N      ], datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_Rho    ], datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_Z      ], datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_Type   ], datasetList[0], True)
    if(bDo_Pixel_SharedFitHits_2D):
        DoPlots( [ pixTk_sharedPixHits_ZVsRho ]       , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_ZVsRho_C ]     , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_ZVsRho_I ]     , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_ZVsRho_F ]     , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_ZVsRho_EtaGE2 ], datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_XVsY ]         , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_XVsY_C ]       , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_XVsY_I ]       , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_XVsY_F ]       , datasetList[0], True)
        DoPlots( [ pixTk_sharedPixHits_XVsY_EtaGE2 ]  , datasetList[0], True)

    ### Pixel Hits (candidates)
    if(bDo_CandPixel_FitHits):
        DoPlots( [ pixTk_candPixHits_N      ], datasetList[0], True)
        #DoPlots( [ pixTk_candPixHits_Rho    ], datasetList[0], True)
        #DoPlots( [ pixTk_candPixHits_Z      ], datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_Type   ], datasetList[0], True)
    if(bDo_CandPixel_FitHits_2D):
        DoPlots( [ pixTk_candPixHits_ZVsRho ]       , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_ZVsRho_C ]     , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_ZVsRho_I ]     , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_ZVsRho_F ]     , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_ZVsRho_EtaGE2 ], datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_XVsY ]         , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_XVsY_C ]       , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_XVsY_I ]       , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_XVsY_F ]       , datasetList[0], True)
        DoPlots( [ pixTk_candPixHits_XVsY_EtaGE2 ]  , datasetList[0], True)
        

############################################################################################################################################

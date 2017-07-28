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
bDoDistributions = False
bDoEfficiencies  = True
tkWP             = "MediumWP" #"MediumWP" #"vLooseWP"

saveFormats = ["png", "pdf"]
datasetList = ["VBF"]
inputPath   = "Macros/TauTrigger/Validation_Histograms_"
savePath    = ""
#inputPath   = "Macros/TauTrigger/results/validation/5FitParams/Validation_Histograms_"
#savePath    = "/Users/attikis/talks/TauTrigger_06March2015/figures/5FitParams/pixVstt/" + datasetList[0] + "/"


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
### Cuts
############################################################### 
if tkWP == "vLooseWP": 
    cut_pt            = [2.0]
    cut_eta           = [-3.0, 3.0]
    cut_pocaz         = [30.0 ]
    cut_chiSq         = [200.0]
    cut_redChiSq      = []
    cut_stubPtCons    = []
    cut_d0Abs         = []
    cut_nStubs        = [4.0]
    cut_nPSStubs      = []
    cutDir_pt         = [[0.0, 2.0, ROOT.kBlack]]
    cutDir_eta        = []
    cutDir_pocaz      = [[-30.5, -30.0, ROOT.kBlack], [+30.0, +30.5, ROOT.kBlack]]
    cutDir_chiSq      = [[+200.0, +300.0, ROOT.kBlack]]
    cutDir_redChiSq   = []
    cutDir_stubPtCons = []
    cutDir_d0Abs      = []
    cutDir_nStubs     = [[0.0, 4.0, ROOT.kBlack]]
    cutDir_nPSStubs   = []
elif (tkWP == "MediumWP"): 
    cut_pt            = [5.0]
    cut_eta           = []
    cut_pocaz         = []
    cut_chiSq         = []
    cut_redChiSq      = [10]
    cut_stubPtCons    = []
    cut_d0Abs         = []
    cut_nStubs        = [5]
    cut_nPSStubs      = [3] #[3, 7]
    cutDir_pt         = [[0.0, 5.0, ROOT.kBlack]]
    cutDir_eta        = []
    cutDir_pocaz      = []
    cutDir_chiSq      = []
    cutDir_redChiSq   = [[10.0, 50.0, ROOT.kBlack]]
    cutDir_stubPtCons = []
    cutDir_d0Abs      = []
    cutDir_nStubs     = [[0.0, 5.0, ROOT.kBlack]]
    #cutDir_nPSStubs   = [[0.0, 3.0, ROOT.kBlack], [7.0, 10.5, ROOT.kBlack]]
    cutDir_nPSStubs   = [[0.0, 3.0, ROOT.kBlack]]
else:
    exit(1)
############################################################### 
### Histogram Options
############################################################### 
yMinRatio  = 0.0
yMaxRatio  = 5.0
bRatio     = False
bInvRatio  = False
ratioLabel = "ratio"
normFactor = "One"

Pt = {
    "xLabel": "p_{T}"           , "xUnits": "GeV/c", "xMin": 0.00 , "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": cut_pt, "xCutBoxes": cutDir_pt, "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": 1e-04 , "yMax": +1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

Eta = {
    "xLabel": "#eta"            , "xUnits": ""     , "xMin": -3.5 , "xMax": +3.5 , "binWidthX": 0.1 , "xCutLines": cut_eta, "xCutBoxes": cutDir_eta, "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1e-04, "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

EtaAbs = {
    "xLabel": "|#eta|"          , "xUnits": ""     , "xMin": +0.0  , "xMax": +3.0 , "binWidthX": 0.1 , "xCutLines": cut_eta, "xCutBoxes": cutDir_eta, "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1e-04 , "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

POCAz = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": 0.5 , "xCutLines": cut_pocaz, "xCutBoxes": cutDir_pocaz, "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1e-04, "yMax":  +1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

POCAzAbs = {
    "xLabel": "|z_{0}|"         , "xUnits": "cm", "xMin": +0.0  , "xMax": +30.0, "binWidthX": 0.5 , "xCutLines": cut_pocaz, "xCutBoxes": cutDir_pocaz, "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1e-04 , "yMax":  +1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

ChiSq = {
    "xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +300.0, "binWidthX": 2.0 , "xCutLines": cut_chiSq, "xCutBoxes": cutDir_chiSq, "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""    , "yMin": 1e-04 , "yMax":  +1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

RedChiSq = {
    "xLabel": "#chi^{2}_{N}"   , "xUnits": "", "xMin": +0.0 , "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": cut_redChiSq, "xCutBoxes": cutDir_redChiSq, "gridX": True, "logX": False,
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-04, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

StubPtCons = {
    "xLabel": "p_{T}^{stub} constistency", "xUnits": "", "xMin": 0.0  , "xMax": 200.0, "binWidthX": 5.0 , "xCutLines": cut_stubPtCons, "xCutBoxes": cutDir_stubPtCons, "gridX": True, "logX": False,
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-04 , "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

d0Abs = {
    "xLabel": "|d_{0}|"        , "xUnits": "cm", "xMin": +0.0 , "xMax": +1.0, "binWidthX": 0.05 , "xCutLines": cut_d0Abs, "xCutBoxes": cutDir_d0Abs, "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1e-04, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

NStubs = {
    "xLabel": "Stubs Multiplicity", "xUnits": "", "xMin": +0.0, "xMax": +10.0, "binWidthX": 1.0 , "xCutLines": cut_nStubs, "xCutBoxes": cutDir_nStubs, "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": 1e-04, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }

NPsStubs = {
    "xLabel": "PS Stubs Multiplicity", "xUnits": "", "xMin": +0.0, "xMax": 10.0, "binWidthX": 1.0 , "xCutLines": cut_nPSStubs, "xCutBoxes": cutDir_nPSStubs, "gridX": True, "logX": False,
    "yLabel": "Entries / %0.0f"      , "yUnits": "", "yMin": 1e-4, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.65, "yLegMax": 0.75 }


############################################################### 
### Group Histograms: L1Tks
############################################################### 
histoFolder    = ""
hL1Tks_Pt                = m_histos.TH1orTH2( histoFolder, "L1Tks_Pt"               , "TTTracks (5)", None, **Pt         )
hL1Tks_Eta               = m_histos.TH1orTH2( histoFolder, "L1Tks_Eta"              , "TTTracks (5)", None, **Eta        )
hL1Tks_EtaAbs            = m_histos.TH1orTH2( histoFolder, "L1Tks_EtaAbs"           , "TTTracks (5)", None, **EtaAbs     )
hL1Tks_POCAz             = m_histos.TH1orTH2( histoFolder, "L1Tks_POCAz"            , "TTTracks (5)", None, **POCAz      )
hL1Tks_POCAzAbs          = m_histos.TH1orTH2( histoFolder, "L1Tks_POCAzAbs"         , "TTTracks (5)", None, **POCAzAbs   )
hL1Tks_ChiSquared        = m_histos.TH1orTH2( histoFolder, "L1Tks_ChiSquared"       , "TTTracks (5)", None, **ChiSq      )
hL1Tks_RedChiSquared     = m_histos.TH1orTH2( histoFolder, "L1Tks_RedChiSquared"    , "TTTracks (5)", None, **RedChiSq   )
hL1Tks_StubPtconsistency = m_histos.TH1orTH2( histoFolder, "L1Tks_StubPtConsistency", "TTTracks (5)", None, **StubPtCons )
hL1Tks_d0Abs             = m_histos.TH1orTH2( histoFolder, "L1Tks_d0Abs"            , "TTTracks (5)", None, **d0Abs      )
hL1Tks_NStubs            = m_histos.TH1orTH2( histoFolder, "L1Tks_NStubs"           , "TTTracks (5)", None, **NStubs     )
hL1Tks_NPsStubs          = m_histos.TH1orTH2( histoFolder, "L1Tks_NPsStubs"         , "TTTracks (5)", None, **NPsStubs   )

############################################################### 
### Group Histograms: L1PixTks
############################################################### 
hL1PixTks_Pt                = m_histos.TH1orTH2( histoFolder, "L1PixTks_Pt"               , "TTPixelTracks", None, **Pt         )
hL1PixTks_Eta               = m_histos.TH1orTH2( histoFolder, "L1PixTks_Eta"              , "TTPixelTracks", None, **Eta        )
hL1PixTks_EtaAbs            = m_histos.TH1orTH2( histoFolder, "L1PixTks_EtaAbs"           , "TTPixelTracks", None, **EtaAbs     )
hL1PixTks_POCAz             = m_histos.TH1orTH2( histoFolder, "L1PixTks_POCAz"            , "TTPixelTracks", None, **POCAz      )
hL1PixTks_POCAzAbs          = m_histos.TH1orTH2( histoFolder, "L1PixTks_POCAzAbs"         , "TTPixelTracks", None, **POCAzAbs   )
hL1PixTks_ChiSquared        = m_histos.TH1orTH2( histoFolder, "L1PixTks_ChiSquared"       , "TTPixelTracks", None, **ChiSq      )
hL1PixTks_RedChiSquared     = m_histos.TH1orTH2( histoFolder, "L1PixTks_RedChiSquared"    , "TTPixelTracks", None, **RedChiSq   )
hL1PixTks_StubPtconsistency = m_histos.TH1orTH2( histoFolder, "L1PixTks_StubPtConsistency", "TTPixelTracks", None, **StubPtCons )
hL1PixTks_d0Abs             = m_histos.TH1orTH2( histoFolder, "L1PixTks_d0Abs"            , "TTPixelTracks", None, **d0Abs      )
hL1PixTks_NStubs            = m_histos.TH1orTH2( histoFolder, "L1PixTks_NStubs"           , "TTPixelTracks", None, **NStubs     )
hL1PixTks_NPsStubs          = m_histos.TH1orTH2( histoFolder, "L1PixTks_NPsStubs"         , "TTPixelTracks", None, **NPsStubs   )

###############################################################
### Main
###############################################################
def DoPlots(hList, datasetList):

    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
        p.AddHisto(hList)
        p.EnableColourPalette(True)
        p.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
        p.GetTLegend().SetHeader( m_datasets.DatasetClass().ConvertDatasetToLatex(dataset) )
        p.SaveHistos(True, savePath, saveFormats, "_" + tkWP)

    return

###############################################################
def DoEfficiency(hList, datasetList, cutDirList, saveExt=""):

    for cutDir in cutDirList:
        
        p = m_plotter.Plotter( Verbose=False, BatchMode=True )
        for dataset in datasetList:
            p.AddDataset(dataset, datasetPaths[dataset])
        p.AddHisto(hList)        
        p.EnableColourPalette(True)
        p.DrawEfficiency(cutDir, "binomial")
        p.GetTLegend().SetHeader( m_datasets.DatasetClass().ConvertDatasetToLatex(dataset) )
        p.SaveHistos(True, savePath, saveFormats, "_" + tkWP)

    return

###############################################################
if __name__ == "__main__":

    if (bDoDistributions):
        DoPlots( [hL1Tks_Pt               , hL1PixTks_Pt]               , datasetList )
        DoPlots( [hL1Tks_Eta              , hL1PixTks_Eta]              , datasetList )
        DoPlots( [hL1Tks_EtaAbs           , hL1PixTks_EtaAbs]           , datasetList )
        DoPlots( [hL1Tks_POCAz            , hL1PixTks_POCAz]            , datasetList )
        DoPlots( [hL1Tks_POCAzAbs         , hL1PixTks_POCAzAbs]         , datasetList )
        DoPlots( [hL1Tks_ChiSquared       , hL1PixTks_ChiSquared]       , datasetList )
        DoPlots( [hL1Tks_RedChiSquared    , hL1PixTks_RedChiSquared]    , datasetList )
        DoPlots( [hL1Tks_StubPtconsistency, hL1PixTks_StubPtconsistency], datasetList )
        DoPlots( [hL1Tks_d0Abs            , hL1PixTks_d0Abs]            , datasetList )
        DoPlots( [hL1Tks_NStubs           , hL1PixTks_NStubs]           , datasetList )
        DoPlots( [hL1Tks_NPsStubs         , hL1PixTks_NPsStubs]         , datasetList )

    if (bDoEfficiencies):
        DoEfficiency( [hL1Tks_Pt               , hL1PixTks_Pt]               , datasetList, [">"] )
        DoEfficiency( [hL1Tks_EtaAbs           , hL1PixTks_EtaAbs]           , datasetList, ["<"] )
        DoEfficiency( [hL1Tks_POCAzAbs         , hL1PixTks_POCAzAbs]         , datasetList, ["<"] )
        DoEfficiency( [hL1Tks_ChiSquared       , hL1PixTks_ChiSquared]       , datasetList, ["<"] )
        DoEfficiency( [hL1Tks_RedChiSquared    , hL1PixTks_RedChiSquared]    , datasetList, ["<"] )
        DoEfficiency( [hL1Tks_StubPtconsistency, hL1PixTks_StubPtconsistency], datasetList, [">"] )
        DoEfficiency( [hL1Tks_d0Abs            , hL1PixTks_d0Abs]            , datasetList, [">"] )
        DoEfficiency( [hL1Tks_NStubs           , hL1PixTks_NStubs]           , datasetList, [">"] )
        DoEfficiency( [hL1Tks_NPsStubs         , hL1PixTks_NPsStubs]         , datasetList, [">"] )
        

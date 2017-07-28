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
track_WP      = "MediumWP"
datasetList   = ["HPlus160"]
inputPath     = "Macros/TauTrigger/Tracks_Histograms_"
#savePath      = "/Users/attikis/talks/TauTrigger_06March2015/figures/5FitParams/pixVstt/" + datasetList[0] + "/"
savePath      = ""
saveFormats   = ["png", "pdf"]


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
bRatio     = True
bInvRatio  = False
yMin       = 1e-5
yMax       = 1.0
yMinRatio  = 0.0
yMaxRatio  = 2.0
normFactor = "One"

NTks = {
    "xLabel": "Multiplicity"    , "xUnits": "" , "xMin": -0.5 , "xMax": 299.5, "binWidthX": 5.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin , "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Pt = {
    "xLabel": "p_{T}"           , "xUnits": "GeV/c", "xMin": 0.00 , "xMax": 100.0, "binWidthX": 2.0 , "xCutLines": [15], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": 1E-05, "yMax": yMax , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Eta = {
    "xLabel": "#eta"            , "xUnits": ""     , "xMin": -3.5 , "xMax": +3.5 , "binWidthX": 0.1 , "xCutLines": [-3.0, 0, +3.0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1E-03, "yMax": 2E-1 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

EtaAbs = {
    "xLabel": "|#eta|"          , "xUnits": ""     , "xMin": 0.0  , "xMax": +3.5 , "binWidthX": 0.1 , "xCutLines": [0, +3.0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1E-03, "yMax": 2E-1 , "binWidthY": None, "yCutLines": []       , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Phi = {
    "xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.4 , "xMax": +3.4 , "binWidthX": 0.1, "xCutLines": [-pi, 0, +pi], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f" , "yUnits": ""       , "yMin": 5E-3 , "yMax": 1E-1 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

POCAz = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": 0.5 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1E-05 , "yMax": 5E-01 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

POCAzAbs = {
    "xLabel": "|z_{0}|"         , "xUnits": "cm", "xMin": 0.0  , "xMax": +30.0, "binWidthX": 0.5 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1E-05, "yMax": 5E-01 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

ChiSq  = {
    "xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +150.0, "binWidthX": 2.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": ""    , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

RedChiSq = {
    "xLabel": "#chi^{2}/d.o.f.", "xUnits": "", "xMin": -0.5, "xMax": +50.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

RInv = {
    "xLabel": "#rho^{-1}"       , "xUnits": "cm^{-1}", "xMin": -0.006, "xMax": +0.006, "binWidthX": 0.0001, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.4f" , "yUnits": ""       , "yMin": 1E-04, "yMax": 5E-01, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

Charge = {
    "xLabel": "Charge (e)"     , "xUnits": "", "xMin": -3.5 , "xMax": 3.5 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin , "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.26, "xLegMax": 0.45, "yLegMin": 0.72, "yLegMax": 0.82 }

d0 = {
    "xLabel": "d_{0}"   , "xUnits": "cm", "xMin": -1.0 , "xMax": +1.0 , "binWidthX": 0.05 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": "", "yMin": 1E-05, "yMax": yMax, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

d0Abs = {
    "xLabel": "|d_{0}|" , "xUnits": "cm", "xMin": +0.0 , "xMax": +1.0, "binWidthX": 0.05 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": "", "yMin": 1E-05, "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

d0Sig = {
    "xLabel": "d_{0}/#sigma(d_{0})", "xUnits": "", "xMin": -20.0 , "xMax": +20.0, "binWidthX": 0.50, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f"    , "yUnits": "", "yMin": 1e-03, "yMax": 1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

NPix = {
    "xLabel": "Pixel Hits", "xUnits": "", "xMin": 1.0, "xMax": 6.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": 1e-2, "yMax":  1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

NStubs = {
    "xLabel": "Stubs Multiplicity", "xUnits": "", "xMin": -0.5, "xMax": 12.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": 1e-4, "yMax":  yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

NPsStubs = {
    "xLabel": "PS-Stubs Multiplicity", "xUnits": "", "xMin": -0.5, "xMax": 12.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": 1e-4, "yMax":  yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

StubPtCons = {
    "xLabel": "p_{T}^{stub} constistency", "xUnits": "", "xMin": 0.0  , "xMax": 50.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-3 , "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": "Ratio", "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

sigmaRInv = {
    "xLabel": "#sigma_{#rho^{-1}}", "xUnits": "", "xMin": 0.0 , "xMax": 2e-04, "binWidthX": None , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f"   , "yUnits": "", "yMin": 1E-05, "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

sigmaPhi0 = {
    "xLabel": "#sigma_{#phi_{0}}" , "xUnits": "", "xMin": +0.0 , "xMax": +0.0003, "binWidthX": None , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f", "yUnits": "", "yMin": 1E-05 , "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

sigmaD0 = {
    "xLabel": "#sigma_{d_{0}}" , "xUnits": "cm", "xMin": +0.0 , "xMax": 0.008, "binWidthX": 5E-05 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f", "yUnits": "", "yMin": 1E-05 , "yMax": 5E-01, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

sigmaT  = {
    "xLabel": "#sigma_{tan(#lambda)}", "xUnits": "", "xMin": 0.0  , "xMax": 1.5e-3, "binWidthX": None , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f"      , "yUnits": "", "yMin": 1E-05, "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

sigmaZ0  = {
    "xLabel": "#sigma_{z_{0}}" , "xUnits": "", "xMin": 0.0 , "xMax": 1.5e-2, "binWidthX": 1e-4 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f", "yUnits": "", "yMin": 1E-05 , "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

############################################################### 
### Group Histogram (of same binning)
############################################################### 
histoFolder    = ""
hL1PixTks_Multiplicity = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_Multiplicity", "TTPixelTracks", None, **NTks )
hL1PixTks_Pt           = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_Pt"          , "TTPixelTracks", None, **Pt )
hL1PixTks_Eta          = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_Eta"         , "TTPixelTracks", None, **Eta )
hL1PixTks_EtaAbs       = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_EtaAbs"      , "TTPixelTracks", None, **EtaAbs )
hL1PixTks_Phi          = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_Phi"         , "TTPixelTracks", None, **Phi )
hL1PixTks_Charge       = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_Charge"      , "TTPixelTracks", None, **Charge )
hL1PixTks_NPixHits     = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_NPixHits"    , "TTPixelTracks", None, **NPix )
hL1PixTks_NStubs       = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_NStubs"      , "TTPixelTracks", None, **NStubs )
hL1PixTks_NPsStubs     = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_NPsStubs"    , "TTPixelTracks", None, **NPsStubs )
hL1PixTks_StubPtCons   = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_StubPtCons"  , "TTPixelTracks", None, **StubPtCons)
hL1PixTks_POCAz        = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_POCAz"       , "TTPixelTracks", None, **POCAz )
hL1PixTks_POCAzAbs     = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_POCAzAbs"    , "TTPixelTracks", None, **POCAzAbs )
hL1PixTks_d0           = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_d0"          , "TTPixelTracks", None, **d0 )
hL1PixTks_d0Abs        = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_d0Abs"       , "TTPixelTracks", None, **d0Abs )
hL1PixTks_d0Sig        = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_d0Sig"       , "TTPixelTracks", None, **d0Sig )
hL1PixTks_ChiSq        = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_ChiSq"       , "TTPixelTracks", None, **ChiSq )
hL1PixTks_RedChiSq     = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_RedChiSq"    , "TTPixelTracks", None, **RedChiSq )
hL1PixTks_RInv         = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_RInv"        , "TTPixelTracks", None, **RInv )
hL1PixTks_SigmaRInv    = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_SigmaRInv"   , "TTPixelTracks", None, **sigmaRInv )
hL1PixTks_SigmaPhi0    = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_SigmaPhi0"   , "TTPixelTracks", None, **sigmaPhi0 )
hL1PixTks_SigmaD0      = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_SigmaD0"     , "TTPixelTracks", None, **sigmaD0 )
hL1PixTks_SigmaT       = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_SigmaT"      , "TTPixelTracks", None, **sigmaT )
hL1PixTks_SigmaZ0      = m_histos.TH1orTH2( histoFolder, "TTPixelTracks_" + track_WP + "_SigmaZ0"     , "TTPixelTracks", None, **sigmaZ0 )

hL1Tks_Multiplicity    = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_Multiplicity", "TTTracks", None, **NTks )
hL1Tks_Pt              = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_Pt"          , "TTTracks", None, **Pt )
hL1Tks_Eta             = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_Eta"         , "TTTracks", None, **Eta )
hL1Tks_EtaAbs          = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_EtaAbs"      , "TTTracks", None, **EtaAbs )
hL1Tks_Phi             = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_Phi"         , "TTTracks", None, **Phi )
hL1Tks_Charge          = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_Charge"      , "TTTracks", None, **Charge )
hL1Tks_NPixHits        = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_NPixHits"    , "TTTracks", None, **NPix )
hL1Tks_NStubs          = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_NStubs"      , "TTTracks", None, **NStubs )
hL1Tks_NPsStubs        = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_NPsStubs"    , "TTTracks", None, **NPsStubs )
hL1Tks_StubPtCons      = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_StubPtCons"  , "TTTracks", None, **StubPtCons)
hL1Tks_POCAz           = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_POCAz"       , "TTTracks", None, **POCAz )
hL1Tks_POCAzAbs        = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_POCAzAbs"    , "TTTracks", None, **POCAzAbs )
hL1Tks_d0              = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_d0"          , "TTTracks", None, **d0 )
hL1Tks_d0Abs           = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_d0Abs"       , "TTTracks", None, **d0Abs )
hL1Tks_d0Sig           = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_d0Sig"       , "TTTracks", None, **d0Sig )
hL1Tks_ChiSq           = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_ChiSq"       , "TTTracks", None, **ChiSq )
hL1Tks_RedChiSq        = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_RedChiSq"    , "TTTracks", None, **RedChiSq )
hL1Tks_RInv            = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_RInv"        , "TTTracks", None, **RInv )
hL1Tks_SigmaRInv       = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_SigmaRInv"   , "TTTracks", None, **sigmaRInv )
hL1Tks_SigmaPhi0       = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_SigmaPhi0"   , "TTTracks", None, **sigmaPhi0 )
hL1Tks_SigmaD0         = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_SigmaD0"     , "TTTracks", None, **sigmaD0 )
hL1Tks_SigmaT          = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_SigmaT"      , "TTTracks", None, **sigmaT )
hL1Tks_SigmaZ0         = m_histos.TH1orTH2( histoFolder, "TTTracks_" + track_WP + "_SigmaZ0"     , "TTTracks", None, **sigmaZ0 )

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
def DoProfile(hList, datasetList, profileAxis="x"):

    saveExt = ""
    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
    p.SetBoolUseDatasetAsLegEntry(True)
    p.AddHisto(hList)
    p.Draw2DProfile(THStackDrawOpt="nostack", bStackInclusive=False, ProfileAxis=profileAxis, firstBin=1, lastBin=-1)
    p.SetTLegendHeader( "TTPixelTracks", "")
    p.SaveHistos(True, savePath, saveFormats, "_Profile" + profileAxis.upper())

    return
###############################################################
def DoCrossPlots(hList, dataset):

    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.AddDataset( dataset, datasetPaths[dataset])
    p.AddDataset( dataset, datasetPaths2[dataset])
    p.AddHisto(hList)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p.SaveHistos(True, savePath, saveFormats)

    return

###############################################################
if __name__ == "__main__":

    # DoPlots( [hL1PixTks_NPixHits ], datasetList, True )
    # DoPlots( [hL1PixTks_SigmaRInv], datasetList, True )
    # DoPlots( [hL1PixTks_SigmaPhi0], datasetList, True )
    # DoPlots( [hL1PixTks_SigmaD0  ], datasetList, True )
    # DoPlots( [hL1PixTks_SigmaT   ], datasetList, True )
    # DoPlots( [hL1PixTks_SigmaZ0  ], datasetList, True )                    
    DoPlots( [hL1Tks_Multiplicity , hL1PixTks_Multiplicity] , datasetList, True )
    DoPlots( [hL1Tks_Pt           , hL1PixTks_Pt]           , datasetList, True )
    DoPlots( [hL1Tks_Eta          , hL1PixTks_Eta]          , datasetList, True )
    DoPlots( [hL1Tks_EtaAbs       , hL1PixTks_EtaAbs]       , datasetList, True )
    DoPlots( [hL1Tks_Phi          , hL1PixTks_Phi]          , datasetList, True )
    DoPlots( [hL1Tks_Charge       , hL1PixTks_Charge]       , datasetList, True )
    DoPlots( [hL1Tks_NStubs       , hL1PixTks_NStubs]       , datasetList, True )
    DoPlots( [hL1Tks_NPsStubs     , hL1PixTks_NPsStubs]     , datasetList, True )
    DoPlots( [hL1Tks_StubPtCons   , hL1PixTks_StubPtCons]   , datasetList, True )
    DoPlots( [hL1Tks_POCAz        , hL1PixTks_POCAz]        , datasetList, True )
    DoPlots( [hL1Tks_POCAzAbs     , hL1PixTks_POCAzAbs]     , datasetList, True )
    DoPlots( [hL1Tks_d0           , hL1PixTks_d0]           , datasetList, True )
    DoPlots( [hL1Tks_d0Abs        , hL1PixTks_d0Abs]        , datasetList, True )
    DoPlots( [hL1Tks_d0Sig        , hL1PixTks_d0Sig]        , datasetList, True )
    DoPlots( [hL1Tks_ChiSq        , hL1PixTks_ChiSq]        , datasetList, True )
    DoPlots( [hL1Tks_RedChiSq     , hL1PixTks_RedChiSq]     , datasetList, True )
    DoPlots( [hL1Tks_RInv         , hL1PixTks_RInv]         , datasetList, True )
    

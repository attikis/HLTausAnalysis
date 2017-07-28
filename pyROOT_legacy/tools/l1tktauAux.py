###############################################################
### Author .........: Alexandros X. Attikis 
### Institute ......: University Of Cyprus (UCY)
### Email ..........: attikis@cern.ch
###############################################################

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
### Histogram Options
###############################################################
yMinRatio   = 0.0
yMaxRatio   = 2.9
bRatio      = True
bInvRatio   = False
ratioLabel  = "Ratio"
normFactor  = "One"
pi          = 4*math.atan(1)

###############################################################
### Histogram Attributes
###############################################################
kwargs   = {}

NTks = {
    "xLabel": "Multiplicity"    , "xUnits": "" , "xMin": -0.5 , "xMax": 299.5, "binWidthX": 5.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": 1e-05 , "yMax": 1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

PtRel = {
    "xLabel": "p_{T}^{rel}"     , "xUnits": "GeV/c", "xMin": 0.00 , "xMax": +2.0, "binWidthX": 0.05 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 1e-04, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

Pt = {
    "xLabel": "p_{T}"           , "xUnits": "GeV/c", "xMin": 0.00 , "xMax": 60.0, "binWidthX": 1.0 , "xCutLines": [5, 10, 15], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": 1e-03, "yMax": 1.0  , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

Eta = {
    "xLabel": "#eta"            , "xUnits": ""     , "xMin": -2.8 , "xMax": +2.8 , "binWidthX": 0.1 , "xCutLines": [-2.5, 0, +2.5], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": 2E-03, "yMax": 2E-1 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

POCAz = {
    "xLabel": "z_{0}"           , "xUnits": "cm", "xMin": -26.0, "xMax": +26.0, "binWidthX": 0.5 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1E-05 , "yMax": 5E-01 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

dPOCAz = {
    "xLabel": "|#Delta z_{0}|"  , "xUnits": "cm", "xMin": 0.00 , "xMax": +1.0 , "binWidthX": 0.05, "xCutLines": [0.1, 0.5], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": 1e-03, "yMax": +1e-0, "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

POCAzSig = {
    "xLabel": "z_{0} / #sigma(z_{0})", "xUnits": "", "xMin": -900.0, "xMax": +900.0, "binWidthX": 50.0 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.0f"      , "yUnits": "", "yMin": +5E-03, "yMax": 1e-01 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True ,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

ChiSq  = {
    "xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +150.0, "binWidthX": 2.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": ""    , "yMin": 1E-05 , "yMax": +1e00, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

RedChiSq = {
    "xLabel": "#chi^{2}_{N}"   , "xUnits": "", "xMin": +0.00, "xMax": +10.0, "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-04, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

Charge = {
    "xLabel": "Charge (e)"     , "xUnits": "", "xMin": -5.0 , "xMax": 5.0 , "binWidthX": 1.0 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-5 , "yMax": 2e0 , "binWidthY": None, "yCutLines": [ ], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

D0 = {
    "xLabel": "d_{0}"          , "xUnits": "cm", "xMin": -1.0 , "xMax": +1.0, "binWidthX": 0.02, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1E-05, "yMax": +2.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

D0Abs = {
    "xLabel": "|d_{0}|"         , "xUnits": "cm", "xMin": 0.0 , "xMax": +0.4, "binWidthX": 0.01, "xCutLines": [4.0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f", "yUnits": ""  , "yMin": 1e-04, "yMax": +2.0, "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

D0Sig = {
    "xLabel": "d_{0}/#sigma(d_{0})", "xUnits": "", "xMin": -20.0 , "xMax": +20.0, "binWidthX": 0.50, "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f"    , "yUnits": "", "yMin": 1e-03, "yMax": 1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

D0SigAbs = {
    "xLabel": "|d_{0}/#sigma(d_{0})|", "xUnits": "", "xMin": 0.0  , "xMax": +20.0, "binWidthX": 0.50, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f"      , "yUnits": "", "yMin": 1e-03, "yMax": 1.0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

SigmaRInv = {
    "xLabel": "#sigma_{#rho^{-1}}", "xUnits": "", "xMin": 0.0  , "xMax": 5e-05, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f"   , "yUnits": "", "yMin": 1E-05, "yMax": 2e+00, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

SigmaPhi0 = {
    "xLabel": "#sigma_{#phi_{0}}", "xUnits": "", "xMin": +0.0001, "xMax": +0.0005, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f"  , "yUnits": "", "yMin": 1E-05  , "yMax": 1.0    , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

SigmaD0 = {
    "xLabel": "#sigma_{d_{0}}" , "xUnits": "cm", "xMin": +0.0015, "xMax": 0.005, "binWidthX": 1E-04 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f", "yUnits": ""  , "yMin": 1E-04  , "yMax": 2E-00, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

SigmaT  = {
    "xLabel": "#sigma_{tan(#lambda)}", "xUnits": "", "xMin": 0.0  , "xMax": 2e-3, "binWidthX": 1e-04, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.4f"      , "yUnits": "", "yMin": 1E-05, "yMax": 1.0   , "binWidthY": None , "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

SigmaZ0  = {
    "xLabel": "#sigma_{z_{0}}" , "xUnits": "", "xMin": 0.0 , "xMax": 1.8e-2, "binWidthX": 5e-4, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.5f", "yUnits": "", "yMin": 1E-05 , "yMax": 1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

StubPtCons = {
    "xLabel": "p_{T}^{stub} constistency", "xUnits": "", "xMin": 0.0  , "xMax": 50.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-3 , "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

NPix = {
    "xLabel": "Pixel Hits", "xUnits": "", "xMin": 1.0, "xMax": 6.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": 1e-2, "yMax":  1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.66, "xLegMax": 0.85, "yLegMin": 0.72, "yLegMax": 0.82 }

NTkTaus = {
    "xLabel": "Multiplicity"    , "xUnits": "", "xMin": +0.0, "xMax": 10.0, "binWidthX": 1.0 , "xCutLines": [2], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f" , "yUnits": "", "yMin": 1e-4, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": True, 
    "ratioLabel": ratioLabel, "ratio": True, "invRatio": False, "yMinRatio": 1e-1, "yMaxRatio": 200, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

NStubs = {
    "xLabel": "Stubs Multiplicity", "xUnits": "", "xMin": +0.0, "xMax": 12.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": 1e-4, "yMax": +2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

NPsStubs = {
    "xLabel": "PS Stubs Multiplicity", "xUnits": "", "xMin": +0.0, "xMax": 12.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"      , "yUnits": "", "yMin": 1e-4, "yMax": +2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

NBStubs = {
    "xLabel": "Barrel Stubs Multiplicity", "xUnits": "", "xMin": +0.0, "xMax": 12.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-3, "yMax":  2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

NEStubs = {
    "xLabel": "Endcap Stubs Multiplicity", "xUnits": "", "xMin": +0.0, "xMax": 12.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": 1e-3, "yMax":  2.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82 }

IsGenuine = {
    "xLabel": "is Genuine"     , "xUnits": "", "xMin": -0.5, "xMax": 1.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-2, "yMax": 1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

IsUnknown = {
    "xLabel": "is Unknown"     , "xUnits": "", "xMin": -0.5, "xMax": 1.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-2, "yMax": 1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

IsCombinatoric = {
    "xLabel": "is Combinatoric", "xUnits": "", "xMin": -0.5, "xMax": 1.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": 1e-2, "yMax": 1.2, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

Rtau = {       
    "xLabel": "R_{#tau} = p_{T}^{ldg}/E_{T}^{calo}", "xUnits": "", "xMin": 0.0 , "xMax": 2.0, "binWidthX": 0.05, "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.2f"                    , "yUnits": "", "yMin": 1e-4, "yMax": 1.0, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": True , "logYRatio": True ,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": 1e-1, "yMaxRatio": 100, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

CHF = {       
    "xLabel": "CHF = #Sigma_{i}^{sig} p_{T}^{i}/E_{T}^{calo}", "xUnits": "", "xMin": 0.0 , "xMax": 10.0, "binWidthX": 0.10, "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f"                        , "yUnits": "", "yMin": 1e-4, "yMax": 1.0, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": True , "logYRatio": True , 
    "ratioLabel": ratioLabel, "ratio": True, "invRatio": False, "yMinRatio": 1e-1, "yMaxRatio": 50, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

NHF = {       
    "xLabel": "NHF = (E_{T}^{calo}-#Sigma_{i}^{sig} p_{T}^{i})/E_{T}^{calo}", "xUnits": "", "xMin": -2.0 , "xMax": +1.0, "binWidthX": 0.05, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f"                        , "yUnits": "", "yMin": 1e-4, "yMax": 1.0, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": True , "logYRatio": True , 
    "ratioLabel": ratioLabel, "ratio": True, "invRatio": False, "yMinRatio": 1e-1, "yMaxRatio": 50, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

NHFAbs = {       
    "xLabel": "|NHF| = |E_{T}^{calo}-#Sigma_{i}^{sig} p_{T}^{i}|/E_{T}^{calo}", "xUnits": "", "xMin": 0.0 , "xMax": +3.0, "binWidthX": 0.05, "xCutLines": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f"                                               , "yUnits": "", "yMin": 1e-5, "yMax": +1e0, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": True , 
    "logYRatio": True , "ratioLabel": ratioLabel, "ratio": True, "invRatio": False, "yMinRatio": 1e-1, "yMaxRatio": 50, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

NSigTks = {
    "xLabel": "Signal Tks Multiplicity"   , "xUnits": "", "xMin": -0.5 , "xMax": 10.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +1e-4, "yMax":  1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

NIsoTks = {
    "xLabel": "Iso Tks Multiplicity"   , "xUnits": "", "xMin": -0.5 , "xMax": 10.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": +1e-4, "yMax":  1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

InvMass = {  
    #"xLabel": " m_{tks}"       , "xUnits": "GeV/c^{2}", "xMin": +0.0 , "xMax": 2.5, "binWidthX": 0.05 , "xCutLines": [1.77], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "xLabel": " m_{tks}"       , "xUnits": "GeV/c^{2}", "xMin": +0.0 , "xMax": 5.0, "binWidthX": 0.05 , "xCutLines": [1.77], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.2f", "yUnits": ""         , "yMin": +1e-4, "yMax": 1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

DeltaR = {
    "xLabel": "#Delta R"       , "xUnits": "", "xMin": +0.0 , "xMax": 0.35, "binWidthX": 0.005, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.3f", "yUnits": "", "yMin": +1e-4, "yMax": 1.00, "binWidthY": None , "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

RelIso = {      
    # "xLabel": "RelIso = #Sigma_{i}^{iso} p_{T}^{i}/p_{T}^{m}", "xUnits": "", "xMin": +0.0 , "xMax": +1.0, "binWidthX": None, "gridX": True, "logX": False, "logXRatio": False,
    "xLabel": "RelIso = #Sigma_{i}^{iso} p_{T}^{i}/p_{T}^{m}", "xUnits": "", "xMin": +0.0 , "xMax": +1.5, "binWidthX": None, "gridX": True, "logX": False, "logXRatio": False,
    "yLabel": "Entries / %0.2f"                              , "yUnits": "", "yMin": +1e-3, "yMax": 1e-1, "binWidthY": None, "gridY": True, "logY": True , "logYRatio": True,
    "ratioLabel": ratioLabel, "ratio": True, "invRatio": False, "yMinRatio": 1e-01, "yMaxRatio": 500, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

VtxIso = {      
    "xLabel": "VtxIso = min(|z_{0}^{tk_{m}} - z_{0}^{tk_{i}}|)", "xUnits": "cm", "xMin": -10.0, "xMax": +10.0, "binWidthX": 0.2 , "xCutLines": [-1.0, +1.0], "gridX": True, "logX": False,
    "yLabel": "Entries / %0.2f"                                , "yUnits": ""  , "yMin": +1e-3, "yMax": +1e0 , "binWidthY": None, "yCutLines": []          , "gridY": True, "logY": True ,
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}

VtxIsoAbs = {      
    "xLabel": "VtxIso = min(|z_{0}^{tk_{m}} - z_{0}^{tk_{i}}|)", "xUnits": "cm", "xMin": +0.0 , "xMax": 2.0, "binWidthX": None, "xCutLines": [+0.4], "xCutBoxes": [], "gridX": True, "logX": False, 
    "yLabel": "Entries / %0.2f"                                , "yUnits": ""  , "yMin": +1e-3, "yMax": 1e0, "binWidthY": None, "yCutLines": []    , "yCutBoxes": [], "gridY": True, "logY": True , 
    "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.65, "xLegMax": 0.90, "yLegMin": 0.68, "yLegMax": 0.82}


###############################################################
### L1TkTaus
###############################################################
folder = ""
hL1TkTau_Multiplicity= m_histos.TH1orTH2( folder, "L1TkTau_Multiplicity", "L1TkTau", None , **NTkTaus   )
hL1TkTau_Rtau        = m_histos.TH1orTH2( folder, "L1TkTau_Rtau"        , "L1TkTau", None , **Rtau      )
hL1TkTau_CHF         = m_histos.TH1orTH2( folder, "L1TkTau_CHF"         , "L1TkTau", None , **CHF       )
hL1TkTau_NHF         = m_histos.TH1orTH2( folder, "L1TkTau_NHF"         , "L1TkTau", None , **NHF       )
hL1TkTau_NHFAbs      = m_histos.TH1orTH2( folder, "L1TkTau_NHFAbs"      , "L1TkTau", None , **NHFAbs    )
hL1TkTau_NSigTks     = m_histos.TH1orTH2( folder, "L1TkTau_NSigTks"     , "L1TkTau", None , **NSigTks   )
hL1TkTau_NIsoTks     = m_histos.TH1orTH2( folder, "L1TkTau_NIsoTks"     , "L1TkTau", None , **NIsoTks   )
hL1TkTau_InvMass     = m_histos.TH1orTH2( folder, "L1TkTau_InvMass"     , "L1TkTau", None , **InvMass   )
hL1TkTau_InvMassIncl = m_histos.TH1orTH2( folder, "L1TkTau_InvMassIncl" , "L1TkTau", None , **InvMass   )
hL1TkTau_SigConeRMin = m_histos.TH1orTH2( folder, "L1TkTau_SigConeRMin" , "L1TkTau", None , **DeltaR    )
hL1TkTau_SigConeRMax = m_histos.TH1orTH2( folder, "L1TkTau_SigConeRMax" , "L1TkTau", None , **DeltaR    )
hL1TkTau_IsoConeRMin = m_histos.TH1orTH2( folder, "L1TkTau_IsoConeRMin" , "L1TkTau", None , **DeltaR    )
hL1TkTau_IsoConeRMax = m_histos.TH1orTH2( folder, "L1TkTau_IsoConeRMax" , "L1TkTau", None , **DeltaR    )
hL1TkTau_Charge      = m_histos.TH1orTH2( folder, "L1TkTau_Charge"      , "L1TkTau", None , **Charge    )
hL1TkTau_RelIso      = m_histos.TH1orTH2( folder, "L1TkTau_RelIso"      , "L1TkTau", None , **RelIso    )
hL1TkTau_VtxIso      = m_histos.TH1orTH2( folder, "L1TkTau_VtxIso"      , "L1TkTau", None , **VtxIso    )
hL1TkTau_VtxIsoAbs   = m_histos.TH1orTH2( folder, "L1TkTau_VtxIsoAbs"   , "L1TkTau", None , **VtxIsoAbs )
hL1TkTau_DeltaRGenP  = m_histos.TH1orTH2( folder, "L1TkTau_DeltaRGenP"  , "L1TkTau", None , **DeltaR    )

###############################################################
### L1TkTaus; SigTks 
###############################################################
folder = ""
hL1TkTau_SigTks_Pt            = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_Pt"           , "SigTks", None, **Pt        )
hL1TkTau_SigTks_Eta           = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_Eta"          , "SigTks", None, **Eta       )
hL1TkTau_SigTks_POCAz         = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_POCAz"        , "SigTks", None, **POCAz     )
hL1TkTau_SigTks_DeltaPOCAz    = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_DeltaPOCAz"   , "SigTks", None, **dPOCAz    )
hL1TkTau_SigTks_d0            = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_d0"           , "SigTks", None, **D0        )
hL1TkTau_SigTks_d0Abs         = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_d0Abs"        , "SigTks", None, **D0Abs     )
hL1TkTau_SigTks_d0Sig         = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_d0Sig"        , "SigTks", None, **D0Sig     )  
hL1TkTau_SigTks_d0SigAbs      = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_d0SigAbs"     , "SigTks", None, **D0SigAbs  )  
hL1TkTau_SigTks_PtRel         = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_PtRel"        , "SigTks", None, **PtRel     )
hL1TkTau_SigTks_StubPtCons    = m_histos.TH1orTH2( folder, "L1TkTau_SigTks_StubPtCons"   , "SigTks", None, **StubPtCons)

###############################################################
### L1TkTaus; IsoTks 
###############################################################
folder = ""
hL1TkTau_IsoTks_Pt            = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_Pt"           , "IsoTks", None, **Pt        )
hL1TkTau_IsoTks_Eta           = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_Eta"          , "IsoTks", None, **Eta       )
hL1TkTau_IsoTks_POCAz         = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_POCAz"        , "IsoTks", None, **POCAz     )
hL1TkTau_IsoTks_DeltaPOCAz    = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_DeltaPOCAz"   , "IsoTks", None, **dPOCAz    )
hL1TkTau_IsoTks_d0            = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_d0"           , "IsoTks", None, **D0        )
hL1TkTau_IsoTks_d0Abs         = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_d0Abs"        , "IsoTks", None, **D0Abs     )
hL1TkTau_IsoTks_d0Sig         = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_d0Sig"        , "IsoTks", None, **D0Sig     )  
hL1TkTau_IsoTks_d0SigAbs      = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_d0SigAbs"     , "IsoTks", None, **D0SigAbs  )  
hL1TkTau_IsoTks_PtRel         = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_PtRel"        , "IsoTks", None, **PtRel     )
hL1TkTau_IsoTks_StubPtCons    = m_histos.TH1orTH2( folder, "L1TkTau_IsoTks_StubPtCons"   , "IsoTks", None, **StubPtCons)

###############################################################
### L1TkTaus: Matching Pixel Track
###############################################################
folder = ""
hL1TkTau_MatchPixTk_DeltaR        = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_DeltaR"       , "Match PixTk", None, **DeltaR    )
hL1TkTau_MatchPixTk_PtRel         = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_PtRel"        , "Match PixTk", None, **PtRel     )
hL1TkTau_MatchPixTk_NPixHits      = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_NPixHits"     , "Match PixTk", None, **NPix      )
hL1TkTau_MatchPixTk_NStubs        = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_NStubs"       , "Match PixTk", None, **NStubs    )
hL1TkTau_MatchPixTk_NPsStubs      = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_NPsStubs"     , "Match PixTk", None, **NPsStubs  )
hL1TkTau_MatchPixTk_NBarrelStubs  = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_NBarrelStubs" , "Match PixTk", None, **NBStubs   )
hL1TkTau_MatchPixTk_NEndcapStubs  = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_NEndcapStubs" , "Match PixTk", None, **NEStubs   )
hL1TkTau_MatchPixTk_StubPtCons    = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_StubPtCons"   , "Match PixTk", None, **StubPtCons)
hL1TkTau_MatchPixTk_Pt            = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_Pt"           , "Match PixTk", None, **Pt        )
hL1TkTau_MatchPixTk_Eta           = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_Eta"          , "Match PixTk", None, **Eta       )
hL1TkTau_MatchPixTk_POCAz         = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_POCAz"        , "Match PixTk", None, **POCAz     )
hL1TkTau_MatchPixTk_POCAzSig      = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_POCAzSig"     , "Match PixTk", None, **POCAzSig  )
hL1TkTau_MatchPixTk_d0            = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_d0"           , "Match PixTk", None, **D0        )
hL1TkTau_MatchPixTk_d0Abs         = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_d0Abs"        , "Match PixTk", None, **D0Abs     )
hL1TkTau_MatchPixTk_d0Sig         = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_d0Sig"        , "Match PixTk", None, **D0Sig     )  
hL1TkTau_MatchPixTk_d0SigAbs      = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_d0SigAbs"     , "Match PixTk", None, **D0SigAbs  )  
hL1TkTau_MatchPixTk_ChiSquared    = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_ChiSquared"   , "Match PixTk", None, **ChiSq     )
hL1TkTau_MatchPixTk_RedChiSquared = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_RedChiSquared", "Match PixTk", None, **RedChiSq  )
hL1TkTau_MatchPixTk_SigmaRInv     = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_SigmaRInv"    , "Match PixTk", None, **SigmaRInv )
hL1TkTau_MatchPixTk_SigmaPhi0     = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_SigmaPhi0"    , "Match PixTk", None, **SigmaPhi0 )
hL1TkTau_MatchPixTk_SigmaD0       = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_SigmaD0"      , "Match PixTk", None, **SigmaD0   )
hL1TkTau_MatchPixTk_SigmaT        = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_SigmaT"       , "Match PixTk", None, **SigmaT    )
hL1TkTau_MatchPixTk_SigmaZ0       = m_histos.TH1orTH2( folder, "L1TkTau_MatchPixTk_SigmaZ0"      , "Match PixTk", None, **SigmaZ0   )

###############################################################
### L1TkTaus: Matching Track
###############################################################
hL1TkTau_MatchTk_DeltaR         = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_DeltaR"        , "Match Tk", None, **DeltaR         )
hL1TkTau_MatchTk_PtRel          = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_PtRel"         , "Match Tk", None, **PtRel          )
hL1TkTau_MatchTk_Pt             = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_Pt"            , "Match Tk", None, **Pt             )
hL1TkTau_MatchTk_Eta            = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_Eta"           , "Match Tk", None, **Eta            )
hL1TkTau_MatchTk_POCAz          = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_POCAz"         , "Match Tk", None, **POCAz          )
hL1TkTau_MatchTk_d0             = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_d0"            , "Match Tk", None, **D0             )
hL1TkTau_MatchTk_d0Abs          = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_d0"            , "Match Tk", None, **D0Abs          )
hL1TkTau_MatchTk_ChiSquared     = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_ChiSquared"    , "Match Tk", None, **ChiSq          )
hL1TkTau_MatchTk_RedChiSquared  = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_RedChiSquared" , "Match Tk", None, **RedChiSq       )
hL1TkTau_MatchTk_NStubs         = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_NStubs"        , "Match Tk", None, **NStubs         )
hL1TkTau_MatchTk_NPsStubs       = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_NPsStubs"      , "Match Tk", None, **NPsStubs       )
hL1TkTau_MatchTk_NBarrelStubs   = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_NBarrelStubs"  , "Match Tk", None, **NBStubs        )
hL1TkTau_MatchTk_NEndcapStubs   = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_NEndcapStubs"  , "Match Tk", None, **NEStubs        )
hL1TkTau_MatchTk_StubPtCons     = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_StubPtCons"    , "Match Tk", None, **StubPtCons     )
hL1TkTau_MatchTk_IsGenuine      = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_IsGenuine"     , "Match Tk", None, **IsGenuine      )
hL1TkTau_MatchTk_IsUnknown      = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_IsUnknown"     , "Match Tk", None, **IsUnknown      )
hL1TkTau_MatchTk_IsCombinatoric = m_histos.TH1orTH2( folder, "L1TkTau_MatchTk_IsCombinatoric", "Match Tk", None, **IsCombinatoric )

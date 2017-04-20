###############################################################
### Author .........: Alexandros X. Attikis 
### Institute ......: University Of Cyprus (UCY)
### Email ..........: attikis@cern.ch
###############################################################

###############################################################
### All imported modules
###############################################################
# System
import os
import sys
import numpy
import math

# User
import HLTausAnalysis.NtupleAnalysis.tools.old.plotter as m_plotter
import HLTausAnalysis.NtupleAnalysis.tools.old.histos as m_histos
import HLTausAnalysis.NtupleAnalysis.tools.old.styles as m_styles
import HLTausAnalysis.NtupleAnalysis.tools.old.aux as m_aux

# ROOT
import ROOT

###############################################################
### Histogram Options
###############################################################
hFolder     = ""
yMin        = 1
yMax        = None
yMinRatio   = 0.0
yMaxRatio   = 4.9
bRatio      = False
bInvRatio   = True
ratioLabel  = "1/Ratio"
normFactor  = None

###############################################################
### Histogram Attributes
###############################################################
kwargs   = {}

Rate = {
    #"xLabel": "E_{T}", "xUnits": "GeV", "xMin": 0.0 , "xMax": +100.0, "binWidthX": None, "xCutLines": []  , "xCutBoxes": [[0.0, 15.0, ROOT.kBlack]], "gridX": True, "logX": False,
    "xLabel": "E_{T}", "xUnits": "GeV", "xMin": 0.0 , "xMax": +120.0, "binWidthX": None, "xCutLines": []  , "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Rate (kHz) / %0.0f", "yUnits": "", "yMin": +1E0, "yMax": 5E+04, "yCutLines": [50], "yCutBoxes": [[50.0, 5E+04, ROOT.kBlack]], "gridY": True, "logY": True ,
    "logYRatio": True, "ratioLabel": ratioLabel, "ratio": False, "invRatio": True, "yMinRatio": 1e-1, "yMaxRatio": 50.0, "normaliseTo": normFactor, "drawOptions": "PF",
    "legOptions": "LP", "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.70, "yLegMax": 0.82 }

DiRate = {
    "xLabel": "E_{T}", "xUnits": "GeV", "xMin": 0.0 , "xMax": +100.0, "binWidthX": None, "xCutLines": []  , "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Rate (kHz) / %0.0f", "yUnits": ""   , "yMin": +1E0, "yMax": 5E+04 , "yCutLines": [50], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": True, 
    "ratioLabel": ratioLabel, "ratio": False, "invRatio": True, "yMinRatio": +1e0, "yMaxRatio": 200.0, "normaliseTo": normFactor, "drawOptions": "PF", "legOptions": "LP",
    "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.70, "yLegMax": 0.82 }

Rate2D = {
    "xLabel": "E_{T}^{Ldg} / %0.0f"   , "xUnits": "GeV", "xMin": 0.0 , "xMax": +100.0, "binWidthX": 1, "xCutLines": [], "xCutBoxes": [], "gridX": False, "logX": False,
    "yLabel": "E_{T}^{subLdg} / %0.0f", "yUnits": "GeV", "yMin": 0.0 , "yMax": +100.0, "binWidthY": 1, "yCutLines": [], "yCutBoxes": [], "gridY": False, "logY": False ,
    "zLabel": "Rate", "zUnits": "kHz", "zMin": 1E+0, "zMax": 5E+04 , "zCutLinesErrors": True, "zCutLines": [50.0], "zCutBoxes": [], "gridZ": False, "logZ": True ,
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.70, "yLegMax": 0.82 }

Eff = {
    "xLabel": "E_{T}"             , "xUnits": "GeV", "xMin": 0.0 , "xMax": +100.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.0f", "yUnits": ""   , "yMin": 0.0 , "yMax": +1.09 , "binWidthY": None, "yCutLines": [], "yCutBoxes":  [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": True, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 1.1, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.70, "yLegMax": 0.82 }

Eff2D = {
    "xLabel": "E_{T}^{Ldg} / %0.0f"   , "xUnits": "GeV", "xMin": 0.0, "xMax": +100.0, "binWidthX": 1, "xCutLines": [], "xCutBoxes": [], "gridX": False, "logX": False,
    "yLabel": "E_{T}^{subLdg} / %0.0f", "yUnits": "GeV", "yMin": 0.0, "yMax": +100.0, "binWidthY": 1, "yCutLines": [], "yCutBoxes": [], "gridY": False, "logY": False ,
    "zLabel": "Efficiency"            , "zUnits": ""   , "zMin": 0.0, "zMax": 1.0   , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.70, "yLegMax": 0.82 }

ROC = {
    "xLabel": "Efficiency", "xUnits": ""   , "xMin": 0.00, "xMax": +0.60 , "binWidthX": None, "xCutLines": [0.51] , "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Rate" , "yUnits": "kHz", "yMin": 1E+0, "yMax": +0.5E+03, "binWidthY": None, "yCutLines": [], "yCutBoxes": [[50.0, 50.0, ROOT.kBlue]], "gridY": True, "logY": True,
    "legOptions": "FL", "drawOptions": "ACE3", "xLegMin": 0.20, "xLegMax": 0.40, "yLegMin": 0.74, "yLegMax": 0.92}

ROC_C = {
    "xLabel": "Efficiency", "xUnits": ""   , "xMin": 0.00, "xMax": +0.45 , "binWidthX": None, "xCutLines": [] , "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Rate" , "yUnits": "kHz", "yMin": 1E+0, "yMax": +0.5E+03, "binWidthY": None, "yCutLines": [], "yCutBoxes": [[50.0, 50.0, ROOT.kBlue]], "gridY": True, "logY": True,
    "legOptions": "FL", "drawOptions": "ACE3", "xLegMin": 0.20, "xLegMax": 0.40, "yLegMin": 0.74, "yLegMax": 0.92}

ROC_I = {
    "xLabel": "Efficiency", "xUnits": ""   , "xMin": 0.00, "xMax": +0.25 , "binWidthX": None, "xCutLines": [] , "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Rate" , "yUnits": "kHz", "yMin": 1E+0, "yMax": +0.5E+03, "binWidthY": None, "yCutLines": [], "yCutBoxes": [[50.0, 50.0, ROOT.kBlue]], "gridY": True, "logY": True,
    "legOptions": "FL", "drawOptions": "ACE3", "xLegMin": 0.20, "xLegMax": 0.40, "yLegMin": 0.74, "yLegMax": 0.92}

ROC_F = {
    "xLabel": "Efficiency", "xUnits": ""   , "xMin": 0.00, "xMax": +0.14 , "binWidthX": None, "xCutLines": [] , "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Rate" , "yUnits": "kHz", "yMin": 1E+0, "yMax": +0.5E+03, "binWidthY": None, "yCutLines": [], "yCutBoxes": [[50.0, 50.0, ROOT.kBlue]], "gridY": True, "logY": True,
    "legOptions": "FL", "drawOptions": "ACE3", "xLegMin": 0.20, "xLegMax": 0.40, "yLegMin": 0.74, "yLegMax": 0.92}

TurnOn = {
    #"xLabel": "#tau_{MC} E_{T}^{vis}", "xUnits": "GeV", "xMin": 0.0 , "xMax": 150.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    #"yLabel": "Signal Efficiency / %0.0f"   , "yUnits": ""   , "yMin": 0.0 , "yMax": +1.17 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": False,
    "xLabel": "E_{T}^{vis}", "xUnits": "GeV", "xMin": 0.0 , "xMax": 150.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.0f"   , "yUnits": ""   , "yMin": 0.0 , "yMax": +1.17 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": 1.1, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.26, "yLegMax": 0.40 }

TurnOn50 = {
    #"xLabel": "#tau_{MC} E_{T}^{vis}", "xUnits": "GeV", "xMin": 0.0 , "xMax": 150.0, "binWidthX": None, "xCutLines": [50], "xCutBoxes": [], "gridX": True, "logX": False,
    #"yLabel": "Signal Efficiency / %0.0f"   , "yUnits": ""   , "yMin": 0.0 , "yMax": +1.17 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": False,
    "xLabel": "E_{T}^{vis}", "xUnits": "GeV", "xMin": 0.0 , "xMax": 150.0, "binWidthX": None, "xCutLines": [50], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.0f"   , "yUnits": ""   , "yMin": 0.0 , "yMax": +1.17 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": 1.1, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.26, "yLegMax": 0.40 }

TurnOn25 = {
    #"xLabel": "#tau_{MC} E_{T}^{vis}", "xUnits": "GeV", "xMin": 0.0 , "xMax": +150.0, "binWidthX": None, "xCutLines": [25], "xCutBoxes": [], "gridX": True, "logX": False,
    #"yLabel": "Signal Efficiency / %0.0f"   , "yUnits": ""   , "yMin": 0.0 , "yMax": +1.17 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": False,
    "xLabel": "E_{T}^{vis}", "xUnits": "GeV", "xMin": 0.0 , "xMax": +150.0, "binWidthX": None, "xCutLines": [25], "xCutBoxes": [], "gridX": True, "logX": False,
    "yLabel": "Efficiency / %0.0f"   , "yUnits": ""   , "yMin": 0.0 , "yMax": +1.17 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": False,
    "ratioLabel": "Ratio", "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": 1.1, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.26, "yLegMax": 0.40 }


VertexZ = {
    "xLabel": "z_{vtx}"         , "xUnits": "cm", "xMin": -20.0, "xMax": +20.0, "binWidthX": 0.5 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
    "yLabel": "Entries / %0.1f" , "yUnits": ""  , "yMin": 1e-3 , "yMax": None , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": False, "logYRatio": False, 
    "ratioLabel": ratioLabel, "ratio": True, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
    "xLegMin": 0.68, "xLegMax": 0.85, "yLegMin": 0.70, "yLegMax": 0.82 }

VertexXY = {
    "xLabel": "x_{vtx}  / %0.3f", "xUnits": "cm", "xMin": -0.0049, "xMax": +0.0049, "binWidthX": 0.0005, "gridX": False, "logX": False, 
    "yLabel": "y_{vtx}  / %0.3f", "yUnits": "cm", "yMin": -0.0049, "yMax": +0.0049, "binWidthY": 0.0005, "gridY": False, "logX": False, 
    "zLabel": "Entries"         , "zUnits": ""   , "zMin": 0.0, "zMax": None  , "zCutLinesErrors": False, "zCutLines": [], "zCutBoxes": [], "gridZ": False, "logZ": False,
    "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.70, "yLegMax": 0.82 }



###############################################################
### Miscellaneous
###############################################################
hHepMCEvt_VtxZ        = m_histos.TH1orTH2( hFolder, "HepMCEvt_VtxZ"     , "Truth", None, **VertexZ  )
hHepMCEvt_VtxX_VtxY   = m_histos.TH1orTH2( hFolder, "HepMCEvt_VtxX_VtxY", "Truth", None, **VertexXY  )

###############################################################
### SingleTau
###############################################################
hCalo_Rate   = m_histos.TH1orTH2( hFolder, "Calo_Rate"  , "Calo", "SingleTau_Rate_Calo"  , **Rate )
hCalo_Rate_C = m_histos.TH1orTH2( hFolder, "Calo_Rate_C", "Calo", "SingleTau_Rate_Calo_C", **Rate )
hCalo_Rate_I = m_histos.TH1orTH2( hFolder, "Calo_Rate_I", "Calo", "SingleTau_Rate_Calo_I", **Rate )
hCalo_Rate_F = m_histos.TH1orTH2( hFolder, "Calo_Rate_F", "Calo", "SingleTau_Rate_Calo_F", **Rate )

hTk_Rate     = m_histos.TH1orTH2( hFolder, "Tk_Rate"    , "+Tk" , "SingleTau_Rate_Tk"    , **Rate )
hTk_Rate_C   = m_histos.TH1orTH2( hFolder, "Tk_Rate_C"  , "+Tk" , "SingleTau_Rate_Tk_C"  , **Rate )
hTk_Rate_I   = m_histos.TH1orTH2( hFolder, "Tk_Rate_I"  , "+Tk" , "SingleTau_Rate_Tk_I"  , **Rate )
hTk_Rate_F   = m_histos.TH1orTH2( hFolder, "Tk_Rate_F"  , "+Tk" , "SingleTau_Rate_Tk_F"  , **Rate )

hVtxIso_Rate   = m_histos.TH1orTH2( hFolder, "VtxIso_Rate"  , "CaloTk", "SingleTau_Rate_VtxIso", **Rate )
hVtxIso_Rate_C = m_histos.TH1orTH2( hFolder, "VtxIso_Rate_C", "CaloTk", "SingleTau_Rate_VtxIso", **Rate )
hVtxIso_Rate_I = m_histos.TH1orTH2( hFolder, "VtxIso_Rate_I", "CaloTk", "SingleTau_Rate_VtxIso", **Rate )
hVtxIso_Rate_F = m_histos.TH1orTH2( hFolder, "VtxIso_Rate_F", "CaloTk", "SingleTau_Rate_VtxIso", **Rate )

hCalo_Eff   = m_histos.TH1orTH2( hFolder, "Calo_Eff"  , "Calo", "SingleTau_Eff_Calo"  , **Eff )
hCalo_Eff_C = m_histos.TH1orTH2( hFolder, "Calo_Eff_C", "Calo", "SingleTau_Eff_Calo_C", **Eff )
hCalo_Eff_I = m_histos.TH1orTH2( hFolder, "Calo_Eff_I", "Calo", "SingleTau_Eff_Calo_I", **Eff )
hCalo_Eff_F = m_histos.TH1orTH2( hFolder, "Calo_Eff_F", "Calo", "SingleTau_Eff_Calo_F", **Eff )

hTk_Eff     = m_histos.TH1orTH2( hFolder, "Tk_Eff"  , "+Tk", "SingleTau_Eff_Tk"  , **Eff )
hTk_Eff_C   = m_histos.TH1orTH2( hFolder, "Tk_Eff_C", "+Tk", "SingleTau_Eff_Tk_C", **Eff )
hTk_Eff_I   = m_histos.TH1orTH2( hFolder, "Tk_Eff_I", "+Tk", "SingleTau_Eff_Tk_I", **Eff )
hTk_Eff_F   = m_histos.TH1orTH2( hFolder, "Tk_Eff_F", "+Tk", "SingleTau_Eff_Tk_F", **Eff )

hVtxIso_Eff   = m_histos.TH1orTH2( hFolder, "VtxIso_Eff"  , "CaloTk" , "SingleTau_Eff_VtxIso"  , **Eff )
hVtxIso_Eff_C = m_histos.TH1orTH2( hFolder, "VtxIso_Eff_C", "CaloTk" , "SingleTau_Eff_VtxIso_C", **Eff )
hVtxIso_Eff_I = m_histos.TH1orTH2( hFolder, "VtxIso_Eff_I", "CaloTk" , "SingleTau_Eff_VtxIso_I", **Eff )
hVtxIso_Eff_F = m_histos.TH1orTH2( hFolder, "VtxIso_Eff_F", "CaloTk" , "SingleTau_Eff_VtxIso_F", **Eff )

hCalo_TurnOn50    = m_histos.TH1orTH2( hFolder, "Calo_TurnOn50"  , "Calo"   , "TurnOn_SingleTau_50GeV"  , **TurnOn50 )
hTk_TurnOn50      = m_histos.TH1orTH2( hFolder, "Tk_TurnOn50"    , "+Tk"    , "TurnOn_SingleTau_50GeV"    , **TurnOn50 )
#hVtxIso_TurnOn50  = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn50", "+VtxIso", "TurnOn_SingleTau_50GeV", **TurnOn50 )
hVtxIso_TurnOn50  = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn50", "Calo+Tks", "TurnOn_SingleTau_50GeV", **TurnOn50 )

hCalo_TurnOn_SingleTau50KHz   = m_histos.TH1orTH2( hFolder, "Calo_TurnOn_SingleTau50KHz"  , "Calo"   , "TurnOn_SingleTau_50kHz", **TurnOn )
hTk_TurnOn_SingleTau50KHz     = m_histos.TH1orTH2( hFolder, "Tk_TurnOn_SingleTau50KHz"    , "+Tk"    , "TurnOn_SingleTau_50kHz", **TurnOn )
#hVtxIso_TurnOn_SingleTau50KHz = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn_SingleTau50KHz", "+VtxIso", "TurnOn_SingleTau_50kHz", **TurnOn )
hVtxIso_TurnOn_SingleTau50KHz = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn_SingleTau50KHz", "Calo + Tks", "TurnOn_SingleTau_50kHz", **TurnOn )

###############################################################
### DiTau (Indistinghuishable)
###############################################################
hDiTau_Rate_Calo   = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Calo"  , "Calo", "DiTau_Rate_Calo"  , **DiRate )
hDiTau_Rate_Calo_C = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Calo_C", "Calo", "DiTau_Rate_Calo_C", **DiRate )
hDiTau_Rate_Calo_I = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Calo_I", "Calo", "DiTau_Rate_Calo_I", **DiRate )
hDiTau_Rate_Calo_F = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Calo_F", "Calo", "DiTau_Rate_Calo_F", **DiRate )

hDiTau_Rate_Tk   = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Tk"  , "+Tk", "DiTau_Rate_Tk"  , **DiRate )
hDiTau_Rate_Tk_C = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Tk_C", "+Tk", "DiTau_Rate_Tk_C", **DiRate )
hDiTau_Rate_Tk_I = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Tk_I", "+Tk", "DiTau_Rate_Tk_I", **DiRate )
hDiTau_Rate_Tk_F = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Tk_F", "+Tk", "DiTau_Rate_Tk_F", **DiRate )

hDiTau_Rate_VtxIso   = m_histos.TH1orTH2( hFolder, "DiTau_Rate_VtxIso"  , "CaloTk", "DiTau_Rate_VtxIso"  , **DiRate )
hDiTau_Rate_VtxIso_C = m_histos.TH1orTH2( hFolder, "DiTau_Rate_VtxIso_C", "CaloTk", "DiTau_Rate_VtxIso_C", **DiRate )
hDiTau_Rate_VtxIso_I = m_histos.TH1orTH2( hFolder, "DiTau_Rate_VtxIso_I", "CaloTk", "DiTau_Rate_VtxIso_I", **DiRate )
hDiTau_Rate_VtxIso_F = m_histos.TH1orTH2( hFolder, "DiTau_Rate_VtxIso_F", "CaloTk", "DiTau_Rate_VtxIso_F", **DiRate )

hDiTau_Eff_Calo   = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Calo"  , "Calo", "DiTau_Eff_Calo"  , **Eff )
hDiTau_Eff_Calo_C = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Calo_C", "Calo", "DiTau_Eff_Calo_C", **Eff )
hDiTau_Eff_Calo_I = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Calo_I", "Calo", "DiTau_Eff_Calo_I", **Eff )
hDiTau_Eff_Calo_F = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Calo_F", "Calo", "DiTau_Eff_Calo_F", **Eff )

hDiTau_Eff_Tk     = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Tk"  , "+Tk", "DiTau_Eff_Tk"  , **Eff )
hDiTau_Eff_Tk_C   = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Tk_C", "+Tk", "DiTau_Eff_Tk_C", **Eff )
hDiTau_Eff_Tk_I   = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Tk_C", "+Tk", "DiTau_Eff_Tk_I", **Eff )
hDiTau_Eff_Tk_F   = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Tk_C", "+Tk", "DiTau_Eff_Tk_F", **Eff )

hDiTau_Eff_VtxIso   = m_histos.TH1orTH2( hFolder, "DiTau_Eff_VtxIso"  , "CaloTk", "DiTau_Eff_VtxIso"  , **Eff )
hDiTau_Eff_VtxIso_C = m_histos.TH1orTH2( hFolder, "DiTau_Eff_VtxIso_C", "CaloTk", "DiTau_Eff_VtxIso_C", **Eff )
hDiTau_Eff_VtxIso_I = m_histos.TH1orTH2( hFolder, "DiTau_Eff_VtxIso_I", "CaloTk", "DiTau_Eff_VtxIso_I", **Eff )
hDiTau_Eff_VtxIso_F = m_histos.TH1orTH2( hFolder, "DiTau_Eff_VtxIso_F", "CaloTk", "DiTau_Eff_VtxIso_F", **Eff )

hCalo_TurnOn25    = m_histos.TH1orTH2( hFolder, "Calo_TurnOn25"  , "Calo"   , "TurnOn_SingleTau_25GeV", **TurnOn25 )
hTk_TurnOn25      = m_histos.TH1orTH2( hFolder, "Tk_TurnOn25"    , "+Tk"    , "TurnOn_SingleTau_25GeV", **TurnOn25 )
#hVtxIso_TurnOn25  = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn25", "+VtxIso", "TurnOn_SingleTau_25GeV", **TurnOn25 )
hVtxIso_TurnOn25  = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn25", "Calo + Tks", "TurnOn_SingleTau_25GeV", **TurnOn25 )

hCalo_TurnOn_DiTau50KHz   = m_histos.TH1orTH2( hFolder, "Calo_TurnOn_DiTau50KHz"  , "Calo"   , "TurnOn_DiTau_50KHz", **TurnOn )
hTk_TurnOn_DiTau50KHz     = m_histos.TH1orTH2( hFolder, "Tk_TurnOn_DiTau50KHz"    , "+Tk"    , "TurnOn_DiTau_50KHz", **TurnOn )
#hVtxIso_TurnOn_DiTau50KHz = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn_DiTau50KHz", "+VtxIso", "TurnOn_DiTau_50KHz", **TurnOn )
hVtxIso_TurnOn_DiTau50KHz = m_histos.TH1orTH2( hFolder, "VtxIso_TurnOn_DiTau50KHz", "Calo + Tks", "TurnOn_DiTau_50KHz", **TurnOn )

###############################################################
### DiTau (Distinghuishable): Calo-Iso
###############################################################
hDiTau_Rate_Calo_Tk     = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Calo_Tk"    , "Calo, Tk"     , "DiTau_Rate_Calo_Tk"    , **Rate2D )
hDiTau_Rate_Calo_VtxIso = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Calo_VtxIso", "Calo, +VtxIso", "DiTau_Rate_Calo_VtxIso", **Rate2D )

hDiTau_Eff_Calo_Tk     = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Calo_Tk"    , "Calo, +Tk"    , "DiTau_Eff_Calo_Tk"    , **Eff2D )
hDiTau_Eff_Calo_VtxIso = m_histos.TH1orTH2( hFolder, "DiTau_Eff_Calo_VtxIso", "Calo, +VtxIso", "DiTau_Eff_Calo_VtxIso", **Eff2D )

###############################################################
### DiTau (Distinghuishable): Calo-Iso
###############################################################
hDiTau_Rate_Tk_VtxIso = m_histos.TH1orTH2( hFolder, "DiTau_Rate_Tk_VtxIso", "+Tk, +VtxIso", "DiTau_Rate_Tk_VtxIso", **Rate2D )
hDiTau_Eff_Tk_VtxIso  = m_histos.TH1orTH2( hFolder,"DiTau_Eff_Tk_VtxIso"  , "+Tk, +VtxIso", "DiTau_Eff_Tk_VtxIso" , **Eff2D  )

############################################################### 
### SingleTau 
############################################################### 
SingleTau_Rate = []
SingleTau_Rate.append(hCalo_Rate)
#SingleTau_Rate.append(hTk_Rate)
SingleTau_Rate.append(hVtxIso_Rate)

SingleTau_Rate_C = []
SingleTau_Rate_C.append(hCalo_Rate_C)
#SingleTau_Rate_C.append(hTk_Rate_C)
SingleTau_Rate_C.append(hVtxIso_Rate_C)

SingleTau_Rate_I = []
SingleTau_Rate_I.append(hCalo_Rate_I)
#SingleTau_Rate_I.append(hTk_Rate_I)
SingleTau_Rate_I.append(hVtxIso_Rate_I)

SingleTau_Rate_F = []
SingleTau_Rate_F.append(hCalo_Rate_F)
#SingleTau_Rate_F.append(hTk_Rate_F)
SingleTau_Rate_F.append(hVtxIso_Rate_F)

SingleTau_Eff = []
SingleTau_Eff.append(hCalo_Eff)
#SingleTau_Eff.append(hTk_Eff)
SingleTau_Eff.append(hVtxIso_Eff)

SingleTau_Eff_C = []
SingleTau_Eff_C.append(hCalo_Eff_C)
#SingleTau_Eff_C.append(hTk_Eff_C)
SingleTau_Eff_C.append(hVtxIso_Eff_C)

SingleTau_Eff_I = []
SingleTau_Eff_I.append(hCalo_Eff_I)
#SingleTau_Eff_I.append(hTk_Eff_I)
SingleTau_Eff_I.append(hVtxIso_Eff_I)

SingleTau_Eff_F = []
SingleTau_Eff_F.append(hCalo_Eff_F)
#SingleTau_Eff_F.append(hTk_Eff_F)
SingleTau_Eff_F.append(hVtxIso_Eff_F)

SingleTau_TurnOn = []
SingleTau_TurnOn.append(hCalo_TurnOn50)
#SingleTau_TurnOn.append(hTk_TurnOn50)
SingleTau_TurnOn.append(hVtxIso_TurnOn50)

SingleTau_TurnOn_50KHz = []
SingleTau_TurnOn_50KHz.append(hCalo_TurnOn_SingleTau50KHz)
#SingleTau_TurnOn_50KHz.append(hTk_TurnOn_SingleTau50KHz)
SingleTau_TurnOn_50KHz.append(hVtxIso_TurnOn_SingleTau50KHz)

SingleTau_ROCs = {}
SingleTau_ROCs[hCalo_Rate]   = hCalo_Eff
#SingleTau_ROCs[hTk_Rate]     = hTk_Eff
SingleTau_ROCs[hVtxIso_Rate] = hVtxIso_Eff

SingleTau_ROCs_C = {}
SingleTau_ROCs_C[hCalo_Rate_C]   = hCalo_Eff_C
#SingleTau_ROCs_C[hTk_Rate_C]     = hTk_Eff_C
SingleTau_ROCs_C[hVtxIso_Rate_C] = hVtxIso_Eff_C

SingleTau_ROCs_I = {}
SingleTau_ROCs_I[hCalo_Rate_I]   = hCalo_Eff_I
#SingleTau_ROCs_I[hTk_Rate_I]     = hTk_Eff_I
SingleTau_ROCs_I[hVtxIso_Rate_I] = hVtxIso_Eff_I

SingleTau_ROCs_F = {}
SingleTau_ROCs_F[hCalo_Rate_F]   = hCalo_Eff_F
#SingleTau_ROCs_F[hTk_Rate_F]     = hTk_Eff_F
SingleTau_ROCs_F[hVtxIso_Rate_F] = hVtxIso_Eff_F

############################################################### 
### DiTau (Indistinguishable)
############################################################### 
DiTau_TurnOn = []
DiTau_TurnOn.append(hCalo_TurnOn25)
#DiTau_TurnOn.append(hTk_TurnOn25)
DiTau_TurnOn.append(hVtxIso_TurnOn25)

DiTau_Rate = []
DiTau_Rate.append(hDiTau_Rate_Calo)
#DiTau_Rate.append(hDiTau_Rate_Tk)
DiTau_Rate.append(hDiTau_Rate_VtxIso)

DiTau_Rate_C = []
DiTau_Rate_C.append(hDiTau_Rate_Calo_C)
#DiTau_Rate_C.append(hDiTau_Rate_Tk_C)
DiTau_Rate_C.append(hDiTau_Rate_VtxIso_C)

DiTau_Rate_I = []
DiTau_Rate_I.append(hDiTau_Rate_Calo_I)
#DiTau_Rate_I.append(hDiTau_Rate_Tk_I)
DiTau_Rate_I.append(hDiTau_Rate_VtxIso_I)

DiTau_Rate_F = []
DiTau_Rate_F.append(hDiTau_Rate_Calo_F)
#DiTau_Rate_F.append(hDiTau_Rate_Tk_F)
DiTau_Rate_F.append(hDiTau_Rate_VtxIso_F)

DiTau_Eff = []
DiTau_Eff.append(hDiTau_Eff_Calo)
#DiTau_Eff.append(hDiTau_Eff_Tk)
DiTau_Eff.append(hDiTau_Eff_VtxIso)

DiTau_Eff_C = []
DiTau_Eff_C.append(hDiTau_Eff_Calo_C)
#DiTau_Eff_C.append(hDiTau_Eff_Tk_C)
DiTau_Eff_C.append(hDiTau_Eff_VtxIso_C)

DiTau_Eff_I = []
DiTau_Eff_I.append(hDiTau_Eff_Calo_I)
#DiTau_Eff_I.append(hDiTau_Eff_Tk_I)
DiTau_Eff_I.append(hDiTau_Eff_VtxIso_I)

DiTau_Eff_F = []
DiTau_Eff_F.append(hDiTau_Eff_Calo_F)
#DiTau_Eff_F.append(hDiTau_Eff_Tk_F)
DiTau_Eff_F.append(hDiTau_Eff_VtxIso_F)

DiTau_ROCs = {}
DiTau_ROCs[hDiTau_Rate_Calo]   = hDiTau_Eff_Calo
#DiTau_ROCs[hDiTau_Rate_Tk]     = hDiTau_Eff_Tk
DiTau_ROCs[hDiTau_Rate_VtxIso] = hDiTau_Eff_VtxIso

DiTau_ROCs_C = {}
DiTau_ROCs_C[hDiTau_Rate_Calo_C]   = hDiTau_Eff_Calo_C
#DiTau_ROCs_C[hDiTau_Rate_Tk_C]     = hDiTau_Eff_Tk_C
DiTau_ROCs_C[hDiTau_Rate_VtxIso_C] = hDiTau_Eff_VtxIso_C

DiTau_ROCs_I = {}
DiTau_ROCs_I[hDiTau_Rate_Calo_I]   = hDiTau_Eff_Calo_I
#DiTau_ROCs_I[hDiTau_Rate_Tk_I]     = hDiTau_Eff_Tk_I
DiTau_ROCs_I[hDiTau_Rate_VtxIso_I] = hDiTau_Eff_VtxIso_I

DiTau_ROCs_F = {}
DiTau_ROCs_F[hDiTau_Rate_Calo_F]   = hDiTau_Eff_Calo_F
#DiTau_ROCs_F[hDiTau_Rate_Tk_F]     = hDiTau_Eff_Tk_F
DiTau_ROCs_F[hDiTau_Rate_VtxIso_F] = hDiTau_Eff_VtxIso_F

############################################################### 
### DiTau (Distinguishable: Calo-Iso
############################################################### 
DiTau_Rate_CaloIso = []
DiTau_Rate_CaloIso.append(hDiTau_Rate_Calo_Tk)
DiTau_Rate_CaloIso.append(hDiTau_Rate_Calo_VtxIso)

DiTau_Eff_CaloIso = []
DiTau_Eff_CaloIso.append(hDiTau_Eff_Calo_Tk)
DiTau_Eff_CaloIso.append(hDiTau_Eff_Calo_VtxIso)

DiTau_ROCs_CaloIso = {}
DiTau_ROCs_CaloIso[hDiTau_Rate_Calo]        = hDiTau_Eff_Calo
DiTau_ROCs_CaloIso[hDiTau_Rate_Tk]          = hDiTau_Eff_Tk
DiTau_ROCs_CaloIso[hDiTau_Rate_Calo_Tk]     = hDiTau_Eff_Calo_Tk
DiTau_ROCs_CaloIso[hDiTau_Rate_Calo_VtxIso] = hDiTau_Eff_Calo_VtxIso

############################################################### 
### DiTau (Distinguishable: Tk-Iso
############################################################### 
DiTau_Rate_TkIso = []
DiTau_Rate_TkIso.append(hDiTau_Rate_Tk_VtxIso)

DiTau_Eff_TkIso = []
DiTau_Eff_TkIso.append(hDiTau_Eff_Tk_VtxIso)

DiTau_ROCs_TkIso = {}
DiTau_ROCs_TkIso[hDiTau_Rate_Calo]      = hDiTau_Eff_Calo
DiTau_ROCs_TkIso[hDiTau_Rate_Tk]        = hDiTau_Eff_Tk
DiTau_ROCs_TkIso[hDiTau_Rate_Tk_VtxIso] = hDiTau_Eff_Tk_VtxIso

DiTau_ROCs_TP = {}
DiTau_ROCs_TP[hDiTau_Rate_Calo]        = hDiTau_Eff_Calo
#DiTau_ROCs_TP[hDiTau_Rate_Tk]          = hDiTau_Eff_Tk
DiTau_ROCs_TP[hDiTau_Rate_VtxIso]      = hDiTau_Eff_VtxIso
#DiTau_ROCs_TP[hDiTau_Rate_Calo_VtxIso] = hDiTau_Eff_Calo_VtxIso

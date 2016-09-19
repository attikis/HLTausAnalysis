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
yMin       = 0.0
yMax       = None
yMinRatio  = 0.5
yMaxRatio  = 2.9
bRatio     = False
bInvRatio  = False
ratioLabel = "Ratio"
normFactor = "One"
bLogY      = False
bLogX      = False

###############################################################
### Histogram Attributes
###############################################################
kwargs   = {}

N = {"xLabel": "Multiplicity"     , "xUnits": ""       , "xMin": -0.5, "xMax": +9.5, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
     "yLabel": "Entries / %0.0f"  , "yUnits": ""       , "yMin": +0.0, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
     "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP", 
     "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

Pt = {"xLabel": "p_{T}"           , "xUnits": "GeV/c", "xMin": +0.00, "xMax": +150.0, "binWidthX": +5.0, "xCutLines": [20.0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
      "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax   , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False,
      "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
      "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

Et = {"xLabel": "E_{T}"           , "xUnits": "GeV"  , "xMin": +0.00, "xMax": +150.0, "binWidthX": None, "xCutLines": [20.0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
      "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin , "yMax": yMax  , "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False,
      "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
      "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

#Eta = {"xLabel": "#eta"            , "xUnits": ""     , "xMin": -3.0 , "xMax": +3.0 , "binWidthX": 0.1, "xCutLines": [-2.3, +2.3], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
Eta = {"xLabel": "#eta"            , "xUnits": ""     , "xMin": -3.0 , "xMax": +3.0 , "binWidthX": 0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax  , "binWidthY": None, "yCutLines": []          , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
       "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
       "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

Phi = {"xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.6, "xMax": +3.6 , "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
       "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
       "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

Mass = {"xLabel": "m_{inv}"        , "xUnits": "GeV/c^{2}" , "xMin": 0.0 , "xMax": +2.0 , "binWidthX": 0.02, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
        "yLabel": "Entries / %0.2f", "yUnits": ""        , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
        "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
        "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

Charge = {"xLabel": "Charge"         , "xUnits": "e", "xMin": -3.0 , "xMax": +3.0 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f", "yUnits": "" , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

Status = {"xLabel": "Status"         , "xUnits": "", "xMin": 0.0 , "xMax": +3.0 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin , "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VertexX = {"xLabel": "x_{vtx}"         , "xUnits": "cm", "xMin": -1.5, "xMax": +1.5, "binWidthX": None , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": yMin, "yMax": yMax, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
           "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
           "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VertexY= {"xLabel": "y_{vtx}"         , "xUnits": "cm", "xMin": -1.5, "xMax": +1.5, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.2f" , "yUnits": ""  , "yMin": yMin, "yMax": yMax, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VertexZ = {"xLabel": "z_{vtx}"         , "xUnits": "cm", "xMin": -25.0, "xMax": +25.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.1f" , "yUnits": ""  , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
           "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
           "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

DecayMode = {"xLabel": "Decay Mode"     , "xUnits": "", "xMin": +0.0, "xMax": +30.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [[1.0, +8.0, ROOT.kGreen+2]], "gridX": True, "logX": False, 
             "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": None, "yMax": None , "binWidthY": None, "yCutLines": [], "xCutBoxes": []                     , "gridY": True, "logY": bLogY, 
             "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": None, "drawOptions": "P", "legOptions": "LP",
             "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

PdgId = {"xLabel": "PdgId"          , "xUnits": "", "xMin": -600.0, "xMax": +600.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
         "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin  , "yMax": yMax  , "binWidthY": None, "yCutLines": [], "xCutBoxes": [], "gridY": True, "logY": bLogY, 
         "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" ,
         "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VisEt = {"xLabel": "E_{T}^{vis}"     , "xUnits": "GeV", "xMin": +0.00, "xMax": +150.0, "binWidthX": +5.0, "xCutLines": [20.0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
         "yLabel": "Entries / %0.0f" , "yUnits": ""   , "yMin": yMin , "yMax": yMax  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False,
         "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
         "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VisEta = {"xLabel": "#eta^{vis}"      , "xUnits": ""     , "xMin": -3.0 , "xMax": +3.0, "binWidthX": 0.04, "xCutLines": [-2.3, +2.3], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax  , "binWidthY": None, "yCutLines": []          , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VisPhi = {"xLabel": "#phi^{vis}"      , "xUnits": "rads" , "xMin": -3.6, "xMax": +3.6 , "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VisMass = {"xLabel": "m_{inv}^{vis}"  , "xUnits": "GeV/c^{2}" , "xMin": 0.0 , "xMax": +2.0 , "binWidthX": 0.02, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.2f", "yUnits": ""        , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
           "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
           "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgPt = {"xLabel": "p_{T}^{ldg}"     , "xUnits": "GeV/c", "xMin": 0.00, "xMax": 100.0, "binWidthX": +5.0, "xCutLines": [15.0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
         "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False,
         "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
         "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgEt = {"xLabel": "E_{T}^{ldg}"     , "xUnits": "GeV", "xMin": 0.00, "xMax": 100.0, "binWidthX": +5.0, "xCutLines": [20.0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
         "yLabel": "Entries / %0.0f" , "yUnits": ""   , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False,
         "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
         "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgEta = {"xLabel": "#eta^{ldg}"      , "xUnits": "", "xMin": -3.0, "xMax": +3.0, "binWidthX": +0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.2f" , "yUnits": "", "yMin": yMin, "yMax": yMax, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgPhi = {"xLabel": "#phi^{ldg}"      , "xUnits": "rads" , "xMin": -3.6, "xMax": +3.6 , "binWidthX": +0.2, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgMass = {"xLabel": "m_{inv}^{ldg}"        , "xUnits": "GeV/c^{2}" , "xMin": 0.0 , "xMax": +2.0 , "binWidthX": 0.02, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.2f", "yUnits": ""        , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
           "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
           "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgCharge = {"xLabel": "Charge^{ldg}"   , "xUnits": "e", "xMin": -3.0 , "xMax": +3.0 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.0f", "yUnits": "" , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
             "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
             "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgPdgId = {"xLabel": "PdgId^{ldg}"    , "xUnits": "", "xMin": -600.0, "xMax": +600.0, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
            "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin  , "yMax": yMax  , "binWidthY": None, "yCutLines": []                         , "xCutBoxes": [], "gridY": True, "logY": bLogY, 
            "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" ,
            "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgStatus = {"xLabel": "Status^{ldg}"   , "xUnits": "", "xMin": +0.0 , "xMax": +3.5 , "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
             "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
             "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgVertexX = {"xLabel": "x_{vtx}^{ldg}"   , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.1f" , "yUnits": ""  , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
              "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
              "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgVertexY = {"xLabel": "y_{vtx}^{ldg}"   , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.1f" , "yUnits": ""  , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
              "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
              "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

LdgVertexZ = {"xLabel": "z_{vtx}^{ldg}"   , "xUnits": "cm", "xMin": -30.0, "xMax": +30.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.1f" , "yUnits": ""  , "yMin": yMin , "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
              "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
              "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

VisPtOverPt = {"xLabel": "p_{T}^{vis}/p_{T}", "xUnits": "", "xMin": +0.00, "xMax": +2.0, "binWidthX": None, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
               "yLabel": "Entries / %0.2f"  , "yUnits": "", "yMin": yMin , "yMax": yMax, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False,
               "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
               "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

Rtau = {"xLabel": "R_{#tau} = p^{ldg} / p^{vis}_{#tau}", "xUnits": "", "xMin": 0.0 , "xMax": +1.1, "binWidthX": +0.05, "xCutLines": [] , "gridX": True, "logX": False, 
        "yLabel": "Entries / %0.2f"             , "yUnits": "", "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [] , "gridY": True, "logY": bLogY, 
        "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": 0.0, "yMaxRatio": 4.5, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
        "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

EcalFrac = {"xLabel": "E_{T}^{#pi^{0}} / E_{T}^{vis}", "xUnits": "", "xMin": 0.0, "xMax": +1.0 , "binWidthX": 0.1, "xCutLines": [] , "gridX": True, "logX": False, 
            "yLabel": "Entries / %0.2f"              , "yUnits": "", "yMin": yMin , "yMax": yMax  , "binWidthY": None , "yCutLines": [] , "gridY": True, "logY": bLogY, 
            "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
            "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

HcalFrac = {"xLabel": "E_{T}^{#pi^{#pm}} / E_{T}^{vis}", "xUnits": "", "xMin": 0.0, "xMax": +1.0 , "binWidthX": 0.1, "xCutLines": [] , "gridX": True, "logX": False, 
            "yLabel": "Entries / %0.2f"                , "yUnits": "", "yMin": yMin, "yMax": yMax  , "binWidthY": None , "yCutLines": [] , "gridY": True, "logY": bLogY, 
            "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": "One", "drawOptions": "P", "legOptions": "LP",
            "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

DeltaR = {"xLabel": "#DeltaR"        , "xUnits": "", "xMin": +0.0, "xMax": +2.0, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
          "yLabel": "Entries / %0.3f", "yUnits": "", "yMin": yMin, "yMax": yMax, "binWidthY": None, "yCutLines": [], "gridY": True, "logY": bLogY, 
          "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

DeltaRMax = {"xLabel": "#DeltaR_{max}(#pi^{#pm}_{ldg}, #pi^{#pm})", "xUnits": "", "xMin": +0.0, "xMax": +0.25, "binWidthX": None, "xCutLines": [], "gridX": True, "logX": False, 
             "yLabel": "Entries / %0.3f"                          , "yUnits": "", "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "gridY": True, "logY": bLogY, 
             "ratioLabel": ratioLabel, "ratio": False, "invRatio": False, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
             "xCutBoxes": [[0.15, 0.25, ROOT.kBlack]], "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }
             #"xCutBoxes": [[0.02, 0.10, ROOT.kBlack]], "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

DiTauPt = {"xLabel": "di-#tau p_{T}"   , "xUnits": "GeV/c", "xMin": +0.0, "xMax": +200.0, "binWidthX": +5.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
           "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax   , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False,
           "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP", 
           "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

DiTauM = {"xLabel": "m_{#tau#tau}"   , "xUnits": "GeV/c^{2}", "xMin": 0.0 , "xMax": 200.0, "binWidthX": +5.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f", "yUnits": ""       , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
          "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

DiTauVisM = {"xLabel": "m_{#tau#tau}^{vis}", "xUnits": "GeV/c^{2}", "xMin": 0.0 , "xMax": 200.0, "binWidthX": +5.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.2f", "yUnits": ""          , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": bLogY, "logYRatio": False, 
             "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP",
             "xLegMin": 0.70, "xLegMax": 0.88, "yLegMin": 0.72, "yLegMax": 0.84 }

### TH2D
DeltaRMax_Pt = {"xLabel": "#DeltaR_{max}(#pi^{#pm}_{ldg}, #pi^{#pm})  / %0.3f", "xUnits": ""    , "xMin": 0.0, "xMax":   0.25, "binWidthX": None , "gridX": False, "logX": False,
                "yLabel": "p_{T}^{ldg} / %0.0f"                                , "yUnits": "GeV", "yMin": 0.0, "yMax": 200.0 , "binWidthY": 2.0 , "gridY": False, "logY": False,
                "zLabel": "", "zUnits": ""    , "zMin": 0.0, "zMax": 0.005 , "zCutLines": []   , "gridZ": False, "logZ": False, 
                "xCutBoxes": [[0.15, 0.25, ROOT.kBlack]], "yCutLines": [15.0],
                "legOptions": "P", "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.78, "yLegMax": 0.84, "normaliseTo": "One"}

DeltaRMax_Et = {"xLabel": "#DeltaR_{max}(#pi^{#pm}_{ldg}, #pi^{#pm})  / %0.3f", "xUnits": ""   , "xMin": 0.0, "xMax":   0.25, "binWidthX": None, "gridX": False, "logX": False, 
                "yLabel": "E_{T} / %0.0f"                                     , "yUnits": "GeV", "yMin": 0.0, "yMax": 200.0 , "binWidthY": 2.0, "gridY": False, "logY": False,
                "zLabel": "", "zUnits": ""   , "zMin": 0.0, "zMax": 0.005 , "zCutLines": []  , "gridZ": False, "logZ": False,
                "xCutBoxes": [[0.15, 0.25, ROOT.kBlack]], "yCutLines": [20.0], #"yCutBoxes": [[20.0, 200.0, ROOT.kRed]], 
                "legOptions": "P", "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.78, "yLegMax": 0.84, "normaliseTo": "One"}

DeltaRMax_VisEt = {"xLabel": "#DeltaR_{max}(#pi^{#pm}_{ldg}, #pi^{#pm})  / %0.3f", "xUnits": ""   , "xMin": 0.0, "xMax":   0.25, "binWidthX": None , "gridX": False, "logX": False, 
                   "yLabel": "E_{T}^{vis} / %0.0f"                               , "yUnits": "GeV", "yMin": 0.0, "yMax": 200.0 , "binWidthY": 2.0 , "gridY": False, "logY": False,
                   "zLabel": "", "zUnits": ""   , "zMin": 0.0, "zMax": 0.005 , "zCutLines": []   , "gridZ": False, "logZ": False,
                   "xCutLines": [0.15], "yCutLines": [20.0], "xCutBoxes": [], 
                   #"xCutLines": [], "yCutLines": [20.0], "xCutBoxes": [[0.15, 0.25, ROOT.kBlack]], 
                   "legOptions": "P", "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.78, "yLegMax": 0.84, "normaliseTo": "One"}

LdgPt_VisEt = {"xLabel": "p_{T}^{ldg} / %0.0f", "xUnits": "GeV/c", "xMin": 0.00, "xMax": 100.0, "binWidthX": +2.0, "xCutLines": [15.0], "xCutBoxes": [], "gridX": False, "logX": False,
               "yLabel": "E_{T}^{vis} / %0.0f", "yUnits": "GeV"  , "yMin": 0.00, "yMax": 200.0, "binWidthY": +2.0, "yCutLines": [20.0], "xCutBoxes": [], "gridX": False, "logY": False,
               "zLabel": "", "zUnits": ""   , "zMin": 0.0, "zMax": None , "zCutLines": []   , "gridZ": False, "logZ": False,
               "legOptions": "P", "drawOptions": "COLZ", "xLegMin": 0.62, "xLegMax": 0.85, "yLegMin": 0.78, "yLegMax": 0.84, "normaliseTo": "One"}

### TF1
f1Args = {"lineColour": ROOT.kBlack, "lineStyle": ROOT.kSolid, "lineWidth": 3, "fillColour": ROOT.kBlack, "fillStyle": 0}
#f1Args = {"lineColour": ROOT.kBlack, "lineStyle": ROOT.kSolid, "lineWidth": 3, "fillColour": ROOT.kBlack, "fillStyle": 3001}

### Efficiencies
LdgPtEff    = {"xLabel": "p_{T}^{ldg}", "xUnits": "GeV", "xMin": 0.00, "xMax": +100.0, "binWidthX": None, "xCutLines": [15.0], "xCutBoxes": [], "gridX": True, "logX": False,
               "yLabel": "Efficiency %0.1f" , "yUnits": "", "yMin": 0.00, "yMax": +  1.0, "binWidthY": None, "yCutLines": []    , "yCutBoxes": [], "gridY": True, "logY": False,
               "legOptions": "LP", "drawOptions": "AP", "xLegMin": 0.60, "xLegMax": 0.84, "yLegMin": 0.20, "yLegMax": 0.46}


###############################################################
### Histograms
###############################################################
histoFolder = ""
legTitle    = " " #"#tau_{h}"

### Tau
HadronicTau_N       = m_histos.TH1orTH2( histoFolder, "HadronicTau_N"      , legTitle , None, **N       )
HadronicTau_Pt      = m_histos.TH1orTH2( histoFolder, "HadronicTau_Pt"     , legTitle , None, **Pt      )
HadronicTau_Eta     = m_histos.TH1orTH2( histoFolder, "HadronicTau_Eta"    , legTitle , None, **Eta     )
HadronicTau_Phi     = m_histos.TH1orTH2( histoFolder, "HadronicTau_Phi"    , legTitle , None, **Phi     )
HadronicTau_Mass    = m_histos.TH1orTH2( histoFolder, "HadronicTau_Mass"   , legTitle , None, **Mass    )
HadronicTau_Charge  = m_histos.TH1orTH2( histoFolder, "HadronicTau_Charge" , legTitle , None, **Charge  )
HadronicTau_PdgId   = m_histos.TH1orTH2( histoFolder, "HadronicTau_PdgId"  , legTitle , None, **PdgId   )
HadronicTau_Status  = m_histos.TH1orTH2( histoFolder, "HadronicTau_Status" , legTitle , None, **Status  )
HadronicTau_VertexX = m_histos.TH1orTH2( histoFolder, "HadronicTau_VertexX", legTitle , None, **VertexX )
HadronicTau_VertexY = m_histos.TH1orTH2( histoFolder, "HadronicTau_VertexY", legTitle , None, **VertexY )
HadronicTau_VertexZ = m_histos.TH1orTH2( histoFolder, "HadronicTau_VertexZ", legTitle , None, **VertexZ )

### Visible-Tau
HadronicTau_VisEt        = m_histos.TH1orTH2( histoFolder, "HadronicTau_VisEt"       , legTitle , None, **VisEt    )
HadronicTau_VisEta       = m_histos.TH1orTH2( histoFolder, "HadronicTau_VisEta"      , legTitle , None, **VisEta   )
HadronicTau_VisPhi       = m_histos.TH1orTH2( histoFolder, "HadronicTau_VisPhi"      , legTitle , None, **VisPhi   )
HadronicTau_VisMass      = m_histos.TH1orTH2( histoFolder, "HadronicTau_VisMass"     , legTitle , None, **VisMass  )
HadronicTau_DecayMode    = m_histos.TH1orTH2( histoFolder, "HadronicTau_DecayMode"   , legTitle , None, **DecayMode)
HadronicTau_EcalFraction = m_histos.TH1orTH2( histoFolder, "HadronicTau_EcalFraction", legTitle , None, **EcalFrac )
HadronicTau_HcalFraction = m_histos.TH1orTH2( histoFolder, "HadronicTau_HcalFraction", legTitle , None, **HcalFrac )

### Tau Charged Pions 
HadronicTau_ChargedPion_N       = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_N"      , legTitle , None, **N       )
HadronicTau_ChargedPion_Pt      = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_Pt"     , legTitle , None, **Pt      )
HadronicTau_ChargedPion_Eta     = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_Eta"    , legTitle , None, **Eta     )
HadronicTau_ChargedPion_Phi     = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_Phi"    , legTitle , None, **Phi     )
HadronicTau_ChargedPion_Mass    = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_Mass"   , legTitle , None, **Mass    )
HadronicTau_ChargedPion_Charge  = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_Charge" , legTitle , None, **Charge  )
HadronicTau_ChargedPion_PdgId   = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_PdgId"  , legTitle , None, **PdgId   )
HadronicTau_ChargedPion_Status  = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_Status" , legTitle , None, **Status  )
HadronicTau_ChargedPion_VertexX = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_VertexX", legTitle , None, **VertexX )
HadronicTau_ChargedPion_VertexY = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_VertexY", legTitle , None, **VertexY )
HadronicTau_ChargedPion_VertexZ = m_histos.TH1orTH2( histoFolder, "HadronicTau_ChargedPion_VertexZ", legTitle , None, **VertexZ )

### Tau Neutral Pions
HadronicTau_NeutralPion_N       = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_N"      , legTitle , None, **N       )
HadronicTau_NeutralPion_Pt      = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_Pt"     , legTitle , None, **Pt      )
HadronicTau_NeutralPion_Eta     = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_Eta"    , legTitle , None, **Eta     )
HadronicTau_NeutralPion_Phi     = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_Phi"    , legTitle , None, **Phi     )
HadronicTau_NeutralPion_Mass    = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_Mass"   , legTitle , None, **Mass    )
HadronicTau_NeutralPion_Charge  = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_Charge" , legTitle , None, **Charge  )
HadronicTau_NeutralPion_PdgId   = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_PdgId"  , legTitle , None, **PdgId   )
HadronicTau_NeutralPion_Status  = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_Status" , legTitle , None, **Status  )
HadronicTau_NeutralPion_VertexX = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_VertexX", legTitle , None, **VertexX )
HadronicTau_NeutralPion_VertexY = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_VertexY", legTitle , None, **VertexY )
HadronicTau_NeutralPion_VertexZ = m_histos.TH1orTH2( histoFolder, "HadronicTau_NeutralPion_VertexZ", legTitle , None, **VertexZ )

### Leading Charged Pion
HadronicTau_LdgChPion_Pt      = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Pt"      , legTitle , None, **LdgPt      )
HadronicTau_LdgChPion_Eta     = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Eta"     , legTitle , None, **LdgEta     )
HadronicTau_LdgChPion_Phi     = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Phi"     , legTitle , None, **LdgPhi     )
HadronicTau_LdgChPion_Mass    = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Mass"    , legTitle , None, **LdgMass    )
HadronicTau_LdgChPion_Charge  = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Charge"  , legTitle , None, **LdgCharge  )
HadronicTau_LdgChPion_PdgId   = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_PdgId"   , legTitle , None, **LdgPdgId   )
HadronicTau_LdgChPion_Status  = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Status"  , legTitle , None, **LdgStatus  )
HadronicTau_LdgChPion_VertexX = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_VertexX" , legTitle , None, **LdgVertexX )
HadronicTau_LdgChPion_VertexY = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_VertexY" , legTitle , None, **LdgVertexY )
HadronicTau_LdgChPion_VertexZ = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_VertexZ" , legTitle , None, **LdgVertexZ )
HadronicTau_LdgChPion_Rtau    = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Rtau"    , legTitle , None, **Rtau       )
HadronicTau_LdgChPion_Pt_VisEt= m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_Pt_VisEt", legTitle , None, **LdgPt_VisEt)

### 3-prong
HadronicTau_LdgChPion_DeltaRMax          = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax"         , legTitle , None, **DeltaRMax       )
HadronicTau_LdgChPion_DeltaRMax_Pt       = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax_Pt"      , legTitle , None, **DeltaRMax_Pt    )
HadronicTau_LdgChPion_DeltaRMax_TauEt    = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax_TauEt"   , legTitle , None, **DeltaRMax_Et    )
HadronicTau_LdgChPion_DeltaRMax_VisTauEt = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax_VisTauEt", legTitle , None, **DeltaRMax_VisEt )

HadronicTau_LdgChPion_DeltaRMax_MinPt          = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax_MinPt"         , legTitle , None, **DeltaRMax       )
HadronicTau_LdgChPion_DeltaRMax_Pt_MinPt       = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax_Pt_MinPt"      , legTitle , None, **DeltaRMax_Pt    )
HadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt    = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt"   , legTitle , None, **DeltaRMax_Et    )
HadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt = m_histos.TH1orTH2( histoFolder, "HadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt", legTitle , None, **DeltaRMax_VisEt )
  
### DiTau System
DiTau_Pt      = m_histos.TH1orTH2( histoFolder, "DiTau_Pt"     , legTitle , None, **DiTauPt   )
DiTau_InvMass = m_histos.TH1orTH2( histoFolder, "DiTau_InvMass", legTitle , None, **DiTauM    )
DiTau_VisMass = m_histos.TH1orTH2( histoFolder, "DiTau_VisMass", legTitle , None, **DiTauVisM )

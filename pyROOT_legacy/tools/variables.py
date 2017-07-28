###############################################################
### Author .........: Alexandros X. Attikis 
### Institute ......: University Of Cyprus (UCY)
### Email ..........: attikis@cern.ch
###############################################################

###############################################################
### All imported modules
###############################################################
### System modules
import sys
import array
import math
import copy
### Other 
import ROOT

###############################################################
### Histogram attributes
###############################################################
yMin     = 1E-05
yMin     = 1.0
ratio    = False
invratio = False
yMinR    = 0.0
yMaxR    = 10.0
norm     = "One"


NTks = {"xLabel": "Multiplicity"     , "xUnits": ""       , "xMin": -0.5, "xMax": 100.5, "binWidthX": 1.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
        "yLabel": "Entries / %0.0f"  , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
        "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

MatchIndex = {"xLabel": "match index"      , "xUnits": ""       , "xMin": -0.5, "xMax": 12.5, "binWidthX": 1.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f"  , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }


NTkTaus = {"xLabel": "Multiplicity"     , "xUnits": ""       , "xMin": -0.5, "xMax": 20.5, "binWidthX": 1.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.0f"  , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
           "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

N = {"xLabel": "Multiplicity"     , "xUnits": ""       , "xMin": -0.5, "xMax": 12.5, "binWidthX": 1.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
     "yLabel": "Entries / %0.0f"  , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
     "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }


Pt = {"xLabel": "p_{T}"           , "xUnits": "(GeV/c)", "xMin": 0.0  , "xMax": 100.0, "binWidthX": 2.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
      "yLabel": "Entries / %0.0f" , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
      "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

ScalarSumPt = {"xLabel": "#scale[0.70]{#sum} p_{T}^{tk}", "xUnits": "(GeV/c)", "xMin": 0.0, "xMax": 100.0, "binWidthX": 5.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, 
               "yLabel": "Entries / %0.0f" , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , 
               "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

CaloEt = {"xLabel": "E_{T}"           , "xUnits": "(GeV)", "xMin": 0.0  , "xMax": 100.0, "binWidthX": 5.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
          "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

Et = {"xLabel": "E_{T}"           , "xUnits": "(GeV)", "xMin": 0.0  , "xMax": 100.0, "binWidthX": 2.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
      "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
      "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

Eta = {"xLabel": "#eta"            , "xUnits": ""     , "xMin": -2.6, "xMax": +2.6 , "binWidthX": 0.2 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

Phi = {"xLabel": "#phi"            , "xUnits": "(rads)" , "xMin": -3.6, "xMax": +3.6 , "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

VisEt = {"xLabel": "E_{T}^{vis}" , "xUnits": "(GeV)", "xMin": 0.0 , "xMax": 200.0, "binWidthX": 5.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
      "yLabel": "Entries / %0.0f", "yUnits": ""     , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
      "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

VisEta = {"xLabel": "#eta^{vis}"     , "xUnits": ""     , "xMin": -2.6, "xMax": +2.6, "binWidthX": 0.2 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.2f", "yUnits": ""     , "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

VisPhi = {"xLabel": "#phi^{vis}"            , "xUnits": "(rads)" , "xMin": -3.6 , "xMax": +3.6 , "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

LdgEt = {"xLabel": "E_{T}^{ldg}"           , "xUnits": "(GeV)", "xMin": 0.0  , "xMax": 200.0, "binWidthX": 5.0 , "xCutLines": [], "xCutLines": [], "gridX": True, "logX": False, "logXRatio": False, 
      "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": [], "gridY": True, "logY": True , "logYRatio": False, 
      "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

LdgPt = {"xLabel": "p_{T}^{ldg}"     , "xUnits": "(GeV/c)", "xMin": 0.0 , "xMax": 100.0, "binWidthX": 5.0 , "xCutLines": [], "xCutLines": [5, 10], "gridX": True, "logX": False, "logXRatio": False, 
         "yLabel": "Entries / %0.0f" , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutLines": []     , "gridY": True, "logY": True , "logYRatio": False, 
         "ratioLabel": "1/ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

LdgEta = {"xLabel": "#eta^{ldg}"  , "xUnits": "", "xMin": -2.6, "xMax": +2.6, "binWidthX": 0.2 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

LdgPhi = {"xLabel": "#phi^{ldg}"            , "xUnits": "(rads)" , "xMin": -3.6 , "xMax": +3.6 , "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""       , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

CaloBx = {"xLabel": "bx"              , "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
          "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

CaloType = {"xLabel": "Type"            , "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
            "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
            "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

MatchFound = {"xLabel": "match found"     , "xUnits": "" , "xMin": -0.5 , "xMax": +1.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NProngs = {"xLabel": "n-Prong"         , "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
           "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NTracks = {"xLabel": "Track Multiplicity", "xUnits": "" , "xMin": -0.5 , "xMax": +8.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.0f"   , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
           "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NTracksIso = {"xLabel": "Tracks in Isolation Annulus", "xUnits": "", "xMin": -0.5, "xMax": +8.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f"            , "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NNeutrals = {"xLabel": "Neutral Multiplicity", "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.0f"      , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
             "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NGammas = {"xLabel": "Gamma Multiplicity", "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
           "yLabel": "Entries / %0.0f"    , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
           "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NStables = {"xLabel": "Stable-Particle Multiplicity", "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
            "yLabel": "Entries / %0.0f"             , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
            "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NNeutrinos = {"xLabel": "Neutrino Multiplicity", "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f"      , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NDaughters = {"xLabel": "Daughter Multiplicity", "xUnits": "" , "xMin": -0.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f"      , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DecayType = {"xLabel": "Decay Type"    , "xUnits": "" , "xMin": -0.5 , "xMax": +8.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
             "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

Status = {"xLabel": "Status"    , "xUnits": "" , "xMin": -0.5 , "xMax": +3.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
             "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

ParticleId = {"xLabel": "particle id"     , "xUnits": "" , "xMin": -50.0, "xMax": +50.0 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

MomId = {"xLabel": "mom id"          , "xUnits": "" , "xMin": -50.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [25], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
         "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
         "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

RealMomId = {"xLabel": "real mom id"     , "xUnits": "" , "xMin": -50.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [25], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
             "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
             "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

VtxZ0  = {"xLabel": "z^{vtx}"         , "xUnits": "(cm)", "xMin": -30.0, "xMax": +30.0, "binWidthX": 1.0 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f" , "yUnits": ""    , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
          "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaVtxZ0Iso = {"xLabel": "|z^{vtx}_{ldg} - z^{vtx}_{iso-tk}|", "xUnits": "(cm)", "xMin": 0.0 ,  "xMax": +5.0, "binWidthX": 0.1 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
                 "yLabel": "Entries / %0.1f"               , "yUnits": ""    , "yMin": yMin, "yMax": +1.0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
                 "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaZ_MCVertex = {"xLabel": "z^{vtx}_{pv} - z^{vtx}_{mc}", "xUnits": "(cm)", "xMin": -10.0 ,  "xMax": +10.0, "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
                   "yLabel": "Entries / %0.1f"            , "yUnits": ""    , "yMin": yMin  , "yMax": +1.0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
                   "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

AbsDeltaZ_MCVertex = {"xLabel": "|z^{vtx}_{pv} - z^{vtx}_{mc}|", "xUnits": "(cm)", "xMin": 0.0 ,  "xMax": +10.0, "binWidthX": 0.2 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
                      "yLabel": "Entries / %0.1f"              , "yUnits": ""    , "yMin": yMin, "yMax": +1.0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
                      "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaVtxZ0Sig = {"xLabel": "|z^{vtx}_{ldg} - z^{vtx}_{sig-tk}|", "xUnits": "(cm)", "xMin": 0.0 ,  "xMax": +5.0, "binWidthX": 0.1 , "xCutLines": [1.0], "xCutBoxes": [], "gridX": True, "logX": False, 
                 "yLabel": "Entries / %0.0f"               , "yUnits": ""    , "yMin": yMin, "yMax": +1.0  , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
                 "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

PVDeltaZ_LdgTk = {"xLabel": "|z^{vtx}_{ldg-tk} - z^{vtx}_{pv}|", "xUnits": "(cm)", "xMin": 0.0 ,  "xMax": +5.0, "binWidthX": 0.1 , "xCutLines": [1.0], "xCutBoxes": [], "gridX": True, "logX": False, 
                  "yLabel": "Entries / %0.1f"                  , "yUnits": ""    , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": True , 
                  "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

PVDeltaZ_SigTk = {"xLabel": "|z^{vtx}_{sig-tk} - z^{vtx}_{pv}|", "xUnits": "(cm)", "xMin": 0.0 ,  "xMax": +5.0, "binWidthX": 0.1 , "xCutLines": [1.0], "xCutBoxes": [], "gridX": True, "logX": False, 
                  "yLabel": "Entries / %0.1f"                  , "yUnits": ""    , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": True , 
                  "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

PVDeltaZ_IsoTk = {"xLabel": "|z^{vtx}_{iso-tk} - z^{vtx}_{pv}|", "xUnits": "(cm)", "xMin": 0.0 ,  "xMax": +5.0, "binWidthX": 0.1 , "xCutLines": [1.0], "xCutBoxes": [], "gridX": True, "logX": False, 
                  "yLabel": "Entries / %0.1f"                  , "yUnits": ""    , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": True , 
                  "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

ChiSquared  = {"xLabel": "#chi^{2}"        , "xUnits": ""    , "xMin": +0.0 , "xMax": +60.0, "binWidthX": 1.0 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
               "yLabel": "Entries / %0.0f" , "yUnits": ""    , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
               "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

RedChiSquared = {"xLabel": "#chi^{2}/d.o.f.", "xUnits": "", "xMin": -0.5, "xMax": +10.5, "binWidthX": 1.0 , "xCutLines": [5,8], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
                 "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []   , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
                 "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

StubPtCons = {"xLabel": "p_{T}^{stub} constistency", "xUnits": "", "xMin": 0.0  , "xMax": 200.0, "binWidthX": 5.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
              "yLabel": "Entries / %0.0f"          , "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
              "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

HasBarrelStub = {"xLabel": "has barrel stub", "xUnits": "", "xMin": -0.5 , "xMax": 1.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
                 "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
                 "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

HasEndcapStub = {"xLabel": "has endcap stub", "xUnits": "", "xMin": -0.5 , "xMax": 1.5 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
                 "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
                 "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

Sector = {"xLabel": "Sector"         , "xUnits": "", "xMin": -1000.0 , "xMax": +1000.0, "binWidthX": 10.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin   , "yMax": +1.0   , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
          "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

Wedge = {"xLabel": "Wedge"          , "xUnits": "", "xMin": -1000.0 , "xMax": +1000.0, "binWidthX": 10.0, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
         "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin   , "yMax": +1.0   , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
         "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

RInv = {"xLabel": "1/r"             , "xUnits": "(cm^{-1})", "xMin": -1E-02  , "xMax": +1E-02 , "binWidthX": 1E-03, "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
         "yLabel": "Entries / %0.0f", "yUnits": ""         , "yMin": yMin   , "yMax": +1.0   , "binWidthY": None , "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
        "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

NStubs = {"xLabel": "Stubs Multiplicity", "xUnits": "", "xMin": -0.5 , "xMax": 10.5, "binWidthX": 1.0 , "xCutLines": [5], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f"   , "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
          "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

iDisk = {"xLabel": "Stub iDisk"     , "xUnits": "", "xMin": -0.5 , "xMax": 5.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
         "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
         "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

iLayer = {"xLabel": "Stub iLayer"    , "xUnits": "", "xMin": -0.5 , "xMax": 7.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
          "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
          "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

iPhi = {"xLabel": "Stub iPhi"      , "xUnits": "", "xMin": -0.5 , "xMax": 84.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
        "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
        "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

iRing = {"xLabel": "Stub iRing"     , "xUnits": "", "xMin": -0.5 , "xMax": 16.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
         "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
         "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

iSide = {"xLabel": "Stub iSide"     , "xUnits": "", "xMin": -0.5 , "xMax": 3.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
         "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
         "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

iZ = {"xLabel": "Stub iZ   "     , "xUnits": "", "xMin": -0.5 , "xMax": 65.5, "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
      "yLabel": "Entries / %0.0f", "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
      "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

RelIso = {"xLabel": "#scale[0.70]{#sum} p_{T}^{tks} / p_{T}^{ldg}", "xUnits": "", "xMin": 0.0 , "xMax": +1.5, "binWidthX": 0.05, "xCutLines": [0.1, 0.2], "xCutBoxes": [], "gridX": True, "logX": False, 
          "yLabel": "Entries / %0.1f"                             , "yUnits": "", "yMin": yMin, "yMax": +1.0, "binWidthY": None, "yCutLines": []        , "yCutBoxes": [], "gridY": True, "logY": True, 
          "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

RelIsoInv = {"xLabel": "p_{T}^{ldg}/#scale[0.70]{#sum} p_{T}^{tks}", "xUnits": "", "xMin": 0.0 , "xMax": +30.0, "binWidthX": 0.5 , "xCutLines": [5, 10], "xCutBoxes": [], "gridX": True, "logX": False, 
             "yLabel": "Entries / %0.1f"                           , "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []     , "yCutBoxes": [], "gridY": True, "logY": True , 
             "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

IsoTracksIndex = {"xLabel": "tk^{iso} index of TkTaus", "xUnits": "", "xMin": 0.0 , "xMax": +10.0 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
                  "yLabel": "Entries / %0.1f"        , "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
                  "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

SignalTracksIndex = {"xLabel": "tk^{sig} index of TkTaus", "xUnits": "", "xMin": 0.0 , "xMax": +10.0 , "binWidthX": 1.0 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
                     "yLabel": "Entries / %0.1f"         , "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
                     "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaEta = {"xLabel": "#Delta#eta(TkTau, Calo)", "xUnits": "", "xMin": 0.0 , "xMax": +0.25, "binWidthX": 0.01 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
            "yLabel": "Entries / %0.2f"        , "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
            "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaPhi = {"xLabel": "#Delta#phi(TkTau, Calo)", "xUnits": "(rads)", "xMin": 0.0 , "xMax": +0.25, "binWidthX": 0.01 , "xCutLines": [], "xCutBoxes": [], "gridX": True, "logX": False, 
            "yLabel": "Entries / %0.2f"        , "yUnits": ""      , "yMin": yMin, "yMax": +1.0 , "binWidthY": None , "yCutLines": [], "yCutBoxes": [], "gridY": True, "logY": True , 
            "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaR = {"xLabel": "#DeltaR(TkTau, Calo)", "xUnits": "", "xMin": 0.0 , "xMax": +0.3 , "binWidthX": 0.01, "xCutLines": [0.15], "xCutBoxes": [], "gridX": True, "logX": False, 
          "yLabel": "Entries / %0.2f"     , "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": []    , "yCutBoxes": [], "gridY": True, "logY": True , 
          "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaPtEt = {"xLabel": "#scale[0.70]{#sum} p_{T}^{tks} - E_{T}^{calo}", "xUnits": "", "xMin": -65.0, "xMax": +65.0, "binWidthX": 5.0 , "xCutLines": [0], "gridX": True, "logX": False, 
             "yLabel": "Entries / %0.0f"                              , "yUnits": ""      , "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "yCutLines": [],  "gridY": True, "logY": True , 
             "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

DeltaPtEtOverEt = {"xLabel": "(#scale[0.70]{#sum} p_{T}^{tks} - E_{T}^{calo})/E_{T}^{calo}", "xUnits": "", "xMin": -3.0, "xMax": +3.0, "binWidthX": 0.2 , "gridX": True, "xCutLines": [0], 
                   "yLabel": "Entries / %0.2f"                                             , "yUnits": ""      , "yMin": yMin, "yMax": +1.0, "binWidthY": None, "gridY": True, "logY": True , 
                   "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

LdgPtOverCaloEt = {"xLabel": "p_{T}^{ldg} / E_{T}^{calo}", "xUnits": "", "xMin": +0.0, "xMax": +10.0, "binWidthX": 0.2 , "gridX": True, "xCutLines": [], 
                   "yLabel": "Entries / %0.2f"                           , "yUnits": "", "yMin": yMin, "yMax": +1.0 , "binWidthY": None, "gridY": True, "logY": True   , 
                   "ratioLabel": "ratio", "ratio": ratio, "invRatio": invratio, "yMinRatio": yMinR, "yMaxRatio": yMaxR, "normaliseTo": norm, "drawOptions": "P", "legOptions": "LP" }

### 2D
LdgPtVsRelIso =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                                , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                  "yLabel": "#scale[0.70]{#sum} p_{T}^{tks} / p_{T}^{ldg} / %0.2f", "yUnits": "", "yMin": 0.0, "yMax": 5.0  , "binWidthY": 0.05, "yCutLines": [0.1]      , "logY": False, 
                  "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : +1E+04, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsRelIsoInv =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                                , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                     "yLabel": "1 / #scale[0.70]{#sum} p_{T}^{tks} / p_{T}^{ldg} / %0.2f", "yUnits": "", "yMin": 0.0, "yMax": 30.0  , "binWidthY": 0.5, "yCutLines": [0.1]      , "logY": False, 
                     "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : +1E+04, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso0p05cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{0.05cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso0p10cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{0.10cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso0p25cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{0.25cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso0p50cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{0.50cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso0p75cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{0.75cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso1p00cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{1.00cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso1p25cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{1.25cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso1p25cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{1.25cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso1p50cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{1.50cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso1p75cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{1.75cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

LdgPtVsVtxIso2p00cm =  {"xLabel": "p_{T}^{ldg}  / %0.0f"                          , "xUnits": "", "xMin": 0.0, "xMax": +50.0, "binWidthX": 1.0 , "xCutLines": [5.0, 10.0], "logX": False,
                        "yLabel": "#scale[0.70]{#sum} tks_{iso}^{2.00cm}  / %0.1f", "yUnits": "", "yMin": 0.0, "yMax": 10.0  , "binWidthY": 1.0, "yCutLines": []         , "logY": False, 
                        "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

SumPt_Vs_ZVertex =  {"xLabel": "#scale[0.70]{#sum} p_{T}^{tk} / %0.0f", "xUnits": "(GeV/c)", "xMin": +0.00, "xMax": 100.0, "binWidthX": 5.0 , "xCutLines": [], "gridX": False, "logX": False,
                     "yLabel": "z^{vtx} /  %0.0f"                     , "yUnits": "(cm)"   , "yMin": -30.0, "yMax": +30.0, "binWidthY": 1.0 , "yCutLines": [], "gridY": False, "logY": False, 
                     "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

SumPt_Vs_AbsDeltaZMC =  {"xLabel": "#scale[0.70]{#sum} p_{T}^{tk} / %0.0f", "xUnits": "(GeV/c)", "xMin": +0.00, "xMax": 100.0, "binWidthX": 5.0 , "xCutLines": [], "gridX": False, "logX": False,
                         "yLabel": "|z^{vtx}_{pv} - z^{vtx}_{mc}| / %0.1f", "yUnits": "(cm)"   , "yMin": +0.00, "yMax": +10.0, "binWidthY": 0.2 , "yCutLines": [], "gridY": False, "logY": False, 
                         "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

SumPt_Vs_DeltaZMC =  {"xLabel": "#scale[0.70]{#sum} p_{T}^{tk} / %0.0f", "xUnits": "(GeV/c)", "xMin": +0.00, "xMax": 100.0, "binWidthX": 5.0 , "xCutLines": [], "gridX": False, "logX": False,
                      "yLabel": "z^{vtx}_{pv} - z^{vtx}_{mc} / %0.1f"  , "yUnits": "(cm)"   , "yMin": -10.0, "yMax": +10.0, "binWidthY": 0.2 , "yCutLines": [], "gridY": False, "logY": False, 
                      "zLabel": "Entries", "zUnits": "", "zMin" : 1.0, "zMax" : 1E4, "zCutLines": [], "zCutLinesErrors": False, "logZ": True, "normaliseTo": None, "drawOptions": "COLZ"}

###############################################################
### Standard cuts-lines definitions
###############################################################
EtaRegionLines  = [-2.4, -1.6, -0.8, +0.8, +1.6, +2.4]
EtaRegionsBoxes = [ [-2.4, -1.6, ROOT.kGreen+2], [+1.6, +2.4, ROOT.kGreen+2], [-1.6, -0.8, ROOT.kBlue], [+0.8, +1.6, ROOT.kBlue], [-0.8, +0.8, ROOT.kRed+2] ]

############################################################### ###############################################################

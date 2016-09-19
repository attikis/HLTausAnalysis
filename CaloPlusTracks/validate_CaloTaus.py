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
import tools.cuts as m_cuts
import tools.aux as m_aux
#from tools.variables import *

###############################################################
### Options here
###############################################################
bBatchMode     = True
bVerbose       = False
bSavePlots     = True

#inputPath      = "Macros/TauTrigger/Validation_Histograms_"
inputPath      = "test/Validation_Histograms_"
mySavePath     = ""
mySaveFormats  = ["png"]
histoFolder    = ""
datasetList    = ["MinBias"]

datasetPaths = {}
datasetPaths["MinBias"]           = inputPath + "nugun.root"
datasetPaths["VBF"]               = inputPath + "vbf.root"
datasetPaths["SinglePionPlus"]    = inputPath + "onepionplus.root"
datasetPaths["SinglePionMinus"]   = inputPath + "onepionminus.root"
datasetPaths["SingleTauOneProng"] = inputPath + "onetaugun1p.root"
datasetPaths["TauThreeProngsEnr"] = inputPath + "twotaugun3p.root"
datasetPaths["TTbar"]             = inputPath + "ttbar.root"
datasetPaths["HPlus160"]          = inputPath + "litehiggs.root"
datasetPaths["HPlus200"]          = inputPath + "heavyhiggs.root"

### Histogram Global Options
yMin       = 1E+00 #1E-04 #1E+00
yMax       = None  #1E+00 #None
yMinRatio  = 0.0
yMaxRatio  = 5.0
bRatio     = False
bInvRatio  = False
ratioLabel = "ratio"
normFactor = None #"One" #None
nCaloTaus  = 3 #12 

############################################################### 
### Histogram Options
############################################################### 
N = {"xLabel": "Multiplicity"     , "xUnits": ""       , "xMin": -0.5 , "xMax": 21.5, "binWidthX": 1.0 , "xCutLines": [20], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
     "yLabel": "Entries / %0.0f"  , "yUnits": ""       , "yMin": yMin, "yMax": yMax, "binWidthY": None, "yCutLines": []  , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
     "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" }

Energy = {"xLabel": "E"              , "xUnits": "GeV", "xMin": 0.00 , "xMax": 100.0, "binWidthX": 2.0 , "xCutLines": [5], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
          "yLabel": "Entries / %0.0f", "yUnits": ""     , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
          "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" }

Et = {"xLabel": "E_{T}"           , "xUnits": "GeV", "xMin": 0.00 , "xMax": 100.0, "binWidthX": 2.0 , "xCutLines": [5], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False,
      "yLabel": "Entries / %0.0f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False,
      "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" }

Eta = {"xLabel": "#eta"            , "xUnits": ""     , "xMin": -2.6 , "xMax": +2.6 , "binWidthX": 0.2 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""     , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" }

Phi = {"xLabel": "#phi"            , "xUnits": "rads" , "xMin": -3.6 , "xMax": +3.6 , "binWidthX": 0.2 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
       "yLabel": "Entries / %0.2f" , "yUnits": ""       , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
       "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" }

Type = {"xLabel": "Type"            , "xUnits": "" , "xMin": -0.5 , "xMax": +9.5 , "binWidthX": 1.0 , "xCutLines": [3], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
        "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
        "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" }

Bx = {"xLabel": "bx"              , "xUnits": "" , "xMin": -5.5 , "xMax": +5.5 , "binWidthX": 1.0 , "xCutLines": [0], "xCutBoxes": [], "gridX": True, "logX": False, "logXRatio": False, 
      "yLabel": "Entries / %0.0f" , "yUnits": "" , "yMin": yMin, "yMax": yMax , "binWidthY": None, "yCutLines": [] , "yCutBoxes": [], "gridY": True, "logY": True , "logYRatio": False, 
      "ratioLabel": ratioLabel, "ratio": bRatio, "invRatio": bInvRatio, "yMinRatio": yMinRatio, "yMaxRatio": yMaxRatio, "normaliseTo": normFactor, "drawOptions": "P", "legOptions": "LP" }

############################################################### 
### Group Histograms
############################################################### 
L1CaloTaus                = []
L1CaloTau_E               = []
L1CaloTau_Et              = []
L1CaloTau_Eta             = []
L1CaloTau_Phi             = []
L1CaloTau_Bx              = []
L1CaloTau_Type            = []
L1CaloTau_Multiplicity    = m_histos.TH1orTH2( histoFolder, "L1CaloTau_Multiplicity", "Calo", None, **N )

### Histograms for eT-sorted L1CaloTaus
for i in range(0, nCaloTaus, 1):
    index = str(i)
    L1CaloTau_E.append              ( m_histos.TH1orTH2( histoFolder, "L1CaloTau_E%s"    % (index), "Calo (%s)" % (index) , None, **Energy ) )
    L1CaloTau_Et.append             ( m_histos.TH1orTH2( histoFolder, "L1CaloTau_Et%s"   % (index), "Calo (%s)" % (index) , None, **Et  ) )
    L1CaloTau_Eta.append            ( m_histos.TH1orTH2( histoFolder, "L1CaloTau_Eta%s"  % (index), "Calo (%s)" % (index) , None, **Eta ) )
    L1CaloTau_Phi.append            ( m_histos.TH1orTH2( histoFolder, "L1CaloTau_Phi%s"  % (index), "Calo (%s)" % (index) , None, **Phi ) )
    L1CaloTau_Bx.append             ( m_histos.TH1orTH2( histoFolder, "L1CaloTau_Bx%s"   % (index), "Calo (%s)" % (index) , None, **Bx  ) )
    L1CaloTau_Type.append           ( m_histos.TH1orTH2( histoFolder, "L1CaloTau_Type%s" % (index), "Calo (%s)" % (index) , None, **Type) )

###############################################################
### Main
###############################################################
def DoPlots(hList, datasetList):

    p = m_plotter.Plotter( bVerbose, bBatchMode )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
    #p.GetDatasets(datasetProd, inputPath, datasetList)
    p.AddHisto(hList)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p.SaveHistos(bSavePlots, mySavePath, mySaveFormats)

    return

###############################################################
if __name__ == "__main__":

    DoPlots( [L1CaloTau_Multiplicity], datasetList )
    DoPlots( L1CaloTau_E             , datasetList )
    DoPlots( L1CaloTau_Et            , datasetList )
    DoPlots( L1CaloTau_Eta           , datasetList )
    DoPlots( L1CaloTau_Phi           , datasetList )
    DoPlots( L1CaloTau_Bx            , datasetList )
    DoPlots( L1CaloTau_Type          , datasetList )

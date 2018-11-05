#!/usr/bin/env python
'''
DESCRIPTION:
Script that plots the ROC curvs (Rate Vs Efficiency) 
for all datases in a given multicrab.


USAGE:
./plotResolutions.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotResolutions.py -m multicrab_CaloTkSkim_v92X_20180801T1203
./plotResolutions.py -m multicrab_CaloTk_v92X_IsoConeRMax0p3_VtxIso0p5_RelIso0p2_14h29m15s_23Aug2018 -e "SingleE"
./plotResolutions.py -e "SingleE" -m multicrab_TkTau_v92X_SeedPt5ChiSq50Stubs5_SigRMax0p25Const2p5Pt2ChiSq50Stub5dZ0p5_IsoRMax0p35Pt2ChiSq50Stubs4_VtxIso0p5_RelIso0p20dZ0p5_VtxIso0p5RelIso0p30_18h08m38s_15Oct2018 -n


LAST USED:
./plotResolutions.py -e "SingleE" -m 

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
import re
import array
from optparse import OptionParser

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

import HLTausAnalysis.NtupleAnalysis.tools.dataset as dataset
import HLTausAnalysis.NtupleAnalysis.tools.histograms as histograms
import HLTausAnalysis.NtupleAnalysis.tools.counter as counter
import HLTausAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HLTausAnalysis.NtupleAnalysis.tools.styles as styles
import HLTausAnalysis.NtupleAnalysis.tools.plots as plots
import HLTausAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HLTausAnalysis.NtupleAnalysis.tools.aux as aux
import HLTausAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HLTausAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

#================================================================================================
# Variable definition
#================================================================================================
ss = ShellStyles.SuccessStyle()
ns = ShellStyles.NormalStyle()
ts = ShellStyles.NoteStyle()
hs = ShellStyles.HighlightAltStyle()
ls = ShellStyles.HighlightStyle()
es = ShellStyles.ErrorStyle()
cs = ShellStyles.CaptionStyle()


#================================================================================================ 
# Function Definition
#================================================================================================ 
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return

def Verbose(msg, printHeader=True, verbose=False):
    if not opts.verbose:
        return
    Print(msg, printHeader)
    return


def GetDatasetsFromDir(opts):
    Verbose("Getting datasets")
    
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode, 
                                                        analysisName=opts.analysis,
                                                        optimizationMode=opts.optMode)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysis,
                                                        includeOnlyTasks=opts.includeOnlyTasks,
                                                        optimizationMode=opts.optMode)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysis,
                                                        excludeTasks=opts.excludeTasks,
                                                        optimizationMode=opts.optMode)
    else:
        raise Exception("This should never be reached")
    return datasets
    
def getAlgos():
    '''
    https://root.cern.ch/doc/master/classTAttText.html
    '''
    #return ["TkTaus", "TkTaus (RelIso)", "TkTaus (VtxIso)", "TkTaus (VtxIso-L)", "TkTaus (VtxIso-T)", "TkTaus (RelIso-L)", "TkTaus (RelIso-T)"]
    return ["TkTaus", "TkTaus #font[72]{RelIso}", "TkTaus #font[72]{VtxIso}", "TkTaus #font[72]{VtxIso-L}", "TkTaus #font[72]{VtxIso-T}", "TkTaus #font[72]{RelIso-L}", "TkTaus #font[72]{RelIso-T}"]

def main(opts):
    
    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)

    optModes = [""]
    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # For-loop: All opt Mode
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        if 0:
            datasetsMgr.loadLuminosities() # from lumi.json

        # Custom filtering of datasets 
        datasetsToRemove = []
        # For-loop: All dsets to be removed
        for i, d in enumerate(datasetsToRemove, 0):
            msg = "Removing dataset %s" % d
            Print(ShellStyles.WarningLabel() + msg + ShellStyles.NormalStyle(), i==0)
            datasetsMgr.remove(filter(lambda name: d in name, datasetsMgr.getAllDatasetNames()))        

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 

        # Print dataset information
        datasetsMgr.PrintInfo()
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
 
        # Define signal and background datasets names
        dsets_signal  = []
        dsets_minBias = []
        for d in datasetsMgr.getAllDatasetNames():
            if "SingleNeutrino" in d:
                dsets_minBias.append(d)
            else:
                dsets_signal.append(d)

        # ROC curve ingredients (histograms)
        resList = [
            ["TkEG_EtResolution" , "TkEG_EtResolution_1pr" , "TkEG_EtResolution_3pr" , "TkEG_EtResolution_withNeutrals" , "TkEG_EtResolution_noNeutrals" ], 
            ["TkEG_EtaResolution", "TkEG_EtaResolution_1pr", "TkEG_EtaResolution_3pr", "TkEG_EtaResolution_withNeutrals", "TkEG_EtaResolution_noNeutrals"], 
            ["TkEG_PhiResolution", "TkEG_PhiResolution_1pr", "TkEG_PhiResolution_3pr", "TkEG_PhiResolution_withNeutrals", "TkEG_PhiResolution_noNeutrals"],
            
            ["TkEG_EtResolution" , "TkEG_EtResolution_C"   , "TkEG_EtResolution_I"   , "TkEG_EtResolution_F" ],
            ["TkEG_EtaResolution", "TkEG_EtaResolution_C"  , "TkEG_EtaResolution_I"  , "TkEG_EtaResolution_F"],
            ["TkEG_PhiResolution", "TkEG_PhiResolution_C"  , "TkEG_PhiResolution_I"  , "TkEG_PhiResolution_F"],

            ["TkEG_EtResolution_F", "TkEG_EtResolution_1pr_F", "TkEG_EtResolution_3pr_F", "TkEG_EtResolution_withNeutrals_F", "TkEG_EtResolution_noNeutrals_F"],
            ["TkEG_EtaResolution_F", "TkEG_EtaResolution_1pr_F", "TkEG_EtaResolution_3pr_F", "TkEG_EtaResolution_withNeutrals_F", "TkEG_EtaResolution_noNeutrals_F"],
            ["TkEG_PhiResolution_F", "TkEG_PhiResolution_1pr_F", "TkEG_PhiResolution_3pr_F", "TkEG_PhiResolution_withNeutrals_F", "TkEG_PhiResolution_noNeutrals_F"],

            ["TkEG_EtResolution_F_withEGs", "TkEG_EtResolution_1pr_F_withEGs", "TkEG_EtResolution_3pr_F_withEGs", "TkEG_EtResolution_withNeutrals_F_withEGs", "TkEG_EtResolution_noNeutrals_F_withEGs"],
            ["TkEG_EtaResolution_F_withEGs", "TkEG_EtaResolution_1pr_F_withEGs", "TkEG_EtaResolution_3pr_F_withEGs", "TkEG_EtaResolution_withNeutrals_F_withEGs", "TkEG_EtaResolution_noNeutrals_F_withEGs"],
            ["TkEG_PhiResolution_F_withEGs", "TkEG_PhiResolution_1pr_F_withEGs", "TkEG_PhiResolution_3pr_F_withEGs", "TkEG_PhiResolution_withNeutrals_F_withEGs", "TkEG_PhiResolution_noNeutrals_F_withEGs"],
            ["TkEG_EtResolution_F_noEGs", "TkEG_EtResolution_1pr_F_noEGs", "TkEG_EtResolution_3pr_F_noEGs", "TkEG_EtResolution_withNeutrals_F_noEGs", "TkEG_EtResolution_noNeutrals_F_noEGs"],
            ["TkEG_EtaResolution_F_noEGs", "TkEG_EtaResolution_1pr_F_noEGs", "TkEG_EtaResolution_3pr_F_noEGs", "TkEG_EtaResolution_withNeutrals_F_noEGs", "TkEG_EtaResolution_noNeutrals_F_noEGs"],
            ["TkEG_PhiResolution_F_noEGs", "TkEG_PhiResolution_1pr_F_noEGs", "TkEG_PhiResolution_3pr_F_noEGs", "TkEG_PhiResolution_withNeutrals_F_noEGs", "TkEG_PhiResolution_noNeutrals_F_noEGs"],

            ["TkEG_EtResolution_F_withEGs_posEta", "TkEG_EtResolution_1pr_F_withEGs_posEta", "TkEG_EtResolution_3pr_F_withEGs_posEta", "TkEG_EtResolution_withNeutrals_F_withEGs_posEta", "TkEG_EtResolution_noNeutrals_F_withEGs_posEta"],
            ["TkEG_EtaResolution_F_withEGs_posEta", "TkEG_EtaResolution_1pr_F_withEGs_posEta", "TkEG_EtaResolution_3pr_F_withEGs_posEta", "TkEG_EtaResolution_withNeutrals_F_withEGs_posEta", "TkEG_EtaResolution_noNeutrals_F_withEGs_posEta"],
            ["TkEG_PhiResolution_F_withEGs_posEta", "TkEG_PhiResolution_1pr_F_withEGs_posEta", "TkEG_PhiResolution_3pr_F_withEGs_posEta", "TkEG_PhiResolution_withNeutrals_F_withEGs_posEta", "TkEG_PhiResolution_noNeutrals_F_withEGs_posEta"],
            ["TkEG_EtResolution_F_noEGs_posEta", "TkEG_EtResolution_1pr_F_noEGs_posEta", "TkEG_EtResolution_3pr_F_noEGs_posEta", "TkEG_EtResolution_withNeutrals_F_noEGs_posEta", "TkEG_EtResolution_noNeutrals_F_noEGs_posEta"],
            ["TkEG_EtaResolution_F_noEGs_posEta", "TkEG_EtaResolution_1pr_F_noEGs_posEta", "TkEG_EtaResolution_3pr_F_noEGs_posEta", "TkEG_EtaResolution_withNeutrals_F_noEGs_posEta", "TkEG_EtaResolution_noNeutrals_F_noEGs_posEta"],
            ["TkEG_PhiResolution_F_noEGs_posEta", "TkEG_PhiResolution_1pr_F_noEGs_posEta", "TkEG_PhiResolution_3pr_F_noEGs_posEta", "TkEG_PhiResolution_withNeutrals_F_noEGs_posEta", "TkEG_PhiResolution_noNeutrals_F_noEGs_posEta"],

            ["TkEG_EtResolution_F_withEGs_negEta", "TkEG_EtResolution_1pr_F_withEGs_negEta", "TkEG_EtResolution_3pr_F_withEGs_negEta", "TkEG_EtResolution_withNeutrals_F_withEGs_negEta", "TkEG_EtResolution_noNeutrals_F_withEGs_negEta"],
            ["TkEG_EtaResolution_F_withEGs_negEta", "TkEG_EtaResolution_1pr_F_withEGs_negEta", "TkEG_EtaResolution_3pr_F_withEGs_negEta", "TkEG_EtaResolution_withNeutrals_F_withEGs_negEta", "TkEG_EtaResolution_noNeutrals_F_withEGs_negEta"],
            ["TkEG_PhiResolution_F_withEGs_negEta", "TkEG_PhiResolution_1pr_F_withEGs_negEta", "TkEG_PhiResolution_3pr_F_withEGs_negEta", "TkEG_PhiResolution_withNeutrals_F_withEGs_negEta", "TkEG_PhiResolution_noNeutrals_F_withEGs_negEta"],
            ["TkEG_EtResolution_F_noEGs_negEta", "TkEG_EtResolution_1pr_F_noEGs_negEta", "TkEG_EtResolution_3pr_F_noEGs_negEta", "TkEG_EtResolution_withNeutrals_F_noEGs_negEta", "TkEG_EtResolution_noNeutrals_F_noEGs_negEta"],
            ["TkEG_EtaResolution_F_noEGs_negEta", "TkEG_EtaResolution_1pr_F_noEGs_negEta", "TkEG_EtaResolution_3pr_F_noEGs_negEta", "TkEG_EtaResolution_withNeutrals_F_noEGs_negEta", "TkEG_EtaResolution_noNeutrals_F_noEGs_negEta"],
            ["TkEG_PhiResolution_F_noEGs_negEta", "TkEG_PhiResolution_1pr_F_noEGs_negEta", "TkEG_PhiResolution_3pr_F_noEGs_negEta", "TkEG_PhiResolution_withNeutrals_F_noEGs_negEta", "TkEG_PhiResolution_noNeutrals_F_noEGs_negEta"],

            ["TkEG_ChargedResolution", "TkEG_ChargedResolution_1pr", "TkEG_ChargedResolution_3pr", "TkEG_ChargedResolution_withNeutrals", "TkEG_ChargedResolution_noNeutrals"],
            ["TkEG_NeutralsResolution", "TkEG_NeutralsResolution_1pr", "TkEG_NeutralsResolution_3pr", "TkEG_NeutralsResolution_withNeutrals", "TkEG_NeutralsResolution_noNeutrals"],

            ]
                       
    
        # For-loop: All signal histos
        for i, s in enumerate(dsets_signal, 1):
            PU = s.split("PU")[1]
            
            # By decay mode
            PlotHistos(datasetsMgr, resList[0], s, PU, "ResolutionEt_%s_all" % (s) )
            PlotHistos(datasetsMgr, resList[1], s, PU, "ResolutionEta_%s_all" % (s) )
            PlotHistos(datasetsMgr, resList[2], s, PU, "ResolutionPhi_%s_all" % (s) )
            
            # By eta region (central, intermediate, forward)
            PlotHistos(datasetsMgr, resList[3], s, PU, "ResolutionEt_%s_cif"  % (s) )
            PlotHistos(datasetsMgr, resList[4], s, PU, "ResolutionEta_%s_cif" % (s) )
            PlotHistos(datasetsMgr, resList[5], s, PU, "ResolutionPhi_%s_cif" % (s) )
            
            PlotHistos(datasetsMgr, resList[6], s, PU, "ResolutionEt_%s_all_F" % (s) )
            PlotHistos(datasetsMgr, resList[7], s, PU, "ResolutionEta_%s_all_F" % (s) )
            PlotHistos(datasetsMgr, resList[8], s, PU, "ResolutionPhi_%s_all_F" % (s) )

            PlotHistos(datasetsMgr, resList[9] , s, PU, "ResolutionEt_%s_all_F_withEGs" % (s) )
            PlotHistos(datasetsMgr, resList[10], s, PU, "ResolutionEta_%s_all_F_withEGs" % (s) )
            PlotHistos(datasetsMgr, resList[11], s, PU, "ResolutionPhi_%s_all_F_withEGs" % (s) )
            
            PlotHistos(datasetsMgr, resList[12], s, PU, "ResolutionEt_%s_all_F_noEGs" % (s) )
            PlotHistos(datasetsMgr, resList[13], s, PU, "ResolutionEta_%s_all_F_noEGs" % (s) )
            PlotHistos(datasetsMgr, resList[14], s, PU, "ResolutionPhi_%s_all_F_noEGs" % (s) )

            PlotHistos(datasetsMgr, resList[15], s, PU, "ResolutionEt_%s_all_F_withEGs_posEta" % (s) )
            PlotHistos(datasetsMgr, resList[16], s, PU, "ResolutionEta_%s_all_F_withEGs_posEta" % (s) )
            PlotHistos(datasetsMgr, resList[17], s, PU, "ResolutionPhi_%s_all_F_withEGs_posEta" % (s) )
            
            PlotHistos(datasetsMgr, resList[18], s, PU, "ResolutionEt_%s_all_F_noEGs_posEta" % (s) )
            PlotHistos(datasetsMgr, resList[19], s, PU, "ResolutionEta_%s_all_F_noEGs_posEta" % (s) )
            PlotHistos(datasetsMgr, resList[20], s, PU, "ResolutionPhi_%s_all_F_noEGs_posEta" % (s) )

            PlotHistos(datasetsMgr, resList[21], s, PU, "ResolutionEt_%s_all_F_withEGs_negEta" % (s) )
            PlotHistos(datasetsMgr, resList[22], s, PU, "ResolutionEta_%s_all_F_withEGs_negEta" % (s) )
            PlotHistos(datasetsMgr, resList[23], s, PU, "ResolutionPhi_%s_all_F_withEGs_negEta" % (s) )
            
            PlotHistos(datasetsMgr, resList[24], s, PU, "ResolutionEt_%s_all_F_noEGs_negEta" % (s) )
            PlotHistos(datasetsMgr, resList[25], s, PU, "ResolutionEta_%s_all_F_noEGs_negEta" % (s) )
            PlotHistos(datasetsMgr, resList[26], s, PU, "ResolutionPhi_%s_all_F_noEGs_negEta" % (s) )

            PlotHistos(datasetsMgr, resList[27], s, PU, "ResolutionCharged_%s_all" % (s) )
            PlotHistos(datasetsMgr, resList[28], s, PU, "ResolutionNeutrals_%s_all" % (s) )

        print

    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return


def PlotHistos(datasetsMgr, histoList, signal, PU, saveName=None):
    
    # Get Histogram name and its kwargs
    kwargs    = GetHistoKwargs(saveName, opts)
    hList     = []
    legDict   = {}
    algos     = getAlgos()
    if "_all"  in saveName:
        algos = ["Inclusive", "1-prong", "3-prong", "#geq 1 #pi^{0}'s", "0 #pi^{0}'s"]
        if "_all_F"  in saveName:
            algos = ["Inclusive (F)", "1-prong (F)", "3-prong (F)", "#geq 1 #pi^{0}'s (F)", "0 #pi^{0}'s (F)"]
    if "_cif"  in saveName:
        #algos = ["Inclusive", "Central", "Intermediate", "Forward"]
        algos = ["Inclusive", "|#eta| < 0.8 (C)", "0.8 < |#eta| < 1.6 (I)", "|#eta| > 1.6 (F)"]

    # For-loop: All tau algorithms
    for l, hName in enumerate(histoList, 0):
        
        algo = hName.split("_")[0]
        msg  = "Resolution for \"%s\" algorithm (%s)" % (algo, signal)
        aux.PrintFlushed(msg, False)
        h = datasetsMgr.getDataset(signal).getDatasetRootHisto(hName).getHistogram()
        h.SetName(hName)
        legDict[hName] = algos[l]
        hList.append(h)

        # Create the rate histograms
        if opts.normalizeToOne:
            p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[])
            norm = True
            for hist in p.histoMgr.getHistos():
                if hist.getRootHisto().Integral() == 0:
                    norm = False
                    break
            if (norm):
                p.histoMgr.forEachHisto(lambda h: h.getRootHisto().Scale(1.0/h.getRootHisto().Integral()) )
            else:
                aux.Print("Cannot normalise empty histo \"%s\" for dataset \"%s\"" % (hName, signal), True)
        else:
            p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[])

    # Set legend labels
    for i, h in enumerate(p.histoMgr.getHistos(), 0):
        hName = h.getName()
        p.histoMgr.forHisto(hName, styles.getCaloStyle(i))

        if "chargedresolution" in hName.lower() or "neutralsresolution" in hName.lower():
            p.histoMgr.setHistoDrawStyle(hName, "HIST")
            p.histoMgr.setHistoLegendStyle(hName, "L")
        else:
            p.histoMgr.setHistoDrawStyle(hName, "AP")
            p.histoMgr.setHistoLegendStyle(hName, "P")
            
    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legDict)

    # Draw and save the plot
    plots.drawPlot(p, saveName, **kwargs)

    # Add additional canvas text
    histograms.addPileupText("PU=%s" % (PU) )
    histograms.addText(0.22, 0.86, plots._legendLabels[signal], 17)

    # Save the plots in custom list of saveFormats
    aux.SavePlot(p, opts.saveDir, saveName, opts.saveFormats, True)
    return

def GetHistoKwargs(h, opts):
    _mvLeg1 = {"dx": -0.15, "dy": -0.00, "dh": -0.0}
    _mvLeg2 = {"dx": -0.00, "dy": -0.00, "dh": -0.0}
    _mvLeg3 = {"dx": -0.05, "dy": -0.00, "dh": -0.0}

    logY    = False
    yMin    = 0.0
    if logY:
        yMin = 1
        yMaxF = 10
    else:
        yMaxF = 1.2

    if opts.normalizeToOne:
        _yLabel = "Arbitrary Units / %.2f"
    else:
        _yLabel = "Events / %.2f"
        
    _kwargs = {
        "xlabel"           : "#delta x / x",
        "ylabel"           : _yLabel,
        "addMCUncertainty" : False, 
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "opts"             : {"xmin": 0.0, "xmax": 0.8, "ymin": yMin, "ymax":1000, "ymaxfactor": yMaxF},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : logY,
        "rebinX"           : 5,
        "moveLegend"       : _mvLeg2,
        "xtitlesize"       : 0.1,#xlabelSize,
        "ytitlesize"       : 0.1,#ylabelSize,
        "cutBox"           : {"cutValue":  0, "fillColor": 16, "box": False, "line": True , "cutGreaterThan": False},
        "cutBoxY"          : {"cutValue": 50, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
        }

    if "resolutionet_" in h.lower():
        #_kwargs["xlabel"]     = "#deltaE_{T} / E_{T}^{vis}"
        #_kwargs["opts"]       = {"xmin": -0.2, "xmax": 0.2, "ymin": 0, "ymaxfactor": yMaxF}
        #if "_all_f" in h.lower() or "_cif" in h.lower():
        _kwargs["xlabel"]     = "#deltaE_{T} / E_{T}^{vis}"
        _kwargs["opts"]       = {"xmin": -1.0, "xmax": 1.0, "ymin": 0, "ymaxfactor": yMaxF}
        _kwargs["rebinX"]     = 50
        #_kwargs["moveLegend"] = _mvLeg1
    if "resolutioneta_" in h.lower():
        _kwargs["xlabel"]     = "#delta#eta / #eta^{vis}"
        #_kwargs["opts"]       = {"xmin": -3.2, "xmax": 3.2, "ymin": 0.8e-4, "ymaxfactor": yMaxF}
        _kwargs["opts"]       = {"xmin": -0.1, "xmax": 0.1, "ymin": 0.8e-4, "ymaxfactor": yMaxF}
        #_kwargs["rebinX"]     = 1
    if "resolutionphi_" in h.lower():
        _kwargs["xlabel"]     = "#delta#phi / #phi^{vis}"
        _kwargs["opts"]       = {"xmin": -0.1, "xmax": 0.1, "ymin": 0.8e-4, "ymaxfactor": yMaxF}
        #_kwargs["opts"]       = {"xmin": -2.2, "xmax": 2.2, "ymin": 0.8e-4, "ymaxfactor": yMaxF}
        #_kwargs["rebinX"]     = 1

    if "resolutionneutrals_" in h.lower():
        _kwargs["xlabel"]     = "#delta#pi^{0}"
        _kwargs["opts"]       = {"xmin": -5, "xmax": 5, "ymin": 0, "ymaxfactor": yMaxF}
        _kwargs["rebinX"]     = 1

    if "resolutioncharged_" in h.lower():
        _kwargs["xlabel"]     = "#delta#pi^{#pm}"
        _kwargs["opts"]       = {"xmin": -5, "xmax": 5, "ymin": 0, "ymaxfactor": yMaxF}
        _kwargs["rebinX"]     = 1
        

        #if "_all_f" in h.lower() or "_cif" in h.lower():
        #    _kwargs["xlabel"]     = "#deltaE_{T} / E_{T}^{vis}"
        #    _kwargs["opts"]       = {"xmin": -1.0, "xmax": 1.0, "ymin": 0, "ymaxfactor": yMaxF}
        #_kwargs["moveLegend"] = _mvLeg1

    if "cif" in h.lower():
        _kwargs["moveLegend"] = _mvLeg1
    if "all_f" in h.lower():
        _kwargs["moveLegend"] = _mvLeg3

    return _kwargs

def getHistos(datasetsMgr, histoName):

    h1 = datasetsMgr.getDataset("Data").getDatasetRootHisto(histoName)
    h1.setName("Data")

    h2 = datasetsMgr.getDataset("EWK").getDatasetRootHisto(histoName)
    h2.setName("EWK")
    return [h1, h2]

#================================================================================================ 
# Main
#================================================================================================ 
if __name__ == "__main__":
    '''
    https://docs.python.org/3/library/argparse.html
 
    name or flags...: Either a name or a list of option strings, e.g. foo or -f, --foo.
    action..........: The basic type of action to be taken when this argument is encountered at the command line.
    nargs...........: The number of command-line arguments that should be consumed.
    const...........: A constant value required by some action and nargs selections.
    default.........: The value produced if the argument is absent from the command line.
    type............: The type to which the command-line argument should be converted.
    choices.........: A container of the allowable values for the argument.
    required........: Whether or not the command-line option may be omitted (optionals only).
    help............: A brief description of what the argument does.
    metavar.........: A name for the argument in usage messages.
    dest............: The name of the attribute to be added to the object returned by parse_args().
    '''
    
    # Default Settings
    ANALYSIS     = "HLTausAnalysis"
    BATCHMODE    = True
    DATAERA      = "TDR2019" #"ID2017" #"TDR2019"
    FOLDER       = ""
    GRIDX        = False
    GRIDY        = False
    OPTMODE      = None
    PRECISION    = 3
    RATIO        = False
    SAVEDIR      = None
    SAVEFORMATS = [".pdf"] #[".C", ".png", ".pdf"]
    SEARCHMODE   = None
    URL          = False
    VERBOSE      = False
    NORMTOONE    = True

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-n", "--normalizeToOne", dest="normalizeToOne", action="store_true", default=NORMTOONE,
                  help="Normalise all histograms to unit area? [default: %s]" % (NORMTOONE) )

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysis", dest="analysis", type="string", default=ANALYSIS,
                      help="Override default analysis name [default: %s]" % ANALYSIS)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX, 
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY, 
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)

    parser.add_option("--formats", dest="formats", default = None,
                      help="Formats in which all plots will be saved in. Provide as list of comma-separated (NO SPACE!) formats. [default: None]")
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    parser.add_option("--ratio", dest="ratio", action="store_true", default = RATIO,
                      help="Draw ratio canvas for Data/MC curves? [default: %s]" % (RATIO) )

    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if opts.mcrab == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)

    # Determine path for saving plots
    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/TkEG/", postfix="Resolutions")
    else:
        print "opts.saveDir = ", opts.saveDir

    # Overwrite default save formats?
    if opts.formats != None:
        opts.saveFormats = opts.formats.split(",")
    else:
        opts.saveFormats = SAVEFORMATS

    # Sanity check
    allowedFolders = [""]

    if opts.folder not in allowedFolders:
        Print("Invalid folder \"%s\"! Please select one of the following:" % (opts.folder), True)
        for m in allowedFolders:
            Print(m, False)
        sys.exit()
    
    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotResolutions.py: Press any key to quit ROOT ...")

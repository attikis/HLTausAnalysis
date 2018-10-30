#!/usr/bin/env python
'''
DESCRIPTION:
Script that plots the ROC curvs (Rate Vs Efficiency) 
for all datases in a given multicrab.


USAGE:
./plotRateVsEff.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotRateVsEff.py -m multicrab_CaloTkSkim_v92X_20180801T1203
./plotRateVsEff.py -m multicrab_CaloTk_v92X_IsoConeRMax0p3_VtxIso0p5_RelIso0p2_14h29m15s_23Aug2018 -e "SingleE"


LAST USED:
./plotRateVsEff.py -e "SingleE" -m 

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
        effLists    = [["Tk_Eff", "RelIso_Eff", "VtxIso_Eff", "VtxIsoLoose_Eff", "VtxIsoTight_Eff", "RelIsoLoose_Eff", "RelIsoTight_Eff"], 
                       ["DiTau_Eff_Tk", "DiTau_Eff_RelIso", "DiTau_Eff_VtxIso", "DiTau_Eff_VtxIsoLoose", "DiTau_Eff_VtxIsoTight", "DiTau_Eff_RelIsoLoose", "DiTau_Eff_RelIsoTight", ]]

        rateLists   = [["Tk_Rate", "RelIso_Rate", "VtxIso_Rate", "VtxIsoLoose_Rate", "VtxIsoTight_Rate", "RelIsoLoose_Rate", "RelIsoTight_Rate"], 
                       ["DiTau_Rate_Tk", "DiTau_Rate_RelIso", "DiTau_Rate_VtxIso", "DiTau_Rate_VtxIsoLoose", "DiTau_Rate_VtxIsoTight", "DiTau_Rate_RelIsoLoose", "DiTau_Rate_RelIsoTight"]]

        turnOnLists = [["Tk_TurnOn25", "RelIso_TurnOn25", "VtxIso_TurnOn25"],# "VtxIsoLoose_TurnOn25", "VtxIsoTight_TurnOn25", "RelIsoLoose_TurnOn25", "RelIsoTight_TurnOn25"], 
                       ["Tk_TurnOn50", "RelIso_TurnOn50", "VtxIso_TurnOn50"]]#, "VtxIsoLoose_TurnOn50", "VtxIsoTight_TurnOn50", "RelIsoLoose_TurnOn50", "RelIsoTight_TurnOn50"]]

        turnOnLists_noNeutrals = [["Tk_TurnOn25_noNeutrals", "RelIso_TurnOn25_noNeutrals", "VtxIso_TurnOn25_noNeutrals"], 
                                  ["Tk_TurnOn50_noNeutrals", "RelIso_TurnOn50_noNeutrals", "VtxIso_TurnOn50_noNeutrals"]]

        turnOnLists_withNeutrals = [["Tk_TurnOn25_withNeutrals", "RelIso_TurnOn25_withNeutrals", "VtxIso_TurnOn25_withNeutrals"], 
                                    ["Tk_TurnOn50_withNeutrals", "RelIso_TurnOn50_withNeutrals", "VtxIso_TurnOn50_withNeutrals"]]

        turnOnLists_1pr = [["Tk_TurnOn25_1pr", "RelIso_TurnOn25_1pr", "VtxIso_TurnOn25_1pr"], 
                           ["Tk_TurnOn50_1pr", "RelIso_TurnOn50_1pr", "VtxIso_TurnOn50_1pr"]]

        turnOnLists_3pr = [["Tk_TurnOn25_3pr", "RelIso_TurnOn25_3pr", "VtxIso_TurnOn25_3pr"], 
                           ["Tk_TurnOn50_3pr", "RelIso_TurnOn50_3pr", "VtxIso_TurnOn50_3pr"]]

        turnOnLists_all = [["VtxIso_TurnOn25", "VtxIso_TurnOn25_1pr", "VtxIso_TurnOn25_3pr", "VtxIso_TurnOn25_withNeutrals", "VtxIso_TurnOn25_noNeutrals"],
                           ["VtxIso_TurnOn50", "VtxIso_TurnOn50_1pr", "VtxIso_TurnOn50_3pr", "VtxIso_TurnOn50_withNeutrals", "VtxIso_TurnOn50_noNeutrals"]]
    
        # For-loop: All background histos (min bias)
        for i, b in enumerate(dsets_minBias, 1):
            bPU = b.split("PU")[1]

            # Create rate plots (SingleTau, DiTau)
            if 1:
                opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="Rates")
                PlotRate(datasetsMgr, rateLists[0], b, bPU)
                PlotRate(datasetsMgr, rateLists[1], b, bPU)
            
            # For-loop: All signal histos
            for j, s in enumerate(dsets_signal, 1):
                sPU = s.split("PU")[1]
             
                # Create rate plots (SingleTau, DiTau)
                if i == 1: # (since inside minBias loop)
                    opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="Efficiencies")
                    PlotEfficiency(datasetsMgr, effLists[0], s, sPU)
                    PlotEfficiency(datasetsMgr, effLists[1], s, sPU)

                # Skip non-matching signal and bkg PU pairs?
                if sPU != bPU and sPU != "":
                    continue
                else:
                    if sPU == "":
                        sPU = "0" #rename before saving

                # For-loop: All triggers
                for k in range(0, len(effLists)):
                    eff  = effLists[k]
                    rate = rateLists[k]
                    Verbose("Bkg = %s, Signal = %s" % (b, s), False)
                    if 1:
                        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="ROC")
                        PlotRateVsEff(datasetsMgr, eff, rate, s, b, sPU, bPU)

        # For-loop: All signal histos
        for i, s in enumerate(dsets_signal, 1):
            PU = s.split("PU")[1]
            
            # Create rate plots (SingleTau, DiTau) 
            if 1: 
                opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="TurnOns")
                PlotTurnOns(datasetsMgr, turnOnLists[0], s, PU, "TurnOns_25GeV_%s_Inclusive" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists[1], s, PU, "TurnOns_50GeV_%s_Inclusive" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_noNeutrals[0], s, PU, "TurnOns_25GeV_%s_noNeutrals" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_noNeutrals[1], s, PU, "TurnOns_50GeV_%s_noNeutrals" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_withNeutrals[0], s, PU, "TurnOns_25GeV_%s_withNeutrals" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_withNeutrals[1], s, PU, "TurnOns_50GeV_%s_withNeutrals" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_1pr[0], s, PU, "TurnOns_25GeV_%s_1pr" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_1pr[1], s, PU, "TurnOns_50GeV_%s_1pr" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_3pr[0], s, PU, "TurnOns_25GeV_%s_3pr" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_3pr[1], s, PU, "TurnOns_50GeV_%s_3pr" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_all[0], s, PU, "TurnOns_25GeV_%s_all" % (s) )
                PlotTurnOns(datasetsMgr, turnOnLists_all[1], s, PU, "TurnOns_50GeV_%s_all" % (s) )
        print

    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return


def PlotRate(datasetsMgr, histoList, bkg, PU):

    # Get Histogram name and its kwargs
    taus     = "SingleTau"
    saveName = "Rate_%s_PU%s" % (taus, PU)
    kwargs   = GetHistoKwargs(saveName, opts)
    hList    = []
    legDict  = {}
    algos    = getAlgos()

    # For-loop: All tau algorithms
    for i, hName in enumerate(histoList, 0):
        algo = hName.split("_")[0]
        if algo == "DiTau":
            taus = "DiTau"
            algo = hName.split("_")[-1]
        
        aux.PrintFlushed("Plotting rate (%s-%s)" % (algo, taus), False)
        h = datasetsMgr.getDataset(bkg).getDatasetRootHisto(hName).getHistogram()
        h.SetName(hName)
        legDict[hName] = algos[i]
        hList.append(h)
        
    # Create the rate histograms
    p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[])

    # Set legend labels
    for i, h in enumerate(p.histoMgr.getHistos(), 0):
        hName = h.getName()
        algo  = h.getName().split("_")[0]
        if algo == "DiTau":
            algo = hName.split("_")[-1]
        p.histoMgr.forHisto(hName, styles.getTauAlgoStyle(algo))
        p.histoMgr.setHistoDrawStyle(hName, "HIST")
        p.histoMgr.setHistoLegendStyle(hName, "L")
    
    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legDict)

    # Draw and save the plot
    saveName = "Rate_%s_PU%s" % (taus, PU)
    plots.drawPlot(p, saveName, **kwargs)

    # Add additional canvas text
    histograms.addPileupText("PU=%s" % (PU) )
    histograms.addText(0.66, 0.86, taus, 17)

    # Save the plots in custom list of saveFormats
    aux.SavePlot(p, opts.saveDir, saveName, opts.saveFormats, True)
    #print
    return

def PlotEfficiency(datasetsMgr, histoList, signal, PU):

    # Get Histogram name and its kwargs
    taus     = "SingleTau"
    saveName = "Efficiency_%s_%s" % (taus, signal)
    kwargs   = GetHistoKwargs(saveName, opts)
    hList    = []
    legDict  = {}
    algos    = getAlgos()

    # For-loop: All tau algorithms
    count = -1
    for i, hName in enumerate(histoList, 0):
        algo = hName.split("_")[0]
        if algo == "DiTau":

            if "TT_" in signal or "GluGlu" in signal:
                pass
            else:
                return

            taus = "DiTau"
            algo = hName.split("_")[-1]
        count+=1  
        aux.PrintFlushed("Plotting efficiency (%s-%s-%s)" % (algo, taus, signal), False) #count==0)
        #aux.Print("Plotting efficiency (%s-%s-%s)" % (algo, taus, signal), False) #count==0)
        h = datasetsMgr.getDataset(signal).getDatasetRootHisto(hName).getHistogram()
        h.SetName(hName)
        legDict[hName] = algos[i]
        hList.append(h)

    # Create the rate histograms
    p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[])

    # Set legend labels
    for h in p.histoMgr.getHistos():
        hName = h.getName()
        algo  = h.getName().split("_")[0]
        if algo == "DiTau":
            algo = hName.split("_")[-1]
        p.histoMgr.forHisto(hName, styles.getTauAlgoStyle(algo))
        p.histoMgr.setHistoDrawStyle(hName, "HIST")
        p.histoMgr.setHistoLegendStyle(hName, "L")
    
    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legDict)

    # Draw and save the plot
    saveName = "Efficiency_%s_%s_%s" % (taus, algo, signal)
    plots.drawPlot(p, saveName, **kwargs)

    # Add additional canvas text
    histograms.addPileupText("PU=%s" % (PU) )
    histograms.addText(0.66, 0.86, plots._legendLabels[signal], 17)

    # Save the plots in custom list of saveFormats
    aux.SavePlot(p, opts.saveDir, saveName, opts.saveFormats, True)
    #print
    return


def PlotTurnOns(datasetsMgr, histoList, signal, PU, saveName=None):
    
    # Get Histogram name and its kwargs
    myRegex   = "(?:TurnOn)(.*)"
    m         = re.search(myRegex, histoList[0])
    threshold = m.group(1)
    if saveName==None:
        #saveName  = "TurnOns_%sGeV_%s" % (threshold, signal)
        saveName  = "TurnOns_%s_%s" % (threshold, signal)
    kwargs    = GetHistoKwargs(saveName, opts)
    hList     = []
    legDict   = {}
    algos     = getAlgos()
    if "_all"  in saveName:
        algos = ["Inclusive", "1-prong", "3-prong", "#geq 1 #pi^{0}'s", "0 #pi^{0}'s"]

    # For-loop: All tau algorithms
    for l, hName in enumerate(histoList, 0):
        
        algo = hName.split("_")[0]
        msg  = "Turn-on for \"%s\" algorithm (%s)" % (algo, signal)
        aux.PrintFlushed(msg, False)
        #aux.Print(msg, True)
        h = datasetsMgr.getDataset(signal).getDatasetRootHisto(hName).getHistogram()
        h.SetName(hName)
        legDict[hName] = algos[l]
        hList.append(h)

    # Create the rate histograms
    p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[])

    # Set legend labels
    for i, h in enumerate(p.histoMgr.getHistos(), 0):
        hName = h.getName()
        algo  = h.getName().split("_")[0]
        if algo == "DiTau":
            algo = hName.split("_")[-1]
        if "_all" in saveName:
            p.histoMgr.forHisto(hName, styles.getCaloStyle(i))
        else:
            p.histoMgr.forHisto(hName, styles.getTauAlgoStyle(algo))
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

def PlotRateVsEff(datasetsMgr, effHistoList, rateHistoList, signal, bkg, sPU, bPU):

    # Definitions
    tgraphs = []
    legDict = {}
    aux.PrintFlushed("Plotting ROC (%s-%s)" % (bkg, signal), False) #count==0)

    # Get Histogram name and its kwargs
    if "ditau" in effHistoList[0].lower():
        if "TT_" in signal or "GluGlu" in signal:
            pass
        else:
            return
        saveName = "RateVsEff_DiTau_%s_PU%s_MBPU%s" % (signal.split("_")[0], sPU, bPU)
    else:
        saveName = "RateVsEff_SingleTau_%s_PU%s_MBPU%s" % (signal.split("_")[0], sPU, bPU)
    kwargs_  = GetHistoKwargs(saveName, opts)

    for i in range (0, len(effHistoList)):
        if (i==0):
            g1 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], signal, bkg)
            g1.SetName("Tk")
        elif (i==1):
            g2 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], signal, bkg)
            g2.SetName("RelIso")
        elif (i==2):
            g3 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], signal, bkg)
            g3.SetName("VtxIso")
        elif (i==3):
            g4 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], signal, bkg)
            g4.SetName("VtxIsoLoose")
        elif (i==4):
            g5 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], signal, bkg)
            g5.SetName("VtxIsoTight")
        elif (i==5):
            g6 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], signal, bkg)
            g6.SetName("RelIsoLoose")
        elif (i==6):
            g7 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], signal, bkg)
            g7.SetName("RelIsoTight")

    # Create the Rate Vs Efficiency TGraphs
    p = plots.ComparisonManyPlot(g1, [g2, g3, g4, g5, g6, g7], saveFormats=[])
    algos = getAlgos()

    # Set individual styles
    for index, h in enumerate(p.histoMgr.getHistos()):
        hName = h.getName()
        legDict[hName] = algos[index] #styles.getCaloLegend(index)
        p.histoMgr.forHisto(hName, styles.getTauAlgoStyle(h.getName())) #styles.getCaloStyle(index))
        p.histoMgr.setHistoDrawStyle(h.getName(), "LX") # "X" = Do not draw error bars
        p.histoMgr.setHistoLegendStyle(h.getName(), "L") #LP

    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legDict)

    # Draw and save the plot
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary

    # Draw Error bands
    if opts.errorBands:
    #for i, g in enumerate([g0, g1, g2, g3, g4]):
        for i, g in enumerate([g1, g2, g3]):
            shapes, min, max = DrawErrorBand(g) 
            for shape in shapes:
                shape.SetFillColor( p.histoMgr.getHistos()[i].getRootHisto().GetFillColor())
                shape.SetFillStyle(3002)
                shape.Draw("f same")
            ROOT.gPad.RedrawAxis()

    histograms.addPileupText("PU=%s" % (bPU) )
    histograms.addText(0.55, 0.48, plots._legendLabels[signal], 18)

    # Save the plots in custom list of saveFormats
    aux.SavePlot(p, opts.saveDir, saveName, opts.saveFormats, True)
    return

def convert2RateVsEffTGraph(datasetsMgr, effHistoName, rateHistoName, signal, bkg):

    hEff  = datasetsMgr.getDataset(signal).getDatasetRootHisto(effHistoName).getHistogram()
    hRate = datasetsMgr.getDataset(bkg).getDatasetRootHisto(rateHistoName).getHistogram()
    
    # Sanity Checks
    if (hEff.GetXaxis().GetBinWidth(0) != hRate.GetXaxis().GetBinWidth(0)):
        Print("Efficiency histogram '%s' and rate histogram '%s' have different binning." % (effHistoName,rateHistoName), True)
        sys.exit()
    if (hEff.GetNbinsX() != hRate.GetNbinsX()):
        Print("Efficiency histogram '%s' and rate histogram %s have different number of bins." % (effHistoName,rateHistoName), True)
        sys.exit()
    
    # Lists for values                                                                                                                                         
    x     = []
    y     = []
    xerrl = []
    xerrh = []
    yerrl = []
    yerrh = []
    
    nBinsX = hEff.GetNbinsX()
    nBins  = nBinsX

    for i in range (0, nBinsX):
        # Get values
        xVal  = hEff.GetBinContent(i)
        xLow  = hEff.GetBinError(i)
        xHigh = xLow
        yVal  = hRate.GetBinContent(i)
        yLow  = hRate.GetBinError(i)
        yHigh = yLow            

        # Force error bars to not be above (belo) 1.0 (0.0)
        if 0:
            if abs(yVal + yHigh) > 1.0:
                yHigh = 1.0-yVal
            if yVal - yLow < 0.0:
                yLow = yVal
        
        # WARNING! Ugly trick so that zero points are not visible on canvas
        if 1:
            if xVal == 0.0:
                nBins=nBins-1
                continue

        # Store values
        x.append(xVal)
        xerrl.append(xLow)
        xerrh.append(xHigh)
        y.append(yVal)
        yerrl.append(yLow)
        yerrh.append(yHigh)
        
    # Create the TGraph with asymmetric errors
    tgraph = ROOT.TGraphAsymmErrors(nBins,
                                    array.array("d",x),
                                    array.array("d",y),
                                    array.array("d",xerrl),
                                    array.array("d",xerrh),
                                    array.array("d",yerrl),
                                    array.array("d",yerrh))
    
    tgraph.GetXaxis().SetLimits(0.0, 1.0)
    # Construct info table (debugging)
    table  = []
    align  = "{0:>6} {1:^10} {2:>10} {3:>10} {4:>10} {5:^3} {6:<10}"
    header = align.format("#", "xLow", "Efficiency", "xUp", "Rate", "+/-", "Error") #Purity = 1-EWK/Data
    hLine  = "="*70
    table.append("")
    table.append(hLine)
    table.append("{0:^70}".format(effHistoName))
    table.append(header)
    table.append(hLine)
    
    # For-loop: All values x-y and their errors
    for i, xV in enumerate(x, 0):
        row = align.format(i+1, "%.4f" % xerrl[i], "%.4f" %  x[i], "%.4f" %  xerrh[i], "%.5f" %  y[i], "+/-", "%.5f" %  yerrh[i])
        table.append(row)
    table.append(hLine)

    if 0:#printValues:
        for i, line in enumerate(table, 1):
            Print(line, False) #i==1)        
    return tgraph

def DrawErrorBand(graph):
    isErrorBand = graph.GetErrorXhigh(0) != -1 and graph.GetErrorXlow(0) != -1
    npoints     = graph.GetN()
 
    if not isErrorBand:
        graph.Draw("l same")
        return
 
    # Declare individual TGraph objects used in drawing error band
    central, min, max = ROOT.TGraph(), ROOT.TGraph(), ROOT.TGraph()
    shapes = []
    for i in range((npoints-1)*4):
        shapes.append(ROOT.TGraph())
 
    # Set ownership of TGraph objects
    ROOT.SetOwnership(central, False)
    ROOT.SetOwnership(    min, False)
    ROOT.SetOwnership(    max, False)
    for shape in shapes:
        ROOT.SetOwnership(shape, False)
 
    # Get data points from TGraphAsymmErrors
    x, y, xmin, xmax = [], [], [], []
    for i in range(npoints):
        tmpX, tmpY = ROOT.Double(0), ROOT.Double(0)
        graph.GetPoint(i, tmpX, tmpY)
        x.append(tmpX)
        y.append(tmpY)
        xmin.append(tmpX - graph.GetErrorXlow(i))
        xmax.append(tmpX + graph.GetErrorXhigh(i))

    # Fill central, min and max graphs
    for i in range(npoints):
        # central.SetPoint(i, x[i], y[i])
        min.SetPoint(i, xmin[i], y[i])
        max.SetPoint(i, xmax[i], y[i])
 
    # Fill shapes which will be shaded to create the error band
    for i in range(npoints-1):
        for version in range(4):
            shapes[i+(npoints-1)*version].SetPoint((version+0)%4, xmax[i],   y[i])
            shapes[i+(npoints-1)*version].SetPoint((version+1)%4, xmax[i+1], y[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+2)%4, xmin[i+1], y[i+1])
            shapes[i+(npoints-1)*version].SetPoint((version+3)%4, xmin[i],   y[i])
 
    # Set attributes to those of input graph
    min.SetLineStyle(graph.GetLineStyle())
    max.SetLineStyle(graph.GetLineStyle())
    return shapes, min, max 

def GetHistoKwargs(h, opts):
    _mvLeg1 = {"dx": -0.20, "dy": -0.45, "dh": -0.02}
    _mvLeg2 = {"dx": -0.10, "dy": -0.07, "dh": -0.00}
    _mvLeg3 = {"dx": -0.15, "dy": -0.45, "dh": -0.02}
    _mvLeg4 = {"dx": -0.52, "dy": -0.07, "dh": -0.05}
    #_mvLeg5 = {"dx": -0.52, "dy": -0.45, "dh": -0.02}
    _mvLeg5 = _mvLeg1

    logY    = True
    yMin    = 0.0
    if logY:
        yMin = 1
        yMaxF = 10
    else:
        yMaxF = 1.2
        
    _kwargs = {
        "xlabel"           : "Efficiency",
        "ylabel"           : "Rate (kHz)", #"Rate (kHz) / %.0f GeV",
        "addMCUncertainty" : False, 
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "opts"             : {"xmin": 0.0, "xmax": 0.8, "ymin": yMin, "ymax":1000, "ymaxfactor": yMaxF},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : logY,
        "moveLegend"       : _mvLeg1,
        "xtitlesize"       : 0.1,#xlabelSize,
        "ytitlesize"       : 0.1,#ylabelSize,
        "cutBoxY"           : {"cutValue": 50, "fillColor": 16, "box": False, "line": True, "cutGreaterThan"   : False}
        }

    if "RateVsEff_" in h:
        _kwargs["moveLegend"] = _mvLeg5
        _kwargs["opts"]       = {"xmin": 0.0, "xmax": 1.0, "ymin": yMin, "ymax":1000, "ymaxfactor": yMaxF}
        if "ditau" in h.lower():
            _kwargs["opts"]   = {"xmin": 0.0, "xmax": 0.6, "ymin": yMin, "ymax":1000, "ymaxfactor": yMaxF}
            #_kwargs["moveLegend"] = _mvLeg1

    if "Rate_" in h:
        _kwargs["xlabel"]     = "E_{T} (GeV)"
        _kwargs["ylabel"]     = "Rate (kHz)"# / %.0f GeV"
        #_kwargs["opts"]       = {"xmin": 0.0, "xmax": 300.0, "ymin": yMin, "ymax":1e5, "ymaxfactor": yMaxF}
        _kwargs["opts"]       = {"xmin": 0.0, "xmax": 100.0, "ymin": 1, "ymax":5e4, "ymaxfactor": yMaxF}
        _kwargs["moveLegend"] = _mvLeg2
        _kwargs["cutBoxY"]    = {"cutValue": 50, "fillColor": 16, "box": False, "line": True, "cutGreaterThan": False}

    if "Efficiency_" in h:
        units = "GeV"
        _kwargs["xlabel"]     = "E_{T} (%s)" % (units)
        _kwargs["ylabel"]     = "Efficiency / %0.0f " + units
        _kwargs["log"]        = False
        _kwargs["opts"]       = {"xmin": 0.0, "xmax": 100.0, "ymin": 0.0, "ymax": 1.0, "ymaxfactor": yMaxF}
        _kwargs["moveLegend"] = _mvLeg2 #_mvLeg3
        _kwargs["cutBoxY"]    = {"cutValue": 50, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
        _kwargs["cutBoxX"]    = {"cutValue": 10, "fillColor": 16, "box": True, "line": True, "cutGreaterThan": False}

    if "TurnOn" in h:
        _units = "GeV"
        _kwargs["xlabel"]     = "#tau_{h} E_{T}^{vis} (%s)" % (_units)
        _kwargs["ylabel"]     = "Efficiency / %0.0f " + _units
        _kwargs["log"]        = False
        _kwargs["rebinX"]     = 1 # do NOT change
        _kwargs["opts"]       = {"xmin": 0.0, "xmax": 200.0, "ymin": 0.0, "ymax": 1.15, "ymaxfactor": yMaxF}
        _kwargs["cutBoxY"]    = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
        _kwargs["moveLegend"] = _mvLeg4
        if "50" in h:
            _kwargs["cutBox"] = {"cutValue": 50.0, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
        if "25" in h:
            _kwargs["cutBox"] = {"cutValue": 25.0, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
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
    EBANDS       = False
    GRIDX        = False
    GRIDY        = False
    OPTMODE      = None
    PRECISION    = 3
    RATIO        = False
    SAVEDIR      = None
    SAVEFORMATS = [".C", ".png", ".pdf"]
    SEARCHMODE   = None
    URL          = False
    VERBOSE      = False

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

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

    parser.add_option("--errorBands", dest="errorBands", action="store_true", default=EBANDS, 
                      help="Enable error bands for ROC curves [default: %s]" % EBANDS)

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
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="ROC")
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
        raw_input("=== plotRateVsEff.py: Press any key to quit ROOT ...")

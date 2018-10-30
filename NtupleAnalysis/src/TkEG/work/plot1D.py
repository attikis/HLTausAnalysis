#!/usr/bin/env python
'''
DESCRIPTION:
Basic plotting script for making plots for CaloTk analyzer.


USAGE:
./plotL1TkTau.py -m <multicrab_directory> [opts]
./plotL1TkTau.py -m multicrab_CaloTkSkim_v92X_20180801T1203/ -i "SingleNeutrino_L1TPU140|TT_TuneCUETP8M2T4_14TeV_L1TnoPU" -n


LAST USED:
./plotL1TkTau.py -m multicrab_CaloTkSkim_v92X_20180801T1203/ -i "SingleNeutrino_L1TPU140|TT_TuneCUETP8M2T4_14TeV_L1TnoPU|TT_TuneCUETP8M2T4_14TeV_L1TPU140" -n

'''
#================================================================================================
# Imports
#================================================================================================
import os
import sys
from optparse import OptionParser
import getpass
import socket
import json
import copy

import HLTausAnalysis.NtupleAnalysis.tools.dataset as dataset
import HLTausAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HLTausAnalysis.NtupleAnalysis.tools.styles as styles
import HLTausAnalysis.NtupleAnalysis.tools.plots as plots
import HLTausAnalysis.NtupleAnalysis.tools.histograms as histograms
import HLTausAnalysis.NtupleAnalysis.tools.aux as aux
import HLTausAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

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
# Main
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

def GetLumi(datasetsMgr):
    Verbose("Determininig Integrated Luminosity")

    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            lumi += d.getLuminosity()
    Verbose("Luminosity = %s (pb)" % (lumi), True )
    return lumi

def GetDatasetsFromDir(opts):
    Verbose("Getting datasets")
    
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=None,
                                                        analysisName=opts.analysis)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=None,
                                                        analysisName=opts.analysis,
                                                        includeOnlyTasks=opts.includeOnlyTasks)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=None,
                                                        analysisName=opts.analysis,
                                                        excludeTasks=opts.excludeTasks)
    else:
        raise Exception("This should never be reached")
    return datasets
    
    
def PlotHisto(datasetsMgr, h):
    dsetsMgr = datasetsMgr.deepCopy()

    if "eff" in h.lower():
        dsetsMgr.remove("SingleNeutrino_L1TPU140", close=False)
        dsetsMgr.remove("SingleNeutrino_L1TPU200", close=False)
        opts.normalizeToOne = False
    elif "resolution" in h.lower():
        dsetsMgr.remove("SingleNeutrino_L1TPU140", close=False)
        dsetsMgr.remove("SingleNeutrino_L1TPU200", close=False)
    elif "rate" in h.lower():
        opts.normalizeToOne = False
        for d in dsetsMgr.getAllDatasetNames():
            if "SingleNeutrino" in d:
                continue
            else:
                dsetsMgr.remove(d, close=False)
    else:
        pass


    
    # Create the MC Plot with selected normalization ("normalizeToOne", "normalizeByCrossSection", "normalizeToLumi")
    kwargs = {}
    '''
    hList  = getHistoList(dsetsMgr, h)

    if opts.normalizeToOne:
        if 1:
            p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[], **kwargs)
            p.histoMgr.forEachHisto(lambda h: h.getRootHisto().Scale(1.0/h.getRootHisto().Integral()) )
        else:
            # p = plots.MCPlot(dsetsMgr, h, normalizeToOne=True, saveFormats=[], **kwargs)
            p = plots.PlotSameBase(dsetsMgr, h, normalizeToOne=True, saveFormats=[], **kwargs)
    else:
        if 1:
            p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[], **kwargs) #FIXME
        else:
            # p = plots.MCPlot(dsetsMgr, h, normalizeToLumi=opts.intLumi, saveFormats=[], **kwargs)
            p = plots.PlotSameBase(dsetsMgr, h, saveFormats=[], **kwargs)
    
    '''
    if opts.normalizeToOne:
        #p = plots.MCPlot(dsetsMgr, h, normalizeToOne=True, saveFormats=[], **kwargs)
        p = plots.PlotSameBase(dsetsMgr, h, normalizeToOne=True, saveFormats=[], **kwargs)
    else:
        #p = plots.MCPlot(dsetsMgr, h, normalizeToLumi=opts.intLumi, saveFormats=[], **kwargs)
        p = plots.PlotSameBase(dsetsMgr, h, saveFormats=[], **kwargs) #def __init__(self, datasetMgr, name, normalizeToOne=False, datasetRootHistoArgs={}, **kwargs):

    

    # Set default styles (Called by default in MCPlot)
    p._setLegendStyles()
    p._setLegendLabels()
    p._setPlotStyles()

    # Customise legend
    # p.histoMgr.setHistoLegendStyleAll("L")
    for d in dsetsMgr.getAllDatasetNames():
        if "SingleNeutrino" in d:
            p.histoMgr.setHistoLegendStyle(d, "F")
        else:
            p.histoMgr.setHistoLegendStyle(d, "L")

    p.setLegend(histograms.createLegend(0.68, 0.85-0.05*len(dsetsMgr.getAllDatasetNames()), 0.77, 0.92))
    #p.histoMgr.setHistoLegendStyle("histo_" + dataset, "LP")

    # Draw a customised plot
    plots.drawPlot(p, h, **GetHistoKwargs(h, opts) )
    if ("photons_egs_matching" in h.lower()):
        p.getFrame().GetXaxis().SetLabelSize(14)
        p.getFrame().GetXaxis().LabelsOption("u")
        

    # Remove legend?
    if 0:
        p.removeLegend()

    # Save in all formats chosen by user
    aux.SavePlot(p, opts.saveDir, h, opts.saveFormats, opts.url)
    return


def getHistoList(datasetsMgr, histoName):
    hList = []
    # For-loop: All dataset names
    for d in datasetsMgr.getAllDatasetNames():
        h = datasetsMgr.getDataset(d).getDatasetRootHisto(histoName)
        hList.append(h)
    return hList


def GetHistoKwargs(h, opts):
    
    hName   = h.lower()
    #_xLabel = ""
    _yLabel = "Arbitrary Units / %.2f "
    _rebinX = 1
    _rebinY = None
    _units  = ""
    _format = "%0.0f " + _units
    _cutBox = {"cutValue": 10.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
    _cutBoxY= {"cutValue": 50.0, "fillColor": 16, "box": False, "line": False, "cutGreaterThan"   : False}
    _leg   = {"dx": -0.5, "dy": -0.3, "dh": -0.4}
    _ratio = True
    _log   = False
    _yMin  = 0.0
    _yMaxF = 1.2
    _xMin  = None
    _xMax  = None
    _yNorm = "Arbitrary Units"
    
    if _log:
        if _yMin == 0.0:
            _yMin = 1e-3
        _yMaxF = 10

    # Efficiency & Rate plots
    if "eff" in hName:
        _units  = "GeV"
        _format = "%0.0f " + _units
        _xLabel = "E_{T} (%s)" % (_units)
        _cutBox = {"cutValue": 20.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        _rebinX = 1
        _yLabel = "Efficiency / %.0f " + _units
    if "counters" in hName:
        _units  = ""
        _format = "%0.0f " + _units
        _xLabel = "counters"
        _rebinX = 1
        #_xMax   = +10.0
        _yLabel = _yNorm + " / " + _format
        _log    = False
    if "rate" in hName:
        _units  = "GeV"
        _format = "%0.0f " + _units
        _xLabel = "E_{T} (%s)" % (_units)
        _cutBox = {"cutValue": 50.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        _cutBoxY= {"cutValue": 50, "fillColor": 16, "box": False, "line": True, "cutGreaterThan"   : False}
        _rebinX = 1
        _yLabel = "Rate (kHz) / %.0f " + _units
        _log    = True
        _yMin   = 1e+0
        _xMax   = 250.0
        
    if "etresolution" in hName:
        _units  = ""
        _format = "%0.2f " + _units
        #_xLabel = "(E_{T}^{calo} - p_{T}^{vis}) / p_{T}^{vis}"
        _xLabel = "#deltaE_{T} / p_{T}^{vis}"
        _rebinX = 1
        _xMin   = -4.0 #-5.5
        _xMax   = +4.0 #+5.5

    if "etaresolution" in hName:
        # _xLabel = "(#eta^{calo} - #eta^{vis}) / #eta^{vis}"
        _xLabel = "#delta#eta / #eta^{vis}"
        _xMin   = -4.0 #-5.5
        _xMax   = +4.0 #+5.5
        _cutBox = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "phiresolution" in hName:
        _units  = ""
        _format = "%0.2f " + _units
        #_xLabel = "#phi^{calo} - #phi^{vis} / #phi^{vis}"
        _xLabel = "#delta#phi / #phi^{vis}"
        _rebinX = 1
        _xMin   = -2.0 #-5.0
        _xMax   = +2.0 #+5.0
        _cutBox = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    _kwargs = {
        #"xlabel"           : _xLabel,
        "ylabel"           : _yLabel,
        "rebinX"           : _rebinX,
        "rebinY"           : _rebinY,
        "ratioYlabel"      : "Ratio",
        #"ratio"            : _ratio,
        "stackMCHistograms": False,
        "ratioInvert"      : False,
        "addMCUncertainty" : True,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "opts"             : {"ymin": _yMin, "ymaxfactor": _yMaxF},
        "log"              : _log,
        "cutBox"           : _cutBox,
        "cutBoxY"          : _cutBoxY,
        "createLegend"     : None #_leg,
        }

    kwargs = copy.deepcopy(_kwargs)

    if "nstubs" in hName:
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 10, "ymin": _yMin, "ymaxfactor": _yMaxF}
                

    if "multiplicity" in h.lower():
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 20, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "multiplicitypercluster" in h.lower():
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 6.5, "ymin": _yMin, "ymax" : 1.2, "ymaxfactor": _yMaxF}

    if "mcmatch_chargeddaugh_n" in h.lower():
        kwargs["xlabel"] = "Number of daughters^{+} of matched gen-#tau"
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 10.5, "ymin": _yMin, "ymax" : 1.2,"ymaxfactor": _yMaxF}

    if "mcmatch_neutraldaugh_n" in h.lower():
        kwargs["xlabel"] = "Number of daughters^{0} of matched gen-#tau"
        kwargs["opts"]   = {"xmin": -0.5, "xmax": 10.5, "ymin": _yMin,"ymax" : 1.2,"ymaxfactor": _yMaxF}

    if "chargeddaugh_totalmass" in h.lower():
        kwargs["log"]  = True
        kwargs["opts"] = {"xmin": 0.0, "xmax": 2.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "neutraldaugh_totalmass" in h.lower():
        kwargs["log"]  = True
        kwargs["opts"] = {"xmin": 0.0, "xmax": 2.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.77, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "chargeddaugh_pt" in h.lower():
        #kwargs["log"]  = True                                                                                                                                           
        kwargs["opts"] = {"xmin": 0.0, "xmax": 120.0, "ymin": _yMin, "ymax": 0.065, "ymaxfactor": _yMaxF}

    if "neutraldaugh_et" in h.lower():
        #kwargs["log"]  = True                                                                                                                                           
        kwargs["opts"] = {"xmin": 0.0, "xmax": 120.0, "ymin": _yMin, "ymax": 0.065, "ymaxfactor": _yMaxF}

    if "photons_dr" in h.lower():
        kwargs["opts"] = {"xmin": 0.0, "xmax": 0.5, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "chi2" in h.lower():
        kwargs["rebinX"] = 2

    if "gentaushadronic_n" in h.lower():
        kwargs["xlabel"] = "No of genuine hadronic taus in event"

    if "trkClusters_M" == h:
        kwargs["log"]  = True
        kwargs["opts"]   = {"xmin": 0, "xmax": 2.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "trkClusters_M_beforeCut" == h:
        kwargs["log"]  = True
        kwargs["opts"]   = {"xmin": 0, "xmax": 2.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "trkClusters_Pt" == h:
        kwargs["opts"]   = {"xmin": 0, "xmax": 50.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["log"]  = True


    if "clustegs_eta" in h.lower():
        kwargs["opts"]   = {"xmin": -2.0, "xmax": 2.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
         
    if "sigcone_deltar" in h.lower():
        kwargs["log"]  = True
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 0.16, "ymin": 0.0001, "ymaxfactor": _yMaxF}
        kwargs["xlabel"] = "R_{max}^{sig}"

    if "tkeg_isotracks_invmass" in h.lower():
        kwargs["log"]  = True
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 0.5, "ymin": 0.0001, "ymaxfactor": _yMaxF}
        
    if "isotracks_n" in h.lower():
        kwargs["opts"]   = {"xmin": 0, "xmax": 10.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        
    if "tkeg_reliso" in h.lower() or "badetresolcand_reliso" in h.lower():
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.2, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.20, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["log"]  = True
        
    if "tkeg_vtxiso" in h.lower() or "badetresolcand_vtxiso" in h.lower():
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 4.0, "ymin": 0.0001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["log"]  = True

    if "EGClusters_M" == h:
        kwargs["opts"]   = {"xmin": 0, "xmax": 2.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "EGClusters_Et" == h:
        kwargs["opts"]   = {"xmin": 0, "xmax": 80.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        
    if "TkEG_ET" == h:
        kwargs["opts"]   = {"xmin": 0, "xmax": 120.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["log"]  = True

    if "chf" in hName or "nhf" in hName:
        kwargs["opts"]   = {"xmin": 0, "xmax": 1.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["log"]  = True

    if "tkeg_invmass" in hName:
        kwargs["opts"]   = {"xmin": 0, "xmax": 2.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["log"]  = True
        kwargs["cutBox"] = {"cutValue": 1.77, "fillColor": 16, "box": False, "line": True, "greaterThan": True}


    if "leadtrks_multiplicity" in h.lower():
        kwargs["opts"]   = {"xmin": 0, "xmax": 10.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "leadtrk_eg_drmin" in h.lower():
        kwargs["log"]  = True
        kwargs["opts"]   = {"xmin": 0, "xmax": 5.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.3, "fillColor": 16, "box": False, "line": True, "greaterThan": True}


    if "leadtrks_pt" in h.lower():
        _yLabel = "Arbitrary Units / %.0f "
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0, "xmax": 50.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["log"]    = True                                                                                                                                         
    
    if "clusttrks_pt" in h.lower():
        _yLabel = "Arbitrary Units / %.0f "
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 50.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        #kwargs["log"]    = True                                                                                                                                         

    if "EGs_Et" == h:
        kwargs["log"]  = True
        kwargs["opts"] = {"xmin": 0.0, "xmax": 100.0, "ymin": 0.001, "ymax": 0.2, "ymaxfactor": _yMaxF}

    if "EGs_MCmatched_Et" == h:
        kwargs["log"]  = True
        kwargs["opts"] = {"xmin": 0.0, "xmax": 100.0, "ymin": 0.001, "ymax": 0.2, "ymaxfactor": _yMaxF}

    if "mcmatch_dr" in h.lower():
        kwargs["log"]  = True
        kwargs["opts"] = {"xmin": 0.0, "xmax": 5.0, "ymin": 0.001, "ymaxfactor": _yMaxF}

    if "dz0" in h.lower():
        kwargs["log"]  = True
        kwargs["opts"] = {"xmin": 0.0, "xmax": 10.0, "ymin": 0.001, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.8, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "etresolution" in hName or "ptresolution" in hName:
        kwargs["opts"] = {"xmin": -1.0, "xmax": 1.0, "ymin": 0.0, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "etaresolution" in hName:
        kwargs["opts"] = {"xmin": -0.3, "xmax": 0.3, "ymin": 0.0, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "phiresolution" in hName:
        kwargs["opts"] = {"xmin": -0.3, "xmax": 0.3, "ymin": 0.0, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        


    return kwargs

def GetBinwidthDecimals(binWidth):                                                                                                                                   
    dec =  " %0.0f"
    if binWidth < 1:
        dec = " %0.1f"
    if binWidth < 0.1:
        dec = " %0.2f"
    if binWidth < 0.01:
        dec =  " %0.3f"
    if binWidth < 0.001:
        dec =  " %0.4f"
    if binWidth < 0.0001:
        dec =  " %0.5f"
    if binWidth < 0.00001:
        dec =  " %0.6f"
    if binWidth < 0.000001:
        dec =  " %0.7f"
    return dec

def ReorderDatasets(datasets):
    newOrder =  datasets
    
    for i, d in enumerate(datasets, 0):
        if "PU200" in d:
            newOrder.remove(d)
            newOrder.insert(0, d)
            #newOrder.insert(0, newOrder.pop(i))
    for j, d in enumerate(datasets, 0):
        if "PU140" in d:
            newOrder.remove(d)
            newOrder.insert(0, d)
    for k, d in enumerate(datasets, 0):
        if "noPU" in d:
            newOrder.remove(d)
            newOrder.insert(0, d)
    
    mb140 = "SingleNeutrino_L1TPU140"
    mb200 = "SingleNeutrino_L1TPU200"
    if mb140 in datasets:
        newOrder.remove(mb140)
        newOrder.insert(len(newOrder), mb140)
    if mb200 in datasets:
        newOrder.remove(mb200)
        newOrder.insert(len(newOrder), mb200)
    return newOrder


def main(opts):
    
    # Set the ROOTeError verbosity
    ROOT.gErrorIgnoreLevel=3000 # kUnset=-1, kPrint=0, kInfo=1000, kWarning=2000, kError=3000, kBreak=4000

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    style.setOptStat(False)


    # Obtain dsetMgrCreator and register it to module selector
    dsetMgrCreator = dataset.readFromMulticrabCfg(directory=opts.mcrab)

    # Setup & configure the dataset manager
    datasetsMgr = GetDatasetsFromDir(opts)
    datasetsMgr.updateNAllEventsToPUWeighted()

    if opts.verbose:
        datasetsMgr.PrintCrossSections()
        datasetsMgr.PrintInfo()

    # Get all the histograms and their paths
    #hList = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(opts.folder)
    #print hList
    #sys.exit()

    # Setup & configure the dataset manager (no collision data => not needed)
    if 0:
        datasetsMgr.loadLuminosities()
        datasetsMgr.updateNAllEventsToPUWeighted()

    # Print information
    if opts.verbose:
        datasetsMgr.PrintCrossSections()
        # datasetsMgr.PrintLuminosities()

    # Print dataset information (before merge)        
    datasetsMgr.PrintInfo()
        
    # Merge histograms (see NtupleAnalysis/python/tools/plots.py)    
    plots.mergeRenameReorderForDataMC(datasetsMgr)

    # Get Luminosity
    if 0:
        intLumi = datasetsMgr.getDataset("Data").getLuminosity()

    # Apply new dataset order?
    newOrder = ReorderDatasets(datasetsMgr.getAllDatasetNames())
    datasetsMgr.selectAndReorder(newOrder)

    # Print dataset information (after merge)
    if 0:
        datasetsMgr.PrintInfo() #Requires python 2.7.6 or 2.6.6

    # Plot Histograms
    histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(opts.folder)
    histoPaths = [os.path.join(opts.folder, h) for h in histoList]
    histoType  = type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram())
    plotCount  = 0
    skipList   = []

    # For-loop: All histos in opts.folder
    for i, h in enumerate(histoPaths, 1):
        
        # Obsolete quantity
        #if h in skipList:
        if "counter" in h.lower():
            continue

        histoType  = str(type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram()))
        if "TH1" not in histoType:
            continue
        
        aux.PrintFlushed(h, plotCount==0)
        plotCount += 1
        
        if "negs" in h.lower() or "ntks" in h.lower():
            ROOT.gStyle.SetNdivisions(8, "X")
        
        PlotHisto(datasetsMgr, h)

    print
    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return

#================================================================================================
# Main
#================================================================================================
if __name__ == "__main__":

    # Default Settings 
    ANALYSIS    = "HLTausAnalysis"
    BATCHMODE   = True
    DATAERA     = "TDR2019"
    FOLDER      = ""
    GRIDX       = False
    GRIDY       = False   
    INTLUMI     = 1.0
    NORMTOONE   = True
    SAVEDIR     = None
    SAVEFORMATS =[".pdf"] # [".C", ".png", ".pdf"]
    VERBOSE     = False

    parser = OptionParser(usage="Usage: %prog [options]" , add_help_option=False,conflict_handler="resolve")

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation  NOT generates a window) [default: %s]" % BATCHMODE)

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", default=None,
                      help="Path to the multicrab directory for input")

    parser.add_option("-n", "--normalizeToOne", dest="normalizeToOne", action="store_true", default=NORMTOONE,
                      help="Normalise all histograms to unit area? [default: %s]" % (NORMTOONE) )

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX,
                      help="Enable x-axis grid? [default: %s]" % (GRIDX) )
    
    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY,
                      help="Enable y-axis grid? [default: %s]" % (GRIDY) )

    parser.add_option("--url", dest="url", action="store_true", default=False,
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % False)

    parser.add_option("--formats", dest="formats", default = None,
                      help="Formats in which all plots will be saved in. Provide as list of comma-separated (NO SPACE!) formats. [default: None]")

    parser.add_option("--analysis", dest="analysis", type="string", default=ANALYSIS,
                      help="Override default analysis [default: %s]" % ANALYSIS)

    parser.add_option("--dataEra", dest="dataEra", default = DATAERA,
                      help="Formats in which all plots will be saved in. Provide as list of comma-separated (NO SPACE!) formats. [default: %s]" % (DATAERA))

    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if opts.mcrab == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        print __doc__
        sys.exit(0)
    
    # Determine path for saving plots
    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/TkEG/", postfix="TH1D")
    else:
        print "opts.saveDir = ", opts.saveDir

    # Overwrite default save formats?
    if opts.formats != None:
        opts.saveFormats = opts.formats.split(",")
    else:
        opts.saveFormats = SAVEFORMATS

    # Inform user of compatibility issues
    pyV1  =  sys.version_info[0]
    pyV2  =  sys.version_info[1]
    pyV3  =  sys.version_info[2]
    pyVer = "%d.%d.%d" % (pyV1, pyV2, pyV3)
    if pyV2 < 7 or pyV3 < 6:
        Print("Recommended %sPython 2.7.6%s or later (using %sPython %s). EXIT!" % (hs, ns, es, pyVer + ns), True)
        #sys.exit()
    else:
        Print("Recommended %sPython 2.7.6%s or later (using %sPython %s)" % (hs, ns, ss, pyVer + ns), True)

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotL1TkTau.py: Press any key to quit ROOT ...")

#!/usr/bin/env python
'''
DESCRIPTION:
Basic plotting script for making plots for CaloTk analyzer.


USAGE:
./plotCounters.py -m <multicrab_directory> [opts]
./plotCounters.py -m multicrab_CaloTkSkim_v92X_20180801T1203/ -i "SingleNeutrino_L1TPU140|TT_TuneCUETP8M2T4_14TeV_L1TnoPU" -n
./plotCounters.py -m multicrab_CaloTkSkim_v92X_20180801T1203/ -i "SingleNeutrino_L1TPU140|TT_TuneCUETP8M2T4_14TeV_L1TnoPU|TT_TuneCUETP8M2T4_14TeV_L1TPU140" -n
./plotCounters.py -m multicrab_CaloTk_v92X_IsoConeRMax0p4_VtxIso1p0_08h09m41s_23Aug2018 -e "TT|Glu|SingleTau|Higgs1000|Higgs500" -n
./plotCounters.py -m multicrab_CaloTk_v92X_IsoConeRMax0p4_VtxIso1p0_08h09m41s_23Aug2018 -e "TT|SingleTau|Higgs" -n 
./plotCounters.py -m multicrab_CaloTk_v92X_IsoConeRMax0p4_VtxIso1p0_08h09m41s_23Aug2018 -e "TT|Glu|Higgs" -n
./plotCounters.py -m multicrab_CaloTk_v92X_IsoConeRMax0p3_VtxIso0p5_RelIso0p2_14h29m15s_23Aug2018 -e "TT|SingleTau|Higgs|SingleE" -n


LAST USED:
./plotCounters.py -n -e "SingleE|Charged|TT|GluGlu" -m 


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
    Verbose("Plotting histogram %s" % (h), True)

    dsetsMgr = datasetsMgr.deepCopy()

    # Create the plot with selected normalization
    kwargs = {}
    hList  = getHistoList(dsetsMgr, h)

    if opts.normalizeToOne:
        p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[], **kwargs)
        
        # For-loop: All histograms in manager
        normValue = 1.0
        for hist in p.histoMgr.getHistos():
            # print "hist.getRootHisto().GetName() = ", hist.getRootHisto().GetName()
            # print "histo.getRootHistos().GetBinContent(1) = ", hist.getRootHisto().GetBinContent(1)
            normValue = hist.getRootHisto().GetBinContent(1)
            hist.getRootHisto().Scale(1.0/normValue)
            
        # Normalise wrt bin #1
        #p.histoMgr.forEachHisto(lambda h: h.getRootHisto().Scale(1.0/normValue) )
    else:
        p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[], **kwargs)
        #p = plots.PlotSameBase(dsetsMgr, h, saveFormats=[], **kwargs)
            
    
    # Set default styles (Called by default in MCPlot)
    p._setLegendStyles()
    p._setLegendLabels()
    p._setPlotStyles()

    # Customise legend
    for d in dsetsMgr.getAllDatasetNames():
        if "SingleNeutrino" in d:
            p.histoMgr.setHistoLegendStyle(d, "F")
        else:
            p.histoMgr.setHistoLegendStyle(d, "L")
            p.histoMgr.setHistoDrawStyle(d, "HIST9")

    # Create legend
    if 0:
        p.setLegend(histograms.createLegend(0.18, 0.86-0.04*len(dsetsMgr.getAllDatasetNames()), 0.42, 0.92))
    else:
        p.setLegend(histograms.createLegend(0.58, 0.86-0.04*len(dsetsMgr.getAllDatasetNames()), 0.92, 0.92))

    # Draw a customised plot
    kwargs = GetHistoKwargs(h, opts) 
    plots.drawPlot(p, h, **kwargs)

    # Remove legend?
    if kwargs["removeLegend"]:
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
    _xLabel = ""
    if opts.normalizeToOne:
        _yNorm = "Arbitrary Units"
    else:
        _yNorm = "Events"
    _yLabel = _yNorm + " / %.2f "
    _rebinX = 1
    _rebinY = None
    _units  = ""
    _format = "%0.0f " + _units
    _cutBox = {"cutValue": 10.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
    _cutBoxY= {"cutValue": 1.00, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
    _leg    = {"dx": -0.4, "dy": -0.3, "dh": -0.4}
    _leg2   = {"dx": -0.75, "dy": -0.3, "dh": -0.4}
    _ratio  = True
    _log    = False
    _yMin   = 0.0
    _yMax   = None
    _yMaxF  = 1.15
    _xMin   = None
    _xMax   = None
    _rmLeg  = False
    _mvLeg  = {"dx": -0.50, "dy": -0.50, "dh": -0.0}

    # Default (tdrstyle.py)
    ROOT.gStyle.SetLabelSize(27, "XYZ")

    if "counters" in hName:
        _units  = ""
        _format = "%0.0f " + _units
        _xLabel = "" #"counters"
        _xMin   = 1 #don't show Loose/Tight WPs
        _xMax   = 8 #don't show Loose/Tight WPs
        _yLabel = _yNorm + " / " + _format
        #_log    = True
        #_rmLeg  = True
        #_ratio  = False
        _cutBox = {"cutValue": 0.15, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        ROOT.gStyle.SetLabelSize(13, "X") #"XY"
                
    if _log:
        if _yMin == 0.0:
            if opts.normalizeToOne:
                _yMin = 0.9e-4
            else:
                _yMin = 1e0
        _yMaxF = 5

    if _yMax == None:
        _opts = {"ymin": _yMin, "ymaxfactor": _yMaxF}
    else:
        _opts = {"ymin": _yMin, "ymax": _yMax}
    if _xMax != None:
        _opts["xmax"] = _xMax
    if _xMin != None:
        _opts["xmin"] = _xMin

    _opts2 = {"ymin": 0.0, "ymax": 2.3}

    _kwargs = {
        "xlabel"           : _xLabel,
        "ylabel"           : _yLabel,
        "rebinX"           : _rebinX,
        "rebinY"           : _rebinY,
        "ratioYlabel"      : "Ratio ", #"1/Ratio "
        "ratio"            : _ratio, # only plots.ComparisonManyPlot(). Eitherwise comment out
        "stackMCHistograms": False,
        "ratioInvert"      : False, #True,
        "addMCUncertainty" : True,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "opts"             : _opts,
        "opts2"            : _opts2,
        "log"              : _log,
        "cutBox"           : _cutBox,
        "cutBoxY"          : _cutBoxY,
        #"createLegend"     : None,
        "moveLegend"       : _mvLeg,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        "removeLegend"     : _rmLeg,
        }
    return _kwargs

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

    # Use tau guns for ratio reference => put very first
    if "SingleTau_L1TPU200" in datasets:
        newOrder.remove("SingleTau_L1TPU200")
        newOrder.insert(0,"SingleTau_L1TPU200")
    if "SingleTau_L1TPU140" in datasets:
        newOrder.remove("SingleTau_L1TPU140")
        newOrder.insert(0,"SingleTau_L1TPU140")
    
    mb140 = "SingleNeutrino_L1TPU140"
    mb200 = "SingleNeutrino_L1TPU200"
    if mb140 in datasets:
        newOrder.remove(mb140)
        newOrder.insert(len(newOrder), mb140)
    if mb200 in datasets:
        newOrder.remove(mb200)
        newOrder.insert(len(newOrder), mb200)

#    if mb140 in datasets and mb200 in datasets:
#        newOrder.remove(mb140)
#        newOrder.insert(-1, mb140)
#        newOrder.remove(mb200)
#        newOrder.insert(-1, mb200)

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
    histoPaths = [os.path.join(opts.folder, h) for h in histoList if "Counters" in h]
    histoType  = type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram())
    done       = []
    plotCount  = 0

    # For-loop: All histos in opts.folder    
    for i, h in enumerate(histoPaths, 1):
        histoType  = str(type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram()))

        if "TH1" not in histoType:
            continue

        if h in done:
            continue

        aux.PrintFlushed(h, plotCount==0)
        plotCount += 1
        PlotHisto(datasetsMgr, h)
        done.append(h)

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
    NORMTOONE   = False
    SAVEDIR     = None
    SAVEFORMATS = [".C", ".png", ".pdf"]
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
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="Counters")
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
        raw_input("=== plotCounters.py: Press any key to quit ROOT ...")

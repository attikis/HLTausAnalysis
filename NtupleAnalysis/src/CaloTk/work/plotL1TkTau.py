#!/usr/bin/env python
'''
DESCRIPTION:
Basic plotting script for making plots for CaloTk analyzer.


USAGE:
./plotL1TkTau.py -m <multicrab_directory> [opts]


LAST USED:
./plotL1TkTau.py -m multicrab_CaloTkSkim_v92X_20180801T1203/ -i "SingleNeutrino_L1TPU140|TT_TuneCUETP8M2T4_14TeV_L1TnoPU" 

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

    # Create the MC Plot with selected normalization ("normalizeToOne", "normalizeByCrossSection", "normalizeToLumi")
    p = plots.MCPlot(datasetsMgr, h, saveFormats=[], normalizeToOne=opts.normalizeToOne)
    # p = plots.ComparisonManyPlot(FakeB_inverted, compareHistos, saveFormats=[])       

    '''
    # For-loop: All histos
    for index, h in enumerate(p.histoMgr.getHistos()):
        if index == 0:
            p.histoMgr.setHistoDrawStyle(h.getName(), "HIST")
            p.histoMgr.setHistoLegendStyle(h.getName(), "L")
            #    continue
        else:
            p.histoMgr.setHistoDrawStyle(h.getName(), "AP")
            p.histoMgr.setHistoLegendStyle(h.getName(), "AP")

    # Set default dataset style to all histos
    for index, h in enumerate(p.histoMgr.getHistos()):
        plots._plotStyles[p.histoMgr.getHistos()[index].getDataset().getName()].apply(p.histoMgr.getHistos()[index].getRootHisto())
       #plots._plotStyles[dataset.getName()].apply(p)
    '''

    # Additional text
    #histograms.addText(0.18, 0.92, "Phase-2 Simulation")
    
    # Draw a customised plot
    kwargs_  = GetHistoKwargs(h, opts)
    plots.drawPlot(p, h, **kwargs_)

    # Remove legend?
    if 0:
        p.removeLegend()

    # Save in all formats chosen by user
    aux.SavePlot(p, opts.saveDir, h, opts.saveFormats, opts.url)
    return


def GetHistoKwargs(h, opts):
    _moveLegend = {"dx": -0.1, "dy": 0.0, "dh": -0.15}
    logY    = False
    _yLabel = "Arbitrary Units / %.2f "
    yMin    = 0.0
    if logY:
        yMaxF = 10
    else:
        yMaxF = 1.2

    _kwargs = {
        #"xlabel"           : "x",
        "ylabel"           : _yLabel,
        "rebinX"           : 1,
        "rebinY"           : None,
        "ratioYlabel"      : "Ratio",
        "ratio"            : False, #opts.ratio,
        "stackMCHistograms": False, #opts.nostack,
        "ratioInvert"      : False,
        "addMCUncertainty" : False,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "opts"             : {"ymin": yMin, "ymaxfactor": yMaxF},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : logY,
        "moveLegend"       : _moveLegend,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        }

    kwargs = copy.deepcopy(_kwargs)
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
    newOrder = ["SingleNeutrino_L1TPU140", "TT_TuneCUETP8M2T4_14TeV_L1TnoPU"]
    if 0:
        datasetsMgr.selectAndReorder(newOrder)

    # Print dataset information (after merge)
    if 0:
        datasetsMgr.PrintInfo() #Requires python 2.7.6 or 2.6.6

    # Plot Histograms
    histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(opts.folder)
    histoPaths = [os.path.join(opts.folder, h) for h in histoList]
    histoType  = type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram())

    # For-loop: All histos in opts.folder
    for i, h in enumerate(histoPaths, 1):
        histoType  = str(type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram()))
        if "TH1" not in histoType:
            continue
        PlotHisto(datasetsMgr, h)

    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return


#================================================================================================
# Main
#================================================================================================
if __name__ == "__main__":

    # Default Settings 
    BATCHMODE   = True
    VERBOSE     = False
    INTLUMI     = 1.0
    SAVEDIR     = None
    SAVEFORMATS = [".png"] #[".C", ".png", ".pdf"]
    GRIDX       = False
    GRIDY       = False    
    ANALYSIS    = "HLTausAnalysis"
    DATAERA     = "ID2017" #"TDR2019"
    FOLDER      = ""
    NORMTOONE   = False

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
    print "opts.mcrab = ", opts.mcrab
    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="L1TkTau")
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

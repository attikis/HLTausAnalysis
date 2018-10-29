#!/usr/bin/env python
'''
DESCRIPTION:
Script that plots Data/MC for all histograms under a given folder (passsed as option to the script)
Good for sanity checks for key points in the cut-flow


USAGE:
./plotTH2.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotTH2.py -m <mcrab> --folder ""
./plotTH2.py -m multicrab_TkTau_v92X_SeedPt5_SigRMax0p15_IsoRMax0p35_VtxIso0p5_RelIso0p5_22h48m02s_09Oct2018 --logZ -i "TT_TuneCUETP8M2T4_14TeV_L1TnoPU" --normalizeToLumi -intLumi 1000
./plotTH2.py -m multicrab_TkTau_v92X_SeedPt5_SigRMax0p15_IsoRMax0p35_VtxIso0p5_RelIso0p5_22h48m02s_09Oct2018 --logZ -i "TT_TuneCUETP8M2T4_14TeV_L1TnoPU" --normalizeByCrossSection
./plotTH2.py -m multicrab_TkTau_v92X_SeedPt5_SigRMax0p15_IsoRMax0p35_VtxIso0p5_RelIso0p5_22h48m02s_09Oct2018 --logZ -i "TT_TuneCUETP8M2T4_14TeV_L1TnoPU" --normalizeToOne


LAST USED:
/plotTH2.py -e "SingleE" --logZ --normalizeToOne -m 

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
import HLTausAnalysis.NtupleAnalysis.tools.aux as aux
import HLTausAnalysis.NtupleAnalysis.tools.counter as counter
import HLTausAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HLTausAnalysis.NtupleAnalysis.tools.styles as styles
import HLTausAnalysis.NtupleAnalysis.tools.plots as plots
import HLTausAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HLTausAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HLTausAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

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

def GetLumi(datasetsMgr):
    Verbose("Determininig Integrated Luminosity")
    
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            lumi += d.getLuminosity()
    Verbose("Luminosity = %s (pb)" % (lumi), True)
    return lumi

def GetDatasetsFromDir(opts):
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode, 
                                                        analysisName=opts.analysisName,
                                                        optimizationMode=opts.optMode)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        includeOnlyTasks=opts.includeOnlyTasks,
                                                        optimizationMode=opts.optMode)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        excludeTasks=opts.excludeTasks,
                                                        optimizationMode=opts.optMode)
    else:
        raise Exception("This should never be reached")
    return datasets
    

def main(opts):

    optModes = [""]

    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    style.setLogX(opts.logX)
    style.setLogY(opts.logY)
    style.setLogZ(opts.logZ)
    #style.setWide(False, 0.15)

    # For-loop: All opt Mode
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        # datasetsMgr.loadLuminosities() # from lumi.json

        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()
            datasetsMgr.PrintInfo()

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 

        # Get Luminosity
        if opts.intLumi < 0:
            intLumi = datasetsMgr.getDataset("Data").getLuminosity()

        # Custom Filtering of datasets 
        datasetsToRemove = []
        for i, d in enumerate(datasetsToRemove, 0):
            msg = "Removing dataset %s" % d
            Print(ShellStyles.WarningLabel() + msg + ShellStyles.NormalStyle(), i==0)
            datasetsMgr.remove(filter(lambda name: d in name, datasetsMgr.getAllDatasetNames()))

        if opts.verbose:
            datasetsMgr.PrintInfo()

        # Re-order datasets (different for inverted than default=baseline)
        newOrder = []
        for d in datasetsMgr.getMCDatasets():
            newOrder.append(d.getName())
            
        # Apply new dataset order!
        datasetsMgr.selectAndReorder(newOrder)
        
        # Print dataset information
        datasetsMgr.PrintInfo()        

        # Plot Histograms
        histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(opts.folder)
        histoPaths = [os.path.join(opts.folder, h) for h in histoList]
        
        # For-loop: All histograms in chosen folder
        counter = 0
        for i, h in enumerate(histoPaths, 1):
            histoType  = str(type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram()))
            if "TH2" not in histoType:
                continue
            if "DiTau" in h:
                continue
            if "L1TkIsoTau" in h:
                continue

            # For-loop: All datasets
            for d in datasetsMgr.getAllDatasetNames():
                counter += 1
                #if d != "SingleNeutrino_L1TPU140":
                #    continue
                #if "reliso" not in h.lower():
                #    continue
                Plot2dHistograms(datasetsMgr, d, h, counter)
    print
    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return

def GetHistoKwargs(h, opts):
    _moveLegend = {"dx": -0.1, "dy": 0.0, "dh": -0.15}
    logY    = False
    yMin    = 0.0
    if logY:
        yMin = 0.001
        yMaxF = 10
    else:
        yMaxF = 1.0
    
    # z-axis settings
    zMin   =  0
    zMax   = None
    zLabel = "z-axis"
    if opts.normalizeToLumi:
        zLabel  = "Events"
        zMin    = 1e0
    elif opts.normalizeByCrossSection:
        zLabel  = "#sigma (pb)"
    elif opts.normalizeToOne:
        zLabel  = "Arbitrary Units"
    else:
        zLabel = "Unknown"

   # Cut lines/boxes                                                                                                                                                                                                              
    _cutBoxX = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
    _cutBoxY = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True,
                "mainCanvas": True, "ratioCanvas": False} # box = True not working       

    _kwargs = {
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": opts.normalizeToLumi,
        "addCmsText"       : True,
        "cmsExtraText"     : " Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "log"              : logY,
        "moveLegend"       : _moveLegend,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        "zmin"             : zMin, 
        "zmax"             : zMax,
        "zlabel"           : zLabel,
        "cutBox"           : _cutBoxX,
        "cutBoxY"          : _cutBoxY,
        "rebinX"           : 1, 
        "rebinY"           : 1,
        }

    if 0:
        ROOT.gStyle.SetNdivisions(10, "X")
        ROOT.gStyle.SetNdivisions(10, "Y")
        ROOT.gStyle.SetNdivisions(10, "Z")

    if "VtxIso_Vs_RelIso" in h:
        _kwargs["xlabel"] = "vertex isolation (cm)"# / %.2f (cm)"
        _kwargs["ylabel"] = "relative isolation"# / %.2f"
        #_kwargs["opts"]   = {"xmin": 0.0, "xmax": 3.2, "ymin": 0.0, "ymax": 5.0}
        #_kwargs["opts"]   = {"xmin": 0.0, "xmax": 0.6, "ymin": 0.0, "ymax": 3.0}
        _kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0, "ymin": 0.0, "ymax": 3.0}
        _kwargs["rebinX"] = 1
        _kwargs["rebinY"] = 1
        _kwargs["cutBox"]  = {"cutValue": 0.50, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        _kwargs["cutBoxY"] = {"cutValue": 0.20, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "DonutRatio" in h:
        _kwargs["opts"]   = {"xmin": 0.1, "xmax": 2.0, "ymin": 0.0, "ymax": 2.0}
        _kwargs["rebinX"] = 1
        _kwargs["rebinY"] = 1

    if "VisEt_Vs_dRMax" in h:
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 0.3, "ymin": 0.0, "ymax": 120.0}
        _kwargs["cutBox"]  = {"cutValue":  0.17, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        _kwargs["cutBoxY"] = {"cutValue": 20.00, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "PtLdg_Vs_dRMax" in h:
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 0.3, "ymin": 0.0, "ymax": 100.0}
        _kwargs["cutBox"]  = {"cutValue": 0.25, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        _kwargs["cutBoxY"] = {"cutValue": 5.00, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        
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


def getHisto(datasetsMgr, datasetName, name):

    h1 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name)
    h1.setName("h0" + "-" + datasetName)
    return h1

def getHistos(datasetsMgr, histoName):

    h1 = datasetsMgr.getDataset("Data").getDatasetRootHisto(histoName)
    h1.setName("Data")

    h2 = datasetsMgr.getDataset("EWK").getDatasetRootHisto(histoName)
    h2.setName("EWK")
    return [h1, h2]

def getCustomTGraph(histoName, constant=0.0, coeff=1.0, xmin = 0, xmax = 200, step = 1):

    # Define variables
    xList = []
    yList = []        

    # For-loop: x-y values
    for x in range(xmin, xmax*step):

        if x == 0.0:
            continue

        # Calculate x and y value
        if "GenP_VisEt_Vs" in histoName or "GenP_PtLdg_Vs" in histoName:
            xVal = float(x)/float(step)
            yVal = (constant)/(coeff*xVal)
        else:
            xVal    = float(x)/float(step)
            bracket = constant+xVal

            if xVal >= 0.5:
                if bracket > 0:
                    yVal = coeff*pow(bracket, 0.15)
                else:
                    yVal = 0.0
            else:
                #yVal = 0.5-2*pow(xVal,2) 
                yVal = 0.2-1*pow(xVal,3)
                #yVal = 0.5-xVal
            
        # Debug?
        if 0:
            print "%.2f : %.2f" % (xVal, yVal)
            
        # Save values in lists (must later convert to arrays)
        if yVal > 0.0:
            yList.append( yVal )
            xList.append( xVal )
        
    if len(xList) > 0:
        gr = ROOT.TGraph(len(xList), array.array('d', xList), array.array('d', yList))
    else:
        gr = ROOT.TGraph()
    return gr

def Plot2dHistograms(datasetsMgr, dsetName, histoName, index):

    msg = "%s%s (%s)%s" % (ShellStyles.SuccessStyle(), histoName, dsetName, ShellStyles.NormalStyle() )
    aux.PrintFlushed(msg, index==1)

    # Custom Filtering of datasets 
    dsetsMgr = datasetsMgr.deepCopy()
    if opts.verbose:
        dsetsMgr.PrintInfo()

    # Get Histogram name and its kwargs
    saveName = histoName.rsplit("/")[-1] + "_" + dsetName.split("_")[0] + dsetName.split("_")[-1] 
    kwargs_  = GetHistoKwargs(saveName, opts)

    for i, d in enumerate(dsetsMgr.getAllDatasetNames(), 0):
        if d == dsetName:
            continue
        else:
            # Remove dataset from manager but do NOT close the file!
            dsetsMgr.remove(d, close=False)

    # Sanity check
    nDatasets = len(dsetsMgr.getAllDatasets())
    if nDatasets > 1:
        raise Exception("More than 1 datasets detected in the dataset manager! Can only support 1 dataset. Please use the -i option to choose exactly 1 dataset'")

    # Get the reference histo and the list of histos to compare
    datasets0 = dsetsMgr.getAllDatasets()[0].getName()
    histoList = [getHisto(dsetsMgr, datasets0, histoName)]

    # Create the 2d plot  
    Verbose("Creating the 2d plot", True)
    if opts.normalizeToLumi:
        p = plots.MCPlot(dsetsMgr, histoName, normalizeToLumi=opts.intLumi, saveFormats=[])
    elif opts.normalizeByCrossSection:
        p = plots.MCPlot(dsetsMgr, histoName, normalizeByCrossSection=True, saveFormats=[], **{})
    elif opts.normalizeToOne:
        p = plots.MCPlot(dsetsMgr, histoName, normalizeToOne=True, saveFormats=[], **{})
    else:
        raise Exception("One of the options --normalizeToOne, --normalizeByCrossSection, --normalizeToLumi must be enabled (set to \"True\").")

    Verbose("Setting universal histo styles", True)
    p.histoMgr.setHistoDrawStyleAll("COLZ")
    #p.histoMgr.setHistoLegendStyleAll("L")

    Verbose("Customising histograms", True)
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetZaxis().SetTitleOffset(1.3)) #fixme

    Verbose("Setting plot styles to histograms", True)
    for index, h in enumerate(p.histoMgr.getHistos()):
        plots._plotStyles[p.histoMgr.getHistos()[index].getDataset().getName()].apply(p.histoMgr.getHistos()[index].getRootHisto())

    Verbose("Drawing the plot", True)
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary
    
    # Add fit line for shrinking cone?
    const= 0.0
    coeff= 0.0
    step = 1    
    xmin = 0
    xmax = 0
    if "GenP_VisEt_Vs" in histoName:
        const=   3.5
        coeff=   1.0        
        step =   100    
        xmin =   0
        xmax =  30
    if "GenP_PtLdg_Vs" in histoName:
        const=   2.5 #2.0
        coeff=   1.0
        step = 100
        xmin =   0
        xmax =  30
    if "VtxIso_Vs_RelIso" in histoName:
        const= -0.5
        coeff=  0.3 #0.4
        step =  100
        xmin =   0
        xmax =   0#3

    if "GenP_VisEt_Vs" in histoName or "GenP_PtLdg_Vs" in histoName:
        gr = getCustomTGraph(histoName, const, coeff, xmin, xmax, step)
        gr.SetLineWidth(3)
        gr.Draw("L same")

    Verbose("Removing the legend", True)
    p.removeLegend()

    Verbose("Adding text on canvas", True)
    histograms.addText(0.22, 0.89, plots._legendLabels[datasets0], 18)
    #histograms.addText(0.5, 0.89, plots._legendLabels[datasets0], 18)

    Verbose("Saving the canvas", True)
    aux.SavePlot(p, opts.saveDir, saveName, opts.saveFormats, opts.url)

    return


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
    VERBOSE      = False
    ANALYSIS     = "HLTausAnalysis"
    SEARCHMODE   = None
    DATAERA      = "TDR2019"
    BATCHMODE    = True
    GRIDX        = False
    GRIDY        = False
    LOGX         = False
    LOGY         = False
    LOGZ         = False
    OPTMODE      = None
    INTLUMI      = 1000.0
    SAVEFORMATS = [".png"] #[".C", ".png", ".pdf"]
    URL          = False
    SAVEDIR      = None
    FOLDER       = ""
    NORM2ONE     = False
    NORM2XSEC    = False
    NORM2LUMI    = False

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysisName", dest="analysisName", type="string", default=ANALYSIS,
                      help="Override default analysisName [default: %s]" % ANALYSIS)

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX, 
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY, 
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX, 
                      help="Set the x-axis to log scale? [default: %s]" % LOGX)

    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY, 
                      help="Set the y-axis to log scale? [default: %s]" % LOGY)

    parser.add_option("--logZ", dest="logZ", action="store_true", default=LOGZ, 
                      help="Set the z-axis to log scale? [default: %s]" % LOGZ)

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

    parser.add_option("--normalizeToOne", dest="normalizeToOne", action="store_true", default=NORM2ONE,
                      help="Normalise plot to one [default: %s]" % NORM2ONE)
    
    parser.add_option("--normalizeByCrossSection", dest="normalizeByCrossSection", action="store_true", default=NORM2XSEC,
                      help="Normalise plot by cross-section [default: %s]" % NORM2XSEC)
    
    parser.add_option("--normalizeToLumi", dest="normalizeToLumi", action="store_true", default=NORM2LUMI,
                      help="Normalise plot to luminosity [default: %s]" % NORM2LUMI)

    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

  # Determine path for saving plots                                                                                                                                                                                              
    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="TH2")
    else:
        print "opts.saveDir = ", opts.saveDir

    if opts.mcrab == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)

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
        raw_input("=== plotTH2.py: Press any key to quit ROOT ...")

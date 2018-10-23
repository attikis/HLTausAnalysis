#!/usr/bin/env python
'''
DESCRIPTION:
Script that plots Data/MC for all histograms under a given folder (passsed as option to the script)
Good for sanity checks for key points in the cut-flow


USAGE:
./plotResolution.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotResolution.py -m multicrab_CaloTkSkim_v92X_20180711T1118/ --gridX --gridY --formats .png,.pdf,.C
./plotResolution.py -m multicrab_CaloTkSkim_v92X_20180711T1118/
./plotResolution.py -i TT_TuneCUETP8M2T4_14TeV_L1TPU140 --url -m multicrab_Tracking_v92X_FitParams4Pt2Eta999Stubs0_14h58m32s_15Oct2018


LAST USED:
./plotResolution.py -i TT_TuneCUETP8M2T4_14TeV_L1TnoPU --gridX --url -m 

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
import re
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
def GetLegendDictionary():
    legendDict = {}
    legendDict["C"] = "|#eta| < 0.8"         # C=Central
    legendDict["I"] = "0.8 < |#eta| < 1.6"   # I=Intermediate
    legendDict["F"] = "|#eta| > 1.6"         # F=Forward
    legendDict["L"] = "p_{T} < 5 GeV/c"      # L=Low pT
    legendDict["M"] = "5 < p_{T} < 15 GeV/c" # M=Medium pT
    legendDict["H"] = "p_{T} > 15 GeV/c"     # H=High pT
    return legendDict

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

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    
    optModes = [""]
    if opts.optMode != None:
        optModes = [opts.optMode]

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
    if 0:
        datasetsMgr.selectAndReorder(newOrder)

    # Sanity check (exactly 1 dataset)
    nDsets = len(datasetsMgr.getAllDatasetNames())
    if nDsets > 1:
        msg = "Must have exactly 1 dataset (got %d). Please use the -i and -e options to choose exactly 1 dataset!" % (nDsets)
        Print(msg, True)
        sys.exit()    

    # Print dataset information (after merge)
    if 0:
        datasetsMgr.PrintInfo() #Requires python 2.7.6 or 2.6.6

    # Plot Histograms
    histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(opts.folder)
    histoPaths = [os.path.join(opts.folder, h) for h in histoList]
    skipList   = ["eff_", "counter", "match_trk", "match_tp", "tp_pt", "tp_eta", "resVsEta_ptRel", "resVsPt_ptRel"]
    histoList  = []
    allowedReg = ["C", "I", "F", "L", "M", "H"]

    # For-loop: All histograms (in region)
    for h in histoPaths:
        skip = False
        for s in skipList:
            if s in h.lower():
                skip = True
        if skip:
            continue

        region = h.split("_")[-1]
        if region in allowedReg:
            hName = h.replace(region, "")
            histoList.append(hName)

    # Plot the resolution histograms
    uniqueList = set(histoList)
    
    # For-loop: All histo (generic) names
    for i, h in enumerate(uniqueList, 1):
        if h.endswith("_"):
            h = h[:-1]

        aux.Print("%d/%d: %s" % (i, len(uniqueList), h), i==1)
        PlotHistograms(datasetsMgr, h)

    #print
    Print("All plots saved under directory %s" % (ts + aux.convertToURL(opts.saveDir, opts.url) + ns), True)
        
    return

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


def GetHistoKwargs(h, opts):

    mvLeg  = {"dx": -0.1, "dy": -0.04, "dh": -0.15}
    logY   = False
    yMin   = 0
    yMaxF  = 1.2
    kwargs = {
        "rebinX"           : 1,
        "rebinY"           : None,
        "ratioYlabel"      : "1/Ratio ",
        "ratio"            : opts.ratio,
        "stackMCHistograms": False,
        "ratioInvert"      : True, 
        "addMCUncertainty" : False, 
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "opts"             : {"ymin": yMin, "ymaxfactor": yMaxF},
        "opts2"            : {"ymin": 0.0, "ymax": 2.1},
        "log"              : logY,
        "moveLegend"       : mvLeg,
        "xtitlesize"       : 0.1,#xlabelSize,
        "ytitlesize"       : 0.1,#ylabelSize,
        "cutBox"           : {"cutValue": 0.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        }

        
    if ("resVsPt_" in h):
        ROOT.gStyle.SetNdivisions(6, "X")
        kwargs["opts"] = {"xmin": 0.0, "xmax": 50.0, "ymin": 0.0, "ymaxfactor": 1.4}
        format = "%.0f GeV/c"
        kwargs["cmsTextPosition"] = "left"
        # The ROOT.TGaxis.SetMaxDigits sets the maximum number of digits permitted for the axis labels above which the notation with 10^N is used. 
        # For example, to accept 6 digits number like 900000 on an axis call TGaxis::SetMaxDigits(6). The default value is 5. fgMaxDigits must be greater than 0.
        ROOT.TGaxis.SetMaxDigits(5) #default

        if ("resVsPt_pt" in h):
            kwargs["ylabel"] = "p_{T} resolution / " + format
            # kwargs["ylabel"] = "(p_{T}^{tk} - p_{T}^{tp})" + format

        if ("resVsPt_ptRel" in h):
            ROOT.TGaxis.SetMaxDigits(2) # force scientific scale        
            kwargs["ylabel"] = "relative p_{T} resolution / " + format
            #kwargs["ylabel"] = "(p_{T}^{tk} - p_{T}^{tp}) / p_{T}^{tp} / " + format
            
        if ("resVsPt_eta" in h):
            ROOT.TGaxis.SetMaxDigits(2) # force scientific scale
            kwargs["ylabel"] = "#eta resolution / " + format
            #kwargs["ylabel"] = "(eta^{tk} - eta^{tp})" + format

        if ("resVsPt_phi" in h):
            ROOT.TGaxis.SetMaxDigits(2) # force scientific scale
            kwargs["ylabel"] = "#phi resolution / " + format

        if ("resVsPt_z0" in h):
            #ROOT.TGaxis.SetMaxDigits(2) # force scientific scale
            kwargs["ylabel"] = "z_{0} resolution / " + format

        if ("resVsPt_d0" in h):
            ROOT.TGaxis.SetMaxDigits(2) # force scientific scale
            kwargs["ylabel"] = "d_{0} resolution / " + format

    if ("resVsEta_" in h):
        ROOT.gStyle.SetNdivisions(6, "X")
        format = "%.1f"
        kwargs["cmsTextPosition"] = "left"
        ROOT.TGaxis.SetMaxDigits(2) #default is 5

        if ("resVsEta_pt" in h):
            kwargs["ylabel"] = "p_{T} resolution / " + format

        if ("resVsEta_ptRel" in h):
            kwargs["ylabel"] = "relative p_{T} resolution / " + format

        if ("resVsEta_eta" in h):
            kwargs["ylabel"] = "#eta resolution / " + format
            kwargs["opts"]["ymaxfactor"] = 1.6

        if ("resVsEta_phi" in h):
            kwargs["ylabel"] = "#phi resolution / " + format
            kwargs["opts"]["ymaxfactor"] = 1.6

        if ("resVsEta_z0" in h):
            kwargs["ylabel"] = "z_{0} resolution / " + format
            kwargs["opts"]["ymaxfactor"] = 1.5

        if ("resVsEta_d0" in h):
            kwargs["ylabel"] = "d_{0} resolution / " + format

    if ("res_" in h):
        ROOT.gStyle.SetNdivisions(6, "X")
        kwargs["log"]    = True
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

        if ("res_pt" in h): # = (tk_Pt - tp_Pt)
            kwargs["ylabel"] = "Arbitrary units / %.2f GeV/c"
            kwargs["opts"] = {"xmin": -3.0, "xmax": 3.0, "ymin": 1e-4, "ymaxfactor": 2.0}

        if ("res_ptRel" in h): # = (tk_Pt - tp_Pt)/(tp_Pt)
            kwargs["ylabel"] = "Arbitrary units / %.3f"
            kwargs["opts"] = {"xmin": -0.2, "xmax": 0.2, "ymin": 1e-4, "ymaxfactor": 2.0}

        if ("res_eta" in h): # = (tk_Eta - tp_Eta)
            kwargs["ylabel"] = "Arbitrary units / %.5f"
            kwargs["opts"]   = {"xmin": -0.02, "xmax": 0.02, "ymin": 1e-4, "ymaxfactor": 2.0}

        if ("res_phi" in h): # = (tk_Phi - tp_Phi)
            kwargs["ylabel"] = "Arbitrary units / %.5f rads"
            kwargs["opts"]   = {"xmin": -0.02, "xmax": 0.02, "ymin": 1e-4, "ymaxfactor": 2.0}

        if ("res_z0" in h): # = (tk_Phi - tp_Phi)
            kwargs["ylabel"] = "Arbitrary units / %.3f cm"
            kwargs["opts"] = {"xmin": -1.0, "xmax": 1.0, "ymin": 1e-4, "ymaxfactor": 2.0}

        if ("res_d0" in h): # = (tk_d0 - tp_d0)
            kwargs["ylabel"] = "Arbitrary units / %.4f cm"
            kwargs["opts"] = {"xmin": -1.0, "xmax": 1.0, "ymin": 1e-4, "ymaxfactor": 2.0}

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

def getResolutionHistos(datasetsMgr, datasetName, name, regions=["_C", "_I", "_F"]):
    
    name_C = name + regions[0]
    h1 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name_C)
    h1.setName(name_C)

    name_I = name + regions[1]
    h2 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name_I)
    h2.setName(name_I)

    name_F = name + regions[2]
    h3 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name_F)
    h3.setName(name_F)
    return [h1, h2, h3]

def getHistos(datasetsMgr, datasetName, histoNames, skipIndex=0):

    histos = []
    for i, name in enumerate(histoNames):
        if i==skipIndex:
            continue
        h = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name)
        hName ="h%s-%s" % (i, name)
        h.setName(hName)
        histos.append(h)
    return histos

def PlotHistograms(datasetsMgr, histoName):

    # Customise histo
    kwargs_= GetHistoKwargs(histoName, opts)

    # Get dataset
    datasets0 = datasetsMgr.getAllDatasets()[0].getName()

    # Get histograms for the ComparisonPlot
    if "resVsEta" in histoName:
        myHistos   = getResolutionHistos(datasetsMgr, datasets0, histoName, regions=["_L", "_M", "_H"] )
    else:
        myHistos   = getResolutionHistos(datasetsMgr, datasets0, histoName, regions=["_C", "_I", "_F"] )

    # Get reference histo and list of comparison histos
    hReference = myHistos[0]
    hCompares  = myHistos[1:]

    # Create the plotting object
    p = plots.ComparisonManyPlot(hReference, hCompares, saveFormats=[])

    # Get legend names
    legDict = GetLegendDictionary()
    legends = {}

    # For-loop: All histograms in manager
    for i, h in enumerate(p.histoMgr.getHistos(), 1):
        region = h.getName().split("_")[-1]
            
        if region not in legDict.keys():
            raise Exception("Invalid region \"%s\". Cannot determine legend name for histogram \"%s\"" % (region, h.getName()) )
        legends[h.getName()] = legDict[region]
        p.histoMgr.forHisto(h.getName(), styles.getRegionStyle(region))
        p.histoMgr.setHistoDrawStyle(h.getName(), "LP")
        p.histoMgr.setHistoLegendStyle(h.getName(), "LP")

    # Normalise resolution plots
    if "res_" in histoName:
        if h.getRootHisto().Integral() > 0:
            aux.Verbose("Normalising \"%s\" to unit area" % (histoName), True)
            p.histoMgr.forEachHisto(lambda h: h.getRootHisto().Scale(1.0/h.getRootHisto().Integral()) )
                       
    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legends)

    # Draw and save the plot
    plots.drawPlot(p, histoName, **kwargs_) #the "**" unpacks the kwargs_ dictionary
    
    # Additional text on canvas
    histograms.addText(0.65, 0.90, plots._legendLabels[datasetsMgr.getAllDatasetNames()[0]], 17)

    # Save the plots in custom list of saveFormats
    aux.SavePlot(p, opts.saveDir, histoName, opts.saveFormats, opts.url)
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
    ANALYSISNAME = None #"FakeBMeasurement"
    SEARCHMODE   = None #"80to1000"
    DATAERA      = None #"ID2017"
    GRIDX        = False
    GRIDY        = False
    OPTMODE      = None
    BATCHMODE    = True
    URL          = False
    SAVEDIR      = None
    VERBOSE      = False
    RATIO        = False
    SAVEFORMATS = [".png"]
    FOLDER      = ""

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysisName", dest="analysisName", type="string", default=ANALYSISNAME,
                      help="Override default analysisName [default: %s]" % ANALYSISNAME)

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
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("--formats", dest="formats", default = None,
                      help="Formats in which all plots will be saved in. Provide as list of comma-separated (NO SPACE!) formats. [default: None]")

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("--ratio", dest="ratio", action="store_true", default = RATIO,
                      help="Draw ratio canvas for Data/MC curves? [default: %s]" % (RATIO) )

    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )
    
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
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="hltaus/", postfix="")
    else:
        print "opts.saveDir = ", opts.saveDir

    # Overwrite default save formats?
    if opts.formats != None:
        opts.saveFormats = opts.formats.split(",")
    else:
        opts.saveFormats = SAVEFORMATS


    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotPtEtaRegions.py: Press any key to quit ROOT ...")

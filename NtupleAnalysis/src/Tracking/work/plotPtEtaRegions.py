#!/usr/bin/env python
'''
Description:
Script that plots Data/MC for all histograms under a given folder (passsed as option to the script)
Good for sanity checks for key points in the cut-flow

Usage:
./plot_Folder.py -m <pseudo_mcrab_directory> [opts]

Examples:
./plot_Folder.py -m <mcrab> --folder ""
./plot_Folder.py -m FakeBMeasurement_GE2Medium_GE1Loose0p80_StdSelections_BDTm0p80_AllSelections_BDT0p90_RandomSort_171120_100657 --url --folder ForFakeBMeasurementEWKFakeB --nostack

Last Used:

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
import HLTausAnalysis.NtupleAnalysis.tools.counter as counter
import HLTausAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HLTausAnalysis.NtupleAnalysis.tools.styles as styles
import HLTausAnalysis.NtupleAnalysis.tools.plots as plots
import HLTausAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HLTausAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HLTausAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles


#================================================================================================
# Initializations
#================================================================================================         
def Init():
    
    global etaRegions, ptRegions, etaRegNames, ptRegNames

    etaRegions = ["I", "C", "F"]
    ptRegions = ["L", "M", "H"]
    
    etaRegNames = ["match_tp_pt_",
                   "match_trk_nstub_",
                   "eff_pt_",
                   "res_pt_",
                   "res_eta_",
                   "res_phi_",
                   "res_z0_",
                   "res_d0_",
                   "resVsPt_pt_",
                   "resVsPt_ptRel_",
                   "resVsPt_eta_",
                   "resVsPt_phi_",
                   "resVsPt_z0_",
                   "resVsPt_d0_"]

    ptRegNames = ["match_tp_eta_", 
                  "match_trk_chi2_", 
                  "eff_eta_", 
                  "resVsEta_pt_",
                  "resVsEta_ptRel_",
                  "resVsEta_eta_",
                  "resVsEta_phi_",
                  "resVsEta_z0_",
                  "resVsEta_d0_",
                  ]


    return


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

def GetListOfEwkDatasets():
    Verbose("Getting list of EWK datasets")
    if 0: # TopSelection
        return  ["TT", "WJetsToQQ_HT_600ToInf", "SingleTop", "DYJetsToQQHT", "TTZToQQ",  "TTWJetsToQQ", "Diboson", "TTTT"]
    else: # TopSelectionBDT
        return  ["TT", "noTop", "SingleTop", "ttX"]


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

    Init()
    #optModes = ["", "OptChiSqrCutValue50", "OptChiSqrCutValue100"]
    optModes = [""]

    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # For-loop: All opt Mode
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted() #marina
        if 0:
            datasetsMgr.loadLuminosities() # from lumi.json

        # Set/Overwrite cross-sections
        datasetsToRemove = []
        for d in datasetsMgr.getAllDatasets():
            datasetsMgr.getDataset(d.getName()).setCrossSection(1.0)

        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()

        # Custom Filtering of datasets 
        for i, d in enumerate(datasetsToRemove, 0):
            msg = "Removing dataset %s" % d
            Print(ShellStyles.WarningLabel() + msg + ShellStyles.NormalStyle(), i==0)
            datasetsMgr.remove(filter(lambda name: d in name, datasetsMgr.getAllDatasetNames()))
        #if opts.verbose:
            #datasetsMgr.PrintInfo()

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 
        
        # Get Luminosity
        if 0:
            intLumi = datasetsMgr.getDataset("Data").getLuminosity()

        # Re-order datasets (different for inverted than default=baseline)
        newOrder = []
        # For-loop: All MC datasets
        for d in datasetsMgr.getMCDatasets():
            newOrder.append(d.getName())
            
        # Apply new dataset order!
        datasetsMgr.selectAndReorder(newOrder)
        
        # Print dataset information
        #datasetsMgr.PrintInfo()
        
        # Apply TDR style
        style = tdrstyle.TDRStyle()
        #style.setOptStat(True)
        style.setGridX(opts.gridX)
        style.setGridY(opts.gridY)

        # Plot Histograms
        folder     = "" #opts.folder
        histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(folder)
        histoPaths = [os.path.join(folder, h) for h in histoList]

        # Get the histo names depending on the region we want to plot 
        if opts.etaReg:
            regHistoNames = etaRegNames
            regions = etaRegions
        elif  opts.ptReg:
            regHistoNames = ptRegNames
            regions = ptRegions
        else:
            raise Exception("No region selected. Options: --pt (p_{T} region), --eta (#eta region)")

        # Loop over all the histo names for the specific region 
        for regHisto in regHistoNames:
            histosToPlot     = []
            for i, h in enumerate(histoPaths, 1):
                histoType  = str(type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram()))
                if "TH1" not in histoType:
                    continue
                elif h.split("_")[-1] not in regions:
                    continue
                if regHisto in h:
                    histosToPlot.append(h)
            PlotHistograms(datasetsMgr, histosToPlot[0:3])
        
        
    return




def GetHistoKwargs(h, opts):
    #_moveLegend = {"dx": -0.1, "dy": -0.05, "dh": -0.15}
    _moveLegend = {"dx": -0.52, "dy": -0.05, "dh": -0.15}
    
    logY    = False
    if ("match_tp_" in h) or ("match_trk_" in h) or ("res_" in h):
        logY = True
        

    yMin    = 0.0
    if logY:
        yMin = 0.0001
        yMaxF = 10.0
    else:
        yMaxF = 1.2
        
    _kwargs = {
        #"ylabel"           : _yLabel,
        "rebinX"           : 1,
        "rebinY"           : None,
        "ratioYlabel"      : "Ratio",
        "ratio"            : opts.ratio,
        "stackMCHistograms": opts.nostack,
        "ratioInvert"      : False, 
        "addMCUncertainty" : False, 
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "opts"             : {"ymin": yMin, "ymaxfactor": yMaxF},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : logY,
        "moveLegend"       : _moveLegend,
        "xtitlesize"       : 0.1,#xlabelSize,
        "ytitlesize"       : 0.1,#ylabelSize,
        }

    kwargs = copy.deepcopy(_kwargs)

    '''

    if "_eta" in h.lower():
        #_yLabel = "Arbitrary Units / %.0f "
        units            = ""
        kwargs["xlabel"] = "#eta" 
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -2.5, "xmax": 2.5, "ymin": yMin, "ymaxfactor": yMaxF}

    if "phi" in h.lower():
        #_yLabel = "Arbitrary Units / %.0f "
        units            = "rad"
        kwargs["xlabel"] = "#phi (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -3.15, "xmax": 3.15, "ymin": yMin, "ymaxfactor": yMaxF}
        '''

    if ("match_trk_nstub_" in h) or ("match_trk_chi2_" in h):
        kwargs["moveLegend"] = {"dx": -0.52, "dy": -0.55, "dh": -0.15}

    if ("match_trk_nstub_" in h):
        kwargs["opts"]   = {"xmin": 0, "xmax": 15, "ymin": yMin, "ymax" : 1, "ymaxfactor": yMaxF}

    if ("match_trk_chi2_" in h):
        kwargs["opts"]   = {"xmin": 0, "xmax": 200, "ymin": yMin, "ymax" : 1, "ymaxfactor": yMaxF}

    if ("res_pt_" in h):
        kwargs["opts"]   = {"xmin": -5.0, "xmax": 5.0, "ymin": yMin, "ymax" : 1, "ymaxfactor": yMaxF}

    if ("res_eta_" in h):
        kwargs["opts"]   = {"xmin": -0.01, "xmax": 0.01, "ymin": yMin,"ymax" : 1, "ymaxfactor": yMaxF}
        
    if ("res_phi_" in h):
        kwargs["opts"]   = {"xmin": -0.005, "xmax": 0.005, "ymin": yMin, "ymax" : 1, "ymaxfactor": yMaxF}

    if ("res_z0_" in h):
        kwargs["opts"]   = {"xmin": -1.0, "xmax": 1.0, "ymin": yMin, "ymax" : 1, "ymaxfactor": yMaxF}

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

def getHisto(datasetsMgr, datasetName, name):

    h1 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name)
    h1.setName("h0" + "-" + name)
    return h1

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


def PlotHistograms(datasetsMgr, histoNames):
    Verbose("Plotting Data-MC Histograms")

    ROOT.gStyle.SetNdivisions(8, "X")
    ROOT.gStyle.SetNdivisions(8, "Y")

    # Get Histogram saveNname and its kwargs
    saveName = histoNames[0].rsplit("/")[-1]
    if opts.etaReg:
        saveName = saveName.replace(saveName.split("_")[-1], "CIF")
    elif opts.ptReg:
        saveName = saveName.replace(saveName.split("_")[-1], "LMH")

    kwargs_  = GetHistoKwargs(saveName, opts)
    kwargs ={}

    # Get dataset
    datasets0 = datasetsMgr.getAllDatasets()[0].getName()

    # Get histograms for the ComparisonPlot
    histoReference = getHisto(datasetsMgr, datasets0, histoNames[0])
    histoCompares  = getHistos(datasetsMgr, datasets0, histoNames)

    if ("res_" in histoNames[0]) or ("match_trk_" in histoNames[0]) or ("match_tp_" in histoNames[0]):
        histoReference.normalizeToOne()
        for h in histoCompares:
            h.normalizeToOne()


    # Create the plotting object
    p = plots.ComparisonManyPlot(histoReference, histoCompares, saveFormats=[])

    # Get Legend Names depending on the regions & customize plot
    legendDict ={}
    for index, h in enumerate(p.histoMgr.getHistos()):
        hName = h.getName()
        if opts.etaReg:
            if "_C" in hName:
                legendDict[hName] = "|#eta| < 0.8"
            elif "_I" in hName:
                legendDict[hName] ="0.8 < |#eta| < 1.6"
            elif "_F" in hName:
                legendDict[hName] ="|#eta| > 1.6"
        elif opts.ptReg:
            if "_L" in hName:
                legendDict[hName] = "p_{T} < 5 GeV/c"
            elif "_M" in hName:
                legendDict[hName] ="5 < p_{T} < 15 GeV/c"
            elif "_H" in hName:
                legendDict[hName] ="p_{T} > 15 GeV/c"


        p.histoMgr.forHisto(hName, styles.getRegionStyle(index))
        p.histoMgr.setHistoDrawStyle(hName, "LP")
        p.histoMgr.setHistoLegendStyle(hName, "LP")
                       
    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legendDict)

    # Get the y axis title
    binWidth = p.histoMgr.getHistos()[0].getRootHisto().GetXaxis().GetBinWidth(0)
    if ("res_" in histoNames[0]) or ("match_trk_" in histoNames[0]) or ("match_tp_" in histoNames[0]):
        kwargs_["ylabel"] = "Arbitrary Units / %s" % (GetBinwidthDecimals(binWidth) % (binWidth))
    else:
        ylabel = p.histoMgr.getHistos()[0].getRootHisto().GetYaxis().GetTitle() + " / %s" % (GetBinwidthDecimals(binWidth) % (binWidth))
        kwargs_["ylabel"] = ylabel
        
                   
    # Draw and save the plot
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary

    # Add extra text (dataset name)
    #histograms.addText(0.65, 0.88, plots._legendLabels[datasets0], 17)
    if ("match_trk_nstub_" in histoNames[0]) or ("match_trk_chi2_" in histoNames[0]):
        histograms.addText(0.23, 0.38, plots._legendLabels[datasets0], 17)
    else:
        histograms.addText(0.23, 0.88, plots._legendLabels[datasets0], 17)
        
    # Save the plots in custom list of saveFormats
    SavePlot(p, saveName, os.path.join(opts.saveDir, opts.optMode, opts.folder), [".pdf",".png"] )
    return


def SavePlot(plot, plotName, saveDir, saveFormats = [".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )

    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_").replace(" ", "").replace("(", "").replace(")", "") )
    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        #saveNameURL = saveNameURL.replace("/publicweb/a/aattikis/", "http://home.fnal.gov/~aattikis/")
        saveNameURL = saveNameURL.replace("/afs/cern.ch/user/m/mtoumazo/public/html/hltaus/", "https://cmsdoc.cern.ch/~mtoumazo/hltaus/")
        if opts.url:
            Print(saveNameURL, 1)
        else:
            Print(saveName + ext, 1)
        plot.saveAs(saveName, formats=saveFormats)
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
    GRIDX        = True
    GRIDY        = True
    OPTMODE      = None
    BATCHMODE    = True
    PRECISION    = 3
    INTLUMI      = -1.0
    SUBCOUNTERS  = False
    LATEX        = False
    URL          = True
    NOERROR      = True
    SAVEDIR      = "/afs/cern.ch/user/m/mtoumazo/public/html/hltaus/Tracking/" #os.getcwd()
    VERBOSE      = False
    FOLDER       = ""
    RATIO        = False
    NOSTACK      = False

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

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
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

    parser.add_option("--nostack", dest="nostack", action="store_true", default = NOSTACK,
                      help="Do not stack MC histograms [default: %s]" % (NOSTACK) )

    parser.add_option("--eta", dest="etaReg", action="store_true", default = False,
                      help="Make #eta Region Plots (Central, Intermediate, Forward) [default: False]")

    parser.add_option("--pt", dest="ptReg", action="store_true", default = False,
                      help="Make p_{T} Region Plots (Low, Middle, High) [default: False]")

    
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
        raw_input("=== plot_Folder.py: Press any key to quit ROOT ...")

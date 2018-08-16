#!/usr/bin/env python
'''
Description:
Script that plots Data/MC for all histograms under a given folder (passsed as option to the script)
Good for sanity checks for key points in the cut-flow

Usage:
./plotRateVsEff.py -m <pseudo_mcrab_directory> [opts]

Examples:
./plotRateVsEff.py -m <mcrab> --folder ""
./plotRateVsEff.py -m FakeBMeasurement_GE2Medium_GE1Loose0p80_StdSelections_BDTm0p80_AllSelections_BDT0p90_RandomSort_171120_100657 --url --folder ForFakeBMeasurementEWKFakeB --nostack

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

        # Custom Filtering of datasets 
        datasetsToRemove = []
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
        #datasetSignal = "GluGluHToTauTau_14TeV_L1TPU140" #"TT_TuneCUETP8M2T4_14TeV_L1TPU140"#"GluGluHToTauTau_14TeV_L1TPU140" #"ChargedHiggs200_14TeV_L1TPU140"
        datasetSignal = "TT_TuneCUETP8M2T4_14TeV_L1TPU140"
        datasetBkg    = "SingleNeutrino_L1TPU140"

        # Plot Histograms
        effHistoLists  = [["Calo_Eff", "Tk_Eff", "VtxIso_Eff", "RelIso_Eff"], ["DiTau_Eff_Calo", "DiTau_Eff_Tk", "DiTau_Eff_VtxIso", "DiTau_Eff_RelIso"]]
        rateHistoLists = [["Calo_Rate", "Tk_Rate", "VtxIso_Rate", "RelIso_Rate"], ["DiTau_Rate_Calo", "DiTau_Rate_Tk", "DiTau_Rate_VtxIso", "DiTau_Rate_RelIso"]]

        # For-loop: All histos
        for i in range(0, len(effHistoLists)):
            eff  = effHistoLists[i]
            rate = rateHistoLists[i]
            PlotRateVsEff(datasetsMgr, eff, rate, datasetSignal, datasetBkg)            

    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return


def PlotRateVsEff(datasetsMgr, effHistoList, rateHistoList, datasetSignal, datasetBkg):
    tgraphs=[]
    legendDict = {}

    # Get Histogram name and its kwargs
    if "ditau" in effHistoList[0].lower():
        saveName = "DiTau_RateVsEff_"+ datasetSignal.split("_")[0]
    else:
        saveName = "SingleTau_RateVsEff_"+ datasetSignal.split("_")[0]
    kwargs_  = GetHistoKwargs(saveName, opts)

    for i in range (0, len(effHistoList)):
        if (i==0) :
            g0 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], datasetSignal, datasetBkg)
            g0.SetName("Calo")
        elif (i==1):
            g1 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], datasetSignal, datasetBkg)
            g1.SetName("Tk")
        elif (i==2):
            g2 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], datasetSignal, datasetBkg)
            g2.SetName("VtxIso")
        elif (i==3):
            g3 = convert2RateVsEffTGraph(datasetsMgr, effHistoList[i], rateHistoList[i], datasetSignal, datasetBkg)
            g3.SetName("RelIso")


    # Create & draw the plot
    p = plots.ComparisonManyPlot(g0, [g1,g2, g3], saveFormats=[])
    
    # Set individual styles
    for index, h in enumerate(p.histoMgr.getHistos()):
        hName = h.getName()
        legendDict[hName] = styles.getCaloLegend(index)
        p.histoMgr.forHisto(hName, styles.getCaloStyle(index))

        p.histoMgr.setHistoDrawStyle(h.getName(), "LX")                                                                                                     
        p.histoMgr.setHistoLegendStyle(h.getName(), "LP")

    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legendDict)

    # Draw and save the plot
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary
    # Draw
    for i, g in enumerate([g0, g1, g2, g3]):
        shapes, min, max = DrawErrorBand(g) 
        for shape in shapes:
            shape.SetFillColor( p.histoMgr.getHistos()[i].getRootHisto().GetFillColor())
            shape.SetFillStyle(3003)
            shape.Draw("f same")
        ROOT.gPad.RedrawAxis()


    # Additional text
    histograms.addText(0.65, 0.38, plots._legendLabels[datasetSignal], 17)
    # if "ditau" in saveName.lower():
    #     histograms.addText(0.77, 0.89,"DoubleTau", 20)
    # else:
    #     histograms.addText(0.77, 0.89,"SingleTau", 20)

    # Save the plots in custom list of saveFormats
    aux.SavePlot(p, opts.saveDir, saveName, opts.saveFormats, True)
    return

def convert2RateVsEffTGraph(datasetsMgr, effHistoName, rateHistoName, datasetSignal, datasetBkg):

    hEff  = datasetsMgr.getDataset(datasetSignal).getDatasetRootHisto(effHistoName).getHistogram()
    hRate = datasetsMgr.getDataset(datasetBkg).getDatasetRootHisto(rateHistoName).getHistogram()
    
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
    _moveLegend = {"dx": -0.1, "dy": -0.55, "dh": -0.15}
    logY    = True
    yMin    = 0.0
    if logY:
        yMin = 1
        yMaxF = 10
    else:
        yMaxF = 1.2
        
    _kwargs = {
        "xlabel"           : "Efficiency",
        "ylabel"           : "Rate (kHz)",
        "addMCUncertainty" : False, 
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Phase-2 Simulation",
        "cmsTextPosition"  : "outframe",
        "opts"             : {"xmin": 0.0, "xmax": 1.0, "ymin": yMin, "ymax":1000, "ymaxfactor": yMaxF},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : logY,
        "moveLegend"       : _moveLegend,
        "xtitlesize"       : 0.1,#xlabelSize,
        "ytitlesize"       : 0.1,#ylabelSize,
        "cutBoxY"           : {"cutValue": 50, "fillColor": 16, "box": False, "line": True, "cutGreaterThan"   : False}
        }

    kwargs = copy.deepcopy(_kwargs)
    
    if "ditau" in h.lower():
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 0.6, "ymin": yMin, "ymax":1000, "ymaxfactor": yMaxF}

    return kwargs


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
    ANALYSIS    = "HLTausAnalysis"
    BATCHMODE    = True
    DATAERA     = "ID2017" #"TDR2019"
    FOLDER       = ""
    GRIDX        = True
    GRIDY        = True
    OPTMODE      = None
    PRECISION    = 3
    RATIO        = False
    SAVEDIR      = None
    SAVEFORMATS = [".png"] #[".C", ".png", ".pdf"]
    SEARCHMODE   = None
    URL          = True
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
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="ROC")
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

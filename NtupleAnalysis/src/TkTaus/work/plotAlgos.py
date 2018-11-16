#!/usr/bin/env python
'''
DESCRIPTION:
Performance Comparison (Rate, Eff, ROC, Turn-on) for different Tau algorithms


USAGE:
./plotAlgos.py -m <pseudo_mcrab> [opts]


EXAMPLES:
./plotAlgos.py -m multicrab_CaloTk_v92X_16h12m34s_15Nov2018,multicrab_TkTaus_v92X_14h39m23s_15Nov2018 --bandValue 2
./plotAlgos.py -m multicrab_CaloTk_v92X_16h12m34s_15Nov2018,multicrab_TkTaus_v92X_17h56m36s_15Nov2018,multicrab_TkEG_v92X_17h04m49s_15Nov2018 --gridX --gridY


LAST USED:
./plotAlgos.py -m multicrab_CaloTk_v92X_16h12m34s_15Nov2018,multicrab_TkTaus_v92X_17h56m36s_15Nov2018,multicrab_TkEG_v92X_21h00m14s_15Nov2018

'''
#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
import datetime
from optparse import OptionParser

import getpass

import ROOT
import array

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

ROOT.gErrorIgnoreLevel = ROOT.kError

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

def GetDatasetsFromDir(opts, i):

    aux.Verbose("multicrab = \"%s\"" % (opts.mcrabs[i]), i==0)

    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrabs[i]],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode, 
                                                        analysisName=opts.analysisName,
                                                        optimizationMode=opts.optMode)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrabs[i]],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        includeOnlyTasks=opts.includeOnlyTasks,
                                                        optimizationMode=opts.optMode)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrabs[i]],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        excludeTasks=opts.excludeTasks,
                                                        optimizationMode=opts.optMode)
    else:
        raise Exception("This should never be reached")
    return datasets
    

def GetHistoKwargs(h, opts):

    # Definitions
    _mvLeg1 = {"dx": -0.14, "dy": -0.02, "dh": -0.1}
    _mvLeg2 = {"dx": -0.14, "dy": -0.02, "dh": -0.1}
    _mvLeg3 = {"dx": -0.50, "dy": -0.02, "dh": -0.1}
    _mvLeg4 = {"dx": -0.14, "dy": -0.50, "dh": -0.1}
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
        "opts"             : {"xmin": 0.0, "xmax": 1.0, "ymin": yMin, "ymax":15000, "ymaxfactor": yMaxF},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : logY,
        "moveLegend"       : _mvLeg1,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        "cutBoxY"           : {"cutValue": 50, "fillColor": 16, "box": False, "line": True, "cutGreaterThan"   : False}
        }

    if h == "RateVsEff":
        #_kwargs["opts"]   = {"xmin": 0.0, "xmax": 1.0, "ymin": yMin, "ymax":5000, "ymaxfactor": yMaxF}
        _kwargs["opts"]   = {"xmin": 0.0, "xmax": 0.8, "ymin": yMin, "ymax":5000, "ymaxfactor": yMaxF}
        _kwargs["moveLegend"] = _mvLeg3

    if "_Rate" in h:
        _kwargs["xlabel"]     = "E_{T} (GeV)"
        _kwargs["ylabel"]     = "Rate (kHz)"# / %.0f GeV"
        #_kwargs["opts"]       = {"xmin": 0.0, "xmax": 200.0, "ymin": 1, "ymax":5e4, "ymaxfactor": yMaxF}
        _kwargs["opts"]       = {"xmin": 0.0, "xmax": 160.0, "ymin": 1, "ymax":5e4, "ymaxfactor": yMaxF}
        _kwargs["moveLegend"] = _mvLeg1
        _kwargs["cutBoxY"]    = {"cutValue": 50, "fillColor": 16, "box": False, "line": True, "cutGreaterThan": False}

    if "_Eff" in h:
        units = "GeV"
        _kwargs["xlabel"]     = "E_{T} (%s)" % (units)
        _kwargs["ylabel"]     = "Efficiency / %0.0f " + units
        _kwargs["log"]        = False
        #_kwargs["opts"]       = {"xmin": 0.0, "xmax": 200.0, "ymin": 0.0, "ymax": 1.0, "ymaxfactor": yMaxF}
        _kwargs["opts"]       = {"xmin": 0.0, "xmax": 160.0, "ymin": 0.0, "ymax": 1.0, "ymaxfactor": yMaxF}
        _kwargs["moveLegend"] = _mvLeg2
        _kwargs["cutBoxY"]    = {"cutValue": 50, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
        _kwargs["cutBoxX"]    = {"cutValue": 10, "fillColor": 16, "box": True, "line": True, "cutGreaterThan": False}
    
    if "TurnOn" in h:
        _units = "GeV"
        _kwargs["xlabel"]     = "#tau_{h} E_{T}^{vis} (%s)" % (_units)
        _kwargs["ylabel"]     = "Efficiency / %0.0f " + _units
        _kwargs["log"]        = False
        _kwargs["rebinX"]     = 1 # do NOT change
        #_kwargs["opts"]       = {"xmin": 0.0, "xmax": 200.0, "ymin": 0.0, "ymax": 1.15, "ymaxfactor": yMaxF}
        _kwargs["opts"]       = {"xmin": 0.0, "xmax": 160.0, "ymin": 0.0, "ymax": 1.15, "ymaxfactor": yMaxF}
        _kwargs["cutBoxY"]    = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
        _kwargs["moveLegend"] = _mvLeg3
        if "50" in h:
            _kwargs["cutBox"] = {"cutValue": 50.0, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}
        if "25" in h:
            _kwargs["cutBox"] = {"cutValue": 25.0, "fillColor": 16, "box": False, "line": False, "cutGreaterThan": False}

    return _kwargs


def getAlgoLabel(algo):
    '''
    https://root.cern.ch/doc/master/classTAttText.html
    '''
    labelDict = {}
    # labelDict["CaloTk"] = "#font[72]{Calo+Tracks}"
    # labelDict["TkTaus"] = "#font[72]{Tracks-only}"
    # labelDict["TkEG"]   = "#font[72]{Tracks+e/#gamma}"
    labelDict["CaloTk"] = "Calo+Tracks"
    labelDict["TkTaus"] = "Tracks-only"
    labelDict["TkEG"]   = "Tracks+EG"
    if algo not in labelDict.keys():
        msg = "Could not find algorith \"%s\" in supported keys" % (algo)
        raise Exception(es + msg + ns)
    return labelDict[algo]

def getAlgos():
    return ["#font[72]{Calo+Tracks}", "#font[72]{Tracks-only}", "#font[72]{Tracks+e#gamma}"]


def main(opts):

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)

    # Setup & configure the dataset manager 
    datasetsMgr = GetDatasetsFromDir(opts, 0)
    datasetsMgr.updateNAllEventsToPUWeighted()
    if 0:
        datasetsMgr.loadLuminosities() # from lumi.json
        
    if opts.verbose:
        datasetsMgr.PrintCrossSections()
        datasetsMgr.PrintLuminosities()

    # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
    plots.mergeRenameReorderForDataMC(datasetsMgr) 

    # Print datasets info summary
    datasetsMgr.PrintInfo()

    # Define the mapping histograms in numerator->denominator pairs
    VariableList = ["L1Taus_SingleTau_Eff", "L1Taus_SingleTau_Rate", "L1Taus_TurnOn25", "L1Taus_TurnOn50"]

    counter =  0
    opts.nDatasets = len(datasetsMgr.getAllDatasets())
    nPlots  = len(VariableList)*opts.nDatasets

    # For-loop: All datasets
    for dataset in datasetsMgr.getAllDatasets():
        
        opts.saveDir = aux.getSaveDirPath(opts.saveDirBase, prefix="hltaus/", postfix="ROC")
        PlotRateVsEff(datasetsMgr, dataset, "SingleTau", "PU140")
        PlotRateVsEff(datasetsMgr, dataset, "SingleTau", "PU200")
        PlotRateVsEff(datasetsMgr, dataset, "DiTau"    , "PU140")
        PlotRateVsEff(datasetsMgr, dataset, "DiTau"    , "PU200")

        # For-looop: All variables
        for hName in VariableList:
            hPath = os.path.join(opts.folder, hName)

            counter+=1
            msg = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Histogram", "%i" % counter, "/", "%s:" % (nPlots), "%s" % (dataset.getName()))
            aux.Print(ShellStyles.SuccessStyle() + msg + ShellStyles.NormalStyle(), counter==1)

            if "neutrino" in dataset.getName().lower():        
                if "rate" in hName.lower():
                    opts.saveDir = aux.getSaveDirPath(opts.saveDirBase, prefix="hltaus/", postfix="Rates")
                    PlotHistos(dataset.getName(), hPath, hName.split("_")[0] + "_")
                else:
                    pass
            else:
                if "rate" in hName.lower():
                    continue
                else:
                    if "eff" in hName.lower():
                        opts.saveDir = aux.getSaveDirPath(opts.saveDirBase, prefix="hltaus/", postfix="Efficiencies")
                        PlotHistos(dataset.getName(), hPath, hName.split("_")[0] + "_" )

                    if "turnon" in hName.lower():
                        opts.saveDir = aux.getSaveDirPath(opts.saveDirBase, prefix="hltaus/", postfix="TurnOns")
                        PlotHistos(dataset.getName(), hPath, hName.split("_")[1] + "_" )

    aux.Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return
            

def PlotHistos(datasetName, hPath, saveName):

    _kwargs = GetHistoKwargs(hPath, opts)
    datasetsMgr_list = []

    j = 0
    histoList = []
    legDict   = {}

    # For-loop: All pseudo-multicrabs
    while j < len(opts.mcrabs):

        algo    = opts.mcrabs[j].split("_")[1]
        dMgr    = GetDatasetsFromDir(opts, j)
        dataset = dMgr.getDataset(datasetName)
        dMgr.updateNAllEventsToPUWeighted()

        # Get Histogram from dataset
        histo = dataset.getDatasetRootHisto(hPath).getHistogram()                

        # Set style
        styles.styles[j].apply(histo) 
        legDict[algo] = getAlgoLabel(algo)

        if (j == 0):
            if "turnon" in saveName.lower():
                refHisto = histograms.Histo(histo, algo, legendStyle = "LP", drawStyle="AP")
            else:
                refHisto = histograms.Histo(histo, algo, legendStyle = "L", drawStyle="HIST")
                refHisto.getRootHisto().SetLineStyle(ROOT.kSolid)
                refHisto.getRootHisto().SetLineWidth(4)
        else:
            if "turnon" in saveName.lower():
                histoList.append(histograms.Histo(histo, algo, legendStyle = "LP", drawStyle="AP"))
            else:
                histo.SetLineWidth(4)
                histoList.append(histograms.Histo(histo, algo, legendStyle = "L", drawStyle="HIST"))
        j = j + 1

    # Create the plotter object 
    p = plots.ComparisonManyPlot(refHisto, histoList, saveFormats=[])
    p.setLegendHeader(plots._legendLabels[datasetName])

    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legDict)

    # Draw the plots
    plots.drawPlot(p, opts.saveDir, **_kwargs)
    
    # Add text
    if 0:
        histograms.addText(0.65, 0.1, "#\pm %.1f %% band" % (opts.bandValue), 20)

    # Save plot in all formats
    saveName += datasetName
    aux.SavePlot(p, opts.saveDir, saveName, opts.saveFormats, True)
    return

def PlotRateVsEff(datasetsMgr, dataset, Type, PU):

    # Sanity checks
    if PU not in dataset.getName():
        return
    if "Neutrino" in dataset.getName():
        return

    # Definitions
    _kwargs = GetHistoKwargs("RateVsEff", opts)
    j       = -1
    gList   = []
    legDict = {}

    # For-loop: All pseudo-multicrabs
    while j < len(opts.mcrabs)-1:
        j = j + 1

        dMgr = GetDatasetsFromDir(opts, j)
        dMgr.updateNAllEventsToPUWeighted()

        # Get dataset
        signal  = dataset.getName()
        bkg     = "SingleNeutrino_L1T%s" % (PU)

        # Get Histogram from dataset
        hEff   = os.path.join(opts.folder, "L1Taus_%s_Eff" % (Type) )
        hRate  = os.path.join(opts.folder, "L1Taus_%s_Rate" % (Type) )
        graph  = convert2RateVsEffTGraph(dMgr, hEff, hRate, signal, bkg)
        algo   = opts.mcrabs[j].split("_")[1]
        graph.SetName(algo)
        styles.styles[j].apply(graph) 
        graph.SetLineWidth(4)
        gList.append(graph)
        legDict[algo] = getAlgoLabel(algo)


    # Create the plotter object 
    p = plots.ComparisonManyPlot(gList[0], gList[1:], saveFormats=[])
    p.setLegendHeader(plots._legendLabels[signal])

    # Set individual styles
    for index, h in enumerate(p.histoMgr.getHistos(), 1):
        hName = h.getName()
        legDict[hName] = getAlgoLabel(hName)
        p.histoMgr.setHistoDrawStyle(h.getName(), "LX") # "X" = Do not draw error bars
        p.histoMgr.setHistoLegendStyle(h.getName(), "L") #LP

    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legDict)

    # Draw the plots
    plots.drawPlot(p, opts.saveDir, **_kwargs)
    
    # Add text
    if 0:
        histograms.addText(0.65, 0.1, "#\pm %.1f %% band" % (opts.bandValue), 20)

    # Save plot in all formats
    saveName = Type + "_" + dataset.getName()
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


def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):

     # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        aux.Verbose(saveNameURL, i==0)
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
    BATCHMODE    = True
    ANALYSIS     = "HLTausAnalysis"
    DATAERA      = "TDR2019" #"ID2017" #"TDR2019"
    FOLDER       = ""
    EBANDS       = False
    BANDVALUE    = 5.0
    GRIDX        = False
    GRIDY        = False
    OPTMODE      = None
    PRECISION    = 3
    RATIO        = False
    SAVEDIR      = None
    SEARCHMODE   = None
    SAVEFORMATS  = [".png"] #[".C", ".png", ".pdf"]
    URL          = False
    VERBOSE      = False
    PREFIX       = ""
    POSTFIX      = ""

        
    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrabs", dest="mcrabs", action="store", 
                      help="Path to the multicrab directories for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysisName", dest="analysisName", type="string", default=ANALYSIS,
                      help="Override default analysisName [default: %s]" % ANALYSIS)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX, 
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY, 
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)
    
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

    parser.add_option("--prefix", dest="prefix", type="string", default = PREFIX,
                      help="Prefix string to be appended on save directory [default: %s]" % (PREFIX) )

    parser.add_option("--postfix", dest="postfix", type="string", default = POSTFIX,
                      help="Postfix string to be appended on save directory [default: %s]" % (POSTFIX) )

    parser.add_option("--bandValue", dest="bandValue", type="float", default=BANDVALUE,
                      help="Add a symmetric band around 1.0. Value passed should be the percentage (e.g 10 or 5)  [default: %s]" % (BANDVALUE) )

    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if len(opts.mcrabs) == 0:
        aux.Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)

    # Store all (comma-separated) pseudomulticrabs in a list
    if opts.mcrabs == None:
        aux.Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        sys.exit(1)
    else:
        if "," in opts.mcrabs:
            opts.mcrabs = opts.mcrabs.split(",")    
            aux.Print("Will use the following pseudo-multicrab directories:", True)
            for d in opts.mcrabs:
                aux.Print(d, False)
        else:
            cwd  = os.getcwd()
            dirs = [ name for name in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, name)) ]
            mcrabs = [d for d in dirs if opts.mcrabs in d and "Top" in d]
            opts.mcrabs = mcrabs


    # For-loop: All pseudomulticrab dirs
    for mcrab in opts.mcrabs:
        if not os.path.exists("%s/multicrab.cfg" % mcrab):
            msg = "No pseudo-multicrab directory found at path '%s'! Please check path or specify it with --mcrab!" % (mcrab)
            raise Exception(ShellStyles.ErrorLabel() + msg + ShellStyles.NormalStyle())
        else:
            msg = "Using pseudo-multicrab directory %s" % (ShellStyles.NoteStyle() + mcrab + ShellStyles.NormalStyle())
            aux.Verbose(msg , True)

    # Determine path for saving plots
    opts.saveDirBase = None
    if opts.saveDir == None:
        dirName = "multicrab" #"_".join(opts.mcrabs[0].split("_")[:2]) 
        for m in opts.mcrabs:
            dirName +=  "_" + m.split("_")[1]
        dirName += "_" + datetime.date.today().strftime("%d%h%Y")
        opts.saveDir = aux.getSaveDirPath(dirName, prefix="hltaus/", postfix="")
        opts.saveDirBase = aux.getSaveDirPath(dirName, prefix="hltaus/", postfix="")
    else:
        print "opts.saveDir = ", opts.saveDir
        print "opts.saveDirBase = ", opts.saveDirBase

    # Overwrite default save formats?
    if opts.formats != None:
        opts.saveFormats = opts.formats.split(",")
    else:
        opts.saveFormats = SAVEFORMATS

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotAlgos.py: Press any key to quit ROOT ...")

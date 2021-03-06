#!/usr/bin/env python
'''
Description:
Script that plots the Cumulative histogram for a given distribution for given background and signal sample 

Usage:
./plotOptimisation.py -m <pseudo_mcrab_directory> [opts]

Examples:
./plotOptimisation.py -m <mcrab> 

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
        histoType  = type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram())
        for i, h in enumerate(histoPaths, 1):
            histoType  = str(type(datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDatasetRootHisto(h).getHistogram()))
            if "TH1" not in histoType:
                continue
            if "TkEG_vtxIso"==h or "leadTrks_Pt"== h:
                GetCumulativePlot(datasetsMgr, h)
                GetSignificancePlot(datasetsMgr, h)
            
                

    return

def GetCumulativePlot(datasetsMgr, histoName):

    # Define cut direction
    cutDirection= ">" 

    Verbose("Calculating the cut-efficiency (%s) for histo with name %s" % (cutDirection, histoName) )

    # Declare variables & options                                                                      
    first  = True
    isData = False
    histosCumul = []

    kwargs_  = GetHistoKwargs(histoName, opts)

    # Get the ROOT histogram                                                                           
    for d in datasetsMgr.getAllDatasetNames():
        rootHisto = datasetsMgr.getDataset(d).getDatasetRootHisto(histoName)

        # Normalise the histogram                                                                          
        rootHisto.normalizeToOne()

        # Get a clone of the wrapped histogram normalized as requested.                                   
        h = rootHisto.getHistogram()
    
        # Set the x-y axis titles
        binWidth = h.GetXaxis().GetBinWidth(0)
        kwargs_["xlabel"] = h.GetXaxis().GetTitle()
        #kwargs_["xlabel"] = "Relative Isolation"
        kwargs_["ylabel"] = "Efficiency / Value (%s) / %s" % (cutDirection, GetBinwidthDecimals(binWidth) % (binWidth) )
    
        # If empty return                                                                                  
        if h.GetEntries() == 0:
            return

        # Create the eventsNo histogram                                                                    
        eventsNo   = h.Clone(d)

        # Reset and Edit the eventsNo histogram                                                            
        eventsNo.Reset()

        # Calculate the instances passing a given cut (all bins)                                           
        nBinsX = h.GetNbinsX()+1
        for iBin in range(1, nBinsX):
            
            nTotal = h.Integral(0, nBinsX)
            
            if cutDirection == ">":
                nPass  = h.Integral(iBin, nBinsX)
            elif cutDirection == "<":
                nPass  = nTotal - h.Integral(iBin+1, nBinsX)
            else:
                raise Exception("Invalid cut direction  \"%s\". Please choose either \">\" or \"<\"" % (cutDirection))

            # Sanity check                                                                                 
            if nPass < 0:
                nPass = 0

            # Fill the eventNo histogram                                                                   
            eventsNo.SetBinContent(iBin, nPass)
            eventsNo.SetBinError(iBin, math.sqrt(nPass)/10)
        
        histosCumul.append(eventsNo)


    histoName = histoName +"_Cumulative"
    PlotHistograms(datasetsMgr,histoName,histosCumul, **kwargs_)

    return


def GetSignificancePlot(datasetsMgr, histoName):

    # Define cut direction
    cutDirection= ">" 
    histosSignif = []

    Verbose("Calculating the cut-efficiency (%s) for histo with name %s" % (cutDirection, histoName) )

    # Declare variables & options                                                                      
    first  = True
    isData = False
    Signal_Cumul = []

    kwargs_  = GetHistoKwargs(histoName, opts)

    # Get the ROOT histogram                                                                           
    for d in datasetsMgr.getAllDatasetNames():
        rootHisto = datasetsMgr.getDataset(d).getDatasetRootHisto(histoName)

        # Normalise the histogram                                                                          
        rootHisto.normalizeToOne()

        # Get a clone of the wrapped histogram normalized as requested.                                   
        h = rootHisto.getHistogram()
    
        # Set the x-y axis titles
        binWidth = h.GetXaxis().GetBinWidth(0)
        kwargs_["xlabel"] = h.GetXaxis().GetTitle()
        #kwargs_["xlabel"] = "Relative Isolation"
        kwargs_["ylabel"] = "Efficiency / Value (%s) / %s" % (cutDirection, GetBinwidthDecimals(binWidth) % (binWidth) )
    
        # If empty return                                                                                  
        if h.GetEntries() == 0:
            return

        # Create the eventsNo histogram                                                                    
        eventsNo   = h.Clone(d)

        # Reset and Edit the eventsNo histogram                                                            
        eventsNo.Reset()

        # Calculate the instances passing a given cut (all bins)                                           
        nBinsX = h.GetNbinsX()+1
        for iBin in range(1, nBinsX):
            
            nTotal = h.Integral(0, nBinsX)
            
            if cutDirection == ">":
                nPass  = h.Integral(iBin, nBinsX)
            elif cutDirection == "<":
                nPass  = nTotal - h.Integral(iBin+1, nBinsX)
            else:
                raise Exception("Invalid cut direction  \"%s\". Please choose either \">\" or \"<\"" % (cutDirection))

            # Sanity check                                                                                 
            if nPass < 0:
                nPass = 0

            # Fill the eventNo histogram                                                                   
            eventsNo.SetBinContent(iBin, nPass)
            eventsNo.SetBinError(iBin, math.sqrt(nPass)/10)
        
        if "SingleNeutrino" in d:
            Background_Cumul = eventsNo
            datasetsMgr.remove(d)
        else:
            Signal_Cumul.append(eventsNo)

    for index, sig in enumerate(Signal_Cumul):
        x_signal     = []                                                                                 
        y_signal     = []                                                                                 
        x_background = []                                                                                 
        y_background = []                                                                                 
        
        # Clone the Cumulative plots for Signal and Bkg                                                   
        h_signal     = sig.Clone()                                                               
        h_background = Background_Cumul.Clone()                                                           
                                                                                                      
        # Get the binNo of the Cumulative plots for Signal and Bkg                                         
        n_signal = h_signal.GetNbinsX()                                                                   
        n_background = h_background.GetNbinsX()                                                           
                                                                                                      
        # Get the bin-width of the plots                                                                  
        binWidth = h_signal.GetBinWidth(0)                                                                
                                                                                                      
        # Sanity Check                                                                                    
        if (n_signal != n_background):                                                                    
            print "Warning: Number of Bins of signal and background are different!"                      


        # Get the x-y values of the Cumulative Histogram                                                  
        for i in range(1, n_signal+1):                                                                    
                                                                                                      
            # Signal (x,y) Values                                                                         
            x_signal.append(h_signal.GetBinLowEdge(i)+0.5*h_signal.GetBinWidth(i))                        
            y_signal.append(h_signal.GetBinContent(i))                                                    
            
            # Background (x,y) values and errors                                                          
            x_background.append(h_background.GetBinLowEdge(i)+0.5*h_background.GetBinWidth(i))            
            y_background.append(h_background.GetBinContent(i))                                            

        # Create the (x-y) list to append the significance value for each x value                         
        x = []                                                                                            
        y = []                                                                                            
                
        '''
        # Create variables for Maximum Significance                                                       
        maxSignX = 0                                                                                      
        maxSignY = 0                                                                                      
        '''
        print "="*100
        # Calculate the Significance of the Signal                                                        
        for i in range(0, n_signal):                                                                      
            if ((float(y_background[i]) <= 0 ) or (float(y_signal[i]) <= 0)):                             
                significance = 0                                                                          
            else:
                significance = float(y_signal[i])/math.sqrt(float(y_background[i])+ float(y_signal[i]))
                #significance = float(y_signal[i])/math.sqrt(float(y_background[i]))
                #significance = float(y_signal[i]) / float(y_background[i])                                
                #print i, y_signal[i], y_background[i], significance
            '''
            elif (kwargs.get("significanceDef") == "S1"):                                                 
                significance = float(y_signal[i])/math.sqrt(float(y_background[i]))                       
            elif (kwargs.get("significanceDef") == "S2"):                                                 
                significance = float(y_signal[i]) / float(y_background[i])                                
            '''
            x.append(x_signal[i])
            y.append(significance)

        # Create the Significance Plot         
        tGraph = ROOT.TGraph(n_signal, array.array("d", x), array.array("d", y))                          
        # tGraph = ROOT.TGraph(n_signal, x, y)                          
                                                                                                         
        # Customize y-axis title                                                                          
        ytitle = "S/ #sqrt{B} "                                                                       
        '''
        if (kwargs.get("significanceDef") == "S1"):                                                       
        ytitle = "S/ #sqrt{B} "                                                                       
        elif (kwargs.get("significanceDef") == "S2"):                                                     
        ytitle = "S/B"                                
        '''
        
        # Customize the Significance Plot                                                                 
        #styleDict[signal_dataset.getName()].apply(tGraph)                                                 
        tGraph.SetName(datasetsMgr.getAllDatasetNames()[index])                                                          
    #tGraph.GetYaxis().SetTitle("\mathcal{S} =(S/sqrt(S+B))")                                         
    #tGraph.GetYaxis().SetTitle( ytitle +" (%s) / %s" % (cutDirection, GetBinwidthDecimals(binWidth) % (binWidth) ))                      
        #tGraph.GetYaxis().SetTitle( ytitle +" (%s) / %s" % (cutDirection, GetBinwidthDecimals(binWidth) % (binWidth) ))                                        
    #tGraph.GetYaxis().SetTitle( ytitle +" (%s) / %s" % (cutDirection,binWidth))                      
        #tGraph.GetXaxis().SetTitle(h_signal.GetXaxis().GetTitle())                                        
        optGraph = histograms.HistoGraph(tGraph, datasetsMgr.getAllDatasetNames()[index], "")
        histosSignif.append(optGraph)
        
    
    histoName = histoName +"_Significance"
    PlotHistograms(datasetsMgr,histoName,histosSignif, **kwargs_)

    return


def GetOptimizationPlot(datasetsMgr, histoName):

    #print "--plotAux.py: ploting Optimization for histogram", histoName                               

    # Get the options                                                                                  
    #cutDirection = kwargs.get("cutDirection")
    
    
    # Get the Cumulative Plots for Signal and Background                                               
    Signal_Cumul     = GetCumulativePlot(signal_datasets, histoName, **kwargs)
    Background_Cumul = GetCumulativePlot(background_dataset, histoName, **kwargs)


    ######
        # Create lists to append the x-y values of the Cumulative Histogram                               
    x_signal     = []                                                                                 
    y_signal     = []                                                                                 
    x_background = []                                                                                 
    y_background = []                                                                                 
                                                                                                      
    # Clone the Cumulative plots for Signal and Bkg                                                   
    h_signal     = Signal_Cumul.Clone()                                                               
    h_background = Background_Cumul.Clone()                                                           
                                                                                                      
    #Get the binNo of the Cumulative plots for Signal and Bkg                                         
    n_signal = h_signal.GetNbinsX()                                                                   
    n_background = h_background.GetNbinsX()                                                           
                                                                                                      
    # Get the bin-width of the plots                                                                  
    binWidth = h_signal.GetBinWidth(0)                                                                
                                                                                                      
    # Sanity Check                                                                                    
    if (n_signal != n_background):                                                                    
        print "Warning: Number of Bins of signal and background are different!"                      


        # Get the x-y values of the Cumulative Histogram                                                  
    for i in range(1, n_signal+1):                                                                    
                                                                                                      
        # Signal (x,y) Values                                                                         
        x_signal.append(h_signal.GetBinLowEdge(i)+0.5*h_signal.GetBinWidth(i))                        
        y_signal.append(h_signal.GetBinContent(i))                                                    
                                                                                                      
        # Background (x,y) values and errors                                                          
        x_background.append(h_background.GetBinLowEdge(i)+0.5*h_background.GetBinWidth(i))            
        y_background.append(h_background.GetBinContent(i))                                            
                                                                                                      
    # Create the (x-y) list to append the significance value for each x value                         
    x = []                                                                                            
    y = []                                                                                            
                                                                                                      
    # Create variables for Maximum Significance                                                       
    maxSignX = 0                                                                                      
    maxSignY = 0                                                                                      
                                                                                                      
    # Calculate the Significance of the Signal                                                        
    for i in range(0, n_signal):                                                                      
        if ((float(y_background[i]) <= 0 ) or (float(y_signal[i]) <= 0)):                             
            significance = 0                                                                          
        elif (kwargs.get("significanceDef") == "S1"):                                                 
            significance = float(y_signal[i])/math.sqrt(float(y_background[i]))                       
        elif (kwargs.get("significanceDef") == "S2"):                                                 
            significance = float(y_signal[i]) / float(y_background[i])                                

        # Create the Significance Plot                                                                    
    tGraph = ROOT.TGraph(n_signal, array.array("d", x), array.array("d", y))                          
                                                                                                      
    # Customize y-axis title                                                                          
    if (kwargs.get("significanceDef") == "S1"):                                                       
        ytitle = "S/ #sqrt{B} "                                                                       
    elif (kwargs.get("significanceDef") == "S2"):                                                     
        ytitle = "S/B"                                

    # Customize the Significance Plot                                                                 
    styleDict[signal_dataset.getName()].apply(tGraph)                                                 
    tGraph.SetName(signal_dataset.getName())                                                          
    #tGraph.GetYaxis().SetTitle("\mathcal{S} =(S/sqrt(S+B))")                                         
    #tGraph.GetYaxis().SetTitle( ytitle +" (%s) / %s" % (cutDirection, GetBinwidthDecimals(binWidth) % (binWidth) ))                      
    tGraph.GetYaxis().SetTitle( ytitle +" (%s) / %s" % (cutDirection, GetBinwidthDecimals(binWidth) % (binWidth) ))                                        
    #tGraph.GetYaxis().SetTitle( ytitle +" (%s) / %s" % (cutDirection,binWidth))                      
    tGraph.GetXaxis().SetTitle(h_signal.GetXaxis().GetTitle())                                        
    optGraph = histograms.HistoGraph(tGraph, legName, legStyle, drawStyle)      


    return optGraph, maxSignX         



def GetHistoKwargs(h, opts):
    _moveLegend = {"dx": -0.10, "dy": -0.55, "dh": -0.15}
    logY    = False
    _yLabel = "Arbitrary Units / %.2f "
    yMin    = 0.0
    if logY:
        yMaxF = 10
    else:
        yMaxF = 1.2
        
    _kwargs = {
        "ylabel"           : _yLabel,
        "rebinX"           : None,
        "rebinY"           : None,
        "ratioYlabel"      : "Ratio",
        "ratio"            : opts.ratio,
        "stackMCHistograms": opts.nostack,
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
        "xtitlesize"       : 0.1,#xlabelSize,
        "ytitlesize"       : 0.1,#ylabelSize,
        }

    kwargs = copy.deepcopy(_kwargs)

    if "leadtrks_pt" in h.lower():
        kwargs["opts"] = {"xmin": 0.0, "xmax": 80.0, "ymin": yMin, "ymaxfactor": yMaxF}
    

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


def getHistos(datasetsMgr, histoName):

    h1 = datasetsMgr.getDataset("Data").getDatasetRootHisto(histoName)
    h1.setName("Data")

    h2 = datasetsMgr.getDataset("EWK").getDatasetRootHisto(histoName)
    h2.setName("EWK")
    return [h1, h2]

def PlotHistograms(datasetsMgr, histoName, histos, **kwargs_):
    Verbose("Plotting Histograms")

    # Set the saveName
    saveName = histoName# +"_Cumulative"
    legendDict = {} 

    # Create & draw the plot                                                                   
    p = plots.ComparisonManyPlot(histos[0], histos[1:], saveFormats=[])

    # Set individual styles                                                                    
    for index, h in enumerate(p.histoMgr.getHistos()):
        p.histoMgr.setHistoDrawStyle(h.getName(), "PX")            
        p.histoMgr.setHistoLegendStyle(h.getName(), "LP")
        legendDict[h.getName()] = plots._legendLabels[h.getName()]
        plots._plotStyles[h.getName()].apply(p.histoMgr.getHistos()[index].getRootHisto())
                
    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany(legendDict)

    # Draw and save the plot                                                                   
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary                                            

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
    SAVEDIR      = "/afs/cern.ch/user/m/mtoumazo/public/html/hltaus/TkCalo/Optimisation/" #os.getcwd()
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
        raw_input("=== plotOptimisation.py: Press any key to quit ROOT ...")

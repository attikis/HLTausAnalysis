#!/usr/bin/env python
'''

Usage (single plot):
./plotEff.py -m <multicrab_directory> <jsonfile> [opts]

Usage (multiple plots):
./plotEff.py -m <pseudo_mcrab_directory> json/*.json
./plotEff.py -m <pseudo_mcrab_directory> json/*.json json/L1TkTau/*.json

Last Used:
./plotEff.py -m multicrab_CaloPlusTracks_v61XSLHC6_20170420T1537/ json/L1TkTau/Multiplicity.json -i "VBF|Neutrino"
./plotEff.py -m multicrab_CaloPlusTracks_v61XSLHC6_20170420T1537/ json/*/*.json -i "VBF|Neutrino"
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

import HLTausAnalysis.NtupleAnalysis.tools.dataset as dataset
import HLTausAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HLTausAnalysis.NtupleAnalysis.tools.styles as styles
import HLTausAnalysis.NtupleAnalysis.tools.plots as plots
import HLTausAnalysis.NtupleAnalysis.tools.histograms as histograms
import HLTausAnalysis.NtupleAnalysis.tools.aux as aux

import ROOT

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


def GetDatasetsFromDir(opts, json):
    Verbose("Getting datasets")
    
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=json["dataEra"],
                                                        searchMode=None, 
                                                        includeOnlyTasks="|".join(json["samples"]),
                                                        analysisName=json["analysis"])
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=json["dataEra"],
                                                        searchMode=None,
                                                        analysisName=json["analysis"],
                                                        includeOnlyTasks=opts.includeOnlyTasks)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=json["dataEra"],
                                                        searchMode=None,
                                                        analysisName=json["analysis"],
                                                        excludeTasks=opts.excludeTasks)
    else:
        raise Exception("This should never be reached")
    return datasets
    
    
def Plot(jsonfile, opts):
    Verbose("Plotting")

    with open(os.path.abspath(jsonfile)) as jfile:
        j = json.load(jfile)

        Verbose("Plotting %s. Will save under \"%s\"" % (j["title"], j["saveDir"]), True)

        # Setup the style
        style = tdrstyle.TDRStyle()
        style.setGridX(j["gridX"]=="True")
        style.setGridY(j["gridY"]=="True")
    
        # Set ROOT batch mode boolean
        ROOT.gROOT.SetBatch(opts.batchMode)

        # Setup & configure the dataset manager
        datasetsMgr = GetDatasetsFromDir(opts, j)
        if 0:
            datasetsMgr.loadLuminosities()
        #datasetsMgr.updateNAllEventsToPUWeighted()

        # Print information
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()

        # Set/Overwrite cross-sections
        for d in datasetsMgr.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0)

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py)    
        plots.mergeRenameReorderForDataMC(datasetsMgr)

        # Print dataset information
        datasetsMgr.PrintInfo()

        # Plot the histogram
        ComparisonPlot(datasetsMgr, j)
        return


def getHistos(datasetsMgr, datasetName, name1, name2):

    h1 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name1)
    h1.setName("h1" + "-" + datasetName)

    h2 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name2)
    h2.setName("h2" + "-" + datasetName)
    return [h1, h2]

def getHisto(datasetsMgr, datasetName, histoName):

    h1 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(histoName)
    h1.setName("h1" + "-" + datasetName)
    return [h1]


def ComparisonPlot(datasetsMgr, json):
    Verbose("Creating MC plot")
        
    # Create the MC Plot with selected normalization ("normalizeToOne", "normalizeByCrossSection", "normalizeToLumi")
    kwargs     = {}
    ylabel_    = json["ylabel"]
    normToOne_ = json["normalizationToOne"]=="True"
    if normToOne_:
        ylabel_ = ylabel_.replace(json["ylabel"].split(" /")[0], "Arbitrary Units")
    #p = plots.PlotSameBase(datasetsMgr, json["histogram"], normalizeToOne=normToOne_, **kwargs) #works
    #p = plots.ComparisonPlot(*getHistos(datasetsMgr, json['samples'][0], json['histogram'], json['histogram'])) #works
    p = plots.ComparisonPlot(*getHistos(datasetsMgr, json['samples'][0], json['histogram'], "Tk_Eff"))
    
    # Customise styling
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetFillStyle(0))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerSize(1.2))
    p.histoMgr.setHistoDrawStyleAll(json["drawStyle"])
    p.histoMgr.setHistoLegendStyleAll(json["legendStyle"])

    # Apply styles
    ##########################################
    p.histoMgr.forHisto("h1-VBF_HToTauTau" , styles.getBaselineStyle() )
    p.histoMgr.forHisto("h2-VBF_HToTauTau" , styles.getInvertedStyle() )

    # Set draw style
    p.histoMgr.setHistoDrawStyle("h1-VBF_HToTauTau", "LP")
    p.histoMgr.setHistoLegendStyle("h2-VBF_HToTauTau", "LP")
    # p.histoMgr.setHistoLegendStyleAll("LP")

    # Set legend labels
    p.histoMgr.setHistoLegendLabelMany({
            "h1-VBF_HToTauTau" : "h1",
            "h2-VBF_HToTauTau" : "h2",
            })
    ##########################################
    
    # Label size (optional. Commonly Used in counters)
    xlabelSize = None
    if "xlabelsize" in json:
        xlabelSize = json["xlabelsize"]
    ylabelSize = None
    if "ylabelsize" in json:
        ylabelSize = json["ylabelsize"]
    # Draw a customised plot
    saveName = os.path.join(json["saveDir"], json["title"])    
    
    plots.drawPlot(p, 
                   saveName,                  
                   xlabel            = json["xlabel"], 
                   ylabel            = ylabel_,
                   rebinX            = json["rebinX"],
                   stackMCHistograms = json["stackMCHistograms"]=="True", 
                   addCmsText        = json["addCmsText"]=="True",
                   cmsExtraText      = json["cmsExtraText"],
                   opts              = json["opts"],
                   log               = json["logY"]=="True", 
                   moveLegend        = json["moveLegend"],
                   cutBox            = {"cutValue": json["cutValue"],
                                        "fillColor": json["cutFillColour"],
                                        "box": json["cutBox"]=="True",
                                        "line": json["cutLine"]=="True",
                                        "greaterThan": json["cutGreaterThan"]=="True"},
                   xlabelsize        = xlabelSize,
                   ylabelsize        = ylabelSize,
                   )
    
    # Remove legend?
    if json["removeLegend"] == "True":
        p.removeLegend()

    # Additional text
    histograms.addText(json["extraText"].get("x"), json["extraText"].get("y"), json["extraText"].get("text"), json["extraText"].get("size") )

    # Save in all formats chosen by user
    saveFormats = json["saveFormats"]
    for i, ext in enumerate(saveFormats):
        Print("%s" % saveName + ext, i==0)
    p.saveAs(saveName, formats=saveFormats)
    return


def main(opts):
    Verbose("main function")

    jsonFiles = []
    # For-loop: All system script arguments
    for arg in sys.argv[1:]:

        # Skip if not a json file
        if ".json" not in arg:
            continue

        # Sanity check - File exists
        if not os.path.exists(arg):
            Print("The JSON file \"%s\" does not seem to be a valid path.. Please check that the file exists. Exit" % (arg), True)
            sys.exit()

        # Load & append json file
        with open(os.path.abspath(arg)) as jsonFile:
            try:
                json.load(jsonFile)
                jsonFiles.append(arg)
            except ValueError, e:
                Print("Problem loading JSON file %s. Please check the file" % (arg))
                sys.exit()
            
    # Sanity check - At least 1 json file found
    if len(jsonFiles) == 0:
        Print("No JSON files found. Please read the script instructions. Exit", True)
        print __doc__
        sys.exit()    

    # For-loop: All json files
    for j in jsonFiles:
        Print("Processing JSON file \"%s\"" % (j), True)
        Plot(j, opts)
    return


#================================================================================================
# Main
#================================================================================================
if __name__ == "__main__":

    # Default Settings 
    global opts
    BATCHMODE = True
    VERBOSE   = False
    #INTLUMI   = 1.0

    parser = OptionParser(usage="Usage: %prog [options]" , add_help_option=False,conflict_handler="resolve")

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation  NOT generates a window) [default: %s]" % BATCHMODE)

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    #parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
    #                  help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")
    
    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if opts.mcrab == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        print __doc__
        sys.exit(0)

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotEff.py: Press any key to quit ROOT ...")

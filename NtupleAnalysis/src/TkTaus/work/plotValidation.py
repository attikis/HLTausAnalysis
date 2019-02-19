#!/usr/bin/env python
'''
DESCRIPTION:
Script for plotting Rates, Efficiency, 
Rate Vs Efficiency (ROC), and Turn-on curves
for 3 Phase-2 L1 Tau algorithms.


USAGE:
./plotValidation.py -m <multicrab_dir>


EXAMPLES:
./plotValidation.py -m <multicrab_dir>


LAST USED:
./plotValidation.py -m <multicrab_dir>

'''

#================================================================================================  
# Imports
#================================================================================================  
import os
import sys
import re
import array 
import copy
import getpass
from optparse import OptionParser

import ROOT

import HLTausAnalysis.NtupleAnalysis.tools.dataset as dataset
import HLTausAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HLTausAnalysis.NtupleAnalysis.tools.styles as styles
import HLTausAnalysis.NtupleAnalysis.tools.plots as plots
import HLTausAnalysis.NtupleAnalysis.tools.histograms as histograms
import HLTausAnalysis.NtupleAnalysis.tools.aux as aux
import HLTausAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

import ROOT
ROOT.gROOT.SetBatch(True)


#================================================================================================
# Shell Styles
#================================================================================================
ws = ShellStyles.WarningStyle()
es = ShellStyles.ErrorStyle()
hs = ShellStyles.HighlightStyle()
cs = ShellStyles.CaptionStyle()
ns = ShellStyles.NormalStyle()
rs = ShellStyles.ResultStyle()
ss = ShellStyles.SuccessStyle()
# s = ShellStyles.NoteStyle()
# s = ShellStyles.AltStyle()
# s = ShellStyles.NoteLabel()
# s = ShellStyles.WarningLabel()
# s = ShellStyles.ErrorLabel()
# s = ShellStyles.HighlightAltStyle()
# s = ShellStyles.SuccessLabel()


#================================================================================================  
# Class definition
#================================================================================================  
def Verbose(msg, printHeader=False):
    '''
    Calls Print() only if verbose options is set to true.
    '''
    if not opts.verbose:
        return
    Print(msg, printHeader)
    return


def Print(msg, printHeader=True):
    '''
    Simple print function. If verbose option is enabled prints, otherwise does nothing.
    '''
    fName = __file__.split("/")[-1]
    if printHeader:
        print "=== ", fName
    print "\t", msg
    return


def getDirectoryContent(rootFile, directory, predicate=None):
    '''
    Get the directory content of a given directory in the ROOT file.
    '''
    return aux.listDirectoryContent(rootFile.Get(directory))
    
class Category:
    def __init__(self, name, opts):
        xMin  = 0
        yMin  = 0.0
        yMaxF = 1.2

        self.name       = name
        self.gOpts      = opts
        self.histonames = []
        self.labels     = {}
        self.colors     = {}
        self.histograms = {}
        #self.opts       = {"xmin": xMin, "xmax" : self.gOpts.xMax, "ymin" : yMin, "ymaxfactor": yMaxF}
        #self.opts2      = {"ymin": 0.3, "ymax": 1.7}
        self.moveLegend = {}
        return

    def getName(self):
        return self.name

    def addData(self,datarootfile):
        self.datafile = datarootfile
        return

    def addHisto(self, histo, legendlabel, color=0):
        Print("Appending histogram \"%s\" with legend label \"%s\"." % (histo, legendlabel), False)
        self.histonames.append(histo)
        self.labels[histo] = legendlabel
        self.colors[histo] = color
        return

    def setMoveLegend(self,move):
        self.moveLegend = move
        return

    def clone(self,name):
        returnCat = Category(name, self.gOpts)
        returnCat.histonames = self.histonames
        returnCat.labels     = self.labels
        returnCat.colors     = self.colors
        returnCat.opts       = self.opts
        returnCat.opts2      = self.opts2
        returnCat.moveLegend = self.moveLegend
        return returnCat

    def fetchHistograms(self, rootFiles, opts):

        hnamesList = copy.deepcopy(self.histonames)
        
        # For-loop: All histogram names
        for i, hName in enumerate(hnamesList, 1):
            # For-loop: All ROOT files
            for j, f in enumerate(rootFiles, 1):
                Verbose("Getting histogram \"%s\" from ROOT file \"%s\"" % (hName, f.GetName()), i==1)

                # Get the histograms
                histo = f.Get(hName)
                hType = str(type(histo))
                
                # Safety net for empty histos (if datasset integral is zero corresponding histo will not exist)
#                if "TH1" in hType:
#                    Verbose("Histogram is of type %s" % (hs + hType + ns), j==1)
#                else: 
#                    msg = "Skipping %s. Not a TH1 type histogram!" % (hName)
#                    Print(es + msg + ns, False)
#                    # Remove the first item from the list whose value is hName. It is an error if there is no such item.
#                    if hName in self.histonames:
#                        self.histonames.remove(hName)
#                    continue

                # Customise histogram
                histo.SetFillColor(self.colors[hName])
                histo.SetLineWidth(3)
                Print("Mapping histogram name (key) %s to histogram with name %s (object)" % (hName, histo.GetName()), i==1 and j==1) #iro - fixme (dataset?)
                self.histograms[hName] = histo # iro - fixme (dataset?)
        return


    def plot(self):
        style = tdrstyle.TDRStyle()
        ROOT.gStyle.SetErrorX(0.5) #required for x-axis error bars! (must be called AFTER tdrstyle.TDRStyle())

        histolist = []        
        # For-loop: All bkg histos
        for i, hName in enumerate(self.histonames, 1):
            myLabel = self.labels[hName]
            myHisto = self.histograms[hName]
            # myHisto = self.histograms[hName.split("/")[-1]] # iro - fixme
            
            Print("Creating histogram %s (Integral = %.1f)" % (hs + hName + ns, myHisto.Integral()), i==1) #iro - fixme
            hh = histograms.Histo(myHisto, hName, legendStyle="F", drawStyle="HIST", legendLabel=myLabel)
            if 0:
                aux.PrintTH1Info(hh.getRootHisto())
            hh.setIsDataMC(isData=False, isMC=True)
            histolist.append(hh)
            
        # Sanity check
        for i, h in enumerate(histolist, 1):
            hName = h.getRootHisto().GetName()
            Verbose("Histogram name is %s" % (hs + hName + ns), i==1)

        # Do the plot
        Verbose("Plotting all histograms for comparison", True)
        p = plots.ComparisonManyPlot(histolist[0], histolist[1:], saveFormats=[])
        p.setDefaultStyles()
        if "Vs" in str(type(histolist[0])):
            p.addMCUncertainty()
        p.setLegendHeader(self.getName())
                
        # Draw the plot
        saveName = histolist[0].getRootHisto().GetName()
        if not os.path.exists(opts.saveDir):
            os.makedirs(opts.saveDir)
        plots.drawPlot(p, os.path.join(opts.saveDir, saveName), **GetHistoKwargs(saveName, self))
        #plots.drawPlot(p, os.path.join(opts.saveDir, "alexandros"), **myParams)

        # Save the plot (not needed - drawPlot saves the canvas already)
        SavePlot(p, saveName, opts.saveDir, saveFormats = opts.saveFormats)
        return


def GetHistoKwargs(saveName, category):
    sName = saveName.lower()
    
    # Define variables
    xlabel  = "x"
    units   = "" #(GeV/c^{2})
    logx    = False
    logy    = False
    ymin    = 0.0
    ymaxf   = 1.2
    xmin    = 0.0
    xmax    = 1.0
    if "numberof" in sName:
        xlabel_ = "multiplicity"
        units   = ""
        xmin    = 0
        xmax    = 20
    elif "eta" in sName:
        xlabel_ = "#eta"
        units   = ""
        xmin    = -2.5
        xmax    = +2.5
    elif "et" in sName:
        xlabel_ = "E_{T}"
        units   = "GeV"
        xmin    =   0
        xmax    = 200
    elif "phi" in sName:
        xlabel_ = "#phi"
        units   = "rads"
        units   = ""
        xmin    = -3.15
        xmax    = +3.15       
    else:
        pass

    ylabel = "Entries (%s)" % (units)
    if units !=  "":
        xlabel = "%s (%s)" % (xlabel_, units)
    else:
        xlabel = xlabel_

    if "EtGenVsEt" in saveName:
        units  = "GeV"
        xlabel = "aaE_{T} (%s)" % (units)
        ylabel = "E_{T} (%s)" % (units)
        ymin   = 0 
        ymax   = 200
        xmin   = 0 
        xmax   = 200

    myParams = {}
    myParams["xlabel"]            = xlabel
    myParams["ylabel"]            = ylabel
    myParams["ratio"]             = False 
    myParams["ratioYlabel"]       = "Ratio "
    myParams["logx"]              = logx
    myParams["log"]               = logy
    myParams["ratioType"]         = "errorScale"
    myParams["ratioErrorOptions"] = {"numeratorStatSyst": False, "denominatorStatSyst": True}
    myParams["divideByBinWidth"]  = False
    myParams["errorBarsX"]        = True
    myParams["xlabelsize"]        = 25
    myParams["ylabelsize"]        = 25
    myParams["addMCUncertainty"]  = True
    myParams["addLuminosityText"] = False
    myParams["opts"]              = {"xmin": xmin, "xmax" : xmax, "ymin" : ymin, "ymaxfactor": ymaxf}
    myParams["opts2"]             = {"ymin": 0.3, "ymax": 1.7} # if ratio is enabled
    myParams["moveLegend"]        = category.moveLegend
    myParams["cmsExtraText"]      = "Phase-2 Simulation"
    myParams["cmsTextPosition"]   = "outframe"
    return myParams
                       
def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )

     # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        Print(saveNameURL, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return


def OpenRootFiles(opts):
    Verbose("Opening all ROOT files from multicrab directory \"%s\"" % (opts.mcrab), True)

    myRootFiles = []
    # For-loop: All datasets 
    for i, dset in enumerate(opts.dsetDict, 1):

        # For-loop: All ROOT files (for given dataset)
        for j, f in enumerate(opts.dsetDict[dset], 1):
            if os.path.isfile(f):
                Print("Opening file %s" % (f), i==1)
                myRootFiles.append(ROOT.TFile.Open(f, "R"))
            else:
                raise Exception("File \"%s\" not found." % (f) )

    Print("%sSuccessfully opened %d ROOT files from multicrab directory %s.%s" % (ss, len(myRootFiles), opts.mcrab, ns), True)
    return myRootFiles


def CloseRootFiles(myRootFiles):
    Verbose("Closing %d ROOT files." % ( len(myRootFiles) ), True)

    # For-loop: All ROOT files    
    for i, f in enumerate(myRootFiles, 1):
        Verbose("Closing file \"%s\"" % (f.GetName()), i==1)
        f.Close()
    Print("%sSuccessfully closed %d ROOT files.%s" % (ss, len(myRootFiles), ns), True)        
    return

def GetFoldersInRootFile(rootFile):
    folders = []
    fKeys   = rootFile.GetListOfKeys()
    # For-loop: All keys
    for k in fKeys:
        folders.append(k.GetName())
    return folders
    
def GetAlgorithmLabel(algo):
    if algo == "CaloTk":
        return "#it{tracks+calo}"
    elif algo == "TkTau":
        return "#it{tracks-only}"
    elif algo == "TkEG":
        return "#it{tracks+e#gamma}"
    else:
        raise Exception("%sUnknown algorithm %s" % (es, algo + ns))

    
def GetAlgorithmNames(rootFile):
    folders     = GetFoldersInRootFile(rootFile)
    algosTmp    = [ f.replace("Eff", "",).replace("Rate", "") for f in folders ]
    algos       = list(set(algosTmp)) #remove duplicate entries
    return algos


def main(opts):

    # Open all ROOT file (Nested for-loops)
    myRootFiles = OpenRootFiles(opts)
    hnameList   = getDirectoryContent(myRootFiles[0], opts.folder, predicate=None)
    hpathList   = [os.path.join(opts.folder, h) for h in hnameList]
    algos       = GetAlgorithmNames(myRootFiles[0])
    cObjects    = []


    for hpath in hpathList:
        
        # List categories and add the histogram names
        for a in algos:
            cObject = Category(GetAlgorithmLabel(a), opts)

        # FIXME - shouldn't these be inside the loop?
        # Customise legend position and size
        cObject.setMoveLegend( {"dx": -0.08, "dy": -0.02, "dh": -0.05 } )

        # Add histogram
        cObject.addHisto(hpath, hpath, color=ROOT.kBlue-1) # iro - fixme

        # Create list of algorithms
        algorithms = [cObject] # iro - fixme

        Verbose("Fetch histograms for %d algorithms" % (len(algorithms)), True)
        # For-loop: All aglorithms
        for i, algo in enumerate(algorithms, 1):
            Verbose("Fetching histograms for algorithm %s" % (hs + algo.getName() + ns), i==1)
            algo.fetchHistograms(myRootFiles, opts)
    
        Verbose("Plotting histograms for all algorithms", True)
        # For-loop: All algorithms
        for i, a in enumerate(algorithms, True):
            Verbose("Plotting histograms for algorithm %s" % (a.getName()), i==1)
            a.plot() 

    # Close all ROOT file (Nested for-loops)
    CloseRootFiles(myRootFiles)
    return


if __name__=="__main__":

    # Default options
    VERBOSE         = False
    #SAVEDIR         = "/afs/cern.ch/user/%s/%s/public/html/FitDiagnostics" % (getpass.getuser()[0], getpass.getuser())
    SAVEDIR         = "tmp" 
    URL             = False
    GRIDX           = False
    GRIDY           = False
    ROOTFILE        = "L1TausPerformance_full.root"
    MCRAB           = None
    FORMATS         = [".png"] #[".png", ".pdf", ".C"]
    FOLDER          = "TkEGRate" #"TkEGEff"
    
    parser = OptionParser(usage="Usage: %prog [options]", add_help_option=True, conflict_handler="resolve")

    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE) )
    
    parser.add_option("--url", dest="url", action="store_true", default=URL,
                      help="Don't print the actual save path the plots are saved, but print the URL instead [default: %s]" % URL)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all plots will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--rootFile", dest="rootFile", type="string", default=ROOTFILE,
                      help="The input ROOT file [default: %s]" % ROOTFILE)

    parser.add_option("--saveFormats", dest="saveFormats", type="string", default=FORMATS,
                      help="List of formats to save the plots [default: %s]" % FORMATS)
    
    parser.add_option("--gridX", dest="gridX", default=GRIDX, action="store_true",
                      help="Enable the grid for the x-axis [default: %s]" % (GRIDX) )

    parser.add_option("--gridY", dest="gridY", default=GRIDY, action="store_true",
                      help="Enable the grid for the y-axis [default: %s]" % (GRIDY) )

    parser.add_option("-m", "--mcrab", dest="mcrab", default=MCRAB, 
                      help="Name of multicrab directory [default: %s]" % (MCRAB) )
    
    parser.add_option("--folder", dest="folder", default=FOLDER, 
                      help="Folder inside the ROOT file containing the histograms to be plotted [default: %s]" % (FOLDER) )

    (opts, args) = parser.parse_args()

    opts.rootFiles = []
    opts.datasets  = []
    opts.dsetDict  = {}    
    # For-loop: All directories under the multicrab
    for crab in os.listdir(opts.mcrab):
        dset = crab.replace("crab_", "")
        if not os.path.isdir(os.path.join(os.getcwd(), opts.mcrab, crab)):
            print crab
            continue
        else:
            opts.datasets.append(dset)
            Verbose("Dataset \"%s\"" % (crab), len(opts.datasets) == 1)
            
        # Store filenames
        allFiles  = os.listdir(os.path.join(opts.mcrab, crab, "results"))
        rootFiles = [os.path.join(opts.mcrab, crab, "results", f) for f in allFiles if ".root" in f]
        opts.dsetDict[dset] = rootFiles
    
        # Sanity check
        if len(rootFiles) < 1:
            raise Exception("Found no ROOT files for dataset %s" % (d))
            
        # For-loop: All ROOT files (for given dataset)
        for i, f in enumerate(rootFiles, 1):
            Verbose("ROOT file \"%s\"" % (f), i==0)

    # Sanity check
    if len(opts.dsetDict) < 1:
        raise Exception("Found no datasets in multicrab dir \"%s\"" % (opts.mcrab))
    else:
        for i, d in enumerate(opts.dsetDict, 1):
            Verbose(d, True)
            for j, f in enumerate(opts.dsetDict[d], 1):
                Verbose(f, j==0)
                
    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="")

    # Histograms to ignore because dataset yield is zero (hack)
    opts.empty = []

    # Main code
    main(opts)

#!/usr/bin/env python
'''
DESCRIPTION:
Script for plotting Rates, Efficiency, 
Rate Vs Efficiency (ROC), and Turn-on curves
for 3 Phase-2 L1 Tau algorithms.


USAGE:
./plotValidation.py -m <multicrab_dir>


EXAMPLES:
./plotValidation.py -m multicrab_L1TkEG_PhaseIIFall17D_MC


LAST USED:
./plotValidation.py -m multicrab_L1TkEG_PhaseIIFall17D_MC_test/ --gridX --gridY --cutY

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
ls = ShellStyles.HighlightStyle()
hs = ShellStyles.HighlightAltStyle()
cs = ShellStyles.CaptionStyle()
ns = ShellStyles.NormalStyle()
rs = ShellStyles.ResultStyle()
ss = ShellStyles.SuccessStyle()
# s = ShellStyles.NoteStyle()
# s = ShellStyles.AltStyle()
# s = ShellStyles.NoteLabel()
# s = ShellStyles.WarningLabel()
# s = ShellStyles.ErrorLabel()
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


def ConverToRateHisto(hh):
    if "EtThreshold" not in hh.getName():
        return
    
    rh = hh.getRootHisto()
    convertTokHz  =  1.0e-03 # 1 Hz -> 1kHz
    crossingRate  = 30.0e+06 # 30MHz
    normFactor    = (crossingRate * convertTokHz) / (rh.GetBinContent(1))
    # normFactor    = (crossingRate * convertTokHz) / (rh.Integral())
    # normFactor    = (crossingRate * convertTokHz) / (rh.GetEntries())
    rh.Scale(normFactor)
    return


class Algorithm:
    def __init__(self, name, opts):
        self.name          = name
        self.label         = self.GetAlgorithmLabel(name)
        self.gOpts         = opts
        self.histonames    = []
        self.histoPaths    = []
        self.labels        = {}
        self.colors        = {}
        self.histograms    = {}
        self.histosSkipped = []
        self.histosCreated = []
        return

    def GetAlgorithmLabel(self, algo):
        if algo == "CaloTk":
            return "#it{tracks+calo}"
        elif algo == "TkTau":
            return "#it{tracks-only}"
        elif algo == "TrkTau":
            return "#it{tracks-only}"
        elif algo == "TkEG":
            return "#it{tracks+e#gamma}"
        else:
            raise Exception("%sUnknown algorithm %s" % (es, algo + ns))

    def getName(self):
        return self.name

    def getHistosCreated(self):
        return self.histosCreated

    def getHistosSkipped(self):
        return self.histosSkipped

    def getNHistosCreated(self):
        return len(self.histosCreated)

    def getNHistosSkipped(self):
        return len(self.histosSkipped)

    def getNHistosTotal(self):
        return self.getNHistosCreated() + self.getNHistosSkipped()

    def getHistonames(self):
        return self.histonames

    def getHistoPaths(self):
        # Ensure list has no duplicates
        return list(set(self.histoPaths))
    
    def getLabel(self):
        return self.label

    def getIndex(self):
        algo = self.name
        if algo == "CaloTk":
            return 0
        elif algo == "TkTau":
            return 1
        elif algo == "TrkTau":
            return 1
        elif algo == "TkEG":
            return 2
        else:
            raise Exception("%sUnknown algorithm %s" % (es, algo + ns))

    def addHisto(self, histo, legendlabel, color=0):
        Verbose("Appending histogram \"%s\" with legend label \"%s\"." % (histo, legendlabel), False)
        self.histonames.append(histo)
        self.labels[histo] = legendlabel
        self.colors[histo] = color
        return

    def clone(self,name):
        returnCat = Algorithm(name, self.gOpts)
        returnCat.histonames = self.histonames
        returnCat.labels     = self.labels
        returnCat.colors     = self.colors
        returnCat.opts       = self.opts
        returnCat.opts2      = self.opts2
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
                if "TH" not in hType:
                    msg = "Skipping non-histogram object \"%s\" ( type = \"%s\")" % (hName, hType)
                    Verbose(es + msg + ns, True)
                    continue
                else:
                    Verbose("Histogram \"%s\" if of type \"%s\"" % (hName, hType), True)
                    
                # Customise histogram                
                histo.SetFillColor(self.colors[hName])
                histo.SetLineWidth(3)
                Verbose("Mapping histogram name (key) %s to histogram with name %s (object)" % (hName, histo.GetName()), i==1 and j==1) #iro - fixme (dataset?)
                print "self.histograms = ", self.histograms
                self.histograms[hName] = histo # iro - fixme (dataset?)
                self.histoPaths.append(hPath)
        return

    def fetchHistogram(self, hPath, rootFiles, opts):

        # For-loop: All ROOT files
        for j, f in enumerate(rootFiles, 1):
            Verbose("Getting histogram \"%s\" from ROOT file \"%s\"" % (hPath, f.GetName()), i==1)
            
            # Get the histograms
            histo = f.Get(hPath)
            hType = str(type(histo))
            if "TH" not in hType:
                msg = "Skipping non-histogram object \"%s\" ( type = \"%s\")" % (hPath, hType)
                Verbose(ls + msg + ns, True)
                continue
            else:
                Verbose("Histogram \"%s\" if of type \"%s\"" % (hPath, hType), True)
                
            # Customise histogram
            histo.SetFillColor(self.colors[hPath])
            histo.SetLineWidth(3)
            Verbose("Mapping histogram name (key) %s to histogram with name %s (object)" % (hPath, histo.GetName()), i==1 and j==1) #iro - fixme (dataset?)
            self.histograms[hPath] = histo
            self.histoPaths.append(hPath)
        return


    def plot(self, hPath, verbose):

        # Sanity check
        if hPath not in self.histograms:
            msg = "Histogram %s was not stored for plotting. Skipping!" % (es + hPath + ns)
            Print(msg, False)
            return
        
        hList   = []
        myLabel  = self.labels[hPath]
        myHisto  = self.histograms[hPath]
        integral = myHisto.Integral()
        
        # Sanity check (Histogram is not empty)
        if integral == 0.0:
            msg = "Skipping histogram %s (Integral = %.1f)" % (ls + hPath + ns, integral)
            Print(msg, False)
            self.histosSkipped.append(hPath)
            return
        else:
            msg = "Creating histogram %s" % (hs + hPath + ns)
            Print(msg, verbose)

        # Create histogram histo
        #hh = histograms.Histo(myHisto, hPath, legendStyle="F", drawStyle="HIST", legendLabel=myLabel)
        hh = histograms.Histo(myHisto, hPath, legendStyle="L", drawStyle="HIST", legendLabel=myLabel)

        # Convert to rate histogram (scale appropriately)
        ConverToRateHisto(hh)

        # Set histogram type (Data or MC)
        hh.setIsDataMC(isData=False, isMC=True)
        hList.append(hh)
        if 0:
            aux.PrintTH1Info(hh.getRootHisto())
            
        # Apply histogram style
        styles.styles[self.getIndex()].apply(hh.getRootHisto())
        self.histosCreated.append(hPath)            
            
        # Sanity check
        if len(hList) < 1:
            return

        # Get the keyword arguments for this specific histogram
        saveName = hPath.replace("/", "_") 
        kwargs   = GetHistoKwargs(saveName, self)
        
        # Make the plot
        Verbose("Plotting all histograms for comparison", True)
        p = plots.ComparisonManyPlot(hList[0], hList[1:], saveFormats=[])
         
        # Apply Univeral style?
        if 0:
            p.histoMgr.setHistoDrawStyleAll("HIST")
            p.histoMgr.setHistoLegendStyleAll("FLP")

        # Set legend labels
        p.histoMgr.setHistoLegendLabelMany({hPath: kwargs["legend"]})

        # Set Legend header?
        if 0:
            p.setLegendHeader(self.getLabel())

        # Add MC uncertainty (if not a TH2)
        if "Vs" in str(type(hList[0])):
            p.addMCUncertainty() # appears to do nothing..
            
        # Draw the plot
        plots.drawPlot(p, hPath, **kwargs)

        # Save the plot
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
    xlabel_ = ""
    ymax    = None
    cutBoxX = None
    cutBoxY = None
    
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
        ymin    = 1e0
        logy    = True
    elif "phi" in sName:
        xlabel_ = "#phi"
        units   = "rads"
        xmin    = -3.15
        xmax    = +3.15       
    else:
        pass

    if units !=  "":
        ylabel = "Entries (%s)" % (units)
        xlabel = "%s (%s)" % (xlabel_, units)
    else:
        ylabel = "Entries"
        xlabel = xlabel_
        
    if "EtGenVsEt" in saveName:
        units  = "GeV"
        xlabel = "E_{T} (%s)" % (units)
        ylabel = "E_{T} (%s)" % (units)
        ymin   = 0 
        ymax   = 200
        xmin   = 0 
        xmax   = 200

    # Rate histogram (
    if "EtThreshold" in saveName:
        ylabel  = "Rate (kHz)"
        units   = ""
        ymin    = 1e0
        ymax    = 5e4
        #ymin    = 1e-4
        #ymax    = 2
        xmin    =   0
        xmax    = 100        
        # 46 GeV threshold line (TkEG threshold for 100 Hz)
        if "TkEG" in saveName:
            cutBoxX = {"cutValue": 50.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        elif "TrkTau" in saveName:
            cutBoxX = {"cutValue": 42.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        elif "CaloTk" in saveName:
            cutBoxX = {"cutValue": 63.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        else:
            pass
        # 50 Hz reference line
        cutBoxY = {"cutValue": 50, "fillColor": 16, "box": False, "line": True, "greaterThan": True, "mainCanvas": True, "ratioCanvas": False}
        
    if ymax == None:
        opts1 = {"xmin": xmin, "xmax" : xmax, "ymin" : ymin, "ymaxfactor": ymaxf}
    else:
        opts1 = {"xmin": xmin, "xmax" : xmax, "ymin" : ymin, "ymax": ymax}

    
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
    myParams["opts"]              = opts1
    myParams["opts2"]             = {"ymin": 0.3, "ymax": 1.7} # if ratio is enabled
    myParams["moveLegend"]        = {"dx": -0.02, "dy": +0.02, "dh": -0.15 }
    myParams["cmsExtraText"]      = "Phase-2 Simulation"
    myParams["cmsTextPosition"]   = "outframe"
    myParams["legend"]            = category.getLabel()
    if opts.cutX:
        myParams["cutBox"]  = cutBoxX
    if opts.cutY:
        myParams["cutBoxY"] = cutBoxY
    return myParams
                       
def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )

     # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats, 0):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        Verbose(saveNameURL, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return


def OpenRootFiles(opts):
    Verbose("Opening all ROOT files from multicrab directory \"%s\"" % (opts.mcrab), True)

    myRootFiles = []
    # For-loop: All datasets
    Print("Found %d datasets:\n\t%s" % (len(opts.dsetDict), "\n\t".join(opts.dsetDict.keys())), True)
    for i, dset in enumerate(opts.dsetDict, 1):
        
        # For-loop: All ROOT files (for given dataset)
        for j, f in enumerate(opts.dsetDict[dset], 1):
            if os.path.isfile(f):
                Verbose("Opening file %s" % (f), i==1)
                myRootFiles.append(ROOT.TFile.Open(f, "R"))
            else:
                msg = "File \"%s\" not found." % (f)
                raise Exception(es + msg + ns)

    if len(myRootFiles) < 1:
        msg = "Could not find any ROOT files under multicrab directory %s" % (opts.mcrab)
        raise Exception(es + msg + ns)
    
    Verbose("%sSuccessfully opened %d ROOT files from multicrab directory %s" % (ss, len(myRootFiles), opts.mcrab + ns), True)
    return myRootFiles


def CloseRootFiles(myRootFiles):
    Verbose("Closing %d ROOT files." % ( len(myRootFiles) ), True)

    # For-loop: All ROOT files    
    for i, f in enumerate(myRootFiles, 1):
        Verbose("Closing file \"%s\"" % (f.GetName()), i==1)
        f.Close()
    Verbose("%sSuccessfully closed %d ROOT files.%s" % (ss, len(myRootFiles), ns), True)        
    return

def GetFoldersInRootFile(rootFile):
    folders = []
    fKeys   = rootFile.GetListOfKeys()
    # For-loop: All keys
    for k in fKeys:
        folders.append(k.GetName())
    return folders
    
    
def GetAlgorithmNames(rootFile):
    folders     = GetFoldersInRootFile(rootFile)
    algosTmp    = [ f.replace("Eff", "",).replace("Rate", "") for f in folders ]
    algos       = list(set(algosTmp)) #remove duplicate entries
    return algos


def GetHistoPaths(rootFile):
    hPathList  = []
    # For-loop: All folders in ROOT file
    for folder in opts.folders:
        hNameList  = []
        hNameList.extend( getDirectoryContent(rootFile, folder, predicate=None) )
        hPathList.extend( [os.path.join(folder, h) for h in hNameList] )

    # Remove duplicates (safety)
    hPathList = list(set(hPathList))
    Print("Found %d histos:\n\t%s" % (len(hPathList), "\n\t".join(hPathList)), True)
    
    # Sort list alphabetically before returning
    hPathList.sort() 

    return hPathList
    
    
def main(opts):

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    style.setOptStat(False)
    ROOT.gStyle.SetErrorX(0.5) #required for x-axis error bars! (must be called AFTER tdrstyle.TDRStyle())
        
    # Open all ROOT file (Nested for-loops)
    myRootFiles = OpenRootFiles(opts)

    # Declare a reference ROOT file
    rootFile = myRootFiles[0]
    
    # Get all folders inside these files
    opts.folders = GetFoldersInRootFile(rootFile)
    Print("Found %d folders:\n\t%s" % (len(opts.folders), "\n\t".join(opts.folders)), True)    

    # Get all the histogram names
    algorithms = []
    algoNames  = []
    hPathList  = GetHistoPaths(rootFile)
    
    # Get algorithms present in ROOT files
    algos = GetAlgorithmNames(myRootFiles[0])
    Print("Found %d algorithms:\n\t%s" % (len(algos), "\n\t".join(algos)), True)

    # For-loop: All histograms inside ROOT files
    #print "\n".join(hPathList)
    #sys.exit()

    for i, hPath in enumerate(hPathList, 1):
        # print "hPath = ", hPath
        algoName = hPath.split("/")[0].replace("Rate", "").replace("Eff", "")
        # algoNames.append(algoName)
        algo = Algorithm(algoName, opts)

        # Add histogram
        Verbose("Adding histogram \"%s\"" % (hPath), i==1)
        algo.addHisto(hPath, os.path.basename(hPath), i+2) # iro - fixme
                
        Verbose("Fetching histograms for algorithm %s" % (hs + algo.getName() + ns), j==1)
        #algo.fetchHistograms(myRootFiles, opts) # f.Get(histoPath)
        algo.fetchHistogram(hPath, myRootFiles, opts) # f.Get(histoPath)

        Verbose("Plotting histogram %s" % (hPath), True)
        algo.plot(hPath, i==1)        
        
        # Create list of algorithms
        algorithms.append(algo)
    
    # Close all ROOT file (Nested for-loops)
    CloseRootFiles(myRootFiles)

    # Inform user of location of files
    Print("%sPlots saved under director \"%s\"%s" % (ss, opts.saveDir, ns), True)

    return


if __name__=="__main__":

    # Default options
    VERBOSE         = False
    #SAVEDIR         = "/afs/cern.ch/user/%s/%s/public/html/FitDiagnostics" % (getpass.getuser()[0], getpass.getuser())
    SAVEDIR         = "tmp" 
    URL             = False
    CUTX            = False
    CUTY            = False
    GRIDX           = False
    GRIDY           = False
    ROOTFILE        = "L1TausPerformance_full.root"
    MCRAB           = None
    FORMATS         = [".png"] #[".png", ".pdf", ".C"]
    
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

    parser.add_option("--cutX", dest="cutX", default=CUTX, action="store_true",
                      help="Enable the cut-line for the x-axis [default: %s]" % (CUTX) )

    parser.add_option("--cutY", dest="cutY", default=CUTY, action="store_true",
                      help="Enable the cut-line for the y-axis [default: %s]" % (CUTY) )

    parser.add_option("-m", "--mcrab", dest="mcrab", default=MCRAB, 
                      help="Name of multicrab directory [default: %s]" % (MCRAB) )
    
    (opts, args) = parser.parse_args()

    opts.rootFiles = []
    opts.datasets  = []
    opts.dsetDict  = {}    
    # For-loop: All directories under the multicrab
    for crab in os.listdir(opts.mcrab):
        dset = crab.replace("crab_", "")
        if not os.path.isdir(os.path.join(os.getcwd(), opts.mcrab, crab)):
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

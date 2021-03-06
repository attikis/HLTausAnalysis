#!/usr/bin/env python
'''
DESCRIPTION:
Script for plotting Rates, Efficiency, 
Rate Vs Efficiency (ROC), and Turn-on curves
for 3 Phase-2 L1 Tau algorithms. The input used
should be a multicrab with the output of the EDAnalyzer
that was written specifically to test the EDProducers of
the 3 L1 Tau algorithms!


USAGE:
./plotValidation.py -m <multicrab_dir>


EXAMPLES:
./plotValidation.py -m multicrab_L1TkEG_PhaseIIFall17D_MC
./plotValidation.py -m multicrab_L1TkEG_PhaseIIFall17D_MC_test/ --gridX --gridY --cutY


LAST USED:
./plotValidation.py -m multicrab_Phase2L1Taus_PhaseIIFall17D_MC --gridX --gridY --cutY

'''

#================================================================================================  
# Imports
#================================================================================================  
import os
import sys
import re
import fnmatch
import array 
import copy
import subprocess
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


def MergeRootFiles(allFiles, opts):
    if any(opts.rootFile in f for f in allFiles):
        pass
    else:
        # Get the number of root files in the crab directory
        nFiles = len(allFiles)
        
        # Merge all histogram root files to one and add it to the list of all files
        dirName = os.path.join(opts.mcrab, crab, "results")
        mergedRootFileName = "%s/%s-%s.root" % (dirName, opts.rootFile, dset)
        Verbose("Merging all the histogram root files for dataset %s using hadd" % (dset), True)
	subprocess.check_output("hadd %s %s/*.root" % (mergedRootFileName, dirName), shell=True)  # example: hadd result.root file1.root file2.root ... filen.root
        allFiles.append(mergedRootFileName)

        # Open the merged root file and add a histogram containing the number of the root files which are now merged
        mergedRootFile = ROOT.TFile(mergedRootFileName,"UPDATE");
        histo  = ROOT.TH1D('nFiles','NumberOfFiles',1,0.0,1.0)
        histo.SetBinContent(1, nFiles)
        mergedRootFile.Write()
	mergedRootFile.Close()
        return
    
def GetDatasetLabel(dsetName):
    if "GluGluHToTauTau" in dsetName:
        return  plots._legendLabels["VBF_HToTauTau"]
    elif "SingleNeutrino" in dsetName:
        return  plots._legendLabels["MinBias"]
    elif "ChargedHiggs200" in dsetName:
        return plots._legendLabels["ChargedHiggs200_14TeV_L1TPU200"]
    else:
        raise Exception("%sUnknown dataset %s" % (es, dataset + ns))

def getDirectoryContent(rootFile, directory, predicate=None):
    '''
    Get the directory content of a given directory in the ROOT file.
    '''
    return aux.listDirectoryContent(rootFile.Get(directory))


def ConvertToEfficiencyHisto(hh, dset, algorithm):
    if "EtEfficiency" not in hh.getName():
        return

    opts.dsetDict[dset] = rootFiles
    
    rh = hh.getRootHisto()
    nFiles = opts.nFiles[dset]
    print "nFiles = ", nFiles #iro
    normFactor    = 1.0/float(nFiles)
    rh.Scale(normFactor)
    return

'''    
   # Get denominator histogram!
    if 1:
        return
    hPath = algorithm.getName() + "Eff/EtEfficiencyNorm"
    if hPath in algorithm.histograms or hPath == algorithm.histograms:
        hNorm = algorithm.histograms[hPath]
    else:
        msg = "Could not find histogram \"%s\"" % (hPath)
        raise Exception(es + msg + ns)
    rh.Scale(1.0/hNorm.GetBinContent(1))
    return
'''    
    
def ConvertToRateHisto(hh):
    if "EtThreshold" not in hh.getName():
        return
    
    rh = hh.getRootHisto()
    convertTokHz  =  1.0e-03 # 1 Hz -> 1kHz
    crossingRate  = 30.0e+06 # 30MHz
    normFactor    = (crossingRate * convertTokHz) / (rh.GetBinContent(1)) # correct
    # normFactor    = (crossingRate * convertTokHz) / (rh.Integral())   # wrong!
    # normFactor    = (crossingRate * convertTokHz) / (rh.GetEntries()) # wrong!
    rh.Scale(normFactor)
    return

def ConvertToTurnOnHisto(hh, algorithm):
    if "GenEtTurnOn" not in hh.getName():
        return

    # Get denominator histogram!
    hPath    = algorithm.getName() + "Eff/GenEtVis"
    fullPath = None
    
    # For-loop: All histograms
    for h in algorithm.getHistograms():
        if hPath in h:
            fullPath = h
            # print "fullPath = ", fullPath
        else:
            pass

    if fullPath == None:
        return
    hDenominator = algorithm.histograms[fullPath]
    hIntegral    = hDenominator.Integral()
    
    # Sanity check
    if hIntegral == 0.0:
        msg = "Cannot create turn-on curve. Histogram \"%s\" has integral %.2f" % (hIntegral)
        raise Exception(es + msg + ns)

    # Get the ROOT histogram
    rh = hh.getRootHisto()

    # Divide by proper denominator to get turn-on curve
    Verbose("Converting histogram \"%s\" to a turn-on curve" % (hh.getName()), True)
    rh.Divide(hDenominator)
    return


#================================================================================================
# Class definition
#================================================================================================  
class Algorithm:
    def __init__(self, name, opts):
        self.name          = name
        self.label         = self.GetAlgorithmLabel(name)
        self.gOpts         = opts
        self.histoNames    = []
        self.histoPaths    = []
        self.histograms    = {}
        self.histos        = {}
        self.histosSkipped = []
        self.histosCreated = []
        self.labels        = {}
        self.colors        = {}        
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

    def getHistograms(self):
        '''
        Return the ROOT histograms
        '''        
        return self.histograms

    def getHistos(self):
        '''
        Return the histos
        '''
        return self.histos
    
    def getHistoNames(self):
        return self.histoNames
    
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

    def clone(self,name):
        returnCat = Algorithm(name, self.gOpts)
        returnCat.histoNames = self.histoNames
        returnCat.histoPaths = self.histoPaths
        returnCat.labels     = self.labels
        returnCat.colors     = self.colors
        returnCat.opts       = self.opts
        returnCat.opts2      = self.opts2
        return returnCat

    def fetchHistogram(self, dset, hPath, rootFiles, opts):
        
        # For-loop: All ROOT files
        for j, f in enumerate(rootFiles, 1):
            Verbose("Getting histogram \"%s\" from ROOT file \"%s\"" % (hPath, f.GetName()), i==1)
            
            # Get the histograms
            histo = f.Get(hPath)
            hType = str(type(histo))
            if "TH" not in hType:
                msg = "Skipping non-histogram object \"%s\" ( type = \"%s\")" % (hPath, hType)
                Print(ls + msg + ns, True)
                continue
            else:
                Verbose("Histogram \"%s\" if of type \"%s\"" % (hPath, hType), True)

            # Ensure that dataset information is contained in the name to avoid duplication
            histoName = os.path.join(dset, hPath)
            histo.SetName(histoName)
            Verbose("Mapping histogram name (key) %s to histogram with name %s (object)" % (hPath, histo.GetName()), i==1 and j==1) #iro - fixme (dataset?)
            # self.histograms[hPath] = copy.deepcopy(histo) # required so that they are available after closing ROOT files
            self.histograms[histoName] = copy.deepcopy(histo) # required so that they are available after closing ROOT files
            self.histoPaths.append(hPath)
            self.histoNames.append(histoName)
        return


    def plot(self, dset, hPath):

        # Sanity check
        if hPath not in self.histoPaths:
            msg = "Histogram %s was not stored for plotting. Skipping!" % (es + hPath + ns)
            Print(msg, False)
            return
        
        hList    = []
        myLabel  = os.path.basename(hPath)
        fullPath = os.path.join(dset, hPath)
        myHisto  = self.histograms[fullPath]
        integral = myHisto.Integral()
        
        # Sanity check (Histogram is not empty)
        if integral == 0.0:
            msg = "Skipping histogram %s (Integral = %.1f)" % (ls + hPath + ns, integral)
            Verbose(msg, False)
            self.histosSkipped.append(hPath)
            #return ls
            return es
        else:
            msg = "Creating histogram %s" % (hs + hPath + ns)
            Verbose(msg, False)

        # Create histogram histo
        if "Vs" in hPath:
            # opts.style.setWide(onoff=True, percIncrease=0.15)
            hh = histograms.Histo(myHisto, hPath, legendStyle="L", drawStyle="COLZ", legendLabel=myLabel)
        elif "genetturnon" in hPath.lower():
            hh = histograms.Histo(myHisto, hPath, legendStyle="LP", drawStyle="AP", legendLabel=myLabel)
        else:
            hh = histograms.Histo(myHisto, hPath, legendStyle="L", drawStyle="HIST", legendLabel=myLabel)

        # Convert to rate histogram (scale appropriately)
        ConvertToRateHisto(hh)

        # Convert to turn-on histogram (EtTurnOn/EtGen)
        ConvertToTurnOnHisto(hh, self)

        # Convert to efficiency histo (requires fix due to use of hadd)
        ConvertToEfficiencyHisto(hh, dset, self)
            
        # Skip meaningless histos
        if "EtThreshold" in hh.getName() and "neutrino" not in dset.lower():
            return es
        
        # Set histogram type (Data or MC)
        hh.setIsDataMC(isData=False, isMC=True)
        hList.append(hh)
        if 0:
            aux.PrintTH1Info(hh.getRootHisto())
            
        # Apply histogram style
        styles.styles[self.getIndex()].apply(hh.getRootHisto())
        self.histosCreated.append(hPath)            

        # Replace original histo with final version of it
        self.histos[fullPath] = hh #iro
        if "GenEtTurnOn" in hh.getName():
            pass
        # aux.PrintTH1Info(hh.getRootHisto())
            #sys.exit()

        # Sanity check
        if len(hList) < 1:
            return es

        # Get the keyword arguments for this specific histogram
        # saveName = hPath.replace("/", "_").replace(self.getName(), self.getName() + "_") #iro
        saveName = hPath.replace("/", "_")
        kwargs   = GetHistoKwargs(saveName, self.getLabel())
        
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
        p.setLegendHeader(GetDatasetLabel(dset)) # p.setLegendHeader(self.getLabel())

        # Add MC uncertainty (if not a TH2)
        if "Vs" in str(type(hList[0])):
            p.addMCUncertainty() # appears to do nothing..
            
        # Draw the plot
        plots.drawPlot(p, dset + ":" + hPath, **kwargs)

        # Save the plot
        SavePlot(p, saveName, os.path.join(opts.saveDir, dset), saveFormats = opts.saveFormats)
        return ss


def GetHistoKwargs(saveName, label):
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
    legend  = {"dx": -0.02, "dy": +0.02, "dh": -0.15 }
    format_ = "%.0f"
            
    if "numberof" in sName or "multiplicity" in sName:
        xlabel_ = "multiplicity"
        units   = ""
        format_ = "%.0f"
        xmin    = 0
        xmax    = 20
    elif "eta" in sName:
        xlabel_ = "#eta"
        format_ = "%.2f"
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
        format_ = "%.2f"
        xmin    = -3.15
        xmax    = +3.15       
    else:
        pass

    if units !=  "":
        ylabel = "Entries / %s %s" % (format_, units)
        xlabel = "%s (%s)" % (xlabel_, units)
    else:
        ylabel = "Entries / %s" % (format_)
        xlabel = xlabel_
        
    if "GenEtVsEt" in saveName or "EtVsGenEt" in saveName:
        units  = "GeV"
        ylabel = "E_{T}^{GEN} (%s)" % (units)
        xlabel = "E_{T}^{L1}  (%s)" % (units)
        ymin   = 0 
        ymax   = 200
        xmin   = 0 
        xmax   = 200
        logx   = False
        logy   = False
        legend = {"dx": -0.1, "dy": -0.62, "dh": -0.15 }        

    # Turn-on curves
    if "genetturnon" in sName:
        units  = "GeV"
        xlabel = "#tau_{h} E_{T}^{vis} (%s)" % (units) # VERIFY WITH MARINA/CODE
        ylabel = "Efficiency / %.0f " + units
        ymin   = 0 
        ymax   = 1.0
        xmin   = 0 
        xmax   = 200 #160
        logy   = False
        
    #  Efficiency histogram
    if "EtEfficiency" in saveName:
        units   = "GeV"
        ylabel  = "Efficiency / %.0f " + units
        ymin    = 0.0
        ymax    = 1.0
        xmin    =   0
        xmax    = 160
        logy   = False
        
    # Rate histogram
    if "EtThreshold" in saveName:
        ylabel  = "Rate (kHz) / %.0f GeV"
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
    myParams["moveLegend"]        = legend
    myParams["cmsExtraText"]      = "Phase-2 Simulation"
    myParams["cmsTextPosition"]   = "outframe"
    myParams["legend"]            = label
    # myParams["rebin"]             =  2
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

def GetFoldersInRootFiles(opts):
    folders  = []
    rDataset = opts.dsetDict[opts.dsetDict.keys()[0]]
    rootFile = rDataset[0]
    fKeys    = rootFile.GetListOfKeys()
    
    # For-loop: All keys
    for k in fKeys:
        if "nFiles" in k.GetName():
            continue
        folders.append(k.GetName())
    return folders
    
    
def GetAlgorithmNames(opts):
    '''
    Get the first dataset found in the dataset dictionary and retrieves the mapped
    list of ROOT files. From that list it selects the first one and then inspects
    the folders to see which EDProducer was called by the EDAnalyzer and saved information
    in the file. The saved folders (containing histograms) contain the EDProducer name 
    and can hence help us determine which algorithms are present in the ROOT files.

    Assumption: All datasets have ROOT files with the exact same composition
    '''
    allowed     = ["TrkTau", "CaloTk", "TkEG"]
    rDataset    = opts.dsetDict[opts.dsetDict.keys()[0]]
    rootFile    = rDataset[0]    
    # folders     = GetFoldersInRootFile(rootFile)
    algosTmp    = [ f.replace("Eff", "",).replace("Rate", "") for f in opts.folders ]
    algos       = list(set(algosTmp)) #remove duplicate entries

    # Finally impose restrictions of "allowed" algos
    algos       = [a for a in algos if a in allowed]
    return algos


def GetHistoPaths(rootFile, algoList):
    hPathList  = []
    # For-loop: All folders in ROOT file
    for folder in opts.folders:
        hNameList  = []
        hNameList.extend( getDirectoryContent(rootFile, folder, predicate=None) )
        hPathList.extend( [os.path.join(folder, h) for h in hNameList] )

    # Remove duplicates (safety)
    hPathList = list(set(hPathList))
    Verbose("Found %d histos:\n\t%s" % (len(hPathList), "\n\t".join(hPathList)), True)

    # Sort list alphabetically before returning
    hPathList.sort()

    # For-loop: All algorithms
    for algo in algoList:
        hPath = "%sEff/GenEtVis" % algo
        # Ensure that GenEt histos (turn-on ingredients curves) are done first!
        if hPath in hPathList:
            hPathList.insert(0, hPathList.pop(hPathList.index(hPath)))
        else:
            msg = "The histogram \"%s\" is missing. Cannot create turn-on curves!" % (hPath)
            Print(ls + msg + ns, True)
        
    return hPathList
    
def PrintRootFileNames(rootFiles):

    # For-loop: All ROOT file objects
    for i, f in enumerate(rootFiles, 1):
        Print("%d) = %s" % (i, f.GetName()), i==1)
    return

    
def main(opts):

    # Apply TDR style
    opts.style = tdrstyle.TDRStyle()
    opts.style.setGridX(opts.gridX)
    opts.style.setGridY(opts.gridY)
    # opts.style.setLogX(opts.logX)
    # opts.style.setLogY(opts.logY)
    opts.style.setOptStat(False)
    ROOT.gStyle.SetErrorX(0.5) #required for x-axis error bars! (overwrites tdrstyle.TDRStyle())    

    
    # Create  dictionary of algorithms (name -> object) for each dataset (?)
    algorithms = {}
    for algoName in opts.algos:
        algorithms[algoName] = Algorithm(algoName, opts)
        
    # Create dictionary mapping each dataset to the number of merged ROOT files
    opts.nFiles = {}
    
    # For-loop: All datasets
    for i, dset in enumerate(opts.dsetDict, 1):
        
        # Get the list with the ROOT files
        rootFiles = opts.dsetDict[dset]
        if opts.verbose:
            PrintRootFileNames(rootFiles)

         # Get the number of the merged ROOT files
        mergedRootFile = ROOT.TFile(rootFiles[0].GetName(), "OLD");
        if (mergedRootFile.GetListOfKeys().Contains("nFiles")):
            opts.nFiles[dset] = mergedRootFile.Get("nFiles").GetBinContent(1)
        else:
            msg = "Could not find the histogram with the number of the merged root files!"
            raise Exception(es + msg + ns)

        
        # Get all the histogram names
        hPathList = GetHistoPaths(rootFiles[0], opts.algos)

        # For-loop: All histograms inside the ROOT files
        for j, hPath in enumerate(hPathList, 1):

            aName = hPath.split("/")[0].replace("Rate", "").replace("Eff", "")
            algo  = algorithms[aName]
                
            Verbose("Fetching histograms for algorithm %s" % (hs + algo.getName() + ns), j==1)
            algo.fetchHistogram(dset, hPath, rootFiles, opts)

            Verbose("Plotting histogram %s" % (hPath), True)
            style = algo.plot(dset, hPath)
            msg   = "{:>10} {:>20} {:<30}".format("%s:" % (algo.getName()), "%s:" % dset, "%s" % (hPath))
            Print(style + msg + ns, j==1)

        # Close all ROOT files for given dataset
        CloseRootFiles(rootFiles)


    # Declare list of histos to be drawn
    turnOns = []
    rates   = []
    effVBF  = []     

    # For-loop: All final algorithm objects
    for i, algoName in enumerate(algorithms, 1):
        algo  = algorithms[algoName]

        hDict = algo.getHistos()
        for j, hName in enumerate(hDict, 1):
            h = hDict[hName]
            if "GenEtTurnOn" in h.getRootHisto().GetName():
                turnOns.append(h)
            elif "EtThreshold" in h.getRootHisto().GetName():
                if "Neutrino" not in h.getRootHisto().GetName():
                    pass
                else:
                    rates.append(h)
            elif "EtEfficiency" in h.getRootHisto().GetName():
                if "HToTauTau" in h.getRootHisto().GetName():
                    effVBF.append(h)
            else:
                pass

    # Make the plot
    Verbose("Plotting turn-on histograms for comparison", True)
    p = plots.ComparisonManyPlot(turnOns[0], turnOns[1:], saveFormats=[])
    kwargs = GetHistoKwargs("GenEtTurnOn", "TEST-1")
    plots.drawPlot(p, "TurnOns", **kwargs)
    SavePlot(p, "TurnOns", os.path.join(opts.saveDir, "Algos"), saveFormats = opts.saveFormats)

    # Make the plot
    Verbose("Plotting rate histograms for comparison", True)
    p = plots.ComparisonManyPlot(rates[0], rates[1:], saveFormats=[])
    kwargs = GetHistoKwargs("EtThreshold", "TEST-2")
    plots.drawPlot(p, "Rates", **kwargs)
    SavePlot(p, "Rates", os.path.join(opts.saveDir, "Algos"), saveFormats = opts.saveFormats)

    # Make the plot
    Verbose("Plotting efficienciesiciency histograms for comparison", True)
    p = plots.ComparisonManyPlot(effVBF[0], effVBF[1:], saveFormats=[])
    kwargs = GetHistoKwargs("EtEfficiency", "TEST-3")
    plots.drawPlot(p, "Efficiency_VBF", **kwargs)
    SavePlot(p, "Efficiency_VBF", os.path.join(opts.saveDir, "Algos"), saveFormats = opts.saveFormats)
    
    # Inform user of location of files
    Print("%sPlots saved under director \"%s\"%s" % (hs, opts.saveDir, ns), True)

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
    ROOTFILE        = "histograms"
    MCRAB           = None
    FORMATS         = [".png"] #[".png", ".pdf", ".C"]
    
    parser = OptionParser(usage="Usage: %prog [options]", add_help_option=True, conflict_handler="resolve")

    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE) )
    
    parser.add_option("--url", dest="url", action="store_true", default=URL,
                      help="Don't print the actual save path the plots are saved, but print the URL instead [default: %s]" % URL)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all plot=s will be saved [default: %s]" % SAVEDIR)

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
    opts.algos     = []
    
    # For-loop: All directories under the multicrab
    for i, crab in enumerate(os.listdir(opts.mcrab), 1):
        dset = crab.replace("crab_", "")
        if not os.path.isdir(os.path.join(os.getcwd(), opts.mcrab, crab)):
            continue
        else:
            opts.datasets.append(dset)
            Verbose("Dataset \"%s\"" % (crab), len(opts.datasets) == 1)

        # Store filenames
        allFiles      = os.listdir(os.path.join(opts.mcrab, crab, "results"))
        MergeRootFiles(allFiles, opts)
        rootFileNames = [os.path.join(opts.mcrab, crab, "results", f) for f in allFiles if opts.rootFile in f]
        rootFiles     = [ROOT.TFile.Open(f, "R") for f in rootFileNames if os.path.isfile(f)]
        opts.dsetDict[dset] = rootFiles
    
        # Sanity check
        if len(rootFiles) < 1:
            msg = "Found no ROOT files for dataset %s" % (dset)
            raise Exception(es + msg + ns)
            
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
                Print(f.GetName(), j==0)

    # Inform user (briefly)
    Verbose("Found %d datasets: %s" % (len(opts.datasets), ", ".join(opts.datasets)), True)

    # Get folders present in ROOT files
    opts.folders = GetFoldersInRootFiles(opts)
    Verbose("Found %d folders: %s" % (len(opts.folders), ", ".join(opts.folders)), True)    

    # Get algorithms present in ROOT files
    opts.algos = GetAlgorithmNames(opts)
    Verbose("Found %d algorithms: %s" % (len(opts.algos), ", ".join(opts.algos)), True)

    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="")

    # Do the plotting
    main(opts)

#!/usr/bin/env python
'''
DESCRIPTION:
This script is used to create a pseudo multicrab directory, using an input ROOT file (raw2TTree.root).
The purpose is primarily to enable the easy testing of local changes to the  code and the 
NtupleAnalysis code using the hplusGenerateDataFormats.py script.
A script execution will thus create an empty multicrab with identical name and structure as those created
by the multicrab.py script. It will contain a single dataset with a single ROOT file under results/ dir, which is 
a mere copy of the file used as input for the script execution (only renamed) to histograms-<dataset>.root.


USAGE:
hltausPseudoMulticrab.py -f test.root [opts]

To add  a dataset to an existing <multicrab_dir>:
hltausPseudoMulticrab.py -f test.root -r <multicrab_dir> 


EXAMPLES:
hltausPseudoMulticrab.py -f histograms-TT_TuneCUETP8M2T4_14TeV_powheg_pythia8_PhaseIIFall17D_L1TPU140_93X.root --dir multicrab_test
hltausPseudoMulticrab.py -f histograms-PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV_PhaseIIFall17D_L1TPU140_93X.root --dir multicrab_test
hltausPseudoMulticrab.py


LAST USED:
hltausPseudoMulticrab.py --dir multicrab_CaloTk_v92X_IsoConeRMax0p4_VtxIso1p0

'''

#================================================================================================
# Import modules
#================================================================================================
import os
import re
import sys
import time
import datetime
import subprocess
from optparse import OptionParser

import ROOT

import HLTausAnalysis.NtupleAnalysis.tools.git as git
import HLTausAnalysis.NtupleAnalysis.tools.aux as aux
from HLTausAnalysis.NtupleAnalysis.tools.datasets import *
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
# Function Definitions
#================================================================================================ 
def CheckRootFile(fileName):

    if not os.path.isfile(fileName):        
        raise Exception("The ROOT file provided (%s) does not exists!" % (opts.rootFile) )

    # Attempty to open ROOT file (determine if null, zombie, etc..)
    rf = ROOT.TFile.Open(fileName)
    if rf == None:
        raise Exception("The file %s is NULL" % (es + fileName + ns) )
    elif rf.IsZombie():
        raise Exception("The ROOT file %s is a Zombie" % (es + fileName + ns) )
    else:
        Verbose("The ROOT file %s is healthy" % (hs + fileName + ns) )
    return


def getDatasetInfo(opts):
    try:
        opts.datasetName = re.search('histograms-(.+?).root', opts.rootFile).group(1)
        Verbose("File \"%s\" corresponds to dataset is \"%s\""% (opts.rootFile, opts.datasetName), False)
    except AttributeError:
        raise Exception("Could not determine the dataset corresponding to the file \"%s\"" % (opts.rootFile) )

    # Get dataset info
    dgroup = DatasetGroup(opts.dataEra)
    opts.dataset = None
    opts.datasetAlias = None
    opts.dataVersion = None    
    opts.pileup = None    

    # For-loop: All datasets objects 
    for i,d in enumerate(dgroup.GetDatasetList(),1):
        Verbose("Attempting to find alias for dataset %s" % d.getName(), i==1)

        name  = opts.datasetName
        alias = d.getAlias()

        # Stupid "hack" for special case (one-off!)
        if "RelVal" in name:
            name  = opts.datasetName.replace("RelValSingleTauFlatPt2To100_pythia8", "SingleTau").replace("2023D17", "L1T")

        if 0:
            print "%s in %s?" % (alias.split("_"), name.split("_"))

        if set(alias.split("_")).issubset(name.split("_")):
            opts.dataset = d
            opts.datasetAlias = alias
            Verbose("Alias for dataset %s is %s" % (d.getName(), alias), True)
            break
        else:
            pass

    if opts.datasetAlias == None:
        raise Exception("%sCould not determine the alias of datasets \"%s\"%s. Check the file \"datasets.py\"" % (es, opts.datasetAlias, ns) )

    # Get data-version, <PU>, and lumi
    opts.dataVersion = d.getDataVersion()
    opts.pileup = d.getPU()
    if opts.pileup == 140:
        opts.lumi = 5 # 5e34 (<PU>=140)
    elif opts.pileup == 200:
        opts.lumi = 7 # 7e34 (<PU>=200)
    else:
        opts.lumi = 0 # 7e34 (<PU>=200)

    return

def writeConfigInfoToRootFile(fileName, opts):
    '''
    Creates configInfo folder where the following objects are written:
    - dataVersion (CMSSW version, followed by "mc" or "data". e.g. 80Xmc
    - codeVersion (git commit id)
    - configinfo (histogram with COM energy, luminosity, pileup, etc..)
    '''
    if os.path.isfile(fileName) == False:
        raise Exception("The file %s does not exist!" % (fileName) )

    # Definitions
    filePath = fileName
    fileMode = "UPDATE"

    # Open the ROOT file
    nAttempts   = 1
    maxAttempts = 3

    fOUT = GetRootFile(fileName, fileMode="UPDATE")

    # Write the configInfo
    fOUT.cd("/")
    myConfigInfoDir = fOUT.mkdir("configInfo")

    # Create configinfo histogram
    hConfigInfo = ROOT.TH1F("configinfo", "configinfo", 6, 0, 4)
    hConfigInfo.GetXaxis().SetBinLabel(1, "control")
    hConfigInfo.SetBinContent(1, 1)
    hConfigInfo.GetXaxis().SetBinLabel(2, "energy")
    hConfigInfo.SetBinContent(2, 14)
    hConfigInfo.GetXaxis().SetBinLabel(3, "isData")
    hConfigInfo.SetBinContent(3, 0)
    hConfigInfo.GetXaxis().SetBinLabel(3, "isPileupReweighted")
    hConfigInfo.SetBinContent(3, 0)
    hConfigInfo.GetXaxis().SetBinLabel(4, "isTopPtReweighted")
    hConfigInfo.SetBinContent(4, 0)
    hConfigInfo.GetXaxis().SetBinLabel(5,"luminosity")
    hConfigInfo.SetBinContent(5, opts.lumi*1e+34)
    hConfigInfo.GetXaxis().SetBinLabel(6,"pileup")
    hConfigInfo.SetBinContent(6, opts.pileup)
    
    # Write a copy of data version and code version
    hConfigInfo.SetDirectory(myConfigInfoDir)
    dataVersion = ROOT.TNamed("dataVersion", str(opts.dataVersion))
    codeVersion = ROOT.TNamed("codeVersionAnalysis", git.getCommitId())
    myConfigInfoDir.Add(dataVersion)
    myConfigInfoDir.Add(codeVersion)
    
    # Write and close the root file
    fOUT.cd("/")
    fOUT.Write()

    # Close the ROOT file
    Verbose("Closing file %s" % (fOUT.GetName()) )
    fOUT.Close()
    return


def writeCounters(fileName, opts):
    '''
    Creates configInfo folder where the following objects are written:
    - dataVersion (CMSSW version, followed by "mc" or "data". e.g. 80Xmc
    - codeVersion (git commit id)
    - configinfo (histogram with COM energy, luminosity, pileup, etc..)
    '''
    if os.path.isfile(fileName) == False:
        raise Exception("The file %s does not exist!" % (fileName) )

    # Definitions
    filePath = fileName
    fileMode = "UPDATE"

    # Open the ROOT file
    nAttempts   = 1
    maxAttempts = 3

    fOUT = GetRootFile(fileName, fileMode="UPDATE")
    
    # Create the  counters folder
    folderName = opts.analysisName + "_" + opts.dataEra
    hCounters = fOUT.Get(folderName + "/Counters")
    

    fOUT.cd(folderName)
    folderName += "/counters"    
    countersDir = fOUT.mkdir(folderName)
    fOUT.cd(folderName)
    # Write counter histogram
    hCounter = ROOT.TH1F("counter", "counter", 2, 0, 2)
    # hCounter.GetXaxis().SetBinLabel(?, "ttree:: skimCounterAll")
    # hCounter.GetXaxis().SetBinLabel(?, "ttree:: skimCounterPassed")
    hCounter.GetXaxis().SetBinLabel(1, "Base::AllEvents")
    hCounter.SetBinContent(1, hCounters.GetBinContent(1))
    hCounter.GetXaxis().SetBinLabel(2, "Base::Events")
    hCounter.SetBinContent(2, hCounters.GetBinContent(2))
        
    # Create the counters/weighted folder
    folderName += "/weighted"
    weightedCountersDir = fOUT.mkdir(folderName)
    fOUT.cd(folderName)
    # Write counter histogram
    hCounterWeighted = ROOT.TH1F("counter", "counter", 2, 0, 2)
    hCounterWeighted.GetXaxis().SetBinLabel(1, "Base::AllEvents")
    hCounterWeighted.SetBinContent(1, hCounters.GetBinContent(1))
    hCounterWeighted.GetXaxis().SetBinLabel(2, "Base::Events")
    hCounterWeighted.SetBinContent(2, hCounters.GetBinContent(2))

    # Write and close the root file
    fOUT.cd("/")
    fOUT.Write()
        
    # Close the ROOT file
    Verbose("Closing file %s" % (fOUT.GetName()) )
    fOUT.Close()
    return

def moveAllHistosIntoAnalysisFolder(fileName, opts):
    '''
    Moves all histograms into an analysis folder
    '''
    if os.path.isfile(fileName) == False:
        raise Exception("The file %s does not exist!" % (fileName) )

    # Go to home directory
    fOUT = GetRootFile(fileName, fileMode="UPDATE")

    # Get list of all histograms in ROOT file
    fOUT.cd("/")
    keyList = fOUT.GetListOfKeys()

    # Create an analysis folder & Move all histograms under it
    folderName = opts.analysisName + "_" + opts.dataEra
    fOUT.mkdir(folderName)
    fOUT.cd(folderName)

    # For-loop: All objects in ROOT ifle
    for k in keyList:
        # print k
        obj = k.ReadObj()
        obj.Write()

    # Get list of all histograms in ROOT file
    keyListNew = fOUT.GetListOfKeys()
    for i,k in enumerate(keyListNew):
        if folderName in k.GetName():
            pass
        elif "configInfo" in k.GetName():
            pass
        else:
            Verbose("Deleting key with name %s" % (k.GetName()), i==0)
            k.Delete()

    # Close the ROOT file
    Verbose("Closing file %s." % (fOUT.GetName()) )
    fOUT.Close()
    return

def GetRootFile(fileName, fileMode="UPDATE"):
    # Definitions
    filePath    = fileName
    fileMode    = "UPDATE"
    nAttempts   = 1
    maxAttempts = 3
    fOUT        = None

    # Attempt to open
    while nAttempts < maxAttempts:
        try:
            Verbose("Attempt #%s: Opening ROOT file %s in %s mode." % (nAttempts, filePath, fileMode) )
            fOUT = ROOT.TFile.Open(filePath, fileMode)
            fOUT.cd()
            break
        except:
            nAttempts += 1
            Print("TFile::Open(\"%s\", \"%s\") failed (%s/%s). Retrying..." % (filePath, fileMode, nAttempts, maxAttempts) )

    # Safety clause
    if fOUT == None:
        raise Exception("TFile::Open(\"%s\", \"%s\") failed" % (filePath, fileMode) )
    else:
        Verbose("Successfully opened %s in %s mode (after %s attempts)" % (filePath, fileMode, nAttempts) )
    return fOUT

def GetSelfName():
    return __file__.split("/")[-1]

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

def AskUser(msg):
    '''
    Prompts user for keyboard feedback to a certain question. 
    Returns true if keystroke is \"y\", false otherwise.
    '''
    keystroke = raw_input("\t" +  msg + " (y/n): ")
    if (keystroke.lower()) == "y":
        return True
    elif (keystroke.lower()) == "n":
        return False
    else:
        AskUser(msg)
    return False

def GetCMSSW(opts):
    '''
    Get a command-line-friendly format of the CMSSW version currently use.
    https://docs.python.org/2/howto/regex.html
    '''
    if opts.cmssw != "":
        return opts.cmssw
            
    # Get the current working directory
    pwd = os.getcwd()

    # Create a compiled regular expression object
    cmssw_re = re.compile("/CMSSW_(?P<version>\S+?)/")

    # Scan through the string 'pwd' & look for any location where the compiled RE 'cmssw_re' matches
    match = cmssw_re.search(pwd)

    # Return the string matched by the RE. Convert to desirable format
    version = ""
    if match:
	version = match.group("version")
	version = version.replace("_","")
	version = version.replace("pre","p")
	version = version.replace("patch","p")
    return version

def GetTaskDirName(analysis):
    '''
    Get the name of the CRAB task directory to be created. For the user's benefit this
    will include the CMSSW version and possibly important information from
    the dataset used, such as the bunch-crossing time.
    '''

    # Constuct basic task directory name
    if opts.dirName == "":
        dirName = "multicrab"
        dirName+= "_" + analysis
        dirName+= "_v" + GetCMSSW(opts)
        dirName+= "_" + opts.datetime.strftime("%Hh%Mm%Ss_%d%h%Y")
    else:
        dirName = opts.dirName + "_" + opts.datetime.strftime("%Hh%Mm%Ss_%d%h%Y")

    # If directory already exists
    if os.path.exists(dirName) and os.path.isdir(dirName):
	Verbose("The directory %s already exists" % (hs+dirName+ns), True)
    return dirName

def CreateTaskDirAndCfg(dirName, cfg):
    '''
    Create the CRAB task directory and copy inside relevant files (e.g. multicrab.cfg)
    '''
    # Create the pseudo-multicrab main directory ?
    if os.path.exists(dirName):
        msg = "Cannot create directory \"%s\". It already exists!" % (dirName)
        #raise Exception(msg)
        Verbose(msg)
    else:
        Verbose("%sCreated pseudo-multicrab directory %s%s" % (hs, dirName, ns))
        os.mkdir(dirName)

    # Create the "multicrab.cfg" file
    CreateCfgFile(cfg)

    # Write the commit id, "git status", "git diff" command output the directory created for the multicrab task
    gitFileList = git.writeCodeGitInfo(dirName, False)    
    Verbose("Copied %s to '%s'." % ("'" + "', '".join(gitFileList) + "'", dirName) )
    return

def CreateCfgFile(filePath="multicrab.cfg"):
    '''
    Alternative way
    cd <pseudo-multicrab>
    find * -maxdepth 0 -type d | awk '{print "["$1"]"}' > multicrab.cfg
    '''

    cmd = "touch %s" % filePath
    if os.path.isfile(filePath):
        Verbose("Cannot created file %s. It already exists" % (filePath), True)
    else:
        Verbose(cmd, True)
        os.system(cmd)
    return

def CreateDatasetDir(dirName, dsetName):

    datasetDir = os.path.join(dirName, dsetName)
    if os.path.exists(datasetDir) and os.path.isdir(datasetDir):
        Verbose("%sCannot create directory \"%s\". It already exists%s" % (es, datasetDir, ns), True)
    else:
        Verbose("%sCreated dataset directory %s%s" % (hs, datasetDir, ns))
        os.mkdir(datasetDir)
    return

def CreateJob(opts, args):
    '''
    Create a pseudo-multicrab directory.
    '''    

    # Get general info
    getDatasetInfo(opts)

    # Determine pseudo-multicrab directory name
    opts.taskDirName = GetTaskDirName(opts.algorithm)

    # Create CRAB task diractory
    cfg = os.path.join(opts.taskDirName, "multicrab.cfg")
    CreateTaskDirAndCfg(opts.taskDirName, cfg)

    # Create the dataset subdirectory
    CreateDatasetDir(opts.taskDirName, opts.datasetAlias)
        
    # Write the new dataset to multicrab.cfg
    multicrabCfg = open(cfg, "a")
    multicrabCfg.write("[" + str(opts.datasetAlias) + "]")
    multicrabCfg.write("\n")
    multicrabCfg.close()

    Verbose("Creating directory structure for dataset with name \"%s\"" % (opts.datasetName))
    dirs = ["results", "inputs"]
    # For-loop: All sub-directories to be created
    for d in dirs:
        newDir = os.path.join(opts.taskDirName, opts.datasetAlias, d)
        if os.path.exists(newDir) and os.path.isdir(newDir):
            #raise Exception("Cannot create directory \"%s\". It already exists!" % (newDir))
            Verbose("Cannot create directory \"%s\". It already exists!" % (newDir), True)
        else:
            os.mkdir(newDir)

        # The "results" directory
        resultsDir  = os.path.join(opts.taskDirName, opts.datasetAlias, "results")
        resultsFile = "histograms-%s.root" % (opts.datasetName)
        resultsPath = os.path.join(resultsDir, resultsFile)
	Verbose("Copying the ROOT file \"%s\" in the directory \"%s\"" % (opts.rootFile, resultsDir))
        cmd = "cp %s %s" % (opts.rootFile, resultsPath)
	os.system(cmd)

        # Save Auxiliary information and count4ers
        writeConfigInfoToRootFile(resultsPath, opts)
        moveAllHistosIntoAnalysisFolder(resultsPath, opts)
        writeCounters(resultsPath, opts) 
    return 0


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

    # Default Values
    VERBOSE      = False
    DIRNAME      = ""
    CMSSW        = "92X"
    ANALYSIS     = "HLTausAnalysis"
    ALGORITHM    = "CaloTk"
    DATAERA      = "TDR2019" #"ID2017"
    DATAVERSION  = None
    ROOTFILE     = None
    DATETIME     = datetime.datetime.now()
    
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))

    parser.add_option("-f", "--rootFile", dest="rootFile", default= ROOTFILE,
                      help="The ROOT file (raw2TTree.root) to be copied inside the multicrab dir[default: %s]" % (ROOTFILE) )

    parser.add_option("--dir", dest="dirName", default=DIRNAME, type="string",
                      help="Custom name for CRAB directory name. Overwrite automatic naming [default: %s]" % (DIRNAME))

    parser.add_option("--cmssw", dest="cmssw" , default=CMSSW, type="string",
                      help="CMSSW version for the CRAB directory name [default: %s]" % (CMSSW))

    parser.add_option("--dataVersion", dest="dataVersion", default=DATAVERSION,
                      help="Data version (<xy> where x=CMSSW release and y=mc or data) [default: %s]" % (DATAVERSION))

    parser.add_option("--analysisName", dest="analysisName", default=ANALYSIS, type="string",
                      help="The analysis name which corresponds to the folder where all histograms will be moved to in the new ROOT file [default: %s]" % (ANALYSIS))

    parser.add_option("--algorithm", dest="algorithm", default=ALGORITHM, type="string",
                      help="The tau algorithm analyzer used to create the output ROOT files (informative) [default: %s]" % (ALGORITHM))

    parser.add_option("--dataEra", dest="dataEra", default=DATAERA, type="string",
                      help="The dataEra name which will be appended to the analysis folder where all histograms will be moved to in the new ROOT file [default: %s]" % (DATAERA))

    parser.add_option("--datetime", dest="datetime", default=DATETIME, type="string",
                      help="Datetime to be appended to the pseudo-multicrab directory name[default: %s]" % (DATETIME))

    (opts, args) = parser.parse_args()

    # Require at least one argument (input ROOT file)
    if len(sys.argv) < 1:
        parser.print_help()
        sys.exit(1)
        
    # The ROOT file options must be provided
    opts.rootFiles = []
    if opts.rootFile == None:
        # For-loop: dir contents
        for dirName, subdirList, fileList in os.walk("."):
            for f in fileList:
                if f.endswith(".root"):
                    opts.rootFiles.append(f)
            # Do it just once
            break
        if len(opts.rootFiles) < 1:
            raise Exception("No specific ROOT file provided and no ROOT files found in current working directory")
        else:
            Print("No specific ROOT file provided. Will use all %d ROOT files in current directory (if any)" % (len(opts.rootFiles)), True)
    else:
        Verbose("The ROOT file(s) to be copied inside the multicrab dir is \"%s\"" % (opts.rootFile), True)

    # Data-era setting
    validEras = ["TP2015", "ID2017", "TDR2019"]
    if opts.dataEra not in validEras:
        Print("Invalid data-eta \"%s\". Please select one of the following:\n\t%s" % (opts.dataEra, ", ".join(validEras)), True)
        sys.exit()


    #  Set ROOT ingore level for verbosity
    ROOT.gErrorIgnoreLevel = ROOT.kFatal #kPrint = 0,  kInfo = 1000, kWarning = 2000, kError = 3000, kBreak = 4000, kSysError = 5000, kFatal = 6000   
    
    # Create the pseudo-multicrab dir (if the ROOT file exists)
    nFiles = len(opts.rootFiles)    
    if opts.rootFile != None:
        Verbose("Single ROOT file mode", True)

        # Check & assign ROOT file 
        CheckRootFile(f)
        opts.rootFile = f
        
        # Create directory and subdirectories
        CreateJob(opts, args)
        aux.Print("%sFile %d/%d: %s%s" % (hs, i, nFiles, f, ns), i==1)
        
        Print("Created pseudo-multicrab directory %s!" % (ss + opts.dirName + ns), True)        
        sys.exit()
    elif nFiles > 0: 
        Verbose("Multiple ROOT files mode", True)

        # For-loop: All ROOT files
        for i, f in enumerate(opts.rootFiles, 1):
            
            # Check & assign ROOT file 
            CheckRootFile(f)
            opts.rootFile = f
            
            # Create directory and subdirectories
            CreateJob(opts, args)
            aux.PrintFlushed("%sFile %d/%d: %s%s" % (hs, i, nFiles, f, ns), i==1)
        print

        Print("Created pseudo-multicrab %s" % (ss + opts.taskDirName + ns), True)
        sys.exit()
    else:
        sys.exit()

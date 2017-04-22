#!/usr/bin/env python
'''
Description:
This script is used to create a pseudo multicrab directory, using an input ROOT file (miniaod2tree.root).
The purpose is primarily to enable the easy testing of local changes to the MiniAOD2TTree code and the 
NtupleAnalysis code using the hplusGenerateDataFormats.py script.
A script execution will thus create an empty multicrab with identical name and structure as those created
by the multicrab.py script. It will contain a single dataset with a single ROOT file under results/ dir, which is 
a mere copy of the file used as input for the script execution (only renamed) to histograms-<dataset>.root.

Usage:
hltausPseudoMulticrab.py -f test.root [opts]

To add  a dataset to an existing <multicrab_dir>:
hltausPseudoMulticrab.py -f test.root -r <multicrab_dir> 

Examples:
hltausPseudoMulticrab.py -f test.root 
hltausPseudoMulticrab.py -f test.root --cmssw "612SLHC6mc"
hltausPseudoMulticrab.py -f test.root --dir multicrab_test -v 
hltausPseudoMulticrab.py -f test.root --dir multicrab_CaloPlusTracks_v_20170419T1926/ --dataset Neutrino_Pt2to20_gun
hltausPseudoMulticrab.py -f results/L1CaloTaus_CaloCorr_TTTracks_Stubs_TTPixelTracks_CandPixHits_TPs_GenPs_v620SLHC12p1_07Nov2016/CaloPlusTracks_Histograms_VBF.root
hltausPseudoMulticrab.py -f results/L1CaloTaus_CaloCorr_TTTracks_Stubs_TTPixelTracks_CandPixHits_TPs_GenPs_v620SLHC12p1_07Nov2016/CaloPlusTracks_Histograms_MinBias.root --dataset Neutrino_Pt2to20_gun --dir multicrab_CaloPlusTracks_v612SLHC6_20170420T1226

Last Used:
hltausPseudoMulticrab.py -f CaloPlusTracks_Histograms_VBF.root
hltausPseudoMulticrab.py -f CaloPlusTracks_Histograms_MinBias.root --dataset Neutrino_Pt2to20_gun --dir multicrab_CaloPlusTracks_v61XSLHC6_20170421T1554
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
from HLTausAnalysis.NtupleAnalysis.tools.datasets import *


#================================================================================================ 
# Function Definitions
#================================================================================================ 
def FileExists(filePath, opts):
    '''
    Checks that a file exists by executing the ls command for its full path, 
    or the corresponding "EOS" command if opts.filesInEOS is enabled.
    '''
    Verbose("FileExists()", False)
    
    if "CRAB3_TransferData" in filePath: # fixme: just a better alternative to "if opts.filesInEOS:"
        cmd = ConvertCommandToEOS("ls", opts) + " " + filePath
        ret = Execute("%s" % (cmd) )
        # If file is not found there won't be a list of files; there will be an error message
        errMsg = ret[0]
        if "Unable to stat" in errMsg:
            return False
        elif errMsg == filePath.split("/")[-1]:
            return True
        else:
            raise Exception("This should not be reached! Execution of command %s returned %s" % (cmd, errMsg))
    else:
        if os.path.isfile(filePath):
            return True
        else:
            return False
    return True


def writeConfigInfoToRootFile(fileName, opts):
    '''
    Creates configInfo folder where the following objects are written:
    - dataVersion (CMSSW version, followed by "mc" or "data". e.g. 80Xmc
    - codeVersion (git commit id)
    - configinfo (histogram with COM energy, luminosity, pileup, etc..)
    '''
    Verbose("writeConfigInfoToRootFile()", True)
    
    if FileExists(fileName, opts ) == False:
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
    Verbose("writeCounters()", True)
    
    if FileExists(fileName, opts ) == False:
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
    Verbose("moveAllHistosIntoAnalysisFolder()", True)
    
    if FileExists(fileName, opts ) == False:
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
    Print("Closing file %s." % (fOUT.GetName()) )
    fOUT.Close()
    return


def GetRootFile(fileName, fileMode="UPDATE"):
    '''
    '''
    Verbose("GetRootFile()", True)

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
    Verbose("GetSelfName()")    
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
    Verbose("AskUser()")
    
    keystroke = raw_input("\t" +  msg + " (y/n): ")
    if (keystroke.lower()) == "y":
        return True
    elif (keystroke.lower()) == "n":
        return False
    else:
        AskUser(msg)
    

def GetCMSSW(opts):
    '''
    Get a command-line-friendly format of the CMSSW version currently use.
    https://docs.python.org/2/howto/regex.html
    '''
    Verbose("GetCMSSW()")

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


def GetAnalysis():
    '''
    Get the analysis type. This will later-on help determine the datasets to be used.
    https://docs.python.org/2/howto/regex.html
    '''
    Verbose("GetAnalysis()")
    
    # Create a compiled regular expression object
    leg_re = re.compile("miniAOD2TTree_(?P<leg>\S+)_cfg.py")

    # Scan through the string 'pwd' & look for any location where the compiled RE 'cmssw_re' matches
    match = leg_re.search(opts.pset)

    # Return the string matched by the RE. Convert to desirable format
    analysis = "DUMMY"
    if match:
	analysis = match.group("leg")
    else:
        raise Exception("Could not determine the analysis type from the PSET \"%s\"" % (opts.pset) )

    return analysis


def AskToContinue(taskDirName, analysis, opts):
    '''
    Inform user of the analysis type and datasets to be user in the multicrab job creation. Offer chance to abort sequence 
    '''
    Verbose("AskToContinue()")

    Print("Creating pseudo-multicrab directory \"%s\" using the file \"%s\" as input" % (taskDirName, opts.rootFile) )
    #DatasetGroup(analysis).PrintDatasets(False)
    
    AbortTask(keystroke="q")
    return


def AbortTask(keystroke):
    '''
    Give user last chance to abort CRAB task creation.
    '''
    Verbose("AbortTask()")
    if not opts.ask:
        return

    message  = "=== %s:\n\tPress \"%s\" to abort, any other key to proceed: " % (GetSelfName(), keystroke)
    response = raw_input(message)
    if (response!= keystroke):
        return
    else:
        print "=== %s:\n\tEXIT" % (GetSelfName())
        sys.exit()
    return


def GetTaskDirName(analysis, version, datasets):
    '''
    Get the name of the CRAB task directory to be created. For the user's benefit this
    will include the CMSSW version and possibly important information from
    the dataset used, such as the bunch-crossing time.
    '''
    Verbose("GetTaskDirName()")

    # Constuct basic task directory name
    dirName = "multicrab"
    dirName+= "_"  + analysis
    dirName+= "_v" + version
    
    # Add dataset-specific info, like bunch-crossing info
    bx_re = re.compile("\S+(?P<bx>\d\dns)_\S+")
    match = bx_re.search(datasets[0].URL)
    if match:
	dirName+= "_"+match.group("bx")

    # Append the creation time to the task directory name    
    # time = datetime.datetime.now().strftime("%d%b%Y_%Hh%Mm%Ss")
    time = datetime.datetime.now().strftime("%Y%m%dT%H%M")
    dirName+= "_" + time

    # If directory already exists (resubmission)
    if os.path.exists(opts.dirName) and os.path.isdir(opts.dirName):
	dirName = opts.dirName
    return dirName


def CreateTaskDir(dirName):
    '''
    Create the CRAB task directory and copy inside it the PSET to be used for the CRAB job.
    '''
    Verbose("CreateTaskDir()")
    
    if  os.path.exists(dirName):
        msg = "Cannot create directory \"%s\". It already exists!" % (dirName)
        #raise Exception(msg)
        Print(msg)
        return

    # Create the pseudo multicrab directory
    Verbose("mkidr %s" % (dirName))
    os.mkdir(dirName)

    # Place an empty PSet file inside the pseudo multicrab directory
    cmd = "touch %s" % dirName + "/" + opts.pset
    Verbose(cmd)
    os.system(cmd)

    # Create the multicrab.cfg
    multicrabCfg = dirName + "/" + "multicrab.cfg"
    cmd = "touch %s" % multicrabCfg
    if not os.path.exists(multicrabCfg):
        Verbose(cmd)
        os.system(cmd)

    # Write the commit id, "git status", "git diff" command output the directory created for the multicrab task
    gitFileList = git.writeCodeGitInfo(dirName, False)    
    Verbose("Copied %s to '%s'." % ("'" + "', '".join(gitFileList) + "'", dirName) )
    return


def GetRequestName(dataset):
    '''
    Return the file name and path to an (empty) crabConfig_*.py file where "*" 
    contains the dataset name and other information such as tune, COM, Run number etc..
    of the Data or MC sample used
    '''
    Verbose("GetRequestName()")
    
    # Create compiled regular expression objects
    datadataset_re = re.compile("^/(?P<name>\S+?)/(?P<run>Run\S+?)/")
    mcdataset_re   = re.compile("^/(?P<name>\S+?)/")
    tune_re        = re.compile("(?P<name>\S+)_Tune")
    tev_re         = re.compile("(?P<name>\S+)_13TeV")
    ext_re         = re.compile("(?P<name>_ext\d+)-")
    runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_")
    # runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_13TeV_PromptReco_Collisions15(?P<BunchSpacing>\S*)_JSON(?P<Silver>(_\S+|))\.")
    # runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_13TeV_PromptReco_Collisions15(?P<BunchSpacing>\S*)_JSON")
    # runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_13TeV_PromptReco_Collisions15_(?P<BunchSpacing>\d+ns)_JSON_v")
    
    # Scan through the string 'dataset.URL' & look for any location where the compiled RE 'mcdataset_re' matches
    match = mcdataset_re.search(dataset.URL)
    if dataset.isData():
	match = datadataset_re.search(dataset.URL)
        
    # Append the dataset name
    if match:
	requestName = match.group("name")

    # Append the Run number (for Data samples only)
    if dataset.isData():
	requestName+= "_"
	requestName+= match.group("run")

    # Append the MC-tune (for MC samples only) 
    tune_match = tune_re.search(requestName)
    if tune_match:
	requestName = tune_match.group("name")

    # Append the COM Energy (for MC samples only) 
    tev_match = tev_re.search(requestName)
    if tev_match:
	requestName = tev_match.group("name")

    # Append the Ext
    ext_match = ext_re.search(dataset.URL)
    if ext_match:
	requestName+=ext_match.group("name")

    # Append the Run Range (for Data samples only)
    if dataset.isData():
	runRangeMatch = runRange_re.search(dataset.lumiMask)
	if runRangeMatch:
	    runRange= runRangeMatch.group("RunRange")
	    runRange = runRange.replace("-","_")
	    #bunchSpace = runRangeMatch.group("BunchSpacing")
	    requestName += "_" + runRange #+ bunchSpace
	    #Ag = runRangeMatch.group("Silver")
	    #if Ag == "_Silver": # Use  chemical element of silver (Ag)
            #    requestName += Ag

    # Finally, replace dashes with underscores    
    requestName = requestName.replace("-","_")

    return requestName


def EnsurePathDoesNotExist(taskDirName, requestName):
    '''
    Ensures that file does not already exist
    '''
    Verbose("EnsurePathDoesNotExist()")
    
    filePath = os.path.join(taskDirName, requestName)
    
    if not os.path.exists(filePath):
	return
    else:
        msg = "File '%s' already exists!" % (filePath)
        if AskUser(msg + " Proceed and overwrite it?"):
            return
	else:
            raise Exception(msg)
    return


def CreateJob(opts, args):
    '''
    Create & submit a CRAB task, using the user-defined PSET and list of datasets.
    '''
    Verbose("CreateJob()")
    
    # Get general info
    version     = GetCMSSW(opts)
    analysis    = GetAnalysis()
    dataset = None
    for d in DatasetGroup("All").GetDatasetList():
        if opts.dataset in d.URL.replace("-", "_"):
            dataset = d
    if dataset == None:
        raise Exception("Could not find dataset object for dataset with name \"%s\"." % (opts.dataset) )
    else:
        datasets= [dataset]

    if opts.dirName == "":
        taskDirName = GetTaskDirName(analysis, version, datasets)
    else:
        taskDirName = opts.dirName

    # Give user last chance to abort
    AskToContinue(taskDirName, analysis, opts)
    
    # Create CRAB task diractory
    if opts.dirName == "":
        CreateTaskDir(taskDirName)

    # Create the "multicrab.cfg" file
    multicrabCfg = open(taskDirName + "/" + "multicrab.cfg", 'a')
        
    # For-loop: All datasets [always 1 in this case]
    for dataset in datasets:        
	Verbose("Determining request name for dataset with URL \"%s\"" % (dataset.URL))
        requestName = GetRequestName(dataset)

	Verbose("Creating directory for dataset with request name \"%s\"" % (requestName))
        datasetDir = os.path.join(taskDirName, requestName)
        if os.path.exists(datasetDir) and os.path.isdir(datasetDir):
            raise Exception("Cannot create directory \"%s\". It already exists!" % (datasetDir))
        else:
            os.mkdir(datasetDir)

        # Write dataset to multicrab.cfg
        multicrabCfg.write("[" + str(requestName) + "]")
        multicrabCfg.write("\n")
        
	Verbose("Creating directory structure for dataset with request name \"%s\"" % (requestName))
        dirs = ["results", "inputs"]
        for d in dirs:
            newDir = os.path.join(datasetDir, d)
            if os.path.exists(newDir) and os.path.isdir(newDir):
                raise Exception("Cannot create directory \"%s\". It already exists!" % (newDir))
            else:
                os.mkdir(newDir)

        resultsDir  = os.path.join(datasetDir, "results")
        resultsFile = "histograms-%s.root" % (requestName)
        resultsPath = os.path.join(resultsDir, resultsFile)
	Verbose("Copying the ROOT file \"%s\" in the directory \"%s\"" % (opts.rootFile, resultsDir))
        cmd = "cp %s %s" % (opts.rootFile, resultsPath)
	os.system(cmd)

        writeConfigInfoToRootFile(resultsPath, opts)
        moveAllHistosIntoAnalysisFolder(resultsPath, opts)
        writeCounters(resultsPath, opts)

    Print("Successfully created pseudo-multicrab directory \"%s\" " % (taskDirName))
    multicrabCfg.close()
    if 0:
        os.system("ls -lt")
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
    ASK          = False
    PSET         = "miniAOD2TTree_CaloPlusTracks_cfg.py"
    DIRNAME      = ""
    DATASET      = ""
    CMSSW        = "61XSLHC6"
    ANALYSIS     = "HLTausAnalysis"
    DATAERA      = "TP2015" #"TDR2019"
    DATAVERSION  = CMSSW + "mc"
    ROOTFILE     = "test.root"
    NOMINALLUMI  = 5 # 5e34 (<PU>=140)
    ULTIMATELUMI = 7 # 7e34 (<PU>=200)

    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))

    parser.add_option("--dataset", dest="dataset", default="VBF_HToTauTau_125_14TeV_powheg_pythia6", 
                      help="Dataset to include in multicrab dir [default: %s]" % (DATASET))

    parser.add_option("-f", "--rootFile", dest="rootFile", default= ROOTFILE,
                      help="The ROOT file (miniaod2tree.root) to be copied inside the multicrab dir[default: %s]" % (ROOTFILE) )

    parser.add_option("-p", "--pset", dest="pset", default=PSET, type="string",
                      help="The python cfg file to be used by cmsRun [default: %s]" % (PSET))

    parser.add_option("--dir", dest="dirName", default=DIRNAME, type="string",
                      help="Custom name for CRAB directory name. Overwrite automatic naming [default: %s]" % (DIRNAME))

    parser.add_option("--cmssw", dest="cmssw" , default=CMSSW, type="string",
                      help="CMSSW version for the CRAB directory name [default: %s]" % (CMSSW))

    parser.add_option("--dataVersion", dest="dataVersion", default=DATAVERSION,
                      help="Data version (<xy> where x=CMSSW release and y=mc or data) [default: %s]" % (DATAVERSION))

    parser.add_option("--ask", dest="ask", default=ASK,
                      help="Ask before proceeding to create the pseudo-multicrab directory [default: %s]" % (ASK))

    parser.add_option("--lumi", dest="lumi", default=NOMINALLUMI, type=int,
                      help="The lumonisoty coefficient \"a\" that corresponds to the luminosity of the samples  \"a\"E+34 [default: %s]" % (NOMINALLUMI))

    parser.add_option("--analysisName", dest="analysisName", default=ANALYSIS, type="string",
                      help="The analysis name which corresponds to the folder where all histograms will be moved to in the new ROOT file [default: %s]" % (ANALYSIS))

    parser.add_option("--dataEra", dest="dataEra", default=DATAERA, type="string",
                      help="The dataEra name which will be appended to the analysis folder where all histograms will be moved to in the new ROOT file [default: %s]" % (DATAERA))


    (opts, args) = parser.parse_args()

    # Require at least one argument (input ROOT file)
    if len(sys.argv) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        pass
    # The ROOT file options must be provided
    if opts.rootFile == None:
        raise Exception("Must provide a ROOT file (miniaod2tree.root) as argument!")

    # Luminosity setting
    validLumis = [NOMINALLUMI, ULTIMATELUMI]
    if opts.lumi not in validLumis:
        Print("Invalid value for lumi (\"%s\"). Luminosity value must correspond to either the HL-LHC nominal luminosity %s (E+34) or the HL-LHC ultimate luminosity %s (E+34)." % (opts.lumi, NOMINALLUMI, ULTIMATELUMI), True)
        sys.exit()

    # Pileup setting        
    if opts.lumi == NOMINALLUMI:
        opts.pileup = 140
    else:
        opts.pileup = 200
    Print("The average pileup for the HL-LHC luminosity %sE+34 is set to <PU>=%s" % (opts.lumi, opts.pileup), True)

    # Dataset check
    validDataset  = True
    validDatasets = []
    for d in DatasetGroup("All").GetDatasetList():
        dName = GetRequestName(d)
        validDatasets.append(dName)
        if opts.dataset == dName:
            validDataset = True
        else:
            pass

    if not validDataset:
        Print("Invalid dataset %s. Please select one of the following:\n\t" % (opts.dataset, ", ".join(validDatasets)), True)
        sys.exit()
    
    # Additional sanity checks
    datasets = ["VBF_HToTauTau_125_14TeV_powheg_pythia6", "Neutrino_Pt2to20_gun", "PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV", "TauThreeProngs"]
    if opts.dataset not in datasets:
        Print("Dataset %s not valid. Please select one of the following:\n\t%s" % (opts.dataset, ", ".join(datasets)))
        sys.exit()

    # Create the pseudo-multicrab dir (if the ROOT file exists)
    if os.path.exists(opts.rootFile):
        sys.exit( CreateJob(opts, args) )
    else:
        raise Exception("The ROOT file provided (%s) does not exists!" % (opts.rootFile) )

#lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt" # ICHEP dataset 271036-276811
lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"


#================================================================================================ 
# Class Definition
#================================================================================================ 
import os

class Dataset:
    def __init__(self, url, dbs="global", dataVersion="80Xmc", dasQuery="", lumiMask=lumiMask, name=""):
        self.URL = url
        self.DBS = dbs
        self.dataVersion = dataVersion
        if not os.path.dirname(lumiMask):
            #lumiMask = os.path.join(os.environ['CMSSW_BASE'],"src/HiggsAnalysis/MiniAOD2TTree/data",lumiMask)
            lumiMask = os.path.join("", "src/HiggsAnalysis/MiniAOD2TTree/data",lumiMask)
        self.lumiMask = lumiMask
        self.DASquery = dasQuery
	self.name = name

    def isData(self):
        if "data" in self.dataVersion:
            return True
        return False

    def getName(self):
	return self.name

    
#================================================================================================ 
# MC-UpgFall13d, BE5D samples (PU = 140) - https://hypernews.cern.ch/HyperNews/CMS/get/slhc-triggersim/285.html
#================================================================================================ 
datasetsVBF_UpgFall13d = []
das = "https://cmsweb.cern.ch/das/request?input=dataset%3D%2FVBF_HToTauTau_125_14TeV_powheg_pythia6%2FUpgFall13d-PU140bx25_POSTLS261_V3-v1%2FGEN-SIM-DIGI-RAW&instance=prod%2Fglobal"
datasetsVBF_UpgFall13d.append(Dataset('/VBF_HToTauTau_125_14TeV_powheg_pythia6/UpgFall13d-PU140bx25_POSTLS261_V3-v1/GEN-SIM-DIGI-RAW', dataVersion="612SLHC6mc", dasQuery=das))

datasetsMinBias_UpgFall13d = []
das = "https://cmsweb.cern.ch/das/request?input=dataset%3D%2FNeutrino_Pt2to20_gun%2FUpgFall13d-PU140bx25_POSTLS261_V3-v1%2FGEN-SIM-DIGI-RAW&instance=prod%2Fglobal"
datasetsMinBias_UpgFall13d.append(Dataset('/Neutrino_Pt2to20_gun/UpgFall13d-PU140bx25_POSTLS261_V3-v1/GEN-SIM-DIGI-RAW', dataVersion="612SLHC6mc", dasQuery=das))

datasetsTTbar_UpgFall13d = []
das = "https://cmsweb.cern.ch/das/request?input=dataset%3D%2FPYTHIA6_Tauola_TTbar_TuneZ2star_14TeV%2FUpgFall13d-PU140bx25_POSTLS261_V3-v1%2FGEN-SIM-DIGI-RAW&instance=prod%2Fglobal"
datasetsTTbar_UpgFall13d.append(Dataset('/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/UpgFall13d-PU140bx25_POSTLS261_V3-v1/GEN-SIM-DIGI-RAW', dataVersion="612SLHC6mc", dasQuery=das))

datasetsTauThreeProngs_UpgFall13d = []
das = "https://cmsweb.cern.ch/das/request?input=dataset%3D%2FTauThreeProngs%2FUpgFall13d-PU140bx25_POSTLS261_V3-v1%2FGEN-SIM-DIGI-RAW&instance=prod%2Fglobal"
datasetsTauThreeProngs_UpgFall13d.append(Dataset('/TauThreeProngs/UpgFall13d-PU140bx25_POSTLS261_V3-v1/GEN-SIM-DIGI-RAW', dataVersion="612SLHC6mc", dasQuery=das))


#================================================================================================ 
# Dataset Grouping
#================================================================================================ 
TP2015Datasets = []
TP2015Datasets.extend(datasetsVBF_UpgFall13d)
TP2015Datasets.extend(datasetsMinBias_UpgFall13d)
TP2015Datasets.extend(datasetsTTbar_UpgFall13d)
TP2015Datasets.extend(datasetsTauThreeProngs_UpgFall13d)

TDR2019Datasets = []


#================================================================================================ 
# Class Definition
#================================================================================================ 
class DatasetGroup:
    def __init__(self, analysis):
        self.verbose   = False
        self.analysis  = analysis
        self.GroupDict = {}
        self.CreateGroups()

    def SetVerbose(verbose):
        self.verbose = verbose
        return


    def Verbose(self, msg, printHeader=False):
        '''
        Simple print function. If verbose option is enabled prints, otherwise does nothing.
        '''
        if not self.verbose:
            return
        self.Print(msg, printHeader)
        return


    def Print(self, msg, printHeader=True):
        '''
        Simple print function. If verbose option is enabled prints, otherwise does nothing.
        '''
        fName = __file__.split("/")[-1]
        cName = self.__class__.__name__
        name  = fName + ": " + cName
        if printHeader:
                print "=== ", name
        print "\t", msg
        return


    def CreateGroups(self):
        '''
        Create dataset grouping in a dictionary for easy access.
        '''

        analyses = ["SignalAnalysis", "Hplus2tbAnalysis", "TauLeg", "METLeg", "L1Study", "All"]
        if self.analysis not in analyses:
            raise Exception("Unknown analysis \"%s\". Please select one of the following: \"%s" % (self.analysis, "\", \"".join(analyses) + "\".") )


        self.GroupDict["TP2015"]  = TP2015Datasets
        self.GroupDict["TDR2019"] = TDR2019Datasets
        self.GroupDict["All"]     = TP2015Datasets + TDR2019Datasets
        return


    def GetDatasetList(self):
        '''
        Return the dataset list according to the analysis name. 
        Uses pre-defined dictionary mapping: analysis->dataset list
        '''
        return self.GroupDict[self.analysis]


    def PrintDatasets(self, printHeader=False):
        '''
        Print all datasets for given analysis
        '''
        datasetList = self.GroupDict[self.analysis]

        if printHeader==True:
            self.Print("The datasets for analysis \"%s\" are:\n\t%s" % (self.analysis, "\n\t".join(str(d.URL) for d in datasetList) ), True)
        else:
            self.Print("\n\t".join(str(d.URL) for d in datasetList), False)
        return

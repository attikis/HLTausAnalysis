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
# PhaseIISpring17D
#================================================================================================ 
dasAll  = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2F*PhaseIISpring17D*%2F*+status%3D*"
dasDone = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2F*PhaseIISpring17D*%2F*"
das     = dasAll

datasetsBs_PhaseIISpring17D = []
datasetsBs_PhaseIISpring17D.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsBs_PhaseIISpring17D.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsBs_PhaseIISpring17D.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))

datasetsSingleE_PhaseIISpring17D = []
datasetsSingleE_PhaseIISpring17D.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleE_PhaseIISpring17D.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleE_PhaseIISpring17D.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSingleMu_PhaseIISpring17D = []
datasetsSingleMu_PhaseIISpring17D.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleMu_PhaseIISpring17D.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleMu_PhaseIISpring17D.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSingleNu_PhaseIISpring17D = []
datasetsSingleNu_PhaseIISpring17D.append(Dataset('/SingleNeutrino/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsSingleNu_PhaseIISpring17D.append(Dataset('/SingleNeutrino/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSinglePhoton_PhaseIISpring17D = []
datasetsSinglePhoton_PhaseIISpring17D.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePhoton_PhaseIISpring17D.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePhoton_PhaseIISpring17D.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSinglePion0_PhaseIISpring17D = []
datasetsSinglePion0_PhaseIISpring17D.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion0_PhaseIISpring17D.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion0_PhaseIISpring17D.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSinglePion_PhaseIISpring17D = []
datasetsSinglePion_PhaseIISpring17D.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion_PhaseIISpring17D.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion_PhaseIISpring17D.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSingleTau_PhaseIISpring17D = []
datasetsSingleTau_PhaseIISpring17D.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsSingleTau_PhaseIISpring17D.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsSingleTau_PhaseIISpring17D.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))

datasetsTTJpsiFilter_PhaseIISpring17D = []
datasetsTTJpsiFilter_PhaseIISpring17D.append(Dataset('/TT_JpsiFilter_TuneCUETP8M1_mtop166_5_14TeV-powheg-tauola-pythia8/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTTJpsiFilter_PhaseIISpring17D.append(Dataset('/TT_JpsiFilter_TuneCUETP8M1_mtop166_5_14TeV-powheg-tauola-pythia8/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsTT_PhaseIISpring17D = []
datasetsTT_PhaseIISpring17D.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTT_PhaseIISpring17D.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsTT_PhaseIISpring17D.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-PU200_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))

datasetsTTbar_PhaseIISpring17D = []
datasetsTTbar_PhaseIISpring17D.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTTbar_PhaseIISpring17D.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTTbar_PhaseIISpring17D.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-PU200_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 


#================================================================================================ 
# Dataset Grouping
#================================================================================================ 
TP2015Datasets = []
TP2015Datasets.extend(datasetsVBF_UpgFall13d)
TP2015Datasets.extend(datasetsMinBias_UpgFall13d)
TP2015Datasets.extend(datasetsTTbar_UpgFall13d)
TP2015Datasets.extend(datasetsTauThreeProngs_UpgFall13d)

ID2017Datasets = []
ID2017Datasets.extend(datasetsBs_PhaseIISpring17D)
ID2017Datasets.extend(datasetsSingleE_PhaseIISpring17D)
ID2017Datasets.extend(datasetsSingleMu_PhaseIISpring17D)
ID2017Datasets.extend(datasetsSingleNu_PhaseIISpring17D)
ID2017Datasets.extend(datasetsSinglePhoton_PhaseIISpring17D)
ID2017Datasets.extend(datasetsSinglePion0_PhaseIISpring17D)
ID2017Datasets.extend(datasetsSinglePion_PhaseIISpring17D)
ID2017Datasets.extend(datasetsSingleTau_PhaseIISpring17D)
ID2017Datasets.extend(datasetsTTJpsiFilter_PhaseIISpring17D)
ID2017Datasets.extend(datasetsTT_PhaseIISpring17D)
ID2017Datasets.extend(datasetsTTbar_PhaseIISpring17D)

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

        analyses = ["TP2015", "ID2017", "TDR2019", "All"]

        if self.analysis not in analyses:
            raise Exception("Unknown analysis \"%s\". Please select one of the following: \"%s" % (self.analysis, "\", \"".join(analyses) + "\".") )


        self.GroupDict["TP2015"]  = TP2015Datasets
        self.GroupDict["ID2017"]  = ID2017Datasets
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

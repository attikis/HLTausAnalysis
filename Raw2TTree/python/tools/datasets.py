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
            lumiMask = os.path.join(os.environ['CMSSW_BASE'],"src/HiggsAnalysis/MiniAOD2TTree/data",lumiMask)
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
# Data
#================================================================================================ 
datasetsTauData = []
das = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FTau%2FRun2016*-PromptReco-v2%2FMINIAOD"
datasetsTauData.append(Dataset('/Tau/Run2016B-03Feb2017_ver2-v2/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_273150-275376_13TeV_PromptReco_Collisions16_JSON_Tau_Run2016B.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016C-03Feb2017-v1/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_275420-276283_13TeV_PromptReco_Collisions16_JSON_Tau_Run2016C.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016D-03Feb2017-v1/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_276315-276437_13TeV_PromptReco_Collisions16_JSON_Tau_Run2016D_MET80.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016D-03Feb2017-v1/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_276453-276811_13TeV_PromptReco_Collisions16_JSON_Tau_Run2016D_noMET80.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016E-03Feb2017-v1/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_276824-277420_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016E.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016F-03Feb2017-v1/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_277816-278800_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016F_HIP.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016F-03Feb2017-v1/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_278801-278808_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016F_HIPfixed.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016G-03Feb2017-v1/MINIAOD', dataVersion="80Xdata", dasQuery=das, lumiMask="Cert_278816-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016G.txt"))
#datasetsTauData.append(Dataset('/Tau/Run2016H-PromptReco-v1/MINIAOD', dataVersion="80Xdata2016H", dasQuery=das, lumiMask="Cert_281010-281202_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016H.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016H-03Feb2017_ver2-v1/MINIAOD', dataVersion="80Xdata2016H", dasQuery=das, lumiMask="Cert_281207-284035_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016H.txt"))
datasetsTauData.append(Dataset('/Tau/Run2016H-03Feb2017_ver3-v1/MINIAOD', dataVersion="80Xdata2016H", dasQuery=das, lumiMask="Cert_284036-284068_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016H.txt"))

#================================================================================================ 
# PhaseIISpring17D
#================================================================================================ 
dasAll  = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2F*PhaseIISpring17D*%2F*+status%3D*"
dasDone = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2F*PhaseIISpring17D*%2F*"
das     = dasAll

datasetsBs = []
datasetsBs.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsBs.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsBs.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))

datasetsSingleE = []
datasetsSingleE.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleE.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleE.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSingleMu = []
datasetsSingleMu.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleMu.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSingleMu.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSingleNu = []
datasetsSingleNu.append(Dataset('/SingleNeutrino/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsSingleNu.append(Dataset('/SingleNeutrino/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSinglePhoton = []
datasetsSinglePhoton.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePhoton.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePhoton.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSinglePion0 = []
datasetsSinglePion0.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion0.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion0.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSinglePion = []
datasetsSinglePion.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsSinglePion.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsSingleTau = []
datasetsSingleTau.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsSingleTau.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))
datasetsSingleTau.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das))

datasetsTTJpsiFilter = []
datasetsTTJpsiFilter.append(Dataset('/TT_JpsiFilter_TuneCUETP8M1_mtop166_5_14TeV-powheg-tauola-pythia8/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTTJpsiFilter.append(Dataset('/TT_JpsiFilter_TuneCUETP8M1_mtop166_5_14TeV-powheg-tauola-pythia8/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsTT = []
datasetsTT.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTT.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTT.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-PU200_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 

datasetsTTbar = []
datasetsTTbar.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTTbar.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 
datasetsTTbar.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-PU200_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc", dasQuery=das)) 


#================================================================================================ 
# Dataset Grouping
#================================================================================================ 
L1TauDatasets = []
#L1TauDatasets.extend(datasetsSinglePion0)
#L1TauDatasets.extend(datasetsSinglePion)
#L1TauDatasets.extend(datasetsSingleTau)
#L1TauDatasets.extend(datasetsSingleNu)
#L1TauDatasets.extend(datasetsTTJpsiFilter)
L1TauDatasets.extend(datasetsTT)
#L1TauDatasets.extend(datasetsTTbar)

AllDatasets = []
AllDatasets.extend(datasetsBs)
AllDatasets.extend(datasetsSingleE)
AllDatasets.extend(datasetsSingleMu)
AllDatasets.extend(datasetsSingleNu)
AllDatasets.extend(datasetsSinglePhoton)
AllDatasets.extend(datasetsSinglePion0)
AllDatasets.extend(datasetsSinglePion)
AllDatasets.extend(datasetsSingleTau)
AllDatasets.extend(datasetsTTJpsiFilter)
AllDatasets.extend(datasetsTT)
AllDatasets.extend(datasetsTTbar)


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
        analyses = ["CaloTk", "TkCalo", "PFTau", "CaloPix", "All"]
        if self.analysis not in analyses:
            raise Exception("Unknown analysis \"%s\". Please select one of the following: \"%s" % (self.analysis, "\", \"".join(analyses) + "\".") )
        self.GroupDict["CaloTk"]  = L1TauDatasets
        self.GroupDict["TkCalo"]  = L1TauDatasets
        self.GroupDict["PFTau"]   = L1TauDatasets
        self.GroupDict["CaloPix"] = L1TauDatasets
        self.GroupDict["All"]     = AllDatasets
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

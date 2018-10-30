#lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt" # ICHEP dataset 271036-276811
lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"


#================================================================================================ 
# Class Definition
#================================================================================================ 
import os

class Dataset:
    def __init__(self, url, dbs="global", dataVersion="80Xmc", dasQuery="", alias="", lumiMask=lumiMask):
        self.URL = url
        self.DBS = dbs
        self.dataVersion = dataVersion
        if not os.path.dirname(lumiMask):
            lumiMask = os.path.join(os.environ['CMSSW_BASE'],"src/HiggsAnalysis/MiniAOD2TTree/data",lumiMask)
        self.lumiMask = lumiMask
        self.DASquery = dasQuery
	self.alias = alias
        return

    def isData(self):
        if "data" in self.dataVersion:
            return True
        return False

    def isMC(self):
        return not self.isData()

    def getAlias(self):
	return self.alias

    def getFriendlyUrl(self):
        if self.URL.startswith("/"):
            friendlyURL = self.URL[1:]
        else:
            friendlyURL = self.URL[1:]

        friendlyURL = friendlyURL.replace("-", "_").replace("/", "_")
	return friendlyURL

    def getName(self):
	return self.getUrl()


    def getUrl(self):
	return self.URL

    def getDBS(self):
	return self.DBS

    def getDataVersion(self):
	return self.dataVersion

    def getDataVersion(self):
	return self.dataVersion

    def getDasQuery(self):
	return self.DASquery

    def getLumiMask(self):
	return self.lumiMask

    def getPU(self):
        return self.getPileup()

    def getPileup(self):
        pu = self.alias.split("L1T")[1]
        if "noPU" in pu:
            return 0
        elif "140" in pu:
            return 140
        elif "200" in pu:
            return 200
        else:
            raise Exception("Could not determine Pileup value for dataset=%s (alias = %s)" % (self.URL, self.alias) )


#================================================================================================ 
# Data
#================================================================================================ 
datasetsData = []
das = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2FTau%2FRun2016*-PromptReco-v2%2FMINIAOD"
#datasetsTauData.append(Dataset('/Tau/Run2016H-03Feb2017_ver3-v1/MINIAOD', dataVersion="80Xdata2016H", dasQuery=das, lumiMask="Cert_284036-284068_13TeV_PromptReco_Collisions16_JSON_NoL1T_Tau_Run2016H.txt"))

#================================================================================================ 
# PhaseIIFall17D
#================================================================================================ 
das = "https://cmsweb.cern.ch/das/request?view=list&limit=100&instance=prod%2Fglobal&input=dataset%3D%2F*%2FPhaseIIFall17D*%2FGEN-SIM-DIGI-RAW"
datasetsSingleNu = []
datasetsSingleNu.append(Dataset('/SingleNeutrino/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="SingleNeutrino_L1TPU140")) 
datasetsSingleNu.append(Dataset('/SingleNeutrino/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="SingleNeutrino_L1TPU200")) 

das = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=%2FRelValSingleTauFlatPt2To100_pythia8%2FCMSSW_9_3_7*D17*%2FGEN-SIM-DIGI-RAW"
datasetsSingleTau = []
datasetsSingleTau.append(Dataset('/RelValSingleTauFlatPt2To100_pythia8/CMSSW_9_3_7-93X_upgrade2023_realistic_v5_2023D17noPU-v2/GEN-SIM-DIGI-RAW'        , dataVersion="93Xmc" , dasQuery=das, alias="SingleTau_L1TnoPU"))
datasetsSingleTau.append(Dataset('/RelValSingleTauFlatPt2To100_pythia8/CMSSW_9_3_7-PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/GEN-SIM-DIGI-RAW', dataVersion= "93Xmc", dasQuery=das, alias="SingleTau_L1TPU200"))

das = "https://cmsweb.cern.ch/das/request?view=list&limit=100&instance=prod%2Fglobal&input=dataset%3D%2F*%2FPhaseIIFall17D*%2FGEN-SIM-DIGI-RAW"
datasetsTT = []
datasetsTT.append(Dataset('/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW' , dataVersion="93Xmc", dasQuery=das, alias="TT_TuneCUETP8M2T4_14TeV_L1TnoPU"))
datasetsTT.append(Dataset('/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="TT_TuneCUETP8M2T4_14TeV_L1TPU140"))
datasetsTT.append(Dataset('/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="TT_TuneCUETP8M2T4_14TeV_L1TPU200"))

das = "https://cmsweb.cern.ch/das/request?view=list&limit=100&instance=prod%2Fglobal&input=dataset%3D%2F*%2FPhaseIIFall17D*%2FGEN-SIM-DIGI-RAW"
datasetsH2tautau = []
datasetsH2tautau.append(Dataset('/GluGluHToTauTau_M125_14TeV_powheg_pythia8/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW' , dataVersion="93Xmc", dasQuery=das, alias="GluGluHToTauTau_14TeV_L1TnoPU"))
datasetsH2tautau.append(Dataset('/GluGluHToTauTau_M125_14TeV_powheg_pythia8/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="GluGluHToTauTau_14TeV_L1TPU140"))
datasetsH2tautau.append(Dataset('/GluGluHToTauTau_M125_14TeV_powheg_pythia8/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="GluGluHToTauTau_14TeV_L1TPU200"))

das = "https://cmsweb.cern.ch/das/request?view=list&limit=100&instance=prod%2Fglobal&input=dataset%3D%2F*%2FPhaseIIFall17D*%2FGEN-SIM-DIGI-RAW"
datasetsHPlus = []
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW'  , dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs200_14TeV_L1TnoPU"))
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW' , dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs200_14TeV_L1TPU140"))
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW' , dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs200_14TeV_L1TPU200"))
#datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW'  , dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs500_14TeV_L1TnoPU"))
#datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW' , dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs500_14TeV_L1TPU140"))
#datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW' , dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs500_14TeV_L1TPU200"))
#datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW' , dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs1000_14TeV_L1TnoPU"))
#datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs1000_14TeV_L1TPU140"))
#datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v2/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="ChargedHiggs1000_14TeV_L1TPU200"))

datasetsSingleE = []
datasetsSingleE.append(Dataset('/SingleE_FlatPt-2to100/PhaseIIFall17D-L1TnoPU_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW',dataVersion="93Xmc", dasQuery=das, alias="SingleE_L1TnoPU"))
datasetsSingleE.append(Dataset('/SingleE_FlatPt-2to100/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="SingleE_L1TPU140"))
datasetsSingleE.append(Dataset('/SingleE_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW', dataVersion="93Xmc", dasQuery=das, alias="SingleE_L1TPU200"))

#================================================================================================ 
# Dataset Grouping
#================================================================================================ 
PhaseIIFall17D = []
PhaseIIFall17D.extend(datasetsSingleNu)
PhaseIIFall17D.extend(datasetsTT)
PhaseIIFall17D.extend(datasetsH2tautau)
PhaseIIFall17D.extend(datasetsHPlus)
PhaseIIFall17D.extend(datasetsSingleTau)
PhaseIIFall17D.extend(datasetsSingleE)

AllDatasets = []
AllDatasets += PhaseIIFall17D

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
        analyses = ["CaloTk", "TkCalo", "PFTau", "CaloPix", "NEW", "All", "TDR2019"]
        if self.analysis not in analyses:
            raise Exception("Unknown analysis \"%s\". Please select one of the following: \"%s" % (self.analysis, "\", \"".join(analyses) + "\".") )
        self.GroupDict["CaloTk"]  = PhaseIIFall17D
        self.GroupDict["TkCalo"]  = PhaseIIFall17D
        self.GroupDict["PFTau"]   = PhaseIIFall17D
        self.GroupDict["CaloPix"] = PhaseIIFall17D
	self.GroupDict["NEW"]	  = PhaseIIFall17D
	self.GroupDict["TDR2019"] = PhaseIIFall17D
        self.GroupDict["All"]     = AllDatasets
        return


    def GetDatasetList(self):
        '''
        Return the dataset list according to the analysis name. 
        Uses pre-defined dictionary mapping: analysis->dataset list
        '''
        return self.GroupDict[self.analysis]

    def GetDatasetListAlt(self):
        '''
        Return the dataset list according to the analysis name. 
        Uses pre-defined dictionary mapping: analysis->dataset list
        '''
        return self.GroupDict["TDR2019"]


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

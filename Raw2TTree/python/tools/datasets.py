#lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt" # ICHEP dataset 271036-276811
lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"

#================================================================================================ 
# Class Definition
#================================================================================================ 
import os

class Dataset:
    def __init__(self, url, dbs="global", dataVersion="80Xmc", dasQuery="", lumiMask=lumiMask, name=None):
        self.URL = url
        self.DBS = dbs
        self.dataVersion = dataVersion
        if not os.path.dirname(lumiMask):
            lumiMask = os.path.join(os.environ['CMSSW_BASE'],"src/HiggsAnalysis/MiniAOD2TTree/data", lumiMask)
        self.lumiMask = lumiMask
        self.DASquery = dasQuery
	if name != None:
            self.name = name
        else:
            self.name = url
            
    def isData(self):
        if "data" in self.dataVersion:
            return True
        return False

    def getName(self):
	return self.name

#================================================================================================ 
# PhaseIISpring17D
#================================================================================================ 
das_PhaseIISpring17D = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset%3D%2F*%2F*PhaseIISpring17D*%2F*"

datasetsBs = []
datasetsBs.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsBs.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsBs.append(Dataset('/BsToPhiPhi_SoftQCDnonD_TuneCUEP8M1_14TeV_Pythia8/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )

datasetsSingleE = []
datasetsSingleE.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW' , dataVersion="90Xmc") ) 
datasetsSingleE.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSingleE.append(Dataset('/SingleE_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsSingleMu = []
datasetsSingleMu.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSingleMu.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSingleMu.append(Dataset('/SingleMu_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsSingleNu = []
datasetsSingleNu.append(Dataset('/SingleNeutrino/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsSingleNu.append(Dataset('/SingleNeutrino/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsSinglePhoton = []
datasetsSinglePhoton.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSinglePhoton.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSinglePhoton.append(Dataset('/SinglePhoton_FlatPt-8to150/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsSinglePion0 = []
datasetsSinglePion0.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSinglePion0.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSinglePion0.append(Dataset('/SinglePion0_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsSinglePion = []
datasetsSinglePion.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSinglePion.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsSinglePion.append(Dataset('/SinglePion_FlatPt-8to100/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsSingleTau = []
datasetsSingleTau.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsSingleTau.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsSingleTau.append(Dataset('/SingleTau_FlatPt-8to150/PhaseIISpring17D-PU200_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )

datasetsSingleTau1pr = []
datasetsSingleTau1pr.append(Dataset('/SingleTauOneProngFlatPt10To100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsSingleTau1pr.append(Dataset('/SingleTauOneProngFlatPt10To100/PhaseIISpring17D-PU140_MB100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsSingleTau1pr.append(Dataset('/SingleTauOneProngFlatPt10To100/PhaseIISpring17D-PU200_MB100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )

datasetsSingleTau3pr = []
datasetsSingleTau3pr.append(Dataset('/TauThreeProngsEnriched/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsSingleTau3pr.append(Dataset('/TauThreeProngsEnriched/PhaseIISpring17D-PU140_MB100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsSingleTau3pr.append(Dataset('/TauThreeProngsEnriched/PhaseIISpring17D-PU200_MB100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )

datasetsTTJpsiFilter = []
datasetsTTJpsiFilter.append(Dataset('/TT_JpsiFilter_TuneCUETP8M1_mtop166_5_14TeV-powheg-tauola-pythia8/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsTTJpsiFilter.append(Dataset('/TT_JpsiFilter_TuneCUETP8M1_mtop166_5_14TeV-powheg-tauola-pythia8/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsTT = []
datasetsTT.append(Dataset('/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsTT.append(Dataset('/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIISpring17D-PU140_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsTT.append(Dataset('/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIISpring17D-PU200_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )

datasetsTT_pilot = []
datasetsTT.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsTT.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )
datasetsTT.append(Dataset('/TT_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIISpring17D-PU200_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") )

datasetsTTbar = []
datasetsTTbar.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-NoPU_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsTTbar.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-PU140_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsTTbar.append(Dataset('/TTbar_14TeV_TuneCUETP8M1_PhaseIIFall16/PhaseIISpring17D-PU200_pilot_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsQCD = []
datasetsQCD.append(Dataset('/QCD_Pt-15to3000_Tune4C_14TeV_pythia8/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsQCD.append(Dataset('/QCD_Pt-15to3000_Tune4C_14TeV_pythia8/PhaseIISpring17D-PU140_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsQCD.append(Dataset('/QCD_Pt-15to3000_Tune4C_14TeV_pythia8/PhaseIISpring17D-PU200_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

datasetsHPlus = []
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV/PhaseIISpring17D-PU140_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs200_14TeV/PhaseIISpring17D-PU200_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV/PhaseIISpring17D-PU140_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs500_14TeV/PhaseIISpring17D-PU200_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV/PhaseIISpring17D-PU140_100M_90X_upgrade2023_realistic_v9-v1/GEN-DIGI-SIM-RAW', dataVersion="90Xmc") ) 
datasetsHPlus.append(Dataset('/PYTHIA_Tauola_TB_ChargedHiggs1000_14TeV/PhaseIISpring17D-PU200_100M_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW', dataVersion="90Xmc") ) 

#================================================================================================ 
# Dataset Grouping
#================================================================================================ 
Test92XDatasets = []
# Test92XDatasets.extend(datasetsSingleE)
Test92XDatasets.extend(datasetsSinglePion)
Test92XDatasets.extend(datasetsSinglePion0)
Test92XDatasets.extend(datasetsSingleTau)
Test92XDatasets.extend(datasetsSingleTau1pr)
Test92XDatasets.extend(datasetsSingleTau3pr)
Test92XDatasets.extend(datasetsSingleNu)
#Test92XDatasets.extend(datasetsTTJpsiFilter)
Test92XDatasets.extend(datasetsTT)
#Test92XDatasets.extend(datasetsHPlus)
#Test92XDatasets.extend(datasetsTT_pilot)
#Test92XDatasets.extend(datasetsTTbar)
#Test92XDatasets.extend(datasetsQCD)

L1TauDatasets = []
L1TauDatasets.extend(datasetsSingleE)
L1TauDatasets.extend(datasetsSinglePion)
L1TauDatasets.extend(datasetsSinglePion0)
L1TauDatasets.extend(datasetsSingleTau)
L1TauDatasets.extend(datasetsSingleTau1pr)
L1TauDatasets.extend(datasetsSingleTau3pr)
L1TauDatasets.extend(datasetsSingleNu)
#L1TauDatasets.extend(datasetsTTJpsiFilter)
L1TauDatasets.extend(datasetsTT)
L1TauDatasets.extend(datasetsHPlus)
#L1TauDatasets.extend(datasetsTT_pilot)
#L1TauDatasets.extend(datasetsTTbar)
#L1TauDatasets.extend(datasetsQCD)

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
        analyses = ["CaloTk", "TkCalo", "PFTau", "CaloPix", "All", "Test_92X"]
        if self.analysis not in analyses:
            raise Exception("Unknown analysis \"%s\". Please select one of the following: \"%s" % (self.analysis, "\", \"".join(analyses) + "\".") )
        self.GroupDict["CaloTk"]   = L1TauDatasets
        self.GroupDict["TkCalo"]   = L1TauDatasets
        self.GroupDict["PFTau"]    = L1TauDatasets
        self.GroupDict["CaloPix"]  = L1TauDatasets
        self.GroupDict["All"]      = AllDatasets
        self.GroupDict["Test_92X"] = Test92XDatasets
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

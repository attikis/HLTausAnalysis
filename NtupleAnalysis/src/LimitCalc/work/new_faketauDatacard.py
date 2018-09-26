import HiggsAnalysis.NtupleAnalysis.tools.systematics as systematics

DataCardName    = 'Default_8TeV'
#Path = "/home/wendland/data/v533/2014-01-29-noMetSF-withL1ETMfix"#2014-01-29-noMetSF-withL1ETMfix"
#Path = "/home/wendland/data/v533/2014_02_14_v3_etacorrected"
#Path = "/home/wendland/data/v533/2014-03-20"
#Path = "/home/wendland/data/v533/2014-03-20_expCtrlPlots"
#Path = "/home/wendland/data/v533/2014-04-14_nominal_norm5GeVLRB"
#Path = "/home/wendland/data/v533/2014-03-20_optTau60Met80_mt20gev"
#Path = "/home/wendland/data/v533/2014-03-20_METprecut30"
#Path = "/home/wendland/data/v533/2014_03_12_metphicorrected"
#Path = "/home/wendland/data/v533/2014_02_14_v3_decaymode1"
#Path            = '/home/wendland/data/v445/met50rtaunprongs'
#Path            = '/mnt/flustre/slehti/hplusAnalysis/QCDInverted/CMSSW_4_4_5/src/HiggsAnalysis/NtupleAnalysis/test/datacardGenerator/TESTDATA/'
Path = "/mnt/flustre/epekkari/FakeTauDatacard"
#Path = "/mnt/flustre/epekkari/NominalDatacard"

LightMassPoints      = [80,90,100,120,140,150,155,160]
#LightMassPoints      = [80,120,160]
#LightMassPoints      = [120]
#LightMassPoints      = []
HeavyMassPoints      = [180,190,200,220,250,300,400,500,600] # mass points 400-600 are not available for 2011 branch
#HeavyMassPoints      = [180,220,300,600]
#HeavyMassPoints      = [180]
#HeavyMassPoints      = []
MassPoints = LightMassPoints[:]+HeavyMassPoints[:]

BlindAnalysis   = False
OptionBlindThreshold = None # If signal exceeds this fraction of expected events, data is blinded; set to None to disable

# Rate counter definitions
SignalRateCounter = "Selected events"
FakeRateCounter = "EWKfaketaus:SelectedEvents"

# Options
OptionMassShape = "TransverseMass"
#OptionMassShape = "FullMass"
#OptionMassShape = "TransverseAndFullMass2D" #FIXME not yet supported!!!

# Choose source of EWK+tt genuine tau background
OptionGenuineTauBackgroundSource = "DataDriven"                          # State-of-the-art: embedded data used (use for optimization and results)
#OptionGenuineTauBackgroundSource = "MC_FakeAndGenuineTauNotSeparated" # MC used, fake taus are not separated from genuine taus
#OptionGenuineTauBackgroundSource = "MC_FullSystematics"               # MC used, fake and genuine taus separated (use for embedding closure test)
#OptionGenuineTauBackgroundSource = "MC_RealisticProjection"            # MC used, fake and genuine taus separated (can be used for optimization)

OptionRealisticEmbeddingWithMC = True # Only relevant for OptionReplaceEmbeddingByMC==True
OptionTreatTriggerUncertaintiesAsAsymmetric = True # Set to true, if you produced multicrabs with doAsymmetricTriggerUncertainties=True
OptionTreatTauIDAndMisIDSystematicsAsShapes = True # Set to true, if you produced multicrabs with doTauIDandMisIDSystematicsAsShapes=True
OptionIncludeSystematics = True # Set to true if you produced multicrabs with doSystematics=True

OptionPurgeReservedLines = True # Makes limit running faster, but cannot combine leptonic datacards
#OptionDoControlPlots = True
OptionDoMergeFakeTauColumns = False
OptionDoControlPlots = True
OptionCombineSingleColumnUncertainties = True # Makes limit running faster
OptionCtrlPlotsAtMt = True
OptionDisplayEventYieldSummary = True
OptionNumberOfDecimalsInSummaries = 1
OptionRemoveHHDataGroup = True
OptionLimitOnSigmaBr = False # Is automatically set to true for heavy H+
OptionDoTBbarForHeavy = False # NOTE: usable only for 2012

# For projections
trg_MET_dataeffScaleFactor = None # Default is None, i.e. 1.0

# Options for reports and article
OptionBr = 0.01  # Br(t->bH+)
OptionSqrtS = 8 # sqrt(s)

# Tolerance for throwing error on luminosity difference (0.01 = 1 percent agreement is required)
ToleranceForLuminosityDifference = 0.05
# Tolerance for almost zero rate (columns with smaller rate are suppressed)
ToleranceForMinimumRate = 1.5
# Minimum stat. uncertainty to set to bins with zero events
MinimumStatUncertainty = 0.5


# Shape histogram definitions
SignalShapeHisto = None
FakeShapeHisto = None
ShapeHistogramsDimensions = None

if OptionMassShape == "TransverseMass":
    SignalShapeHisto = "shapeTransverseMass"
    FakeShapeTTbarHisto = "shapeEWKFakeTausTransverseMass"
    FakeShapeOtherHisto = "shapeProbabilisticBtagEWKFakeTausTransverseMass"
elif OptionMassShape == "FullMass":
    SignalShapeHisto = "shapeInvariantMass"
    FakeShapeOtherHisto = "shapeEWKFakeTausInvariantMass"
    FakeShapeTTbarHisto = FakeShapeOtherHisto
elif OptionMassShape == "TransverseAndFullMass2D": # FIXME: preparing to add support, not yet working
    SignalShapeHisto = "shapetransverseAndFullMass2D" # FIXME: Not yet implemented to signal analysis, etc.
    FakeShapeOtherHisto = "shapeEWKFakeTausTransverseAndFullMass2D" # FIXME: Not yet implemented to signal analysis, etc.
    FakeShapeTTbarHisto = FakeShapeOtherHisto
ShapeHistogramsDimensions = systematics.getBinningForPlot(SignalShapeHisto)

DataCardName += "_"+OptionMassShape.replace("TransverseMass","mT").replace("FullMass","invMass")

##############################################################################
# Observation definition (how to retrieve number of observed events)
#
from HiggsAnalysis.LimitCalc.InputClasses import ObservationInput
Observation = ObservationInput(datasetDefinition="Data",
                               shapeHisto=SignalShapeHisto)
#Observation.setPaths(signalPath,signalDataPaths)

##############################################################################
# Systematics lists
myTrgShapeSystematics = []
if OptionTreatTriggerUncertaintiesAsAsymmetric:
    myTrgShapeSystematics = ["trg_tau_dataeff","trg_tau_MCeff","trg_L1ETM_dataeff","trg_L1ETM_MCeff"] # Variation done separately for data and MC efficiencies
else:
    myTrgShapeSystematics = ["trg_tau","trg_L1ETM"] # Variation of trg scale factors

myTauIDShapeSystematics = []
if OptionTreatTauIDAndMisIDSystematicsAsShapes:
    myTauIDShapeSystematics = ["tau_ID_shape","tau_ID_eToTauBarrel_shape","tau_ID_eToTauEndcap_shape","tau_ID_muToTau_shape","tau_ID_jetToTau_shape"] # tau ID and mis-ID systematics done with shape variation
else:
    myTauIDShapeSystematics = ["tau_ID"] # tau ID and mis-ID systematics done with constants

myShapeSystematics = []
myShapeSystematics.extend(myTrgShapeSystematics)
myShapeSystematics.extend(myTauIDShapeSystematics)
myShapeSystematics.extend(["ES_taus","ES_jets","JER","ES_METunclustered","pileup"]) # btag is not added, because it has the tag and mistag categories

myEmbeddingMETUncert = "trg_L1ETM"
if OptionTreatTriggerUncertaintiesAsAsymmetric:
    myEmbeddingMETUncert += "_dataeff"
    myEmbeddingShapeSystematics = ["trg_tau_dataeff",myEmbeddingMETUncert,"trg_muon_dataeff","ES_taus","Emb_mu_ID","Emb_WtauTomu"]
else:
    myEmbeddingShapeSystematics = ["trg_tau",myEmbeddingMETUncert,"trg_muon","ES_taus","Emb_mu_ID","Emb_WtauTomu"]
# Add tau ID uncert. to embedding either as a shape or as a constant
if "tau_ID_shape" in myTauIDShapeSystematics:
    myEmbeddingShapeSystematics.append("tau_ID_shape")
else:
    myEmbeddingShapeSystematics.append("tau_ID")

myFakeShapeSystematics = []
for item in myShapeSystematics:
    if item == "tau_ID":
        myFakeShapeSystematics.append("tau_misID")
    else:
        if not item == "tau_ID_shape":
            myFakeShapeSystematics.append(item)

##############################################################################
# DataGroup (i.e. columns in datacard) definitions
#
from HiggsAnalysis.LimitCalc.InputClasses import DataGroup
DataGroups = []
EmbeddingIdList = []
EWKFakeIdList = []

signalTemplate = DataGroup(datasetType="Signal",
                           shapeHisto=SignalShapeHisto)

for mass in LightMassPoints:
    myMassList = [mass]
    if not OptionRemoveHHDataGroup:
        hhx = signalTemplate.clone()
        hhx.setLabel("HH"+str(mass)+"_a")
        hhx.setLandSProcess(-1)
        hhx.setValidMassPoints(myMassList)
        hhx.setNuisances(myShapeSystematics[:]+["e_mu_veto","b_tag","xsect_tt_8TeV","lumi"])
        hhx.setDatasetDefinition("TTToHplusBHminusB_M"+str(mass))
        DataGroups.append(hhx)

    hwx = signalTemplate.clone()
    hwx.setLabel("HW"+str(mass)+"_a")
    hwx.setLandSProcess(0)
    hwx.setValidMassPoints(myMassList)
    hwx.setNuisances(myShapeSystematics[:]+["e_mu_veto","b_tag","xsect_tt_8TeV","lumi"])
    hwx.setDatasetDefinition("TTToHplusBWB_M"+str(mass))
    DataGroups.append(hwx)

for mass in HeavyMassPoints:
    myMassList = [mass]
    hx = signalTemplate.clone()
    hx.setLabel("Hp"+str(mass)+"_a")
    hx.setLandSProcess(0)
    hx.setValidMassPoints(myMassList)
    hx.setNuisances(myShapeSystematics[:]+["e_mu_veto","b_tag","lumi"])
    if not OptionDoTBbarForHeavy:
        hx.setDatasetDefinition("HplusTB_M"+str(mass))
    else:
        hx.setDatasetDefinition("HplusToTBbar_M"+str(mass))
    DataGroups.append(hx)

myQCDShapeSystematics = myShapeSystematics[:]
#for i in range(0,len(myQCDShapeSystematics)):
    #if myQCDShapeSystematics[i].startswith("trg_CaloMET") and not "forQCD" in myQCDShapeSystematics[i]:
    #    myQCDShapeSystematics[i] = myQCDShapeSystematics[i]+"_forQCD"

myQCDFact = DataGroup(
    label        = "QCDfact",
    landsProcess = 3,
    validMassPoints = MassPoints,
    datasetType  = "QCD factorised",
    datasetDefinition = "QCDfactorisedmt",
    nuisances    = myQCDShapeSystematics[:]+["b_tag","top_pt","QCD_metshape","xsect_tt_8TeV_forQCD"],
    shapeHisto   = SignalShapeHisto,
)

myQCDInv = DataGroup(
    label        = "QCDinv",
    landsProcess = 3,
    validMassPoints = MassPoints,
    datasetType  = "QCD inverted",
    datasetDefinition = "QCDinvertedmt",
    nuisances    = myQCDShapeSystematics[:]+["b_tag","top_pt","QCD_metshape","xsect_tt_8TeV_forQCD","QCDinvTemplateFit"],
    shapeHisto   = SignalShapeHisto,
)

if OptionMassShape == "TransverseMass":
    myQCDFact.setDatasetDefinition("QCDfactorisedmt")
    myQCDInv.setDatasetDefinition("QCDinvertedmt")
elif OptionMassShape == "FullMass":
    myQCDFact.setDatasetDefinition("QCDfactorisedinvmass")
    myQCDInv.setDatasetDefinition("QCDinvertedinvmass")

DataGroups.append(myQCDFact)
DataGroups.append(myQCDInv)

if OptionGenuineTauBackgroundSource == "DataDriven":
    # EWK + ttbar with genuine taus
    EmbeddingIdList = [3]
    DataGroups.append(DataGroup(
        label        = "EWK_Tau",
        landsProcess = 1,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        #datasetDefinition   = ["SingleMu"],
        datasetDefinition   = "Data",
        validMassPoints = MassPoints,
        #additionalNormalisation = 0.25, # not needed anymore
        nuisances    = myEmbeddingShapeSystematics[:]+["Emb_QCDcontam","Emb_hybridCaloMET"]
        #nuisances    = ["trg_tau_embedding","tau_ID","ES_taus","Emb_QCDcontam","Emb_WtauTomu","Emb_musel_ditau_mutrg","stat_Emb"]
    ))
    EWKFakeIdList = [] #this is used

elif OptionGenuineTauBackgroundSource == "MC_FullSystematics" or OptionGenuineTauBackgroundSource == "MC_RealisticProjection":
    # Mimic embedding with MC analysis (introduces double counting of EWK fakes, but that should be small effect)
    EmbeddingIdList = [4]
    myEmbeddingShapeSystematics = []
    if OptionGenuineTauBackgroundSource == "MC_RealisticProjection":
        # Mimic with uncertainties the outcome of data-driven embedding
        if OptionTreatTriggerUncertaintiesAsAsymmetric:
            myEmbeddingShapeSystematics.append("trg_tau_dataeff")
            myEmbeddingShapeSystematics.append("trg_L1ETM_dataeff")
        else:
            myEmbeddingShapeSystematics.append("trg_tau")
            myEmbeddingShapeSystematics.append("trg_L1ETM")
        myEmbeddingShapeSystematics.append("ES_taus")
        if OptionTreatTauIDAndMisIDSystematicsAsShapes:
            myEmbeddingShapeSystematics.append("tau_ID_shape")
            myFakeShapeSystematics.append("tau_ID_shape")
        else:
            myEmbeddingShapeSystematics.append("tau_ID")
            myFakeShapeSystematics.append("tau_ID")
        myEmbeddingShapeSystematics.extend(["Emb_QCDcontam","Emb_hybridCaloMET","Emb_rest"])
    elif OptionGenuineTauBackgroundSource == "MC_FullSystematics":
        # Use full MC systematics; approximate xsect uncertainty with ttbar xsect unsertainty
        myEmbeddingShapeSystematics = myShapeSystematics[:]+["top_pt","e_mu_veto","b_tag","xsect_tt_8TeV","lumi"]
    DataGroups.append(DataGroup(
        label        = "pseudo_emb_TTJets_MC",
        landsProcess = 1,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "TTJets",
        validMassPoints = MassPoints,
        nuisances    = myEmbeddingShapeSystematics[:],
    ))
    DataGroups.append(DataGroup(
        label        = "pseudo_emb_Wjets_MC",
        landsProcess = None,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "WJets",
        validMassPoints = MassPoints,
        nuisances    = myEmbeddingShapeSystematics,
    ))
    DataGroups.append(DataGroup(
        label        = "pseudo_emb_t_MC",
        landsProcess = None,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "SingleTop",
        validMassPoints = MassPoints,
        nuisances    = myEmbeddingShapeSystematics,
    ))
    DataGroups.append(DataGroup(
        label        = "pseudo_emb_DY_MC",
        landsProcess = None,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition   = "DYJetsToLL",
        validMassPoints = MassPoints,
        nuisances    = myEmbeddingShapeSystematics,
    ))
    DataGroups.append(DataGroup(
        label        = "pseudo_emb_VV_MC",
        landsProcess = None,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition   = "Diboson",
        validMassPoints = MassPoints,
        nuisances    = myEmbeddingShapeSystematics,
    ))
    EWKFakeIdList = [1,5,6,7,8]
    DataGroups.append(DataGroup(
        label        = "tt_EWK_faketau",
        landsProcess = 1,
        shapeHisto   = FakeShapeTTbarHisto,
        datasetType  = "EWKfake",
        datasetDefinition = "TTJets",
        validMassPoints = MassPoints,
        nuisances    = myFakeShapeSystematics[:]+["e_mu_veto_fakes","b_tag_fakes","top_pt","xsect_tt_8TeV","lumi"],
    ))
    DataGroups.append(DataGroup(
        label        = "W_EWK_faketau",
        landsProcess = 5,
        shapeHisto   = FakeShapeOtherHisto,
        datasetType  = "EWKfake",
        datasetDefinition = "WJets",
        validMassPoints = MassPoints,
        nuisances    = myFakeShapeSystematics[:]+["e_mu_veto_fakes","b_mistag_fakes","xsect_Wjets","lumi","probBtag"],
    ))
    DataGroups.append(DataGroup(
        label        = "t_EWK_faketau",
        landsProcess = 6,
        shapeHisto   = FakeShapeOtherHisto,
        datasetType  = "EWKfake",
        datasetDefinition = "SingleTop",
        validMassPoints = MassPoints,
        nuisances    = myFakeShapeSystematics[:]+["e_mu_veto_fakes","b_tag_fakes","xsect_singleTop","lumi","probBtag"],
    ))
    DataGroups.append(DataGroup(
        label        = "DY_EWK_faketau",
        landsProcess = 7,
        shapeHisto   = FakeShapeOtherHisto,
        datasetType  = "EWKfake",
        datasetDefinition   = "DYJetsToLL",
        validMassPoints = MassPoints,
        nuisances    = myFakeShapeSystematics[:]+["e_mu_veto_fakes","b_mistag_fakes","xsect_DYtoll","lumi","probBtag"],
    ))
    DataGroups.append(DataGroup(
        label        = "VV_EWK_faketau",
        landsProcess = 8,
        shapeHisto   = FakeShapeOtherHisto,
        datasetType  = "EWKfake",
        datasetDefinition   = "Diboson",
        validMassPoints = MassPoints,
        nuisances    = myFakeShapeSystematics[:]+["e_mu_veto_fakes","b_mistag_fakes","xsect_VV","lumi","probBtag"],
    ))
elif OptionGenuineTauBackgroundSource == "MC_FakeAndGenuineTauNotSeparated":
    # Replace embedding and fakes with MC
    EmbeddingIdList = [1,4,5,6,7]
    DataGroups.append(DataGroup(
        label        = "ttbar_MC",
        landsProcess = 1,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "TTJets",
        validMassPoints = MassPoints,
        nuisances    = myShapeSystematics[:]+["e_mu_veto","b_tag","top_pt","xsect_tt_8TeV","lumi"],
    ))
    DataGroups.append(DataGroup(
        label        = "Wjets_MC",
        landsProcess = 4,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "WJets",
        validMassPoints = MassPoints,
        nuisances    = myShapeSystematics[:]+["e_mu_veto","b_mistag","xsect_Wjets","lumi"],
    ))
    DataGroups.append(DataGroup(
        label        = "t_MC",
        landsProcess = 5,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "SingleTop",
        validMassPoints = MassPoints,
        nuisances    = myShapeSystematics[:]+["e_mu_veto","b_tag","xsect_singleTop","lumi"],
    ))
    DataGroups.append(DataGroup(
        label        = "DY_MC",
        landsProcess = 6,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "DYJetsToLL",
        validMassPoints = MassPoints,
        nuisances    = myShapeSystematics[:]+["e_mu_veto","b_mistag","xsect_DYtoll","lumi"],
    ))
    DataGroups.append(DataGroup(
        label        = "VV_MC",
        landsProcess = 7,
        shapeHisto   = SignalShapeHisto,
        datasetType  = "Embedding",
        datasetDefinition = "Diboson",
        validMassPoints = MassPoints,
        nuisances    = myShapeSystematics[:]+["e_mu_veto","b_mistag","xsect_VV","lumi"],
    ))
else:
    raise Exception("Error: unknown value for flag OptionGenuineTauBackgroundSource!")

# Reserve column 2
DataGroups.append(DataGroup(
    label        = "res.",
    landsProcess = 2,
    datasetType  = "None",
    validMassPoints = MassPoints,
))


##############################################################################
# Definition of nuisance parameters
#
# Note: Remember to include 'stat.' into the label of nuistances of statistical nature
#
from HiggsAnalysis.LimitCalc.InputClasses import Nuisance
ReservedNuisances = []
ReservedNuisances.append(["05", "reserved for leptonic"])
ReservedNuisances.append(["06", "reserved for leptonic"])
ReservedNuisances.append(["08", "reserved for leptonic"])
ReservedNuisances.append(["20", "reserved for leptonic"])
ReservedNuisances.append(["21", "reserved for leptonic"])
ReservedNuisances.append(["23", "reserved for leptonic"])
if OptionPurgeReservedLines:
    ReservedNuisances = []

Nuisances = []

if "trg_tau" in myShapeSystematics:
    Nuisances.append(Nuisance(
        id            = "trg_tau",
        label         = "tau+MET trg tau part",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "TauTrgSF",
    ))
else:
    Nuisances.append(Nuisance(
        id            = "trg_tau_dataeff",
        label         = "tau+MET trg tau part data eff.",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "TauTrgDataEff",
    ))

    Nuisances.append(Nuisance(
        id            = "trg_tau_MCeff",
        label         = "tau+MET trg tau part MC eff.",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "TauTrgMCEff",
    ))

if OptionIncludeSystematics:
    if OptionTreatTriggerUncertaintiesAsAsymmetric:
        Nuisances.append(Nuisance(
            id            = "trg_L1ETM_dataeff",
            label         = "tau+MET trg L1ETM data eff.",
            distr         = "shapeQ",
            function      = "ShapeVariation",
            systVariation = "L1ETMDataEff",
        ))
        Nuisances.append(Nuisance(
            id            = "trg_L1ETM_MCeff",
            label         = "tau+MET trg L1ETM MC eff.",
            distr         = "shapeQ",
            function      = "ShapeVariation",
            systVariation = "L1ETMMCEff",
        ))
    else:
        Nuisances.append(Nuisance(
            id            = "trg_L1ETM",
            label         = "tau+MET trg L1ETM",
            distr         = "shapeQ",
            function      = "ShapeVariation",
            systVariation = "L1ETMSF",
        ))

if OptionGenuineTauBackgroundSource == "DataDriven":
    if OptionTreatTriggerUncertaintiesAsAsymmetric:
        Nuisances.append(Nuisance(
            id            = "trg_muon_dataeff",
            label         = "SingleMu trg data eff.",
            distr         = "shapeQ",
            function      = "ShapeVariation",
            systVariation = "MuonTrgDataEff",
        ))
    else:
        Nuisances.append(Nuisance(
            id            = "trg_muon",
            label         = "SingleMu trg data eff.",
            distr         = "lnN",
            function      = "Constant",
            value         = 0.02,
        ))

if not "tau_ID_shape" in myShapeSystematics:
    Nuisances.append(Nuisance(
        id            = "tau_ID",
        label         = "tau-jet ID (no Rtau)",
        distr         = "lnN",
        function      = "Constant",
        value         = systematics.getTauIDUncertainty(isGenuineTau=True)
    ))

    Nuisances.append(Nuisance(
        id            = "tau_misID",
        label         = "tau-jet mis ID (no Rtau)",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.15, # FIXME
    ))

Nuisances.append(Nuisance(
    id            = "tau_ID_constShape",
    label         = "tau-jet ID (no Rtau)",
    distr         = "shapeQ",
    function      = "ConstantToShape",
    value         = systematics.getTauIDUncertainty(isGenuineTau=True)
))

if "tau_ID_shape" in myShapeSystematics:
    Nuisances.append(Nuisance(
        id            = "tau_ID_shape",
        label         = "tau-jet ID (no Rtau) genuine taus",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "GenuineTau",
    ))

if "tau_ID_eToTauBarrel_shape" in myShapeSystematics:
    Nuisances.append(Nuisance(
        id            = "tau_ID_eToTauBarrel_shape",
        label         = "tau-jet ID (no Rtau) e->tau (barrel)",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "FakeTauBarrelElectron",
    ))

if "tau_ID_eToTauEndcap_shape" in myShapeSystematics:
    Nuisances.append(Nuisance(
        id            = "tau_ID_eToTauEndcap_shape",
        label         = "tau-jet ID (no Rtau) e->tau (endcap)",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "FakeTauEndcapElectron",
    ))

if "tau_ID_muToTau_shape" in myShapeSystematics:
    Nuisances.append(Nuisance(
        id            = "tau_ID_muToTau_shape",
        label         = "tau-jet ID (no Rtau) mu->tau",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "FakeTauMuon",
    ))

if "tau_ID_jetToTau_shape" in myShapeSystematics:
    Nuisances.append(Nuisance(
        id            = "tau_ID_jetToTau_shape",
        label         = "tau-jet ID (no Rtau) jet->tau",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "FakeTauJet",
    ))

if OptionIncludeSystematics:
    Nuisances.append(Nuisance(
        id            = "ES_taus",
        label         = "TES bin-by-bin uncertainty",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "TES",
    ))
    Nuisances.append(Nuisance(
        id            = "ES_jets",
        label         = "JES bin-by-bin uncertainty",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "JES",
    ))
    Nuisances.append(Nuisance(
        id            = "JER",
        label         = "Jet energy resolution",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "JER",
    ))
    Nuisances.append(Nuisance(
        id            = "ES_METunclustered",
        label         = "MET unclustered scale bin-by-bin uncertainty",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "MET",
    ))
    Nuisances.append(Nuisance(
        id            = "b_tag",
        label         = "btagging",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "BTagSF",
    ))
    Nuisances.append(Nuisance(
        id            = "b_tag_fakes",
        label         = "btagging for EWK fake taus",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "BTagSF",
    ))
    Nuisances.append(Nuisance(
        id            = "b_mistag",
        label         = "mistagging",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "BTagSF",
    ))
    Nuisances.append(Nuisance(
        id            = "b_mistag_fakes",
        label         = "mistagging EWK fake taus",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "BTagSF",
    ))
    Nuisances.append(Nuisance(
        id            = "top_pt",
        label         = "top pT reweighting",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "TopPtWeight",
    ))
else:
    Nuisances.append(Nuisance(
        id            = "ES_taus",
        label         = "NON-EXACT VALUE for JES/JER/MET/Rtau effect on mT shape",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.03,
    ))
    Nuisances.append(Nuisance(
        id            = "ES_jets",
        label         = "NON-EXACT VALUE for JES/JER/MET/Rtau effect on mT shape",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.03,
    ))
    Nuisances.append(Nuisance(
        id            = "ES_METunclustered",
        label         = "NON-EXACT VALUE for JES/JER/MET/Rtau effect on mT shape",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.01,
    ))
    Nuisances.append(Nuisance(
        id            = "JER",
        label         = "NON-EXACT VALUE for Jet energy resolution",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.02,
    ))
    Nuisances.append(Nuisance(
        id            = "b_tag",
        label         = "NON-EXACT VALUE for btagging",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05,
    ))
    Nuisances.append(Nuisance(
        id            = "b_tag_fakes",
        label         = "NON-EXACT VALUE for btagging for EWK fake taus",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05,
    ))
    Nuisances.append(Nuisance(
        id            = "b_mistag",
        label         = "NON-EXACT VALUE for mistagging",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05,
    ))
    Nuisances.append(Nuisance(
        id            = "b_mistag_fakes",
        label         = "NON-EXACT VALUE for mistagging EWK fake taus",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05,
    ))
    Nuisances.append(Nuisance(
        id            = "top_pt",
        label         = "NON-EXACT VALUE for top pT reweighting",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.15,
    ))

if OptionGenuineTauBackgroundSource == "DataDriven":
    Nuisances.append(Nuisance(
        id            = "Emb_mu_ID",
        label         = "Muon ID for embedding",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "MuonIdDataEff",
    ))

Nuisances.append(Nuisance(
    id            = "e_mu_veto",
    label         = "lepton veto",
    distr         = "lnN",
    function      = "Ratio",
    numerator     = "muon veto", # main counter name after electron and muon veto
    denominator   = "tau trigger scale factor", # main counter name before electron and muon veto
    scaling       = 0.02
))

Nuisances.append(Nuisance(
    id            = "e_mu_veto_fakes",
    label         = "lepton veto",
    distr         = "lnN",
    function      = "Ratio",
    numerator     = "EWKfaketaus:muon veto", # main counter name after electron and muon veto
    denominator   = "EWKfaketaus:taus == 1", # main counter name before electron and muon veto # the name is misleading, it is actually after tau trg scale factor
    scaling       = 0.02
))

Nuisances.append(Nuisance(
    id            = "QCD_metshape",
    label         = "QCD met shape syst.",
    distr         = "shapeQ",
    function      = "QCDShapeVariation",
    systVariation = "QCDNormSource",
))

if OptionGenuineTauBackgroundSource == "DataDriven" or OptionGenuineTauBackgroundSource == "MC_RealisticProjection":
    Nuisances.append(Nuisance(
        id            = "Emb_QCDcontam",
        label         = "EWK with taus QCD contamination",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.020 #FIXME
    ))
    Nuisances.append(Nuisance(
        id            = "Emb_hybridCaloMET",
        label         = "EWK with taus hybrid calo MET and L1ETM",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.12 #FIXME
    ))


if OptionGenuineTauBackgroundSource == "DataDriven":
    if "Emb_WtauTomu" in myEmbeddingShapeSystematics:
        Nuisances.append(Nuisance(
            id            = "Emb_WtauTomu",
            label         = "EWK with taus W->tau->mu",
            distr         = "shapeQ",
            function      = "ShapeVariation",
            systVariation = "WTauMu",
        ))
    else:
        Nuisances.append(Nuisance(
            id            = "Emb_WtauTomu",
            label         = "EWK with taus W->tau->mu",
            distr         = "lnN",
            function      = "Constant",
            value         = 0.007
        ))

if OptionGenuineTauBackgroundSource == "MC_RealisticProjection":
    Nuisances.append(Nuisance(
        id            = "Emb_rest",
        label         = "EWK with taus W->tau->mu",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.03
    ))

#Nuisances.append(Nuisance(
    #id            = "Emb_musel_ditau_mutrg",
    #label         = "EWK with taus muon selection+ditau+mu trg",
    #distr         = "lnN",
    #function      = "Constant",
    #value         = 0.031 #FIXME
#))

Nuisances.append(Nuisance(
    id            = "xsect_tt_8TeV",
    label         = "ttbar cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = systematics.getCrossSectionUncertainty("TTJets").getUncertaintyDown(),
    upperValue    = systematics.getCrossSectionUncertainty("TTJets").getUncertaintyUp(),
))

Nuisances.append(Nuisance(
    id            = "xsect_tt_8TeV_forQCD",
    label         = "ttbar cross section",
    distr         = "lnN",
    function      = "ConstantForQCD",
    value         = systematics.getCrossSectionUncertainty("TTJets").getUncertaintyDown(),
    upperValue    = systematics.getCrossSectionUncertainty("TTJets").getUncertaintyUp(),
))

Nuisances.append(Nuisance(
    id            = "xsect_Wjets",
    label         = "W+jets cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = systematics.getCrossSectionUncertainty("WJets").getUncertaintyDown()
))

Nuisances.append(Nuisance(
    id            = "xsect_singleTop",
    label         = "single top cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = systematics.getCrossSectionUncertainty("SingleTop").getUncertaintyDown()
))

Nuisances.append(Nuisance(
    id            = "xsect_DYtoll",
    label         = "Z->ll cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = systematics.getCrossSectionUncertainty("DYJetsToLL").getUncertaintyDown()
))

Nuisances.append(Nuisance(
    id            = "xsect_VV",
    label         = "diboson cross section",
    distr         = "lnN",
    function      = "Constant",
    value         = systematics.getCrossSectionUncertainty("Diboson").getUncertaintyDown()
))

Nuisances.append(Nuisance(
    id            = "lumi",
    label         = "luminosity",
    distr         = "lnN",
    function      = "Constant",
    value         = systematics.getLuminosityUncertainty()
))

if OptionIncludeSystematics:
    Nuisances.append(Nuisance(
        id            = "pileup",
        label         = "pileup",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "PUWeight",
    ))

    Nuisances.append(Nuisance(
        id            = "pileup_fakes",
        label         = "pileup",
        distr         = "shapeQ",
        function      = "ShapeVariation",
        systVariation = "PUWeight",
    ))
else:
    Nuisances.append(Nuisance(
        id            = "pileup",
        label         = "FAKE pileup",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05
    ))
    Nuisances.append(Nuisance(
        id            = "pileup_fakes",
        label         = "FAKE pileup",
        distr         = "lnN",
        function      = "Constant",
        value         = 0.05
    ))

Nuisances.append(Nuisance(
    id            = "probBtag",
    label         = "probabilistic btag", 
    distr         = "lnN",
    function      = "Constant", 
    value         = 0.5
))

Nuisances.append(Nuisance(
    id            = "QCDinvTemplateFit",
    label         = "QCDInv: fit", 
    distr         = "lnN",
    function      = "Constant", 
    value         = 0.03,
))

MergeNuisances = []
if "tau_ID_constShape" in myEmbeddingShapeSystematics:
    MergeNuisances.append(["tau_ID_shape", "tau_ID_constShape"])
#if OptionTreatTriggerUncertaintiesAsAsymmetric:
#    MergeNuisances.append(["trg_CaloMET_dataeff", "trg_CaloMET_dataeff_forQCD"])
#    MergeNuisances.append(["trg_CaloMET_MCeff", "trg_CaloMET_MCeff_forQCD"])
#else:
#    MergeNuisances.append(["trg_CaloMET", "trg_CaloMET_forQCD"])
#MergeNuisances.append(["ES_taus","ES_taus_fakes","ES_taus_tempForEmbedding"])
#MergeNuisances.append(["ES_jets","ES_jets_fakes"])
#MergeNuisances.append(["JER","JER_fakes"])
#MergeNuisances.append(["ES_METunclustered","ES_METunclustered_fakes"])
MergeNuisances.append(["e_mu_veto","e_mu_veto_fakes"])
MergeNuisances.append(["b_tag","b_tag_fakes"])
MergeNuisances.append(["b_mistag","b_mistag_fakes"])
MergeNuisances.append(["pileup","pileup_fakes"])
MergeNuisances.append(["xsect_tt_8TeV", "xsect_tt_8TeV_forQCD"])


# Control plots
from HiggsAnalysis.LimitCalc.InputClasses import ControlPlotInput
ControlPlots = []

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_pT_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_pT_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_pT_AfterStandardSelections",
    details          = { "xlabel": "Selected #tau p_{T}",
                         "ylabel": "Events/#Deltap_{T}",
                         "divideByBinWidth": True,
                         "unit": "GeV/c",
                         "log": True,
                         "opts": {"ymin": 0.0009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_p_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_p_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_p_AfterStandardSelections",
    details          = { "xlabel": "Selected #tau p",
                         "ylabel": "Events/#Deltap",
                         "divideByBinWidth": True,
                         "unit": "GeV/c",
                         "log": True,
                         "opts": {"ymin": 0.009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_eta_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_eta_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_eta_AfterStandardSelections",
    details          = { "xlabel": "Selected #tau #eta",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "legendPosition": "SW",
                         "opts": {"ymin": 0.009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_phi_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_phi_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_phi_AfterStandardSelections",
    details          = { "xlabel": "Selected #tau #phi",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "^{o}",
                         "log": True,
                         "legendPosition": "SW",
                         "opts": {"ymin": 0.009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_LeadingTrackPt_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_LeadingTrackPt_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_LeadingTrackPt_AfterStandardSelections",
    details          = { "xlabel": "#tau leading track p_{T}",
                         "ylabel": "Events/#Deltap_{T}",
                         "divideByBinWidth": True,
                         "unit": "GeV/c",
                         "log": True,
                         "ratioLegendPosition": "right",
                         "opts": {"ymin": 0.0009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_LeadingTrackP_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_LeadingTrackP_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_LeadingTrackP_AfterStandardSelections",
    details          = { "xlabel": "#tau leading track p",
                         "ylabel": "Events/#Deltap",
                         "divideByBinWidth": True,
                         "unit": "GeV/c",
                         "log": True,
                         "ratioLegendPosition": "right",
                         "opts": {"ymin": 0.0009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_Rtau_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_Rtau_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_Rtau_AfterStandardSelections",
    details          = { "xlabel": "Selected #tau R_{#tau}",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "legendPosition": "SE",
                         "opts": {"ymin": 0.009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "SelectedTau_DecayMode_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "SelectedTau_DecayMode_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "SelectedTau_DecayMode_AfterStandardSelections",
    details          = { "xlabel": "Selected #tau Decay mode",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "ratioLegendPosition": "right",
                         "opts": {"ymin": 0.9} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "Njets_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "Njets_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "Njets_AfterStandardSelections",
    details          = { "xlabel": "Number of selected jets",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "opts": {"ymin": 0.9} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "#tau_{h}+#geq3j", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "JetPt_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "JetPt_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "JetPt_AfterStandardSelections",
    details          = { "xlabel": "jet p_{T}",
                         "ylabel": "Events/#Deltap_{T}",
                         "divideByBinWidth": True,
                         "unit": "GeV/c",
                         "log": True,
                         "opts": {"ymin": 0.009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "JetEta_AfterStandardSelections",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "JetEta_AfterStandardSelections",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "JetEta_AfterStandardSelections",
    details          = { "xlabel": "jet #eta",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "legendPosition": "SW",
                         "opts": {"ymin": 0.09} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "CollinearTailKillerMinimum",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "ImprovedDeltaPhiCutsCollinearMinimum",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "ImprovedDeltaPhiCutsCollinearMinimum",
    details          = { "xlabel": "R_{coll}^{min}",
        #"xlabel": "min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{o}-#Delta#phi(jet_{1..3},MET))^{2}})",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "^{o}",
                         "log": True,
                         "legendPosition": "SE",
                         "opts": {"ymin": 0.09} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "R_{coll}^{min}", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "BJetSelection",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "NBjets",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "NBjets",
    details          = { "xlabel": "Number of selected b jets",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "opts": {"ymin": 0.09} },
    blindedRange=[],
    #blindedRange     = [1.5,10], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "#geq1 b tag", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "BtagDiscriminator",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "BtagDiscriminator",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "BtagDiscriminator",
    details          = { "xlabel": "b tag discriminator",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "legendPosition": "NE",
                         "opts": {"ymin": 0.9} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "BJetPt",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "BJetPt",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "BJetPt",
    details          = { "xlabel": "b jet p_{T}",
                         "ylabel": "Events/#Deltap_{T}",
                         "divideByBinWidth": True,
                         "unit": "GeV/c",
                         "log": True,
                         "opts": {"ymin": 0.0009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "BJetEta",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "BJetEta",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "BJetEta",
    details          = { "xlabel": "b jet #eta",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "",
                         "log": True,
                         "legendPosition": "SW",
                         "opts": {"ymin": 0.09} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))


ControlPlots.append(ControlPlotInput(
    title            = "MET",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "MET",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "MET",
    details          = { "xlabel": "E_{T}^{miss}",
                         "ylabel": "Events/#DeltaE_{T}^{miss}",
                         "divideByBinWidth": True,
                         "unit": "GeV",
                         "log": True,
                         "opts": {"ymin": 0.0009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "E_{T}^{miss}", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "METPhi",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "METPhi",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "METPhi",
    details          = { "xlabel": "E_{T}^{miss} #phi",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "^{o}",
                         "log": True,
                         "legendPosition": "SW",
                         "opts": {"ymin": 0.09} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

ControlPlots.append(ControlPlotInput(
    title            = "TauPlusMETPt",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "TauPlusMETPt",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "TauPlusMETPt",
    details          = { "xlabel": "p_{T}(#tau + E_{T}^{miss})",
                         "ylabel": "Events/#Deltap_{T}",
                         "divideByBinWidth": True,
                         "unit": "GeV",
                         "log": True,
                         "opts": {"ymin": 0.0009} },
    blindedRange     = [], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
))

#TODO: add as preselection for all ctrl plots in signal analysis MET30 and/or collinear tail killer and/or full tail killer
#TODO: Add to signal analysis ctrl plots tail killer plots

#ControlPlots.append(ControlPlotInput(
    #title            = "DeltaPhi",
    #signalHistoPath  = "",
    #signalHistoName  = "deltaPhi",
    #EWKfakeHistoPath  = "",
    #EWKfakeHistoName  = "EWKFakeTausDeltaPhi",
    #details          = { "bins": 11,
                         #"rangeMin": 0.0,
                         #"rangeMax": 180.0,
                         #"variableBinSizeLowEdges": [0., 10., 20., 30., 40., 60., 80., 100., 120., 140., 160.], # if an empty list is given, then uniform bin width is used
                         #"binLabels": [], # leave empty to disable bin labels
                         #"xlabel": "#Delta#phi(#tau_{h},E_{T}^{miss})",
                         #"ylabel": "Events",
                         #"unit": "^{o}",
                         #"log": True,
                         #"DeltaRatio": 0.5,
                         #"ymin": 0.9,
                         #"ymax": -1},
    #blindedRange     = [-1, 300], # specify range min,max if blinding applies to this control plot
    #evaluationRange  = [], # specify range to be evaluated and saved into a file
    #flowPlotCaption  = "N_{b jets}", # Leave blank if you don't want to include the item to the selection flow plot
#))

#ControlPlots.append(ControlPlotInput(
    #title            = "MaxDeltaPhi",
    #signalHistoPath  = "",
    #signalHistoName  = "maxDeltaPhiJetMet",
    #details          = { "bins": 18,
                         #"rangeMin": 0.0,
                         #"rangeMax": 180.0,
                         #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         #"binLabels": [], # leave empty to disable bin labels
                         #"xlabel": "max(#Delta#phi(jet,E_{T}^{miss})",
                         #"ylabel": "Events",
                         #"unit": "^{o}",
                         #"log": True,
                         #"DeltaRatio": 0.5,
                         #"ymin": 0.9,
                         #"ymax": -1},
    #blindedRange     = [-1, 300], # specify range min,max if blinding applies to this control plot
    #evaluationRange  = [], # specify range to be evaluated and saved into a file
    #flowPlotCaption  = "#Delta#phi(#tau_{h},E_{T}^{miss})", # Leave blank if you don't want to include the item to the selection flow plot
#))

#ControlPlots.append(ControlPlotInput(
    #title            = "WMass",
    #signalHistoPath  = "TopChiSelection",
    #signalHistoName  = "WMass",
    #details          = { "bins": 20,
                         #"rangeMin": 0.0,
                         #"rangeMax": 200.0,
                         #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         #"binLabels": [], # leave empty to disable bin labels
                         #"xlabel": "m_{jj}",
                         #"ylabel": "Events",
                         #"unit": "GeV/c^{2}",
                         #"log": True,
                         #"DeltaRatio": 0.5,
                         #"ymin": 0.9,
                         #"ymax": -1},
    #blindedRange     = [-1, 400], # specify range min,max if blinding applies to this control plot
    #evaluationRange  = [], # specify range to be evaluated and saved into a file
    #flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
#))

#ControlPlots.append(ControlPlotInput(
    #title            = "TopMass",
    #signalHistoPath  = "TopChiSelection",
    #signalHistoName  = "TopMass",
    #details          = { "bins": 20,
                         #"rangeMin": 0.0,
                         #"rangeMax": 400.0,
                         #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                         #"binLabels": [], # leave empty to disable bin labels
                         #"xlabel": "m_{bjj}",
                         #"ylabel": "Events",
                         #"unit": "GeV/c^{2}",
                         #"log": True,
                         #"DeltaRatio": 0.5,
                         #"ymin": 0.9,
                         #"ymax": -1},
    #blindedRange     = [-1, 400], # specify range min,max if blinding applies to this control plot
    #evaluationRange  = [], # specify range to be evaluated and saved into a file
    #flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
#))


ControlPlots.append(ControlPlotInput(
    title            = "BackToBackTailKillerMinimum",
    signalHistoPath  = "ForDataDrivenCtrlPlots",
    signalHistoName  = "ImprovedDeltaPhiCutsBackToBackMinimum",
    EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
    EWKfakeHistoName  = "ImprovedDeltaPhiCutsBackToBackMinimum",
    details          = { "xlabel": "R_{bb}^{min}",
    #"xlabel": "min(#sqrt{(180^{o}-#Delta#phi(#tau,MET))^{2}+#Delta#phi(jet_{1..3},MET)^{2}})",
                         "ylabel": "Events",
                         "divideByBinWidth": False,
                         "unit": "^{o}",
                         "log": True,
                         "legendPosition": "SE",
                         "opts": {"ymin": 0.09} },
    blindedRange     = [81,159], # specify range min,max if blinding applies to this control plot
    evaluationRange  = [], # specify range to be evaluated and saved into a file
    flowPlotCaption  = "R_{bb}^{min}", # Leave blank if you don't want to include the item to the selection flow plot
))

if OptionMassShape == "TransverseMass":
    ControlPlots.append(ControlPlotInput(
        title            = "TransverseMass",
        signalHistoPath  = "",
        signalHistoName  = "shapeTransverseMass",
        EWKfakeHistoPath  = "",
        EWKfakeHistoName  = "shapeEWKFakeTausTransverseMass",
        details          = { "xlabel": "m_{T}(#tau_{h},E_{T}^{miss})",
                         "ylabel": "Events/#Deltam_{T}",
                         "divideByBinWidth": True,
                         "unit": "GeV",
                         "log": False,
                         "opts": {"ymin": 0.0}},
        blindedRange     = [-1, 1000], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [60, 180], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "final", # Leave blank if you don't want to include the item to the selection flow plot
    ))
    ControlPlots.append(ControlPlotInput(
        title            = "TransverseMassLog",
        signalHistoPath  = "",
        signalHistoName  = "shapeTransverseMass",
        EWKfakeHistoPath  = "",
        EWKfakeHistoName  = "shapeEWKFakeTausTransverseMass",
        details          = { "xlabel": "m_{T}(#tau_{h},E_{T}^{miss})",
                         "ylabel": "Events/#Deltam_{T}",
                         "divideByBinWidth": True,
                         "unit": "GeV",
                         "log": True,
                         "opts": {"ymin": 0.009}},
        blindedRange     = [-1, 1000], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))
elif OptionMassShape == "FullMass":
    ControlPlots.append(ControlPlotInput(
        title            = "FullMass",
        signalHistoPath  = "",
        signalHistoName  = "shapeInvariantMass",
        EWKfakeHistoPath  = "",
        EWKfakeHistoName  = "shapeEWKFakeTausInvariantMass",
        details          = { "xlabel": "m(#tau_{h},E_{T}^{miss})",
                             "ylabel": "Events/#Deltam",
                             "divideByBinWidth": True,
                             "unit": "GeV",
                             "log": False,
                             "opts": {"ymin": 0.0},
                             "opts": {"ymin": 1e-5} },
        blindedRange     = [-1, 1000], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [80, 180], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "final", # Leave blank if you don't want to include the item to the selection flow plot
    ))

if OptionCtrlPlotsAtMt:
    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_pT_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_pT_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_pT_AfterMtSelections",
        details          = { "xlabel": "Selected #tau p_{T}",
                             "ylabel": "Events/#Deltap_{T}",
                             "divideByBinWidth": True,
                             "unit": "GeV/c",
                             "log": True,
                             "opts": {"ymin": 0.0009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_p_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_p_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_p_AfterMtSelections",
        details          = { "xlabel": "Selected #tau p",
                             "ylabel": "Events/#Deltap",
                             "divideByBinWidth": True,
                             "unit": "GeV/c",
                             "log": True,
                             "opts": {"ymin": 0.009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_eta_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_eta_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_eta_AfterMtSelections",
        details          = { "xlabel": "Selected #tau #eta",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "legendPosition": "SW",
                             "opts": {"ymin": 0.009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_phi_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_phi_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_phi_AfterMtSelections",
        details          = { "xlabel": "Selected #tau #phi",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "^{o}",
                             "log": True,
                             "legendPosition": "SW",
                             "opts": {"ymin": 0.09} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_LeadingTrackPt_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_LeadingTrackPt_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_LeadingTrackPt_AfterMtSelections",
        details          = { "xlabel": "#tau leading track p_{T}",
                             "ylabel": "Events/#Deltap_{T}",
                             "divideByBinWidth": True,
                             "unit": "GeV/c",
                             "log": True,
                             "ratioLegendPosition": "right",
                             "opts": {"ymin": 0.0009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_LeadingTrackP_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_LeadingTrackP_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_LeadingTrackP_AfterMtSelections",
        details          = { "xlabel": "#tau leading track p",
                             "ylabel": "Events/#Deltap",
                             "divideByBinWidth": True,
                             "unit": "GeV/c",
                             "log": True,
                             "ratioLegendPosition": "right",
                             "opts": {"ymin": 0.0009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_Rtau_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_Rtau_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_Rtau_AfterMtSelections",
        details          = { "xlabel": "Selected #tau R_{#tau}",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "legendPosition": "SE",
                             "opts": {"ymin": 0.009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "SelectedTau_DecayMode_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "SelectedTau_DecayMode_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "SelectedTau_DecayMode_AfterMtSelections",
        details          = { "xlabel": "Selected #tau Decay mode",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "ratioLegendPosition": "right",
                             "opts": {"ymin": 0.9} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "Njets_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "Njets_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "Njets_AfterMtSelections",
        details          = { "xlabel": "Number of selected jets",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "opts": {"ymin": 0.9} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "JetPt_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "JetPt_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "JetPt_AfterMtSelections",
        details          = { "xlabel": "jet p_{T}",
                             "ylabel": "Events/Deltap_{T}",
                             "divideByBinWidth": True,
                             "unit": "GeV/c",
                             "log": True,
                             "opts": {"ymin": 0.009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "JetEta_AfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "JetEta_AfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "JetEta_AfterMtSelections",
        details          = { "xlabel": "jet #eta",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "legendPosition": "SW",
                             "opts": {"ymin": 0.09} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "CollinearTailKillerMinimumAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "ImprovedDeltaPhiCutsCollinearMinimumAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "ImprovedDeltaPhiCutsCollinearMinimumAfterMtSelections",
        details          = { "xlabel": "R_{coll}^{min}",
        #"xlabel": "min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{o}-#Delta#phi(jet_{1..3},MET))^{2}})",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "^{o}",
                             "log": True,
                             "legendPosition": "SE",
                             "opts": {"ymin": 0.09} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "BJetSelectionAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "NBjetsAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "NBjetsAfterMtSelections",
        details          = { "xlabel": "Number of selected b jets",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "opts": {"ymin": 0.09} },
        blindedRange=[],
        #blindedRange     = [1.5,10], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "BtagDiscriminatorAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "BtagDiscriminatorAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "BtagDiscriminatorAfterMtSelections",
        details          = { "xlabel": "b tag discriminator",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "legendPosition": "NE",
                             "opts": {"ymin": 0.9} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "BJetPtAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "BJetPtAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "BJetPtAfterMtSelections",
        details          = { "xlabel": "b jet p_{T}",
                             "ylabel": "Events/#Deltap_{T}",
                             "divideByBinWidth": True,
                             "unit": "GeV/c",
                             "log": True,
                             "opts": {"ymin": 0.0009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "BJetEtaAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "BJetEtaAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "BJetEtaAfterMtSelections",
        details          = { "xlabel": "b jet #eta",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "",
                             "log": True,
                             "legendPosition": "SW",
                             "opts": {"ymin": 0.09} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "METAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "METAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "METAfterMtSelections",
        details          = { "xlabel": "E_{T}^{miss}",
                             "ylabel": "Events/#DeltaE_{T}^{miss}",
                             "divideByBinWidth": True,
                             "unit": "GeV",
                             "log": True,
                             "opts": {"ymin": 0.0009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "METPhiAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "METPhiAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "METPhiAfterMtSelections",
        details          = { "xlabel": "E_{T}^{miss} #phi",
                             "ylabel": "Events/#DeltaE_{T}^{miss}#phi",
                             "divideByBinWidth": True,
                             "unit": "^{o}",
                             "log": True,
                             "legendPosition": "SW",
                             "opts": {"ymin": 0.09} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    ControlPlots.append(ControlPlotInput(
        title            = "TauPlusMETPtAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "TauPlusMETPtAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "TauPlusMETPtAfterMtSelections",
        details          = { "xlabel": "p_{T}(#tau + E_{T}^{miss})",
                             "ylabel": "Events/#Deltap_{T}",
                             "divideByBinWidth": True,
                             "unit": "GeV",
                             "log": True,
                             "opts": {"ymin": 0.0009} },
        blindedRange     = [], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

    #ControlPlots.append(ControlPlotInput(
        #title            = "DeltaPhiAfterMtSelections",
        #signalHistoPath  = "",
        #signalHistoName  = "deltaPhiAfterMtSelections",
        #EWKfakeHistoPath  = "",
        #EWKfakeHistoName  = "EWKFakeTausDeltaPhiAfterMtSelections",
        #details          = { "bins": 11,
                             #"rangeMin": 0.0,
                             #"rangeMax": 180.0,
                             #"variableBinSizeLowEdges": [0., 10., 20., 30., 40., 60., 80., 100., 120., 140., 160.], # if an empty list is given, then uniform bin width is used
                             #"binLabels": [], # leave empty to disable bin labels
                             #"xlabel": "#Delta#phi(#tau_{h},E_{T}^{miss})",
                             #"ylabel": "Events",
                             #"unit": "^{o}",
                             #"log": True,
                             #"DeltaRatio": 0.5,
                             #"ymin": 0.9,
                             #"ymax": -1},
        #blindedRange     = [-1, 300], # specify range min,max if blinding applies to this control plot
        #evaluationRange  = [], # specify range to be evaluated and saved into a file
        #flowPlotCaption  = "N_{b jets}", # Leave blank if you don't want to include the item to the selection flow plot
    #))

    #ControlPlots.append(ControlPlotInput(
        #title            = "MaxDeltaPhi",
        #signalHistoPath  = "",
        #signalHistoName  = "maxDeltaPhiJetMet",
        #details          = { "bins": 18,
                             #"rangeMin": 0.0,
                             #"rangeMax": 180.0,
                             #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                             #"binLabels": [], # leave empty to disable bin labels
                             #"xlabel": "max(#Delta#phi(jet,E_{T}^{miss})",
                             #"ylabel": "Events",
                             #"unit": "^{o}",
                             #"log": True,
                             #"DeltaRatio": 0.5,
                             #"ymin": 0.9,
                             #"ymax": -1},
        #blindedRange     = [-1, 300], # specify range min,max if blinding applies to this control plot
        #evaluationRange  = [], # specify range to be evaluated and saved into a file
        #flowPlotCaption  = "#Delta#phi(#tau_{h},E_{T}^{miss})", # Leave blank if you don't want to include the item to the selection flow plot
    #))

    #ControlPlots.append(ControlPlotInput(
        #title            = "WMass",
        #signalHistoPath  = "TopChiSelection",
        #signalHistoName  = "WMass",
        #details          = { "bins": 20,
                             #"rangeMin": 0.0,
                             #"rangeMax": 200.0,
                             #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                             #"binLabels": [], # leave empty to disable bin labels
                             #"xlabel": "m_{jj}",
                             #"ylabel": "Events",
                             #"unit": "GeV/c^{2}",
                             #"log": True,
                             #"DeltaRatio": 0.5,
                             #"ymin": 0.9,
                             #"ymax": -1},
        #blindedRange     = [-1, 400], # specify range min,max if blinding applies to this control plot
        #evaluationRange  = [], # specify range to be evaluated and saved into a file
        #flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    #))

    #ControlPlots.append(ControlPlotInput(
        #title            = "TopMass",
        #signalHistoPath  = "TopChiSelection",
        #signalHistoName  = "TopMass",
        #details          = { "bins": 20,
                             #"rangeMin": 0.0,
                             #"rangeMax": 400.0,
                             #"variableBinSizeLowEdges": [], # if an empty list is given, then uniform bin width is used
                             #"binLabels": [], # leave empty to disable bin labels
                             #"xlabel": "m_{bjj}",
                             #"ylabel": "Events",
                             #"unit": "GeV/c^{2}",
                             #"log": True,
                             #"DeltaRatio": 0.5,
                             #"ymin": 0.9,
                             #"ymax": -1},
        #blindedRange     = [-1, 400], # specify range min,max if blinding applies to this control plot
        #evaluationRange  = [], # specify range to be evaluated and saved into a file
        #flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    #))

    ControlPlots.append(ControlPlotInput(
        title            = "BackToBackTailKillerMinimumAfterMtSelections",
        signalHistoPath  = "ForDataDrivenCtrlPlots",
        signalHistoName  = "ImprovedDeltaPhiCutsBackToBackMinimumAfterMtSelections",
        EWKfakeHistoPath  = "ForDataDrivenCtrlPlotsEWKFakeTaus",
        EWKfakeHistoName  = "ImprovedDeltaPhiCutsBackToBackMinimumAfterMtSelections",
        details          = { #"xlabel": "min(#sqrt{(180^{o}-#Delta#phi(#tau,MET))^{2}+#Delta#phi(jet_{1..3},MET)^{2}})",
                             "xlabel": "R_{bb}^{min}",
                             "ylabel": "Events",
                             "divideByBinWidth": False,
                             "unit": "^{o}",
                             "log": True,
                             "legendPosition": "SE",
                             "opts": {"ymin": 0.09} },
        blindedRange     = [81,159], # specify range min,max if blinding applies to this control plot
        evaluationRange  = [], # specify range to be evaluated and saved into a file
        flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    ))

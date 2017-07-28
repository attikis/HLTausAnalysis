#!/usr/bin/env python
###############################################################
### All imported modules
###############################################################
import ROOT
import os
import sys
import numpy
import math
import tools.plotter as m_plotter
import tools.histos as m_histos
from tools.l1tktauAux import *
import tools.styles as m_styles
import tools.aux as m_aux

###############################################################
### Options here
###############################################################
bDoL1TkTau       = True
bDoL1TkTauExtra  = False
bDoMatchTk       = False
bDoSigTks        = True
bDoIsoTks        = False
bDoEfficiencies  = True
bTTPixelTracks   = True
bTTTracks        = not bTTPixelTracks
nFitParams       = "5FitParams/"
datasetList      = ["MinBias", "VBF"] #["MinBias", "VBF", "HPlus160"]
saveFormats      = ["png", "pdf"]
if bTTPixelTracks:
    tkType = "pixTks"
else:
    tkType = "ttTks"

inputPath    = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/"
#inputPath    = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/" + nFitParams + tkType + "/default/"
#inputPath    = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/"
#inputPath   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/trial5/"
#inputPath   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/trial6/"
#inputPath   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/finalTrial/Tk/"
#inputPath   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/"
#inputPath   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/tauTrigger/5FitParams/" + tkType + "/optimisation/new/"
savePath     = ""
#savePath     = "/Users/attikis/talks/TauTrigger_13March2015/figures/5FitParams/" + tkType + "/"
#savePath    = "/Users/attikis/talks/Phase2_Trigger_Workshop_17March2015/figures/5FitParams/" + tkType + "/"

datasetPaths   = {}
datasetPaths["MinBias"]           = inputPath + "TauTrigger_Histograms_MinBias.root"
datasetPaths["PiPlus"]            = inputPath + "TauTrigger_Histograms_PiPlus.root"
datasetPaths["PiMinus"]           = inputPath + "TauTrigger_Histograms_PiMinu.root"
datasetPaths["SingleTauGun1p"]    = inputPath + "TauTrigger_Histograms_SingleTauGun1p.root"
datasetPaths["DiTauGun3p"]        = inputPath + "TauTrigger_Histograms_DiTauGun3p.root"
datasetPaths["VBF"]               = inputPath + "TauTrigger_Histograms_VBF.root"
datasetPaths["TTbar"]             = inputPath + "TauTrigger_Histograms_TTBar.root"
datasetPaths["HPlus200"]          = inputPath + "TauTrigger_Histograms_HPlus200.root"
datasetPaths["HPlus160"]          = inputPath + "TauTrigger_Histograms_HPlus160.root"

###############################################################
### Main
###############################################################
def PrintTreeCreationParameters():
    
    inputPath = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v3/TTracks_TTPixelTracks_CaloTausCorr_FitParam5_PtMin2_TTI2023Upg14D_TestVersionWithPixelBugs_28Jan2015_144735/"
    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    p.GetDatasets("TTI2023Upg14D", inputPath, ["minbias"])
    p.PrintPSet("TkTauFromCaloNTupleMaker/configInfo/parameterSet")
    
    return

###############################################################
def DoPlots(hList, datasetList, saveExt=""):

    p = m_plotter.Plotter( Verbose=False, BatchMode=True )
    for dataset in datasetList:
        p.AddDataset(dataset, datasetPaths[dataset])
    p.SetBoolUseDatasetAsLegEntry(True)
    p.AddHisto(hList)
    p.EnableColourPalette(False)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False, bAddReferenceHisto = False)
    #p.GetTLegend().SetHeader( "" )
    p.SaveHistos(True, savePath, saveFormats, saveExt)

    return

###############################################################
def DoEfficiency(hList, datasetList, cutDirList, saveExt=""):

    for cutDir in cutDirList:
        
        p = m_plotter.Plotter( Verbose=False, BatchMode=True )
        for dataset in datasetList:
            p.AddDataset(dataset, datasetPaths[dataset])
        p.SetBoolUseDatasetAsLegEntry(True)
        p.AddHisto(hList)
        p.EnableColourPalette(False)
        p.DrawEfficiency(cutDir, "binomial")
        p.GetTLegend().SetHeader( "L1TkTau" )
        p.SaveHistos(True, savePath, saveFormats, saveExt)
        
    return

    ###############################################################
if __name__ == "__main__":

    #PrintTreeCreationParameters()
    if bDoL1TkTau:
        # DoPlots( hL1TkTau_Multiplicity, datasetList) 
        DoPlots( hL1TkTau_CHF        , datasetList) 
        DoPlots( hL1TkTau_NHF        , datasetList) 
        DoPlots( hL1TkTau_NHFAbs     , datasetList) 
        DoPlots( hL1TkTau_InvMass    , datasetList) 
        DoPlots( hL1TkTau_Rtau       , datasetList) 
        # DoPlots( hL1TkTau_RelIso     , datasetList) 
        # DoPlots( hL1TkTau_VtxIso     , datasetList) 
        # DoPlots( hL1TkTau_VtxIsoAbs  , datasetList) 
        # DoPlots( hL1TkTau_InvMassIncl, datasetList) 


    if bDoL1TkTauExtra:
        DoPlots( hL1TkTau_SigConeRMin, datasetList) 
        DoPlots( hL1TkTau_SigConeRMax, datasetList) 
        DoPlots( hL1TkTau_IsoConeRMin, datasetList) 
        DoPlots( hL1TkTau_IsoConeRMax, datasetList) 
        DoPlots( hL1TkTau_DeltaRGenP , datasetList) 
        DoPlots( hL1TkTau_Charge     , datasetList) 


    if bDoEfficiencies:
        DoEfficiency( hL1TkTau_CHF                  , datasetList, ["<"] )
        DoEfficiency( hL1TkTau_NHF                  , datasetList, ["<", ">"] )
        DoEfficiency( hL1TkTau_NHFAbs               , datasetList, ["<", ">"] )
        # DoEfficiency( hL1TkTau_InvMass              , datasetList, ["<"] )
        # DoEfficiency( hL1TkTau_Rtau                 , datasetList, [">", "<"] )
        # DoEfficiency( hL1TkTau_RelIso               , datasetList, ["<"] )
        # DoEfficiency( hL1TkTau_VtxIsoAbs            , datasetList, [">"] )


    if bDoMatchTk:
        if bTTPixelTracks:
            # DoPlots( hL1TkTau_MatchPixTk_Pt           , datasetList)
            # DoPlots( hL1TkTau_MatchPixTk_PtRel        , datasetList)
            # DoPlots( hL1TkTau_MatchPixTk_ChiSquared   , datasetList)
            # DoPlots( hL1TkTau_MatchPixTk_RedChiSquared, datasetList)
            # DoPlots( hL1TkTau_MatchPixTk_d0Abs        , datasetList)
            # DoPlots( hL1TkTau_MatchPixTk_d0SigAbs     , datasetList)
            # DoPlots( hL1TkTau_MatchPixTk_DeltaR       , datasetList)
            # DoPlots( hL1TkTau_MatchPixTk_NStubs       , datasetList)
            DoPlots( hL1TkTau_MatchPixTk_StubPtCons   , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_Eta          , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_NPixHits     , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_NPsStubs     , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_NBarrelStubs , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_NEndcapStubs , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_d0           , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_d0Sig        , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_POCAz        , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_POCAzSig     , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_SigmaRInv    , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_SigmaPhi0    , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_SigmaD0      , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_SigmaT       , datasetList)
            ### DoPlots( hL1TkTau_MatchPixTk_SigmaZ0      , datasetList)
        if bTTTracks:
            DoPlots( hL1TkTau_MatchTk_PtRel         , datasetList)
            DoPlots( hL1TkTau_MatchTk_DeltaR        , datasetList)
            DoPlots( hL1TkTau_MatchTk_Pt            , datasetList)
            DoPlots( hL1TkTau_MatchTk_d0Abs         , datasetList)
            DoPlots( hL1TkTau_MatchTk_ChiSquared    , datasetList)
            DoPlots( hL1TkTau_MatchTk_RedChiSquared , datasetList)
            DoPlots( hL1TkTau_MatchTk_NStubs        , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_NPsStubs      , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_NBarrelStubs  , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_NEndcapStubs  , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_Eta           , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_POCAz         , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_d0            , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_IsGenuine     , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_IsUnknown     , datasetList)
            ### DoPlots( hL1TkTau_MatchTk_IsCombinatoric, datasetList)


    if bDoSigTks:
        DoPlots( hL1TkTau_NSigTks          , datasetList)
        # DoPlots( hL1TkTau_SigTks_Pt        , datasetList)
        # DoPlots( hL1TkTau_SigTks_Eta       , datasetList)
        DoPlots( hL1TkTau_SigTks_d0        , datasetList)
        DoPlots( hL1TkTau_SigTks_d0Abs     , datasetList)
        DoPlots( hL1TkTau_SigTks_d0Sig     , datasetList)
        DoPlots( hL1TkTau_SigTks_d0SigAbs  , datasetList)
        DoPlots( hL1TkTau_SigTks_PtRel     , datasetList)
        # DoPlots( hL1TkTau_SigTks_POCAz     , datasetList)
        DoPlots( hL1TkTau_SigTks_DeltaPOCAz, datasetList)
        DoPlots( hL1TkTau_SigTks_StubPtCons, datasetList)


    if bDoIsoTks:
        DoPlots( hL1TkTau_NIsoTks          , datasetList)
        # DoPlots( hL1TkTau_IsoTks_Pt        , datasetList)
        # DoPlots( hL1TkTau_IsoTks_Eta       , datasetList)
        DoPlots( hL1TkTau_IsoTks_d0        , datasetList)
        DoPlots( hL1TkTau_IsoTks_d0Abs     , datasetList)
        DoPlots( hL1TkTau_IsoTks_d0Sig     , datasetList)
        DoPlots( hL1TkTau_IsoTks_d0SigAbs  , datasetList)
        DoPlots( hL1TkTau_IsoTks_PtRel     , datasetList)
        # DoPlots( hL1TkTau_IsoTks_POCAz     , datasetList)
        DoPlots( hL1TkTau_IsoTks_DeltaPOCAz, datasetList)
        DoPlots( hL1TkTau_IsoTks_StubPtCons, datasetList)



###############################################################     

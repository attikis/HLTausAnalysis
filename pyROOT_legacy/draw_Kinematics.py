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
from tools.kinematicsAux import *
import tools.styles as m_styles
import tools.aux as m_aux

###############################################################
### Options here
###############################################################
bBatchMode  = True
bVerbose    = False
bSavePlots  = True
bSingleTau  = True
saveFormats = ["png", "pdf"]
datasetList = ["VBF"] #SingleTauOneProng , TauThreeProngsEnr, VBF , TTbar , HPlus160 , HPlus200
visEt       = "20" #"20"
#inputPath   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/"
inputPath   = "/Users/attikis/hltaus/pyROOT/Macros/TauTrigger/results/kinematics/visEt" + visEt + "/"
savePath    = "/Users/attikis/talks/TauTrigger_13March2015/figures/kinematics/" + datasetList[0] + "/"
#savePath    = ""

datasetPaths_all = {}
tauDecayMode     = "_3p" #"_all" #"_1p" #"_3p"
datasetPaths_all["MinBias"]           = inputPath + "TauTriggerKinematics_Histograms_nugun%s.root"       % (tauDecayMode)
datasetPaths_all["VBF"]               = inputPath + "TauTriggerKinematics_Histograms_vbf%s.root"         % (tauDecayMode)
datasetPaths_all["SingleTauOneProng"] = inputPath + "TauTriggerKinematics_Histograms_onetaugun1p%s.root" % (tauDecayMode)
datasetPaths_all["TauThreeProngsEnr"] = inputPath + "TauTriggerKinematics_Histograms_twotaugun3p%s.root" % (tauDecayMode)
datasetPaths_all["TTbar"]             = inputPath + "TauTriggerKinematics_Histograms_ttbar%s.root"       % (tauDecayMode)
datasetPaths_all["HPlus160"]          = inputPath + "TauTriggerKinematics_Histograms_litehiggs%s.root"   % (tauDecayMode)
datasetPaths_all["HPlus200"]          = inputPath + "TauTriggerKinematics_Histograms_heavyhiggs%s.root"  % (tauDecayMode)

datasetPaths_1p = {}
tauDecayMode    = "_1p"
datasetPaths_1p["MinBias"]           = inputPath + "TauTriggerKinematics_Histograms_nugun%s.root"       % (tauDecayMode)
datasetPaths_1p["VBF"]               = inputPath + "TauTriggerKinematics_Histograms_vbf%s.root"         % (tauDecayMode)
datasetPaths_1p["SingleTauOneProng"] = inputPath + "TauTriggerKinematics_Histograms_onetaugun1p%s.root" % (tauDecayMode)
datasetPaths_1p["TauThreeProngsEnr"] = inputPath + "TauTriggerKinematics_Histograms_twotaugun3p%s.root" % (tauDecayMode)
datasetPaths_1p["TTbar"]             = inputPath + "TauTriggerKinematics_Histograms_ttbar%s.root"       % (tauDecayMode)
datasetPaths_1p["HPlus160"]          = inputPath + "TauTriggerKinematics_Histograms_litehiggs%s.root"   % (tauDecayMode)
datasetPaths_1p["HPlus200"]          = inputPath + "TauTriggerKinematics_Histograms_heavyhiggs%s.root"  % (tauDecayMode)

datasetPaths_3p = {}
tauDecayMode    = "_3p"
datasetPaths_3p["MinBias"]           = inputPath + "TauTriggerKinematics_Histograms_nugun%s.root"       % (tauDecayMode)
datasetPaths_3p["VBF"]               = inputPath + "TauTriggerKinematics_Histograms_vbf%s.root"         % (tauDecayMode)
datasetPaths_3p["SingleTauOneProng"] = inputPath + "TauTriggerKinematics_Histograms_onetaugun1p%s.root" % (tauDecayMode)
datasetPaths_3p["TauThreeProngsEnr"] = inputPath + "TauTriggerKinematics_Histograms_twotaugun3p%s.root" % (tauDecayMode)
datasetPaths_3p["TTbar"]             = inputPath + "TauTriggerKinematics_Histograms_ttbar%s.root"       % (tauDecayMode)
datasetPaths_3p["HPlus160"]          = inputPath + "TauTriggerKinematics_Histograms_litehiggs%s.root"   % (tauDecayMode)
datasetPaths_3p["HPlus200"]          = inputPath + "TauTriggerKinematics_Histograms_heavyhiggs%s.root"  % (tauDecayMode)

###############################################################
### Main
###############################################################
def PrintTreeCreationParameters():
    
    inputPath = "/Users/attikis/hltaus/rootFiles/TTrees/CMSSW_6_2_0_SLHC12_patch1/TkTauFromCaloAnalyzer_v2/SigTk5_R01to04_TTI2023Upg14D_DefaultTrackingSequence_08Aug2014_164642/"
    p = m_plotter.Plotter( bVerbose, bBatchMode )
    p.GetDatasets("TTI2023Upg14D", inputPath, ["minbias"])
    p.PrintPSet("TkTauFromCaloAnalyzer/configInfo/parameterSet")

    return

###############################################################
def DoPlots(hList, datasetList, tauDecayMode=""):

    saveExt = ""
    p = m_plotter.Plotter( bVerbose, bBatchMode )
    for dataset in datasetList:
        if tauDecayMode=="":
            p.AddDataset(dataset + ":all" , datasetPaths_all[dataset])
            p.AddDataset(dataset + ":1pr", datasetPaths_1p[dataset])
            p.AddDataset(dataset + ":3pr", datasetPaths_3p[dataset])
        else:
            saveExt = "_" + tauDecayMode
            p.AddDataset(dataset + ":" + tauDecayMode, datasetPaths_all[dataset])
    p.SetBoolUseDatasetAsLegEntry(True)
    p.AddHisto(hList)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p.SetTLegendHeader( dataset, tauDecayMode )
    p.SaveHistos(bSavePlots, savePath, saveFormats, saveExt)

    return


###############################################################
def DoPlotsWithTF1(hList, datasetList, myF1Exp, xMin, xMax, tauDecayMode="", kwargs=None):

    saveExt = ""
    p = m_plotter.Plotter( bVerbose, bBatchMode )
    for dataset in datasetList:
        if tauDecayMode=="":
            p.AddDataset(dataset + ":all" , datasetPaths_all[dataset])
            p.AddDataset(dataset + ":1pr", datasetPaths_1p[dataset])
            p.AddDataset(dataset + ":3pr", datasetPaths_3p[dataset])
        else:
            saveExt = "_" + tauDecayMode
            p.AddDataset(dataset + ":" + tauDecayMode, datasetPaths_all[dataset])
    p.SetBoolUseDatasetAsLegEntry(True)
    p.AddHisto(hList)
    p.AddTF1( myF1Exp, xMin, xMax, kwargs)
    p.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
    p.SetTLegendHeader( dataset, tauDecayMode )
    p.SaveHistos(bSavePlots, savePath, saveFormats,  saveExt)
    
    return

###############################################################
def DoEfficiency(hList, datasetList, cutDirList, tauDecayMode=""):

    saveExt = ""
    for histo in hList:
        for cutDir in cutDirList:
                
            p = m_plotter.Plotter( bVerbose, bBatchMode )
            for dataset in datasetList:        
                if tauDecayMode=="":
                    p.AddDataset(dataset + ":all" , datasetPaths_all[dataset])
                    p.AddDataset(dataset + ":1pr", datasetPaths_1p[dataset])
                    p.AddDataset(dataset + ":3pr", datasetPaths_3p[dataset])
                else:
                    saveExt = "_" + tauDecayMode
                    p.AddDataset(dataset + ":" + tauDecayMode, datasetPaths_all[dataset])
            p.AddHisto([histo])
            p.SetBoolUseDatasetAsLegEntry(True)
            p.EnableColourPalette(True)
            p.DrawEfficiency(cutDir, "binomial")
            p.SetTLegendHeader( dataset, tauDecayMode )
            p.SaveHistos(bSavePlots, savePath, saveFormats, saveExt)
        
    return

###############################################################
def DoProfile(hList, datasetList, tauDecayMode, profileAxis="x"):

    saveExt = ""
    p = m_plotter.Plotter( bVerbose, bBatchMode )
    for dataset in datasetList:
        if tauDecayMode=="":
            p.AddDataset(dataset + ":all" , datasetPaths_all[dataset])
            p.AddDataset(dataset + ":1pr", datasetPaths_1p[dataset])
            p.AddDataset(dataset + ":3pr", datasetPaths_3p[dataset])
        else:
            saveExt = "_" + tauDecayMode
            p.AddDataset(dataset + ":" + tauDecayMode, datasetPaths_all[dataset])
    p.SetBoolUseDatasetAsLegEntry(True)
    p.AddHisto(hList)
    p.Draw2DProfile(THStackDrawOpt="nostack", bStackInclusive=False, ProfileAxis=profileAxis, firstBin=1, lastBin=-1)
    p.SetTLegendHeader( dataset, "")
    p.SaveHistos(bSavePlots, savePath, saveFormats, saveExt)

    return

###############################################################
if __name__ == "__main__":
    
    '''
    #PrintTreeCreationParameters()
    ### Tau
    #DoPlots( [HadronicTau_N]      , datasetList)
    DoPlots( [HadronicTau_Pt]     , datasetList)
    #DoPlots( [HadronicTau_Eta]    , datasetList)
    #DoPlots( [HadronicTau_Phi]    , datasetList)
    #DoPlots( [HadronicTau_Mass]   , datasetList)
    #DoPlots( [HadronicTau_Charge] , datasetList)
    #DoPlots( [HadronicTau_PdgId]  , datasetList)
    #DoPlots( [HadronicTau_Status] , datasetList)
    #DoPlots( [HadronicTau_VertexX], datasetList)
    #DoPlots( [HadronicTau_VertexY], datasetList)
    #DoPlots( [HadronicTau_VertexZ], datasetList)

    ### Visible-Tau
    DoPlots( [HadronicTau_VisEt]       , datasetList)
    #DoPlots( [HadronicTau_VisEta]      , datasetList)
    #DoPlots( [HadronicTau_VisPhi]      , datasetList)
    DoPlots( [HadronicTau_VisMass]     , datasetList)
    #DoPlots( [HadronicTau_DecayMode]   , datasetList)
    DoPlots( [HadronicTau_EcalFraction], datasetList)
    DoPlots( [HadronicTau_HcalFraction], datasetList)

    ### Tau Charged Pions 
    #DoPlots( [HadronicTau_ChargedPion_N]      , datasetList)
    DoPlots( [HadronicTau_ChargedPion_Pt]     , datasetList)
    #DoPlots( [HadronicTau_ChargedPion_Eta]    , datasetList)
    #DoPlots( [HadronicTau_ChargedPion_Phi]    , datasetList)
    #DoPlots( [HadronicTau_ChargedPion_Mass]   , datasetList)
    #DoPlots( [HadronicTau_ChargedPion_Charge] , datasetList)
    #DoPlots( [HadronicTau_ChargedPion_PdgId]  , datasetList)
    #DoPlots( [HadronicTau_ChargedPion_Status] , datasetList)
    #DoPlots( [HadronicTau_ChargedPion_VertexX], datasetList)
    #DoPlots( [HadronicTau_ChargedPion_VertexY], datasetList)
    #DoPlots( [HadronicTau_ChargedPion_VertexZ], datasetList)

    ### Tau Neutral Pions
    #DoPlots( [HadronicTau_NeutralPion_N]      , datasetList)
    DoPlots( [HadronicTau_NeutralPion_Pt]     , datasetList)
    #DoPlots( [HadronicTau_NeutralPion_Eta]    , datasetList)
    #DoPlots( [HadronicTau_NeutralPion_Phi]    , datasetList)
    #DoPlots( [HadronicTau_NeutralPion_Mass]   , datasetList)
    #DoPlots( [HadronicTau_NeutralPion_Charge] , datasetList)
    #DoPlots( [HadronicTau_NeutralPion_PdgId]  , datasetList)
    #DoPlots( [HadronicTau_NeutralPion_Status] , datasetList)
    #DoPlots( [HadronicTau_NeutralPion_VertexX], datasetList)
    #DoPlots( [HadronicTau_NeutralPion_VertexY], datasetList)
    #DoPlots( [HadronicTau_NeutralPion_VertexZ], datasetList)

    ### Leading Charged Pion
    DoPlots(   [HadronicTau_LdgChPion_Pt]     , datasetList)
    DoPlots(   [HadronicTau_LdgChPion_Eta]    , datasetList)
    #DoPlots(   [HadronicTau_LdgChPion_Phi]    , datasetList)
    #DoPlots(   [HadronicTau_LdgChPion_Mass]   , datasetList)
    #DoPlots(   [HadronicTau_LdgChPion_Charge] , datasetList)
    #DoPlots(   [HadronicTau_LdgChPion_PdgId]  , datasetList)
    #DoPlots(   [HadronicTau_LdgChPion_Status] , datasetList)
    #DoPlots(   [HadronicTau_LdgChPion_VertexX], datasetList)
    #DoPlots(   [HadronicTau_LdgChPion_VertexY], datasetList)
    DoPlots(   [HadronicTau_LdgChPion_VertexZ], datasetList)
    DoPlots(   [HadronicTau_LdgChPion_Rtau]   , datasetList)
    DoPlots(   [HadronicTau_LdgChPion_Pt_VisEt], datasetList, "all")
    DoProfile( [HadronicTau_LdgChPion_Pt_VisEt], datasetList, "all", "x")
    
    ### DiTau System
    DoPlots( [DiTau_Pt]     , datasetList)
    #DoPlots( [DiTau_InvMass], datasetList)
    DoPlots( [DiTau_VisMass], datasetList)


    ### Efficiencies
    DoEfficiency( [HadronicTau_Pt]                       , datasetList, [">"] )
    DoEfficiency( [HadronicTau_VisEt]                    , datasetList, [">"] )
    DoEfficiency( [HadronicTau_EcalFraction]             , datasetList, ["<"])
    DoEfficiency( [HadronicTau_HcalFraction]             , datasetList, [">"] )
    #DoEfficiency( [HadronicTau_LdgChPion_Pt]             , datasetList, [">"])
    #DoEfficiency( [HadronicTau_LdgChPion_DeltaRMax]      , datasetList, ["<"], "3pr")
    #DoEfficiency( [HadronicTau_LdgChPion_DeltaRMax_MinPt], datasetList, ["<"], "3pr")
    #DoEfficiency( [HadronicTau_LdgChPion_Rtau]           , datasetList, [">"] )

    '''
    ### Signal Cone Opening
    #DoPlots( [HadronicTau_LdgChPion_DeltaRMax]      , datasetList, "3pr")
    #DoPlots( [HadronicTau_LdgChPion_DeltaRMax_MinPt], datasetList, "3pr")

#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_Pt]      , datasetList, "2.2/x", 0.01, 0.15, "3pr", f1Args)
#    DoProfile(      [HadronicTau_LdgChPion_DeltaRMax_Pt]      , datasetList, "3pr", "x")
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_Pt_MinPt], datasetList, "2.2/x", 0.01, 0.15, "3pr", f1Args)
#    DoProfile(      [HadronicTau_LdgChPion_DeltaRMax_Pt_MinPt], datasetList, "3pr", "x")
    
#    DoProfile(      [HadronicTau_LdgChPion_DeltaRMax_TauEt]      , datasetList, "Profile", "x")
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_TauEt]      , datasetList, "3.5/x", 0.01, 0.15, "3p5GeV", f1Args)
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_TauEt]      , datasetList, "3.0/x", 0.01, 0.15, "3p0GeV", f1Args)
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_TauEt]      , datasetList, "2.5/x", 0.01, 0.15, "2p5GeV", f1Args)

#    DoProfile(      [HadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt]      , datasetList, "Profile", "x")
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt]      , datasetList, "3.5/x", 0.01, 0.15, "3p5GeV", f1Args)
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt]      , datasetList, "3.0/x", 0.01, 0.15, "3p0GeV", f1Args)
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_TauEt_MinPt]      , datasetList, "2.5/x", 0.01, 0.15, "2p5GeV", f1Args)

    DoProfile(      [HadronicTau_LdgChPion_DeltaRMax_VisTauEt]      , datasetList, "Profile", "x")
    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt]      , datasetList, "3.5/x", 0.01, 0.15, "3p5GeV", f1Args)
    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt]      , datasetList, "3.0/x", 0.01, 0.15, "3p0GeV", f1Args)
    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt]      , datasetList, "2.5/x", 0.01, 0.15, "2p5GeV", f1Args)
    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt]      , datasetList, "1.5/x", 0.001, 0.1, "1p5GeV", f1Args)
    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt]      , datasetList, "1.0/x", 0.001, 0.1, "1p0GeV", f1Args)
    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt]      , datasetList, "0.5/x", 0.001, 0.1, "0p5GeV", f1Args)

#    DoProfile(      [HadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt], datasetList, "Profile", "x")
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt], datasetList, "3.5/x", 0.01, 0.15, "3p5GeV", f1Args)
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt], datasetList, "3.0/x", 0.01, 0.15, "3p0GeV", f1Args)
#    DoPlotsWithTF1( [HadronicTau_LdgChPion_DeltaRMax_VisTauEt_MinPt], datasetList, "2.5/x", 0.01, 0.15, "2p5GeV", f1Args)
    

###############################################################

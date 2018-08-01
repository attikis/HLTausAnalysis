#!/usr/bin/env python

# a script for plotting any histogram in the analysis results
# input:
#     - result dir containing at least signalAnalysis and QCDAnalysis 
#       pseudomulticrabs (same input as for the datacard generator)
#     - json file containing descriptions for the plotting
# 02122016/S.Lehti

import os
import sys
import re
import json

import HiggsAnalysis.LimitCalc.MulticrabPathFinder as PathFinder
import HiggsAnalysis.NtupleAnalysis.tools.dataset as dataset
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms

import ROOT
ROOT.gROOT.SetBatch(True)

formats = [".pdf",".png"]
plotDir = "Plots"

def usage():
    print
    print "### Usage:  ",sys.argv[0],"<resultdir containing at least signalAnalysis and QCD pseudomulticrab> <histogram-json>"
    print
    sys.exit()


def plot(resultdir,jsonfile):
    with open(os.path.abspath(jsonfile)) as jfile:
        j = json.load(jfile)
        print "Plotting",j["title"],"in",resultdir

        if "outputdir" in j:
            global plotDir
            plotDir = j["outputdir"]
        multicrabPaths = PathFinder.MulticrabPathFinder(resultdir)
        print "multicrabPaths = ", multicrabPaths.

        paths = []
        if os.path.exists(multicrabPaths.getSignalPath()):
            paths.append(multicrabPaths.getSignalPath())
        if os.path.exists(multicrabPaths.getQCDInvertedPath()):
            paths.append(multicrabPaths.getQCDInvertedPath())
        if os.path.exists(multicrabPaths.getEWKPath()):
            paths.append(multicrabPaths.getEWKPath())

        datasets = dataset.getDatasetsFromMulticrabDirs(paths)

        #datasets.loadLuminosities()
        style = tdrstyle.TDRStyle()
        plots.mergeRenameReorderForDataMC(datasets)

        alldsets = datasets.getAllDatasets()
        print "Merged datasets"
        for d in alldsets:
            print "       ",d.getName()

        # lumi = 0.0
        # for d in datasets.getDataDatasets():
        #     print "luminosity",d.getName(),d.getLuminosity()
        #     lumi += d.getLuminosity()
        # print "luminosity, sum",lumi

        if len(j["samples"])>0:
           for s in j["samples"]:
               #print s
               h = datasets.getDataset(s).getDatasetRootHisto(j["histogram"]).getHistogram()
               name = j["histogram"]+s
               plotgraph([h],lumi,j,name)

def plotgraph(histolist,lumi,j,name=""):

    name = os.path.basename(name)

    p = plots.DataMCPlot2(histolist)
    opts = {}
    if "opts" in j:
        opts = j["opts"]
    opts2 = {}
    if "opts2" in j:
        opts2 = j["opts2"]
    p.createFrame(os.path.join(plotDir, name), opts=opts, opts2=opts2)
    if "xlabel" in j:
        p.getFrame().GetXaxis().SetTitle(j["xlabel"])
    if "ylabel" in j:
        p.getFrame().GetYaxis().SetTitle(j["ylabel"])
#    p.setDrawOptions("COLZ")
    if "drawStyle" in j:
        p.histoMgr.setHistoDrawStyleAll(j["drawStyle"])
    if "rebinx" in j:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().RebinX(j["rebinx"]))
    if "rebiny" in j:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().RebinY(j["rebiny"]))

    p.draw()

    histograms.addStandardTexts(lumi=lumi)

    if not os.path.exists(plotDir):
        os.mkdir(plotDir)
    p.save(formats)
    print "Saved plot",os.path.join(plotDir, name)




"""
        signalname_re = re.compile("Hplus\S+_(M_|M)(?P<mass>\d+)$")
        histolist = []  
#        histolist.append(dataset_signal.getDataset("Data").getDatasetRootHisto(j["histogram"]).getHistogram())
        if not dataset_genuinetau == None:
            histolist.append(dataset_genuinetau.getDataset("Embedding").getDatasetRootHisto(j["histogram"]).getHistogram())
        else:
            for d in dataset_signal.getMCDatasets():
                match = signalname_re.search(d.getName())
                if not match:
                    histolist.append(d.getDatasetRootHisto(j["histogram"]).getHistogram())

        alldsets = dataset_faketau.getAllDatasets()
        for d in alldsets:
            print "Dataset",d.getName()
        h = dataset_faketau.getDataset("QCDMeasurementMT").getDatasetRootHisto(j["histogram"]).getHistogram()
        h.GetEntries()
        sys.exit()
        histolist.append(dataset_faketau.getDataset("QCDMeasurementMT").getDatasetRootHisto(j["histogram"]).getHistogram())





        datasets.extend(dataset_faketau)
        if not dataset_genuinetau == None:
            datasets.extend(dataset_genuinetau)


        plots.mergeRenameReorderForDataMC(datasets)

        hplusMass = j["hplusMass"]
        print hplusMass

        names = datasets.getMCDatasetNames()
        signalname_re = re.compile("Hplus\S+_(M_|M)(?P<mass>\d+)$")
        for n in names:
            match = signalname_re.search(n)
            if match:
                mass = match.group("mass")
                if mass == hplusMass:
                    continue
                datasets.remove(n)

        alldsets = datasets.getAllDatasets()
        for d in alldsets:
            print "Dataset",d.getName()


        
#        h_data     = datasets.getDataset("Data").getDatasetRootHisto(j["histogram"]).getHistogram()
#        h_data    = dataset_signal.getDataDatasets().getDatasetRootHisto(j["histogram"]).getHistogram()
#        h_faketau = dataset_faketau.getDatasetRootHisto(j["histogram"]).getHistogram()
#        h_signal  = 
#        h_ewk     = dataset_signal

#        h_backgr = h_data

#        p = plots.ComparisonPlot(h_data,h_backgr)

        p = plots.DataMCPlot2(histolist)
#datasets, j["histogram"])

        p.draw()

        lumi = 0.0
        for d in dataset_signal.getDataDatasets():
            print "luminosity",d.getName(),d.getLuminosity()
            lumi += d.getLuminosity()
        print "luminosity, sum",lumi
        histograms.addStandardTexts(lumi=lumi)

        if not os.path.exists(plotDir):
            os.mkdir(plotDir)
        p.save(formats)
"""

def main():

    resultdir = os.getcwd()
    jsonfiles = []
    for arg in sys.argv[1:]:
        if not os.path.exists(arg):
            continue
        if os.path.isdir(arg):
            resultdir = arg
            continue
        with open(os.path.abspath(arg)) as jsonfile:
            try:
                json.load(jsonfile)
                jsonfiles.append(arg)
            except ValueError, e:
                print "Problem loading json file",arg,", please check the file"
                sys.exit()

    if len(jsonfiles) == 0:
        usage()

    for j in jsonfiles:
        plot(resultdir,j)

if __name__ == "__main__":
    main()

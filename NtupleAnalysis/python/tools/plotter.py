###############################################################
### Author .........: Alexandros X. Attikis 
### Institute ......: University Of Cyprus (UCY)
### Email ..........: attikis@cern.ch
###############################################################
### Normal assignment operations (a = b) will simply point the new variable towards the existing object.
### A deep copy ( copy.deepcopy() ) constructs a new compound object and then, recursively, inserts copies into it of the objects found in the original.
### A shallow copy ( copy.copy() ) constructs a new compound object and then (to the extent possible) inserts references into it to the objects found in the original.

###############################################################
### All imported modules
###############################################################
### System
import os, sys
import array
import math
import copy
import inspect
import glob
from optparse import OptionParser
import itertools
import time
import numpy
from array import array
import re
import collections #for ordered dictionaries

### User
import HLTausAnalysis.NtupleAnalysis.tools.datasets as m_datasets
import HLTausAnalysis.NtupleAnalysis.tools.tdrstyle as m_tdrstyle
import HLTausAnalysis.NtupleAnalysis.tools.text as m_text
import HLTausAnalysis.NtupleAnalysis.tools.aux as m_aux
import HLTausAnalysis.NtupleAnalysis.tools.styles as m_styles
import HLTausAnalysis.NtupleAnalysis.tools.histos as m_histos

### ROOT
import ROOT
from ROOT import std

###############################################################
### Define class
###############################################################
class Plotter(object): 
    def __init__(self, Verbose=False, BatchMode=True):
        self.bVerbose          = Verbose
        self.BatchMode         = BatchMode
        self.SetupROOT()       
        self.TDRStyleObject    = m_tdrstyle.TDRStyle()
        self.StyleObject       = m_styles.StyleClass(verbose = self.bVerbose)
        self.TextObject        = m_text.TextClass(verbose=self.bVerbose)
        self.AuxObject         = m_aux.AuxClass(verbose=self.bVerbose)
        self.DatasetObject     = m_datasets.DatasetClass(verbose=self.bVerbose)
        self.CanvasFactor      = 1.25
        self.DatasetToRootFileMap  = collections.OrderedDict()
        self.DatasetToHistoMap     = collections.OrderedDict()
        self.DatasetToLatexNameMap = collections.OrderedDict()
        self.DivisionPoint     = 1-1/self.CanvasFactor
        self.DrawObjectList    = []
        self.DrawObjectListR   = []
        self.HistoObjectList   = []
        self.TCanvas           = None
        self.TLegend           = None
        self.PadCover          = None
        self.PadPlot           = None
        self.PadRatio          = None        
        self.TBoxList          = []
        self.xTLineList        = []
        self.yTLineList        = []
        self.THRatio           = None
        self.THDumbie          = None
        self.THStack           = ROOT.THStack("THStack" + "@" + str(time.time()), "Stack for PadPlot Histograms")
        self.THStackHistoList  = [] #needed because looping over self.THSTack.GetHists() crashes!
        self.THStackRatio      = ROOT.THStack("THStackRatio" + "@" + str(time.time()), "Stack for PadRatio Histograms")
        self.bInvPadRatio      = False
        self.bPadRatio         = False
        self.RatioErrorType    = None
        self.bStackInclusive   = False
        self.TMultigraph       = ROOT.TMultiGraph("TMultigraph"  + "@" + str(time.time()), "ROOT.TMultiGraph holding various ROOT.TGraphs")
        self.startTime         = time.time()        
        self.MsgCounter        = 0
        self.CutOption         = None
        self.IsTH1             = False
        self.IsTH2             = False
        self.CutLineColour     = ROOT.kBlack #ROOT.kGray
        self.SaveContourPoints = False
        self.MaxDecimals       = None
        self.UseDatasetAsLegEntry = False
        self.PileUp            = None
        
    def Verbose(self, messageList=None):
        '''
        Custome made verbose system. Will print all messages in the messageList
        only if the verbosity boolean is set to true.
        '''
        if self.bVerbose == False:
            return
        
        print "%s:" % (self.__class__.__name__ + "." + sys._getframe(1).f_code.co_name + "()")
        if messageList==None:
            return
        else:
            for message in messageList:
                print "\t", message
        return


    def Print(self, messageList=[""]):
        '''
        Custome made print system. Will print all messages in the messageList
        even if the verbosity boolean is set to false.
        '''

        self.MsgCounter = self.MsgCounter  + 1
        print "[%s] %s:" % (self.MsgCounter, self.__class__.__name__ + "." + sys._getframe(1).f_code.co_name + "()")
        for message in messageList:
            print "\t", message
        return


    def SetupROOT(self):
        '''
        Setup ROOT before doing anything else. Reset ROOT settings, disable statistics box,
        apply Techical Design Report (TDR) style.
        '''
        self.Verbose()
        
        ROOT.gROOT.Reset()
        ROOT.gStyle.SetOptStat(0)
        ROOT.TGaxis.SetMaxDigits(4) #18 July 2015 (4)
        #m_tdrstyle.TDRStyle() # fotis
        ROOT.gROOT.SetBatch(self.BatchMode)
        self.SetupROOTErrorIgnoreLevel(2000)
        ROOT.gStyle.SetNumberContours(999)
        return


    def SetupStatsBox(self, xPos=0.94, yPos=0.84, width=0.20, height=0.12, options = 0):
        '''
        The parameter mode can be = ksiourmen  (default = 000001111)
        k = 1 (2); kurtosis printed (kurtosis and kurtosis error printed)
        s = 1 (2); skewness printed (skewness and skewness error printed)
        i = 1 (2); integral of bins printed (integral of bins with option "width" printed)
        o = 1;     number of overflows printed
        u = 1;     number of underflows printed
        r = 1 (2); rms printed (rms and rms error printed)
        m = 1 (2); mean value printed (mean and mean error values printed)
        e = 1;     number of entries printed
        n = 1;     name of histogram is printed
        
        Example: gStyle->SetOptStat(11);
        print only name of histogram and number of entries.
        '''
        self.Verbose()
        self.PadCover

        ### Beautifications/Styling
        ROOT.gStyle.SetStatBorderSize(0)
        ROOT.gStyle.SetStatColor(ROOT.kWhite)
        ROOT.gStyle.SetStatStyle(3001) #3001
        ROOT.gStyle.SetStatTextColor(ROOT.kBlack)
        #ROOT.gStyle.SetStatFont(62)
        ROOT.gStyle.SetStatFontSize(15)
        # ROOT.gStyle.SetStatFormat(const char* format = "6.4g")
        
        ### Dimensions
        ROOT.gStyle.SetStatY(yPos)
        ROOT.gStyle.SetStatX(xPos)
        ROOT.gStyle.SetStatW(width)
        ROOT.gStyle.SetStatH(height)
        ROOT.gStyle.SetOptStat(options)

        ### Workaround for not being able to not remove it
        if yPos == -1 and xPos == -1 and width==-1 and height==-1:
            ROOT.gStyle.SetStatTextColor(ROOT.kWhite)
            ROOT.gStyle.SetStatFontSize(0)
            ROOT.gStyle.SetStatStyle(0)
            
        return


    def SetupROOTErrorIgnoreLevel(self, level):
        '''
        Available options (const Int_t):
        kUnset    =  -1
        kPrint    =   0
        kInfo     =   1000
        kWarning  =   2000
        kError    =   3000
        kBreak    =   4000
        kSysError =   5000
        kFatal    =   6000
        '''
        self.Verbose()

        ROOT.gErrorIgnoreLevel = level
        return
    

    def AddDataset(self, datasetName, rootFilePath):
        '''
        Add a new dataset and associate it to a root file.
        '''
        self.Verbose(["From file '%s', adding dataset '%s'." % (rootFilePath, datasetName)])

        ### Try to get Pile Up
        self.GetPileUp(datasetName.rsplit(":", 1)[0].lower())

        ### Distinguish latex names and datasets
        datasetName      = datasetName.rsplit(":", 1)[-1]
        datasetLatexName = datasetName.rsplit(":", 1)[0]
                
        ### Ensure that dataset is unique (dictionaries must have unique keys)
        if datasetName in self.DatasetToRootFileMap.keys():
            self.Print(["Dataset '%s' already exists, associated with file '%s'" % (datasetName, self.DatasetToRootFileMap[datasetName]), "Please ensure that dataset names used are unique. EXIT"])
            sys.exit()
        elif datasetName in self.DatasetToLatexNameMap.keys():
            self.Print(["Dataset LaTeX name '%s' already exists, associated with dataset '%s'" % (datasetName), "Please ensure that dataset LaTeX names are unique. EXIT"])
            sys.exit()
        else:
            pass

        ### Map dataset to ROOT file, append dataset (name) to datasetlist, append ROOT file (name) to TFile list
        self.DatasetToRootFileMap[datasetName]  = self.GetRootFile(rootFilePath)
        
        ### Map dataset to Latex name
        if datasetLatexName != datasetName:
            self.DatasetToLatexNameMap[datasetName] = datasetLatexName
        else:
            self.DatasetToLatexNameMap[datasetName] = self.DatasetObject.ConvertDatasetToLatex(datasetName)
        return


    def GetPileUp(self, dataset):
        '''
        '''

        if self.PileUp == None:
            self.PileUp = self.DatasetObject.DatasetToPileUpMap[dataset]
            return
        elif self.PileUp != self.DatasetObject.DatasetToPileUpMap[dataset]:
            self.Print(["Datasets with different PU values ('%s' and '%s'). Will display biggest one" % (self.DatasetObject.DatasetToPileUpMap[dataset] , self.PileUp)])
            if self.DatasetObject.DatasetToPileUpMap[dataset] > self.PileUp:
                self.PileUp = self.DatasetObject.DatasetToPileUpMap[dataset]
            return
        else:
            #self.Print(["ERROR!", "This should not print!"])
            return
        return
    
    
    def GetDatasets(self, McProduction, MulticrabDirPath, datasetList=[""]):
        '''
        OBSOLETE
        '''
        self.Verbose()

        ### First setup the datasets object
        self.DatasetObject.SetMulticrabDirPath(MulticrabDirPath)
        self.DatasetObject.SetMcProduction(McProduction)

        ### Get the dataset to file maps
        for dataset in datasetList:
            dataset = dataset.lower()
            if dataset in self.DatasetObject.GetDatasetToFileMap().keys():
                self.DatasetToRootFileMap[dataset] = self.GetRootFile(self.DatasetObject.GetFile(dataset))                
            else:
                raise Exception("The dataset provided ('%s') is not valid. Please select one of the following options:\n\t%s" % (dataset, "\n\t".join(self.DatasetObject.GetDatasetToFileMap().keys())) )
        return


    def SetVerbose(self, verbose):
        '''
        Manually enable/disable verbosity.
        '''
        self.Verbose()

        self.bVerbose = verbose
        self.Verbose(["Verbose mode = ", self.bVerbose])
        return


    def GetRootFile(self, filePath, mode = "READ"):
        '''
        Open a TFile from a given filePath in a given mode (default = "READ") and return it.
        Ensure that file does exist.
        '''

        self.Verbose()

        ### Check that the path actually exists
        if os.path.exists(filePath) == False:
            raise Exception("File '%s' does not exist! Please make sure the provided path for the file is correct." % (filePath) )

        ### Open file
        fileName = filePath.rsplit("/", 1)[-1]
        rootFile = ROOT.TFile.Open( filePath, mode, fileName, 1, 0)

        ### Ensure that file being accessed is indeed a ROOT file. If so return it, else raise exception
        if isinstance(rootFile, ROOT.TFile) == False:
            raise Exception("The file '%s' exists but it is not a ROOT file!" % (rootFile.GetName()) )
        else:
            return rootFile


    def _CreateCanvas(self):
        '''
        Create a name for a TCanvas and then create it. 
        This name will later on be used to save the canvas under a user specific format ("png", "pdf", "eps", etc..)
        '''                    
        self.Verbose()

        ### Normal assignment operations (a = b) will simply point the new variable towards the existing object.
        histo = self.THDumbie 

        ### Create a name for the TCanvas
        if histo.saveName == None:
            canvasName = histo.name
        else:
            canvasName = histo.saveName

        ### Create TCanvas and activate it
        self.Verbose(["Creating canvas with name '%s'" % ( canvasName)])

        ### Add the time in the canvas name to avoid memory duplicates when handling the same histos. Truncate "@time" when saving canvas
        canvasName =  canvasName + "@" + str(time.time())
        self.TCanvas = ROOT.TCanvas( canvasName, canvasName, 1)
        time.sleep(0.2) #to avoid getting the same time (canvas/histo overwrite)

        self._SetLogAxesCanvas()
        self.TCanvas.cd()
        return


    def _Create2PadCanvas(self):
        '''
        Create a 2-pad TCanvas to accomodate a ratio plot and customise it.
        '''                    
        self.Verbose()

        ### Normal assignment operations (a = b) will simply point the new variable towards the existing object.
        histo = self.THDumbie 

        ### Create the THRatio histogram (PadRatio equivalent of PadPlot's THDumbie)
        self.THRatio = copy.deepcopy(histo)
        self.THRatio.TH1orTH2.SetName("THRatio")

        ### Create TCanvas name that included the dataset name
        if histo.saveName == None:
            canvasName = histo.name
        else:
            canvasName = histo.saveName

        ### Create TCanvas and divide it into two pads: one for plot pad, one for ratio pad
        self.Verbose(["Creating canvas with name '%s'" % ( canvasName)])
        canvasName =  canvasName + "@" + str(time.time())
        self.TCanvas = ROOT.TCanvas( canvasName, canvasName, ROOT.gStyle.GetCanvasDefW(), int(ROOT.gStyle.GetCanvasDefH()*self.CanvasFactor))
        self.TCanvas.Divide(1,2)

        ### Remove x-axis title and labels from the THDumbie to avoid overlap with those THRatio
        self._RemoveXaxisBinLabels(histo.TH1orTH2)
        self._RemoveXaxisTitle(histo.TH1orTH2)

        ### Create plot, ratio and cover pads
        self._CreatePads()
        
        ### Take care of log-axis settings
        self._SetLogAxes2PadCanvas()

        ### Update canvas and change back to the PadPlot
        self.TCanvas.Update()
        self.PadPlot.cd()

        return


    def _CreatePads(self):
        '''
        Create the plot, ratio and cover pads.
        '''
        self.Verbose()

        self._CreatePadPlot()
        self._CreatePadRatio()
        self._CreatePadCover()
        return

    
    def _CreatePadCover(self, xMin=0.09, yMin=0.285, xMax=0.16, yMax=0.32):
        '''
        Creates a cover pad to cover the overlap of the y-axis divisions between the PadPlot and the PadRatio.
        '''
        self.Verbose()
        
        self.TCanvas.cd()
        self.PadCover = ROOT.TPad("PadCover", "PadCover", xMin, yMin, xMax, yMax)
        self.PadCover.SetName("PadCover")
        self.PadCover.SetBorderMode(0)
        self.PadCover.SetFillStyle(1001)
        self.PadCover.SetFillColor(ROOT.kWhite) #ROOT.kRed
        self.PadCover.Draw()
        self.PadRatio.Draw() # Re-draw PadRatio to put back the covered y-axis numbers
        self.TCanvas.Update()
        return


    def _CreatePadPlot(self):
        '''
        Creates a plot pad to draw the histogram stack.
        '''
        self.Verbose()

        self.PadPlot  = self.TCanvas.cd(1)
        self.PadPlot.SetName("PadPlot")
        (xlow, ylow, xup, yup) = [ROOT.Double(x) for x in [0.0]*4]
        self.PadPlot.GetPadPar(xlow, ylow, xup, yup)
        self.PadPlot.SetPad(xlow, self.DivisionPoint, xup, yup)
        self.PadPlot.Draw()
        self.TCanvas.Update()
        return


    def _CreatePadRatio(self):
        '''
        Creates a ratio pad to draw the histogram ratio stack.
        '''
        self.Verbose()
        
        CanvasHeightCorr = 0.022

        self.PadRatio = self.TCanvas.cd(2)
        self.PadRatio.SetName("PadRatio")
        (xlow, ylow, xup, yup) = [ROOT.Double(x) for x in [0.0]*4]
        self.PadRatio.GetPadPar(xlow, ylow, xup, yup)
        self.PadRatio.SetPad(xlow, ylow, xup, self.DivisionPoint + ROOT.gStyle.GetPadBottomMargin() - ROOT.gStyle.GetPadTopMargin() + CanvasHeightCorr)
        self.PadRatio.SetFillStyle(4000)
        self.PadRatio.SetTopMargin(0.0)
        self.PadRatio.SetBottomMargin(self.PadRatio.GetBottomMargin()+0.16)
        self.PadRatio.Draw()
        self.TCanvas.Update()
        return


    def _RemoveXaxisBinLabels(self, histo):
        '''
        Removes all the x-axis labels from the histogram pass as argument.
        '''
        self.Verbose()

        bIsTH1 = isinstance(histo, ROOT.TH1)
        if bIsTH1 == False:
            self.Print(["The histogram '%s' is not a ROOT.TH1 instance. Doing nothing" % histo.name])
            return
        histo.GetXaxis().SetLabelSize(0)
        return


    def _RemoveXaxisTitle(self, histo):
        '''
        Removes the x-axis title from the histogram passed as argument.
        '''
        self.Verbose()

        bIsTH1 = isinstance(histo, ROOT.TH1)
        if bIsTH1 == False:
            self.Print(["The histogram '%s' is not a ROOT.TH1 instance. Doing nothing" % histo.name])
            return
        histo.GetXaxis().SetTitleSize(0)
        return


    def _SetLogAxesCanvas(self):
        '''
        Apply axes customisations to a TCanvas.
        '''
        self.Verbose()

        ### Determine whether to set log for y- and x- axis.
        self._SetLogY(self.THDumbie, self.TCanvas)
        self._SetLogX(self.THDumbie, self.TCanvas)
        self._SetLogZ(self.THDumbie, self.TCanvas)
        return


    def _SetLogX(self, histo, PadOrCancas):
        '''
        Determine whether to set log for x-axis.
        '''
        self.Verbose()

        if histo.logX==False:
            return

        ### Set log-scale for x-axis 
        if histo.xMin == None:
            histo.xMin = histo.TH1orTH2.GetXaxis().GetXmin()

        if histo.xMin > 0:
            PadOrCancas.SetLogx(True)
        else:
            raise Exception("Request for TCanvas::SetLogx(True) rejected. The '%s' minimum x-value is '%s'." % (histo.name, histo.xMin))
        return


    def _SetLogY(self, histo, PadOrCancas):
        '''
        Determine whether to set log for y-axis.
        '''
        self.Verbose()

        if histo.logY==False:
            return

        if histo.yMin == None:
            histo.yMin = histo.TH1orTH2.GetMinimum()

        if histo.yMin > 0:
            PadOrCancas.SetLogy(True)
        else:
            raise Exception("Request for TCanvas::SetLogy(True) rejected. The '%s' minimum y-value is '%s'." % (histo.name, histo.yMin))
        return


    def _SetLogZ(self, histo, PadOrCancas):
        '''
        Determine whether to set log for z-axis.
        '''
        self.Verbose()

        if histo.logZ==False or isinstance(histo.TH1orTH2, ROOT.TH2) == False:
            return

        ### Set log-scale for z-axis 
        if histo.zMin == None:
            histo.zMin = histo.TH1orTH2.GetZaxis().GetXmin()

        if histo.zMin > 0:
            PadOrCancas.SetLogz(True)
        else:
            raise Exception("Request for TCanvas::SetLogz(True) rejected. The '%s' minimum z-value is '%s'." % (histo.name, histo.zMin))
        return


    def _SetLogXYRatio(self, histo, PadOrCancas):
        '''
        Determine whether to set log for x-axis.
        '''
        self.Verbose()

        ### Set log-scale for y-axis 
        if histo.logYRatio==True:
            if histo.yMin == None:
                histo.yMin = histo.TH1orTH2.GetMinimum()
            
            if histo.yMin > 0:
                PadOrCancas.SetLogy(True)
            else:
                raise Exception("Request for TCanvas::SetLogy(True) rejected. The '%s' minimum y-value is '%s'." % (histo.name, histo.yMin))

        ### Set log-scale for x-axis 
        if histo.logXRatio==True:
            if histo.xMin == None:
                histo.xMin = histo.TH1orTH2.GetXaxis().GetXmin()

            if histo.xMin > 0:
                PadOrCancas.SetLogx(True)
            else:
                raise Exception("Request for TCanvas::SetLogx(True) rejected. The '%s' minimum x-value is '%s'." % (histo.name, histo.xMin))
        return


    def _SetLogAxes2PadCanvas(self):
        '''
        Apply customisations to a 2-pad TCanvas.
        '''
        self.Verbose()

        ### Determine whether to set log for y- and x- axis.
        self._SetLogY(self.THDumbie, self.PadPlot)
        self._SetLogX(self.THDumbie, self.PadPlot)
        self._SetLogZ(self.THDumbie, self.PadPlot)
        self._SetLogXYRatio(self.THRatio, self.PadRatio)
        return


    def _CreateLegend(self):
        '''
        Create a TLegend, customise it and return it.
        '''
        self.Verbose()
        
        if isinstance(self.TLegend, ROOT.TLegend) == True:
            #raise Exception("A TLegend has already been created. Only 1 allowed per Plotter object.")
            return
        else:
            histo = self.THDumbie
            self.TLegend = ROOT.TLegend(histo.xLegMin, histo.yLegMin, histo.xLegMax, histo.yLegMax, self._GetLegendHeader(), "brNDC")
            self._CustomiseLegend()
            self.DrawObjectList.append( self.TLegend )
        return


    def _CustomiseLegend(self):
        '''
        Customise a TLegend.
        '''
        self.Verbose()
        self.TLegend.SetName("legend_" + str(time.time()) )
        self.TLegend.SetFillStyle(0)
        self.TLegend.SetLineColor(ROOT.kBlack)
        self.TLegend.SetLineWidth(1)
        self.TLegend.SetBorderSize(0)
        self.TLegend.SetShadowColor(ROOT.kWhite)
        self.TLegend.SetTextSize(0.03)
        self.TLegend.SetTextFont(62)
        return

    

    def AddHisto(self, histoObject):
        '''
        '''
        self.Verbose()
        
        if type(histoObject) == list:
            for h in histoObject:
                self._AddHistoToQueue(h)
        elif isinstance(histoObject, m_histos.TH1orTH2):
            self._AddHistoToQueue(histoObject)
        else:
            self.Print(["ERROR!", "Expected object of type '%s' but got '%s' instead" % ( m_histos.TH1orTH2, type(histoObject) ), "EXIT"])
            sys.exit()


    def _AddHistogramsToStack2D(self, ProfileAxis, firstBin, LastBin):
        '''
        Add all histograms (except Dumbie) to a THStack. For each histogram add a TLegend entry
        and automatically extend the size of the TLegend to accomodate the next entry.
        '''
        self.Verbose()

        entryLabel = ""
        if  ProfileAxis =="x" or  ProfileAxis =="y":
            bAddLegendEntries = True
        else:
            bAddLegendEntries = False
        self.TCanvas.SetName( self.TCanvas.GetName() + "_Profile%s" % (ProfileAxis.upper()) )

        ### For-loop: Histos
        for h in self.DatasetToHistoMap.values():
            
            self.Verbose( ["Creating Profile%s for histogram '%s'  and adding to THStack." % ( ProfileAxis.upper(), h.name)] )
            if ProfileAxis == "x":
                hProfileX_Name = h.name +"_ProfileX"
                hProfileX      = h.TH1orTH2.ProfileX(hProfileX_Name, firstBin, LastBin)
                self.THStack.Add( hProfileX )
            else:
                #raise Exception("Although this works, some validation would have to be carried out with a simple a well undestood 2D histo")
                self.Print(["WARNING! Although this works, some validation would have to be carried out with a simple a well undestood 2D histo"])
                hProfileY_Name = h.name +"_ProfileY"
                hProfileY      = h.TH1orTH2.ProfileY(hProfileY_Name, firstBin, LastBin)
                self.THStack.Add( hProfileY )
                
            ### Add legend entries for THStack?
            if bAddLegendEntries == True:
                self.TLegend.AddEntry( h.TH1orTH2, self._GetLegEntryLabel(h), h.legOptions)
                self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)
            else:
                pass
        return


    def _AddHistoToQueue(self, histoObject):
        '''
        '''
        self.Verbose()

        ### Sanity check. At least one dataset is present 
        if len(self.DatasetToRootFileMap.keys()) < 1:
            self.Print(["ERROR!", "No datasets found! Exit"])
            sys.exit()
        else:
            pass        

        ### Ensure that the pass argument is a valid histo object
        self.IsValidHistoObject(histoObject)

        ### Loop over all datasets
        for dataset in self.DatasetToRootFileMap.keys():

            ### Get the ROOT file for given dataset
            f = self.DatasetToRootFileMap[dataset]

            ### Ensure that histogram exists in the file
            self.CheckThatHistoExists(f, histoObject)

            ### Map histo to dataset
            mapKey = dataset + ":" + histoObject.name
            self.Verbose(["Mapping '%s' to key '%s'." % ( histoObject.name, mapKey)])
            self.DatasetToHistoMap[mapKey] = copy.deepcopy(histoObject)
                
            ### Assign user-defined attributes to histo object
            histoObject = self.DatasetToHistoMap[mapKey] 
            self._AssignHistoObjectAttributes( f, dataset, histoObject)

            ### Print extensive histogram information
            self.PrintHistoInfo(histoObject, False)

            ### Append to histoObject list
            self.HistoObjectList.append(histoObject)

        return

    
    def _AssignHistoObjectAttributes(self, rootFile, dataset, histoObject):
        '''
        '''
        self.Verbose()
        
        ### Declare shorter name references
        h = histoObject
        f = rootFile

        ### Assign attributes
        h.TH1orTH2      = self.GetHistoFromFile(f, h)
        h.TFileName     = f.GetName()
        h.dataset       = dataset
        h.rangeIntegral = h.TH1orTH2.Integral()

        ### Determine the histogram integral
        if self.IsTH2:
            self.TDRStyleObject.setWide(True)
            self._CheckThatNoTH2WithMoreThanOneDatasets(h)
            h.integral  = h.TH1orTH2.Integral(0, h.TH1orTH2.GetNbinsX()+1, 0, h.TH1orTH2.GetNbinsY()+1)
        else:
            h.integral  = h.TH1orTH2.Integral(0, h.TH1orTH2.GetNbinsX()+1)
            
        ### Assign global values
        self.bPadRatio      = h.ratio
        self.bInvPadRatio   = h.invRatio
        self.RatioErrorType = h.ratioErrorType
            
        return
        


    def GetHistoFromFile(self, f, histo):
        '''
        '''
        self.Verbose()
        
        histoPath = ""
        if histo.path == "" or histo.path==None:
            histoPath = histo.name
        else:
            histoPath = histo.path + "/" + histo.name

        return f.Get(histoPath)


    def PrintHistoInfo(self, histo, bVerbose=False):
        '''
        '''
        self.Verbose()
                
        if histo.TH1orTH2.Integral() == 0:
            self.Print(["WARNING!", "File: '%s'" % (histo.TFileName), "Dataset: '%s'" % (histo.dataset), "HistoPath: '%s'" % (histo.path), "HistoName: '%s'" % (histo.name), 
                        "Integral(): '%s'" % (histo.rangeIntegral), "Integral(0, nBins+1): '%s'" % (histo.integral), "normaliseTo: '%s'" % (histo.normaliseTo)])
        elif bVerbose == True:
            self.Print(["File: '%s'" % (histo.TFileName), "Dataset: '%s'" % (histo.dataset), "HistoPath: '%s'" % (histo.path), "HistoName: '%s'" % (histo.name), 
                        "Integral(): '%s'" % (histo.rangeIntegral), "Integral(0, nBins+1): '%s'" % (histo.integral), "normaliseTo: '%s'" % (histo.normaliseTo)])
        else:
            self.Verbose(["File: '%s'" % (histo.TFileName), "Dataset: '%s'" % (histo.dataset), "HistoPath: '%s'" % (histo.path), "HistoName: '%s'" % (histo.name), 
                        "Integral(): '%s'" % (histo.rangeIntegral), "Integral(0, nBins+1): '%s'" % (histo.integral), "normaliseTo: '%s'" % (histo.normaliseTo)])
        return


    def ConvertToOneMinusCumulativeHistos(self):
        '''
        This method converts all histograms into a (1-cumulative integral) histograms.
        '''
        self.Verbose()
        
        ### For-loop: All datasets
        for dataset in self.DatasetToHistoMap.keys():
            h = self.DatasetToHistoMap[dataset]
            self._ConvertToOneMinusCumulativeHisto(h)
        return

    
    def GetROCHistoLists(self, rateToEffMap):
        '''
        '''
        self.Verbose()

        rateHistoList = []
        effHistoList  = []

        ### Get match for rate and eff histos
        for key in sorted(rateToEffMap):
            rateHistoList.append(key)
            effHistoList.append(rateToEffMap[key])
            trigger1 = key.name.rsplit("_", -11)[0]
            trigger2 = rateToEffMap[key].name.rsplit("_", -1)[0]

            ### Make sure you match rate and eff histos of same trigger
            if trigger1 != trigger2:
                raise Exception("Trigger names differ. This should not happen!")

        ### Sanity Check
        if len(rateHistoList) != len(effHistoList):
            raise Exception("Length of RateHistoList ('%s') and  EffHistoList ('%s') is not the same! Cannot produce ROC curves." % (len(RateHistoList), len(EffHistoList)))
        return rateHistoList, effHistoList


    def ConvertToROC(self, rateHistoList, effHistoList, decimals=3, bDrawZValues=False, **kwargs):
        '''
        '''
        self.Verbose()

        for hRate, hEff in zip(rateHistoList, effHistoList):
            self.ConvertHistosToROC( hRate, hEff, bDrawZValues, decimals, **kwargs)
        return


    def ConvertHistosToROC(self, rateHisto, effHisto, bDrawZValues=False, decimals=3, **kwargs):
        '''
        Will create TGraphsErrors.
        '''
        self.Verbose()

        ### Enable colour palette (change colour for same datasample, don't just use shades of same colour)
        self.StyleObject.EnableColourPalette(True)
        
        ### Definitions
        et       = []
        etUp     = []
        etLow    = []
        rate     = []
        rateUp   = []
        rateLow  = []
        eff      = []
        effUp    = []
        effLow   = []

        ### Sanity check
        trigger1 = rateHisto.name.rsplit("_", -1)[0]
        trigger2 = effHisto.name.rsplit("_", -1)[0]
        if trigger1 != trigger2:
            raise Exception("Trigger names differ ('%s' != '%s'). This should not happen!" % (trigger1, trigger2) )

        ### For-loop: Rate histos
        self.PrintHistoInfo(rateHisto, False)
        self.PrintHistoInfo(effHisto, False)
                
        ### Tread 2d ROCs differently
        hRate = None
        hEff  = None
        if ("TH2" in str( type(rateHisto.TH1orTH2)) ):
            hRate = self.GetDiagonalOf2DHistoAndFill1DHisto(rateHisto)
        else:
            hRate = rateHisto.TH1orTH2

        if ("TH2" in str( type(effHisto.TH1orTH2)) ):
            hEff = self.GetDiagonalOf2DHistoAndFill1DHisto(effHisto)
        else:
            hEff = effHisto.TH1orTH2

        ### For-loop: Target rates
        for target_rate in self.GetRateValuesList():
            iBin, calo_et, mb_rate, mb_rateErr  = self.FindFirstBinBelowValue(hRate, target_rate)
            if (iBin == -1):
                continue
            et.append     ( round(calo_et, decimals) )
            etUp.append   ( round(0.0, decimals) )
            etLow.append  ( round(0.0, decimals) )
            rate.append   ( round(mb_rate, decimals) )
            rateUp.append ( round(mb_rateErr, decimals) )
            rateLow.append( round(mb_rateErr, decimals) )
            #rateUp.append ( round(math.sqrt(mb_rate), decimals) ) # conservative approach
            #rateLow.append( round(math.sqrt(mb_rate), decimals) ) # conservative approach
            self.Verbose([ "Target (kHz): '%s'" % (mb_rate), "Rate: '%s + %s - %s (kHz)'" % (mb_rate, math.sqrt(mb_rate), math.sqrt(mb_rate)), "Et: '%s'" % (et) ])
            

        ### Sanity Check
        if (len(et) == len(etUp) == len(etLow) == len(rate) == len(rateUp) == len(rateLow) ):
            pass
        else:
            self.Print(["ERROR! Arrays have different length! EXIT"])
            sys.exit()
        
        ### For-loop: Et (corresponding to target rates)
        for calo_et in et:
            etBin   = hEff.FindBin(calo_et) 
            eff.append   ( round(hEff.GetBinContent(etBin), decimals) )
            effUp.append ( round(hEff.GetBinError(etBin), decimals) )
            effLow.append( round(hEff.GetBinError(etBin), decimals) )
            self.Verbose([ "Dataset:'%s'" % (effHisto.dataset), "HistoName : '%s'" % (effHisto.name), "Trigger: '%s'" % (trigger1), 
                         "Et  (GeV)  : '%s'" % (et), 
                         "Et+ (GeV)  : '%s'" % (etUp), 
                         "Et- (GeV)  : '%s'" % (etLow), 
                         "Rate  (kHz): '%s'" % (rate),
                         "Rate+ (kHz): '%s'" % (rateUp), 
                         "Rate- (kHz): '%s'" % (rateLow), 
                         "Efficiency : '%s'" % (eff), 
                         "Efficiency+: '%s'" % (effUp),  
                         "Efficiency-: '%s'" % (effLow)])

        ### Sanity check
        if (len(eff) == len(effUp) == len(effLow) ):
            pass
        else:
            self.Print(["ERROR! Arrays have different length! EXIT"])
            sys.exit()

        ### Use values to create a TGraphErrors
        self.AddTGraphErrors(effHisto, eff, effUp, effLow, rate, rateUp, rateLow, et, etUp, etLow, bDrawZValues, **kwargs)
        return


    def GetDiagonalOf2DHistoAndFill1DHisto(self, histo):
        '''
        '''
        self.Verbose()

        nBinsX = histo.TH1orTH2.GetNbinsX()+1

        ### Reset Integral, Contents, Errors and Statistics
        tmpHisto = copy.deepcopy( histo.TH1orTH2.ProjectionX(histo.name + "_ProjectionX", 0, -1, "") )
        tmpHisto.Reset("ICES")

        ### For-loop: x-axis (get y=x bins only)
        for bx in range(0, nBinsX):
            tmpHisto.SetBinContent( bx, histo.TH1orTH2.GetBinContent(bx, bx) )
            tmpHisto.SetBinError( bx, histo.TH1orTH2.GetBinError(bx, bx) )
        
        return tmpHisto


    def DrawEfficiency(self, cutDirection=">", errType="binomial"):
        '''
        '''
        self.Verbose()
        
        self._CustomiseHistograms()

        ### Change some histoObject attributes
        histo = self.HistoObjectList[0]
        saveName = histo.saveName
        kwargs   = histo.kwargs
        binWidthX = histo.binWidthX
        if histo.binWidthX == None:
            binWidthX = histo.TH1orTH2.GetBinWidth(0)
        
        kwargs["yUnits"]      = ""
        kwargs["logX"]        = False
        kwargs["logY"]        = False
        kwargs["yMin"]        = 0.0
        kwargs["yMax"]        = 1.15
        kwargs["normaliseTo"] = ""
        yTitleOld             = kwargs["yLabel"].rsplit("/", 1)[0]
        kwargs["yLabel"]      = kwargs["yLabel"].replace(yTitleOld, "Efficiency (" + cutDirection +  ") ")
        kwargs["yLabel"]      = kwargs["yLabel"]  % (binWidthX)  + " " + histo.xUnits
        kwargs["drawOptions"] = "AP" #"ACE3"
        kwargs["legOptions"]  = "LP" #"FL"

        ### Change histo saveName according to cut-direction
        if cutDirection == ">":
            saveName = saveName + "_GreaterThan"
        elif cutDirection == "<":
            saveName= saveName + "_LessThan"
        else:
            raise Exception("Invalid cut-direction '%s' selected for efficiency plot. Please select either '<' or '>'." % (cutDirection) )

        self._ConvertToEfficiencyHistos(cutDirection, errType, **kwargs)
        self.DrawMultigraph(saveName, **kwargs)
        return


    def _ConvertToEfficiencyHistos(self, cutDirection, errType="binomial", **kwargs):
        '''
        '''
        self.Verbose()
        
        for dataset in self.DatasetToHistoMap:
            h = self.DatasetToHistoMap[dataset]
            self._ConvertToEfficiencyHisto(h, cutDirection, errType, **kwargs)
        return


    def _ConvertToEfficiencyHisto(self, histo, cutDirection, errType="binomial", **kwargs):
        '''
        Replaces bin content with the efficiency of the given bin. Cut direction can be chosen.
        '''
        self.Verbose()
        
        ### Declare lists
        xVals   = []
        xLow    = []
        xUp     = []
        effVals = []
        effLow  = []
        effUp   = []

        ### For-loop: Histo Bins
        nBinsX  = histo.TH1orTH2.GetNbinsX()+1
        for b in range(0, nBinsX+1):

            binWidth   = histo.TH1orTH2.GetBinWidth(b)
            binCenter  = histo.TH1orTH2.GetBinCenter(b)
            binLowEdge = histo.TH1orTH2.GetBinLowEdge(b)
            binUpEdge  = binCenter + binWidth/2

            nPass      = histo.TH1orTH2.Integral(b+1, nBinsX) #events that pass the up-edge of the bin
            nTotal     = histo.TH1orTH2.Integral( 0, nBinsX )
            eff        = -1.0

            ### Calculate the efficiency and its error
            eff, err = self.AuxObject.Efficiency(nPass, nTotal, errType)
            if (cutDirection == ">"):
                pass
            elif (cutDirection == "<"):
                eff = 1-eff
            else:
                self.Print(["ERROR! Illegal logic operator ('%s') for cut-direction. EXIT" % (cutDirection)] )
                sys.exit()

            self.Verbose(["bin = %s, x %s %s,  eff = %s +/- %s" % (b, cutDirection, binUpEdge, eff, err)])
            
            ### Save into lists
            xVals.append(binUpEdge)
            xUp.append(0.0)
            xLow.append(0.0)

            effVals.append(eff)
            effUp.append(err)
            effLow.append(err)


        ### Use values to create a TGraphErrors
        histo.TH1orTH2.SetMaximum(1.0)
        self.AddTGraphErrors(histo, xVals, xUp, xLow, effVals, effUp, effLow, None, None, None, False, **kwargs)
        return


    def _ConvertToOneMinusCumulativeHisto(self, histo):
        '''
        This method convert a given histogram to a (1-cumulative integral) histogram. 
        '''
        self.Verbose()

        ### Sanity check: Ensure that no bin-resining has been requested
        if (histo.binWidthX != None):
            self.Print(["ERROR! The requested x-axis bin width ('%s') is invalid. Please set both axes bin-widths to 'None'. EXIT" % histo.binWidthX] )
            sys.exit()
        elif (histo.binWidthY != None):
            self.Print(["ERROR! The requested y-axis bin width ('%s') is invalid. Please set both axes bin-widths to 'None'. EXIT" % histo.binWidthY] )
            sys.exit()
        else:
            pass
        
        ### Loop over all histogram bins
        h = histo.TH1orTH2
        self.Print( ["Converting histo '%s' into a (1-cumulative integral) histogram" % ( h.GetName() )] )

        if isinstance(h, ROOT.TH1D) == True:
            nBins = h.GetNbinsX()+1
            for b in range(0, nBins+1):
                value = h.Integral(b, nBins) 
                error = math.sqrt(value)
                h.SetBinContent( b, value )
                h.SetBinError  ( b, error )
            return
        elif isinstance(h, ROOT.TH2D) == True:
            nBinsX = h.GetNbinsX()+1
            nBinsY = h.GetNbinsY()+1
            for bX in range(0, nBinsX+1):
                for bY in range(0, nBinsY+1):
                    # Draw only y<= x values
                    if bY >= bX and self.SaveContourPoints == False: 
                        break
                    else: #if you want to calculate thresholds proceed as normal
                        pass
                    value = h.Integral(bX, nBinsX, bY, nBinsY)
                    error = math.sqrt(value)
                    h.SetBinContent( bX, bY, value)
                    h.SetBinError  ( bX, bY, error)
            return
        else:
            raise Exception("Cannot only convert a TH1D or a TH2D into a rate histogram, but instead got an object of type '%s'." % (type(h) ) )
        return


    def _CheckThatNoTH2WithMoreThanOneDatasets(self, histoObject):
        ''''
        Ensure that no TH2 is drawn where more than 1 dataset is created. Very difficult to distinguish so we need to take care.
        '''
        self.Verbose()
        
        hType = type(histoObject.TH1orTH2)
        
        if len(self.DatasetToRootFileMap.keys())>1 and "TH2" in str(hType):
            raise Exception("Cannot draw a TH2 while more than 1 datasets are present.")
        else:
            return


    def IsTypeTH2(self, histo):
        '''
        Check if histoObject is of type ROOT.TH2
        '''

        hType = type(histo)

        if "TH2" in str(hType):
            return True
        else:
            return False


    def IsValidHistoObject(self, histoObject):
        ''''
        Ensure that the histoObject is of valid type (m_histos.TH1 or m_histos.TH2). Raise an exception otherwise.
        '''
        self.Verbose()

        if isinstance(histoObject, m_histos.TH1orTH2):
            return
        else:
            self.Print(["ERROR!", "Unknown histo type. Please make sure the histo object '%s' (type = '%s') is either a TH1 or a TH2" % (histoObject, type(histoObject)), "EXIT"])
            sys.exit()            


    def CheckThatHistoExists(self, rootFile, hObject):
        '''
        Ensure that the histogram you are trying to get from a TFile really exists.
        '''
        self.Verbose()
        
        hPath = ""
        if (hObject.path == "") or (hObject.path == None):
            hPath = hObject.name
        else:
            hPath = hObject.path + "/" + hObject.name

        self.IsTH1 = isinstance( rootFile.Get(hPath), ROOT.TH1D)
        self.IsTH2 = isinstance( rootFile.Get(hPath), ROOT.TH2D)

        if not self.IsTH1 and not self.IsTH2:
            raise Exception( "Could not find histo object '%s' in TFile '%s' under folder '%s'." % (hObject.name, rootFile.GetName(), hObject.path) )
        else:
            self.Verbose(["File: '%s'" % (rootFile.GetName()), "Histo: '%s'" % (hPath), "IsTH1: '%s'" % (self.IsTH1), "IsTH2: '%s'" % (self.IsTH2)])

        return


    def CreateDumbie(self, THDumbie=None):
        '''
        Create a dumbie histogram that will be the first to be drawn on each canvas. 
        This should have zero entries but have exactly the same attribues  (binning, axes titles etc..) as the ones to be drawn.
        '''
        self.Verbose()

        myMax = -1E10
        myMin = +1E10
        ### Loop over all TH1's and determine dataset histo with yMax
        for dataset in self.DatasetToHistoMap.keys():
            h = self.DatasetToHistoMap[dataset]
            tmpMax =  h.TH1orTH2.GetMaximum()
            if tmpMax > myMax:
                myMax = h.yMax
                myMin = h.yMin
                self.THDumbie = copy.deepcopy(self.DatasetToHistoMap[dataset])                
                self.THDumbie.TH1orTH2.SetName("THDumbie")
            else:
                continue

        ### Reset only Integral, Contents, Errors and Statistics (not Minimum and Maximum)
        self.THDumbie.TH1orTH2.Reset("ICES")

        ### Set custom range for x- and y- axis and pad margins            
        self.THDumbie.TH1orTH2.GetYaxis().SetRangeUser(myMin, myMax)
        self.THDumbie.TH1orTH2.GetXaxis().SetRangeUser(h.xMin, h.xMax) #does nothing if xMax > max x-value when histogram was created

        ### Set Number of divisions! 
        if (self.IsTH2):
            self.THDumbie.TH1orTH2.GetXaxis().SetNdivisions(510) 
        else:
            #self.THDumbie.TH1orTH2.GetXaxis().SetNdivisions(510) #default
            self.THDumbie.TH1orTH2.GetXaxis().SetNdivisions(505)

        ### Set Line Colour and Width
        self.THDumbie.TH1orTH2.SetLineColor(ROOT.kBlack)
        self.THDumbie.TH1orTH2.SetLineWidth(1)

        ### Increase right pad margin to accomodate z-axis scale and title
        if isinstance(self.THDumbie.TH1orTH2, ROOT.TH2) == True:
            ROOT.gStyle.SetPadRightMargin(0.15)
            #ROOT.gStyle.SetPadRightMargin(0.15)
        return
                

    def AppendToDrawObjectList(self, objectToBeDrawn):
        '''
        Append a drawable object of any type (TCanvas, TLegend, TLine, TBox, etc..) to a list.
        This list will be used later on to draw all objects.
        '''
        self.Verbose()

        self.DrawObjectList.append(objectToBeDrawn)
        self.DrawObjectListR.append(copy.deepcopy(objectToBeDrawn))
        return 

        
    def EnableColourPalette(self, bEnable=False):
        '''
        Changes colour for each histogram within a given dataset only if 1 dataset is present.
        '''
        self.Verbose()
        
        self.StyleObject.EnableColourPalette(bEnable)
        return


    def _CustomiseHistograms(self):
        '''
        Customise all histograms found inside the DatasetToHistoMap.
        '''
        self.Verbose()
        
        ### For-loop: All datasets
        for histo in self.HistoObjectList:
            histo.ApplyStyles( self.StyleObject, type(histo.TH1orTH2))
        return

    
    def _NormaliseHistograms(self):
        '''
        Normalise all histograms found inside the DatasetToHistoMap.
        '''
        self.Verbose()
        
        ### Loop over all histograms and normalise according to user options
        for mapKey in self.DatasetToHistoMap:
            h = self.DatasetToHistoMap[mapKey]
            self._NormaliseHisto(h)
        return


    def _NormaliseHisto(self, h):
        '''
        Normalise the histoObject passed to this function according to user-specified criteria. 
        '''
        self.Verbose()
        
        if h.normaliseTo==None:
            return

        scaleFactor = 1
        if h.TH1orTH2.GetEntries() == 0:
            self.Print(["WARNING! Cannot normalise histogram.", "HistoName: '%s'" % (h.name), "Entries: '%s'" % (h.TH1orTH2.GetEntries()), "TFile: '%s'" % (h.TFileName)])
            return
        
        if h.normaliseTo == "One":
            scaleFactor = h.integral #Note: Using h.rangeIntegral is wrong, as it might depend on histogram binning and maximum of x-axis!
            if scaleFactor!=0:
                h.scaleFactor = 1.0/scaleFactor
                h.TH1orTH2.Scale(h.scaleFactor)
            else:
                self.Print(["WARNING! Cannot normalise histogram. Will do nothing.", "HistoName: '%s'" % (h.name), "TFile: '%s'" % (h.TFileName), "ScaleFactor: '%s'" % (scaleFactor)])
                return
        elif h.normaliseTo == "MB@40MHz":
            ConvertRateTokHz  = 1.0E-3
            CrossingRate      = 40.0E+6 #Hz
            MCEventsTotal     = self.DatasetObject.GetEvents(h.dataset)
            scaleFactor       = (CrossingRate)/(MCEventsTotal)*ConvertRateTokHz
            h.scaleFactor     = scaleFactor
            h.TH1orTH2.Scale(scaleFactor)
        elif h.normaliseTo == "MB@30MHz":
            ConvertRateTokHz  = 1.0E-3
            CrossingRate      = 30.0E+6 #Hz
            MCEventsTotal     = self.DatasetObject.GetEvents(h.dataset)
            scaleFactor       = (CrossingRate)/(MCEventsTotal)*ConvertRateTokHz
            h.scaleFactor     = scaleFactor
            h.TH1orTH2.Scale(scaleFactor)
        elif type(h.normaliseTo) == float:
            h.scaleFactor     = float(h.normaliseTo)
            h.TH1orTH2.Scale(h.scaleFactor)
        else:
            raise Exception("Unknown histoObject normalisation option '%s'.!" % (h.normaliseTo))
        
        self.Verbose(["File: '%s'" % (h.TFileName), "Dataset: '%s'" % (h.dataset), "HistoPath: '%s'" % (h.path), "HistoName: '%s'" % (h.name), 
                      "Integral(): '%s'" % (h.rangeIntegral), "Integral(0, nBins+1): '%s'" % (h.integral), "normaliseTo: '%s'" % (h.normaliseTo)])
        return


    def Draw(self, THStackDrawOpt="nostack", bStackInclusive=False, bAddReferenceHisto=True):
        '''
        Draw all necessary histograms for all datasets.
        '''
        self.Verbose()

        self.bStackInclusive = bStackInclusive
        self._NormaliseHistograms()
        self._CustomiseHistograms()
        self._CreateCanvasAndLegendAndDumbie()
        self._CheckHistogramBinning()
        self._AddHistogramsToStack()
        self.Verbose(["Drawing histo '%s'" % (self.THDumbie.name), "THStackDrawOpt: '%s'" % (THStackDrawOpt), "bStackInclusive: '%s'" % (bStackInclusive)])
        self._DrawHistograms(THStackDrawOpt)
        self._DrawRatioHistograms(bAddReferenceHisto)
        self._DrawNonHistoObjects()
        self._CustomiseStack()
        #self.THStack.Draw("same")            #new: needed when drawing cut-boxes
        self.TLegend.Draw("same")            #new: needed when drawing cut-boxes
        #self.AddDefaultCmsText()             #new: needed when drawing cut-boxes
        self.THDumbie.TH1orTH2.Draw("same")
        return


    def ConvertHistosToEfficiency(self, cutDirection=">",  errType = "binomial",  **kwargs):
        '''
        Draw all necessary histograms for all datasets.
        '''
        self.Verbose()

        self._ConvertToEfficiencyHistos(cutDirection, errType, **kwargs)
        return


    def DrawSame(self, HistoObjectList, TLegendHeader=""):
        '''
        This was designed to be used  in conjuction with GetHistos(). 
        For example:

        p1 = m_plotter.Plotter( bVerbose, bBatchMode )
        for dataset in datasetList:
           p1.AddDataset(dataset, datasetPaths[dataset])
        p1.AddHisto(hList1)
        p1.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
        
        p2 = m_plotter.Plotter( bVerbose, bBatchMode )
        for dataset in datasetList:
           p2.AddDataset(dataset, datasetPaths_D[dataset])
        p2.AddHisto(hList2)
        p2.Draw(THStackDrawOpt="nostack", bStackInclusive = False)
        p2.DrawSame(p1.GetHistos())
        p2.SaveHistos(bSavePlots, mySavePath, mySaveFormats)
        '''
        self.Verbose()

        for h in HistoObjectList:
            self.IsValidHistoObject(h)
            h.ApplyStyles( self.StyleObject, type(h.TH1orTH2))
            h.TH1orTH2.Draw(h.drawOptions + ",9same,")
            self.TLegend.AddEntry( h.TH1orTH2, h.legTitle, self._GetLegEntryOptions(h) )
            self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)

        self.TLegend.SetHeader(TLegendHeader)
        self.THDumbie.TH1orTH2.Draw("same")
        return


    def Draw2DProfile(self, THStackDrawOpt="nostack", bStackInclusive=False, ProfileAxis=None, firstBin=1, lastBin=-1):
        '''
        Draw a ProfileX or ProfileY of a TH2D. Basically, this will plot a weighted 2D histo with a single entry replacing 
        all entries in X (or Y axis). The entry that replaces all other entries in the Profile direction is the average.

        Profile histograms are used to display the mean value of Y and its error for each bin in X. The displayed error is by default the
        standard error on the mean (i.e. the standard deviation divided by the sqrt(n) ). Profile histograms are in many cases an elegant 
        replacement of two-dimensional histograms : the inter-relation of two histogram or scatter-plot; its representation on the line-printer 
        is not particularly satisfactory, except for sparse data. If Y is an unknown (but single-valued) approximate function of X, this function
        is displayed by a profile histogram with much better precision than by a scatter-plot.
        See: http://root.cern.ch/root/html/TProfile.html
        '''
        self.Verbose()

        allowedValues = [None, "x", "y"]
        if ProfileAxis not in allowedValues:
            raise Exception("Invalid ProfileAxis option selected ('%s'). You need to speficy the axis of the Profile (x or y). Available options are 'x' and 'y'." % (ProfileAxis) )

        self.bStackInclusive = bStackInclusive
        self.EnableColourPalette(True)
        self._NormaliseHistograms()
        self._CustomiseHistograms()
        self._CreateCanvasAndLegendAndDumbie()
        self._CheckHistogramBinning()
        self._AddHistogramsToStack2D(ProfileAxis, firstBin, lastBin)
        self.Verbose(["Drawing histo '%s'" % (self.THDumbie.name), "THStackDrawOpt: '%s'" % (THStackDrawOpt), "bStackInclusive: '%s'" % (bStackInclusive)])
        self._DrawHistograms(THStackDrawOpt)
        self._DrawRatioHistograms()
        self._DrawNonHistoObjects()
        self._CustomiseStack()
        self.THDumbie.TH1orTH2.Draw("same")
        return


    def AddTF1(self, myFunction, xMin, xMax, kwargs={}):
        '''
        '''
        self.Verbose()

        f1 = ROOT.TF1("f1", myFunction, xMin, xMax)        

        ### Customise Line Style
        if kwargs.get("fillColour"):
            f1.SetFillColor( kwargs.get("fillColour") )
        if kwargs.get("fillStyle"):
            f1.SetFillStyle( kwargs.get("fillStyle") )
        if kwargs.get("lineColour"):
            f1.SetLineColor( kwargs.get("lineColour") )
        if kwargs.get("lineStyle"):
            f1.SetLineStyle( kwargs.get("lineStyle") )
        if kwargs.get("lineWidth"):
            f1.SetLineWidth( kwargs.get("lineWidth") )

        self.AppendToDrawObjectList(f1)
        return


    def _AddHistogramsToStack(self):
        '''
        Add all histograms (except Dumbie) to a THStack. For each histogram add a TLegend entry
        and automatically extend the size of the TLegend to accomodate the next entry.
        '''
        self.Verbose()

        bAddLegendEntries = "TH1" in str(type(self.THDumbie.TH1orTH2)) and isinstance(self.TLegend, ROOT.TLegend)
        entryLabel = ""
        
        ### Loop over all histo objects
        for histo in self.HistoObjectList:
            
            dataset = histo.dataset + ":" + histo.name
            
            h = self.DatasetToHistoMap[dataset]
            self.Verbose(["Adding histogram '%s' to the THStack" % (histo.name), 
                          "Dataset: '%s'" % (dataset), "Integral(): '%s'" % (histo.rangeIntegral), "Integral(0, nBins+1): '%s'" % (histo.integral), "normaliseTo: '%s'" % (histo.normaliseTo)])

            self.THStack.Add(h.TH1orTH2)
            self.THStackHistoList.append(h.TH1orTH2)
            
            ### Add legend entries for THStack
            if bAddLegendEntries == False:
                continue

            ### Add legend entries only for TH1 type histos        
            if h.legTitle != None:
                self.TLegend.AddEntry( h.TH1orTH2, self._GetLegEntryLabel(h), self._GetLegEntryOptions(h) )
                self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)
        return


    def _CustomiseStack(self):
        '''
        Customise the THStack. Apply x- and y-axis range. Can also implement the inclusive error bar (if histos stacks) in the future,
        by cloning the stack and drawing only the errors with "E1" option?
        '''
        self.Verbose()
        self.THStack.GetYaxis().SetRangeUser(self.THDumbie.yMin, self.THDumbie.yMax)
        self.THStack.GetXaxis().SetRangeUser(self.THDumbie.xMin, self.THDumbie.xMax)
        return


    def _CreateCanvasAndLegendAndDumbie(self):
        '''
        Create a TCanvas, a TLegend, and a dubmie TH1.
        '''
        self.Verbose()

        self.CreateDumbie()
        if self.bPadRatio == True or self.bInvPadRatio == True:
            self._Create2PadCanvas()
        else:
            self._CreateCanvas()

        self._CreateLegend()
        return

        
    def CreateCutLines(self):
        '''
        Create TLines for each cut-line defined by the user when creating a histo instance. 
        Append them to the DrawObjectList so that they can be drawn later on.
        '''
        self.Verbose(["Creating cut-lines."])

        ### Create the TLines for each axis
        self._AppendXYCutLinesToTLineList("x")
        self._AppendXYCutLinesToTLineList("y")
        self._DoZCutLines()

        ### Extend the DrawObjectList with the TLineList
        self.DrawObjectList.extend(self.xTLineList)
        self.DrawObjectList.extend(self.yTLineList)

        self.DrawObjectListR.extend( copy.deepcopy(self.xTLineList) )
        if (self.THDumbie.yCutLinesRatioPad == True):
            self.DrawObjectListR.extend( copy.deepcopy(self.yTLineList) )
        return


    def _AppendXYCutLinesToTLineList(self, axis):
        '''
        Create, customise and append x- or y-axis cut lines to the TLineList. Also add entry to TLegend and provide extra space for another TLegened entry.
        '''
        self.Verbose()

        bXaxisCut = False
        bYaxisCut = False

        if axis == "x":
            bXaxisCut = True
            cLines = self.THDumbie.xCutLines
        elif axis == "y":
            bYaxisCut = True
            cLines = self.THDumbie.yCutLines
        else:
            raise Exception("The option 'axis' can either be \"x\" or \"y\". Passed option was \"%s\"." % (axis) )

        if cLines == None:
            return
                
        ### Loop over all x-axis cut values
        for value in cLines:
            if bXaxisCut == True:
                xMin = value
                xMax = value
            else:
                xMin = self.THDumbie.xMin
                xMax = self.THDumbie.xMax

            if bYaxisCut == True:
                yMin = value
                yMax = value
            else:
                yMin = self.THDumbie.yMin
                yMax = self.THDumbie.yMax

            line = ROOT.TLine(xMin, yMin, xMax, yMax)
            self._CustomiseTLine(line, lineColour=self.CutLineColour, lineWidth=3, lineStyle=ROOT.kDashed) #ROOT.kDashDotted
            self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)
            self.AppendToTLineList( line, axis)

        return


    def _DoZCutLines(self):
        '''
        Create, customise and append z-axis cut lines to the TLineList. Also add entry to TLegend and provide extra space for another TLegened entry.
        "CONT"	 Draw a contour plot (same as CONT0).
        "CONT0"	 Draw a contour plot using surface colors to distinguish contours.
        "CONT1"	 Draw a contour plot using line styles to distinguish contours.
        "CONT2"	 Draw a contour plot using the same line style for all contours.
        "CONT3"	 Draw a contour plot using fill area colors.
        "CONT4"	 Draw a contour plot using surface colors (SURF option at theta = 0).
        "CONT5"	 (TGraph2D only) Draw a contour plot using Delaunay triangles.
        See http://root.cern.ch/root/html/THistPainter.html for contour draw options!
        '''
        self.Verbose()

        zCutLines = self.THDumbie.zCutLines
        if zCutLines == None or len(zCutLines)==0:
            return
        
        ### Copy list into a numpy array, otherwise ROOT will complain
        contours = numpy.asarray(self.THDumbie.zCutLines)
        if (self.THDumbie.zCutLinesErrors == True):
            self._DrawZCutLinesWithErrors(zCutLines)
        else:
            self.DrawZCutLines(zCutLines)
        return


    def _DrawZCutLinesWithErrors(self, zCutLines):
        '''
        Draw z-cut lines defined by the user at histogram creation time.
        Add error bands around the nominal contour that correspond to +/- 1 sigma.
        For information regarding the extraction of the contours as TGraphs see:
        http://root.cern.ch/root/html/THistPainter.html
        [ Look for Draw("CONT LIST") option ]
        '''
        ### Copy list into a numpy array, otherwise ROOT will complain
        contours = numpy.asarray(self.THDumbie.zCutLines)

        ### Set Contour bands
        contourBands   = [1.0, 1.0/self.THDumbie.scaleFactor, 1.0/self.THDumbie.scaleFactor]
        lineStyles     = [ROOT.kSolid, ROOT.kDashDotted, ROOT.kDashDotted] #[ROOT.kSolid, ROOT.kSolid, ROOT.kSolid] #,ROOT.kDashed, ROOT.kDashed]
        lineColours    = [self.CutLineColour, ROOT.kRed-4, ROOT.kAzure+6, ROOT.kSpring+2, ROOT.kMagenta-2, ROOT.kYellow-4, self.CutLineColour, ROOT.kOrange+8, 
                          ROOT.kBlue-4, ROOT.kGreen-2, ROOT.kViolet-5, ROOT.kPink-8, ROOT.kTeal-1, ROOT.kCyan-7, 
                          self.CutLineColour, ROOT.kRed-4, ROOT.kAzure+6, ROOT.kSpring+2, ROOT.kMagenta-2, ROOT.kYellow-4, self.CutLineColour, ROOT.kOrange+8,
                          ROOT.kBlue-4, ROOT.kGreen-2, ROOT.kViolet-5, ROOT.kPink-8, ROOT.kTeal-1, ROOT.kCyan-7]
        lineShades     = [+0, -3, -3]
        lineWidths     = [+3, +1, +1]
        
        ### Get Last histogram in the stack and customise it
        hStackList = []
        total = len(zCutLines)*len(contourBands)
        for i in range(0, total ):
            hStackList.append(copy.deepcopy(self.THStack.GetStack()))

        ### Declare the arrays where the contours values (nominal, Up, Down) will be saved
        tmpHistoDict = {}

        ### for-Loop: Over all histogram types (Nominal, FluctuateUp, FluctuateDown)
        for iBand in range(0, len(contourBands)): 

            ### Get the total histogram
            h = hStackList[iBand].Last()

            ### Default, Fluctuate up or Fluctuate down
            if iBand == 0:
                histoType = "Nominal"
            elif iBand == 1: # +1sigma
                self.FluctuateHisto(h, contourBands[iBand], True)
                histoType = "FluctuateUp"
            elif iBand == 2: # -1sigma
                self.FluctuateHisto(h, contourBands[iBand], False)
                histoType = "FluctuateDown"
            else:
                pass

            h.SetName(h.GetName() + "_" + histoType )
            tmpHistoList = []

            ### for-Loop: Over all z-cuts, for the given iBand(0 = Nominal, 1 = +1sigma, 2 = -1sigma)
            for i in range(0, len(zCutLines)):
                self.Verbose( ["iHistoType: '%s'/'%s' (Nominal, FluctuateUp, FluctuateDown)" % (iBand+1, len(contourBands) ), "zCutNumber: '%s'/'%s'" % (i+1, len(zCutLines)) ] )
                h.SetLineColor( lineColours[i] )
                h.SetFillColor( lineColours[i] + lineShades[iBand] )
                h.SetLineStyle( lineStyles[iBand] )
                h.SetLineWidth( lineWidths[iBand] )
                h.SetContour(i+1, contours) #should start at 1

                ### Before drawing, save the histogram to some list
                if self.SaveContourPoints == True:
                    tmpHistoList.append(h)

                ### Draw the histogram
                h.Draw("cont3same9")

                if iBand == 0:
                    self.TLegend.AddEntry( h, "%s %s #pm 1#sigma" % (zCutLines[i], self.THDumbie.zUnits), "L" )
                    self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)

            ### Save the Nominal, FluctuateUp, FluctuateDown histos for the specific cut value
            tmpHistoDict[ histoType ] = tmpHistoList

        self.SaveContourPointsToFile(tmpHistoDict)
        return


    def DrawZCutLines(self, zCutLines):
        '''
        Draw z-cut lines defined by the user at histogram creation time.
        '''
        ### Copy list into a numpy array, otherwise ROOT will complain
        contours = numpy.asarray(self.THDumbie.zCutLines)

        ### Set Contour bands
        lineColours    = [ROOT.kWhite, ROOT.kYellow, ROOT.kRed, ROOT.kBlack]
        
        ### Get Last histogram in the stack and customise it
        hStackList = []
        total = len(zCutLines)
        for i in range(0, total ):
            hStackList.append(copy.deepcopy(self.THStack.GetStack()))
        
        ### Loop over all z-cuts
        for i in range(0, len(zCutLines)):
            ### Get the total histogram
            h = hStackList[i].Last()
            h.SetLineColor( lineColours[i] )
            h.SetFillColor( lineColours[i] )
            h.SetLineStyle( ROOT.kSolid )
            h.SetLineWidth( 3 )
            h.SetContour(i+1, contours) #should start at 1
            h.Draw("cont3same9")
            self.TLegend.AddEntry( h, "%s %s" % (zCutLines[i], self.THDumbie.zUnits), "L" )
            self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)
        return



    def FluctuateHisto(self, histo, scaleFactor, bFluctuateUp):
        '''
        Statistical fluctuation of histograms to enable contour error bands.
        '''
        ### Reverse scaleFactor
        histo.Scale( scaleFactor )
        zeroEntriesPoissonError = 1.7 #As instructed by Ptochos

        if "TH1" in str(type(histo)):
            ### Fluctate with bin content with sqrt(N)
            for x in range(0, histo.GetNbinsX()):
                binContent = histo.GetBinContent(x)
                if (bFluctuateUp == True):
                    if binContent < 1.0:
                        histo.SetBinContent(x, histo.GetBinContent(x) + scaleFactor*zeroEntriesPoissonError )
                    else:
                        histo.SetBinContent(x, histo.GetBinContent(x) + scaleFactor*math.sqrt(histo.GetBinContent(x)) )
                else:
                    if binContent > 0.0:
                        histo.SetBinContent(x, histo.GetBinContent(x ) - scaleFactor*math.sqrt(histo.GetBinContent(x)) )
                    else:
                        histo.SetBinContent(x, binContent) # do nothing
        elif "TH2" in str(type(histo)):
            ### Fluctate with bin content with sqrt(N)
            for x in range(0, histo.GetNbinsX()):
                for y in range(0, histo.GetNbinsY()):
                    binContent = histo.GetBinContent(x, y)
                    if (bFluctuateUp == True):
                        if binContent < 1.0:
                            histo.SetBinContent(x, y, histo.GetBinContent(x, y) + scaleFactor*zeroEntriesPoissonError )
                        else:
                            histo.SetBinContent(x, y, histo.GetBinContent(x, y) + scaleFactor*math.sqrt(histo.GetBinContent(x, y)) )
                    elif (bFluctuateUp == False):
                        if binContent > 0.0:
                            histo.SetBinContent(x, y, histo.GetBinContent(x, y) - scaleFactor*math.sqrt(histo.GetBinContent(x, y)) )
                        else:
                            histo.SetBinContent(x, y, binContent) # do nothing
        else:
            raise Exception("Unsupported histo type '%s'. Cannot proceed." % (type(histo) ) )

        ### Scale back to normal
        histo.Scale( 1.0/scaleFactor )
        return
        

    def CreateCutBoxes(self):
        '''
        Create TBoxes with associated TLines (custom colour) for each list of cut-range defined by the user when creating a histo instance. 
        Append them to the DrawObjectList so that they can be drawn later on.
        '''
        self.Verbose()
        
        ### Loop over list of xMin-xMax-colour pairs (also a list)
        self._AppendXYCutBoxesToTBoxList("x")
        self._AppendXYCutBoxesToTBoxList("y")

        ### Extend the DrawObjectList with the TLineList and TBoxList
        self.DrawObjectList.extend(self.xTLineList)
        self.DrawObjectList.extend(self.yTLineList)
        self.DrawObjectList.extend(self.TBoxList)

        self.DrawObjectListR.extend(copy.deepcopy(self.xTLineList))
        if (self.THDumbie.yCutLinesRatioPad == True):
            self.DrawObjectListR.extend(copy.deepcopy(self.yTLineList))
            self.DrawObjectListR.extend(copy.deepcopy(self.TBoxList))
        return


    def _AppendXYCutBoxesToTBoxList(self, axis):
        '''
        Create, customise and append x- or y-axis cut boxes to the TBoxList. Also add entry to TLegend and provide extra space for another TLegened entry.
        '''
        self.Verbose()


        bXaxisCut = False
        bYaxisCut = False
        if axis == "x":
            bXaxisCut = True
            cBoxes = self.THDumbie.xCutBoxes
        elif axis == "y":
            bYaxisCut = True
            cBoxes = self.THDumbie.yCutBoxes
        else:
            raise Exception("The option 'axis' can either be \"x\" or \"y\". Passed option was \"%s\"." % (axis) )

        if cBoxes == None:
            return

        for v in cBoxes:
            if bXaxisCut == True:
                xMin = v[0]
                xMax = v[1]
            else:
                xMin = self.THDumbie.xMin
                xMax = self.THDumbie.xMax

            if bYaxisCut == True:
                yMin = v[0]
                yMax = v[1]
            else:
                yMin = self.THDumbie.yMin
                yMax = self.THDumbie.yMax

            cutBox  = ROOT.TBox( xMin, yMin, xMax, yMax)
            if bXaxisCut == True:
                cLine1  = ROOT.TLine(xMin, yMin, xMin, yMax)
                cLine2  = ROOT.TLine(xMax, yMin, xMax, yMax)
            else:
                cLine1  = ROOT.TLine(xMin, yMin, xMax, yMin)
                cLine2  = ROOT.TLine(xMin, yMax, xMax, yMax)
            cLine1.SetLineColor(v[2])
            cLine2.SetLineColor(v[2])
            cLine1.SetLineWidth(1)
            cLine2.SetLineWidth(1)
            self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)
            self.AppendToTLineList( cLine1, axis)
            self.AppendToTLineList( cLine2, axis )
            self.AppendToTBoxList(cutBox, boxColour= v[2] )
        return


    def _DrawHistograms(self, THStackDrawOpt):
        '''
        Draw the THDumbie. Draw the THStack . Create the CutBoxes and CutLines. 
        Re-draw the TH1Dubmie to unhide the hidden tickmards (true when drawing histograms 
        will fill style 1001. Draw all objects in the DrawObjectList.
        For drawing options (TH1, TH2, THStack etc..) see:
        http://root.cern.ch/root/html/THistPainter.html
        '''
        self.Verbose()

        self.THDumbie.TH1orTH2.Draw(self.THDumbie.drawOptions)
        self.DrawStackInclusive()
        self.THStack.Draw(THStackDrawOpt + "," + self.THDumbie.drawOptions + "," +  "9same") #"PADS"
        ROOT.gPad.RedrawAxis() #the histo fill area may hide the axis tick marks. Force a redraw of the axis over all the histograms.
        self.AddDefaultCmsText()
        self.TCanvas.Update()
        self.TCanvas.SetGridx(self.THDumbie.gridX)
        self.TCanvas.SetGridy(self.THDumbie.gridY)
        return


    def _DrawRatioHistograms(self, bAddReferenceHisto=True):
        '''
        Draw all plots on the PadRatio (if applicable).
        For efficiencies and errors see:
        http://steve.cooleysekula.net/goingupalleys/2011/08/09/python-and-root-tricks-efficiency-graphs/
        '''
        self.Verbose()

        if self.bPadRatio == False and self.bInvPadRatio == False:
            return

        self.PadRatio.cd()        

        ### Create the histogram that will divide all other histograms in the THStackRatio (Normalisation histogram)
        UnityTH1 = self.GetUnityTH1()
        hDenominator = copy.deepcopy( self.THStackHistoList[0] ) 
        self.Verbose(["Using histogram '%s' as denominator for ratio plots! " % (hDenominator.GetName())])

        ### Add the reference ratio plot (to enable the identification of the histogram used for the normalisation)
        ### Note: Do not add the hReference histogram  before calling the function self.CustomiseTHRatio(). 
        if bAddReferenceHisto:
            hReference = copy.deepcopy(hDenominator)
            hReference.Divide(hReference)
            hReference.SetMarkerSize(0)
            hReference.SetLineStyle(ROOT.kSolid)            
            for iBin in range (1, hReference.GetNbinsX()+1):
                hReference.SetBinError(iBin, 0.00000000001)
            self.THStackRatio.Add(hReference)

        ### Loop over all histograms in the THStack and create the THStackRatio.
        for h in self.THStackHistoList:
            if h == self.THStackHistoList[0]:                
                continue

            ### Copy Active histogram
            hRatio     = copy.deepcopy(h)
            hNumerator = copy.deepcopy(h)
            
            ### Divide the active histogram with the normalisation histogram
            hRatio.Divide(hNumerator, hDenominator, 1.0, 1.0, self.RatioErrorType)

            ### Inverts ratio histogram if requested (i.e. each bin has content 1/bin)
            if self.bInvPadRatio == True:
                hRatio.Divide(UnityTH1, hRatio) 

            ### Save histogram values to a txt file (for later processing if needed)
            #self.SaveHistoAsTxtFile(hRatio)

            ### Finally, add this ratio histogram to the THStackRatio
            self.THStackRatio.Add(hRatio)

        ### Customise axes and titles
        self.CustomiseTHRatio()

        ### Draw the Ratio Stack with "nostack" option
        self.THRatio.TH1orTH2.Draw()
        self.THStackRatio.Draw("nostack9sameAP")
        
#        ### Finally add the reference ratio plot (to enable the identification of the histogram used for the normalisation)
#        ### Note: Do not add the hReference histogram  before calling the function self.CustomiseTHRatio(). 
#        if bAddReferenceHisto:
#            hReference = copy.deepcopy(hDenominator)
#            hReference.Divide(hReference)
#            hReference.SetMarkerSize(0)
#            hReference.SetLineStyle(ROOT.kSolid)            
#            for iBin in range (1, hReference.GetNbinsX()+1):
#                hReference.SetBinError(iBin, 0.00000000001)
#            self.THStackRatio.Add(hReference)

        ### Switch back to the PadPlot (necessary)
        self.PadPlot.cd()
        return

    
    def SaveHistoAsTxtFile(self, h):
        '''
        '''
        self.Verbose()
        print "Fix me. Saving as .txt not needed anymore. EXIT"
        sys.exit()

        xVals     = []
        xErrDown  = []
        xErrUp    = []
        yVals     = []
        yErrDown  = []
        yErrUp    = []

        ### Loop over all histogram bins (except under/over bins)
        binZeroWidth = h.GetXaxis().GetBinWidth(0)
        for iBin in range (1, h.GetNbinsX()+1):

            xVals.append( h.GetBinCenter(iBin) ) #h.GetBinCenter(iBin)
            xErrUp.append( binZeroWidth/2.0 )
            xErrDown.append( binZeroWidth/2.0 )

            yVals.append( h.GetBinContent(iBin) )
            yErrUp.append( h.GetBinErrorUp(iBin) )
            yErrDown.append( h.GetBinErrorLow(iBin) )
        self.SaveTGraphArraysDataToFile( h.GetName() + ".txt", xVals, xErrUp, xErrDown, yVals, yErrUp, yErrDown, 3)
        return



    def CustomiseTHRatio(self):
        '''
        Apply all necessary customisations to self.THRatio histogram.
        '''

        ### Customise axes and titles
        if self.THRatio.yMinRatio == None:
            self.THRatio.yMinRatio = self.THStackRatio.GetMinimum("nostack")*self.THRatio.GetYMinFactor(self.THRatio.logYRatio)
        if self.THRatio.yMaxRatio == None:
            self.THRatio.yMaxRatio = self.THStackRatio.GetMaximum("nostack")*self.THRatio.GetYMaxFactor(self.THRatio.logYRatio)

        ### Customise the title
        self.THRatio.TH1orTH2.GetXaxis().SetTitleOffset(2.8)
        if self.THRatio.ratioLabel == None:
            if self.bInvPadRatio == False:
                self.THRatio.ratioLabel = "Ratio"
            else:
                self.THRatio.ratioLabel = "1/Ratio"
        self.THRatio.TH1orTH2.GetYaxis().SetTitle(self.THRatio.ratioLabel)

        ### Customise the y-axis
        self.THRatio.yMax = self.THRatio.yMinRatio
        self.THRatio.yMin = self.THRatio.yMaxRatio
        self.THRatio.TH1orTH2.GetYaxis().SetNdivisions(505)
        self.THRatio.TH1orTH2.GetYaxis().SetRangeUser(self.THRatio.yMinRatio, self.THRatio.yMaxRatio)
        self.THRatio.TH1orTH2.GetYaxis().SetTitleOffset(1.8) 
        self.THDumbie.TH1orTH2.GetYaxis().SetTitleOffset(1.8)
        
        ### Enable grid to easy readout of histo
        self.PadRatio.SetGridx(self.THDumbie.gridX)
        self.PadRatio.SetGridy(self.THDumbie.gridY)
        return


    def DrawStackInclusive(self):
        '''
        The GetStack function returns a TObjArray* of TH1* where the TH1 at index i is the sum of histograms 0->i.
        TObjArray::Last() returns the last TH1 in the list, hence the sum of all TH1.
        For help see: http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=12138
        '''
        self.Verbose()
        if self.bStackInclusive==False:
            return

        inclusive = self.THStack.GetStack().Last()
        inclusive.SetFillStyle(0)
        inclusive.SetFillColor(ROOT.kGray)
        inclusive.SetLineColor(ROOT.kGray)
        inclusive.SetLineStyle(ROOT.kSolid)
        inclusive.SetLineWidth(3)

        if self.THDumbie.yMax < inclusive.GetMaximum():
            yMaxNew = inclusive.GetMaximum()
            h       = self.THDumbie
            h.yMax  = yMaxNew*h.GetYMaxFactor(self.THDumbie.logY)
            h.TH1orTH2.GetYaxis().SetRangeUser(h.yMin, h.yMax)        
        else:
            pass
        inclusive.Draw("HISTsame9")
        
        ### Add histogram entry to the legend
        self.TLegend.AddEntry( inclusive, "inclusive", "L" )
        self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)

        return

    

    def GetNumberOfHistosPerDataset(self):
        '''
        '''
        self.Verbose()
        
        numerator    = len(self.HistoObjectList)
        denominator  = len(self.DatasetToHistoMap.keys())
        return float(numerator)/float(denominator)


    def GetNumberOfDatasets(self):
        '''
        '''
        self.Verbose()
        
        return len(self.DatasetToHistoMap.keys())


    def _GetLegEntryLabel(self, histoObject):
        '''
        Determine and return the TLegenend entry label for this instance of histogram object.
        '''
        self.Verbose()
        
        entryLabel = "empty"
        if histoObject.legTitle == None:            
            return  entryLabel

        if self.UseDatasetAsLegEntry==True:
            entryLabel   = self.DatasetToLatexNameMap[histoObject.dataset]
        else:
            entryLabel = histoObject.legTitle

        if isinstance(entryLabel, str) == True:
            self.Verbose(["The TLegend entry name is '%s'." % (entryLabel)])
            return entryLabel
        else:
            raise Exception("The TLegend entry label cannot be returned as it is not of type string but instead of type '%s'." % ( type(entryLabel) ) )


    def _GetLegEntryOptions(self, histoObject):
        '''
        Determine the draw options for all histograms in the THStack by examining the TLegend entry styles:
        "L": draw line associated with TAttLine if obj inherits from TAttLine
        "P": draw polymarker associated with TAttMarker if obj inherits from TAttMarker
        "F": draw a box with fill associated wit TAttFill if obj inherits TAttFill
        "E": draw vertical error bar
        '''
        self.Verbose()
        
        options = histoObject.legOptions
        #if self.bPadRatio == False and  self.bInvPadRatio == False:
        #    options = histoObject.legOptions
        #else:
        #    options = histoObject.legOptions + "P"
        return options


    def AppendToTLineList(self, line, axis):
        '''
        Append a TLine to the TLineList. 
        '''
        self.Verbose()

        if axis == "x":
            self.xTLineList.append(line)
        elif axis == "y":
            self.yTLineList.append(line)
        else:
            raise Exception("The option 'axis' can either be \"x\" or \"y\". Passed option was \"%s\"." % (axis) )

        return


    def _CustomiseTLine(self, line, lineColour=ROOT.kBlack, lineWidth=3, lineStyle=ROOT.kSolid):
        '''
        '''
        self.Verbose()

        line.SetLineWidth(lineWidth)
        line.SetLineStyle(lineStyle)
        line.SetLineColor(lineColour)
        line.Draw()
        return


    def AppendToTBoxList(self, box, boxColour=18):
        '''
        Append a TBox to the TBoxList. 
        '''
        self.Verbose()
        box.SetFillStyle(3003) #3003
        box.SetFillColor(boxColour)
        self.TBoxList.append(box)
        return


    def _DrawNonHistoObjects(self):
        '''
        Draw all drawable objects found in self.DrawObjectList.
        '''
        self.Verbose()                

        if (self.bPadRatio == True or self.bInvPadRatio == True):
            self._DrawNonHistoObjectsWithPadRatio()
        else:
            self._DrawNonHistoObjectsNoPadRatio()        
        return


    def _DrawNonHistoObjectsNoPadRatio(self):
        '''
        Draw all drawable objects found in self.DrawObjectList.
        '''
        self.Verbose()                

        ### First create the draw objects (TLines, TBoxes etc..)
        self.CreateCutBoxes()
        self.CreateCutLines()
        
        ### Draw all objects on the PadPlot
        for o in self.DrawObjectList:
            o.Draw("same")
        return


    def _DrawNonHistoObjectsWithPadRatio(self):
        '''
        Draw all drawable objects found in self.DrawObjectList.
        '''
        self.Verbose()                

        ### First create the draw objects (TLines, TBoxes etc..)
        self.CreateCutBoxes()
        self.CreateCutLines()
        
        ### Draw all objects on the PadPlot
        self.PadPlot.cd()
        for o in self.DrawObjectList:
            o.Draw("same")

        ### Update modified canvas and re-draw the PadPlot axes
        self.PadPlot.Modified()
        self.PadPlot.RedrawAxis()
        self.PadPlot.SetGridx(self.THDumbie.gridX)
        self.PadPlot.SetGridy(self.THDumbie.gridY)

        self.PadRatio.Modified()
        self.PadRatio.RedrawAxis()


        ### Draw all objects on the PadRatio
        self.PadRatio.cd()
        for o in self.DrawObjectListR:
            if isinstance(o, ROOT.TLegend) == True:
                continue

            if ( ( o.GetY1() != o.GetY2() ) and (o.GetX1() == self.THRatio.xMin ) and ( self.THRatio.xMax == o.GetX2() ) ):
                continue
            elif ( (o.GetX1() == self.THRatio.xMin ) and ( self.THRatio.xMax == o.GetX2() ) ):
                continue               
            elif( o.GetX1() == o.GetX2() ):
                o.SetY1(self.THDumbie.yMinRatio)
                o.SetY2(self.THDumbie.yMaxRatio)
            else:
                #print "o.GetX1() = %s, o.GetX2() = %s, o.GetY1() = %s, o.GetY2() = %s" % ( o.GetX1(), o.GetX2(),o.GetY1(), o.GetY2())
                o.SetY1(self.THDumbie.yMinRatio)
                o.SetY2(self.THDumbie.yMaxRatio)
            o.Draw("same")

        ### Update modified canvas and re-draw the PadPlot axes
        self.PadRatio.Modified()
        self.PadRatio.RedrawAxis()
        
        self.PadPlot.cd()
        return

    
    def _CheckHistoBinning(self, histoObject):
        '''
        Ensure that the histoObject has exactly the same binning as the TH1Dubmie.
        '''
        self.Verbose()

        binningIsOk       = False
        binWidthX         = self.THDumbie.binWidthX
        binZeroWidth      = self.THDumbie.TH1orTH2.GetXaxis().GetBinWidth(0)
        tmpBinWidthX      = histoObject.binWidthX
        tmpBinZeroWidth   = histoObject.TH1orTH2.GetXaxis().GetBinWidth(0)
        if (tmpBinWidthX != binWidthX or tmpBinZeroWidth!=binZeroWidth):
            raise Exception("At least one of the histogram in the plotting queue has a different x-axis binning! Please make sure all your histogram bins are identical.")
        return 


    def _CheckHistogramBinning(self):
        '''
        Ensure that all histoObjects have exactly the same binning as the TH1Dubmie.
        '''
        self.Verbose()

        binningIsOk  = False
        binWidthX    = self.THDumbie.binWidthX
        binZeroWidth = self.THDumbie.TH1orTH2.GetXaxis().GetBinWidth(0)
        
        ### Check that the x-axis bin is same for all histograms
        for dataset in self.DatasetToHistoMap.keys():
            h = self.DatasetToHistoMap[dataset]
            self._CheckHistoBinning(h)
        return 


    def AddDefaultCmsText(self):
        '''
        Add the default CMS text on the canvas. Several defaults are available. 
        For available options see the class TextClass(object) under tools/text.py.
        '''
        self.Verbose()
        
        self.TCanvas.cd()
        self.TextObject.AddEnergyText( text = " " )
        self.TextObject.AddCmsSimulationPhaseTwoText(self.IsTH2, self.PileUp) 
        self.TextObject.AddIntLumiText( text = "" )

        self.TCanvas.Update()
        return
    

    def SaveHistos(self, bSave=False, savePath=os.getcwd() + "/", saveFormats=[".png", ".C", ".eps", ".pdf"], saveExtension=""):
        '''
        Save all canvases to specified path and with the desirable format.
        '''
        self.Verbose()
        
        if bSave == False:
            return
            
        ### Sanity checks
        if savePath == "" or savePath == None:
            savePath = os.getcwd() + "/"            
        if os.path.exists(savePath) == False: 
            raise Exception("The path '%s' does not exist! Please make sure the provided path is correct." % (savePath) )

        ### Define path and save
        saveName = savePath + self.TCanvas.GetName().rsplit('@', 1)[0] + saveExtension
        #creationTime =  self.TCanvas.GetName().rsplit('@', 1)[1]
                                                              

        self.Print(["SaveName: '%s' " % (saveName), "Format(s): '%s' " % ( ", ".join(saveFormats) ) ])
        for ext in saveFormats:
            self.TCanvas.Update()
            self.TCanvas.SaveAs( saveName + "." + ext )

        #ROOT.gDirectory.ls()
        return


    def GetCanvasName(self):
        '''
        Return canvas name upon which the plotter objects histograms will be saved on.
        '''
        self.Verbose()
        
        return self.TCanvas.GetName()


    def SetCanvasName(self, canvasName):
        '''
        '''
        self.Verbose()
        
        self.TCanvas.SetName(canvasName)
        return


    def AppendToCanvasName(self, canvasNameExt):
        '''
        '''
        self.Verbose()
        
        self.TCanvas.SetName( self.TCanvas.GetName() + canvasNameExt)
        return


    def GetUnityTH1(self):
        '''
        Returns a TH1 with identical attributes to those of self.TH1Dubmie. But, all its bins
        are filled with 1. So you have a flat distribution histogram at y=1, over the entire x-axis range.
        '''
        self.Verbose()
        hUnity = copy.deepcopy(self.THRatio)
        hUnity.TH1orTH2.Reset()
        hUnity.TH1orTH2.SetName("hUnity")
        hUnity.name = "hUnity"

        nBins = hUnity.TH1orTH2.GetNbinsX()
        for i in range (0, nBins+1):
            #hUnity.TH1orTH2.Fill(i, 1) # error bars are quite large. need to investigate if they are correct.
            ### In the meantime, set the bin-error to zero. I think this is the correct way to do it
            hUnity.TH1orTH2.SetBinContent(i, 1)
            hUnity.TH1orTH2.SetBinError(i, 0)
            
        return hUnity.TH1orTH2


    def PrintPSet(self, psetFolderPath, rootFileName = None):
        ''' 
        Print the PSet used when generating a particular ROOT file.
        '''

        self.Verbose()

        ### Loop over all datasets
        if rootFileName == None:

            if len(self.DatasetToRootFileMap.keys()) < 1:
                raise Exception("Cannot find PSets as no datasets are available. Either add a dataset or provide explicitly the ROOT file you want to read the PSets from." )
            else:
                for dataset in self.DatasetToRootFileMap.keys():
                    f     = self.DatasetToRootFileMap[dataset]
                    named = f.Get(psetFolderPath)
                    rootFileName = f.GetName()
                    break
        else:
            named = self.GetRootFile(rootFileName).Get(psetFolderPath)

        self.Print(["Printing PSets for ROOT file '%s':" % (rootFileName), named.GetTitle() ])
        return


    def CustomiseTGraph(self, graph, histoObject, dataset, kwargs):
        '''
        '''
        self.Verbose()

        graph.SetName( histoObject.name + "_TGraph")

        ### Once drawn, the TGraph can be customised
        if kwargs.get("xMin") != None and kwargs.get("xMax") != None:
            graph.GetXaxis().SetRangeUser( kwargs.get("xMin"), kwargs.get("xMax") ) 
        if kwargs.get("yMin") != None and kwargs.get("yMax") != None:
            graph.GetYaxis().SetRangeUser( kwargs.get("yMin"), kwargs.get("yMax") ) 
            graph.GetYaxis().SetLimits( kwargs.get("yMin"), kwargs.get("yMax") )         
            
        ### Customise TGraph (colours/styles)
        if dataset == None:
            raise Exception("Could not determine dataset('%s'). Cannot continue." % (dataset) )
        elif ":" in dataset:
            styleType = dataset.rsplit(":", 1)[0].lower()    
        else:
            styleType = dataset

        (fillColour, lineColour, markerStyle, markerSize, lineWidth, lineStyle, fillStyle, drawOptions, legOptions) = self.StyleObject.GetTGraphStyles( styleType )
        graph.SetFillColor(fillColour)
        graph.SetFillStyle(fillStyle)
        graph.SetLineColor(lineColour)
        graph.SetLineStyle(lineStyle)
        graph.SetLineWidth(lineWidth)
        graph.SetMarkerColor(fillColour)
        graph.SetMarkerStyle(markerStyle)
        graph.SetMarkerSize(markerSize)
        
        ### Create legend if there isn't one already
        if self.TLegend == None:
            if ":" in dataset:
                datasetLatex = self.DatasetObject.ConvertDatasetToLatex( dataset.rsplit(":", 1)[0] )
            else:
                datasetLatex = self.DatasetObject.ConvertDatasetToLatex( dataset )
            self.TLegend = ROOT.TLegend( kwargs.get("xLegMin"), kwargs.get("yLegMin"), kwargs.get("xLegMax"), kwargs.get("yLegMax"), datasetLatex, "brNDC")
            self._CustomiseLegend()
            
        ### Add TGraph entry to legend
        legOpts = kwargs.get("legOptions") 
        if legOpts == None:
            legOpts = "LFP"

        self.TLegend.AddEntry( graph, self._GetLegEntryLabel(histoObject), legOpts)
        self.TLegend.SetY1( self.TLegend.GetY1() - 0.02)
        return        
        

    def DrawMultigraph(self, saveName, **kwargs):
        '''
        Draws a ROOT.TMultiGraph object that has been created before-hand, using the options provided with the keywords arguments.
        '''
        self.Verbose()
        
        ### Create THDumbie (needed to apply kwargs options)
        self.THDumbie = m_histos.TH1orTH2( path="", name="MGraphDumbie", legTitle="", saveName=saveName, **kwargs)
        
        ### Create a canvas
        self._CreateCanvas()

        ### Draw the multigraph and customise it 
        self.Verbose(["Drawing ROOT.TMultiGraph '%s'" % self.TMultigraph.GetName() ])
        self.TMultigraph.Draw(self.THDumbie.drawOptions)
        self.CustomiseMultigraph(kwargs)
        self._DrawNonHistoObjectsNoPadRatio()
        return


    def CustomiseMultigraph(self, kwargs):
        '''
        '''
        self.Verbose()

        ### Customise x-axis
        xLabel = ""
        xUnits = kwargs.get("xUnits", "")
        xMin   = kwargs.get("xMin", None)
        xMax   = kwargs.get("xMax", None)

        if xUnits == "":
            xLabel = kwargs.get("xLabel", "x-label")
        else:
            xLabel = kwargs.get("xLabel", "x-label") + " (" + xUnits + ")"

        ### Customise y-axis
        yLabel = ""
        yUnits = kwargs.get("yUnits", "")
        yMin   = kwargs.get("yMin", None)
        yMax   = kwargs.get("yMax", None)
        if yUnits == "":
            yLabel = kwargs.get("yLabel", "y-label")
        else:
            yLabel = kwargs.get("yLabel", "y-label") + " (" + yUnits + ")"
            
        ### Apply options
        self.TMultigraph.GetXaxis().SetTitle( xLabel )
        self.TMultigraph.GetYaxis().SetTitle( yLabel )
        self.TMultigraph.GetXaxis().SetLimits(xMin, xMax)
        self.TMultigraph.GetYaxis().SetLimits(yMin, yMax)
        self.TMultigraph.GetYaxis().SetRangeUser( yMin, yMax) 

        ### Customise/Draw TLegend
        self.TLegend.SetX1( kwargs.get("xLegMin") )
        self.TLegend.SetX2( kwargs.get("xLegMax") )
        self.TLegend.SetY1( kwargs.get("yLegMin") )
        self.TLegend.SetY2( kwargs.get("yLegMax") )
        self.TLegend.Draw()

        ### Default Text
        self.AddDefaultCmsText()
        self.TCanvas.Update()

        ### Logarithmic axes
        self.TCanvas.SetLogx( kwargs.get("logX") )
        self.TCanvas.SetLogy( kwargs.get("logY") )
        self.TCanvas.Update()

        ### Grids
        self.TCanvas.SetGridx( kwargs.get("gridX") )
        self.TCanvas.SetGridy( kwargs.get("gridY") )
        self.TCanvas.Update()
        return        
        

    def SaveTGraphArraysDataToFile(self, fileName, xVals, xErrUp, xErrDown, yVals, yErrUp, yErrDown, maxDecimals = 1, bWriteHeader=True):
        '''
        Reading and Writing Files in Python. The modes can be:
        
        'r' when the file will only be read
        
        'w' for only writing (an existing file with the same name will be erased)
        
        'a' opens the file for appending; any data written to the file is automatically
        added to the end. 
        
        'r+' opens the file for both reading and writing.
        '''
        self.Verbose()

        ### Loop over all lines in myFile. Append info to file
        self.Verbose( ["SaveName: '%s' " % (fileName), "MaxDecimals: '%s'" % (maxDecimals)])
        
        ### Open file in append ("a") mode
        f = open(fileName, "a")

        ### Define the column widths: 10 each
        template = "{0:10} : {1:10} : {2:10} : {3:10} : {4:10} : {5:10}"

        ### Create the table headers
        header = template.format("x", "xErrUp", "xErrDown", "y", "yErrUp", "yErrDown")
        if bWriteHeader == True:
            f.write( header + "\n" )
     
        ### Now write the actual values
        for x, xUp, xDown, y, yUp, yDown in  zip(xVals, xErrUp, xErrDown, yVals, yErrUp, yErrDown):
            column1  = str( round( x    , maxDecimals) )
            column2  = str( round( xUp  , maxDecimals) ) 
            column3  = str( round( xDown, maxDecimals) )
        
            column4  = str( round( y    , maxDecimals) )
            column5  = str( round( yUp  , maxDecimals) ) 
            column6  = str( round( yDown, maxDecimals) )

            row = column1 + column2 + column3 + column4 + column5 + column6

            line = template.format(column1, column2, column3, column4, column5, column6) + "\n"
            self.Verbose( ["Writing line: %s " % (line) ] )
            f.write( line  )

        f.close()
        return



    def SaveTGraphArraysDataToFileWithZ(self, fileName, xVals, xErrUp, xErrDown, yVals, yErrUp, yErrDown, zVals, zErrUp, zErrDown, maxDecimals = 1, bWriteHeader=True):
        '''
        Reading and Writing Files in Python. The modes can be:
        
        'r' when the file will only be read
        
        'w' for only writing (an existing file with the same name will be erased)
        
        'a' opens the file for appending; any data written to the file is automatically
        added to the end. 
        
        'r+' opens the file for both reading and writing.
        '''
        self.Verbose()

        ### Sanity check: All value lists have the same number of entries
        length = len(xVals)
        if all(len(lst) != length for lst in [xErrUp, xErrDown, yVals, yErrDown, zVals, zErrUp, zErrDown]):
            self.Print( ["len(xVals) = '%s'"    % (len(xVals)), "len(xErrUp) = '%s'"   % (len(xErrUp)), "len(xErrDown) = '%s'" % (len(xErrDown)), 
                         "len(yVals) = '%s'"    % (len(yVals)), "len(yErrUp) = '%s'"   % (len(yErrUp)), "len(yErrDown) = '%s'" % (len(yErrDown)), 
                         "len(zVals) = '%s'"    % (len(zVals)), "len(zErrUp) = '%s'"   % (len(zErrUp)), "len(zErrDown) = '%s'" % (len(zErrDown)), "EXIT"] )

        ### Loop over all lines in myFile. Append info to file
        self.Verbose( ["SaveName: '%s' " % (fileName), "MaxDecimals: '%s'" % (maxDecimals)])
        
        ### Open file in append ("a") mode
        f = open(fileName, "a")

        ### Define the column widths: 10 each
        template = "{0:10} : {1:10} : {2:10} : {3:10} : {4:10} : {5:10} : {6:10} : {7:10} : {8:10}"

        ### Create the table headers
        header = template.format("x", "xErrUp", "xErrDown", "y", "yErrUp", "yErrDown", "z", "zErrUp", "zErrDown")
        if bWriteHeader == True:
            f.write( header + "\n" )

        self.Verbose( ["(x , y , z) = (%s , %s , %s)" % (xVals , yVals , zVals)]  )
        self.Verbose( ["(x+ , y+ , z+) = (%s , %s , %s)" % (xErrUp , yErrUp , zErrUp)]  )
        self.Verbose( ["(x- , y- , z-) = (%s , %s , %s)" % (xErrDown , yErrDown , zErrDown)]  )
        ### Now write the actual values
        for x, xUp, xDown, y, yUp, yDown, z, zUp, zDown in  zip(xVals, xErrUp, xErrDown, yVals, yErrUp, yErrDown, zVals, zErrUp, zErrDown):
            column1  = str( round( x    , maxDecimals) )
            column2  = str( round( xUp  , maxDecimals) ) 
            column3  = str( round( xDown, maxDecimals) )
        
            column4  = str( round( y    , maxDecimals) )
            column5  = str( round( yUp  , maxDecimals) ) 
            column6  = str( round( yDown, maxDecimals) )

            column7  = str( round( z    , maxDecimals) )
            column8  = str( round( zUp  , maxDecimals) ) 
            column9  = str( round( zDown, maxDecimals) )

            row = column1 + column2 + column3 + column4 + column5 + column6 + column7 + column8 + column9

            line = template.format(column1, column2, column3, column4, column5, column6, column7, column8, column9) + "\n"
            self.Verbose( ["Writing line: %s " % (line) ] )
            f.write( line  )

        f.close()
        return


    def GetTGraphArraysDataFromFile(self, inputFilePath):
        '''
        '''
        self.Verbose()

        if os.path.exists(inputFilePath) == False:
            raise Exception("File '%s' does not exist! Please make sure the provided path for the file is correct." % (inputFilePath) )
        else:
            pass

        ### Variable declaration            
        xVals    = []
        xErrUp   = []
        xErrDown = []
        yVals    = []
        yErrUp   = []
        yErrDown = []

        ### Loop over all lines in myFile
        with open(inputFilePath) as f:
            ### Skip the header row and start reading a file from line2
            next(f)
            ### Loop over all lines
            for line in f:
                
                x     = ( line.rsplit(":", 5)[0] ).replace(" ", "")
                xUp   = ( line.rsplit(":", 5)[1] ).replace(" ", "")
                xDown = ( line.rsplit(":", 5)[2] ).replace(" ", "")
                y     = ( line.rsplit(":", 5)[3] ).replace(" ", "")
                yUp   = ( line.rsplit(":", 5)[4] ).replace(" ", "")
                yDown = ( line.rsplit(":", 5)[5] ).replace(" ", "").replace("\n", "")

                xVals   .append( x     )
                xErrUp  .append( xUp   )
                xErrDown.append( xDown )
                yVals   .append( y     )
                yErrUp  .append( yUp   )
                yErrDown.append( yDown )

        return xVals, xErrUp, xErrDown, yVals, yErrUp, yErrDown


    def GetTGraphArraysDataFromFileWithZ(self, inputFilePath):
        '''
       '''
        self.Verbose()

        if os.path.exists(inputFilePath) == False:
            raise Exception("File '%s' does not exist!" % (inputFilePath) )

        ### Variable declaration            
        xVals    = []
        xErrUp   = []
        xErrDown = []

        yVals    = []
        yErrUp   = []
        yErrDown = []

        zVals    = []
        zErrUp   = []
        zErrDown = []
        
        nColumns = 8

        ### Loop over all lines in myFile
        with open(inputFilePath) as f:
            
            ### Skip the header row and start reading a file from line2
            #next(f)
            
            ### Loop over all lines
            for line in f:

                ### Skip headers column titles (quick and dirty fix)
                if "x" in line:
                    continue

                x     = ( line.rsplit(":", nColumns)[0] ).replace(" ", "")
                xUp   = ( line.rsplit(":", nColumns)[1] ).replace(" ", "")
                xDown = ( line.rsplit(":", nColumns)[2] ).replace(" ", "")
         
                y     = ( line.rsplit(":", nColumns)[3] ).replace(" ", "")
                yUp   = ( line.rsplit(":", nColumns)[4] ).replace(" ", "")
                yDown = ( line.rsplit(":", nColumns)[5] ).replace(" ", "")

                z     = ( line.rsplit(":", nColumns)[6] ).replace(" ", "")
                zUp   = ( line.rsplit(":", nColumns)[7] ).replace(" ", "")
                zDown = ( line.rsplit(":", nColumns)[8] ).replace(" ", "").replace("\n", "")

                xVals   .append( x     )
                xErrUp  .append( xUp   )
                xErrDown.append( xDown )

                yVals   .append( y     )
                yErrUp  .append( yUp   )
                yErrDown.append( yDown )

                zVals   .append( z     )
                zErrUp  .append( zUp   )
                zErrDown.append( zDown )

        return xVals, xErrUp, xErrDown, yVals, yErrUp, yErrDown, zVals, zErrUp, zErrDown


    def AddTGraphErrors(self, histoObject, x, xUp, xLow, y, yUp, yLow, z=None, zUp=None, zLow=None, bDrawZValues=False, **kwargs):
        '''
        '''
        self.Verbose()
        self.PrintHistoInfo(histoObject, False)

        ### Convert array items from string to float!
        x     = [ float(item) for item in x     ]
        xUp   = [ float(item) for item in xUp   ]
        xLow  = [ float(item) for item in xLow ]

        y     = [ float(item) for item in y     ]
        yUp   = [ float(item) for item in yUp   ]
        yLow  = [ float(item) for item in yLow ]

        ### Create the numpy arrays as input to TGraph
        xVals    = numpy.asarray( x    )
        xErrsUp  = numpy.asarray( xUp  )
        xErrsLow = numpy.asarray( xLow )

        yVals    = numpy.asarray( y    )
        yErrsUp  = numpy.asarray( yUp  )
        yErrsLow = numpy.asarray( yLow )
        
        if z!=None:
            zVals    = numpy.asarray( z    )
            zErrsUp  = numpy.asarray( zUp  )
            zErrsLow = numpy.asarray( zLow )

            ### Add z-value as text  at each (x,y) point
            if (bDrawZValues):
                for i in range(0, len(x) ):
                    zName = str(int( z[i]) )
                    xPos  = x[i]
                    yPos  = y[i]
                    if i != 48: # with current settings this means 50kHz (older values: 12)
                        continue
                    self.TextObject.AddText(xPos, yPos, zName )
                    l = ROOT.TLatex( xPos, yPos, zName + " GeV")
                    l.SetTextSize(16)  #10
                    l.SetTextAngle(45) #60
                    self.DrawObjectList.append(l)
            else:
                pass

        ### Create TGraphAsymmErrors
        histoObject.dataset = histoObject.dataset
        
        ### Sanity checks
        if (len(xVals) == 0):
            self.Print(["WARNING! Arrays have zero length! Skipping TGraph of histoObject '%s'" % (histoObject.name)])
            return
        elif (len(xVals) == len(yVals) == len(xErrsUp) == len(xErrsLow) == len(yErrsUp) == len(yErrsLow) ):
            pass
        else:
            self.Print(["ERROR! Arrays have different length! EXIT"])
            sys.exit()
        g = ROOT.TGraphAsymmErrors( len(xVals), xVals, yVals, xErrsUp, xErrsLow, yErrsUp, yErrsLow)

        ### Draw graph to enable its customisation
        g.Draw()

        ### Once drawn, customise the TGraph 
        datasetName, datasetLatexName = self.GetDatasetAndLatexName(histoObject.dataset)
        histoObject.datasetName = datasetLatexName
        self.CustomiseTGraph(g, histoObject, datasetLatexName, kwargs)

        ### Finally, add the TGraph to the TMultiGraph
        self.Verbose(["Adding TGraphAsymmErrors '%s' to TMultigraph." % ( g.GetName() )])
        self.TMultigraph.Add(g)

        return


    def SetBoolUseDatasetAsLegEntry(self, myBool):
        '''
        '''
        self.Verbose()
        self.UseDatasetAsLegEntry = myBool
        return


    def AddTGraphErrorsFromFile(self, histoObject, x, xUp, xDown, y, yUp, yDown, **kwargs):
        '''
        '''
        self.Verbose()

        ### Convert array items from string to float!
        x     = [ float(item) for item in x     ]
        xUp   = [ float(item) for item in xUp   ]
        xDown = [ float(item) for item in xDown ]
        y     = [ float(item) for item in y     ]
        yUp   = [ float(item) for item in yUp   ]
        yDown = [ float(item) for item in yDown ]

        ### Create the numpy arrays as input to TGraph
        xVals     = numpy.asarray( x     )
        xErrsUp   = numpy.asarray( xUp   )
        xErrsDown = numpy.asarray( xDown )

        yVals     = numpy.asarray( y     )
        yErrsUp   = numpy.asarray( yUp   )
        yErrsDown = numpy.asarray( yDown )
        
        ### Create TGraphAsymmErrors
        histoObject.dataset = kwargs.get("dataset")
        g = ROOT.TGraphAsymmErrors( len(xVals), xVals, yVals, xErrsUp, xErrsDown, yErrsUp, yErrsDown)

        ### Draw graph to enable its customisation
        g.Draw()

        ### Once drawn, customise the TGraph 
        datasetName, datasetLatexName = self.GetDatasetAndLatexName(histoObject.dataset)
        histoObject.datasetName = datasetLatexName        
        self.CustomiseTGraph(g, histoObject, datasetLatexName, histoObject.kwargs)

        ### Finally, add the TGraph to the TMultiGraph
        self.TMultigraph.Add(g)
        return


    def SaveTGraphErrorsDataToFile(self, fileName, dataset, effVals, effErrUp, effErrDown, effEtThresholds, rateVals, rateValsUp, rateValsDown):
        '''
        Save all the x-axis and y-axis values (and their corresponding errors) to a txt file. 
        The txt file's name is given as the parameter 'fileName' and also contains the dataset information.
        '''
        self.Verbose()

        ### Declarations
        datasetName = dataset
        maxDecimals = 1
        fileName    = fileName + "_"  + datasetName.replace(":", "_") + ".txt"
        f           = open( fileName, 'w')

        self.Print( ["Dataset: '%s'" % (datasetName), "SaveName: '%s' " % (fileName), "MaxDecimals: '%s'" % (maxDecimals)])

        ### Define the column widths: 10 each
        template = "{0:20} & {1:20} & {2:20}" + r"\\"

        ### Create the table headers
        header = template.format("Efficiency (\%)", "Rate (kHZ)", "Ldg CaloTau E_{T} (GeV)")
        f.write( r"\begin{tabular}{ %s }" % ( " c " * 3 ) + "\n" )
        f.write( r"\hline"  + "\n" )
        f.write( header + "\n" )
        f.write( r"\hline"  + "\n" )
     
        ### Now write the actual values
        for effVal, effErrUp, effErrDown, rateVal, rateValUp, rateValDown, etThreshold in zip(effVals, effErrUp, effErrDown, rateVals, rateValsUp, rateValsDown, effEtThresholds):
            column1  = str( round( effVal*100  , maxDecimals) ) + "\pm " + str( round( effErrDown*100  , maxDecimals) )
            column2  = str( round( rateVal , maxDecimals) ) + " + " + str( round( rateValUp , maxDecimals) ) + " - " + str( round( rateValDown , maxDecimals) )
            column3  = str( etThreshold.rsplit(">=", 1)[-1] )
            row      = column1 + column2 + column3
            f.write( template.format(column1, column2, column3) + "\n" )

        f.write( r"\hline"  + "\n" )
        f.write( r"\end{tabular}" )
        f.close()
        return


    def PrintElapsedTime(self, units = "seconds", marker = ""):
        '''
        Print the time elapses since the creation of the plotter object.
        
        '''
        self.Verbose()
        
        deltaT = time.time() - self.startTime
        if units == "seconds":
            pass
        elif units == "minutes":
            deltaT = deltaT/60.0
        elif units == "hours":
            deltaT = deltaT/(60.0*60)
        else:
            raise Exception("Unsupported units of time. Please choose from 'seconds', 'minutes' and 'hours'.")
            
        self.Print(["Elapsed time%s: '%s' %s" % (marker, deltaT, units) ])
        return


    def FindFirstBinBelowValue(self, histo, targetValue, axis=1):
        '''
        FindLastBinAbove(targetValue, axis=1): 
        find last bin with content > threshold for axis (1=x, 2=y, 3=z)
        if no bins with content > threshold is found the function returns -1.
        '''
        self.Verbose()

        iBin = histo.FindLastBinAbove(targetValue, axis)
        if iBin == -1:
            self.Verbose(["WARNING! Could not find target value '%s'." % (targetValue)])
        binCenter      = histo.GetBinCenter ( iBin   )
        binCenterUp    = histo.GetBinCenter ( iBin+1 )
        binCenterDown  = histo.GetBinCenter ( iBin-1 )
        binContent     = histo.GetBinContent( iBin   )
        binContentUp   = histo.GetBinContent( iBin+1 )
        binContentDown = histo.GetBinContent( iBin-1 )
        binError       = histo.GetBinError( iBin )

        ### Sanity check
        if binContent != 0:
            percentageOffset = ((binContent - targetValue)/binContent )*100
        else:
            percentageOffset = 99999.9

        if abs(percentageOffset) > 50:
            self.Verbose( ["Target: '%s'" % (targetValue) , "BinContent: '%s'" % (binContent), "BinContent (+1): '%s'" % (binContentUp), "BinContent (-1): '%s'" % (binContentDown), 
                         "BinCenter: '%s'" % (binCenter), "BinCenter (+1): '%s'" % (binCenterUp), "BinCenter (-1): '%s'" % (binCenterDown), "Rate Offset (%%): '%f'" % (percentageOffset)] )
        return iBin, binCenter, binContent, binError


    def FindFirstBinAboveValue(self, histo, targetValue):
        '''
        '''
        self.Verbose()
        
        iBin = histo.FindFirstBinAbove(targetValue)
        if iBin == -1:
            raise Exception("Could not find target value '%s'!" % (targetValue))

        ### Get actual values
        binCenter  = histo.GetBinCenter(  iBin )
        binContent = histo.GetBinContent( iBin)

        ### Sanity check
        if binContent != 0:
            percentageOffset = ((binContent - targetValue)/binContent )*100
        else:
            percentageOffset = 99999.9

        if abs(percentageOffset) > 50:
            self.Verbose( ["Target: '%s'" % (targetValue) , "BinContent: '%s'" % (binContent), "BinCenter: '%s'" % (binCenter), "Offset (%%): '%f'" % (percentageOffset)] )
        return iBin, binCenter, binContent


    def GetTGraphValuesAsLists(self, myTGraph):
        '''
        See: 
        http://root.cern.ch/phpBB3/viewtopic.php?t=2499

        d1, d2 = ROOT.Long(0), ROOT.Long(0)
        myGraph.GetPoint( 41, d1, d2 )        
        '''
        self.Verbose()
        
        ### Get ploted array dimension
        nPoints = myTGraph.GetN()
        xVals = [None]*nPoints
        yVals = [None]*nPoints

        for p in range(0, nPoints):
            x = ROOT.Double(0) 
            y = ROOT.Double(0)                
	    myTGraph.GetPoint(p, x, y)
            xVals[p] = x
            yVals[p] = y

        ### Filter the x- and y-value lists. Keep only entries that are non-None
        xVals = filter(lambda a: a != None, xVals)
        yVals = filter(lambda a: a != None, yVals)
        return xVals, yVals        

        
    def SaveContourPointsToFile(self, tmpHistoDict):
        '''
        When option "LIST" is specified together with option "CONT", the points used to draw the contours are saved in TGraph objects:
        
        h->Draw("CONT LIST");
        gPad->Update();
    
        The contour are saved in TGraph objects once the pad is painted. Therefore to use this functionnality in a macro, gPad->Update() should be performed after the histogram drawing. 
        Once the list is built, the contours are accessible in the following way:

        TObjArray *contours = gROOT->GetListOfSpecials()->FindObject("contours")
        Int_t ncontours     = contours->GetSize();
        TList *list         = (TList*)contours->At(i);
        Where "i" is a "contour number", and "list" contains a "list of TGraph objects".     
        For one given contour, more than one disjoint polyline may be generated. The number of TGraphs per contour is given by:
        
        list->GetSize();
        
        To access the first graph in the list one should do:
        
        TGraph *gr1 = (TGraph*)list->First();       
        '''
        self.Verbose()

        if self.SaveContourPoints == False:
            return
    
        ### Variable declaration
        bWriteHeader = True

        ### x = LdgJet ET
        xVals        = []
        xErrUp       = []
        xErrDown     = []

        ### y = SubLdgJet ET
        yVals        = []
        yErrUp       = []
        yErrDown     = []
        
        ### z = Rate (kHz)
        zVals        = []
        zErrUp       = []
        zErrDown     = []
        
        ### 
        failList     = []
        failTypeList = []

        ### for-Loop: Over all histo types (Nominal, FluctuateUp, FluctuateDown) 
        for typeCounter, histoType in enumerate(tmpHistoDict):

            ### Sanity check
            hList       = tmpHistoDict[histoType]
            nHistoTypes = len(tmpHistoDict)
            if nHistoTypes != 3:
                raise Exception("Expected the list to have size 3 (Nominal, FluctuateUp, FluctuateDown) but is '%s' instead." % (nHistoTypes) )

            bSkipThisCut = False
            ### for-Loop: Over all zCutValues
            for nZCutCounter, h in enumerate(hList):

                if nZCutCounter == len(hList):
                    break

                h.Draw("cont list")
                ROOT.gPad.Update()
                ### Save Contours values by first converting them to TGraphs...
                contours  = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
                nContoursPerHisto = contours.GetSize()
                self.Verbose( ["HistoNumber: '%s'/'%s'" % (nZCutCounter, len(hList)-1), "HistoType (Counter): '%s' ('%s'/'%s')" % (histoType, typeCounter, len(tmpHistoDict)-1), 
                               "HistoName: '%s'" % (h.GetName()), "CutValue (Counter): '%s' ('%s'/'%s')" % (self.THDumbie.zCutLines[nZCutCounter], nZCutCounter, len(self.THDumbie.zCutLines)-1) ] )

                ### Access the first graph in the list
                list = contours.At(nZCutCounter)
                myTGraph = list.First()
                if not myTGraph:
                    bSkipThisCut = True
                    continue
                    
                myTGraph.SetName( h.GetName() + "_zCut" + str(self.THDumbie.zCutLines[nZCutCounter]) )
                x, y = self.GetTGraphValuesAsListsContours( myTGraph )
                zVal = []
                zErr = []
                ### In case multiple values are saved from the contour (advise against this)
                for val in range( 0, len(x) ):
                    if bSkipThisCut==True:
                        break
                    zVal.append( self.THDumbie.zCutLines[nZCutCounter] )
                    zErr.append( h.GetBinContent( h.GetXaxis().FindBin(x[val]) , h.GetYaxis().FindBin(y[val]) ) )

                if histoType == "Nominal":
                    xVals.append(x)
                    yVals.append(y)
                    zVals.append( zVal )
                elif histoType == "FluctuateUp":
                    xErrUp.append(x)
                    yErrUp.append(y)
                    zErrUp.append( zErr )
                elif histoType == "FluctuateDown":
                    xErrDown.append(x)
                    yErrDown.append(y)
                    zErrDown.append( zErr )
                else:
                    raise Exception("Unexpected histoType '%s'. Expected histoType to be one of 'Nominal', 'FluctuateUp' or 'FluctuateDown'." % ( histoType ) )

        ### Sanity check
        nValues = len(xVals)
        nDiff   = 0
        if all(len(lst) != nValues for lst in [yVals, zVals]):
            self.Verbose( ["xVals ('%s') = '%s'" % (len(xVals), xVals), "yVals ('%s') = '%s'" % (len(yVals), yVals), "zVals ('%s') = '%s'" % (len(zVals), zVals)] )
            self.Verbose( ["xErrUp ('%s') = '%s'" % (len(xErrUp), xErrUp), "yErrUp ('%s') = '%s'" % (len(yErrUp), yErrUp), "zErrUp ('%s') = '%s'" % (len(zErrUp), zErrUp)] )
            self.Verbose( ["xErrDown ('%s') = '%s'" % (len(xErrDown), xErrDown), "yErrDown ('%s') = '%s'" % (len(yErrDown), yErrDown), "zErrDown ('%s') = '%s'" % (len(zErrDown), zErrDown)])
            nDiff = ( nValues - len(yVals) ) + ( nValues - len(zVals) )
        elif all(len(lst) != nValues for lst in [xErrUp, yErrUp]):
            self.Verbose( ["xVals ('%s') = '%s'" % (len(xVals), xVals), "yVals ('%s') = '%s'" % (len(yVals), yVals), "zVals ('%s') = '%s'" % (len(zVals), zVals)] )
            self.Verbose( ["xErrUp ('%s') = '%s'" % (len(xErrUp), xErrUp), "yErrUp ('%s') = '%s'" % (len(yErrUp), yErrUp), "zErrUp ('%s') = '%s'" % (len(zErrUp), zErrUp)] )
            self.Verbose( ["xErrDown ('%s') = '%s'" % (len(xErrDown), xErrDown), "yErrDown ('%s') = '%s'" % (len(yErrDown), yErrDown), "zErrDown ('%s') = '%s'" % (len(zErrDown), zErrDown)])
            nDiff = ( nValues - len(xErrUp) ) + ( nValues - len(yErrUp) ) + ( nValues - len(zErrUp) )
        elif all(len(lst) != nValues for lst in [xErrDown, yErrDown]):
            self.Verbose( ["xVals ('%s') = '%s'" % (len(xVals), xVals), "yVals ('%s') = '%s'" % (len(yVals), yVals), "zVals ('%s') = '%s'" % (len(zVals), zVals)] )
            self.Verbose( ["xErrUp ('%s') = '%s'" % (len(xErrUp), xErrUp), "yErrUp ('%s') = '%s'" % (len(yErrUp), yErrUp), "zErrUp ('%s') = '%s'" % (len(zErrUp), zErrUp)] )
            self.Verbose( ["xErrDown ('%s') = '%s'" % (len(xErrDown), xErrDown), "yErrDown ('%s') = '%s'" % (len(yErrDown), yErrDown), "zErrDown ('%s') = '%s'" % (len(zErrDown), zErrDown)])
            nDiff = ( nValues - len(xErrDown) ) + ( nValues - len(yErrDown) ) + ( nValues - len(zErrDown) )
        else:
            pass
        

        ### Wait untill all three are filled (x, xErrUp, xErrDown  [ditto for y, yErrUp, yErrDown] )  and then save to file
        nPoints = nValues - abs(nDiff)
        for i in range( 0, nPoints, +1 ):
            
            ### Convert error values to absolute errors
            xErrUpAbs   = []
            xErrDownAbs = []
            yErrUpAbs   = []
            yErrDownAbs = []
            zErrUpAbs   = []
            zErrDownAbs = []

            for j in range( 0, len(xVals[i]), +1 ):

                xErrUpAbs  .append( round(abs( xVals[i][j] - xErrUp  [i][j] ), 2) )
                xErrDownAbs.append( round(abs( xVals[i][j] - xErrDown[i][j] ), 2) )

                yErrUpAbs  .append( round(abs( yVals[i][j] - yErrUp  [i][j] ), 2) )
                yErrDownAbs.append( round(abs( yVals[i][j] - yErrDown[i][j] ), 2) )

                zErrUpAbs  .append( round(abs( zVals[i][j] - zErrUp  [i][j] ), 2) )
                zErrDownAbs.append( round(abs( zVals[i][j] - zErrDown[i][j] ), 2) )
            
            self.Verbose( ["HistoType: '%s'" % (histoType), "CutValue (%s): '%s'" % (i, self.THDumbie.zCutLines[i]), 
                         "xVals[%s]: '%s'" % (i, xVals[i]), "xErrUp: '%s'" % (xErrUpAbs), "xErrDown: '%s'" % (xErrDownAbs), 
                         "yVals[%s]: '%s'" % (i, yVals[i]), "yErrUp: '%s'" % (yErrUpAbs), "xErrDown: '%s'" % (yErrDownAbs),
                         "zVals[%s]: '%s'" % (i, zVals[i]), "zErrUp: '%s'" % (zErrUpAbs), "zErrDown: '%s'" % (zErrDownAbs) ] )

            self.SaveTGraphArraysDataToFileWithZ( self.THDumbie.name + "_ContourValues.txt", 
                                             xVals[i], xErrUpAbs, xErrDownAbs, 
                                             yVals[i], yErrUpAbs, yErrDownAbs, 
                                             zVals[i], zErrUpAbs, zErrDownAbs, 
                                             3, bWriteHeader )
            bWriteHeader = False
        return


    def GetRateValuesList(self):
        '''
        In kHz.
        '''
        self.Verbose()

        rateVals = []
        for i in range(1, 501, 1):
            rateVals.append(i)
            
        return rateVals


    def _GetLegendHeader(self):
        '''
        Create a TLegend header and return it.
        '''
        self.Verbose()

        histo     = self.THDumbie         
        header    = ""
        nDatasets = len(self.DatasetToHistoMap.keys())
        legHeader = "empty"

        ### Create TLegend name
        if ( nDatasets == 1 ) :
            legHeader = self.DatasetObject.ConvertDatasetToLatex(histo.dataset)
        else:
            legHeader = histo.legTitle

        return legHeader


    def GetHistos(self):
        '''
        '''
        self.Verbose()
        return self.HistoObjectList


    def GetTLegend(self):
        '''
        '''
        self.Verbose()
        return self.TLegend



    def SetTLegendHeader(self, primaryHeader="", secondaryHeader=""):
        '''
        '''
        self.Verbose()
        
        if secondaryHeader=="" or secondaryHeader==None:
            self.TLegend.SetHeader( self.DatasetObject.ConvertDatasetToLatex(primaryHeader) )
        else:
            self.TLegend.SetHeader( self.DatasetObject.ConvertDatasetToLatex(primaryHeader) + ": " + secondaryHeader )
        return


    def GetTHStackHistoList(self):
        '''
        '''
        self.Verbose()
        return self.THStackHistoList


    def GetTHStack(self):
        '''
        '''
        self.Verbose()
        return self.THStack


    def GetTHDumbie(self):
        '''
        '''
        self.Verbose()
        return self.THDumbie


    def GetDatasetList(self):
        '''
        '''
        self.Verbose()
        return self.DatasetToRootFileMap.keys()


    def GetDatasetAndLatexName(self, dataset):
        '''
        '''
        self.Verbose()
        
        datasetName      = dataset.rsplit(":", 1)[-1]
        datasetLatexName = dataset.rsplit(":", 1)[0]
        return datasetName, datasetLatexName

###############################################################
### Author .........: Alexandros X. Attikis 
### Institute ......: University Of Cyprus (UCY)
### Email ..........: attikis@cern.ch
###############################################################

###############################################################
### All imported modules
###############################################################
### System modules
import os, sys
import array
import math
import copy
import inspect
from optparse import OptionParser
### Other
import ROOT
import styles as m_styles

###############################################################
### Class definition here
###############################################################
class TH1orTH2:    
    def __init__(self, path, name, legTitle, saveName, **kwargs):
        self.path            = path
        self.name            = name
        self.legTitle        = legTitle
        self.xUnits          = kwargs.get("xUnits", "")
        self.xMin            = kwargs.get("xMin", None)
        self.xMax            = kwargs.get("xMax", None)
        self.yUnits          = kwargs.get("yUnits", "")
        self.yMin            = kwargs.get("yMin", None)
        self.yMax            = kwargs.get("yMax", None)
        self.yMinRatio       = kwargs.get("yMinRatio", None)
        self.yMaxRatio       = kwargs.get("yMaxRatio", None)
        self.zUnits          = kwargs.get("zUnits", "")
        self.zMin            = kwargs.get("zMin", None)
        self.zMax            = kwargs.get("zMax", None)
        self.xLegMin         = kwargs.get("xLegMin", 0.62)
        self.xLegMax         = kwargs.get("xLegMax", 0.92)
        self.yLegMin         = kwargs.get("yLegMin", 0.82)
        self.yLegMax         = kwargs.get("yLegMax", 0.92)
        self.xCutLines       = kwargs.get("xCutLines", None)
        self.xCutBoxes       = kwargs.get("xCutBoxes", None)
        self.yCutLines       = kwargs.get("yCutLines", None)
        self.yCutBoxes       = kwargs.get("yCutBoxes", None)
        self.zCutLines       = kwargs.get("zCutLines", None)
        self.yCutLinesRatioPad = kwargs.get("yCutLinesRatioPad", True)
        self.zCutLinesErrors = kwargs.get("zCutLinesErrors", True)
        self.normaliseTo     = kwargs.get("normaliseTo", None)
        self.ratio           = kwargs.get("ratio", False)
        self.ratioErrorType  = kwargs.get("ratioErrorType", "") #"B" = Binomial
        self.invRatio        = kwargs.get("invRatio", False) #inverse ratio = 1/ratio
        self.logX            = kwargs.get("logX", False)
        self.logY            = kwargs.get("logY", False)
        self.logZ            = kwargs.get("logZ", False)
        self.logXRatio       = kwargs.get("logXRatio", False)
        self.logYRatio       = kwargs.get("logYRatio", False)
        self.gridX           = kwargs.get("gridX", False)
        self.gridY           = kwargs.get("gridY", False)
        self.binWidthX       = kwargs.get("binWidthX", None)
        self.binWidthY       = kwargs.get("binWidthY", None)
        self.verbose         = kwargs.get("verbose", False)
        self.styleType       = kwargs.get("styleType", None)
        self.saveName        = saveName
        self.dataset         = kwargs.get("dataset", None)
        self.datasetLatex    = kwargs.get("datasetLatex", None)
        self.zLabel          = kwargs.get("zLabel", "z-label")
        if self.xUnits == "" or self.xUnits == None:
            self.xLabel      = kwargs.get("xLabel", "x-label")
        else:
            self.xLabel      = kwargs.get("xLabel", "x-label") + " (" + self.xUnits + ")"

        if self.yUnits == "" or self.yUnits == None:
            self.yLabel      = kwargs.get("yLabel", "y-label")
        else:
            self.yLabel      = kwargs.get("yLabel", "y-label") + " (" + self.yUnits + ")"

        if self.zUnits == "" or self.zUnits == None:
            self.zLabel      = kwargs.get("zLabel", "z-label")
        else:
            self.zLabel      = kwargs.get("zLabel", "z-label") + " (" + self.zUnits + ")"
        self.ratioLabel      = kwargs.get("ratioLabel", None)
        self.drawOptions     = kwargs.get("drawOptions", None)
        self.legOptions      = kwargs.get("legOptions", None)
        self.scaleFactor     = +1.0
        self.rangeIntegral   = -1.0
        self.integral        = -1.0
        self.TH1orTH2        = None
        self.TFileName       = None
        self.treeVarExp      = name
        self.kwargs          = kwargs
        if self.saveName == None or self.saveName == "":
            self.saveName = self.name
        self.MsgCounter      = 0


    def Verbose(self, messageList=None):
        '''
        Custome made verbose system. Will print all messages in the messageList
        only if the verbosity boolean is set to true.
        '''
        if self.verbose == False:
            return

        self.MsgCounter = self.MsgCounter  + 1            
        print "[%s] %s:" % (self.MsgCounter, self.__class__.__name__ + "." + sys._getframe(1).f_code.co_name + "()")
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


    def SetVerbose(self, verbose):
        '''
        Manually enable/disable verbosity.
        '''
        self.Verbose()
        
        self.verbose = verbose
        self.Verbose(["Verbose mode = ", self.verbose])
        return


    def PrintAttributes(self):
        '''
        Call this function to print all histogram attributes.
        '''
        self.Print(["Attributes: %s" % (self.__dict__)])
        return


    def ApplyStyles(self, styleObject, histoObjectType):
        '''
        Takes a style object as input to customise the histogram (self). First the rebinning is done
        to the user-defined x-axis bin width. Then the fill/line/marker styles are applied. Then 
        the default axes styles are also configured.
        '''
        self.Verbose()

        self._RebinXToWidth(histoObjectType)
        self._RebinYToWidth(histoObjectType)
        self._SetFillLineMarkerStyles(styleObject, histoObjectType)
        self._SetAxisStyle(histoObjectType)
        return


    def _SetAxisStyle(self, histoObjectType):
        '''
        Sets the x- and y-axis defaults, like offset, label size, label font and yMax.
        '''
        self.Verbose()

        ### Set histogram axis labels
        self.TH1orTH2.SetTitle("")
        self.binWidthX = self.TH1orTH2.GetXaxis().GetBinWidth(0)

        if "%" not in self.yLabel:
            self.Print(["WARNING! No provision for y-units provided for '%s' in yLabel(='%s'). " % (self.TH1orTH2.GetName(), self.yLabel) ])

        if "TH1" in str(histoObjectType):
            self.TH1orTH2.GetXaxis().SetTitle( self.xLabel )
            self.TH1orTH2.GetYaxis().SetTitle( self.yLabel % (self.binWidthX) + " " + self.xUnits )
        elif "TH2" in str(histoObjectType):
            self.yLabel = self.yLabel# + " " + self.yUnits
            if "%" not in self.xLabel:
                self.Print(["WARNING! No provision for x-units provided in xLabel(='%s'). " % (self.xLabel) ])
            self.binWidthY = self.TH1orTH2.GetYaxis().GetBinWidth(0)
            self.TH1orTH2.GetXaxis().SetTitle( self.xLabel % (self.binWidthX) )
            self.TH1orTH2.GetYaxis().SetTitle( self.yLabel % (self.binWidthY) )
        else:
            raise Exception("The type of histoObject passed is not a ROOT.TH1 or a ROOT.TH2 (type = '%s')." % (type(self.TH1orTH2)) )
    
        ### Customise x- and y-axis title font, size, offset
        #for i in range(1, self.TH1orTH2.GetNbinsX()) :
        #    self.TH1orTH2.GetXaxis().SetBinLabel(i, "")
        self.TH1orTH2.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize("Z"))
        self.TH1orTH2.GetXaxis().SetTitleFont(ROOT.gStyle.GetTitleFont("Z"))
        self.TH1orTH2.GetXaxis().SetTitleOffset(1.0)
        
        #for j in range(1, self.TH1orTH2.GetNbinsY()) :
        #    self.TH1orTH2.GetYaxis().SetBinLabel(j, "")
        self.TH1orTH2.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize("Z"))
        self.TH1orTH2.GetYaxis().SetTitleFont(ROOT.gStyle.GetTitleFont("Z"))
        #self.TH1orTH2.GetYaxis().SetTitleOffset(1.35)
        self.TH1orTH2.GetYaxis().SetTitleOffset(1.40)

        ### Customise x- and y-axis label font, size, offset
        self.TH1orTH2.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Z"))
        self.TH1orTH2.GetXaxis().SetLabelFont(ROOT.gStyle.GetLabelFont("Z"))
        self.TH1orTH2.GetXaxis().SetLabelOffset(ROOT.gStyle.GetLabelOffset("Z"))
        #
        self.TH1orTH2.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Z"))
        self.TH1orTH2.GetYaxis().SetLabelFont(ROOT.gStyle.GetLabelFont("Z"))
        self.TH1orTH2.GetYaxis().SetLabelOffset(ROOT.gStyle.GetLabelOffset("Z"))
        
        ### Customise x- and y-axis label font, size, offset
        self.TH1orTH2.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Z"))
        self.TH1orTH2.GetXaxis().SetLabelFont(ROOT.gStyle.GetLabelFont("Z"))
        self.TH1orTH2.GetXaxis().SetLabelOffset(ROOT.gStyle.GetLabelOffset("Z"))
        #
        self.TH1orTH2.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Z"))
        self.TH1orTH2.GetYaxis().SetLabelFont(ROOT.gStyle.GetLabelFont("Z"))
        self.TH1orTH2.GetYaxis().SetLabelOffset(ROOT.gStyle.GetLabelOffset("Z"))
        
        xMin = None
        xMax = None
        yMin = None
        yMax = None

        ### Now for the range of the axes
        if (self.xMin != None):
            xMin = self.xMin
        else:
            self.xMin = self.TH1orTH2.GetXaxis().GetXmin()
            
        if (self.xMax != None):
            xMax = self.xMax
        else:
            self.xMax = self.TH1orTH2.GetXaxis().GetXmax()

        if (self.yMin != None):
            yMin = self.yMin 
        else:
            if "TH1" in str(histoObjectType):
                self.yMin = self.TH1orTH2.GetMinimum()
            else:
                self.yMin = self.TH1orTH2.GetYaxis().GetXmin()
                
        if (self.yMax != None):
            yMax = self.yMax
        else:
            if "TH1" in str(histoObjectType):
                self.yMax = self.TH1orTH2.GetMaximum()*self.GetYMaxFactor(self.logY)
            else:
                self.yMax = self.TH1orTH2.GetYaxis().GetXmax()

        self.TH1orTH2.GetYaxis().SetRangeUser(self.yMin, self.yMax)
        self.TH1orTH2.GetXaxis().SetRangeUser(self.xMin, self.xMax) #Only works if xMin (xMax) is greater (smaller) at the histogram creation time
        
        ### Take care of z-axis range (only applicable for ROOT.TH2's)
        if "TH2" in str(histoObjectType):
            self.TH1orTH2.GetZaxis().SetTitle( self.zLabel )
            #self.TH1orTH2.GetZaxis().SetTitleSize( self.TH1orTH2.GetZaxis().GetTitleSize()*0.8 )
            self.TH1orTH2.GetZaxis().SetTitleOffset(1.30)
            if (self.zMax != None):
                zMax = self.zMax
                self.TH1orTH2.GetZaxis().SetRangeUser(self.zMin, self.zMax) 
            else:
                pass
        else:
            return


        return


    def GetYMaxFactor(self, bLogY):
        '''
        Returns a factor with which we multiply the y-axis max to extend it.
        '''
        self.Verbose()

        yMaxFactor = None
        if bLogY == False:
            yMaxFactor = 1.4 #1.25
        else:
            yMaxFactor = 10.0 #50.0
        return yMaxFactor


    def GetYMinFactor(self, bLogY):
        '''
        Returns a factor with which we multiply the y-axis min to extend it.
        '''
        self.Verbose()

        return 1.0/self.GetYMaxFactor(bLogY)


    def _SetFillLineMarkerStyles(self, styleObject, histoObjectType):
        '''
        This function customises all the histogram-related styles (fill, marker, line). It uses a style object as input to determine all these according to 
        either the dataset name (if "styleType": None) or the actual user-defined styleType.
        '''
        self.Verbose(["For help see: 'http://root.cern.ch/root/html/TAttMarker.html' and 'http://root.cern.ch/root/html/TAttLine.html'."])

        if "TH1" in str(histoObjectType):
            (fillColour, lineColour, markerStyle, markerSize, lineWidth, lineStyle, fillStyle, drawOptions, legOptions) = styleObject.GetTH1Styles(self)
        elif "TH2" in str(histoObjectType):
            (fillColour, lineColour, markerStyle, markerSize, lineWidth, lineStyle, fillStyle, drawOptions, legOptions) = styleObject.GetTH2Styles(self)
        else:
            raise Exception("The type of histoObject passed  is not a ROOT.TH1 or a ROOT.TH2 (type = '%s')." % (type(histoObjectType)) )

        ### Apply colours/styles
        self.TH1orTH2.SetFillColor(fillColour)
        self.TH1orTH2.SetFillStyle(fillStyle)

        self.TH1orTH2.SetLineColor(lineColour)
        self.TH1orTH2.SetLineStyle(lineStyle)
        self.TH1orTH2.SetLineWidth(lineWidth)
        
        self.TH1orTH2.SetMarkerColor(fillColour)
        self.TH1orTH2.SetMarkerStyle(markerStyle)
        self.TH1orTH2.SetMarkerSize(markerSize)
        if self.drawOptions == None:
            self.drawOptions = drawOptions
        if self.legOptions == None:
            self.legOptions  = legOptions
        return


    def _RebinXToWidth(self, histoObjectType):
        '''
        Rebin a histogram x-axis according to the user-defined bin width. 
        '''
        self.Verbose()

        if self.binWidthX==None:
            return


        hName             = self.TH1orTH2.GetName()
        originalBinWidthX = self.TH1orTH2.GetXaxis().GetBinWidth(0)
        originalNBinsX    = self.TH1orTH2.GetNbinsX()

        ### Exact float comparison is tricky in python
        if ( abs(originalBinWidthX - self.binWidthX) < 1e-10): 
            self.Verbose(["Requested binWidthX '%f' is same as original bin-width size ('%f'). Doing nothing." % ( self.binWidthX, originalBinWidthX )])
            return

        self.Verbose(["Rebinning histogram '%s' of original bin size '%s'." % ( hName, originalBinWidthX )])
        ### Calculate the number of bins that correspond to the new bin width. Convert number to an integer
        xMin     = self.TH1orTH2.GetXaxis().GetXmin()
        xMax     = self.TH1orTH2.GetXaxis().GetXmax()
        nBinsX   = (xMax-xMin)/self.binWidthX
        intBinsX = int(nBinsX+0.5)

        ### Check that the user-requested binning makes sense
        if intBinsX !=0:
            remainderX = self.TH1orTH2.GetNbinsX() % intBinsX
        else:
            self.Print( ["Cannot achieve requested binning. Integer modulo by zero (intBinsX = %s). Skipping this histo." % (intBinsX)] )
            return
            
        self.Verbose(["remainderX = %s %s %s = %s" % (self.TH1orTH2.GetNbinsX(), "%", intBinsX, remainderX)])
        
        if remainderX != 0:
            self.Print(["WARNING! Trying to rebin histogram '%s' of x-axis bin-width '%s' to new bin-width of '%s'. The xMin is '%g' and xMax '%g' => number of bins would be '%g', which is not divisor of the number of bins '%d', remainder is '%d'. Will do nothing." % (hName, originalBinWidthX, self.binWidthX, xMin, xMax, nBinsX, originalNBinsX, remainderX)])
            return
            
        rebinNBinsToOne = self.TH1orTH2.GetNbinsX()/intBinsX
        if "TH1" in str(histoObjectType):
            self.TH1orTH2.Rebin(rebinNBinsToOne)
        elif "TH2" in str(histoObjectType):
            self.TH1orTH2.RebinX(rebinNBinsToOne)
        else:
            raise Exception("Something went wrong. This should not be printed.")

        ### Send a warning message if the user-defined binWidthX could not be achieved exactly.
        if self.TH1orTH2.GetXaxis().GetBinWidth(0) != self.binWidthX:
            self.Print(["WARNING! Could not achieve bin-width of '%f' for x-axis of hist '%s'. Actual bin-width is '%f'" % ( self.binWidthX,  self.name, self.TH1orTH2.GetXaxis().GetBinWidth(0))])


        return


    def _RebinYToWidth(self, histoObjectType):
        '''
        Rebin a histogram y-axis according to the user-defined bin width. 
        '''
        self.Verbose()


        if "TH1" in str(histoObjectType) == True:
            return
            
        hName             = self.TH1orTH2.GetName()
        originalBinWidthY = self.TH1orTH2.GetYaxis().GetBinWidth(0)
        originalNBinsY    = self.TH1orTH2.GetNbinsY()

        if self.binWidthY==None:
            return
        else:
            self.Verbose(["Rebinning histogram '%s' of original bin size '%s'." % ( hName, originalBinWidthY )])

        ### Calculate the number of bins that correspond to the new bin width. Convert number to an integer
        yMin     = self.TH1orTH2.GetYaxis().GetXmin()
        yMax     = self.TH1orTH2.GetYaxis().GetXmax()
        nBinsY   = (yMax-yMin)/self.binWidthY
        intBinsY = int(nBinsY+0.5)

        ### Check that the user-requested binning makes sense
        remainderY = self.TH1orTH2.GetNbinsY() % intBinsY

        self.Verbose(["remainderX = %s %s %s = %s" % (self.TH1orTH2.GetNbinsY(), "%", intBinsY, remainderY)])
        
        if remainderY != 0:
            self.Print(["WARNING! Trying to rebin histogram '%s' of y-axis bin-width '%s' to new bin-width of '%s'. The yMin is '%g' and yMax '%g' => number of bins would be '%g', which is not divisor of the number of bins '%d', remainder is '%d'. Will do nothing." % (hName, originalBinWidthY, self.binWidthY, yMin, yMax, nBinsY, originalNBinsY, remainderY)])
            return
        else:    
            rebinNBinsToOne = self.TH1orTH2.GetNbinsY()/intBinsY
            self.TH1orTH2.RebinY(rebinNBinsToOne)
        ### Send a warning message if the user-defined binWidthX could not be achieved exactly.
        if self.TH1orTH2.GetYaxis().GetBinWidth(0)!=self.binWidthY:
            self.Print(["WARNING! Could not exactly achieve a new bin-width of '%s' for y-axis. The new bin-width will instead be '%s'." % ( self.binWidthY, self.TH1orTH2.GetYaxis().GetBinWidth(0) )])

        return

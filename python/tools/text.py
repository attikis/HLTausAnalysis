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
import glob
from optparse import OptionParser
### Other
import ROOT

###############################################################
### Define class here
###############################################################
class TextClass(object):

    def __init__(self, xPos=0.0, yPos=0.0, text="", size=None, bold=False, align="left", color=ROOT.kBlack, verbose = False):
        self.verbose    = verbose
        self.MsgCounter = 0
        self.xPos       = xPos
        self.yPos       = yPos
        self.text       = text
        self.tlatex     = ROOT.TLatex()
        self.tlatex.SetNDC()
        if not bold:
            self.tlatex.SetTextFont(self.tlatex.GetTextFont()-20) # bold -> normal
        if size != None:
            self.tlatex.SetTextSize(size)
        if align.lower() == "left":
            self.tlatex.SetTextAlign(11)
        elif align.lower() == "center":
            self.tlatex.SetTextAlign(21)
        elif align.lower() == "right":
            self.tlatex.SetTextAlign(31)
        else:
            raise Exception("Error: Invalid option '%s' for text alignment! Options are: 'left', 'center', 'right'."%align)
        self.tlatex.SetTextColor(color)
        self._SetDefaults("cmsPreliminary", xPos=0.62, yPos=0.96, text = "CMS Preliminary")
        self._SetDefaults("cmsSimulation", xPos=0.62, yPos=0.96, text = "CMS Simulation")
        self._SetDefaults("cmsSimulationPhaseTwo", xPos=0.46, yPos=0.96, text = "CMS Simulation, Phase-2")
        self._SetDefaults("PU140", xPos=0.2, yPos=0.96, text = "<PU>=140")
        self._SetDefaults("cmsPreliminarySimulation", xPos=0.62, yPos=0.96, text = "CMS Preliminary Simulation")
        self._SetDefaults("cmsPreliminarySimulationPhaseTwo", xPos=0.25, yPos=0.96, text = "CMS Preliminary Simulation, Phase-2")
        self._SetDefaults("cms", xPos=0.62, yPos=0.96, text = "CMS")
        self._SetDefaults("comEnergy", xPos=0.16, yPos=0.96, text = "#sqrt{s}=14 TeV") #HPlus: x=0.19
        self._SetDefaults("intLumi", xPos=0.43, yPos=0.96, text = "L=10 fb^{-1} ")
        self._SetDefaults("COM_PU200-TP", xPos=0.63, yPos=0.96, text = "PU=200, 14 TeV")
        self._SetDefaults("COM_PU140-TP", xPos=0.63, yPos=0.96, text = "PU=140, 14 TeV")
        self._SetDefaults("COM_PU0-TP", xPos=0.63, yPos=0.96, text = "PU=0, 14 TeV")
        self._SetDefaults("cmsSimPhase2-TP", xPos=0.45, yPos=0.90, text = "#font[62]{CMS PhaseII Simulation}")  #tools/tdrstyle.py: self.rightMargin = 0.05
        self._SetDefaults("Preliminary-TP", xPos=0.71, yPos=0.86, text = "#font[52]{Preliminary}")  #tools/tdrstyle.py: self.rightMargin = 0.05
        #self._SetDefaults("cmsSimPhase2-TP", xPos=0.41, yPos=0.90, text = "#font[62]{CMS PhaseII Simulation}")  #tools/tdrstyle.py: self.rightMargin = 0.1
        #self._SetDefaults("Preliminary-TP", xPos=0.67, yPos=0.86, text = "#font[52]{Preliminary}")  #tools/tdrstyle.py: self.rightMargin = 0.1

        self.Verbose()


    def Verbose(self, messageList=None):
        '''
        Custome made verbose system. Will print all messages in the messageList
        only if the verbosity boolean is set to true.
        '''
        if self.verbose == False:
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


    def SetVerbose(self, verbose):
        self.verbose = verbose
        self.Verbose(["self.verbose = %s" % (self.verbose)])
        return


    def _SetDefaults(self, name, **kwargs):
        self.Verbose()

        ### Set all arguments and their values
        for argument, value in kwargs.iteritems():
            ### Set values to all the arguments
            setattr(self, name + "_" + argument, value)
        return


    def _GetValues(self, name, xPos, yPos, text):
        self.Verbose()

        if xPos == None:
            xPos = getattr(self, name + "_xPos" )
        if yPos == None:
            yPos = getattr(self, name + "_yPos" )
        if text == None:
            text = getattr(self, name + "_text" )
        return (xPos, yPos, text)


    def Draw(self, options=None):
        self.Verbose()

        self.tlatex.DrawLatex(self.xPos, self.yPos, str(self.text) )
        return


    def AddEnergyText(self, xPos=None, yPos=None, text=None):
        self.Verbose()

        (_x, _y, _text) = self._GetValues("comEnergy", xPos, yPos, text)
        self.AddText(_x, _y, _text)
        return


    def AddCmsPreliminaryText(self, xPos=None, yPos=None, text=None):
        self.Verbose()

        (_x, _y, _text) = self._GetValues("cmsPreliminary", xPos, yPos, text)
        self.AddText(_x, _y, _text)
        return


    def AddCmsText(self, xPos=None, yPos=None, text=None):
        self.Verbose()

        (_x, _y, _text) = self._GetValues("cms", xPos, yPos, text)
        self.AddText(_x, _y, _text)
        return


    def AddCmsSimulationText(self, xPos=None, yPos=None, text=None):
        self.Verbose()

        (_x, _y, _text) = self._GetValues("cmsSimulation", xPos, yPos, text)
        self.AddText(_x, _y, _text)
        return


    def AddCmsPreliminarySimulationText(self, xPos=None, yPos=None, text=None):
        self.Verbose()

        (_x, _y, _text) = self._GetValues("cmsPreliminarySimulation", xPos, yPos, text)
        self.AddText(_x, _y, _text)
        return


    def AddCmsSimulationPhaseTwoText(self, bIsTH2, PileUp=0, xPos=None, yPos=None, text=None):
        self.Verbose()
    
        if bIsTH2 == False:
            (_x1, _y1, _text1) = self._GetValues("cmsSimPhase2-TP", xPos, yPos, text)
        else:
            (_x1, _y1, _text1) = self._GetValues("cmsSimPhase2-TP", 0.4, 0.9, text)
            
        #(_x2, _y2, _text2) = self._GetValues("COM_PU140-TP", xPos, yPos, text) #xenios
        if str(PileUp) == "None":
            self.Print(["WARNING! PileUp variable is not a string! (%s)" % (PileUp)])
            (_x2, _y2, _text2) = self._GetValues("COM_PU140-TP", xPos, yPos, text)
        else:
            (_x2, _y2, _text2) = self._GetValues("COM_PU%0.0f-TP" %(PileUp), xPos, yPos, text)

        
        if bIsTH2 == False:
            (_x3, _y3, _text3) = self._GetValues("Preliminary-TP", xPos, yPos, text)
        else:
            (_x3, _y3, _text3) = self._GetValues("Preliminary-TP", 0.64, 0.86, text)

        self.AddText(_x1, _y1, _text1)
        self.AddText(_x2, _y2, _text2)
        self.AddText(_x3, _y3, _text3)
        return


    def AddCmsPreliminarySimulationPhaseTwoText(self, xPos=None, yPos=None, text=None):
        self.Verbose()

        (_x, _y, _text) = self._GetValues("cmsPreliminarySimulationPhaseTwo", xPos, yPos, text)
        self.AddText(_x, _y, _text)
        return


    def AddIntLumiText(self, xPos=None, yPos=None, text=None):
        self.Verbose()

        (_x, _y, _text) = self._GetValues("intLumi", xPos, yPos, text)
        self.AddText(_x, _y, _text)
        return


    def AddText(self,xPos, yPos, text, *args, **kwargs):
        self.Verbose()

        t = TextClass(xPos, yPos, text, *args, **kwargs)
        t.Draw()
        return

    def ConvertLatexToText(self, myString):
        self.Verbose()
        #return myString.replace("-", "").replace("#", "").replace("{", "").replace("}", "").replace("rightarrow", "to").replace(" ", "")
        return myString.replace("-", "").replace("#", "").replace("{", "").replace("}", "").replace("rightarrow", "to").replace(" ", "").replace("()","").replace("TMath::", "").replace("Sqrt(", "").replace("(", "").replace(")", "").replace(".", "_")

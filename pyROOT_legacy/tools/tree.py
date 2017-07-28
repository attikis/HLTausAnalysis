###############################################################
### Author .........: Alexandros X. Attikis 
### Institute ......: University Of Cyprus (UCY)
### Email ..........: attikis@cern.ch
###############################################################

###############################################################
### All imported modules
###############################################################
### System modules
import sys
import array
import math
import copy
import re
### Other 
import ROOT

###############################################################
### Define class here
###############################################################
class TreeClass(object): 
    def __init__(self, verbose=True):
        self.bVerbose          = verbose
        self.MsgCounter        = 0
        self.SupportedFormulae = ["Max$", "SecondMax$", "Min$", "Sum$"]

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


    def GetCurrentValue(self, event, branch):
        '''
        Returns the value of a given ROOT.TBranch (branch) from a given event (e.g. of a ROOT.TTree)
        '''
        self.Verbose()
        return getattr(event, branch)

    
    def TreePassFailCut(self, myTree, preSelections, mySelections):
        '''
        Determines the number of events passing a set of selections (TCuts) on a TTree. 
        The TTree could be pre-filtered with a set of pre-selections.

        About TTree::GetEntries() from: http://root.cern.ch/root/html/TTree.html#TTree:GetEntries@1
        Return the number of entries matching the selection.
        Return -1 in case of errors.

        If the selection uses any arrays or containers, we return the number
        of entries where at least one element match the selection.
        GetEntries is implemented using the selector class TSelectorEntries,
        which can be used directly (see code in TTreePlayer::GetEntries) for
        additional option.
        If SetEventList was used on the TTree or TChain, only that subset
        of entries will be considered.
        '''
        self.Verbose()
        #self.Print(["WARNING! Using TTree::GetEntries(string selection) with EXTREME caution. Read method documentation before proceeding"])

        ### Variable definition
        nPass   = -1
        nTotal  = -1

        ### Combine pre-selections and selections
        if preSelections != "":
            mySelections = preSelections + " && " + mySelections

        ### Get the number of events passing pre-selections and selections
        nTotal = myTree.GetEntries(preSelections)
        nPass  = myTree.GetEntries(mySelections)
        self.Verbose(["Pre-Selections: '%s'" % (preSelections), "Selections: '%s'" % (mySelections), "Pass: '%s'" % (nPass), "Total: '%s'" % (nTotal)])
        return nPass, nTotal


    def GetCopyOfTree(self, myTree, mySelections):
        '''
        Return a copy of a given ROOT.TTree and returns another that passes the selections passed as arguments
        '''
        self.Verbose()

        ### I must open a new ROOT file in order to copy my tree with some selections!
        output = ROOT.TFile.Open(myTree.GetName() + ".root","RECREATE");

        self.Print(["Copying TTree '%s' to a new one surviving the selections '%s'"  % ( myTree.GetName()), mySelections ])
        newTree  = myTree.CopyTree( mySelections )
        #output.Close()
        return newTree


    def LoopOverTree(self, myTree, myBranch, threshold, selections):
        '''
        Still under construction
        '''
        self.Verbose()

        nPass  = 0
        nTotal = 0

        ### I must open a new ROOT file in order to copy my tree with some selections!
        output = ROOT.TFile.Open("tmp.root","RECREATE");
        myTreeNew = myTree.CopyTree( selections )

        ### For-Loop: All events
        for event in myTreeNew:

            val = self.GetCurrentValue(event, myBranch)
        
            ### Nested for-loop: All objects in container
            for value in val:
                nTotal = nTotal + 1
                bIsAboveThreshold = value.Pt() > threshold
                if bIsAboveThreshold == True:
                    nPass = nPass + 1
                else:
                    pass

        return nPass


    def GetBranchNameAndFormulaAndFunctionFromVarExp(self, varexp):
        '''
        '''
        self.Verbose()

        reObject    = re.compile("\[.\]")
        arrayIndex     = None #Should be non-None only if the user asked for a specific element to be processes
        branchName     = None
        branchFunction = None
        branchFormula  = None
        functionName   = None

        ### De-couple branch and function from variable expression
        if "." in varexp:
            branchName     = varexp.rsplit(".", 1)[0]
            branchFunction = varexp.rsplit(".", 1)[-1].replace(" ", "")
            matchObject = re.search(reObject, branchName)
            ### Strip vector element index (if present) to get bare branch name
            if matchObject != None:
                arrayIndex = matchObject.group()
                branchName = branchName.replace(arrayIndex, "")
                arrayIndex = arrayIndex.replace("[", "").replace("]", "")
            else:
                pass
            ### Get Function name of TBranch (if any). e.g. for TLorentzVectors
            functionName = branchFunction
        else:
            branchName   = varexp
            functionName = ""

        if "$" in varexp:
            branchName     = varexp.rsplit("$", 1)[-1]
            branchFormula  = varexp.rsplit("$", 1)[0]

        ### Final touches
        if branchFunction != None:
            branchName     = branchName.replace(branchFunction, "").replace("(", "").replace(")","").replace(".", "")
            branchFunction = branchFunction.replace(")", "").replace("(", "") + "()"

        self.Verbose( ["BranchName: '%s'" % (branchName), "BranchFormula: '%s'" % (branchFormula), "BranchFunction: '%s'" % (branchFunction), "ArrayIndex: '%s'" % (arrayIndex) ] )
        return branchName, branchFormula, branchFunction, arrayIndex



    def SplitStringIntoTwoUsingCharacterList(self, myString, characterList):
        '''
        '''
        self.Verbose()

        left      = myString
        right     = ""
        character = ""
        for c in characterList:
            if c in myString:
                left      = myString.rsplit(c, 1)[0]
                right     = myString.rsplit(c, 1)[-1]
                character = c
                break
            else:
                continue
        return character, left, right


    def GetCutListFromSelectionCuts(self, selectionCuts):
        '''
        '''
        self.Verbose()
        
        ### Some definitions
        selectionCutsListRaw         = selectionCuts.split("&&")
        #selectionOptionalCutsListRaw = selectionCuts.split("||")
        logicSymbols         = ["==", "!=", "<>", ">", "<", ">=", "<="] #http://www.tutorialspoint.com/python/python_basic_operators.htm
        cutList              = []
        splitList            = []
        cutDictionary        = {}
        logicSymbol          = ""

        if selectionCuts == "" or selectionCuts == None:
            return cutList, cutDictionary


        ### Loop over all selection cuts
        for cut in selectionCutsListRaw:
            if "=<" in cut or "=>" in cut:
                raise Exception("Unsopported logic symbol provided in cut ('%s'). Please select one of the following:\n'%s'." % ( cut, logicSymbols) )
            bLogicFound = False
            cutList.append(cut)    

            ### Look for logic symbol. If found strip expression to get variable expression
            for l in logicSymbols:
                regex = re.compile(l)
                if bool(regex.search(cut)) == True:
                    bLogicFound = True
                    logicSymbol = l
                else:
                    continue

            ### Sanity check once loop over logic symbols is done
            if bLogicFound == False:
                raise Exception("Unsopported logic symbol provided ('%s'). Please select one of the following:\n'%s'." % ( logicSymbol, logicSymbols) )

            splitList = None
            formula   = None
            if any(f in cut for f in self.SupportedFormulae):
                for f in self.SupportedFormulae:
                    if f in cut:
                        formula = f
                    else:
                        pass

                #print "Formula = ", formula
                ### Find any epxression within curly braces
                m = re.search(r"\{(.*)\}", cut)
                if logicSymbol in m.group(1):
                    cutExp        = m.group(1).rsplit(logicSymbol, 1)[0]
                    cutVal        = m.group(1).rsplit(logicSymbol, 1)[1]
                    cutExpression = formula + "(" + cutExp + ")" + logicSymbol + cutVal
                    cutValue      = cut.rsplit("}", 1)[1]
                else:
                    cutExp        = m.group(1).rsplit(logicSymbol, 1)[0]
                    cutExpression = formula + "(" + cutExp + ")" + logicSymbol
                    cutValue      = cut.rsplit("}", 1)[1]
            else:
                splitList = cut.rsplit(logicSymbol, 1)
                if len(splitList) < 2:
                    raise Exception("Unexpected split list using logic symbol. Please investigate further")
                else:
                    cutExpression = splitList[0]
                    cutValue      = logicSymbol + splitList[-1]
                    cutDictionary[cutExpression] = cutValue
        
        self.Verbose(["CutExpression: '%s'" % (cutExpression), "CutValue: '%s'" % (cutValue)])
        return cutList, cutDictionary
        

    def ApplyTreeFormulaToListOfValues(self, treeFormula, valuesList):
        '''
        '''
        self.Verbose()        

        if treeFormula == "Max":
            return max(valuesList)
        elif treeFormula == "Min":
            return min(valuesList)
        elif treeFormula == "Sum":
            return len(valuesList)
        elif treeFormula == "SecondMax":
            return self.FindSubLeadingValueInList(valuesList)
        else:
            raise Exception("Unsupported formula '%s' provided. Please select one of the following options:\n\t%s" % (treeFormula, "\n\t".join(self.SupportedFormulae) ) )

    
    def FindSubLeadingValueInList(self, valuesList):
        '''
        '''
        self.Verbose()
        
        max1, max2 = None, None
        for x in valuesList:
            if x >= max1:
                max1, max2 = x, max1
            elif x > max2:
                max2 = x
        return max2


    def EnableOnlyRequiredBranches(self, myTree, myBranchList):
        '''
        if  GetEntry() is not called before setting the branch status it is (supposedly) slower
        Not setting the branch status is also (supposedly) much slower.
        '''
        self.Verbose()

        myTree.SetBranchStatus("*", False)
        for branchName in myBranchList:
            myTree.SetBranchStatus(branchName, True)
            
        print "WARNING! This does not work yet! For some reason disabling all branches causes the read branches to display wrong information. Fix it some time soon to speed up Tree loops."
        sys.exit() 
        return

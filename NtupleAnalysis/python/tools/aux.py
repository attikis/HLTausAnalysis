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
import itertools
import progressbar
from itertools import tee, izip

### ROOT
import ROOT

###############################################################
### Define class here
###############################################################
class AuxClass(object): 
    def __init__(self, verbose=False):
        self.bVerbose = verbose

    def Verbose(self, messageList=None):
        '''
        Custome made verbose system. Will print all messages in the messageList
        only if the verbosity boolean is set to true.
        '''
        if self.bVerbose == True:
            print "*** %s:" % (self.__class__.__name__ + "." + sys._getframe(1).f_code.co_name + "()")
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
        print "*** %s:" % (self.__class__.__name__ + "." + sys._getframe(1).f_code.co_name + "()")
        for message in messageList:
            print "\t", message
        return


    def GetKeyFromDictValue(self, value, dictionary):
        '''
        Return the key from a dictionary given its value
        '''
        self.Verbose()
        
        myKey = None
        for key in dictionary:
            if dictionary[key] == value:
                myKey = key
            else:
                pass
        if myKey == None:
            raise Exception( "Could not find key='%s' for value='%s' in dictionary='%s'. Make sure these exist." % (myKey, value, dictionary) )
        else:
            pass
        return myKey


    def Divide(self, numerator, denominator):
        '''
        Safely divide two numbers without worrying about division by zero.
        '''
        self.Verbose()
        
        if (denominator <= 0):
            self.Print(["WARNING! Could not proceed with the division '%s/%s'. Returning zero." % ( numerator, denominator)])
            return 0
        else:
            return float(numerator)/float(denominator)


    def Efficiency(self, nPass, nTotal, errType = "binomial"):
        '''
        Get efficiency given the number of events passing and the total number of events.
        Return efficiency and associated error.
        '''
        self.Verbose()
        
        if nTotal == 0:
            eff = 0
            err = 0
            return eff, err
        
        eff = self.Divide(nPass, nTotal)
        if errType == "binomial":
            #errSquared = nTotal*eff*(1-eff) #eff*(1-eff)/nTotal #bug?
            #err        = math.sqrt(errSquared)
            variance = eff*(1-eff)/nTotal
            err      = math.sqrt(variance)
        else:
            raise Exception( "Only the 'binomial' error type is supported at the moment.")
        self.Verbose(["eff = '%s' +/- '%s'" % (eff, err)])

        return eff, err
            

    def StartProgressBar(self, maxValue, widgets = ["", progressbar.Percentage(), " " , progressbar.RotatingMarker(),  " (" , progressbar.ETA() , ")" ]):
        '''
        widgets = [progressbar.FormatLabel(''), ' ', progressbar.Percentage(), ' ', progressbar.Bar('/'), ' ', progressbar.RotatingMarker()]
        '''
        self.Verbose()

        if maxValue <= 0:
            self.Verbose( ["WARNING: Maximum value for progress bar is '%s'. Setting to '1' to avoid crash." % (maxValue)] )
            maxValue = 1

        pBar = progressbar.ProgressBar(widgets=widgets, maxval=maxValue)
        if pBar.start_time is None:
            pBar.start()
        else:
            raise Exception( "Progress bar already started!")

        return pBar


    def UpdateProgressBar(self, pBar, currentValue):
        '''
        '''
        self.Verbose()

        pBar.update(currentValue)
        return


    def StopProgressBar(self, pBar):
        '''
        This function tells the progress bar when to stop. Call this function right after executing the potentially long function/module.
        '''
        self.Verbose()

        pBar.finish()
        return pBar


    def SaveListsToFile(self, filePath, columnWidths, columnTitles, columnValues, nDecimals=3, delimiter = "~", align = " c ", bLatexTable = True):
        '''
        '''
        self.Verbose()
        
        ### Definitions
        hLine = r"\hline" + "\n"

        ### Open file in write ("w") mode
        if os.path.exists(filePath) == True:
            self.Print( ["WARNING! File '%s' already exists!" % (filePath)] )
        f = open(filePath, "w")
        
        nColumnWidths = len(columnWidths)
        nColumnTitles = len(columnTitles)
        nColumnValues = len(columnValues)
        if( (nColumnWidths != nColumnTitles) or ( nColumnWidths!= nColumnValues) or (nColumnTitles != nColumnValues) ):
            raise Exception("Number of column-widths ('%s'), column-titles ('%s') and number of column-values ('%s') differ!" % (nColumnWidths, nColumnTitles, nColumnValues) )

        nColumns = nColumnValues
        nRows    = len(columnValues[0])
            
        ### Define the template formate
        template = ""
        for i in range(0, nColumns):
            if i != nColumns-1:
                template += "{" + str(i) + ":" + str(columnWidths[i]) + "} %s " % (delimiter)
            else:
                template += "{" + str(i) + ":" + str(columnWidths[i]) + "}"
                
        ### Create the header of the file
        header = template.format(*columnTitles)

        ### Writ eteh header
        if (bLatexTable):
            f.write(r"\begin{table}" + "\n")
            f.write(r"\centering" + "\n")
            f.write(r"\begin{tabular}{%s}" % (align * nColumns) + "\n")
            f.write(hLine)
            f.write(header + r"\\" + "\n")
            f.write(hLine)
        else:
            f.write(header + "\n")
        
        ### Loop over all rows
        for r in range(0, nRows):
            row = []

            ### Loop over all columns
            for c in range(0, nColumns):
                valuesList = columnValues[c]
                rowText    = valuesList[r]
                
                ### Round floats to a max of nDecimals
                if type(rowText) == float:
                    rowText = str( round(valuesList[r], nDecimals) )
                else:
                    rowText = str( valuesList[r] )
                
                ### Append row to a list
                row.append(rowText)

            ### Create the line for this row and format it as desired
            rowLine = template.format(*row)

            self.Verbose( ["Writing line: %s " % (rowLine) ] )
            if (bLatexTable):
                f.write(rowLine + r"\\" + "\n")
            else:
                f.write(rowLine + "\n")

        ### Close file and return
        if (bLatexTable):
            f.write(hLine)
            f.write(r"\end{tabular}" + "\n")
            f.write(r"\end{table}")
        f.close()

        return



    def ReadListsFromFile(self, filePath, delimiter = "~"):
        '''
        '''
        self.Verbose()

        if os.path.exists(filePath) == False:
            raise Exception("File '%s' does not exist! Please make sure the provided path for the file is correct." % (filePath) )

        ### Declarations
        ignoreList = ["begin{table}", "centering", "begin{tabular}", "hline", "end{tabular}", "end{table}"]
        removeList = ["\\", "\n", " "]
        lineList   = []
        nColumns   = 0

        self.Print(["FilePath: '%s'" % (filePath), "Delimiter: '%s'" % (delimiter)])
        ### Loop over all lines in myFile
        with open(filePath) as f:

            ### Loop over all lines
            for line in f:

                ### Skip latex-table related stuff
                if any(ext in line for ext in ignoreList):
                    continue
                else:
                    pass

                ### Save this line to a list (after removing unwanted characters)
                for char in removeList:
                    line = line.replace(char, "")
                lineList.append(line)

                ### Determine the number of columns
                if nColumns == 0:
                    nColumns = len(line.split(delimiter))

        ### Create list of columnLists
        listOfLists = []

        ### Loop over all columns
        for c in range(0, nColumns):
            columnList  = []
            
            ### Loop over all lines (rows)
            for l in lineList:
                columns = l.split(delimiter, nColumns)
                column  = columns[c]
                columnList.append(column)

            ### Append column to listOfLists
            listOfLists.append(columnList)

        return listOfLists



    def SaveDictOfListsToFile(self, mySavePath, fileName, dictOfLists, maxDecimals = 3):
        '''
        '''
        self.Verbose()

        self.Print(["FileName: '%s'" % (fileName), "MaxDecimals: '%s'" % (maxDecimals)])

        ### Declarations
        myFile = open( mySavePath+fileName, 'w')
        template  = ""
        iRow      = ""
        iColumn   = ""

        ### Write to file
        for iKey in dictOfLists:
            iList = dictOfLists[iKey]
            
            ### Treat differently if all elements in the list are float or strings
            myFile.write(iKey)

            if all(isinstance(item, str) for item in iList): 
                myFile.write( "".join(iList) + "\n")
            else:
                formatted = [format(iEntry, maxDecimals) for iEntry in iList]
                myFile.write(str(formatted) + "\n")
        myFile.close()
        return


    def ConvertStringToList(self, string):
        '''
        '''
        self.Verbose()
        theList = string.replace("\n", "").replace("[", "").replace("]","").replace("'","").replace(" ","").split(",")

        return theList


    def ConvertListElementsToFloat(self, myList):
        '''
        '''
        self.Verbose()
        theList = [ float(item) for item in myList ]

        return theList


    def AdjacentPairForLoop(self, iterable):
        '''
        iterate through pairs of items in python list:
        s -> (s0,s1), (s1,s2), (s2, s3), ...

        Example use:
        import aux as m_aux

        AuxObject = m_aux.AuxClass(verbose=True)
        for x1, x2 in AuxObject.AdjacentPairForLoop(xVals):
                diff = abs( x1 - x2)
        '''
        self.Verbose()
        
        a, b = tee(iterable)
        next(b, None)
        return izip(a, b)

def format(value, maxDecimals):
    if type(value)==float:
        return str( round( value , maxDecimals) )
    elif type(value)==int:
        return str( round( float(value) , maxDecimals) )
    else:
        raise Exception("ERROR! A float or int is required but got a '%s' instead with value '%s'" % (type(value), value))
        return str(value)
